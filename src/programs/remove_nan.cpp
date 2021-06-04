#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include "Error.h"
#include <string.h>
#include "PDSLabel.h"
#include "Utils.h"

using std::cerr;
using std::cout;
using std::endl;

#define BAD_VALUE_HEX 0xff7ffffb

using std::set_unexpected;
using std::terminate;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  terminate();
  }


int main(int argc, char* argv[])

{

set_unexpected(myUnexpected);
 string input;
try
  {
    int v=BAD_VALUE_HEX;
    float* p=(float*)&v;
    float BAD_VALUE =*p;

    //------------------------------
    // Parse Command Line
    //------------------------------
    int clidx=0;
    string command=argv[clidx++];
    if(argc!=3 ){
      cerr << "Usage: " << command << " infile outfile"
      << endl;
      exit(1);
    }
    
    
    string infile=argv[clidx++];
    string outfile=argv[clidx++];

    // open input file and header output file
    FileMgr ifm(infile,"r");

    // open output file
    FileMgr ofm(outfile,"w");

    // read label from input file
    PDSLabel label(ifm);
    
    // compute file size info (ASSUMES FLOATING POINT FORMAT)
    int nrows=label.fileRecords()-1;
    int ncols=label.recordLength()/sizeof(float);


    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << infile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }

    // write labels to output files
    ofm.write(label.label());

    // close files
    ofm.close();
    ifm.close();

    // reopen files as standard FILE* files for fast read/write
    FILE* ifp=fopen(infile,"r");
    FILE* ofp=fopen(outfile,"r+");

    
    // skip headers in BIDR (assumes 1 line in header)
    fseek(ifp,ncols*sizeof(float),SEEK_SET);
    fseek(ofp,ncols*sizeof(float),SEEK_SET);
 

     // allocate buffer for reading file
    float* value=(float*)malloc(sizeof(float)*ncols);
    

    // write first n/2 bad_value rows 
    for(int r=0;r<nrows;r++){

      if((int)fread(&value[0],sizeof(float),ncols,ifp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
		       +infile);
	e.throwMe();
      }
      
      // loop over columns
      for(int c=0;c<ncols;c++){
        if(isnan(value[c])) value[c]=BAD_VALUE;
        if(!finite(value[c])) value[c]=BAD_VALUE;
      }
      
      // write output rows

      if((int)fwrite(&value[0],sizeof(float),ncols,ofp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(r+1)+" to file "
		     +outfile);
	e.throwMe();
      }

      cout << "finished row " << r <<endl;
    } // end row loop


    // close files
    fclose(ifp);
    fclose(ofp);


    // free arrays
    free(value);


  }
catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
    cerr << "Non Error Message Exception caught" << endl;
  }
}
