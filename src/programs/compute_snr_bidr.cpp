#include <iomanip>
#include <sstream>
#include <stdlib.h>
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
    if(argc!=4 ){
      cerr << "Usage: " << command << " s0file noises0file snrfile"
      << endl;
      exit(1);
    }
    
    string snfile=argv[clidx++];
    string nfile=argv[clidx++];
    string snrfile=argv[clidx++];

    // open files
    FileMgr snfm(snfile,"r");
    FileMgr nfm(nfile,"r");
    FileMgr snrfm(snrfile,"w");


    // read label from input files
    PDSLabel snlabel(snfm), nlabel(nfm);
    
    // compute file size info (ASSUMES FLOATING POINT FORMAT)
    int nrows=nlabel.fileRecords()-1;
    int ncols=nlabel.recordLength()/sizeof(float);
    int snrows=snlabel.fileRecords()-1;
    int sncols=snlabel.recordLength()/sizeof(float);

    // if file is not PDS
    if(snlabel.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << snfile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }

    // check for file size mismatch
    if(nrows!=snrows || ncols!=sncols){
      ErrorMessage e("Error: Mismatch between input file sizes.");
      e.throwMe();
    }

    // write label to snr file
    snrfm.write(snlabel.label());


    // close files
    snfm.close();
    nfm.close();
    snrfm.close();

    // reopen files as standard FILE* files for fast read/write
    FILE* snfp=fopen(snfile,"r");
    FILE* nfp=fopen(nfile,"r");
    FILE* snrfp=fopen(snrfile,"r+");


    
    // skip headers in BIDR (assumes 1 line in header)
    fseek(snfp,ncols*sizeof(float),SEEK_SET);
    fseek(nfp,ncols*sizeof(float),SEEK_SET);
    fseek(snrfp,ncols*sizeof(float),SEEK_SET);


    // allocate buffers
    float* sn=(float*)malloc(sizeof(float)*ncols);
    float* n=(float*)malloc(sizeof(float)*ncols);
    float* snr=(float*)malloc(sizeof(float)*ncols);

    // loop over rows 
    for(int r=0;r<nrows;r++){
  
      // read row of data from each input file
      if((int)fread(&sn[0],sizeof(float),ncols,snfp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
		       +snfile);
	e.throwMe();
      }

      if((int)fread(&n[0],sizeof(float),ncols,nfp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "+
		       nfile);
	e.throwMe();
      }

      

      // loop over columns
      for(int c=0;c<ncols;c++){
        if(sn[c]==BAD_VALUE || n[c]==BAD_VALUE) snr[c]=BAD_VALUE;
	else snr[c]= sn[c]/n[c] -1;
      } // end column loop

      // write output rows

      if((int)fwrite(&snr[0],sizeof(float),ncols,snrfp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(r+1)+" to file "
		     +snrfile);
	e.throwMe();
      }

      cout << "finished row " << r <<endl;
    } // end row loop


    // close files
    fclose(snfp);
    fclose(nfp);
    fclose(snrfp);

    // free arrays
    free(n);
    free(sn);
    free(snr);

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
