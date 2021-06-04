#include <iomanip>
#include <sstream>
#include "Error.h"
#include <string.h>
#include "PDSLabel.h"
#include "Utils.h"
#include <stdlib.h>
using std::cerr;
using std::cout;
using std::endl;


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

    //------------------------------
    // Parse Command Line
    //------------------------------
    int clidx=0;
    string command=argv[clidx++];
    if(argc!=3 && argc!=4){
      cerr << "Usage: " << command << " infile header_file [body_file]"
      << endl;
      exit(1);
    }
    
    string infile=argv[clidx++];
    string hfile=argv[clidx++];
    string bfile;
    bool output_body=false;

    if(argc==4){
      output_body=true;
      bfile=argv[clidx++];
    }

    // open input file and header output file
    FileMgr ifm(infile,"r");
    FileMgr hfm(hfile,"w");



    // read label from input file and output to header file
    PDSLabel label(ifm);

    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << infile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }

    hfm.write(label.label());

    // close header and input file
    hfm.close();
    ifm.close();



    // output records in loop
    if(output_body){

      // open input file boay file as FILE* in order to perform fast I/O
      FILE* ifp=fopen(infile,"r");
      FILE* bfp=fopen(bfile,"w");
      int records_left=label.fileRecords();
      int num_data_records=label.fileRecords()-label.labelRecords();
      int record_size=label.recordLength();
      char* buffer=(char*)malloc(record_size);
      while(records_left>0){
	// read record
        fread(buffer,1,record_size,ifp);
        
	// if body record then write
        if(records_left<=num_data_records){
	  fwrite(buffer,1,record_size,bfp);
	}
	records_left--;
      }
      free(buffer);
      fclose(ifp);
      fclose(bfp);
    }

    

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
