#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include "Error.h"
#include <string.h>
#include "PDSLabel.h"
#include "Utils.h"

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
    if(argc!=4){
      cerr << "Usage: " << command << " row col infile"
      << endl;
      exit(1);
    }
    
    int row=atoi(argv[clidx++]);
    int col=atoi(argv[clidx++]);
    string infile=argv[clidx++];

    // open input file and header output file
    FileMgr ifm(infile,"r");


    // read label from input file
    PDSLabel label(ifm);

    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << infile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }


    // skip to correct position in file
    int value=0;
    int nrows=label.fileRecords();
    int ncols=label.recordLength()/sizeof(int);
    int headrecs=label.labelRecords();
    if(row>nrows || row < 1 ){
      cerr << command << ":Row Number " << row << " is outside of [ 1, "
	   << nrows << " ]" << endl;
      exit(1);
    }
    if(col>ncols || col < 1 ){
      cerr << command << ":Column Number " << col << " is outside of [ 1, "
	   << ncols << " ]" << endl;
      exit(1);
    }
    int offset = sizeof(int)*(col-1) + (row-1)*label.recordLength(); 
    offset+=headrecs*label.recordLength();
    ifm.setPosition(offset);
    ifm.read(value);
    cout << (int)value << endl;

    // close input file
    ifm.close();
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
