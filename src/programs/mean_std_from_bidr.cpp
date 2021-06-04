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
    if(argc!=5){
      cerr << "Usage: " << command << " row col neigborhood_size infile"
      << endl;
      exit(1);
    }
    
    int row=atoi(argv[clidx++]);
    int col=atoi(argv[clidx++]);
    int n=atoi(argv[clidx++]);
    if(n%2==0 || n<1){
      cerr << command << 
	": Neighborhood_size should be an odd positive integer" << endl;
      exit(1);
    }
    n=n/2;
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
    float value=0;
    int nrows=label.fileRecords();
    int ncols=label.recordLength()/sizeof(int);
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
    double mean=0.0;
    double std=0.0;
    double min=9999999999999999.99;
    double max=-9999999999999999.99;
    int m=0;
    int pos = ifm.getPosition();  

    for(int r=row-n;r<=row+n;r++){
        if(r<=0 || r> nrows) continue;
      for(int c=col-n;c<=col+n;c++){
        if(c<=0 || c> ncols) continue;
	int offset = sizeof(float)*(c-1) + (r-1)*label.recordLength(); 
	ifm.setPosition(pos+offset);
	ifm.read(value);
	mean+=value;
	std+=value*value;
        if(min>value)min=value;
        if(max<value)max=value;
	m++;
      }
    }
    mean=mean/m;
    std=std-mean*mean*m;
    std=sqrt(std/(m-1));
    cout << mean << " " << std << " " << min << " " << max << endl;

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
