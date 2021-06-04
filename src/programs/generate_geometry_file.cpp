#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"
#include "L1I.h"
#include "L1B.h"
#include "Frame.h"
#include "DebugInfo.h"

using std::cout;
using std::cerr;
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

  //------------------------
  // Parse the command line
  //------------------------  

    const char* command = argv[0];
    if (argc != 2 && argc!=3 )
      {
	cerr << "Usage: " << command << " config_file old_geom_file" << endl;
	exit(-1);
      }

    int clidx = 1;
    const char* config_file = argv[clidx++];
    string cmp_file = "";
    if(argc==3) cmp_file=argv[clidx++];


    Uvar::setMode("no_unit_support");

    //------------------------------
    // Load configuration parameters
    //------------------------------

    Config cfg(config_file);
    DebugInfo::config(cfg);

    DebugInfo dbg("main");

    //---------------------------------
    // Load spice files and setup geometry handling
    //---------------------------------

    Frame::setGeometryMode(Frame::PRECOMPUTE);
    Frame::config(cfg);
    

    Frame::writeGeometryToFile();
    if(dbg.level!=0) Frame::testGeometryArrays();
    if(cmp_file!="") Frame::compareGeometryArraysToFile(cmp_file);
    Frame::cleanUp();
  }

catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
    cerr << "Non Error Message Exception caught" << endl;
  }

 return(0);
 
}
