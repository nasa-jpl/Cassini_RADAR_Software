//-----------------------------------------------------------------------------
// doppler_plot
//
// Analyze doppler in a ckernel for radar use.
// This program provides usage information if invoked without any arguments.
// The -h option will provide a more detailed help message.
//
// Most of the inputs are provided in an external ascii config file which
// consists of keyword,value pairs (one per line).
//-----------------------------------------------------------------------------


#include <string.h>
#include <iostream>
#include <vector>
#include <getopt.h>
#include "Config.h"


using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::terminate;


void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  terminate();
  }

int main(int argc, char* argv[])

{

std::set_unexpected(myUnexpected);

try
  {



  //------------------------
  // Parse the command line
  //------------------------

  if(argc!=2) cerr << "Usage:"<<argv[0] << " config_file" << endl;
  
  int clidx=1;
  string cfg_filename=argv[clidx++];

  //------------------------------
  // Load config file (rewrites on output)
  //------------------------------


  Config cfg(cfg_filename);  

   
  }
catch(ErrorMessage& e)
  {
  cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
  cerr << "Exception caught" << endl;
  }

return(0);

}

