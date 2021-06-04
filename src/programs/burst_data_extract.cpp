#include <string.h>
#include <strings.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <getopt.h>

#include "Plot.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"
#include "BurstData.h"
#include "L1B.h"
#include "config_keywords.h"

// Define command line option flags
#define OPTSTRING "hc:r:t:"
extern int optind;
using std::cout;
using std::cerr;
using std::endl;
using std::set_unexpected;
using std::terminate;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << std::endl;
  terminate();
  }

int main(int argc, char* argv[])

{
typedef std::istringstream ISTRINGSTREAM;

std::set_unexpected(myUnexpected);

try{

  //------------------------
  // Parse the command line
  //------------------------

  char* cfg_filename = NULL;
  double range_low;
  double range_hi;
  bool range_flag = false;
  char* type_string = NULL;
  while (1)
    {
    int c = getopt(argc,argv,OPTSTRING);
    if (c == 'h')
      {
      cout << "Option -h prints this message" << endl;
      cout << "Option -c identifies the config file to read from" << endl;
      cout << "Option -t identifies the data type" << endl;
      cout <<
        "Option -r identifies the optional parameter range for the 1st param";
        cout << endl;
      cout << "paramlist is a list of field names from the following list:";
        cout << endl;
      cout << "cnt_rl cnt_nd cnt_radio rad csr r_mode bem sclk brst" << endl;
      cout << "record_count hip cip rip bpd pul pul tro rwd pri" << endl;
      }
    else if (c == 'c')
      cfg_filename = optarg;
    else if (c == 't')
      type_string=optarg;
    else if (c == 'r')
      {
      range_flag = true;
      char c1 = ' ';
      ISTRINGSTREAM is(optarg);
      is >> range_low;
      is >> c1;
      is >> range_hi;
      if (c1 != ':' || !is.eof())
        {
        cout << "Error determining parameter range: " << optarg << endl;
        return(-1);
        }
      }
    else if (c == -1) break;
    }

    // check to make sure all necessary command line info is available
    if (cfg_filename == NULL || type_string == NULL )
    {
    cerr << "Usage: burst_data_extract -h -c cfgfile -r lo:hi -t data_type paramlist > outfile";
      cerr << endl;
    return(-1);
    }

    //------------------------------
    // Load configuration parameters
    //------------------------------

    Config cfg(cfg_filename);

    Frame::config(cfg);
    BurstData::config(cfg);

    // Set up a data object according to data type argument

    string types[]={"L1BP","L1BA"};
    string filetypes[]={"passive","active"};
    string cfg_tag[]={L1B_PASSIVE_MODE_FILENAME, 
		    L1B_ACTIVE_MODE_FILENAME};

    BurstData* data;
    bool data_type_found=false;
    for(int c=0;c<10;c++){
      if(strcasecmp(type_string,types[c].c_str())==0){
	string filename=cfg.str(cfg_tag[c]);
	if(c>1){ // L1I 
	//  data=new L1I(filename,"rb",filetypes[c]);
	}
	else{ // L1B
	  data=new L1B(filename,"rb",filetypes[c]);
	}
	data_type_found=true;	
	break;
      }
    }
    if(!data_type_found){
      cerr << "burst_data_extract: Invalid data type " << type_string << endl;
      cerr << "valid data types are:";
      for (int c=0;c<10;c++){
        cerr << types[c] << endl;
      }
      return(-1);
    }

    // Read Header

    data->readHeader();
    while(! data->eof()){
      data->extractParams(argv+optind,argc-optind,range_flag,range_low,range_hi);
    }
    return(0);
}

catch(ErrorMessage& e)
  {
  cerr << "Error: " << e.msg << std::endl;
  }
catch(...)
  {
  cerr << "Exception caught" << std::endl;
  }

return(0);
}







