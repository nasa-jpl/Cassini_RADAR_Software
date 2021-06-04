#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <getopt.h>
#include "Array.h"
#include "Plot.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "RadiometerData.h"
#include "Frame.h"
#include "Time.h"
#include "TargetGeom.h"
#include "Utils.h"
#include "config_keywords.h"
#include "Constants.h"
#include "Flyby.h"
#include "SARProcParams.h"
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::terminate;
using std::ofstream;
using std::setw;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  terminate();
  }
int main(int argc, char* argv[]){

std::set_unexpected(myUnexpected);

try
  {
    if(argc!=2){
      cerr << "getClosestApproach config_file > outputfile" << endl;
      exit(1);
    }
    set_unit_option("auto");
    Config cfg(argv[1]);
    Frame::config(cfg);
    SARProcParams spp(cfg,1);
    Flyby flyby(cfg);
    Time epoch_time =flyby.epochTime();
    StateVector s=spp.getClosestApproachState();
    PositionVector pos=s.position();
    FloatVector vel=s.velocity();
    double time_in_s;
    string utc;
    epoch_time.getEt(time_in_s);
    utc=epoch_time.utc("ISOD");
    cout.precision(20);
    cout << time_in_s << " Seconds since J2000, "<<utc<<endl;
    cout << "Position Vector = " << pos << endl;
    cout << "Velocity Vector = " << vel << endl;
  }
catch(ErrorMessage& e)
  {
  cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
  cerr << "Exception caught" << endl;
  }

}
