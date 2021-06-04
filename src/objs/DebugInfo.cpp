#include"DebugInfo.h"
#include "Error.h"
#include "Frame.h"
#include "Utils.h"

using std::cout;
using std::cerr;
using std::endl;
using std::hex;
using std::dec;

// initialize static members
DebugInfo::integers DebugInfo::level_map_;
ofstream DebugInfo::file;
int DebugInfo::global_level;
bool DebugInfo::allWarnings=true;

DebugInfo::DebugInfo(const string& s)
{

  integers::const_iterator p=level_map_.find(s);
  if(p==level_map_.end()){
    level=global_level; 
      // if not in list then use global_level
  }
  else{
    level=p->second;
  }
  
}
 
void DebugInfo::setLevel(const string& routine_name, int l){
  level_map_[routine_name]=l;
}

void DebugInfo::config(Config& cfg)
{
    unsigned int i = 1;
    global_level=0;
    string filename=cfg.str("debug_filename");
    allWarnings=(bool)cfg.getInt("ENABLE_ALL_WARNINGS");
    file.open(filename.c_str());
    if(!file){
      ErrorMessage e("Cannot create debug file " + filename);
      e.throwMe();
    }
    while (1)
      {
      if (cfg.keywordExists("debug_routine_name_" + toStr(i)))
        {
	// read ith debug routine name and debug level for that routine
        string rname = cfg.str("debug_routine_name_" + toStr(i));
        int  l = cfg.getInt("debug_level_" + toStr(i));
        if(rname == "GLOBAL" || rname=="Global" || rname=="global"){
	  global_level=l;
	}
        else{
         level_map_[rname]=l;
	}
        }
      else
        {
        break;  // only load a sequential set of names
        }
      ++i;
      }

}

void
DebugInfo::report( ofstream& fs){
  DebugInfo dbgsp("Special_Calibration_Debug");
  if(dbgsp.level){
    fs << endl << "Special_Calibration_Debug Enabled. Flag=" << dbgsp.level
       << "(" << (int*) dbgsp.level << ")" << endl;

    fs << "This debugging option sets several calibration parameters in the";
    fs << endl << "sar processor and point target simulator to unity." 
       << endl;
    fs << "Which parameters are disabled (set to unity) is determined by the"
       << endl;
    fs << "debug level (flag) for the routine name, Special_Calibration_Debug" << endl;
    fs << "See comments in PointTargetSim::addEcho or L1I::calibrate for details." << endl << endl;

    fs << "Parameters settings are as follows: 1=disabled, 0=enabled" 
       << endl << endl;

    int bit[19];
    int flag=dbgsp.level;
    for(int i=0;i<19;i++){
      bit[i]=flag & 0x00001;
      flag = flag >> 1;
    }
    fs << "L1I:: calibrate: ";
    fs << "Xcal(" << bit[0] << ") ";
    fs << "g2(" << bit[1] << ") ";
    if(bit[1]) fs << " [Set to 1 inside 3 dB contour, 0 outside]" << endl;
    fs << "area(" << bit[2] << ") ";
    fs << "range^4(" << bit[3] << ") ";
    fs << "cag(" << bit[4] << ") ";
    fs << "norm_value(" << bit[5] << ") " << endl;

    fs << "SARProcParams:: ";
    fs << "Xconstant_(" << bit[6] << ") ";
    fs << "Xbeam_(" << bit[7] << ") ";
    fs << "Xtimevar_(" << bit[8] << ") ";
    fs << "Xmode_(" << bit[9] << ") " << endl;
    fs << "PointTargetSim::addEcho: ";
    fs << "scale(" << bit[10] << ") ";
    fs << "area_correction(" << bit[11] << ") ";
    fs << "sigma0_value(" << bit[12] << ") "<<endl;
    fs << "PointTargetSim::convertToByte: ";
    fs << "gain(" << bit[13] << ") "<<endl;
    fs << "PointTargetSim::getDelayDopplerAndScale: ";
    fs << "target_ampl(" << bit[14] << ") ";
    fs << "Xconstant(" << bit[15] << ") ";
    fs << "g2(" << bit[16] << ") ";
    if(bit[16]) fs << " [Set to 1 inside 3 dB contour, 0 outside]" << endl;
    fs << "range^4(" << bit[17] << ") ";
    fs << "2*sar_cal_coeff(" << bit[18] << ") " << endl;
  }

  DebugInfo dbgsp2("Special_Radar_Parameters_Debug");
  if(dbgsp2.level){
    fs << endl << "Special_Radar_Parameter_Debug Enabled. Flag=" 
       << dbgsp2.level
       << "(" << hex <<  dbgsp2.level << ")" << endl;

    fs << dec;
    fs << "This debugging option sets several radar parameters in the";
    fs << endl << "point target simulator to config file values rather than" 
       << endl;
    fs << "IEB parameters when the simulator is run in IEB mode."
       << endl;
    fs << "The debug level (flag) for the routine name, Special_Radar_Parameters_Debug" << endl;
    fs << "See comments in PointTargetSim::specialRadarParamsDebug." 
       << endl << endl;

    fs << "Parameters settings are as follows: 1=from config, 0= from usual source" 
       << endl << endl;

    int bit[32];
    int flag=dbgsp.level;
    for(int i=0;i<1;i++){
      bit[i]=flag & 0x00001;
      flag = flag >> 1;
    }
    fs << "PulseDutyCycle(" << bit[0] << ") " << endl;
  }
}

void
DebugInfo::report(){
  report(file);
}
