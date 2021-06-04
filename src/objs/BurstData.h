//=============================================================================
// BurstData.h
//
// This file contains the BurstData, Parameter, and ParameterList
// class declarations.  
//
// The BurstData class 
// is a base class for the L1A, L1B, L2x classes which service
// Cassini RADAR level data files. 
//
// Parameter is a class which encodes information about a data 
// object in a file
// including offset from start of data record, length in bytes, element size
// and number of elements data type and unit string (if datatype is UVAR)
//
// ParameterList is a class that contains a map of parameter names to 
// Parameters and method for accessing this map
//
// A BurstData object is designed to support either reading or writing; the
// constructor sets the type at the beginning.
//
// Interface summary:
//
// BurstData methods:
//
// Code change log:
// Y. Gim
// Date: March 17, 2005
// Add two paramters: noise_dioded_delta_tau_
//                    resistive_load_delta_tau_
//=============================================================================

#ifndef BURST_DATA_H
#define BURST_DATA_H

// FILE SIZE DEFINITIONS
#define MAX_FILE_LENGTH 2000000000
#define MAX_NUM_FILES 10

// Other definitions
#define UTC_STRING_SIZE_ISOD 21
#define UTC_STRING_PAD_ISOD 3
#define TIME_LENGTH_IN_FILE  (UTC_STRING_SIZE_ISOD+UTC_STRING_PAD_ISOD+UTC_STRING_SIZE_ISOC+UTC_STRING_PAD_ISOC +16)
#define UTC_STRING_SIZE_ISOC 23
#define UTC_STRING_PAD_ISOC 1 
#define TARGET_NAME_STRING_SIZE 16
#define TARGET_FRAME_NAME_STRING_SIZE 24
#define SYNC_VAL 0x77746B6A

#include <string>
#include <vector>
#include <map>
#include <list>
#include "Units.h"
#include "Error.h"
#include "Array.h"
#include "GPHData.h"
#include "Frame.h"
#include "Time.h"
#include "CAGTable.h"
#include "TopoMap.h"

class Parameter{
 public:
  enum DataType {ULONG, FLOAT, DOUBLE, FDATA, UVAR, UVAR_DOUBLE, FULL, TIME, 
		 STATEVEC, DIRVEC, ROTVEL, UVEC, DVEC, UMAT, FMAT, STRING,UVAR_WEST_LON};

  Parameter();
  ~Parameter();
  Parameter(string s, 
	    unsigned int off, unsigned int len, DataType type, 
	    void* output_ptr, Frame* fptr=NULL, string ustr="");
  Parameter(string s, unsigned int off, Uvar& u, string ustr);
 
  Parameter(const Parameter&);
  Parameter& operator=(const Parameter&);
  bool operator!=(const Parameter&);
  void read(FileMgr& f) const;
  void readTime(FileMgr& f) const;
  void readStateVector(FileMgr& f) const;
  void readRotationalVelocity(FileMgr& f) const;
  void readDirectionVector(FileMgr& f) const;
  void readWestLon(FileMgr& f) const;
  void show();
  unsigned int offset;
  unsigned int length;
  DataType datatype;
  void* outpt_ptr;
  bool enabled;
  string name;
  string units;  
  Frame* frame_ptr;
};

class ParameterList{

 public:
  ParameterList();
  ~ParameterList();
  void appendParameter(string s, Parameter p);
  void appendParameter(string s, Uvar& u, string ustr="");
  void appendParameter(string s, unsigned int& uli);
  void appendParameter(string s, unsigned int length, 
		       Parameter::DataType t,void* ptr,Frame* f_ptr=NULL, 
		       string ustr="");
 
  void skipBytes(unsigned length);
  bool isParam(const std::string& param_name) const; // No exceptions
  void enable(const std::string& param_name);
  void disable(const std::string& param_name);
  Parameter* getParam(const std::string& param_name);
  typedef map<string,Parameter*> PARAM_MAP;
  PARAM_MAP param_map;
  typedef list<Parameter> PARAM_LIST;
  PARAM_LIST param_list;
  unsigned int offset;
  
 
};

class BurstData
  {
  public:

  typedef Array1D<float> fdata;

  //---------------------
  // Error handling class
  //---------------------

  enum errorE {unspecified, read_error, write_error, internal_error, 
	       unknown_parameter, invalid_time};
  enum filetypeE{active, passive};
  class BurstDataError;
  
  //--------------
  // Constructors
  //--------------

  BurstData(const std::string& filename, const std::string& mode,
    const std::string& filetype);
  virtual ~BurstData();
  //---
  // configuration
  //---
  static void config(Config& cfg);

 



  //-----
  // I/O
  //-----

  virtual void readHeader();
  virtual void writeHeader();
  virtual void rewriteHeader();
  virtual void readRecord();
  void readSclkOnly(); 
  Time getStartTime();
  Time getEndTime();
  void gotoFirstRecord();
  void skipRecord(int n=1);
  int maxSabCounter(); // returns the sab_counter of the last record
  bool isSAR() const;
  bool isPassive();
  bool isCal() const;
  void createNextFile();
  void openNextFile();
  void openLastFile();
  void showPassiveDataOnScreen();


  //-----------
  // Predicates
  //-----------

  bool eof(); 
  bool goodGeometry();

  //-------------------
  // Parameter Access
  //-------------------

  void getParamNameList(list<string>& param_name_list);
  int returnParams(double* p2d, char* param_list[], int nparams,
		   bool range_flag, double p1min, double p1max ,
		   Dvec* p=NULL);
  void extractParams(char* param_list[], int nparams, bool range_flag,
		     double p1min, double p1max, Dvec* p=NULL);
  
  void enableParams(char* param_list[], bool* special, int nparams);
  virtual void enableSpecialParam(const char* name);
  void disableParams(char* param_list[], int nparams);
  virtual void disableSpecialParam(const char* name);
  double getParam(const std::string& param_name);
  Parameter* getParameter(const std::string& param_name);
  virtual void readParams();

  void readParameter(const string& param_name);
  void writeParameter(const string& param_name,
		      Uvar value, const string& ustr);
  void writeParameter(const string& param_name, const unsigned int& i);
  void writeParameter(const string& param_name, const float& x);
  unsigned int getSabCounter();
  virtual double getSpecialParam(const char* name);

  //-------------------
  // Other Methods
  //-------------------

  // This routine is public so that it can be used by the point target
  // simulator, but care should be taken in calling this routine.
  void computeMeasurementGeometry(const Umat& azim_1way3dB_ellipse_fit,
				  const Umat& elev_1way3dB_ellipse_fit,
				  const Umat& azim_2way3dB_ellipse_fit,
				  const Umat& elev_2way3dB_ellipse_fit);

  float computeCalibratedAttenuatorGain();

  unsigned int recordCount();

  unsigned int sclkToRecord(unsigned int sclk);

  void setFileNumber(int fno);

  //-----------------
  // Public Sab data
  //-----------------

  unsigned int  sync;
  unsigned int  sclk;  // actual telemetry value
  unsigned int  record_id; 
  Uvar  scpr;
  Uvar  brst;
  Uvar  header_tfi;
  unsigned int  header_tnc,header_typ,header_tca,header_tcb,header_tcc;
  unsigned int  pwri, vicc,vimc,tail_len,tail_id,sab_counter,sab_len;
  unsigned int  fswm,fswc,ctbc,ctrx,ctps,ctbe,ctps_ctbe,header_end;

  //slow field variables
  Uvar slow_tfi;
  unsigned int dtn;
  unsigned int slow_typ;
  unsigned int csr,r_mode,slow_instruction_number,bem,baq_mode;
  Uvar tro;
  Uvar rc_bw,adc;
  Uvar  at1,at3,at4;
  unsigned int at1_db,at3_db,at4_db;
  unsigned int at1_each, at3_each,at4_each;
  Uvar rip,csd;
  unsigned int rad,csq;
  Uvar chirp_length,slow_cfs;

  //fast field variables
  Uvar fast_tfi;
  unsigned int fin,fast_typ,pul,bii;
  Uvar bpd,pri,rwd,fast_csf;



  //footer1 variables  
  unsigned int  iebtth, iebttl;
  unsigned int  bgcalls,delvmn,delvda,delvyr;
  unsigned int  cnt_rl,cnt_radio,cnt_nd;
  unsigned int   eout,subr,space_craft_time;
  Uvar hip,cip;
  
  
  //sub_comm_data variables 
  Uvar fwdtmp,be1tmp,be2tmp,be3tmp,be4tmp,be5tmp;
  Uvar diptmp,rlotmp,tadcal1,nsdtmp,lnatmp,evdtmp;
  Uvar mratmp, mruttm,dcgttm,cucttm,twttmp,epctmp;
  Uvar tw1ttm,ep1ttm,p_stmp,p_sttm,fguttm,tadcal4;
  Uvar esstmp,wgb1t1,wgb3t1,wgb3t2,wgb3t3,wgb5t1;
  Uvar pcutmp,adctmp,tadcal2,ecltmp,cputmp,memtmp;
  Uvar sadctmp;
  Uvar tadcal3,frwdpw;
  Uvar dcgmon;
  Uvar lpltlm_db;
  Uvar nsdcur,hpapsm,catcur,p_smon,svlsta,usotmp;
  Uvar cpbnkv, essvlt,tadcal5,pcu5v_72,pcu5i_73,pcu5v_74;
  Uvar pcuti_75,pcu15v_76,pcu15i_77,pcu15v_78,pcu15i_79,pcu12v_80;
  Uvar pcu12i_81,pcucur,pllmon,ctu5i;
  Uvar tadcal6;
  Uvar pcu9v_86;
  Uvar pcu9i_87,pcu9v_88,pcu9i_89, tadcl7,shpttm;
  unsigned int num_bursts_in_flight;
  unsigned int Nradar_data;
  float rms_radar_data;   
 
  // geometry and time fields; computed from ephemeris
  
  Time t;                            // time

 

  Uvar T_scwg;                       // temperatures
  Uvar T_feed;
  Uvar T_hga;

  unsigned int quality_flag;               // quality flag
                                            // for burst
                                            // 0 = good
                                            // 1 = ckernel error
                                            // 255 = misc. geometry error
 
  Uvar transmit_time_offset;
  Uvar time_from_closest_approach;
  Uvar time_from_epoch;
  
  char target_name[TARGET_NAME_STRING_SIZE];
  char tbf_frame_name[TARGET_FRAME_NAME_STRING_SIZE];


  Uvar pole_right_ascension;
  Uvar pole_declination;
  Uvar target_rotation_rate;
  Uvar target_rotation_angle;

  unsigned int beam_number;     // beam number

  Uvar sc_pos_J2000_x;
  StateVector sc_state_J2000;       // spacecraft state vector
  StateVector sc_state_target;      


  DirectionVector sc_X_J2000;        // spacecraft frame definition
  DirectionVector sc_Y_J2000;
  DirectionVector sc_Z_J2000;
  DirectionVector sc_X_target;
  DirectionVector sc_Y_target;
  DirectionVector sc_Z_target;

  RotationalVelocity rot_vel_J2000;  
  RotationalVelocity rot_vel_target;

  Uvar norm_cnt_rl;
  Uvar norm_cnt_nd;
  Uvar norm_cnt_radio;


  //SBDR Science Field
  unsigned science_qual_flag;
  Uvar system_gain;
  Uvar antenna_brightness_temp;
  Uvar system_noise_temp;
  Uvar abt_std;
  Uvar pass_geom_time_offset; 
  Uvar pass_pol_angle;
  Uvar pass_emission_angle;
  Uvar pass_azimuth_angle;
  Uvar pass_centroid_lon;
  Uvar pass_centroid_lat;
  Uvar pass_major_width;
  Uvar pass_minor_width;
  Uvar pass_ellipse_pt1_lon;
  Uvar pass_ellipse_pt2_lon;
  Uvar pass_ellipse_pt3_lon;
  Uvar pass_ellipse_pt4_lon;
  Uvar pass_ellipse_pt1_lat;
  Uvar pass_ellipse_pt2_lat;
  Uvar pass_ellipse_pt3_lat;
  Uvar pass_ellipse_pt4_lat;
  unsigned int num_pulses_received;
  Uvar total_echo_energy;
  Uvar noise_echo_energy;
  Uvar x_factor;
  Uvar sigma0_uncorrected;
  Uvar sigma0_corrected;
  Uvar sigma0_uncorrected_std;
  Uvar altitude_means;
  Uvar altitude_means_std;

  Uvar act_geom_time_offset;
  Uvar act_pol_angle;
  Uvar act_incidence_angle;
  Uvar act_azimuth_angle;
  Uvar act_centroid_lon;
  Uvar act_centroid_lat;
  Uvar act_major_width;
  Uvar act_minor_width;

  Uvar act_ellipse_pt1_lon;
  Uvar act_ellipse_pt2_lon;
  Uvar act_ellipse_pt3_lon;
  Uvar act_ellipse_pt4_lon;
  Uvar act_ellipse_pt1_lat;
  Uvar act_ellipse_pt2_lat;
  Uvar act_ellipse_pt3_lat;
  Uvar act_ellipse_pt4_lat;
  Uvar altimeter_profile_range_start;
  Uvar altimeter_profile_range_step;
  unsigned int altimeter_profile_length;
  Uvar sar_azimuth_res;
  Uvar sar_range_res;
  Uvar sar_centroid_bidr_lon;
  Uvar sar_centroid_bidr_lat;


  unsigned int mw_source_id;

  // public topomap stuff
  TopoMap tmap;
  bool orthorectify;
  int ortho_max_iter;
  float ortho_tolerance;

  protected:

  //--------------------
  // static variables
  //--------------------
 
  static int radiometer_offset_; 
  static Uvar radiometer_delta_tau_;
  static Uvar noise_diode_delta_tau_;
  static Uvar resistive_load_delta_tau_;
  static bool burst_data_configured_;
  static unsigned int data_take_id_;
  static CAGTable cag_table_;
  static string target_name_;
  static Uvar target_radius_;
  static Time epoch_;
  static Time closest_approach_time_;
  static Uvar pole_right_ascension_epoch_;
  static Uvar pole_declination_epoch_;
  static Uvar target_rotation_rate_epoch_;
  static Uvar target_rotation_angle_epoch_;
  static int artificial_sclk_offset_;
  static bool use_artificial_sclk_offset_;
  static Uvar zero_range_time_delay_scatt_; 
  static Uvar zero_range_time_delay_alth_; 
  static Uvar zero_range_time_delay_sarl_; 
  static Uvar zero_range_time_delay_sarh_; 
  // static char target_name_[TARGET_NAME_STRING_SIZE];
  
  //---------------------
  // Internal I/O routines
  //---------------------
  void readPassiveSABData(FileMgr* f=NULL);
  void writePassiveSABData(int abs_record_no=-1);
  void checkEndedness(); //  called prior to read header


  void readGeometry(FileMgr* f=NULL);
  void writeGeometry();
  void readTime(FileMgr* f=NULL);
  void writeTime();
  void computeTime();
  void computeNormalizedRadiometerData();
  void computeGeometry();

  void readMeasurementGeometry(FileMgr* f=NULL);
  void writeMeasurementGeometry();
  int copy(BurstData& bd);
  int copyTo(BurstData& bd);
  void setQualityFlag();
  void editTime(const Time& t0); 
        // edits time in files and recomputes geometry
  void pretendDataRead();


  //-------------------
  // Internal variables
  //-------------------
  std::string mode_;
  std::string filename_;
  std::string base_filename_;
  FileMgr file_;
  filetypeE ft_;
  bool header_handled_;
  unsigned int record_count_;
  bool records_counted_;
  unsigned int record_size_;
  unsigned int header_size_;
  int file_number_;
 
 

  bool data_read_;
  bool T_feed_valid_;
  bool T_hga_valid_;
  bool T_scwg_valid_;

  static Uvar interpolation_valid_time_;
  static bool T_feed_available_;
  static bool T_hga_available_;
  static bool T_scwg_available_;
  

  Frame fsc_;
  static Frame j2000_;
  static Frame j2000_at_cassini_;
  

  // static internal variables 
  static  GPHData T_scwg_ephem_;
  static  GPHData T_feed_ephem_;
  static  GPHData T_hga_ephem_;
 
  static Frame ftarget_;
  static bool target_specified_;
  static bool mw_source_scan_;

  
  static Uvar pi_angle_;
  static Uvar two_pi_angle_;
  static Uvar zero_angle_;
 
  static unsigned int Num_mw_sources_;
  static vector<int> source_id_;
  static vector<Uvar> source_dec_;
  static vector<Uvar>  source_ra_;
  static vector<string> source_name_;

  // container : for simplicity of east and west conversion
  vector<Uvar> east_west_;
  //-----------------------
  // constant class  data 
  //-----------------------

  // Parameter name to variables (in *this) mapping.
  ParameterList parameters_;

  };

//------------------------------------------------------
// Error handling class for BurstData (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class BurstData::BurstDataError : public ErrorMessage
  {
  public:

  // Constructors

    BurstDataError(errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified BurstData Error";
    else if (error_type == read_error)
      msg = "BurstData read error";
    else if (error_type == write_error)
      msg = "BurstData write error";
    else if (error_type == internal_error)
      msg = "BurstData internal error";
    else if (error_type == unknown_parameter)
      msg = "BurstData unknown parameter name";
    }

  BurstDataError(const std::string& emsg, errorE err_type = unspecified) 
    // no exceptions
    : error_type(err_type)
    {
    msg = emsg;
    }
  void throwMe(){throw *this;}
  // Public type flag
  BurstData::errorE error_type;
  };

#endif









