//--------------------------------------
// BurstData , Parameter, and ParameterList class routines
//--------------------------------------

#include "BurstData.h"
#include "Error.h"
#include "Constants.h"
#include "math.h"
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include "config_keywords.h"
#include "TargetGeom.h"
#include "Frame.h"
#include "Flyby.h"
using std::cout;
using std::cerr;
using std::endl;
using std::cin;

//--------------------------//
// Methods for Parameter    //
//--------------------------//

Parameter::Parameter()
  : offset(0), length(0), datatype(ULONG), outpt_ptr(NULL),
    enabled(false),frame_ptr(NULL)
{
}
Parameter::Parameter(string s, 
	    unsigned int off, unsigned int len, DataType type, 
	    void* output_ptr, Frame* fptr, string ustr)
  : offset(off), length(len), datatype(type), outpt_ptr(output_ptr),
    enabled(false),name(s),units(ustr),frame_ptr(fptr)
{
}

Parameter::Parameter(string s, unsigned int off, Uvar& u, string ustr) 
  : offset(off), length(sizeof(float)), datatype(UVAR), outpt_ptr((void*)&u),
    enabled(false),name(s),units(ustr),frame_ptr(NULL)
{
}


Parameter::~Parameter(){}
Parameter::Parameter(const Parameter& p)
{
  *this=p;
}

Parameter& 
Parameter::operator=(const Parameter& p){
  offset=p.offset;
  length=p.length;
  datatype=p.datatype;
  outpt_ptr=p.outpt_ptr;
  enabled=p.enabled;
  name=p.name;
  units=p.units;
  frame_ptr=p.frame_ptr;
  return(*this);
}

bool Parameter::operator!=(const Parameter& p2)  
  {
  if (offset != p2.offset || length != p2.length || datatype != p2.datatype)
    {
    return(true);
    }
  if (datatype == ULONG)
    {
    return(*((unsigned int*)outpt_ptr) != *((unsigned int*)p2.outpt_ptr));
    }
  else if (datatype == FLOAT)
    {
    return(*((float*)outpt_ptr) != *((float*)p2.outpt_ptr));
    }
  else if (datatype == DOUBLE)
    {
    return(*((double*)outpt_ptr) != *((double*)p2.outpt_ptr));
    }
  else if (datatype == UVAR || datatype == UVAR_DOUBLE || datatype == TIME)
    {
    return(*((Uvar*)outpt_ptr) != *((Uvar*)p2.outpt_ptr));
    }
  else
    {
    return(false);
    }
    
  }

void
Parameter::read(FileMgr& f) const {
  unsigned int start_position=f.getPosition();
  if(datatype==ULONG){
    unsigned int* ptr=(unsigned int*)outpt_ptr;
    f.read(ptr[0]);
  }
  else if (datatype==FLOAT){
    float* ptr=(float*) outpt_ptr;
    f.read(ptr[0]);
  }
  else if (datatype==DOUBLE){
     double* ptr=(double*) outpt_ptr;
    f.read(ptr[0]);
  }  
  else if (datatype==UVAR){
    float value;
    f.read(value);
    Uvar* ptr=(Uvar*) outpt_ptr;
    *ptr=coerce_base(value,units);
  }

  else if (datatype==UVAR_DOUBLE){
    double value;
    f.read(value);
    Uvar* ptr=(Uvar*) outpt_ptr;
    *ptr=coerce_base(value,units);
  }

  else if (datatype==FDATA){
    BurstData::fdata* ptr= (BurstData::fdata*) outpt_ptr;
    unsigned int nelem= length/sizeof(float);
    for (unsigned int i = 0; i < nelem; i++){f.read((*ptr)(i));}
  }
  else if (datatype==TIME){
    readTime(f);
  }
  else if(datatype==STRING){
    char* ptr=(char*) outpt_ptr;
    ptr=new char[length];
    for(unsigned int c=0;c<length;++c) f.read(ptr[c]);
  }
  else if (datatype==STATEVEC){
    readStateVector(f);
  }
  else if (datatype==DIRVEC){
    readDirectionVector(f);
  }
  else if (datatype==ROTVEL){
    readRotationalVelocity(f);
  }
  else if (datatype==UVAR_WEST_LON){
    readWestLon(f);
  }
  else{
    BurstData::BurstDataError e("Parameter DataType Not Yet Available");
    e.throwMe();
  }
  if(length!=f.getPosition()-start_position){
    
    BurstData::BurstDataError e("Parameter Length for " + name + " is incorrect");
    e.throwMe();
  }
}

void Parameter::readTime(FileMgr& f) const
{
  double t_sclk;
  f.read(t_sclk);
  Uvar t_et;
  t_et.read(f,0,"s");
  Time* t=(Time*) outpt_ptr;
  t->setEt(t_et);
  string utc_isoc;
  utc_isoc.resize(UTC_STRING_SIZE_ISOC+UTC_STRING_PAD_ISOC);
  f.read(utc_isoc);
  string utc_isod;
  utc_isod.resize(UTC_STRING_SIZE_ISOD+UTC_STRING_PAD_ISOD);
  f.read(utc_isod);
}

void Parameter::readStateVector(FileMgr& f) const
{
  double x,y,z;
  f.readXYZ(x,y,z);
  PositionVector scpos("scpos",*frame_ptr,Time(0.0),Uvar(x,"km"),Uvar(y,"km"),
		       Uvar(z,"km"));
  f.readXYZ(x,y,z);
  FloatVector scvel("scvel",*frame_ptr,Time(0.0),Uvar(x,"km/s"),Uvar(y,"km/s"),
		       Uvar(z,"km/s"));
  StateVector* s=(StateVector*) outpt_ptr;
  s->setState(scpos,scvel);

}

void Parameter::readRotationalVelocity(FileMgr& f) const
{
  double x,y,z;
  f.readXYZ(x,y,z);
  RotationalVelocity* r=(RotationalVelocity*) outpt_ptr;
  *r=RotationalVelocity("",*frame_ptr,Time(0.0),Uvar(x,units),
			Uvar(x,units), Uvar(z,units));
}

void Parameter::readDirectionVector(FileMgr& f) const
{
  double x,y,z;
  f.readXYZ(x,y,z);
  DirectionVector* d=(DirectionVector*) outpt_ptr;
  *d=DirectionVector("",*frame_ptr,Time(0.0),x,y,z);
}

void Parameter::readWestLon(FileMgr& f) const
{
  float value;
  f.read(value);
  if(value >=0.0 && value<=180.0) value = - value;
  else if( value >=180.0 && value <=360.0) value = 360 - value;
  else {
    BurstData::BurstDataError e("Longitude value is out of range" + toStr(value));
    e.throwMe();
  }
  Uvar* ptr=(Uvar*) outpt_ptr;
  *ptr=coerce_base(value,units);
}

//----------------------------------------------------
//
// show()
//
// Print out the parameter according to its type.
//----------------------------------------------------

void Parameter::show()
  {
  if (datatype == ULONG)
    {
    cout << *((unsigned int*)outpt_ptr);
    }
  else if (datatype == FLOAT)
    {
    cout << *((float*)outpt_ptr);
    }
  else if (datatype == DOUBLE)
    {
    cout << *((double*)outpt_ptr);
    }
  else if (datatype == UVAR || datatype == UVAR_DOUBLE || datatype == TIME)
    {
    cout << *((Uvar*)outpt_ptr);
    }
  else
    {
    cout << "unknown parameter type";
    }
  }

//--------------------------//
// Methods for ParameterList //
//--------------------------//

//----------------------------------------------------------------------
// isParam
//
// Check to see if supplied string matches any of the data field names.
//----------------------------------------------------------------------

bool ParameterList::isParam(const std::string& param_name) const // no exceptions
  {
  PARAM_MAP::const_iterator p = param_map.find(param_name);
  return(p != param_map.end());
  }
Parameter* ParameterList::getParam(const std::string& param_name) 
{
  if(!isParam(param_name)){
    BurstData::BurstDataError e("Could not find parameter "+param_name,
				BurstData::unknown_parameter);
    e.throwMe();
  }
  return(param_map[param_name]);
}
ParameterList::ParameterList()
 :offset(0)
{}

ParameterList::~ParameterList()
{}

void
ParameterList::appendParameter(string s, Parameter p){
  param_list.push_back(p);
  param_map[s]=&(param_list.back());
}

void ParameterList::skipBytes(unsigned int length)
{
  offset+=length;
}

void
ParameterList::appendParameter(string s, unsigned int& uli)
{
  Parameter p(s,offset,sizeof(unsigned int), Parameter::ULONG,
	      (void*)&uli);
  offset+=p.length;
  appendParameter(s,p);
}

void
ParameterList::appendParameter(string s, Uvar& u, string ustr)
{
  Parameter p(s,offset,u,ustr);
  offset+=p.length;
  appendParameter(s,p);
}

void
ParameterList::appendParameter(string s, unsigned int length,
			       Parameter::DataType t, void* ptr, Frame* f_ptr,
			       string ustr)
{
  Parameter p(s,offset,length,t,ptr,f_ptr,ustr);
  offset+=p.length;
  appendParameter(s,p);
}



void
ParameterList::enable(const std::string& param_name){
  if(!isParam(param_name)){
    BurstData::BurstDataError e("Could not find parameter "+param_name,BurstData::unknown_parameter);
    e.throwMe();
  }
  param_map[param_name]->enabled=true;
}

void
ParameterList::disable(const std::string& param_name){
  if(!isParam(param_name)){
    BurstData::BurstDataError e("Could not find parameter "+param_name,BurstData::unknown_parameter);
    e.throwMe();
  }
  param_map[param_name]->enabled=false;
}
//-----------------
// Methods for BurstData
//-----------------

// static member initialization
int BurstData::radiometer_offset_; 
Uvar BurstData::radiometer_delta_tau_;
Uvar BurstData::noise_diode_delta_tau_;
Uvar BurstData::resistive_load_delta_tau_;
bool BurstData::burst_data_configured_=false;
unsigned int BurstData::data_take_id_;
Uvar  BurstData::interpolation_valid_time_;
GPHData BurstData::T_scwg_ephem_;
GPHData BurstData::T_feed_ephem_;
GPHData BurstData::T_hga_ephem_;
Frame BurstData::ftarget_;
Frame BurstData::j2000_;
Frame BurstData::j2000_at_cassini_;
bool BurstData::target_specified_=false;
bool BurstData::mw_source_scan_=false;
bool BurstData::T_feed_available_=false;
bool BurstData::T_hga_available_=false;
bool BurstData::T_scwg_available_=false;
CAGTable BurstData::cag_table_;
string BurstData::target_name_;
Uvar BurstData::target_radius_;
Time  BurstData::epoch_;
Time  BurstData::closest_approach_time_;
Uvar  BurstData::pole_right_ascension_epoch_;
Uvar  BurstData::pole_declination_epoch_;
Uvar  BurstData::target_rotation_rate_epoch_;
Uvar  BurstData::target_rotation_angle_epoch_;
Uvar BurstData::pi_angle_;
Uvar BurstData::two_pi_angle_;
Uvar BurstData::zero_angle_;
unsigned int BurstData::Num_mw_sources_;
int BurstData::artificial_sclk_offset_;
bool BurstData::use_artificial_sclk_offset_;
Uvar BurstData::zero_range_time_delay_scatt_;
Uvar BurstData::zero_range_time_delay_alth_;
Uvar BurstData::zero_range_time_delay_sarl_;
Uvar BurstData::zero_range_time_delay_sarh_;
vector<int> BurstData::source_id_;
vector<Uvar> BurstData::source_dec_;
vector<Uvar> BurstData::source_ra_;
vector<string> BurstData::source_name_;
//-------------------
// Virtual methods
//-------------------
BurstData::~BurstData(){}

void BurstData::readHeader(){
  BurstDataError("BurstData::readHeader should Not be called").throwMe();
}

void BurstData::writeHeader(){
  BurstDataError("BurstData::writeHeader should Not be called").throwMe();
}

void BurstData::rewriteHeader(){
  BurstDataError("BurstData::rewriteHeader should Not be called").throwMe();
}

void BurstData::readRecord(){
  BurstDataError("BurstData::readRecord should Not be called").throwMe();
}

//-------------
// configuration
//-------------

void
BurstData::config(Config& cfg)
{
  radiometer_delta_tau_=cfg[RADIOMETER_DELTA_TAU_KEYWORD];
  radiometer_offset_=cfg.getInt(RADIOMETER_OFFSET_KEYWORD);
  noise_diode_delta_tau_=cfg["noise_diode_delta_tau"];
  resistive_load_delta_tau_=cfg["resistive_load_delta_tau"];
  data_take_id_=cfg.getInt("data_take_number");
  if(data_take_id_ < 32) {
    cout<<"------------------------------"<<endl;
    cout<<"Data take number in the config file will be used for LBDR/SBDR "<<endl;
    cout<<"Data take number  is "<< data_take_id_<<endl;
    cout<<"------------------------------"<<endl;
  }
  if(cfg.keywordExists("Reset_data_take_number")){
    cout<<"------------------------------------ "<<endl;
    cout<<"Will reset data take number to the one "<<endl;
    cout<<" specified in the config file "<<endl;
    data_take_id_=cfg.getInt("Reset_data_take_number");
    cout<<"New data take number "<< data_take_id_<<endl;
  }
  string  T_scwg_filename = cfg.str(T_SCWG_FILENAME);
  string  T_feed_filename = cfg.str(T_FEED_FILENAME);
  string  T_hga_filename = cfg.str(T_HGA_FILENAME);

  //This routine checks the existence of temperature telemetry files
  if(T_scwg_filename!="NONE"){
     T_scwg_ephem_.read(T_scwg_filename);
     T_scwg_available_=true;
  }
  if(T_feed_filename!="NONE"){
    T_feed_ephem_.read(T_feed_filename);
    T_feed_available_=true;
  }
  if(T_hga_filename!="NONE"){
    T_hga_ephem_.read(T_hga_filename);
    T_hga_available_=true;
  }

  //-----------------------
  //basic angle declaration
  //------------------------
  pi_angle_ =Uvar(pi,"rad");
  two_pi_angle_ = Uvar(2.0*pi,"rad");
  zero_angle_ =Uvar(0,"rad");


  //--------------
  //J2000 centered at cassini
  // j2000_: its center depends on target
  // j2000_at_cassini_: target independent origin
  //---------------
  j2000_at_cassini_=Frame("J2000","Cassini");

  //---------------
  //target name : planetary target, microwave source, or no-target
  //----------------
  target_name_= cfg.str("target"); 
  if(strcasecmp(target_name_.c_str(),"none")!=0 &&
     strcasecmp(target_name_.c_str(),"source")!=0){
    ftarget_=Frame(default_target_frame_spice_id,default_target_spice_id);
    j2000_=Frame(j2000_frame_spice_id,default_target_spice_id);
   
    
    //----
    //1.  try to compute epoch time
    //----
    Flyby flyby(cfg);
    try{
      epoch_=flyby.epochTime();
    }
    catch(ErrorMessage& e){
      cout<<"Error message "<<e.msg<<endl;
      cout<<"failed to obtain epoch time from flyby object  " <<endl;
      ErrorMessage(e).throwMe();//program stops here
    }
   
    //---------------
    //2. compute closest approach time
    // Sun: special treatment
    //---------------
    if(strcasecmp(target_name_.c_str(),"sun")!=0){
      try{
	closest_approach_time_ = flyby.lowestAltitudeTime();
      }
      catch(ErrorMessage& e){
	cout<<"Error message "<<endl;
	cout<<"failed to compute closest approach time "<<endl;
	ErrorMessage(e).throwMe();//program stops here
      }
    }
    else{
      //exception: when target is Sun, do not compute closest approach time
      cout<<endl;
      cout<<"----- Warning ---------------------"<<endl;
      cout<<"Because target is Sun, "<<endl;
      cout<<" closest approach time is replaced with epoch time "<<endl;
      cout<<"------ Work only for C27_b cruise sequence "<<endl;
      cout<<endl;
      closest_approach_time_ = epoch_;
    }
    
    getTBC_Frame_RA_DEC_ROTATION_in_J2000(epoch_,
					  target_name_,
					  pole_right_ascension_epoch_,
					  pole_declination_epoch_,
					  target_rotation_angle_epoch_);
    target_rotation_rate_epoch_ 
      = getTBC_Frame_Rotation_Rate_in_J2000(epoch_, target_name_,Uvar(10,"s"));
    TargetGeom tg;
    tg.setTarget(target_name_);
    target_radius_ = tg.radius();
    target_specified_=true;
    mw_source_scan_ = false;
  }
  //------------
  //microwave source scan
  //--------------
  else if( strcasecmp(target_name_.c_str(),"source")==0){
    string epoch_string=cfg.str("epoch_time");
    epoch_.setUtc(epoch_string);    
    closest_approach_time_=epoch_;
    j2000_=Frame(j2000_frame_spice_id,earth_spice_id);//centered at cassini
    
    //---------------
    //Read microwave sources' RA and DEC
    //---------------
    Num_mw_sources_ = (unsigned int) cfg.getInt("NUM_SOURCES");
    if(Num_mw_sources_==0) ErrorMessage("No number of microwave sources is specified").throwMe();
    for(unsigned int i=0;i<Num_mw_sources_;i++)
      {
	source_id_.push_back(cfg.getInt("SOURCE_ID_"+toStr(i+1)));
	source_dec_.push_back(cfg["SOURCE_DEC_"+toStr(i+1)]);
	source_ra_.push_back(cfg["SOURCE_RA_"+toStr(i+1)]);
	source_name_.push_back(cfg.str("SOURCE_NAME_"+toStr(i+1)));
      }
    target_specified_=false;
    mw_source_scan_=true;//microwave source scan
    target_name_="none";//do not use source as target name 
  }
  //--------------
  //no target
  //---------------
  else{
    //when target is not set, use epoch and closest_approach_time as one
    // specified in config file 
    string epoch_string=cfg.str("epoch_time");
    epoch_.setUtc(epoch_string);    
    closest_approach_time_=epoch_;
    j2000_=Frame(j2000_frame_spice_id,earth_spice_id);
    target_specified_=false;
    mw_source_scan_=false;
  }

  // Configure artificial time offset for use with C43/TA comparisons
  // Aggressively Warn when this keyword is set.
  if(cfg.keywordExists("TA_C43_ARTIFICIAL_SCLK_OFFSET")){
    cerr << "Warning TA/C43 Artificial time offset in use." << endl;
    cerr << "For normal operations YOU MUST remove the " << endl <<
      "TA_C43_ARTIFICIAL_SCLK_OFFSET keyword from the config file" << endl;
    artificial_sclk_offset_=cfg.getInt("TA_C43_ARTIFICIAL_SCLK_OFFSET");
    use_artificial_sclk_offset_=true;
  }
  else{
    artificial_sclk_offset_=0;
    use_artificial_sclk_offset_=false;
  }

  //---------------------------
  //read zero-range time offset
  //-----------------------------
  zero_range_time_delay_scatt_=cfg["zero_range_time_delay_scatt"];
  if(zero_range_time_delay_scatt_ < Uvar(0,"s")) ErrorMessage("invalid delay, should be positive").throwMe();
  zero_range_time_delay_scatt_ *=-1.0;//positive delay will effectively reduce rwd
  
  zero_range_time_delay_alth_=cfg["zero_range_time_delay_alth"];
  if(zero_range_time_delay_alth_ < Uvar(0,"s")) ErrorMessage("invalid delay, should be positive").throwMe();
  zero_range_time_delay_alth_ *=-1.0;//positive delay will effectively reduce rwd

  zero_range_time_delay_sarl_=cfg["zero_range_time_delay_sarl"];
  if(zero_range_time_delay_sarl_ < Uvar(0,"s")) ErrorMessage("invalid delay, should be positive").throwMe();
  zero_range_time_delay_sarl_ *=-1.0;//positive delay will effectively reduce rwd
  
  zero_range_time_delay_sarh_=cfg["zero_range_time_delay_sarh"];
  if(zero_range_time_delay_sarh_ < Uvar(0,"s")) ErrorMessage("invalid delay, should be positive").throwMe();
  zero_range_time_delay_sarh_ *=-1.0;//positive delay will effectively reduce rwd


 



  //-------------------
  //Interpolation valid time
  //-------------------
  interpolation_valid_time_=Uvar(300,"s");//max interpolation interval
  if(cfg.keywordExists("Max_Temperature_Telemetry_interpolation_time")){
    interpolation_valid_time_=cfg["Max_Temperature_Telemetry_interpolation_time"];
    if(interpolation_valid_time_ < Uvar(0,"s"))
      ErrorMessage("BurstData.cpp::cfg: interpolation valid time should be positive").throwMe();
  }
 
  cerr<<"Max_Temperature_Telemetry_interpolation_time for computing feed/scwg/hga temperatures is "<< interpolation_valid_time_.getInUnits("s")<<" s"<<endl;
  

  burst_data_configured_=true;
}

//--------------
// Constructors
//--------------

BurstData::BurstData(const std::string& filename, const std::string& mode,const std::string& filetype)
  : mode_(mode), filename_(filename),  base_filename_(filename),
    file_(filename,mode), 
    header_handled_(false), record_count_(0), records_counted_(false),
    record_size_(0), header_size_(0), file_number_(0),
    data_read_(false),T_feed_valid_(false),T_hga_valid_(false),
    T_scwg_valid_(false),
    fsc_("CASSINI_SC_COORD","Cassini")
   
  {
    if(!burst_data_configured_){
      BurstDataError e("BurstData must be configured before constructing BurstData object");
      e.throwMe();
    }
    if(filetype=="active"){ 
      ft_=active;
    }
    else if (filetype=="passive"){
      ft_=passive;
    }
    else
    { 
     ErrorMessage e("no proper filetype: either active or passive"); 
     e.throwMe(); 
    }
    

    //---------------------
    //initiallize target name 
    //--------------------
    for(int c=0;c<TARGET_NAME_STRING_SIZE;++c) target_name[c]=' ';
    for(int c=0;c<TARGET_FRAME_NAME_STRING_SIZE;++c) tbf_frame_name[c]=' ';
   

    //--------------
    //default value for number of bursts in flight=1
    //---------------
    num_bursts_in_flight=1;
    
    // Assign pointers into name map (this has to be done manually)
    // for Uvar variables
    // format header = 4 + 6 + 9 + 4 varialbes
    // slow field = 9 + 2 + 2(+2 for db) + 1 + 1(+1 for db)+2+2+1
    // fast field = 6 + 3
    // footer = 8 + 6
    // sub comm data = 6 + 6 + 6 + 6 + 6 + 6 + 1 (subr = 0 .. 6)
    // sub comm data = 4 + 6 + 6 + 6 + 6 + 3 + 1 + 1 (subr = 7..14)
    
 
    parameters_.appendParameter("full_record",0,Parameter::FULL,NULL);

   
    
    parameters_.appendParameter("sync",sync);
    parameters_.appendParameter("sclk",sclk);
    parameters_.appendParameter("record_id",record_id);
    parameters_.appendParameter("scpr",scpr,"Hz");
    parameters_.appendParameter("brst",brst,"s");
    
    parameters_.appendParameter("header_tfi",header_tfi,"s");
    parameters_.appendParameter("header_tnc",header_tnc);
    parameters_.appendParameter("header_typ",header_typ);
    parameters_.appendParameter("header_tca",header_tca);
    parameters_.appendParameter("header_tcb",header_tcb);
    parameters_.appendParameter("header_tcc",header_tcc);
    
    parameters_.appendParameter("pwri",pwri);
    parameters_.appendParameter("vicc",vicc);
    parameters_.appendParameter("vimc",vimc);
    parameters_.appendParameter("tail_len",tail_len);
    parameters_.appendParameter("tail_id",tail_id);
    parameters_.appendParameter("sab_counter",sab_counter);
    parameters_.appendParameter("sab_len",sab_len);
    parameters_.appendParameter("fswm",fswm);
    parameters_.appendParameter("fswc",fswc);
    parameters_.appendParameter("ctbc",ctbc);
    parameters_.appendParameter("ctrx",ctrx);
    
    parameters_.appendParameter("ctps",ctps);
    parameters_.appendParameter("ctbe",ctbe);
    parameters_.appendParameter("ctps_ctbe",ctps_ctbe);
    parameters_.appendParameter("header_end",header_end);
    
    
    //slow field variables
    parameters_.appendParameter("slow_tfi",slow_tfi,"s")  ;
    parameters_.appendParameter("dtn",dtn) ;
    parameters_.appendParameter("slow_typ",slow_typ);
    parameters_.appendParameter("csr",csr);
    parameters_.appendParameter("r_mode",r_mode);
    parameters_.appendParameter("sin",slow_instruction_number);
    parameters_.appendParameter("bem",bem);
    parameters_.appendParameter("baq_mode",baq_mode);
    parameters_.appendParameter("tro",tro);
    
    parameters_.appendParameter("rc_bw",rc_bw,"Hz")  ;
    parameters_.appendParameter("adc",adc,"Hz");
    
    
    parameters_.appendParameter("at1_tot",at1);
    parameters_.appendParameter("at3_tot",at3);
    parameters_.appendParameter("at4_tot",at4);
    
    parameters_.appendParameter("at1_each",at1_each);
    parameters_.appendParameter("at3_each",at3_each);
    parameters_.appendParameter("at4_each",at4_each);
    
    parameters_.appendParameter("rip",rip,"s");
    
    parameters_.appendParameter("csd",csd,"s");
    parameters_.appendParameter("rad",rad);
    
    parameters_.appendParameter("csq",csq);
    parameters_.appendParameter("chirp_length",chirp_length,"s");
    
    parameters_.appendParameter("slow_cfs",slow_cfs,"Hz");
    
    
    //fast field variables
    parameters_.appendParameter("fast_tfi",fast_tfi,"s");
    parameters_.appendParameter("fin",fin);
    parameters_.appendParameter("fast_typ",fast_typ);
    parameters_.appendParameter("pul",pul);
    parameters_.appendParameter("bii",bii);
    parameters_.appendParameter("bpd",bpd,"s");
    
    parameters_.appendParameter("pri",pri,"s");
    parameters_.appendParameter("rwd",rwd,"s");
    parameters_.appendParameter("fast_csf",fast_csf,"Hz");
    
    
    
    
    //footer1 variables  
    parameters_.appendParameter("iebtth",iebtth)  ; 
    parameters_.appendParameter("iebttl",iebttl)  ; 
    parameters_.appendParameter("bgcalls",bgcalls)  ; 
    parameters_.appendParameter("delvmn",delvmn)  ; 
    parameters_.appendParameter("delvda",delvda)  ; 
    parameters_.appendParameter("delvyr",delvyr)  ; 


    parameters_.appendParameter("cnt_rl",cnt_rl)  ; 
    parameters_.appendParameter("cnt_radio",cnt_radio)  ; 
    parameters_.appendParameter("cnt_nd",cnt_nd)  ; 
    parameters_.appendParameter("eout",eout)  ; 
    parameters_.appendParameter("subr",subr)  ;
    parameters_.appendParameter("space_craft_time",space_craft_time)  ;
    parameters_.appendParameter("hip",hip,"s")  ;
    parameters_.appendParameter("cip",cip,"s");
    
  
    
    
    
    
    //sub_comm_data variables 
    //subr = 0
    parameters_.appendParameter("fwdtmp",fwdtmp,"K")  ;
    parameters_.appendParameter("be1tmp",be1tmp,"K")  ;
    parameters_.appendParameter("be2tmp",be2tmp,"K")  ;
    parameters_.appendParameter("be3tmp",be3tmp,"K")  ;
    parameters_.appendParameter("be4tmp",be4tmp,"K")  ;
    parameters_.appendParameter("be5tmp",be5tmp,"K")  ;
    
    //subr = 1
    parameters_.appendParameter("diptmp",diptmp,"K") ;
    parameters_.appendParameter("rlotmp",rlotmp,"K");
    parameters_.appendParameter("tadcal1",tadcal1,"K") ;
    parameters_.appendParameter("nsdtmp",nsdtmp,"K") ;
    parameters_.appendParameter("lnatmp",lnatmp,"K") ;
    parameters_.appendParameter("evdtmp",evdtmp,"K") ;
    
    //subr = 2
    parameters_.appendParameter("mratmp",mratmp,"K")  ;
    parameters_.appendParameter("mruttm",mruttm,"K")  ;
    parameters_.appendParameter("dcgttm",dcgttm,"K")  ;
    parameters_.appendParameter("cucttm",cucttm,"K")  ;
    parameters_.appendParameter("twttmp",twttmp,"K")  ;
    parameters_.appendParameter("epctmp",epctmp,"K")  ;
    
    //subr = 3
    parameters_.appendParameter("tw1ttm",tw1ttm,"K")  ;
    parameters_.appendParameter("ep1ttm",ep1ttm,"K")  ;
    parameters_.appendParameter("p_stmp",p_stmp,"K")  ;
    parameters_.appendParameter("p_sttm",p_sttm,"K")  ;
    parameters_.appendParameter("fguttm",fguttm,"K")  ;
    parameters_.appendParameter("tadcal4",tadcal4,"K")  ;
    
    //subr = 4
    parameters_.appendParameter("esstmp",esstmp,"K")   ;
    parameters_.appendParameter("wgb1t1",wgb1t1,"K")   ;
    parameters_.appendParameter("wgb3t1",wgb3t1,"K")   ;
    parameters_.appendParameter("wgb3t2",wgb3t2,"K")   ;
    parameters_.appendParameter("wgb3t3",wgb3t3,"K")   ;
    parameters_.appendParameter("wgb5t1",wgb5t1,"K")   ;
    
    //subr = 5
    parameters_.appendParameter("pcutmp",pcutmp,"K")  ;
    parameters_.appendParameter("adctmp",adctmp,"K") ;
    parameters_.appendParameter("tadcal2",tadcal2,"K")  ;
    parameters_.appendParameter("ecltmp",ecltmp,"K")  ;
    parameters_.appendParameter("cputmp",cputmp,"K")  ;
    parameters_.appendParameter("memtmp",memtmp,"K")  ;
    
    //subr = 6
    parameters_.appendParameter("sadctmp",sadctmp,"K")  ;
    
    //subr = 7
    parameters_.appendParameter("tadcal3",tadcal3,"K")  ;
    parameters_.appendParameter("frwdpw",frwdpw,"W")  ;
    parameters_.appendParameter("dcgmon",dcgmon)  ;
    parameters_.appendParameter("lpltlm_db",lpltlm_db)  ;
    
    //subr = 8
    parameters_.appendParameter("nsdcur",nsdcur,"A")  ;
    parameters_.appendParameter("hpapsm",hpapsm,"A")  ;
    parameters_.appendParameter("catcur",catcur,"A")  ;
    parameters_.appendParameter("p_smon",p_smon,"A")  ;
    parameters_.appendParameter("svlsta",svlsta,"V")  ;
    parameters_.appendParameter("usotmp",usotmp,"K")  ;
    
    //subr = 9
    parameters_.appendParameter("cpbnkv",cpbnkv,"V");
    parameters_.appendParameter("essvlt",essvlt,"V");
    parameters_.appendParameter("tadcal5",tadcal5,"V");
    parameters_.appendParameter("pcu5v_72",pcu5v_72,"V");
    parameters_.appendParameter("pcu5i_73",pcu5i_73,"A");
    parameters_.appendParameter("pcu5v_74",pcu5v_74,"V");
    
    //subr = 10
    parameters_.appendParameter("pcuti_75",pcuti_75,"A")  ;//typo?
    parameters_.appendParameter("pcu15v_76",pcu15v_76,"V")  ;
    parameters_.appendParameter("pcu15i_77",pcu15i_77,"A")  ;
    parameters_.appendParameter("pcu15v_78",pcu15v_78,"V")  ;
    parameters_.appendParameter("pcu15i_79",pcu15i_79,"A") ;
    parameters_.appendParameter("pcu12v_80",pcu12v_80,"V")  ;
    
    //subr = 11
    parameters_.appendParameter("pcu12i_81",pcu12i_81,"A")  ;
    parameters_.appendParameter("pcucur",pcucur,"A")   ;
    parameters_.appendParameter("pllmon",pllmon,"Hz")   ;
    parameters_.appendParameter("ctu5i",ctu5i,"A")   ;
    parameters_.appendParameter("tadcal6",tadcal6)   ;
    parameters_.appendParameter("pcu9v_86",pcu9v_86,"V")   ;
    
    //subr = 12
    parameters_.appendParameter("pcu9i_87",pcu9i_87,"A")  ;
    parameters_.appendParameter("pcu9v_88",pcu9v_88,"V")   ;
    parameters_.appendParameter("pcu9i_89",pcu9i_89,"A")   ;
    
    //subr = 13
    parameters_.appendParameter("tadcl7",tadcl7,"V");
    
    //subr = 14
    parameters_.appendParameter("shpttm",shpttm,"K");
    
    //number of bursts in flight
    parameters_.appendParameter("num_bursts_in_flight",num_bursts_in_flight);
   
    parameters_.appendParameter("raw_active_mode_length",Nradar_data);
    parameters_.appendParameter("raw_active_mode_rms",sizeof(float),Parameter::FLOAT,(void*)(&rms_radar_data));

    //debug
    //cout<<"after raw active mode rms "<< parameters_.offset<<endl;
   


    // computed geometry params
    parameters_.appendParameter("engineer_qual_flag",quality_flag);

    //debug
    //cout<<"before time field "<< parameters_.offset<<endl;

    parameters_.appendParameter("t",TIME_LENGTH_IN_FILE,Parameter::TIME,
				(void*)(&t),NULL,"s");

    //debug 
    //cout<<"after time field "<< parameters_.offset<<endl;

    parameters_.appendParameter("transmit_time_offset",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&transmit_time_offset),NULL,"s");
    parameters_.appendParameter("time_from_closest_approach",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&time_from_closest_approach),NULL,"s");
    parameters_.appendParameter("time_from_epoch",sizeof(double),Parameter::UVAR_DOUBLE,
				(void*)(&time_from_epoch),NULL,"s");
    
    //debug
    //cout<<"index before target name "<< parameters_.offset<<endl;

    //target name and target frame
    parameters_.appendParameter("target_name", TARGET_NAME_STRING_SIZE, Parameter::STRING,
				(void*)(&target_name),NULL,"");
    
    //debug
    //cout<<"index after target name "<< parameters_.offset<<endl;
    parameters_.appendParameter("tbf_frame_name",TARGET_FRAME_NAME_STRING_SIZE,Parameter::STRING,(void*)(&tbf_frame_name),NULL,"");

    //debug
    //    cout<<"index after tbf frame name "<< parameters_.offset<<endl;


    //target specific angles
    parameters_.appendParameter("pole_right_ascension",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&pole_right_ascension),NULL,"deg");
    parameters_.appendParameter("pole_declination",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&pole_declination),NULL,"deg");
    parameters_.appendParameter("target_rotation_rate",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&target_rotation_rate),NULL,"deg/s");
    parameters_.appendParameter("target_rotation_angle",sizeof(double),
				Parameter::UVAR_DOUBLE,(void*)(&target_rotation_angle),NULL,"deg");


    parameters_.appendParameter("scwg_tmp",T_scwg,"K");
    parameters_.appendParameter("feed_tmp",T_feed,"K");
    parameters_.appendParameter("hga_tmp",T_hga,"K");

    //debug
    //cout<<"before beam number "<< parameters_.offset<<endl;

    parameters_.appendParameter("beam_number",beam_number);
    //state vector in J2000


   
    //state vector in target frame
    parameters_.appendParameter("sc_state_J2000",6*sizeof(double),
				Parameter::STATEVEC,(void*)(&sc_state_J2000),
				&j2000_);
    //state vector in target frame
    parameters_.appendParameter("sc_state_target",6*sizeof(double),
				Parameter::STATEVEC,(void*)(&sc_state_target),&ftarget_);
    //direction in J2000
    parameters_.appendParameter("sc_X_J2000",3*sizeof(double),Parameter::DIRVEC,
				(void*)(&sc_X_J2000),&j2000_);
    parameters_.appendParameter("sc_Y_J2000",3*sizeof(double),Parameter::DIRVEC,
				(void*)(&sc_Y_J2000),&j2000_);
    parameters_.appendParameter("sc_Z_J2000",3*sizeof(double),Parameter::DIRVEC,
				(void*)(&sc_Z_J2000),&j2000_);
    //direction in target frame
    parameters_.appendParameter("sc_X_target",3*sizeof(double),
				Parameter::DIRVEC,
				(void*)(&sc_X_target),&ftarget_);
    parameters_.appendParameter("sc_Y_target",3*sizeof(double),
				Parameter::DIRVEC,
				(void*)(&sc_Y_target),&ftarget_);
    parameters_.appendParameter("sc_Z_target",3*sizeof(double),
				Parameter::DIRVEC,
				(void*)(&sc_Z_target),&ftarget_);
    
    //rotation in j2000
    parameters_.appendParameter("rot_vel_J2000",3*sizeof(double),
				Parameter::ROTVEL,(void*)(&rot_vel_J2000),
				&j2000_,"deg/s");
    //rotation in target frame
    parameters_.appendParameter("rot_vel_target",3*sizeof(double),
				Parameter::ROTVEL,(void*)(&rot_vel_target),
				&ftarget_,"deg/s");  
    parameters_.appendParameter("norm_cnt_rl",norm_cnt_rl);
    parameters_.appendParameter("norm_cnt_nd",norm_cnt_nd);
    parameters_.appendParameter("norm_cnt_radio",norm_cnt_radio);
    

    //debug
    //cout<<"before sbdr write "<< parameters_.offset<<endl;

    //SBDR Science Data Field
    parameters_.appendParameter("science_qual_flag",science_qual_flag);
    parameters_.appendParameter("system_gain",system_gain);
    parameters_.appendParameter("antenna_brightness_temp",antenna_brightness_temp,"K");
    parameters_.appendParameter("system_noise_temp",system_noise_temp,"K");
    parameters_.appendParameter("abt_std",abt_std);
    
    
    parameters_.appendParameter("pass_geom_time_offset",pass_geom_time_offset);

    parameters_.appendParameter("pass_pol_angle",pass_pol_angle,"deg");
    parameters_.appendParameter("pass_emission_angle",pass_emission_angle,"deg");
    parameters_.appendParameter("pass_azimuth_angle",pass_azimuth_angle,"deg");
   

    parameters_.appendParameter("pass_centroid_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*) (&pass_centroid_lon),NULL,"deg");
    parameters_.appendParameter("pass_centroid_lat",pass_centroid_lat,"deg");

    parameters_.appendParameter("pass_major_width",pass_major_width);
    parameters_.appendParameter("pass_minor_width",pass_minor_width);
    
    parameters_.appendParameter("pass_ellipse_pt1_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*) (&pass_ellipse_pt1_lon),NULL,"deg");
    parameters_.appendParameter("pass_ellipse_pt2_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*) (&pass_ellipse_pt2_lon),NULL,"deg");
    parameters_.appendParameter("pass_ellipse_pt3_lon",sizeof(float), Parameter::UVAR_WEST_LON,
				(void*) (&pass_ellipse_pt3_lon),NULL,"deg");
    parameters_.appendParameter("pass_ellipse_pt4_lon",sizeof(float), Parameter::UVAR_WEST_LON,
				(void*)(&pass_ellipse_pt4_lon),NULL,"deg");

    parameters_.appendParameter("pass_ellipse_pt1_lat",pass_ellipse_pt1_lat,"deg");   
    parameters_.appendParameter("pass_ellipse_pt2_lat",pass_ellipse_pt2_lat,"deg");
    parameters_.appendParameter("pass_ellipse_pt3_lat",pass_ellipse_pt3_lat,"deg");
    parameters_.appendParameter("pass_ellipse_pt4_lat",pass_ellipse_pt4_lat,"deg");
    parameters_.appendParameter("num_pulses_received",num_pulses_received);

    parameters_.appendParameter("total_echo_energy",total_echo_energy,"J");  
    parameters_.appendParameter("noise_echo_energy",noise_echo_energy,"J");
    parameters_.appendParameter("x_factor",x_factor);
    parameters_.appendParameter("sigma0_uncorrected",sigma0_uncorrected);
    parameters_.appendParameter("sigma0_corrected",sigma0_corrected);
    parameters_.appendParameter("sigma0_uncorrected_std",sigma0_uncorrected_std);
    
    //debug
    //cout<<"before altitude means "<< parameters_.offset<<endl;

    parameters_.appendParameter("altitude_means",altitude_means);
    parameters_.appendParameter("altitude_means_std",altitude_means_std);


    parameters_.appendParameter("act_geom_time_offset",act_geom_time_offset);

    parameters_.appendParameter("act_pol_angle",act_pol_angle,"deg");
    parameters_.appendParameter("act_incidence_angle",act_incidence_angle,"deg");
    parameters_.appendParameter("act_azimuth_angle",act_azimuth_angle,"deg");
    
    parameters_.appendParameter("act_centroid_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*) (&act_centroid_lon),NULL,"deg");
    parameters_.appendParameter("act_centroid_lat",act_centroid_lat,"deg");

    parameters_.appendParameter("act_major_width",act_major_width,"deg");
    parameters_.appendParameter("act_minor_width",act_minor_width,"deg");

    parameters_.appendParameter("act_ellipse_pt1_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*)(&act_ellipse_pt1_lon),NULL,"deg");
    parameters_.appendParameter("act_ellipse_pt2_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*)(&act_ellipse_pt2_lon),NULL,"deg");
    parameters_.appendParameter("act_ellipse_pt3_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*)(&act_ellipse_pt3_lon),NULL,"deg");   
    parameters_.appendParameter("act_ellipse_pt4_lon",sizeof(float),Parameter::UVAR_WEST_LON,
				(void*)(&act_ellipse_pt4_lon),NULL,"deg");


    parameters_.appendParameter("act_ellipse_pt1_lat",act_ellipse_pt1_lat,"deg");
    parameters_.appendParameter("act_ellipse_pt2_lat",act_ellipse_pt2_lat,"deg");
    parameters_.appendParameter("act_ellipse_pt3_lat",act_ellipse_pt3_lat,"deg");
    parameters_.appendParameter("act_ellipse_pt4_lat",act_ellipse_pt4_lat,"deg");
    
    parameters_.appendParameter("altimeter_profile_range_start",altimeter_profile_range_start);
    parameters_.appendParameter("altimeter_profile_range_step",altimeter_profile_range_step);
    parameters_.appendParameter("altimeter_profile_length",altimeter_profile_length);
    parameters_.appendParameter("sar_azimuth_res",sar_azimuth_res);
    parameters_.appendParameter("sar_range_res",sar_range_res);
    parameters_.appendParameter("sar_centroid_bidr_lon",sar_centroid_bidr_lon,"deg");
    parameters_.appendParameter("sar_centroid_bidr_lat",sar_centroid_bidr_lat,"deg");

    //degub
    //cout<<"total bytes "<< parameters_.offset<<endl;    
  }

//--------------
// I/O
//--------------

void
BurstData::writePassiveSABData(int abs_record_no)
{
  //updated on Oct-29-01
  //Write out all the available public variable including science data
  // header slow fast

  if(abs_record_no>=0)
    record_id=data_take_id_*1000000+abs_record_no;

  file_.write(sync); 
  file_.write(sclk);
  file_.write(record_id);

  scpr.writeFloat(file_,0,"Hz");
  brst.writeFloat(file_,0,"s");

  header_tfi.writeFloat(file_,0,"s");
  file_.write(header_tnc);
  file_.write(header_typ);
  file_.write(header_tca);
  file_.write(header_tcb);
  file_.write(header_tcc);

  file_.write(pwri);
  file_.write(vicc);
  file_.write(vimc);
  file_.write(tail_len);
  file_.write(tail_id);
  file_.write(sab_counter);
  file_.write(sab_len);
  file_.write(fswm);
  file_.write(fswc);
  file_.write(ctbc);
  file_.write(ctrx);
 
  file_.write(ctps);
  file_.write(ctbe);
  file_.write(ctps_ctbe);
  file_.write(header_end);
  
  //slow field variables
  slow_tfi.writeFloat(file_,0,"s");
  file_.write(dtn);
  file_.write(slow_typ);
  file_.write(csr);
  file_.write(r_mode);
  file_.write(slow_instruction_number);
  file_.write(bem);
  file_.write(baq_mode);
  tro.writeFloat(file_,0,"s");
  
 
  rc_bw.writeFloat(file_,0,"Hz");
  adc.writeFloat(file_,0,"Hz");

  at1.writeFloat(file_,0);
  at3.writeFloat(file_,0);
  at4.writeFloat(file_,0);

  file_.write(at1_each);
  file_.write(at3_each);
  file_.write(at4_each);

  rip.writeFloat(file_,0,"s");
  csd.writeFloat(file_,0,"s");
  file_.write(rad);

  file_.write(csq);
  chirp_length.writeFloat(file_,0,"s");
  slow_cfs.writeFloat(file_,0,"Hz");


  //fast field variables
  fast_tfi.writeFloat(file_,0,"s");
  file_.write(fin);
  file_.write(fast_typ);
  file_.write(pul);
  file_.write(bii);
  bpd.writeFloat(file_,0,"s");

  pri.writeFloat(file_,0,"s");
  rwd.writeFloat(file_,0,"s");
  fast_csf.writeFloat(file_,0,"Hz");

 //footer variables

  file_.write(iebtth);
  file_.write(iebttl);
  file_.write(bgcalls);
  file_.write(delvmn);
  file_.write(delvda);
  file_.write(delvyr);


  file_.write(cnt_rl);
  file_.write(cnt_radio);
  file_.write(cnt_nd);
  file_.write(eout);
  file_.write(subr);
  file_.write(space_craft_time);
  hip.writeFloat(file_,0,"s");
  cip.writeFloat(file_,0,"s");
    


 
 //sub_comm_data variables 
  
 fwdtmp.writeFloat(file_,0,"K");
 be1tmp.writeFloat(file_,0,"K");
 be2tmp.writeFloat(file_,0,"K");
 be3tmp.writeFloat(file_,0,"K");
 be4tmp.writeFloat(file_,0,"K");
 be5tmp.writeFloat(file_,0,"K");


 diptmp.writeFloat(file_,0,"K");
 rlotmp.writeFloat(file_,0,"K");
 tadcal1.writeFloat(file_,0,"K");
 nsdtmp.writeFloat(file_,0,"K");
 lnatmp.writeFloat(file_,0,"K");
 evdtmp.writeFloat(file_,0,"K");
 
 mratmp.writeFloat(file_,0,"K");
 mruttm.writeFloat(file_,0,"K");
 dcgttm.writeFloat(file_,0,"K");
 cucttm.writeFloat(file_,0,"K");
 twttmp.writeFloat(file_,0,"K");
 epctmp.writeFloat(file_,0,"K"); 

 tw1ttm.writeFloat(file_,0,"K");
 ep1ttm.writeFloat(file_,0,"K");
 p_stmp.writeFloat(file_,0,"K");
 p_sttm.writeFloat(file_,0,"K");
 fguttm.writeFloat(file_,0,"K");
 tadcal4.writeFloat(file_,0,"K"); 

 esstmp.writeFloat(file_,0,"K");
 wgb1t1.writeFloat(file_,0,"K");
 wgb3t1.writeFloat(file_,0,"K");
 wgb3t2.writeFloat(file_,0,"K");
 wgb3t3.writeFloat(file_,0,"K");
 wgb5t1.writeFloat(file_,0,"K");

 pcutmp.writeFloat(file_,0,"K");
 adctmp.writeFloat(file_,0,"K");
 tadcal2.writeFloat(file_,0,"K");
 ecltmp.writeFloat(file_,0,"K");
 cputmp.writeFloat(file_,0,"K");
 memtmp.writeFloat(file_,0,"K");

 sadctmp.writeFloat(file_,0,"K");

 tadcal3.writeFloat(file_,0,"K");
 file_.writeFloatInUnits(frwdpw,"W");
 dcgmon.writeFloat(file_,0);
 lpltlm_db.writeFloat(file_,0);

 nsdcur.writeFloat(file_,0,"A");
 hpapsm.writeFloat(file_,0,"A");
 catcur.writeFloat(file_,0,"A");
 p_smon.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(svlsta,"V");
 usotmp.writeFloat(file_,0,"K");

 file_.writeFloatInUnits(cpbnkv,"V");
 file_.writeFloatInUnits(essvlt,"V");
 file_.writeFloatInUnits(tadcal5,"V");
 file_.writeFloatInUnits(pcu5v_72,"V");
 pcu5i_73.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(pcu5v_74,"V");

 pcuti_75.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(pcu15v_76,"V");
 pcu15i_77.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(pcu15v_78,"V");
 pcu15i_79.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(pcu12v_80,"V");

 pcu12i_81.writeFloat(file_,0,"A");
 pcucur.writeFloat(file_,0,"A");
 pllmon.writeFloat(file_,0,"Hz");
 ctu5i.writeFloat(file_,0,"A");
 tadcal6.writeFloat(file_,0);
 file_.writeFloatInUnits(pcu9v_86,"V");

 pcu9i_87.writeFloat(file_,0,"A");
 file_.writeFloatInUnits(pcu9v_88,"V");
 pcu9i_89.writeFloat(file_,0,"A");

 file_.writeFloatInUnits(tadcl7,"V");

 shpttm.writeFloat(file_,0,"K");


  //num in bursts in flight and data record length
  file_.write(num_bursts_in_flight);
  file_.write(Nradar_data); 
  file_.write(rms_radar_data);
}


void
BurstData::readPassiveSABData(FileMgr* f)
{
  bool data_take_id_checked=false;
  if(f==NULL) f= &file_;
  f->read(sync);
  if(sync!=SYNC_VAL){
    BurstDataError e("Bad sync");
    e.throwMe();
  }
  f->read(sclk);
  f->read(record_id);
  if(!data_take_id_checked){
    if(record_id/1000000 != data_take_id_){
      cerr <<"Data Take ID mismatch between ID:"
	+toStr(data_take_id_)+" and file:"+f->name() << endl;
      cerr <<"Using value in burst data file "+toStr(record_id/1000000) << endl;
      data_take_id_=record_id/1000000;
    }
    data_take_id_checked=true;
  }
  scpr.readFloat(*f,0,"Hz");
  brst.readFloat(*f,0,"s");

  header_tfi.readFloat(*f,0,"s");
  f->read(header_tnc);
  f->read(header_typ);
  f->read(header_tca);
  f->read(header_tcb);
  f->read(header_tcc);

  f->read(pwri);
  f->read(vicc);
  f->read(vimc);
  f->read(tail_len);
  f->read(tail_id);
  f->read(sab_counter);
  f->read(sab_len);
  f->read(fswm);
  f->read(fswc);
  f->read(ctbc);
  f->read(ctrx);

  f->read(ctps);
  f->read(ctbe);
  f->read(ctps_ctbe);
  f->read(header_end);
  
  //slow field variables
  slow_tfi.readFloat(*f,0,"s");
  f->read(dtn);
  f->read(slow_typ);
  f->read(csr);
  f->read(r_mode);
  f->read(slow_instruction_number);
  f->read(bem);
  f->read(baq_mode);
  tro.readFloat(*f,0,"s");

  rc_bw.readFloat(*f,0,"Hz");
  adc.readFloat(*f,0,"Hz");

  at1.readFloat(*f,0);
  at3.readFloat(*f,0);
  at4.readFloat(*f,0);
  f->read(at1_each);
  f->read(at3_each);
  f->read(at4_each);
 
  rip.readFloat(*f,0,"s");
  csd.readFloat(*f,0,"s");
  f->read(rad);

  f->read(csq);
  chirp_length.readFloat(*f,0,"s");
  slow_cfs.readFloat(*f,0,"Hz");


  //fast field variables
  fast_tfi.readFloat(*f,0,"s");
  f->read(fin);
  f->read(fast_typ);
  f->read(pul);
  f->read(bii);
  bpd.readFloat(*f,0,"s");

  pri.readFloat(*f,0,"s");
  rwd.readFloat(*f,0,"s");
  fast_csf.readFloat(*f,0,"Hz");

  //footer variables
  f->read(iebtth);
  f->read(iebttl);
  f->read(bgcalls);
  f->read(delvmn);
  f->read(delvda);
  f->read(delvyr);



  f->read(cnt_rl);
  f->read(cnt_radio);
  f->read(cnt_nd);
  f->read(eout);
  f->read(subr);
  f->read(space_craft_time);
  hip.readFloat(*f,0,"s");
  cip.readFloat(*f,0,"s");
 

 
    
 //sub_comm_data variables 
  fwdtmp.readFloat(*f,0,"K");
  be1tmp.readFloat(*f,0,"K");
  be2tmp.readFloat(*f,0,"K");
  be3tmp.readFloat(*f,0,"K");
  be4tmp.readFloat(*f,0,"K");
  be5tmp.readFloat(*f,0,"K");

  diptmp.readFloat(*f,0,"K");
  rlotmp.readFloat(*f,0,"K");
  tadcal1.readFloat(*f,0,"K");
  nsdtmp.readFloat(*f,0,"K");
  lnatmp.readFloat(*f,0,"K");
  evdtmp.readFloat(*f,0,"K");

  mratmp.readFloat(*f,0,"K");
  mruttm.readFloat(*f,0,"K");
  dcgttm.readFloat(*f,0,"K");
  cucttm.readFloat(*f,0,"K");
  twttmp.readFloat(*f,0,"K");
  epctmp.readFloat(*f,0,"K"); 

  tw1ttm.readFloat(*f,0,"K");
  ep1ttm.readFloat(*f,0,"K");
  p_stmp.readFloat(*f,0,"K");
  p_sttm.readFloat(*f,0,"K");
  fguttm.readFloat(*f,0,"K");
  tadcal4.readFloat(*f,0,"K"); 

  esstmp.readFloat(*f,0,"K");
  wgb1t1.readFloat(*f,0,"K");
  wgb3t1.readFloat(*f,0,"K");
  wgb3t2.readFloat(*f,0,"K");
  wgb3t3.readFloat(*f,0,"K");
  wgb5t1.readFloat(*f,0,"K");

  pcutmp.readFloat(*f,0,"K");
  adctmp.readFloat(*f,0,"K");
  tadcal2.readFloat(*f,0,"K");
  ecltmp.readFloat(*f,0,"K");
  cputmp.readFloat(*f,0,"K");
  memtmp.readFloat(*f,0,"K");

  sadctmp.readFloat(*f,0,"K");

  tadcal3.readFloat(*f,0,"K");
  f->readFloatInUnits(frwdpw,WtoMW,"MW");
  dcgmon.readFloat(*f,0);
  lpltlm_db.readFloat(*f,0);

  nsdcur.readFloat(*f,0,"A");
  hpapsm.readFloat(*f,0,"A");
  catcur.readFloat(*f,0,"A");
  p_smon.readFloat(*f,0,"A");
  f->readFloatInUnits(svlsta,VtoMV,"MV");
  usotmp.readFloat(*f,0,"K");

  f->readFloatInUnits(cpbnkv,VtoMV,"MV");
  f->readFloatInUnits(essvlt,VtoMV,"MV");
  f->readFloatInUnits(tadcal5,VtoMV,"MV");
  f->readFloatInUnits(pcu5v_72,VtoMV,"MV");
  pcu5i_73.readFloat(*f,0,"A");
  f->readFloatInUnits(pcu5v_74,VtoMV,"MV");

  pcuti_75.readFloat(*f,0,"A");
  f->readFloatInUnits(pcu15v_76,VtoMV,"MV");
  pcu15i_77.readFloat(*f,0,"A");
  f->readFloatInUnits(pcu15v_78,VtoMV,"MV");
  pcu15i_79.readFloat(*f,0,"A");
  f->readFloatInUnits(pcu12v_80,VtoMV,"MV");

  pcu12i_81.readFloat(*f,0,"A");
  pcucur.readFloat(*f,0,"A");
  pllmon.readFloat(*f,0,"Hz");
  ctu5i.readFloat(*f,0,"A");
  tadcal6.readFloat(*f,0);
  f->readFloatInUnits(pcu9v_86,VtoMV,"MV");

  pcu9i_87.readFloat(*f,0,"A");
  f->readFloatInUnits(pcu9v_88,VtoMV,"MV");
  pcu9i_89.readFloat(*f,0,"A");

  f->readFloatInUnits(tadcl7,VtoMV,"MV");

  shpttm.readFloat(*f,0,"K");
  f->read(num_bursts_in_flight);
  f->read(Nradar_data); 
  f->read(rms_radar_data);

}

void
BurstData::checkEndedness()
{
  unsigned int fp=file_.getPosition();
  if(fp!=0){
    BurstDataError e("CheckEndedness must be called prior to reading file");
    e.throwMe();
  }
  Parameter* p=parameters_.getParam("sync");
  file_.setPosition(header_size_ + p->offset);
  file_.readKnown(SYNC_VAL);
  file_.rewind();
}
//------------------------------------------------------------------
// readSclkOnly()
//
// Read sclk from a BurstData file and skip the rest
//------------------------------------------------------------------

void BurstData::readSclkOnly()
{
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't read from output file " + filename_,
	       BurstData::read_error);
      e.throwMe();
    }
  if (!header_handled_)
    {
      BurstDataError e("Can't read BurstData record until header is read",
	       BurstData::read_error);
      e.throwMe();
    }

  // This check has the side effect of moving to the next file if necessary
  if(eof())
    {
      BurstDataError e("Unexpected EOF in BurstData file "+filename_,
	       BurstData::read_error);
      e.throwMe();
    } 

  int start_position=file_.getPosition();  
  file_.read(sync); 
  file_.read(sclk);
  file_.setPosition(start_position+record_size_);
  
}

//------------------------------------------------------------------
// Check to see if current record is SAR after it has been read
//
// 
//------------------------------------------------------------------
bool  BurstData::isSAR()  const {
  if(adc>Uvar(900000,"Hz") && adc<Uvar(2100000,"Hz")){
    return(true);
  }
  return(false);
}

//------------------------------------------------------------------
// Check to see if current record is a Cal after it has been read
//
// 
//------------------------------------------------------------------
bool  BurstData::isCal() const {
  if(csr!=0 && csr!=8 && csr!=6){
    return(true);
  }
  return(false);
}


//------------------------------------------------------------------
// Check to see if current record is passive
// (radiometer only) after it has been read
//
// 
//------------------------------------------------------------------
bool  BurstData::isPassive(){
  if(csr==6){
    return(true);
  }
  return(false);
}

//------------------------------------------------------------------
// skipRecord(int n=1)
//
// Skips n records, default is n=1
// Throws an exception on eof
//------------------------------------------------------------------

void BurstData::skipRecord(int n){
  for(int c=0;c<n;c++){

    // side effect of skipping to next file if necessary
    if(eof()){
      BurstDataError e("EOF during skipRecord in "+filename_,
		     BurstData::read_error);
      e.throwMe();
    }
    int start_position=file_.getPosition();
    file_.setPosition(start_position+record_size_);
  }
    
}

int BurstData::maxSabCounter(){
  // save start_position and current sab_counter value
  int old_sab_counter=sab_counter;
  int start_position=file_.getPosition();
  int start_file_no=file_number_;

  // goto to beginning of last record of last file
  openLastFile();
  file_.setPosition(header_size_+record_size_*(record_count_-1));

  // read sab counter of last record
  readParameter("sab_counter");
  int retval=sab_counter;
  
  // reset file number, file position and sab_counter
  setFileNumber(start_file_no);
  file_.setPosition(start_position);
  sab_counter=old_sab_counter;

   
  return(retval);
}

//-------------------------------------------
// go to start of First Record of current file
// for multiple files go to first record of first file
//--------------------------------------------
void BurstData::gotoFirstRecord()
{
  if (file_number_!=0){
    //   BurstData::BurstDataError 
    //  e("Cannot go to first record of multiple file data set", 
    //	BurstData::unspecified);
    // e.throwMe();
    setFileNumber(0);
  }
  file_.setPosition(header_size_);
}

//------------------------------------------------------------------
// createNextFile()
//------------------------------------------------------------------

void BurstData::createNextFile(){
  if (mode_ == "r" || mode_ == "rb")
    {
      BurstDataError e("Can't create next input file " + filename_, 
		       BurstData::write_error);
      e.throwMe();
    }

  // Put correct number of records in header
  rewriteHeader();

  // change file number and name
  file_number_++;
  if(file_number_>=MAX_NUM_FILES){
    BurstDataError e("Maximum number of BurstData files " + base_filename_ + 
		     ".x exceeded", read_error);
    e.throwMe();
  }
  char suffix[4];
  sprintf(suffix,".%d",file_number_+1);
  filename_=base_filename_+suffix;

  // open  new file 
  file_.reopen(filename_,mode_);
  record_count_=0;
  records_counted_=false;
  header_handled_=false;

  // write header
  writeHeader();
}

//------------------------------------------------------------------
// openNextFile()
//------------------------------------------------------------------

void BurstData::openNextFile(){
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't open next output file for reading" + filename_, 
		       BurstData::read_error);
      e.throwMe();
    }

  // change file number and name
  file_number_++;
  if(file_number_>=MAX_NUM_FILES){
    BurstDataError e("Maximum number of BurstData files " + base_filename_ + 
		     ".X exceeded", read_error);
    e.throwMe();
  }
  if(file_number_==0) filename_=base_filename_;
  else {
    char suffix[4];
    sprintf(suffix,".%d",file_number_+1);
    filename_=base_filename_+suffix;
  }

  // open  new file 
  file_.reopen(filename_,mode_);
  record_count_=0;
  records_counted_=false;
  header_handled_=false;

  // read header
  readHeader();
}

//--------------
// Predicates
//--------------

//----------------------------------------------------------------
// eof
//
// Check for end of file for the current read file.
// If current file is opened for writing, an exception is thrown.
// if continuation files exist it opens the next one
// and returns false
//----------------------------------------------------------------

bool BurstData::eof()  
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't check eof for output file " + filename_,
	       BurstData::read_error);
      e.throwMe();
    }
  int position=file_.getPosition();
  if(file_.eof() && position > MAX_FILE_LENGTH - (int)record_size_){
    try{
      openNextFile();
    }
    // returns true if openNextFile fails
    catch(FileMgr::IoError e){
      return(true);
    }
  }
  return(file_.eof());
  }

void BurstData::openLastFile(){
  while(1){
    try{
      openNextFile();
    }
    catch(ErrorMessage e){
      break;
    }
  }
  file_number_-=2;
  openNextFile();
}


//-----------------
// Parameter Access
//-----------------

//------------------------------------------------------------------------
// getParamNameList(name_list)
//
// Return a vector of strings containing the names of all the parameters
// in this BurstData object.
//------------------------------------------------------------------------

void BurstData::getParamNameList(list<string>& name_list)
  {
  name_list.clear();
  for (list<Parameter>::const_iterator p = parameters_.param_list.begin();
       p != parameters_.param_list.end(); ++p)
    {
    name_list.push_back(p->name);
    }
  }

int BurstData::returnParams(double* p2d,char* param_list[], int nparams,
			    bool range_flag, double p1min, double p1max,
			    Dvec* p){
  bool special_case=false; // check to see if range checking on record num
                           // is employed
  static bool initialized=false;
  double num;
  double rval=0.0 ;
  if(range_flag && initialized){ // first record has to be read to initialize
    if (strcmp(param_list[0],"record_num")==0){
      num=getSpecialParam("record_num");
      if(num < p1min){ 
	skipRecord();
	return(2);
      }
      if(num > p1max){
        fflush(stdout);
	cerr << "Data extraction completed successfully." << endl;
	exit(0);
      }
      special_case=true;
    }
  }
  bool *special=new bool[nparams];
  enableParams(param_list,special,nparams);
  if(special_case) special[0]=false; // don't get record_num twice
  for (int i=0;i<nparams;i++){
    if (special[i]) enableSpecialParam(param_list[i]);
  }
  readParams();
  disableParams(param_list,nparams);
  for (int i=0;i<nparams;i++){
    if (special[i]) disableSpecialParam(param_list[i]);
  }
  initialized=true;
  for (int i=0;i<nparams;i++){
    if (special[i]){
      rval=getSpecialParam(param_list[i]);
    }
    else if(special_case && i==0) rval=num; // don't get record_num twice
    else rval=getParam(param_list[i]);  // <-- change here
    
    if(i==0 && range_flag && (rval<p1min || rval>p1max)){
      delete[] special;
      return(1);
    }
    //    cout.precision(10);
    //   cout << value << " " ;

    *(p2d+i)=rval;
    if(p) (*p)(i)=rval;

  }
  //  cout << endl;

  delete[] special;
  return(0);
}


void
BurstData::extractParams(char* param_list[], int nparams, bool range_flag, 
			 double p1min, double p1max, Dvec* p){
  bool special_case=false; // check to see if range checking on record num
                           // is employed
  static bool initialized=false;
  double num;
  if(range_flag && initialized){ // first record has to be read to initialize
    if (strcmp(param_list[0],"record_num")==0){
      num=getSpecialParam("record_num");
      if(num < p1min){ 
	skipRecord();
	return;
      }
      if(num > p1max){
        fflush(stdout);
	cerr << "Data extraction completed successfully." << endl;
	exit(0);
      }
      special_case=true;
    }
  }
  bool *special=new bool[nparams];
  enableParams(param_list,special,nparams);
  if(special_case) special[0]=false; // don't get record_num twice
  for (int i=0;i<nparams;i++){
    if (special[i]) enableSpecialParam(param_list[i]);
  }
  readParams();
  disableParams(param_list,nparams);
  for (int i=0;i<nparams;i++){
    if (special[i]) disableSpecialParam(param_list[i]);
  }
  initialized=true;
  for (int i=0;i<nparams;i++){
    double value;
    if (special[i]){
      value=getSpecialParam(param_list[i]);
    }
    else if(special_case && i==0) value=num; // don't get record_num twice
    else value=getParam(param_list[i]);
    
    if(i==0 && range_flag && (value<p1min || value>p1max)){
      delete[] special;
      return;
    }
    cout.precision(10);
    cout << value << " " ;
    if(p) (*p)(i)=value;
  }
  cout << endl;
  delete[] special;
}



void BurstData::enableParams(char* param_list[], bool* special, int nparams){
  for (int i=0;i<nparams;i++){
    special[i]=false;
    if(parameters_.isParam(param_list[i]))
      parameters_.enable(param_list[i]);
    else special[i]=true;
  }
}

// Virtual method
void BurstData::enableSpecialParam(const char* name){
  if(strcmp(name,"rel_time")==0){
    parameters_.enable("sclk");
  }
  else if(strcmp(name,"dist_from_earth")==0){
    parameters_.enable("sc_state_J2000");
  }
  else if(strcmp(name,"dist_to_target")==0){
    parameters_.enable("sc_state_target");
  }
  else if(strcmp(name,"speed_from_earth")==0){
    parameters_.enable("sc_state_J2000");
  }
  else if(strcmp(name,"speed_to_target")==0){
    parameters_.enable("sc_state_target");
  }
  else if(strcmp(name,"angle_from_sun")==0){
    parameters_.enable("sc_Z_J2000");
    parameters_.enable("sc_state_J2000");
    parameters_.enable("t");
    parameters_.enable("engineer_qual_flag");
  }
  else if(strcmp(name,"sc_in_target_frame_position_x")==0||
	  strcmp(name,"sc_in_target_frame_position_y")==0||
	  strcmp(name,"sc_in_target_frame_position_z")==0||
	  strcmp(name,"sc_in_target_frame_velocity_x")==0||
	  strcmp(name,"sc_in_target_frame_velocity_y")==0||
	  strcmp(name,"sc_in_target_frame_velocity_z")==0){
    parameters_.enable("sc_state_target");

  }
  else if(strcmp(name,"sc_X_target_x")==0||
	  strcmp(name,"sc_X_target_y")==0||
	  strcmp(name,"sc_X_target_z")==0){
    parameters_.enable("sc_X_target");
  }
  else if(strcmp(name,"sc_Y_target_x")==0||
	  strcmp(name,"sc_Y_target_y")==0||
	  strcmp(name,"sc_Y_target_z")==0){
    parameters_.enable("sc_Y_target");
  }
  else if(strcmp(name,"sc_Z_target_x")==0||
	  strcmp(name,"sc_Z_target_y")==0||
	  strcmp(name,"sc_Z_target_z")==0){
    parameters_.enable("sc_Z_target");
  }

  else if(strcmp(name,"record_num")!=0){
    BurstDataError e("No special parameter called "+string(name),unknown_parameter);
    e.throwMe();
  }
}

void BurstData::disableParams(char* param_list[], int nparams){
  for (int i=0;i<nparams;i++){
    if(parameters_.isParam(param_list[i]))
      parameters_.disable(param_list[i]);
  }
}

// virtual method
void BurstData::disableSpecialParam(const char* name){
  if(strcmp(name,"rel_time")==0){
     parameters_.disable("sclk");
  }
  else if(strcmp(name,"dist_from_earth")==0){
    parameters_.disable("sc_state_J2000");
  }
  else if(strcmp(name,"dist_to_target")==0){
    parameters_.disable("sc_state_target");
  }
  else if(strcmp(name,"speed_from_earth")==0){
    parameters_.disable("sc_state_J2000");
  }
  else if(strcmp(name,"speed_to_target")==0){
    parameters_.disable("sc_state_target");
  }
  else if(strcmp(name,"angle_from_sun")==0){
    parameters_.disable("sc_Z_J2000");
    parameters_.disable("sc_state_J2000");
    parameters_.disable("t");
    parameters_.disable("engineer_qual_flag");
  }
  else if(strcmp(name,"sc_in_target_frame_position_x")==0||
	  strcmp(name,"sc_in_target_frame_position_y")==0||
	  strcmp(name,"sc_in_target_frame_position_z")==0||
	  strcmp(name,"sc_in_target_frame_velocity_x")==0||
	  strcmp(name,"sc_in_target_frame_velocity_y")==0||
	  strcmp(name,"sc_in_target_frame_velocity_z")==0){
    parameters_.disable("sc_state_target");
  }
  else if(strcmp(name,"sc_X_target_x")==0||
	  strcmp(name,"sc_X_target_y")==0||
	  strcmp(name,"sc_X_target_z")==0){
    parameters_.disable("sc_X_target");
  }
  else if(strcmp(name,"sc_Y_target_x")==0||
	  strcmp(name,"sc_Y_target_y")==0||
	  strcmp(name,"sc_Y_target_z")==0){
    parameters_.disable("sc_Y_target");
  }
  else if(strcmp(name,"sc_Z_target_x")==0||
	  strcmp(name,"sc_Z_target_y")==0||
	  strcmp(name,"sc_Z_target_z")==0){
    parameters_.disable("sc_Z_target");
  } 
  else if(strcmp(name,"record_num")!=0){
    BurstDataError e("No special parameter called "+string(name),unknown_parameter);
    e.throwMe();
  }
}

double BurstData::getParam(const std::string& param_name) 
{
  if(!parameters_.isParam(param_name)){
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
  }
  Parameter* p=parameters_.getParam(param_name);
  double retval;
  Uvar* u;
  switch (p->datatype){
  case Parameter::ULONG:
    retval=(double)(*((unsigned int*)p->outpt_ptr));
    break;
 
  case Parameter::FLOAT:
    retval=(double)(*((float*)p->outpt_ptr));
    break;
 
  case Parameter::DOUBLE: 
    retval= *((double*)p->outpt_ptr);
    break;

  case Parameter::UVAR:
  case Parameter::UVAR_DOUBLE:
  case Parameter::UVAR_WEST_LON:
    u=(Uvar*)p->outpt_ptr;
    retval=u->getInUnits(p->units);
    break;
  
  default: 
    BurstDataError e("GetParam only works with scalar parameters");
    e.throwMe();
    
  }
  return(retval);
}

// p = getParameter(param_name)
//
// Return a pointer to the Parameter object for the indicated parameter
// in this BurstData object.
// Throws an error if param name is invalid.

Parameter* BurstData::getParameter(const std::string& param_name) 
  {
  if(!parameters_.isParam(param_name))
    {
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
    }
  return(parameters_.getParam(param_name));
  }

void BurstData::readParams(){
  int fp=file_.getPosition();
  
  for (ParameterList::PARAM_LIST::const_iterator p=parameters_.param_list.begin();
       p != parameters_.param_list.end(); p++)
    {
      if(p->enabled){
	if(p->datatype==Parameter::FULL){
	  readRecord();
	  break;
	}
	file_.setPosition(fp+p->offset);
	p->read(file_);
      }
    }
  fp=fp+record_size_;
  file_.setPosition(fp);
  Parameter* p=parameters_.getParam("t");
  if(p->enabled){
    sc_state_J2000.setTime(t);
    sc_state_target.setTime(t);
    
    sc_X_J2000.setTime(t);
    sc_Y_J2000.setTime(t);
    sc_Z_J2000.setTime(t);

    sc_X_target.setTime(t);
    sc_Y_target.setTime(t);
    sc_Z_target.setTime(t);

    rot_vel_J2000.setTime(t);
    rot_vel_target.setTime(t);

  }  
}


//-------------------------------------
// Read a single parameter back from the BurstData file
// and move file pointer to beginning of record
//---------------------------------------
void 
BurstData::readParameter(const string& param_name){
  int fp=file_.getPosition();//assume pointer at the beginning of the record
  if((fp-header_size_)%record_size_!=0) {
    ErrorMessage e("BurstData::readParameter pointer not at the beginning of the record");
    e.throwMe();
  }
  if(!parameters_.isParam(param_name)){
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
  }
  Parameter* p=parameters_.getParam(param_name);
  file_.setPosition(fp+p->offset);
  p->read(file_); // Stores parameter from file in appropriate place in 
                  // BurstData object
  file_.setPosition(fp);
}

//-------------------------------------
//Write a single parameter back into BurstData
// and return file pointer at the beginning of record
//---------------------------------------
void 
BurstData::writeParameter(const string& param_name,
			  Uvar value, const string& ustr){
  int fp=file_.getPosition();//assume pointer at the beginning of the record
  if((fp-header_size_)%record_size_!=0) {
    ErrorMessage e("BurstData::writeParameter pointer not at the beginning of the record");
    e.throwMe();
  }
  if(!parameters_.isParam(param_name)){
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
  }
  Parameter* p=parameters_.getParam(param_name);

  //as of now, only work with Uvar parameter type
  // does not work with Uvar_double or Uvar_west_lon data type
  if(p->datatype==Parameter::UVAR){
    file_.setPosition(fp+p->offset);
    file_.writeFloatInUnits(value,ustr);
  }
  else{
    BurstDataError e("no writing method provided");
    e.throwMe();
  }
  file_.setPosition(fp);
}
void 
BurstData::writeParameter(const string& param_name, const unsigned int& i){
  int fp=file_.getPosition();//assume pointer at the beginning of the record
  if((fp-header_size_)%record_size_!=0) {
    ErrorMessage e("BurstData::writeParameter pointer not at the beginning of the record");
    e.throwMe();
  }
  if(!parameters_.isParam(param_name)){
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
  }
  Parameter* p=parameters_.getParam(param_name);
  if(p->datatype !=Parameter::ULONG && p->length != sizeof(unsigned int)){
    BurstDataError e("variable type is not unsigned int or  size is not that of unsigned int ");
    e.throwMe();
  }
  file_.setPosition(fp+p->offset);
  file_.write(i);
  file_.setPosition(fp);
}

void
BurstData::writeParameter(const string& param_name, const float& x){
  int fp=file_.getPosition();//assume pointer at the beginning of the record
  if((fp-header_size_)%record_size_!=0) {
    ErrorMessage e("BurstData::writeParameter pointer not at the beginning of the record");
    e.throwMe();
  }
  if(!parameters_.isParam(param_name)){
    BurstDataError e("Cannot find parameter "+param_name,unknown_parameter);
    e.throwMe();
  }
  Parameter* p=parameters_.getParam(param_name);
  if(p->datatype !=Parameter::FLOAT && p->length != sizeof(float)){
    BurstDataError e("variable type is not float ");
    e.throwMe();
  }
  file_.setPosition(fp+p->offset);
  file_.write(x);
  file_.setPosition(fp);
}

//-------------------------------------
//Read sab_counter, which is unique for each burst,
// and set the file pointer at the beginning of the record
//---------------------------------------
unsigned int
BurstData::getSabCounter(){
  int fp=file_.getPosition();//assume pointer at the beginning of the record
  if((fp-header_size_)%record_size_!=0) {
    ErrorMessage e("BurstData::writeParameter pointer not at the beginning of the record");
    e.throwMe();
  }
  Parameter *p = parameters_.getParam("sab_counter");
  file_.setPosition(fp+p->offset);
  unsigned int retval;
  file_.read(retval);
  file_.setPosition(fp);
  return(retval);
}


// virtual method
double
BurstData::getSpecialParam(const char* name){
  static double record_num=0;
  static unsigned int sclk_start=sclk;
  double retval;
  if(strcmp(name,"rel_time")==0){
    retval=double(sclk-sclk_start);
  }
  else if(strcmp(name,"record_num")==0){
    record_num++;
    retval=record_num;
  }
  else  if(strcmp(name,"dist_from_earth")==0){
    PositionVector pos=sc_state_J2000.position();
    Uvar mag=pos.magnitude();
    retval=mag.getInUnits("km");
  }
  else if(strcmp(name,"dist_to_target")==0){
    PositionVector pos=sc_state_target.position();
    Uvar mag=pos.magnitude();
    retval=mag.getInUnits("km");
  }
  else if(strcmp(name,"speed_from_earth")==0){
    FloatVector vel=sc_state_J2000.velocity();
    PositionVector pos=sc_state_J2000.position();
    Uvar speed=dot(vel,pos)/pos.magnitude();
    retval=speed.getInUnits("km/s");
  }  
  else if(strcmp(name,"speed_to_target")==0){
    FloatVector vel=sc_state_target.velocity();
    PositionVector pos=sc_state_target.position();
    Uvar speed=-(dot(vel,pos)/pos.magnitude());
    retval=speed.getInUnits("km/s");
  }  
  else if(strcmp(name,"angle_from_sun")==0){
    if(goodGeometry()){
      PositionVector pos=sc_state_J2000.position();
      Frame sun_frame("J2000","Sun");
      PositionVector sun_pos= PositionVector("sun_pos",sun_frame,t,0,0,0);
      sun_pos=sun_pos.representIn(fsc_);
      DirectionVector z_axis=sc_Z_J2000.representIn(fsc_);
      Uvar angle=sun_pos.angle(-z_axis);
      retval=angle.getInUnits("deg");
    }
    else retval=-1000.0;
  }
  else if (strcmp(name,"sc_in_target_frame_position_x")==0) {
    retval=sc_state_target.position()[PositionVector::X].getInUnits("km");
  }
  else if (strcmp(name,"sc_in_target_frame_position_y")==0) {
    retval=sc_state_target.position()[PositionVector::Y].getInUnits("km");
  }
  else if (strcmp(name,"sc_in_target_frame_position_z")==0) {
    retval=sc_state_target.position()[PositionVector::Z].getInUnits("km");
  }
  else if (strcmp(name,"sc_in_target_frame_velocity_x")==0) {
    retval=sc_state_target.velocity()[FloatVector::X].getInUnits("km/s");
  }
  else if (strcmp(name,"sc_in_target_frame_velocity_y")==0) {
    retval=sc_state_target.velocity()[FloatVector::Y].getInUnits("km/s");
  }
  else if (strcmp(name,"sc_in_target_frame_velocity_z")==0) {
    retval=sc_state_target.velocity()[FloatVector::Z].getInUnits("km/s");
  }
  else if (strcmp(name,"sc_X_target_x")==0) {
    retval=sc_X_target[DirectionVector::X];
  }
  else if (strcmp(name,"sc_X_target_y")==0) {
    retval=sc_X_target[DirectionVector::Y];
  }
  else if (strcmp(name,"sc_X_target_z")==0) {
    retval=sc_X_target[DirectionVector::Z];
  }
  else if (strcmp(name,"sc_Y_target_x")==0) {
    retval=sc_Y_target[DirectionVector::X];
  }
  else if (strcmp(name,"sc_Y_target_y")==0) {
    retval=sc_Y_target[DirectionVector::Y];
  }
  else if (strcmp(name,"sc_Y_target_z")==0) {
    retval=sc_Y_target[DirectionVector::Z];
  }
  else if (strcmp(name,"sc_Z_target_x")==0) {
    retval=sc_Z_target[DirectionVector::X];
  }
  else if (strcmp(name,"sc_Z_target_y")==0) {
    retval=sc_Z_target[DirectionVector::Y];
  }
  else if (strcmp(name,"sc_Z_target_z")==0) {
    retval=sc_Z_target[DirectionVector::Z];
  }


  else{
    BurstDataError e("No special parameter called "+string(name),unknown_parameter);
    e.throwMe();
  }
  return(retval);
}

//----------------
// Other Methods
//-----------------

//----------------------------------------------------------------------
// recordCount()
//
// Return the number of records in the associated BurstData file.
// Currently, this method overwrites the public BurstData data fields and is
// not a const method even though it could be if the current data were
// saved and restored.
// Throw an exception if a file access error occurs.
//----------------------------------------------------------------------

unsigned int BurstData::recordCount()
  {
  if (! header_handled_) readHeader();
  if (!records_counted_)
    {
      BurstDataError e("Records not yet counted for "+ filename_, 
		       BurstData::read_error);     
      e.throwMe();
    }

  return(record_count_);
  }

//----------------------------------------
// Return Start Time
//----------------------------------------
Time BurstData::getStartTime(){
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't perform slckToRecord for output file " 
		       + filename_,
		       BurstData::read_error);
      e.throwMe();
    }
  int start_position = file_.getPosition();  // remember current position
  int start_file_number = file_number_;     // remember current file number
  unsigned int start_sclk=sclk;        // remember current sclk

  setFileNumber(0);                         // get base file
  file_.rewind();                           // position at file beginning
    
  header_handled_=false;               // read header
  readHeader();             

  readSclkOnly();
  
  Time retval("Cassini",sclk);  

  setFileNumber(start_file_number);
  file_.setPosition(start_position);         // restore original position
  sclk=start_sclk;   // restore original sclk
  return(retval);
}


//----------------------------------------
// Return End Time
//----------------------------------------
Time BurstData::getEndTime(){
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't perform slckToRecord for output file " 
		       + filename_,
		       BurstData::read_error);
      e.throwMe();
    }
  int start_position = file_.getPosition();  // remember current position
  int start_file_number = file_number_;     // remember current file number
  unsigned int start_sclk=sclk;        // remember current sclk

  setFileNumber(0);                         // get base file
  file_.rewind();                           // position at file beginning
    
  header_handled_=false;               // read header
  readHeader();             

  while(!eof()){
    readSclkOnly();  
  }

  Time retval("Cassini",sclk);  

  setFileNumber(start_file_number);
  file_.setPosition(start_position);         // restore original position
  sclk=start_sclk;   // restore original sclk
  return(retval);
}

//----------------------------------------------------------------------
// sclkToRecord(unsigned int sclk_val, bool reset_position=true)
//
// Return the number of the first record in the associated BurstData file with a 
// sclk which is greater than sclk_val.
//
// Resets file number and position in file
//
// Currently, this method overwrites the public BurstData data fields and is
// not a const method.
// Throw an exception if a file access error occurs.
//----------------------------------------------------------------------

unsigned int BurstData::sclkToRecord(unsigned int sclk_val)
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      BurstDataError e("Can't perform slckToRecord for output file " 
		       + filename_,
		       BurstData::read_error);
      e.throwMe();
    }
  int start_position = file_.getPosition();  // remember current position
  int start_file_number = file_number_;     // remember current file number
  setFileNumber(0);                         // get base file
  file_.rewind();                           // position at file beginning
    
  header_handled_=false;               // read header
  readHeader();             
           
  unsigned int count = 0;
  
  while (!eof())
    {  // read and count records
      readSclkOnly();
      if(sclk>sclk_val) break;
      count++;
    }
  setFileNumber(start_file_number);
  file_.setPosition(start_position);         // restore original position

  return(count);
  }

//---------------------------
// setFileNumber
//
// opens file number fno
// to beginning of file
// and reads header
//---------------------------

void BurstData::setFileNumber(int fno)
{
  file_number_=fno-1;
  openNextFile();
}

void BurstData::editTime(const Time& t0){
  t=t0;
  sclk=t.sclk("Cassini");
  
  // get offset between burst start time and sclk
  Time sclk_time("Cassini",sclk);
  brst=t.et()-sclk_time.et();
  computeGeometry();
}

//--------------------------
// Write Geometry
//--------------------------

void BurstData::writeGeometry()
{
  file_.write(quality_flag);
  // Time
  writeTime();

  // transmit time offset: write method OK for base unit
  transmit_time_offset.write(file_,0,"s");
  time_from_closest_approach.write(file_,0,"s");
  time_from_epoch.write(file_,0,"s");

  //target name and target body fixe frame
  file_.write(target_name);
  file_.write(tbf_frame_name);

  //pol and rotation angle
  
  file_.writeDoubleInUnits(pole_right_ascension,"deg");
  file_.writeDoubleInUnits(pole_declination ,"deg");
  file_.writeDoubleInUnits(target_rotation_rate,"deg/s");
  file_.writeDoubleInUnits(target_rotation_angle,"deg");


  // External Temperature telemetry
  T_scwg.writeFloat(file_,0,"K");
  T_feed.writeFloat(file_,0,"K");
  T_hga.writeFloat(file_,0,"K");

  // beam numbers
  file_.write(beam_number);

  // spacecraft state vectors
  double x,y,z;
  PositionVector scpos=sc_state_J2000.position();
  x=scpos[PositionVector::X].getInUnits("km");
  y=scpos[PositionVector::Y].getInUnits("km");
  z=scpos[PositionVector::Z].getInUnits("km");
  file_.writeXYZ(x,y,z);
  FloatVector scvel=sc_state_J2000.velocity();
  x=scvel[FloatVector::X].getInUnits("km/s");
  y=scvel[FloatVector::Y].getInUnits("km/s");
  z=scvel[FloatVector::Z].getInUnits("km/s");
  file_.writeXYZ(x,y,z);  

  PositionVector scpos2=sc_state_target.position();
  x=scpos2[PositionVector::X].getInUnits("km");
  y=scpos2[PositionVector::Y].getInUnits("km");
  z=scpos2[PositionVector::Z].getInUnits("km");
  file_.writeXYZ(x,y,z);
  FloatVector scvel2=sc_state_target.velocity();
  x=scvel2[FloatVector::X].getInUnits("km/s");
  y=scvel2[FloatVector::Y].getInUnits("km/s");
  z=scvel2[FloatVector::Z].getInUnits("km/s");
  file_.writeXYZ(x,y,z);  
  
  // spacecraft coordinate axes
  x=sc_X_J2000[DirectionVector::X];
  y=sc_X_J2000[DirectionVector::Y];
  z=sc_X_J2000[DirectionVector::Z];
  file_.writeXYZ(x,y,z);
  x=sc_Y_J2000[DirectionVector::X];
  y=sc_Y_J2000[DirectionVector::Y];
  z=sc_Y_J2000[DirectionVector::Z];
  file_.writeXYZ(x,y,z);
  x=sc_Z_J2000[DirectionVector::X];
  y=sc_Z_J2000[DirectionVector::Y];
  z=sc_Z_J2000[DirectionVector::Z];
  file_.writeXYZ(x,y,z);

  x=sc_X_target[DirectionVector::X];
  y=sc_X_target[DirectionVector::Y];
  z=sc_X_target[DirectionVector::Z];
  file_.writeXYZ(x,y,z);

  x=sc_Y_target[DirectionVector::X];
  y=sc_Y_target[DirectionVector::Y];
  z=sc_Y_target[DirectionVector::Z];
  file_.writeXYZ(x,y,z);

  x=sc_Z_target[DirectionVector::X];
  y=sc_Z_target[DirectionVector::Y];
  z=sc_Z_target[DirectionVector::Z];
  file_.writeXYZ(x,y,z);

  //rotational velocities
  x=rot_vel_J2000[RotationalVelocity::X].getInUnits("deg/s");
  y=rot_vel_J2000[RotationalVelocity::Y].getInUnits("deg/s");
  z=rot_vel_J2000[RotationalVelocity::Z].getInUnits("deg/s");
  file_.writeXYZ(x,y,z);

 
  x=rot_vel_target[RotationalVelocity::X].getInUnits("deg/s");
  y=rot_vel_target[RotationalVelocity::Y].getInUnits("deg/s");
  z=rot_vel_target[RotationalVelocity::Z].getInUnits("deg/s");
  file_.writeXYZ(x,y,z);

  norm_cnt_rl.writeFloat(file_,0,"1/s");
  norm_cnt_nd.writeFloat(file_,0,"1/s");
  norm_cnt_radio.writeFloat(file_,0,"1/s");

}


//--------------------------
// Read Geometry
//--------------------------

void BurstData::readGeometry(FileMgr* f)
{

  if(f==NULL) f=&file_;
  f->read(quality_flag);

  // Time
  readTime(f);

 

  //other times
  transmit_time_offset.read(*f,0,"s");
  time_from_closest_approach.read(*f,0,"s");
  time_from_epoch.read(*f,0,"s");

  //target name and target body fixe frame
  f->read(target_name);
  f->read(tbf_frame_name);

  //pol and rotation angle
  f->readDoubleInUnits(pole_right_ascension,degtorad,"rad");
  f->readDoubleInUnits(pole_declination, degtorad,"rad");
  f->readDoubleInUnits(target_rotation_rate,degtorad,"rad/s");
  f->readDoubleInUnits(target_rotation_angle,degtorad,"rad");
 
  // External Temperature telemetry
  T_scwg.readFloat(*f,0,"K");
  T_feed.readFloat(*f,0,"K");
  T_hga.readFloat(*f,0,"K");

 // beam numbers
  f->read(beam_number);

  // spacecraft state vectors
  double x,y,z;
  f->readXYZ(x,y,z);
  PositionVector scpos("scpos",j2000_,t,Uvar(x,"km"),Uvar(y,"km"),
		       Uvar(z,"km"));
  f->readXYZ(x,y,z);
  FloatVector scvel("scvel",j2000_,t,Uvar(x,"km/s"),Uvar(y,"km/s"),
		       Uvar(z,"km/s"));
  sc_state_J2000.setState(scpos,scvel);


  f->readXYZ(x,y,z);
  PositionVector scpos2("scpos2",ftarget_,t,Uvar(x,"km"),Uvar(y,"km"),
		       Uvar(z,"km"));
  f->readXYZ(x,y,z);
  FloatVector scvel2("scvel2",ftarget_,t,Uvar(x,"km/s"),Uvar(y,"km/s"),
		       Uvar(z,"km/s"));
  if(target_specified_) sc_state_target.setState(scpos2,scvel2);
  
  // spacecraft coordinate axes
  f->readXYZ(x,y,z);
  sc_X_J2000=DirectionVector("",j2000_,t,x,y,z);
  f->readXYZ(x,y,z);
  sc_Y_J2000=DirectionVector("",j2000_,t,x,y,z);
  f->readXYZ(x,y,z);
  sc_Z_J2000=DirectionVector("",j2000_,t,x,y,z);

  f->readXYZ(x,y,z);
  sc_X_target=DirectionVector("",ftarget_,t,x,y,z);
  f->readXYZ(x,y,z);
  sc_Y_target=DirectionVector("",ftarget_,t,x,y,z);
  f->readXYZ(x,y,z);
  sc_Z_target=DirectionVector("",ftarget_,t,x,y,z);

  //
  // Rotational Velocities
  //data have been written into deg
  // need to convert to rad by multiplying degtorad
  //
  f->readXYZ(x,y,z);
  rot_vel_J2000=RotationalVelocity("",j2000_,t,
				   Uvar(x*degtorad,"rad/s"),
				   Uvar(x*degtorad,"rad/s"), 
				   Uvar(z*degtorad,"rad/s"));
  f->readXYZ(x,y,z);
  rot_vel_target=RotationalVelocity("",ftarget_,t,
				    Uvar(x*degtorad,"rad/s"),
				    Uvar(x*degtorad,"rad/s"), 
				    Uvar(z*degtorad,"rad/s"));

  norm_cnt_rl.readFloat(*f,0,"1/s");
  norm_cnt_nd.readFloat(*f,0,"1/s");
  norm_cnt_radio.readFloat(*f,0,"1/s");

  
}



// read Time
void
BurstData::readTime(FileMgr* f){
  if(f==NULL) f=&file_;
  double t_sclk;
  //(0) read sclk
  f->read(t_sclk);
  Uvar t_et;
  //(1) ephemeris time
  t_et.read(*f,0,"s");
  t.setEt(t_et);
  //(2) ISOC format
  string utc_isoc;
  utc_isoc.resize(UTC_STRING_SIZE_ISOC+UTC_STRING_PAD_ISOC);
  f->read(utc_isoc);
  //(3) ISOD
  string utc_isod;
  utc_isod.resize(UTC_STRING_SIZE_ISOD+UTC_STRING_PAD_ISOD);
  f->read(utc_isod);
}

// writeTime
void
BurstData::writeTime(){  
  //(0) sclk
  Uvar t_sclk(t.encodedSclk("Cassini"),"s");//continuous sclk time
  t_sclk.write(file_,0,"s");
  //cout<<"t_sclk "<< t.encodedSclk("Cassini")<<endl;
  //(1) ephemeris time
  Uvar t_et=t.et();
  t_et.write(file_,0,"s");
  //(2) ISOC format + pad
  file_.write(t.utc("ISOC"));
  string pad;
  pad.resize(UTC_STRING_PAD_ISOC);
  for(int c=0;c<UTC_STRING_PAD_ISOC;c++) pad[c]=' ';
  file_.write(pad);
  //(3) ISOD format + pad
  file_.write(t.utc("ISOD"));
  pad.resize(UTC_STRING_PAD_ISOD);
  for(int c=0;c<UTC_STRING_PAD_ISOD;c++) pad[c]=' ';
  file_.write(pad);
}

// copies record from file of this object to bd object
int
BurstData::copyTo(BurstData& bd){
  return(bd.copy(*this));
}

// copies record from file of bd object to this object
int
BurstData::copy(BurstData& bd){
 if(bd.ft_ != ft_){
    // do nothing-- This should be OK
    // copying the BurstData part of an sdbr to an lbdr or vice versa
    // should be fine. I want the point target simulator to do this
    // to save space when it is in L1BFILE mode
 }
 
 if (bd.mode_ == "w" || bd.mode_ == "wb")
    {
      BurstData::BurstDataError e("Can't read from output file " + filename_,
	       read_error);
      e.throwMe();
    }  
 if (!bd.header_handled_)
    {
      BurstData::BurstDataError e("Can't read record until header is read",
	      BurstData::read_error);
      e.throwMe();
    } 

  // This check has the side effect of moving to the next  file if necessary
 if(bd.eof())
    {
      BurstData::BurstDataError e("Unexpected EOF in file "+bd.filename_,
	       read_error);
      e.throwMe();
    } 


  int start_position=bd.file_.getPosition();

  // Read in all public variables
  //convert this part to read
  //header slow fast data footer
 
  
  readPassiveSABData(&(bd.file_));  

  readGeometry(&(bd.file_)); 
  readMeasurementGeometry(&(bd.file_));

  data_read_=true;
  return(start_position);
} 

void BurstData::setQualityFlag(){
  if(!T_scwg_valid_) quality_flag+=4;
  if(!T_feed_valid_) quality_flag+=8;
  if(!T_hga_valid_) quality_flag+=16;
}

bool BurstData::goodGeometry(){
  return(quality_flag%4==0);
}


//----------------------------------
// Uses SAB Header/Footer Data to
// compute time of start of burst
//-------------------------------

void
BurstData::computeTime()
{
  if(! data_read_){
    BurstDataError e("computeTime: Sclk data not read");
    e.throwMe();
  }
  
  // Time reference is burst start time sclk+brst
  

  // Aggressively Warn when this keyword is set.
  if(use_artificial_sclk_offset_){
    cerr << "Warning TA/C43 Artificial time offset in use." << endl;
    cerr << "For normal operations YOU MUST remove the " << endl <<
      "TA_C43_ARTIFICIAL_SCLK_OFFSET keyword from the config file" << endl;
    cout << "Warning TA/C43 Artificial time offset in use." << endl;
    cout << "For normal operations YOU MUST remove the " << endl <<
      "TA_C43_ARTIFICIAL_SCLK_OFFSET keyword from the config file" << endl;
    sclk+=artificial_sclk_offset_;
  }
  t.setSclk("Cassini",sclk);
  t+=brst;


  
}

// Compute normalized radiometer data
void 
BurstData::computeNormalizedRadiometerData()
{

  //Y. Gim
  // cip, hip, and rip are reset by using offset values
  //problem with such scheme is that cip, hip,rip can be
  // negative when those values are zero filled due to 
  // BAD sab
  //So, when cip, hip, and rip have negative values,
  // I will assume it is due to zero-filled integration time
  // and therefore use zero for normalized count
  if(cip<Uvar(0,"s"))
    norm_cnt_rl=Uvar(0.0,"1/s");
  else
    norm_cnt_rl=(cnt_rl+radiometer_offset_)/cip;
   
  
  if(hip<Uvar(0,"s"))
    norm_cnt_nd=Uvar(0.0,"1/s");
  else
    norm_cnt_nd=(cnt_nd+radiometer_offset_)/hip;
 
  
  if(rip<Uvar(0,"s"))
    norm_cnt_radio=Uvar(0,"1/s");
  else
    norm_cnt_radio=(cnt_radio+radiometer_offset_*rad)/rip/rad;
	
  // old codes
  //reason for update
  // hip, cip, and rip have been adjusted already
  //norm_cnt_nd=(cnt_nd+radiometer_offset_)/(hip+radiometer_delta_tau_);
  //norm_cnt_rl=(cnt_rl+radiometer_offset_)/(cip+radiometer_delta_tau_);
  //norm_cnt_radio=((cnt_radio+radiometer_offset_*rad)/
  //	    (rip+radiometer_delta_tau_))/rad;
}

// Pretend SAB data was read
void
BurstData::pretendDataRead(){ data_read_=true; }



// Compute BurstData specific data at time t
void
BurstData::computeGeometry()
{ 
 if(! data_read_){
    BurstDataError e("ComputeGeometry: SAB data not read");
    e.throwMe();
  }
 if(! t.valid()){
    BurstDataError e("ComputeGeometry: time not computed", invalid_time);
    e.throwMe();
  }

 //----------------------
 //transmit time offset: 1ms +1PRI
 //---------------------
 transmit_time_offset=Uvar(1*mstos,"s")+pri;



 //--------------------------------
 // get beam number 
 //-------------------------------
 unsigned int i_ctbe= ctbe;
 unsigned int j=0;
 while(i_ctbe){
   i_ctbe>>=1;
   j++;
 }
 beam_number=j;



 //---------------------------------
 //microwave source scan default id
 // once microwave source scan is detected 
 // this value will set set to a target whose beam
 // distance from boresight is a minimum  in beam frame
 //----------------------------------
 mw_source_id=0;




 //--------------------
 //Compute epoch time, target name, target frame
 //-----------------------
 if(target_specified_){
   
   //--------------
   //time from closest approach
   //-------------
   time_from_closest_approach= t - closest_approach_time_;
   //----------------
   //time from epoch
   //--------------
   time_from_epoch = t - epoch_;
   //------------
   //target name
   //-------------
   unsigned int N=target_name_.size();
   if(N>= TARGET_NAME_STRING_SIZE){
     cout<<"ButstData.cpp::computeGeometry(): Target name is too long  and therefore its size ";
     cout<<"will be reduced"<<endl;
   }
   for(unsigned int c=0;c<TARGET_NAME_STRING_SIZE;++c)
     if(c<N) target_name[c]= target_name_[c];
   //------------
   //frame name
   //-------------
   string frame_name="IAU_"+target_name_;
   N= frame_name.size();
   if(N>=TARGET_FRAME_NAME_STRING_SIZE) {
     cout<<"ButstData.cpp::computeGeometry(): frame name is too long";
     cout<<"and its size will be reduced"<<endl;
   }
   for(unsigned int c=0;c<TARGET_FRAME_NAME_STRING_SIZE;++c)
     if(c<N) tbf_frame_name[c]=frame_name[c];
   //--------------
   //target ra, dec, rotation angle and rate
   //-------------
   pole_right_ascension= pole_right_ascension_epoch_;
   pole_declination = pole_declination_epoch_;
   target_rotation_angle = target_rotation_angle_epoch_;
   target_rotation_rate= target_rotation_rate_epoch_;
 }
 else {
   //--------------
   //time from closest approach
   //-------------
   time_from_closest_approach= t - epoch_;
   //----------------
   //time from epoch
   //--------------
   time_from_epoch = t - epoch_;  

   //--------------
   //target name specification
   //-------------- 
   for(unsigned int c=0;c<TARGET_NAME_STRING_SIZE;++c)
     target_name[c]=' ';
   for(unsigned int c=0;c<TARGET_FRAME_NAME_STRING_SIZE;++c)
     tbf_frame_name[c]=' ';
   pole_right_ascension=0.0;
   pole_declination = 0.0;
   target_rotation_angle =0.0;
   target_rotation_rate= 0.0;

   //------------------------------
   //if target is microwave source 
   //----------------------------
   if(mw_source_scan_ && beam_number!=0){
     vector<Uvar> beam_distance; beam_distance.clear();
     Uvar mw_azi, mw_elev;
     Frame fbeam=Frame("CASSINI_RADAR_"+toStr(beam_number),"Cassini");
     for(unsigned int i=0;i<Num_mw_sources_;++i){
       DirectionVector mw_source_dir("mc_source",j2000_,t,0,0,1);// J2000 centered at earth
       mw_source_dir.setRADEC(source_ra_[i],source_dec_[i]);
       mw_source_dir.representIn(fbeam);
       mw_source_dir.getAzimuthElevation(mw_azi,mw_elev);
       beam_distance.push_back( sqrt(mw_azi*mw_azi + mw_elev*mw_elev));
     }
     if(beam_distance.size()!=Num_mw_sources_) ErrorMessage("BurstData:cppbeam distance container size mismatch").throwMe();
     //find nearest target
     mw_source_id=0;
     Uvar min=beam_distance[0];
     for(unsigned int i=0;i<Num_mw_sources_;++i){
       if(min>= beam_distance[i]){
	 min=beam_distance[i];
	 mw_source_id=i;
       }//find min beam distance
     }//loop over microwave sources
  
   
     //------------
     //microwave source target: mw_source_id
     //-----------
     //cout<<"burst number "<< sab_counter<<endl;
     //cout<<"mw source target name "<< source_name_[mw_source_id]<<endl;
     
     //----------------
     //set  microwave source name as target name
     //----------------  
     string mw_source_name=source_name_[mw_source_id];
     unsigned int N=mw_source_name.size();
     if(N>= TARGET_NAME_STRING_SIZE) cout<<"ButstData.cpp::computeGeometry():microwave source name is too long  and therefore its size will be reduced"<<endl;
     for(unsigned int c=0;c<TARGET_NAME_STRING_SIZE;++c)
       if(c<N) target_name[c]= mw_source_name[c];
     
     //------------
     //frame name
     //-------------
     string frame_name="J2000_cassini";
     N= frame_name.size();
     if(N>=TARGET_FRAME_NAME_STRING_SIZE) cout<<"ButstData.cpp::computeGeometry(): frame name for microwave source scan is too long and its size will be reduced"<<endl;
     for(unsigned int c=0;c<TARGET_FRAME_NAME_STRING_SIZE;++c)
       if(c<N) tbf_frame_name[c]=frame_name[c];
   }//microwave scan and beam_number>0
 }//target is not a planetary object: either no target or microwave source

 
  

 // compute temperatures
 // Reset temperatures to 0
 T_scwg=Uvar(0,"K");
 T_feed=Uvar(0,"K");
 T_hga=Uvar(0,"K");
 
 if(T_scwg_available_){
   T_scwg_valid_=true;
   try{
     T_scwg=T_scwg_ephem_.interpolate(t,interpolation_valid_time_);
   }
   catch(ErrorMessage& e){
     cerr<<"Fail to retrieve scwg temp "<<e.msg<<endl;
     T_scwg_valid_=false;
   }
 }
 else{
   T_scwg_valid_=false;
 }
 
 if(T_feed_available_){
   T_feed_valid_=true;
   try{
     T_feed=T_feed_ephem_.interpolate(t, interpolation_valid_time_);
   }
   catch(ErrorMessage& e){
     cerr<<"Fail to retrieve feed temp "<<e.msg<<endl;
     T_feed_valid_=false;
   }
 }
 else{
   T_feed_valid_=false;
 }
 
 if(T_hga_available_){
   T_hga_valid_=true;
   try{
     T_hga=T_hga_ephem_.interpolate(t, interpolation_valid_time_);
   }
   catch(ErrorMessage& e){
     cerr<<"Fail to retrieve hga tmp "<<e.msg<<endl;
     T_hga_valid_=false;
   }
 }
 else{
   T_hga_valid_=false;
 }



 //--------------------------------------
 // compute spacecraft state vector, coordinate axes and
 // rotationalvelocity
 //---------------------------------------

 // in j2000
 double et;
 t.getEt(et);
 j2000_.ephemeris(sc_state_J2000,cassini_spice_id,et);
 sc_X_J2000=DirectionVector("sc_X_J2000",fsc_,t,1,0,0).representIn(j2000_);
 sc_Y_J2000=DirectionVector("sc_X_J2000",fsc_,t,0,1,0).representIn(j2000_);
 sc_Z_J2000=DirectionVector("sc_X_J2000",fsc_,t,0,0,1).representIn(j2000_);

 //compute rotational velocity
 rot_vel_J2000=RotationalVelocity(fsc_,j2000_,t);
 
 // If target is a planetary object
 if(target_specified_){
   ftarget_.ephemeris(sc_state_target,cassini_spice_id,et);
   sc_X_target=DirectionVector("sc_X_target",fsc_,t,1,0,0).representIn(ftarget_);
   sc_Y_target=DirectionVector("sc_Y_target",fsc_,t,0,1,0).representIn(ftarget_);
   sc_Z_target=DirectionVector("sc_Z_target",fsc_,t,0,0,1).representIn(ftarget_);
   rot_vel_target=RotationalVelocity(fsc_,ftarget_,t);  
 }
 //if target is microwave source
 if(mw_source_scan_ && beam_number!=0){
   DirectionVector mw_source_dir("mc_source",j2000_,t,0,0,1);// J2000 centered at earth
   mw_source_dir.setRADEC(source_ra_[mw_source_id],source_dec_[mw_source_id]);
   mw_source_dir.representIn(j2000_at_cassini_);//J2000 centered at Cassini
   PositionVector target_position= -mw_source_dir;// unit vector without distance
   StateVector sc;
   j2000_at_cassini_.ephemeris(sc,cassini_spice_id,et);//state vector in j2000_at_cassini_ frame
   sc_state_target.setState(target_position,sc.velocity());
 }

 computeNormalizedRadiometerData();
}


//compute measurement geometry for science data field 
void
BurstData::computeMeasurementGeometry(const Umat& azim_1way3dB_ellipse_fit,
				      const Umat& elev_1way3dB_ellipse_fit,
				      const Umat& azim_2way3dB_ellipse_fit,
				      const Umat& elev_2way3dB_ellipse_fit)
  {
  if(! data_read_){
    BurstDataError e("ComputeMeasurementGeometry: SAB data not read");
    e.throwMe();
  }
 if(! t.valid()){
    BurstDataError e("ComputeMeasurementGeometry: time not computed", invalid_time);
    e.throwMe();
  }
 if(azim_1way3dB_ellipse_fit.rows()!=5 || azim_1way3dB_ellipse_fit.cols()!=4){
   BurstDataError e("ComputeMeasurementGeometry: azim_1way3dB matrix is not 5x4", invalid_time);
   e.throwMe();
}
 if(azim_2way3dB_ellipse_fit.rows()!=5 || azim_2way3dB_ellipse_fit.cols()!=4){
   BurstDataError e("ComputeMeasurementGeometry: azim_2way3dB matrix is not 5x4", invalid_time);
   e.throwMe();
}
 if(elev_1way3dB_ellipse_fit.rows()!=5 || elev_1way3dB_ellipse_fit.cols()!=4){
   BurstDataError e("ComputeMeasurementGeometry: elev_1way3dB matrix is not 5x4", invalid_time);
   e.throwMe();
}
 if(elev_2way3dB_ellipse_fit.rows()!=5 || elev_2way3dB_ellipse_fit.cols()!=4){
   BurstDataError e("ComputeMeasurementGeometry: elev_2way3dB matrix is not 5x4", invalid_time);
   e.throwMe();
}



  pass_geom_time_offset=Uvar(0,"s");
  pass_pol_angle = Uvar(0,"rad");//default
  pass_emission_angle=Uvar(0,"rad");
  pass_azimuth_angle = Uvar(0,"rad");
  pass_centroid_lon=Uvar(0,"rad");
  pass_centroid_lat=Uvar(0,"rad");
  
  pass_major_width=Uvar(0,"km");
  pass_minor_width=Uvar(0,"km");

  pass_ellipse_pt1_lon=Uvar(0,"rad");
  pass_ellipse_pt2_lon=Uvar(0,"rad");
  pass_ellipse_pt3_lon=Uvar(0,"rad");
  pass_ellipse_pt4_lon=Uvar(0,"rad");
  pass_ellipse_pt1_lat=Uvar(0,"rad");
  pass_ellipse_pt2_lat=Uvar(0,"rad");
  pass_ellipse_pt3_lat=Uvar(0,"rad");
  pass_ellipse_pt4_lat=Uvar(0,"rad");
  
  act_geom_time_offset=Uvar(0,"s");
  act_pol_angle=Uvar(0,"rad");//default
  act_incidence_angle=Uvar(0,"rad");
  act_azimuth_angle=Uvar(0,"rad");
  act_centroid_lon=Uvar(0,"rad");
  act_centroid_lat=Uvar(0,"rad");

  act_major_width=Uvar(0,"km");
  act_minor_width=Uvar(0,"km");

  act_ellipse_pt1_lon=Uvar(0,"rad");
  act_ellipse_pt2_lon=Uvar(0,"rad");
  act_ellipse_pt3_lon=Uvar(0,"rad");
  act_ellipse_pt4_lon=Uvar(0,"rad");
  act_ellipse_pt1_lat=Uvar(0,"rad");
  act_ellipse_pt2_lat=Uvar(0,"rad");
  act_ellipse_pt3_lat=Uvar(0,"rad");
  act_ellipse_pt4_lat=Uvar(0,"rad");
  
  
  altimeter_profile_length=0;//default
  num_pulses_received=0;//default
  sar_azimuth_res=Uvar(0,"km");
  sar_range_res=Uvar(0,"km");
  sar_centroid_bidr_lon=Uvar(0,"rad");
  sar_centroid_bidr_lat=Uvar(0,"rad");

  //set science flag
  science_qual_flag=0;//start from 0
  bitset(science_qual_flag,2,1);//all altimeter fields invalid
  bitset(science_qual_flag,3,1);//all scatt mode fields invalid
  bitset(science_qual_flag,4,1);//all rad science (not geom) fields invalid
  bitset(science_qual_flag,9,1);//all sar fields invalid
  if(num_bursts_in_flight>1) bitset(science_qual_flag,10,1);//act geom invalid
  if(!target_specified_ || beam_number == 0){
    bitset(science_qual_flag,0,1);//passive mode  geom invalid
    bitset(science_qual_flag,1,1);//active mode geom invalid
    return;//no target or beam is used
  }

  //--------------------------------------------------
  //electric field direction: spacecraft x direction
  //-----------------------------------------------
 
  //time at which  polarization angle is computed 
  Time t_pol;
  
  //polarization direction
  DirectionVector pol;
  Frame fbeam("CASSINI_RADAR_"+toStr(beam_number),"Cassini");
  StateVector sc_state_pol;
  DirectionVector bore, look;
  DirectionVector dir_h, dir_v,dir_z;
  Uvar theta, phi;
  DirectionVector east, north;
  DirectionVector to_bore;
  DirectionVector major1, major2;
  DirectionVector minor1, minor2;
  Uvar angle_bet_directions;
  vector<Uvar> ranges_variation_inside_2way3dB;
  ranges_variation_inside_2way3dB.clear();
 
  // active mode polarization angle
  if(!(csr==0 ||csr==8)) {
    bitset(science_qual_flag,1,1);//active mode geom invalid
  }
  else{ //not in cal mode
    if(r_mode>=4 && r_mode <=7) 
      bitset(science_qual_flag,1,1);//active geom invalid
    else if(r_mode<4 || ( r_mode>7 && r_mode<12))
      {//active mode: not in radiometer,IGO , Earth viewing, or Bistatic
	
	//time
	Time t_pol = t;
	
	//mid point of transmission/2  and receiver/2
	act_geom_time_offset=(double(pul)*pri/2.0+rwd+(double(pul)*pri+tro)/2.0)/2.0;
	//time at which act measurement geometries are computed
	t_pol +=act_geom_time_offset;
	pol=DirectionVector("pol",fsc_,t_pol,1,0,0);//sc x-direction
	bore=DirectionVector("bore",fbeam,t_pol,0,0,1);
	ftarget_.ephemeris(sc_state_pol,"Cassini",t_pol,"NONE");
	TargetGeom tg(t_pol);
	tg.setState(sc_state_pol);
	tg.setTarget();
	tg.setLookDirection(bore);
	if(!tg.foundSurfaceIntercept()) 
	  bitset(science_qual_flag,7,1);//no boresight surface intercept
	else
	  {//if there is a surface intercept point for boresight
	    
	    //keep the range
	    ranges_variation_inside_2way3dB.push_back(tg.range());	    
	    //incidence angle
	    act_incidence_angle = tg.incidenceAngle();
	    
	    //polarization angle
	    dir_z = tg.surfaceIntercept();
	    dir_h = cross(bore, dir_z);
	    dir_v = cross(dir_h,bore);
	    act_pol_angle=Uvar(atan2(dot(dir_v,pol),dot(dir_h,pol)),"rad");
	    if(act_pol_angle<Uvar(0,"rad")) act_pol_angle += Uvar(2.0*pi,"rad");
	    
	    //azimuth angle
	    dir_z.getSpherical(theta,phi);
	    //----------------------------------------------
	    //East direction vector (-sin(phi),cos(phi),0)
	    //Target direction (cos(phi)sin(theta), sin(phi)sin(theat),cos(theta)
	    //------------------------------------------------
	    east=DirectionVector("east",ftarget_,t_pol,-sin(phi), cos(phi),0);
	    to_bore=DirectionVector("tobore",ftarget_,t_pol,cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
	    north = cross(to_bore, east);
	   
	    act_azimuth_angle = Uvar(atan2(dot(bore,north),dot(bore,east)),"rad");
	    if(act_azimuth_angle<Uvar(0,"rad")) act_azimuth_angle+=Uvar(2*pi,"rad");
	  
	    //lat lon
	    act_centroid_lon = tg.lon();
	    act_centroid_lat = tg.lat();
	  }
	//work on other than boresight: 2way 3dB
	for(unsigned int j=0;j<4;++j)
	  {
	    TargetGeom tg_look(t_pol);
	    tg_look.setState(sc_state_pol);
	    tg_look.setTarget();
	    look=DirectionVector("",fbeam,t_pol,0,0,1);
	    look.setAzimuthElevation(azim_2way3dB_ellipse_fit(beam_number-1,j),
				     elev_2way3dB_ellipse_fit(beam_number-1,j));
	    tg_look.setLookDirection(look);
	    if(!tg_look.foundSurfaceIntercept())
	      bitset(science_qual_flag,8,1);//no intercept
	    else {
	      //keep the range !!
	      ranges_variation_inside_2way3dB.push_back(tg.range());
	      switch(j)
		{
		case 0:
		  act_ellipse_pt1_lon= tg_look.lon();
		  act_ellipse_pt1_lat = tg_look.lat();
		  break;
		case 1:
		  act_ellipse_pt2_lon = tg_look.lon();
		  act_ellipse_pt2_lat = tg_look.lat();
		  break;
		case 2:
		  act_ellipse_pt3_lon = tg_look.lon();
		  act_ellipse_pt3_lat = tg_look.lat();
		  break;		    
		case 3:
		  act_ellipse_pt4_lon = tg_look.lon();
		  act_ellipse_pt4_lat = tg_look.lat();
		  break;
		default://should not happen
		  break;
		}
	    }//if there is a surface intercept point for 2way 3dB
	  }//four different look directions
	//compute act major width and minor width only when all the suraface
	// intercept points are valid
	if(bitget(science_qual_flag,8,8)==0){
	  major1=DirectionVector("",ftarget_,t_pol,0,0,1);
	  major2=DirectionVector("",ftarget_,t_pol,0,0,1);
	  minor1=DirectionVector("",ftarget_,t_pol,0,0,1);
	  minor2=DirectionVector("",ftarget_,t_pol,0,0,1);
	  
	  //set lat lon
	  major1.setPlanetodetic(act_ellipse_pt1_lat, act_ellipse_pt1_lon);
	  major2.setPlanetodetic(act_ellipse_pt2_lat, act_ellipse_pt2_lon);
	  angle_bet_directions=major1.angle(major2);
	  act_major_width= angle_bet_directions.getInUnits("rad")* target_radius_;
	  //set lat lon
	  minor1.setPlanetodetic(act_ellipse_pt3_lat, act_ellipse_pt3_lon);
	  minor2.setPlanetodetic(act_ellipse_pt4_lat, act_ellipse_pt4_lon);
	  angle_bet_directions=minor1.angle(minor2);
	  act_minor_width= angle_bet_directions.getInUnits("rad")*target_radius_;
	}
	//compute num_pulses_received using range variations
	// when multiple bursts in fligh, no computation, set to -1
	if(ranges_variation_inside_2way3dB.size()==0
	   || num_bursts_in_flight!=1)
	  num_pulses_received=0;//no interception, no echo
	else{
	  Uvar min_range, max_range;
	  min_max(min_range, max_range, ranges_variation_inside_2way3dB);//compute min max
	  Uvar time_min = 2.0* min_range/speed_light;
	  Uvar time_max = 2.0* max_range/speed_light + pri*double(pul);
	  int loss_beginning=0;
	  int loss_end=0;
	  if(time_min <rwd)
	    loss_beginning =(int)round_double(((rwd - time_min)/pri).getInUnits(""));
	  if(time_max >(rwd + pri*double(pul)+tro))
	    loss_end=(int) round_double(((time_max-rwd-pri*double(pul)-tro)/pri).getInUnits(""));
	  num_pulses_received = (int) pul - (loss_beginning+loss_end);
	  if( (loss_end+loss_beginning) > (int) pul){
	    cout<<"---- Warning: check here: too many pulse losses "<<endl;
	    cout<<"No echo received "<<endl;
	    cout<<"lost pulses "<< loss_end + loss_beginning<<endl;
	    cout<<"loss at beginning and end "<< loss_beginning<<endl;
	    cout<<"loss at end "<< loss_end<<endl;
	    cout<<"first pulse "<< time_min<<endl;
	    cout<<"last pulse "<< time_max<<endl;
	    cout<<"time offset "<<	act_geom_time_offset<<endl;
	    cout<<"bdp "<< bpd<<endl;
	    cout<<"rwd "<< rwd<<endl; 
	    num_pulses_received=0;
	  }
	 
	}
      }//active, pulse is transmitted
    else
      {//rmode is larger than 12
	ErrorMessage("BurstData.cpp:rmode >= 12").throwMe();
      }
   


  }//not in calibration mode
  
  // passive mode polarization angle
  if(rad==0) bitset(science_qual_flag,0,1);//rad mode geom fields invalid
  else 
    {//whenever there is a radiometer count window
      t_pol = t;
      pass_geom_time_offset= rwd +double(pul)*pri+tro+ hip + cip;//right before radiometer measurements
      pass_geom_time_offset =  rad * rip/2.0;//half way
      t_pol +=pass_geom_time_offset;
      pol=DirectionVector("pol",fsc_,t_pol,1,0,0);
      bore=DirectionVector("bore",fbeam,t_pol,0,0,1);
      ftarget_.ephemeris(sc_state_pol,"Cassini",t_pol,"NONE");
      TargetGeom tg(t_pol);
      tg.setState(sc_state_pol);
      tg.setTarget();
      tg.setLookDirection(bore);
      if(!tg.foundSurfaceIntercept()) 
	bitset(science_qual_flag,5,1);//no boresight interception
      else 
	{
	  //emission angle
	  pass_emission_angle = tg.incidenceAngle();
	  
	  //polarization angle
	  dir_z = tg.surfaceIntercept();
	  dir_h = cross(bore, dir_z);
	  dir_v = cross(dir_h,bore);
	  pass_pol_angle=Uvar(atan2(dot(dir_v,pol),dot(dir_h,pol)),"rad");
	  if(pass_pol_angle <Uvar(0,"rad")) pass_pol_angle +=Uvar(2.0*pi,"rad");
	  //cout<<"pass pol angle "<< pass_pol_angle.getInUnits("deg")<<endl;
	  
	  //azimuth angle
	  dir_z.getSpherical(theta,phi);
	  east=DirectionVector("east",ftarget_,t_pol,-sin(phi), cos(phi),0);
	  to_bore=DirectionVector("tobore",ftarget_,t_pol,cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
	  north = cross(to_bore, east);
	   
	 
	  pass_azimuth_angle = Uvar(atan2(dot(bore,north),dot(bore,east)),"rad");
	  if(pass_azimuth_angle<Uvar(0,"rad")) pass_azimuth_angle+=Uvar(2*pi,"rad");

	  //cout<<"pass east "<< east<<endl;
	  //cout<<"pass north "<<north<<endl;
	  //cout<<"look direction "<< bore.representIn(ftarget_)<<endl;
	  //cout<<"dot(look north "<< dot(bore, north)<<endl;
	  //cout<<"dot(look east "<< dot(bore,east)<<endl;
	  //cout<<"pass azimuth "<< pass_azimuth_angle.getInUnits("deg")<<endl;

	  //lat lon
	  pass_centroid_lon = tg.lon();
	  pass_centroid_lat = tg.lat();
	}
      //work on other than boresight: 1way 3dB
      for(unsigned int j=0;j<4;++j)
	{
	  TargetGeom tg_look(t_pol);
	  tg_look.setState(sc_state_pol);
	  tg_look.setTarget();
	  look=DirectionVector("",fbeam,t_pol,0,0,1);
	  look.setAzimuthElevation(azim_1way3dB_ellipse_fit(beam_number-1,j),
				   elev_1way3dB_ellipse_fit(beam_number-1,j));
	  tg_look.setLookDirection(look);
	  if(!tg_look.foundSurfaceIntercept())
	    bitset(science_qual_flag,6,1);//no intercept
	  else{
	    switch(j)
	      {
	      case 0:
		pass_ellipse_pt1_lon= tg_look.lon();
		pass_ellipse_pt1_lat = tg_look.lat();
		break;
	      case 1:
		pass_ellipse_pt2_lon = tg_look.lon();
		pass_ellipse_pt2_lat = tg_look.lat();
		break;
	      case 2:
		pass_ellipse_pt3_lon = tg_look.lon();
		pass_ellipse_pt3_lat = tg_look.lat();
		break;		    
	      case 3:
		pass_ellipse_pt4_lon = tg_look.lon();
		pass_ellipse_pt4_lat = tg_look.lat();
		break;
	      default://should not happen
		break;
	      }
	  }//if there is a surface intercept point for 1way 3dB
	}//four different look directions
      if(bitget(science_qual_flag,6,6)==0){
	//compute act major width and minor widht
	//compute act major width and minor widht
	major1=DirectionVector("",ftarget_,t_pol,0,0,1);
	major2=DirectionVector("",ftarget_,t_pol,0,0,1);
	minor1=DirectionVector("",ftarget_,t_pol,0,0,1);
	minor2=DirectionVector("",ftarget_,t_pol,0,0,1);
	
	//set lat lon
	major1.setPlanetodetic(pass_ellipse_pt1_lat, pass_ellipse_pt1_lon);
	major2.setPlanetodetic(pass_ellipse_pt2_lat, pass_ellipse_pt2_lon);
	angle_bet_directions=major1.angle(major2);
	pass_major_width= angle_bet_directions.getInUnits("rad")* target_radius_;
	//set lat lon
	minor1.setPlanetodetic(pass_ellipse_pt3_lat, pass_ellipse_pt3_lon);
	minor2.setPlanetodetic(pass_ellipse_pt4_lat, pass_ellipse_pt4_lon);
	angle_bet_directions=minor1.angle(minor2);
	pass_minor_width= angle_bet_directions.getInUnits("rad")*target_radius_;
      }
    }//when rad != 0
  }

//read measurement geometry
void 
BurstData:: readMeasurementGeometry(FileMgr* f){

  if(f==NULL) f=&file_;

  //reset longitude  container
  east_west_.resize(5);   //reset

  f->read(science_qual_flag);
  system_gain.readFloat(*f,0,"1/(s K)");
  antenna_brightness_temp.readFloat(*f,0,"K");
  system_noise_temp.readFloat(*f,0,"K");
  abt_std.readFloat(*f,0,"K");
  
  pass_geom_time_offset.readFloat(*f,0,"s");

  f->readFloatInUnits(pass_pol_angle,degtorad,"rad");
  f->readFloatInUnits(pass_emission_angle,degtorad,"rad");
  f->readFloatInUnits(pass_azimuth_angle,degtorad,"rad");
  f->readFloatInUnits(east_west_[0],degtorad,"rad");
  f->readFloatInUnits(pass_centroid_lat,degtorad,"rad");

  
  pass_major_width.readFloat(*f,0,"km");
  pass_minor_width.readFloat(*f,0,"km");
  
  f->readFloatInUnits(east_west_[1],degtorad,"rad");
  f->readFloatInUnits(east_west_[2],degtorad,"rad");
  f->readFloatInUnits(east_west_[3],degtorad,"rad");
  f->readFloatInUnits(east_west_[4],degtorad,"rad");
  
  //------------------------------------------  
  //Convert West bound lon to East bound lon
  // from 0 to 180: from -0 to -180
  // from 180 to 360: from 180 to 0 
  //---------------------------------------
  for(unsigned int i=0;i<5;++i){
    if(east_west_[i]>=zero_angle_ && east_west_[i] <= pi_angle_)
      east_west_[i] = -east_west_[i];
    else if(east_west_[i] > pi_angle_ && east_west_[i] <= two_pi_angle_)
      east_west_[i] = two_pi_angle_ -east_west_[i];
    else{
      BurstDataError e("passive lon value is out of range");
      e.throwMe();
    }
  }
  pass_centroid_lon = east_west_[0];
  pass_ellipse_pt1_lon = east_west_[1];
  pass_ellipse_pt2_lon = east_west_[2];
  pass_ellipse_pt3_lon = east_west_[3];
  pass_ellipse_pt4_lon = east_west_[4];

  f->readFloatInUnits(pass_ellipse_pt1_lat,degtorad,"rad");  
  f->readFloatInUnits(pass_ellipse_pt2_lat,degtorad,"rad");
  f->readFloatInUnits(pass_ellipse_pt3_lat,degtorad,"rad");
  f->readFloatInUnits(pass_ellipse_pt4_lat,degtorad,"rad");

  
  f->read(num_pulses_received);
  
  f->readFloatInUnits(total_echo_energy,WtoMW,"MJ");
  
  f->readFloatInUnits(noise_echo_energy, WtoMW,"MJ");
  f->readFloatInUnits(x_factor,WtoMW,"MW");

  sigma0_uncorrected.readFloat(*f,0,"");
  sigma0_corrected.readFloat(*f,0,"");
  sigma0_uncorrected_std.readFloat(*f,0,"");

  altitude_means.readFloat(*f,0,"km");
  altitude_means_std.readFloat(*f,0,"km");
  
  act_geom_time_offset.readFloat(*f,0,"s");

  //reset east west container
  east_west_.resize(5);//reset
  f->readFloatInUnits(  act_pol_angle,degtorad,"rad");
  f->readFloatInUnits(act_incidence_angle,degtorad,"rad");
  f->readFloatInUnits(act_azimuth_angle,degtorad,"rad");  
  f->readFloatInUnits(east_west_[0],degtorad,"rad");
  f->readFloatInUnits(act_centroid_lat,degtorad,"rad");

  act_major_width.readFloat(*f,0,"km");
  act_minor_width.readFloat(*f,0,"km");
  
  
  f->readFloatInUnits(east_west_[1],degtorad,"rad");
  f->readFloatInUnits(east_west_[2],degtorad,"rad");
  f->readFloatInUnits(east_west_[3],degtorad,"rad");
  f->readFloatInUnits(east_west_[4],degtorad,"rad");
  //-----------------------------------------
  //Convert West bound lon to East bound lon
  // from 0 to 180: from -0 to -180
  // from 180 to 360: from 180 to 0 
  //--------------------------------------
  for(unsigned int i=0;i<5;++i){
    if(east_west_[i]>=zero_angle_ && east_west_[i] <= pi_angle_)
      east_west_[i] = -east_west_[i];
    else if(east_west_[i] > pi_angle_ && east_west_[i] <= two_pi_angle_)
      east_west_[i] = two_pi_angle_ -east_west_[i];
    else{
      BurstDataError e("active lon value is out of range");
      e.throwMe();
    }
  }
  act_centroid_lon = east_west_[0];
  act_ellipse_pt1_lon = east_west_[1];
  act_ellipse_pt2_lon = east_west_[2];
  act_ellipse_pt3_lon = east_west_[3];
  act_ellipse_pt4_lon = east_west_[4];

  f->readFloatInUnits(act_ellipse_pt1_lat,degtorad,"rad"); 
  f->readFloatInUnits(act_ellipse_pt2_lat,degtorad,"rad");
  f->readFloatInUnits(act_ellipse_pt3_lat,degtorad,"rad");
  f->readFloatInUnits(act_ellipse_pt4_lat,degtorad,"rad");




  altimeter_profile_range_start.readFloat(*f,0,"km");
  altimeter_profile_range_step.readFloat(*f,0,"km");
  f->read(altimeter_profile_length);
  
  sar_azimuth_res.readFloat(*f,0,"km");
  sar_range_res.readFloat(*f,0,"km");
  
  f->readFloatInUnits(sar_centroid_bidr_lon,degtorad,"rad");
  f->readFloatInUnits(sar_centroid_bidr_lat,degtorad,"rad");
}

//write measurement geometry
void 
BurstData::writeMeasurementGeometry(){
 
  file_.write(science_qual_flag);
  system_gain.writeFloat(file_,0,"1/(s K)");
  antenna_brightness_temp.writeFloat(file_,0,"K");
  system_noise_temp.writeFloat(file_,0,"K");
  abt_std.writeFloat(file_,0,"K");
  
  pass_geom_time_offset.writeFloat(file_,0,"s");

  file_.writeFloatInUnits(  pass_pol_angle,"deg");
  file_.writeFloatInUnits(  pass_emission_angle,"deg");
  file_.writeFloatInUnits(  pass_azimuth_angle,"deg");

 

  //-------------------
  //convert west bound longitude
  //from 0 to 180: from 360 to 180
  //from 0 to -180: from 0 to 180
  //---------------------
  //reset longitude  container
  east_west_.resize(5);
  east_west_[0]= pass_centroid_lon ;
  east_west_[1]= pass_ellipse_pt1_lon ;
  east_west_[2]= pass_ellipse_pt2_lon ;
  east_west_[3]= pass_ellipse_pt3_lon ;
  east_west_[4]= pass_ellipse_pt4_lon ;
  for(unsigned int i=0;i<5;++i){
    if(east_west_[i]<=zero_angle_ && east_west_[i] >= -pi_angle_) 
      east_west_[i]=-east_west_[i];
    else if(east_west_[i] >= zero_angle_ && east_west_[i] <= pi_angle_)
      east_west_[i]= two_pi_angle_ - east_west_[i];
    else{
      BurstDataError e("passive lon angle is out of range ");
      e.throwMe();
    }
  }

  file_.writeFloatInUnits(  east_west_[0],"deg");
  file_.writeFloatInUnits(  pass_centroid_lat,"deg");

  pass_major_width.writeFloat(file_,0,"km");
  pass_minor_width.writeFloat(file_,0,"km");

  //east to west lon #2-5
  file_.writeFloatInUnits( east_west_[1],"deg");
  file_.writeFloatInUnits( east_west_[2],"deg");
  file_.writeFloatInUnits( east_west_[3],"deg");
  file_.writeFloatInUnits( east_west_[4],"deg");
 
  
  file_.writeFloatInUnits( pass_ellipse_pt1_lat,"deg");
  file_.writeFloatInUnits( pass_ellipse_pt2_lat,"deg");
  file_.writeFloatInUnits( pass_ellipse_pt3_lat,"deg");
  file_.writeFloatInUnits( pass_ellipse_pt4_lat,"deg");

  file_.write(num_pulses_received);
  
  file_.writeFloatInUnits(total_echo_energy,"J");
  file_.writeFloatInUnits(noise_echo_energy,"J");
  file_.writeFloatInUnits(x_factor,"W");

  sigma0_uncorrected.writeFloat(file_,0,"");
  sigma0_corrected.writeFloat(file_,0,"");
  sigma0_uncorrected_std.writeFloat(file_,0,"");

  altitude_means.writeFloat(file_,0,"km");
  altitude_means_std.writeFloat(file_,0,"km");
  
  act_geom_time_offset.writeFloat(file_,0,"s");

  

  //-------------------
  //convert west bound longitude
  //from 0 to 180: from 360 to 180
  //from 0 to -180: from 0 to 180
  //---------------------
  east_west_.resize(5);//reset
  east_west_[0]= act_centroid_lon;
  east_west_[1]= act_ellipse_pt1_lon ;
  east_west_[2]= act_ellipse_pt2_lon;
  east_west_[3]= act_ellipse_pt3_lon;
  east_west_[4]= act_ellipse_pt4_lon ;
  for(unsigned int i=0;i<5;++i){
    if(east_west_[i]<=zero_angle_ && east_west_[i] >= -pi_angle_) 
      east_west_[i]=-east_west_[i];
    else if(east_west_[i] >= zero_angle_ && east_west_[i] <= pi_angle_)
      east_west_[i]= two_pi_angle_ - east_west_[i];
    else{
      BurstDataError e("passive lon angle is out of range ");
      e.throwMe();
    }
  }


  file_.writeFloatInUnits( act_pol_angle,"deg");
  file_.writeFloatInUnits(act_incidence_angle,"deg");
  file_.writeFloatInUnits(act_azimuth_angle,"deg");
 
  file_.writeFloatInUnits(east_west_[0],"deg");
  file_.writeFloatInUnits(act_centroid_lat,"deg");
  
  act_major_width.writeFloat(file_,0,"km");
  act_minor_width.writeFloat(file_,0,"km");

  //east to west lon #7-10
  file_.writeFloatInUnits(east_west_[1],"deg");
  file_.writeFloatInUnits(east_west_[2],"deg");
  file_.writeFloatInUnits(east_west_[3],"deg");
  file_.writeFloatInUnits(east_west_[4],"deg");
 
  
  file_.writeFloatInUnits(act_ellipse_pt1_lat,"deg"); 
  file_.writeFloatInUnits(act_ellipse_pt2_lat,"deg");
  file_.writeFloatInUnits(act_ellipse_pt3_lat,"deg");
  file_.writeFloatInUnits(act_ellipse_pt4_lat,"deg");
  
  altimeter_profile_range_start.writeFloat(file_,0,"km");
  altimeter_profile_range_step.writeFloat(file_,0,"km");
  file_.write(altimeter_profile_length);

  sar_azimuth_res.writeFloat(file_,0,"km");
  sar_range_res.writeFloat(file_,0,"km");

 
  file_.writeFloatInUnits( sar_centroid_bidr_lon,"deg");
  file_.writeFloatInUnits( sar_centroid_bidr_lat,"deg");
}



// This simple routine merely selects the calibrated attenuation value
// for the beam number. The nominal to calibrated mapping must have already
// been performed (i.e., during SAB::decode by preprocessor or 
// for the PointTargetSim in getRadarParams or getRadarParamsFromIeb
float BurstData::computeCalibratedAttenuatorGain()
{
  Uvar retval;
  switch(beam_number){
  case 1:
  case 2:
    retval=at1;
    break;
  case 3:
    retval=at3;
    break;
  case 4:
  case 5:
    retval=at4;
    break;
  default:
    ErrorMessage e("BurstData:::computeCalibratedAttenuatorGain bad beam_num");
    e.throwMe();
  }
 
  // returned value is between 0 and 1
  return(retval.getInUnits(""));
}

//-----------------------------------------------------
//
//
// this function is used to print read-in data on the screen
//
//
//-----------------------------------------------------

void BurstData::showPassiveDataOnScreen() 
{
 
 
 
  string input; 

  cout<<"sync"<< sync<<endl;
  cout<<"sclk"<<sclk<<endl; 
  cout<<"scpr"<<scpr<<endl;
  cout<<"brst:"<<brst <<endl;  
  cout <<"header tfi tnc typ: "<< header_tfi<< " "<<header_tnc << " "<< header_typ<<endl;
  cout <<"header tca tcb tcc:"<<header_tca<<" "<<header_tcb <<" "<<header_tcc <<endl;
  cout<<"pwri vicc vimc tail_id sab_counter: "<<pwri<<" "<<vicc<<" "<<vimc<<" "<<tail_id<<" "<<sab_counter<<endl;
  
  cout <<"fswm ctbc ctrx ctps ctbe header_end:"<<fswm<<" "<<ctbc<<" "<<ctrx<<" "<<ctps<<" "<<ctbe<<" "<<header_end<<endl;

getline(cin, input); if (input=="q"){return;}
 
  //-----------------------------------------
  // Extract data fields from uniform arrays.
  //-----------------------------------------

  //slow field
  //for variable names, see blue book  12-2
  cout <<"slow tfi dtn slow type csr r_mode sin: "<<" "<<
    slow_tfi <<" "<<dtn <<" "<<slow_typ <<" "<<csr <<" "<<
    r_mode <<" "<<slow_instruction_number <<endl;
  cout << "bem baq_mode tro rc_bw adc: "<<bem <<" "<<baq_mode <<" "<<tro <<" "<<
    rc_bw<< " "<< adc <<endl;
  cout <<"rip: "<<rip<<endl;
  cout << "at1db  at3_db  at4_db csd rad: "<<" "
    <<at1_db<<" "<<at3_db<<" "<<at4_db<<" "<<csd
       <<" "<<rad<<endl;
 
  cout <<"csq chirp length csf:"<<csq<<" "<<chirp_length<<" "<<slow_cfs<<endl;
  getline(cin, input); if (input == "q"){return;}

  
  //
  //fast field
  //
  cout << "fast tfi fin fast type pul bii bpd rwd pri chirp start frequency:"<<fast_tfi<<" "<<
    fin<<" "<<fast_typ<<" "<<pul<<" "<<bii<<" "<<bpd<<" "<<rwd<<" "<<pri<<" "<<fast_csf
    <<endl;


    
  
  //decoding footer1 depending on radar mode
  //r_mode <4 or r_mode>7 : ALT or Image mode --> BAQ
  //r_mode= 5,6,or 7: Inter-Galatic Object Calibration --> ENG Temp  
 
    cout <<"eng temp"<<endl;
    cout<<  fwdtmp <<endl;
   cout<< be3tmp <<endl;
   cout<< diptmp <<endl;
  cout<<   rlotmp<<endl;
   cout<< nsdtmp <<endl;
  cout<<  lnatmp <<endl;
   cout<< mratmp<<endl;
   cout<< mruttm<<endl;
  cout<<  wgb3t1<<endl;
  cout<<  wgb3t2<<endl;
   cout<< wgb3t3<<endl;
   cout<< nsdcur<<endl;
 
 

  cout <<"cntrl crnradio cntnd eout: "<<cnt_rl<<" "<<cnt_radio<<" "<<cnt_nd<<
    " "<<eout<<endl;

 
  getline(cin, input); if (input=="q"){return;}
  
  if (subr == 15)
  {  // hip and cip invalid for subr == 15, use previous values instead.
    cout <<" hip cip iebtth iebttl bgcalls devlmn devlda devlyr"<< hip<<endl<<
    cip <<endl<<
    iebtth<<endl<< 
    iebttl<<endl<<
    bgcalls<<endl<<
    delvmn<<endl<<
    delvda<<endl<<
    delvyr<<endl;
  }
  else
  {
    cout<< "hip and cip hol"<<endl<<
    hip<<endl<<
    cip<<endl;     
    
  }  
  cout <<"space craft time"<<space_craft_time<<endl;
  cout <<"sub comm data out"<<subr<<endl;

 if (subr == 0)
 {cout <<
   fwdtmp <<endl<<
   be1tmp <<endl<<
   be2tmp<<endl<<
   be3tmp <<endl<<
   be4tmp<<endl<<
   be5tmp<<endl
; 
 }
 else if (subr == 1)
 {   cout <<
   diptmp<<endl<<
   rlotmp <<endl<<
   tadcal1<<endl<<
   nsdtmp <<endl<<
   lnatmp<<endl<<
   evdtmp<<endl;
 }
 else if (subr == 2)
   { cout <<
       mratmp<<endl<<
   mruttm<<endl<<
   dcgttm<<endl<<
   cucttm<<endl<<
   twttmp<<endl<<
       epctmp<<endl;
 }
 else if (subr == 3)
   { cout <<
       tw1ttm <<endl<<
   ep1ttm<<endl<<
   p_stmp<<endl<<
   p_sttm<<endl<<
   fguttm<<endl<<
       tadcal4<<endl;
 }
 else if (subr == 4)
   {cout <<
   esstmp<<endl<<
   wgb1t1<<endl<<
   wgb3t1<<endl<<
   wgb3t2<<endl<<
   wgb3t3<<endl<<
      wgb5t1<<endl;
 }
 else if (subr == 5)
   {cout<<
   pcutmp<<endl<<
   adctmp<<endl<<
   tadcal2<<endl<<
   ecltmp<<endl<<
   cputmp<<endl<<
      memtmp<<endl;
 
 }
 else if (subr == 6)
   {cout <<
      sadctmp <<endl;
 }
 else if (subr == 7)
   {cout<<
      tadcal3<<endl<<
    frwdpw <<endl<<
    dcgmon<<endl<<
      lpltlm_db<<endl;
 }
 else if (subr == 8)
   {     cout <<
	   nsdcur<< endl<<
   hpapsm<< endl<<
   catcur<< endl<<
   p_smon << endl<<
   svlsta<< endl<<
	   usotmp<< endl;
 }
 else if (subr == 9)
   {cout<<
      cpbnkv<<endl<<
   essvlt<<endl<<
   tadcal5<<endl<<
   pcu5v_72 <<endl<<
   pcu5i_73<<endl<<
      pcu5v_74<<endl;
 }
 else if (subr == 10)
   {cout<<
      pcuti_75<<endl<<
   pcu15v_76<<endl<<
   pcu15i_77<<endl<<
   pcu15v_78<<endl<<
   pcu15i_79<<endl<<
      pcu12v_80<<endl;
 }
 else if (subr == 11)
   {cout<<
      pcu12i_81<<endl<<
   pcucur<<endl<<
   pllmon<<endl<<
   ctu5i <<endl<<
   tadcal6<<endl<<
      pcu9v_86<<endl;
 }
 else if (subr == 12)
   {cout<<
      pcu9i_87 <<endl<<
      pcu9v_88 <<endl<<
      pcu9i_89 <<endl;
 }
 else if (subr == 13)
   { cout<<
       tadcl7 <<endl;
 }
 else if (subr == 14)
 { cout <<
       shpttm<<endl;
 }
 else
   {cout <<"subr is 15"<<endl;
   } 

}



