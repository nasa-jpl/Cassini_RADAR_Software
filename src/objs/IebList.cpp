//----------------------------------------------------------------------------
//
// IebList.cpp
//
// This main fuction of this class is to generates idealized Ieb sequence 
// for a given time span and to provide the file I/O of Ieb file.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------


//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>


//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <fstream>
#include "Error.h"
#include "Utils.h"
#include "Units.h"
#include "Ieb.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "Ieb.h"
#include "IebList.h"
#include "Ambiguity.h"
#include "Plot.h"
#include "IebProfile.h"
#include "Flyby.h"
using std::cout;
using std::endl;
using std::cerr;

//----------------------------------------
// Constructors: contruct with cfg only
//-------------------------------------
IebList::IebList(Config& cfg)
  :allbeams_set_(false),
   set_csf_times_pri_integer_(false),
   cfg_(cfg),
   azimuth_3dB_("azimuth 3dB at 0 elevation"),
   elevation_3dB_("elevation 3dB at 0 azimuth"),
   azimuth_5dB_("azimuth 5dB at 0 elevation"),
   elevation_5dB_("elevation 5dB at 0 azimuth"),
   range_3dB_bore_("range 3db bore"),
   doppler_3dB_bore_("doppler 3db bore"),
   range_5dB_bore_("range 5db bore"),
   doppler_5dB_bore_("doppler 5db bore")
  {
  slowfield_.resize(Nslowfield);
  fastfield_.resize(Nfastfield);
  powerfield_.resize(Npowerfield);
  tncfield_.resize(Ntncfield);
  ieb_list_.clear();
  Bvector_.clear();
  Fvector_.clear();
  azimuth_3dB_.resize(5,2);
  elevation_3dB_.resize(5,2);
  azimuth_5dB_.resize(5,2);
  elevation_5dB_.resize(5,2);
  range_3dB_bore_.resize(5,5);
  doppler_3dB_bore_.resize(5,5);
  range_5dB_bore_.resize(5,5);
  doppler_5dB_bore_.resize(5,5);
  config();//read parameters for constructing beams
  }

//----------------------------------------
//Construction of IebList with time span and cfg
// : make a list of Ieb and altitude 
// both start and end time are epoch relative, which
// means I need to get an access to flyby object
//-----------------------------------------
IebList::IebList(const Uvar& start_time, 
		 const Uvar& end_time,
		 Config& cfg)
  :allbeams_set_(false),
   set_csf_times_pri_integer_(false),
   cfg_(cfg),
   azimuth_3dB_("azimuth 3dB at 0 elevation"),
   elevation_3dB_("elevation 3dB at 0 azimuth"),
   azimuth_5dB_("azimuth 5dB at 0 elevation"),
   elevation_5dB_("elevation 5dB at 0 azimuth"),
   range_3dB_bore_("range 3db bore"),
   doppler_3dB_bore_("doppler 3db bore"),
   range_5dB_bore_("range 5db bore"),
   doppler_5dB_bore_("doppler 5db bore")  
  {
  ieb_list_.clear();
  start_time_ = start_time;
  end_time_ = end_time;
  slowfield_.resize(Nslowfield);
  fastfield_.resize(Nfastfield);  
  powerfield_.resize(Npowerfield);
  tncfield_.resize(Ntncfield);
  Bvector_.clear();
  Fvector_.clear();
  azimuth_3dB_.resize(5,2);
  elevation_3dB_.resize(5,2);
  azimuth_5dB_.resize(5,2);
  elevation_5dB_.resize(5,2);
  range_3dB_bore_.resize(5,5);
  doppler_3dB_bore_.resize(5,5);
  range_5dB_bore_.resize(5,5);
  doppler_5dB_bore_.resize(5,5);
  Flyby flyby(cfg_);
  epoch_time_ = flyby.epochTime();
  config();//read parameters for constructing beams
  }

//-----------------------
//destructor
//----------------------
IebList::~IebList()
  {
  }


//------------------------
//Config setup
//------------------------
void IebList::config()
  { 

    //set up default frame
    target_name_=cfg_.str("target");
    if(target_name_=="none" || target_name_=="NONE" || target_name_=="None" ||
       target_name_=="source" || target_name_=="Source" || target_name_=="SOURCE")
      {
	target_name_="none";
      }
    else
      {
	target_frame_=Frame("IAU_"+target_name_,target_name_);
      }
    //usable area conditions  
    amb_ratio_=cfg_["signal_to_amb_ratiodB"].getInUnits("");
    noise_equiv_sigma0_=cfg_["noise_equivalent_sigma0dB"].getInUnits("");
    min_beam_gain_=cfg_["min_oneway_gaindB_wrt_peak"].getInUnits("");
    
    //cout<<"Beam objects are being constructed "<<endl;
    string beam_pattern_source;
    beam_pattern_source=cfg_.str("beam_pattern_source");
    bool beam_pattern = false;
    if (beam_pattern_source=="file") 
      {
	beam_pattern = true;
	//cout<<"beam files are being used"<<endl;
      }
    else if (beam_pattern_source=="sinc_model") 
      {
	beam_pattern = false;
	//cout<<"sinc pattern will be used "<<endl;
      }
    else    ErrorMessage("No beam pattern is set").throwMe();
    
    //make all 5 beams and store
    for (unsigned int i = 0; i < 5;++i)
      {
	Beam beam(i+1,cfg_,beam_pattern);
	Bvector_.push_back(beam); 
	Frame frame("CASSINI_RADAR_"+toStr(i+1),"Cassini");
	Fvector_.push_back(frame);
      }
    allbeams_set_ = true;
    computeBeam3dB_5dBpoints();
    
    //----------------------------
    //Radar System parameter:System noise temperature, Transmit power, 
    // and attenuation calibration factor
    //----------------------------
    data_take_number_ = (unsigned int) cfg_["data_take_number"].getInUnits("");
    Tsys_ = cfg_["Tsys"];
    Pt_ = cfg_["Pt"];
    
    Pn_cal_altl_ = boltzmann_constant * Tsys_*altl_rcv_bandwidth;
    radar_buffer_std_deviation_altl_=cfg_["cal_radar_buffer_std_deviation_altl"].getInUnits("");
    cal_gain_altl_=cfg_["cal_gaindB_altl"].getInUnits("");
    cal_gain_altl_ = pow(10.0,-cal_gain_altl_/10.0);
    squared_deviation_of_system_noise_cal_altl_=cfg_["squared_deviation_of_system_noise_input_at_ALTL"];
    squared_deviation_of_system_noise_cal_alth_=cfg_["squared_deviation_of_system_noise_input_at_ALTH"];
    squared_deviation_of_system_noise_cal_sarl_=cfg_["squared_deviation_of_system_noise_input_at_SARL"];
    squared_deviation_of_system_noise_cal_sarh_=cfg_["squared_deviation_of_system_noise_input_at_SARH"];
    
    
    pbw_ratio_=cfg_["pulse_bandwidth_ratio"].getInUnits("");
    //Muhleman_backscattering coeff: Need to calculate Radar Echo Power
    Uvar keyvalue_reading;
    keyvalue_reading = cfg_["Muhleman_backscatt_model_k1"];
    muhleman_k1_ = keyvalue_reading.getInUnits("");
    keyvalue_reading = cfg_["Muhleman_backscatt_model_k2"];
    muhleman_k2_ = keyvalue_reading.getInUnits("");

    //----------------------
    //csf selection for tone data
    //------------------------
    string csf_selection=cfg_.str("Make_csf_integer_multiple_of_prf");
    if(csf_selection=="Yes" || csf_selection=="yes"){
      set_csf_times_pri_integer_=true;
    }
    else set_csf_times_pri_integer_=false;
    


    if(set_csf_times_pri_integer_){
      min_frequency_of_returned_cw_echo_=cfg_["Min_frequency_of_returned_CW_echo"];
      max_frequency_of_returned_cw_echo_=cfg_["Max_frequency_of_returned_CW_echo"];
      if(min_frequency_of_returned_cw_echo_ > max_frequency_of_returned_cw_echo_) ErrorMessage("IebList.cpp:: cfg(): mininum return frequency of CW echo is larger than max return frequency of CW echo").throwMe();

      if(min_frequency_of_returned_cw_echo_ >Uvar(0,"Hz") ||
	 max_frequency_of_returned_cw_echo_ >Uvar(0,"Hz") ||
	 min_frequency_of_returned_cw_echo_ <-altl_adc/2) ErrorMessage("IebList.cpp:: cfg():return echo of CW is out of band").throwMe();

    }
  }
//-----------------------------------------------------------------------------
//Compute Beam3dB and 5dB  points: along azimuth (elev=0) and elevation(azi=0)
//-----------------------------------------------------------------------------
void IebList::computeBeam3dB_5dBpoints()
  {
    //find 3dB and 5dB  points along azimuth and elevation axis
    if(!allbeams_set_) ErrorMessage("IebList::computeBeam3dBpoints: Not all 5 beams are set up").throwMe();
    
    //cout<<"Find Beam 3dB points along azimuth and elevation axes"<<endl;
    for (unsigned int i = 0; i < 5;++i)
      {
	Uvar zero_degree("zero degree",0,"rad");
	double target_gain = -3.0;//-3.0 dB
	Uvec azimuth_tmp("azimuth 3db ");
	Uvec elevation_tmp("elevation 3db");
	
	Uvar min,max;
	
	elevation_tmp = Bvector_[i].azimAndGainToElev(zero_degree,target_gain);
	min_max(min,max,elevation_tmp);
	elevation_3dB_(i,0) = min;
	elevation_3dB_(i,1) = max;
    
	//cout<<"at 0 azimuth "<< min.getInUnits("deg")<<" "
	//  <<max.getInUnits("deg")<<endl;
    
	azimuth_tmp=Bvector_[i].elevAndGainToAzim(zero_degree,target_gain);
	min_max(min,max,azimuth_tmp);
	azimuth_3dB_(i,0) = min;
	azimuth_3dB_(i,1) = max;
	//cout<<"at 0 elevation "<< min.getInUnits("deg")<<" "
	//  <<max.getInUnits("deg")<<endl;
      }

    //cout<<"Find Beam 5dB points along azimuth and elevation axes"<<endl;
    for (unsigned int i = 0; i < 5;++i)
      {
	Uvar zero_degree("zero degree",0,"rad");
	double target_gain = -5.0;//-5.0 dB
	Uvec azimuth_tmp("azimuth 5db ");
	Uvec elevation_tmp("elevation 5db");
	
	Uvar min,max;
	elevation_tmp = Bvector_[i].azimAndGainToElev(zero_degree,target_gain);
	min_max(min,max,elevation_tmp);
	elevation_5dB_(i,0) = min;
	elevation_5dB_(i,1) = max;
    
	//cout<<"at 0 azimuth "<< min.getInUnits("deg")<<" "
	//<<max.getInUnits("deg")<<endl;
    
	azimuth_tmp=Bvector_[i].elevAndGainToAzim(zero_degree,target_gain);
	min_max(min,max,azimuth_tmp);
	azimuth_5dB_(i,0) = min;
	azimuth_5dB_(i,1) = max;
	//cout<<"at 0 elevation "<< min.getInUnits("deg")<<" "
	//<<max.getInUnits("deg")<<endl;
      }
  }

//--------------------------------------------------------------
//Compute Range and Doppler values on 3dB and 5dB points
// For each beam, 5 points(2 azimuth, 2 elevation's 3dB and 5dB, and
// plus boresight) range and doppler values will be  computed
// range_5dB_bore(beam_number, 0) at (azi,elev) = (-5dB point, 0)
// range_5dB_bore(beam_number, 1) at (azi,elev) = (5dB point, 0)
// range_5dB_bore(beam_number, 2) at (azi,elev) = (0, -5dB point)
// range_5dB_bore(beam_number, 3) at (azi,elev) = (0, 5dB point)
// range_5dB_bore(beam_number,4) at boresight
//----------------------------------------------------------------
void IebList::computeRangeDoppler_3dB_5dB_Borepoints(const Time& t)
  {
    Uvar zero_degree = Uvar(0,"rad");
    StateVector sc_state;
    target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
    for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
      {
      DirectionVector look("look",Fvector_[i_beam],t,0,0,1);
      TargetGeom tg_bore(t);
      tg_bore.setState(sc_state);
      tg_bore.setTarget(target_name_);
      tg_bore.setLookDirection(look);
      range_3dB_bore_(i_beam,4) = tg_bore.range();
      doppler_3dB_bore_(i_beam,4) = tg_bore.doppler(speed_light/carrier_frequency);
      range_5dB_bore_(i_beam,4) = range_3dB_bore_(i_beam,4);
      doppler_5dB_bore_(i_beam,4) = doppler_3dB_bore_(i_beam,4);
      
      TargetGeom tg(t);
      for (unsigned int i = 0; i < 2;++i)
	{
	look.setAzimuthElevation(azimuth_3dB_(i_beam,i),zero_degree);
	tg.setState(sc_state);
	tg.setTarget(target_name_);
	tg.setLookDirection(look);
	range_3dB_bore_(i_beam,i) =tg.range();
	doppler_3dB_bore_(i_beam,i)=tg.doppler(speed_light/carrier_frequency);
	//cout<<"doppler 3db"<<doppler_3dB_bore_(i_beam,i)<<endl;
	tg.reset(t);
	}
      for (unsigned int i = 0; i < 2;++i)
	{
	look.setAzimuthElevation(zero_degree,elevation_3dB_(i_beam,i));
	tg.setState(sc_state);
	tg.setTarget(target_name_);
	tg.setLookDirection(look);
	range_3dB_bore_(i_beam,2+i) = tg.range();
	doppler_3dB_bore_(i_beam,2+i)=tg.doppler(speed_light/carrier_frequency);
	tg.reset(t);
	}
      
      //5db points
      for (unsigned int i = 0; i < 2;++i)
	{
	look.setAzimuthElevation(azimuth_5dB_(i_beam,i),zero_degree);
	tg.setState(sc_state);
	tg.setTarget(target_name_);
	tg.setLookDirection(look);
	range_5dB_bore_(i_beam,i) =tg.range();
	doppler_5dB_bore_(i_beam,i)=tg.doppler(speed_light/carrier_frequency);
	tg.reset(t);
	}
      for (unsigned int i = 0; i < 2;++i)
	{
	look.setAzimuthElevation(zero_degree,elevation_5dB_(i_beam,i));
	tg.setState(sc_state);
	tg.setTarget(target_name_);
	tg.setLookDirection(look);
	range_5dB_bore_(i_beam,2+i) = tg.range();
	doppler_5dB_bore_(i_beam,2+i)=tg.doppler(speed_light/carrier_frequency);
	//cout<<"doppler 5db "<<doppler_5dB_bore_(i_beam,i)<<endl;
	tg.reset(t);
	}
      }//loop over 5 beams
  }
//----------------------------------------
// {get} two azimuth points that produce -3dB oneway gain 
// when elevation =0
//-----------------------------------------
void IebList::getBeam3dB_azimuth_points_at_0_elevation(Umat& azimuth_3dB)
  {
    if(!allbeams_set_) ErrorMessage("Beams are not set up in config method").throwMe();
    azimuth_3dB.resize(azimuth_3dB_.rows()
		       ,azimuth_3dB_.cols());
    azimuth_3dB = azimuth_3dB_;
  }
//----------------------------------------
// {get} two elevation points that produce -3dB oneway gain 
// when azimuth =0
//-----------------------------------------
void IebList::getBeam3dB_elevation_points_at_0_azimuth(Umat& elevation_3dB)
  {
    if(!allbeams_set_) ErrorMessage("Beams are not set up in config method").throwMe();
    elevation_3dB.resize(elevation_3dB_.rows()
			 ,elevation_3dB_.cols()); 
    elevation_3dB = elevation_3dB_;    
  }

//----------------------------------------
// Edited by Chandini
// {get} two elevation points that produce -3dB oneway gain 
// when azimuth =0
//-----------------------------------------
Umat IebList::getBeam3dB_azimuth()
  {  
    return(azimuth_3dB_);   
  }

//----------------------------
//Return the time of the first record
//-----------------------------
Time IebList::getStartTime()
  {
    if (ieb_list_.size() == 0) ErrorMessage("IebList::getStartTime: No record avail").throwMe();
    return( (*ieb_list_.begin()).getTime());
  }

//-----------------------------------
//return the time of the last record
//-----------------------------------
Time IebList:: getEndTime()
  {
  if (ieb_list_.size() == 0) ErrorMessage("IebList::getEndTime: No record avail").throwMe();
   return( (*(ieb_list_.end()-1)).getTime());
  }

//------------------------------
//get first valid slow fast field time
//------------------------------
Time IebList::getFirstValidSlowFastFieldTime()
  {
    if(ieb_list_.size()== 0) ErrorMessage("ieb list container has 0 size").throwMe();
    Ieb_ptr i_ieb=ieb_list_.begin();
    bool slowfield_set = false;
    bool fastfield_set = false;
    Time t;
    for(i_ieb=ieb_list_.begin();i_ieb != ieb_list_.end();++i_ieb)
      {
	if(i_ieb->getIebMod() >1 ) fastfield_set= true;
	if(i_ieb->getIebMod() > 2) slowfield_set = true;
	if(slowfield_set && fastfield_set)
	  {
	    t=i_ieb->getTime();
	    break;
	  }
      }
    if(i_ieb==ieb_list_.end()) ErrorMessage("No ieb record has valid slow and fast field:: your ieblist is not valid for running radar").throwMe();
    return(t);
  }
//------------------------
//return record length
//----------------------
unsigned int IebList::getNumberofIebRecords()
  {
    return(ieb_list_.size());
  }

//--------------------------------
//Index
//--------------------------------
Ieb& IebList::operator[] (const unsigned int& i)
{//linear search
    if(ieb_list_.size()==0) ErrorMessage("IebList.cpp::operator[] zero ieb record length").throwMe();
    if(i > ieb_list_.size()) ErrorMessage("IebList.cpp::operator[] requested index is larger than the size of ieb record").throwMe();
    return(ieb_list_[i]);
  }
//--------------------------
//{get} ieb at time t
//---------------------------
Ieb IebList::getIeb_v0(const Time& t)
  {    
  bool slowfield_set= false;
  bool fastfield_set= false;
  bool powerfield_set=false;
  bool tncfield_set=false;
  


  Ieb ieb(t);
  if(ieb_list_.size()==0 || t < ieb_list_.begin()->getTime())
    ErrorMessage("no record avail ").throwMe();
  else if(ieb_list_.size()==1 ){
    if(t < ieb_list_.begin()->getTime())  ErrorMessage("requested time is earlier than first record ").throwMe();
    else return(ieb_list_[0]);
  }
  else{
    ///---------------------------------------------------
    //if the request time is between the first and last, 
    //find it
    //
    //---------------------------------------------------- 
    Ieb_ptr  ieb_reference= ieb_list_.begin();
    for ( Ieb_ptr ieb_itr = ieb_list_.begin();ieb_itr < ieb_list_.end()-1; ++ieb_itr)
      {
	//linear search
	if ( t >= ieb_itr->getTime() && t < (ieb_itr+1)->getTime())
	  {
	    ieb_reference = ieb_itr;
	    break;
	  }  
      }
    
    //special case when t >= ieb_list_.end()-1: not taken care of yet
    if(t >= (ieb_list_.end()-1)->getTime()) ieb_reference = ieb_list_.end()-1;
    if(ieb_reference == ieb_list_.end()) ErrorMessage("IebList.cpp::request time is between beginning and end of ieblist , but could not find it:: Bug!!").throwMe();
    
    //cout<<"IEBLIST: getIeb:request time and list time "
    //<< t.utc("ISOD")<<" "<<ieb_reference->getTime().utc("ISOD")<<":MOD:"
    //<<ieb_reference->getIebMod() << endl;
    
    
    Ieb_ptr ieb_slow=ieb_reference;
    Ieb_ptr ieb_fast = ieb_reference;
    Ieb_ptr ieb_power=ieb_reference;
    Ieb_ptr ieb_tnc=ieb_reference;
    
    //find  valid slow field    
    for (ieb_slow=ieb_reference;ieb_slow !=ieb_list_.begin();ieb_slow--)
      {
	if(ieb_slow->getIebMod() > 2)
	  { 
	    slowfield_= ieb_slow->getSlowfield(); 
	    slowfield_set = true;
	    break;
	  }
      }
    if(ieb_slow == ieb_list_.begin())
      {
	if( ieb_slow->getIebMod()> 2)
	  { 
	    slowfield_= ieb_slow->getSlowfield(); 
	    slowfield_set = true;
	  }
	else
	  slowfield_set= false;
      }
    
    //find valid fast field
    for (ieb_fast = ieb_reference;ieb_fast!=ieb_list_.begin();ieb_fast--)
      {
	if(ieb_fast->getIebMod()>1)
	  {
	    fastfield_= ieb_fast->getFastfield(); 
	    fastfield_set = true;
	    break;
	  }
      }     
    if(ieb_fast == ieb_list_.begin())
      {
	if(ieb_fast->getIebMod()>1)
	  {
	    fastfield_= ieb_fast->getFastfield(); 
	    fastfield_set = true;
	  }
	else
	  fastfield_set = false;
      }
    
    //find valid power
    for(ieb_power=ieb_reference;ieb_power != ieb_list_.begin();--ieb_power)
      {
	if(ieb_power->getIebMod()==0 || ieb_power->getIebMod()==4)
	  {
	    powerfield_=ieb_power->getPowerfield();
	    powerfield_set = true;
	    break;
	  }
      }
    if(ieb_power==ieb_list_.begin())
      {
	if (ieb_power->getIebMod()==0 ||ieb_power->getIebMod()==4)
	  {
	    powerfield_=ieb_power->getPowerfield();
	    powerfield_set = true;
	  }
	else
	  powerfield_set = false;
      }
    
    //find valid tnc
    for(ieb_tnc=ieb_reference;ieb_tnc!=ieb_list_.begin();--ieb_tnc)
      {
	
	if(ieb_tnc->getIebMod()==1 || ieb_tnc->getIebMod()==4)
	  {
	    tncfield_ = ieb_tnc->getTncfield();
	    tncfield_set = true;
	    break;
	  }
      }
    if(ieb_tnc==ieb_list_.begin())
      {
	if(ieb_tnc->getIebMod()==1 ||ieb_tnc->getIebMod()==4)
	  {
	    tncfield_ = ieb_tnc->getTncfield();
	    tncfield_set = true;
	  }
	else
	  tncfield_set= false;
      }
  }
  unsigned int ieb_mode=1;//default(test and calibration mode)
  if(powerfield_set)  ieb.loadPowerfield(powerfield_);
  if(tncfield_set)  ieb.loadTncfield(tncfield_);
  if(slowfield_set && fastfield_set) ieb.loadSlowFastfield(slowfield_,fastfield_);

  if(powerfield_set&& !tncfield_set && !slowfield_set && !fastfield_set) 
    ieb_mode= 0; 
    
  if(!powerfield_set && tncfield_set && !slowfield_set && !fastfield_set) 
    ieb_mode =1;
 
  if(!powerfield_set && !tncfield_set && !slowfield_set && fastfield_set) 
    ieb_mode = 2;
  
  if(slowfield_set && fastfield_set) 
    ieb_mode = 3;
  
  if(powerfield_set&& tncfield_set && slowfield_set && fastfield_set) 
    ieb_mode = 4;

  switch (ieb_mode){
  case (0):
    //cout << "IebGet Power " << endl; 
    ieb.decodePowerfield();
    break;
  case (1):
    //cout << "IebGet TNC " << endl; 
    ieb.decodeTncfield();
    break;
  case (2):
    //cout << "IebGet Fast " << endl; 
    ieb.decodeSlowFastfield(); // fast by itself???
    break;
  case (3):
    //cout << "IebGet Slow/Fast " << endl; 
    ieb.decodeSlowFastfield();
    break;
  case (4):
    //cout << "IebGet S/F/P/T (4) " << endl;
    ieb.decodeSlowFastfield(); 
    ieb.decodeTncfield(); 
    ieb.decodePowerfield();
    break;
  default:
      ErrorMessage("GetIeb::no iebmode found?").throwMe();

 }

  ieb.setIebMod(ieb_mode);
  return(ieb);
  }

//--------------------
//speed up method
//----------------------
Ieb IebList::getIeb(const Time& t)
{    //bisection search
  bool slowfield_set= false;
  bool fastfield_set= false;
  bool powerfield_set=false;
  bool tncfield_set=false;
  Ieb ieb(t);
  unsigned int Nmark;

  // case(0): no element or request time earlier than first element
  if(ieb_list_.size()==0 || t<ieb_list_.begin()->getTime()){
    ErrorMessage("no record at the requested time").throwMe(); 
  }
  // case (1) only one element
  else if(ieb_list_.size()==1){
    if(t<ieb_list_.begin()->getTime())
      ErrorMessage("requested time is earlier than the first record time").throwMe(); 
    else return(ieb_list_[0]);
  }
  else
    {
      //---------------------------------------------------
      //if the request time is between the first and last, 
      //find it
      //
      //---------------------------------------------------- 
     
      //bisection search
      unsigned int Nstart = 0;
      unsigned int Nend = ieb_list_.size()-1;
      unsigned int Nmiddle=(Nstart+Nend)/2;
      unsigned int Ndiff= Nend - Nstart;
      unsigned int Nbisection = (unsigned int) (log(double(Ndiff))/log(2)+1);
      if(t >= ieb_list_[Nmiddle].getTime()) Nstart= Nmiddle;
      else Nend = Nmiddle;
      for(unsigned int jj=0;jj<=Nbisection;++jj){
	if(t >= ieb_list_[Nstart].getTime() && t< ieb_list_[Nstart+1].getTime() ) 
	  {
	    Nmark=Nstart;
	    break;
	  }
	Ndiff = Nend-Nstart;
	//cout<<" ndiff "<<jj<<" "<<Nstart <<" "<<Nend<<" "<<Ndiff<<endl;
	//cout<<"t "<<t.utc("ISOD")<<endl;
	//cout<<ieb_list_[Nstart].getTime().utc("ISOD")<<endl;
	//cout<<ieb_list_[Nend].getTime().utc("ISOD")<<endl;
	
	if(Ndiff>1)
	  {//there is one ieb between
	    Nmiddle = (Nstart+Nend)/2;
	    if(  t >= ieb_list_[Nmiddle].getTime() )   Nstart = Nmiddle;
	    else Nend = Nmiddle;
	  }
	else 
	  {//there is 0 or 1 element
	    if(ieb_list_[Nstart].getTime() == t)
	      {
		Nmark = Nstart;
		break;
	      } 
	    else if(ieb_list_[Nend].getTime()==t)
	      {
		Nmark = Nend;
		break;
	      }
	    else if(ieb_list_[Nstart].getTime()>= t && ieb_list_[Nend].getTime()< t)
	      {
		Nmark = Nstart;
		break;
	      }
	    else ErrorMessage("ieblist:: could not locate1 t").throwMe();
	  }
      }

      //special case when t == ieb_list_.end()-1: not taken care of yet
      if(t >= (ieb_list_.end()-1)->getTime()) Nstart = ieb_list_.size()-1;
     

      //cout<<"IEBLIST: getIeb:request time and list time "
      //<< t.utc("ISOD")<<" "<<ieb_reference->getTime().utc("ISOD")<<":MOD:"
      //<<ieb_reference->getIebMod() << endl;
     
      
      
      unsigned int  ieb_slow = Nmark;
      unsigned int  ieb_fast = Nmark;
      unsigned int  ieb_power= Nmark;
      unsigned int  ieb_tnc= Nmark;

      //find  valid slow field    
      for (ieb_slow=Nmark;ieb_slow !=0;ieb_slow--)
	{
	  if(ieb_list_[ieb_slow].getIebMod() > 2)
	    { 
	      slowfield_= ieb_list_[ieb_slow].getSlowfield(); 
	      slowfield_set = true;
	      break;
	    }
	}
      if(ieb_slow == 0)
	{
	  if( ieb_list_[0].getIebMod()> 2)
	    { 
	      slowfield_= ieb_list_[0].getSlowfield(); 
	      slowfield_set = true;
	    }
	  else
	    slowfield_set= false;
	}

      //find valid fast field
      for (ieb_fast = Nmark;ieb_fast!=0;ieb_fast--)
	{
	  if(ieb_list_[ieb_fast].getIebMod()>1)
	    {
	      fastfield_= ieb_list_[ieb_fast].getFastfield(); 
	      fastfield_set = true;
	      break;
	    }
	}     
      if(ieb_fast == 0)
	{
	  if(ieb_list_[ieb_fast].getIebMod()>1)
	    {
	      fastfield_= ieb_list_[ieb_fast].getFastfield(); 
	      fastfield_set = true;
	    }
	  else
	    fastfield_set = false;
	}

      //find valid power
      for(ieb_power=Nmark;ieb_power != 0;--ieb_power)
	{
	  if(ieb_list_[ieb_power].getIebMod()==0 || ieb_list_[ieb_power].getIebMod()==4)
	    {
	      powerfield_=ieb_list_[ieb_power].getPowerfield();
	      powerfield_set = true;
	      break;
	    }
	}
      if(ieb_power==0)
	{
	  if (ieb_list_[ieb_power].getIebMod()==0 ||ieb_list_[ieb_power].getIebMod()==4)
	    {
	      powerfield_=ieb_list_[ieb_power].getPowerfield();
	      powerfield_set = true;
	    }
	  else
	    powerfield_set = false;
	}

      //find valid tnc
      for(ieb_tnc=Nmark;ieb_tnc!=0;--ieb_tnc)
	{
	 
	  if(ieb_list_[ieb_tnc].getIebMod()==1 || ieb_list_[ieb_tnc].getIebMod()==4)
	    {
	      tncfield_ =ieb_list_[ieb_tnc].getTncfield();
	      tncfield_set = true;
	      break;
	    }
	}
      if(ieb_tnc==0)
	{
	  if(ieb_list_[ieb_tnc].getIebMod()==1 ||ieb_list_[ieb_tnc].getIebMod()==4)
	    {
	      tncfield_ = ieb_list_[ieb_tnc].getTncfield();
	      tncfield_set = true;
	    }
	  else
	    tncfield_set= false;
	}
    }
  unsigned int ieb_mode=1;//default(test and calibration mode)
  if(powerfield_set)  ieb.loadPowerfield(powerfield_);
  if(tncfield_set)  ieb.loadTncfield(tncfield_);
  if(slowfield_set && fastfield_set) ieb.loadSlowFastfield(slowfield_,fastfield_);

  if(powerfield_set&& !tncfield_set && !slowfield_set && !fastfield_set) 
    ieb_mode= 0; 
    
  if(!powerfield_set && tncfield_set && !slowfield_set && !fastfield_set) 
    ieb_mode =1;
 
  if(!powerfield_set && !tncfield_set && !slowfield_set && fastfield_set) 
    ieb_mode = 2;
  
  if(slowfield_set && fastfield_set) 
    ieb_mode = 3;
  
  if(powerfield_set&& tncfield_set && slowfield_set && fastfield_set) 
    ieb_mode = 4;

  switch (ieb_mode){
  case (0):
    //cout << "IebGet Power " << endl; 
    ieb.decodePowerfield();
    break;
  case (1):
    //cout << "IebGet TNC " << endl; 
    ieb.decodeTncfield();
    break;
  case (2):
    //cout << "IebGet Fast " << endl; 
    ieb.decodeSlowFastfield(); // fast by itself???
    break;
  case (3):
    //cout << "IebGet Slow/Fast " << endl; 
    ieb.decodeSlowFastfield();
    break;
  case (4):
    //cout << "IebGet S/F/P/T (4) " << endl;
    ieb.decodeSlowFastfield(); 
    ieb.decodeTncfield(); 
    ieb.decodePowerfield();
    break;
  default:
      ErrorMessage("GetIeb::no iebmode found?").throwMe();

 }

  ieb.setIebMod(ieb_mode);
  return(ieb);
  }

//------------------------------------------
//get live updated ieb
//--------------------------------------------
Ieb IebList::getLiveIvpUpdateIeb(const Time& t, const Uvar& delta_trigger)
  {
    //-----------------------------------------
    //here is a key point
    //if delta trigger is positive, live update epoch occurs later time.
    // 
    //However, to simulate trigger update with old set of ieblist
    // one needs to ask ieb at delta_trigger earlier
    //----------------------------------------

    //if you are confused, let's say
    //  old epoch at 0 min and trigger time -5 min 
    //  new epoch at 1 min and delta trigger +1 min in order to
    // reset trigger time -4 min 
    //  so, after live update, 0 min ieb will be shifted to +1 min ieb
    // in order to have the same ieb setting at epoch
    // as a result, at epoch (+1 min), one ends up using ieb at 0 min
   // after live update
    //From the above statement, 
    // it is obvious to ask -delta_trigger earlier ieb from old
    // ieblist to have correct ieb at new epoch (+1 min)
    return(getIeb(t-delta_trigger));
  }
//-------------------------------------------------------
//get valid ieb at the requested time
// if time does not match, no ieb will be returned
//-----------------------------------------------------------
Ieb IebList::getIebOn(const Time& t)
  {
    Ieb ieb;
    if (t < ieb_list_.begin()->getTime()) ErrorMessage("IebList.cpp::getIebAt: requested time is earlier than first ieb time").throwMe();    
    if (t > (ieb_list_.end()-1)->getTime()) ErrorMessage("IebList.cpp::getIebAt: requested time is later than last ieb time").throwMe();    
   
    bool found_ieb= false;
    for (unsigned int i=0;i<ieb_list_.size();++i)
      {
	if(t==ieb_list_[i].getTime()) 
	  {
	    ieb=ieb_list_[i];
	    found_ieb=true;
	    break;
	  }
      }
    if(!found_ieb) ErrorMessage("IebList.cpp::getIebOn(): No IEB record at the requested time").throwMe();
    return(ieb);
  }
//------------------------------------
//{add} a single ieb to exiting ieblist
//------------------------------------
void IebList::addIeb(const Ieb& ieb)
  {   
   
  if (ieb_list_.size() == 0)
    {
      ieb_list_.push_back(ieb);
    }
  else if (ieb.getTime() < ieb_list_.begin()->getTime())
    {
    //when ieb.getTime is earlier than ieb_list time
      ieb_list_.insert(ieb_list_.begin(),ieb);//insert one before begin
    }
  else if (ieb.getTime() > (ieb_list_.end()-1)->getTime())
    { //when ieb.getTime is later than ieb_list time
    ieb_list_.push_back(ieb);
    }
  else if (ieb.getTime()> ieb_list_.begin()->getTime() && 
	   ieb.getTime() < (ieb_list_.end()-1)->getTime())
    {
      //when ieb.getTime is between record [0] and [N-1]
     
      for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
	   ieb_pt < (ieb_list_.end()-1);
	   ++ieb_pt)
	{
	  //check whether there is a record with the same time
	  //if it is so, throw an error
	  if (ieb.getTime() == ieb_pt->getTime()) 
	    ErrorMessage("There is an ieb record at the same time").throwMe();
	  if (ieb.getTime() == (ieb_pt+1)->getTime()) 
	    ErrorMessage("There is an ieb record at the same time").throwMe();
	  
	  if (ieb.getTime() > ieb_pt->getTime() && 
	      ieb.getTime() < (ieb_pt+1)->getTime())
	    {
	      //add a record before ieb_add   	
	      ieb_list_.insert(ieb_pt+1,ieb);
	      break;
	    }
	}
    }
  else
    {
    ErrorMessage("IebList::addIeb; ieb was not inserted").throwMe();
    }   
  }

//-----------------------------------------
//back_insert(const Ieb& ieb)
// add an IEB at the end of ieb list
// on the assumption that the inserted one is the latest one
// if not, it will be thrown out
//------------------------------------------
void IebList::back_insert(const Ieb& ieb, bool& valid )
  {
    if(ieb_list_.size()==0)
      {
	ieb_list_.push_back(ieb);
	valid = true;
      }
    else
      {
	if(ieb_list_.back().getTime() < ieb.getTime())
	  {
	    ieb_list_.push_back(ieb);
	    valid = true;
	  }
	else
	  {//in principle, this should throw error
	    // but, as of now, dump warning on the screen
	    valid = false;
	    Time t1,t2;
	    t1= ieb_list_.back().getTime();
	    t2=ieb.getTime();
	    cout<<"--------------warning----------------"<<endl;
	    cout<<"Back inserted ieb is not the latest one "<<endl;
	    cout<<"Current ieb will not be added into ieblist "<<endl;
	    cout<<"Last time record in ieblist "<< t1.utc("ISOD")<<endl;
	    cout<<"Latest Sab ieb time  "<<t2.utc("ISOD")<<endl;
	    cout<<"-------------------------------------"<<endl;
	    cout<<endl;
	  }
      }
  }

//-----------------------
//{delete} a single record of  ieb
//------------------------
void IebList::deleteIeb(const Ieb& ieb)
  {    
    if(ieb_list_.size()==0) ErrorMessage("IebList::delete: no record to delete").throwMe();
    //check there is a record at ieb.getTime();
    unsigned int record_found = 0;
   
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin();
	 ieb_pt < ieb_list_.end()  ;
	 ++ieb_pt)
      {
	if (ieb.getTime() == ieb_pt->getTime())
	  {
	    record_found = 1;
	    ieb_list_.erase(ieb_pt);
	    break;
	  }
      }
    if (record_found != 1) cout<<"IebList.cpp::there is no record to delete "<<endl;
  }

//-----------------------
//{delete} a multiple records of  ieb
// send it a good time and delete all IEB from that time forward until duration
// the good time should not be deleted
//------------------------
void IebList::deleteIebBlock(const Time& goodtime, const unsigned int& duration)
  {    
    //cout << "\tGH:deleteIEBBlock:" << goodtime.getCurrentLocalTime() <<"::SCLK)"<< goodtime.sclk("Cassini")<<" Dur:"<<duration  << endl;

    if(ieb_list_.size()==0) ErrorMessage("IebList::delete: no record to delete").throwMe();
    

    //check there is a record at ieb.goodtime();
    unsigned int record_found = 0;
    Time endtime = goodtime ;
    //Time endtime ( ("Cassini"), goodtime.sclk("Cassini") ) ;

    // Ieb ieb(goodtime);
    //int the_mode = ieb.getMod();
    unsigned int deletecount = 0;
    bool skip_next = false ;

    endtime = endtime + Uvar (duration, "s") ;
    //cout << "goodtime:" << endtime.sclk("Cassini")<<"DUR:"<<duration<< endl;

  
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin();
	 ieb_pt < ieb_list_.end()  ;
	 ++ieb_pt)
      {
	//cout << "Every IEB: " << ieb_pt->getTime() <<":"<< ieb_pt->getFastTfi() << endl;
	if (ieb_pt->getTime() > goodtime  && ieb_pt->getTime() <= endtime ) { //goodtime
	  cout <<"block HERE HERE HERE:" << ieb_pt->getTime().sclk("Cassini") <<":"<<ieb_pt->getFastTfi()<<":"<<endtime.sclk("Cassini") << "MOD:"<< ieb_pt->getIebMod() << endl;

	  if (ieb_pt->getIebMod() == 2 || ieb_pt->getIebMod() == 3 ) // delete slow and slow/fast NOT power or TNC
	    {

	      if (skip_next == false ) { //skip if
		//cout <<"DELETE HERE HERE HERE:" << ieb_pt->getTime().sclk("Cassini") <<":"<<ieb_pt->getFastTfi()<<":"<<endtime.sclk("Cassini") << endl << endl;
	      record_found = 1;
	      deletecount++;
	      ieb_list_.erase(ieb_pt--);  // is it ieb_list_.erase(ieb_pt); or (ieb_pt--)???
	      } else { // skip if (skip was true) 
		skip_next = false;
	      }
 
	    } else skip_next = true;  // This is to not delete right after a power or TNC

	} // end goodtime
	//	 if (ieb.getTime() >  ieb_pt->getTime()+Uvar (duration, "s") ) break;
      }
    if (record_found != 1) ErrorMessage("There is no record to delete");
    cout << "For time: "<< goodtime << "delete iebs = " << deletecount << endl;
  }

//-----------------------------------------------
//Generate Ideal Ieb Sequence
//-------------------------------------------------
void IebList::generateIdealIebSequence()
  {
    IebProfile iebprofile(cfg_);//set up iebprofile object
    Time trigger_time=iebprofile.getTriggerTime();
    Uvar first_div_start_time=iebprofile.getFirstDivisionStartTime();
    Uvar last_div_end_time=iebprofile.getLastDivisionEndTime();
    
    //cout<<"first time "<< first_div_start_time<<endl;
    //cout<<"last time"<<last_div_end_time<<endl;

    Uvar ND = iebprofile.getND();
    Uvar RL = iebprofile.getRL();
    Uvar Max_ND = iebprofile.getMaxND();
    Uvar Max_RL = iebprofile.getMaxRL();
    ieb_list_.clear();//clear ieb container
    
    Time t = epoch_time_;
    Uvar epoch_relative_time = start_time_;
   
    string div_name, radar_mode, autorad;
    unsigned int csr,baq,beam_mask,pulse, bursts_in_flight;
    double  noise_bit;
    Uvar bpd,prf,time_step,rip;
    double dutycycle,data_rate;
    unsigned int percentOfBW = 100;
    int tro=0;
    int ping_pong_counter =1;

    
    // cout<<"ping_pong_counter: "<<ping_pong_counter<<endl;
    while(epoch_relative_time < end_time_ &&
	  epoch_relative_time > first_div_start_time &&
	  epoch_relative_time < last_div_end_time)
      {
	//debug:: count time
	//int tick= int((epoch_relative_time - start_time_).getInUnits(""));
	//if(tick%10==0) cout<<" elapsed time "<< epoch_relative_time<<endl;
	//(0) Find division
	//When there is no division, skip and go next time
	while(!iebprofile.foundDivision(epoch_relative_time))
	  {
	    epoch_relative_time += Uvar(1,"s");
	  }
	div_name=iebprofile.getDivisionName(epoch_relative_time);
	radar_mode= iebprofile.getDivisionMode(div_name) ;
	cout<<"Radar mode: "<<radar_mode<<endl;
	cout<<"last_div_end_time: "<<last_div_end_time<<endl;
	//(1) set absolute time and find surface intercept points
	t = epoch_time_ + epoch_relative_time;
	if (radar_mode != "radiometer") computeRangeDoppler_3dB_5dB_Borepoints(t);//geom cal

	//(2) For scatterometer or scatterometer_compressed,
	//if there is no surface intercept point for beam3  boresight
	// and there is no previous ieb record,
	// skip and go to next ime    
      
	while((radar_mode=="scatterometer" || radar_mode=="scatterometer_compressed")
	      && ieb_list_.size()==0 
	      && range_3dB_bore_(2,4) == Uvar(0,"km"))
	  {
	    epoch_relative_time +=  Uvar(1,"s");
	    while(!iebprofile.foundDivision(epoch_relative_time))
	      {//No division?
		epoch_relative_time += Uvar(1,"s");//then  add 1 s
	      }	  
	    div_name=iebprofile.getDivisionName(epoch_relative_time);
	    radar_mode= iebprofile.getDivisionMode(div_name) ;
	    //find surface point AGAIN
	    t = epoch_time_ + epoch_relative_time;
	    computeRangeDoppler_3dB_5dB_Borepoints(t);//geom cal
	  }


	//(3) construct ieb
	Ieb ieb;
	ieb.setTime(t);//set time
	ieb.setIebMod(4);//every field is valid
	//----------------------------
	//(3-1)Power mode
	//-----------------------------
	Uvar tfi = t - trigger_time;
	if(tfi < Uvar(0,"s"))
	{
	  cout<<"div name "<<div_name<<endl;
	  cout<<"time "<<t.utc("ISOD")<<endl;
	  cout<<"epoch relative time "<<epoch_relative_time<<endl;
	  cout<<"trigger time "<<trigger_time.utc("ISOD")<<endl;
	  ErrorMessage("tfi is less than 0 s").throwMe();
	}
	ieb.setPowerTfi(tfi);
	ieb.setPmd(17);//normal mode operation

	//------------------
	//(3-2)TNC mode
	//---------------------------
	ieb.setTncTfi(tfi);
    
	//---------------------------
	//(3-3)Slow Fast field
	//---------------------------
	ieb.setFastTfi(tfi);
	ieb.setSlowTfi(tfi);
	ieb.setFin(0);
	ieb.setSin(0);
	ieb.setDtn(data_take_number_);



	//(4) pull out division config parameters
	time_step=iebprofile.getDivisionTimestep(div_name);//add time

	setBeamMask(beam_mask,iebprofile.getDivisionBem(div_name));
	autorad=iebprofile.getDivisionAutorad(div_name);
	baq=iebprofile.getDivisionBaq(div_name);
	noise_bit=iebprofile.getDivisionNoiseBitSetting(div_name);
	//bpd=iebprofile.getDivisionBpd(div_name);
	dutycycle=iebprofile.getDivisionDutycycle(div_name);
	prf=iebprofile.getDivisionPrf(div_name);
	pulse=iebprofile.getDivisionNpulse(div_name);
	bursts_in_flight=iebprofile.getDivisionN_bursts_in_flight(div_name);
	percentOfBW = iebprofile.getDivisionPercentBW(div_name);
    	rip=iebprofile.getDivisionRip(div_name);
	data_rate=iebprofile.getDivisionDatarate(div_name);//kbits /sec
	data_rate *= 1000.0;
	csr=iebprofile.getDivisionCsr(div_name);
        tro=iebprofile.getDivisionTro(div_name);
	//cout<<"time "<<epoch_relative_time<<endl;
	//cout<<"time step "<<time_step<<endl;
	//cout<<"mode "<<radar_mode<<endl;



          double data_rate_factor;
          double noise_bit_factor;


	//(5) set autorad 
	if(autorad=="on" || autorad=="On")
	  {//autorad off
	    ieb.setRL_ND_AutoRad(RL,Max_RL,  ND, Max_ND); 
	  }
	else
	  {//autorad off
	    ieb.setRL_ND(RL, ND); 
	  } 

	// Gh added (03/05/2004): busrts in flight must be 1 for all modes
	// except scatterometery it can be ...anything
/*  RW removed 3/29/2017 to allow BIF for ring obs
	if( (radar_mode!="scatterometer"&&radar_mode!="scatterometer_compressed") && bursts_in_flight !=1 ) {
	  cout << "ERROR: mode is: " << radar_mode ;
	  cout << ".  bursts_in_flight is: " << bursts_in_flight ;
	  cout << ".  Burst_in_flight MUST be 1 if not is scattermeter mode or scatterometer_compressed." << endl;
	  ErrorMessage("Error: (from cfg) bursts_in_flight not 1").throwMe();
	}
*/


	//(6) Let's make ieb corresponding to each mode
	if(radar_mode=="radiometer")
	  {//radiometer
	    if(beam_mask!=4) 
	      {
		cout<<"------------------------------------ "<<endl;
		cout<<"Radiometer not using beam 3 only "<<endl;
		cout<<"------------------------------------ "<<endl;
	      }
	    if(csr!=6) ErrorMessage("Radiometer CSR is not 6").throwMe();
	    ieb.computeMod(altl_adc,csr,baq,beam_mask,percentOfBW);
	    ieb.setPri(1/prf);
	    ieb.setDutycycle(dutycycle);
	    ieb.computeWaveform();
	    ieb.setNumberOfPulses(0);      
	    ieb.setMaxDataRate(data_rate);
	  }
	else if(radar_mode=="scatterometer_compressed")
	  {//scat compressed
	    if(beam_mask!=4) 
	      {
		cout<<"------------------------------------ "<<endl;
		cout<<"Scatterometer not using beam 3 only "<<endl;
		cout<<"------------------------------------ "<<endl;
	      }
	    if(baq!=3) ErrorMessage("Wrong baq mode for scatterometer compression mode").throwMe();
	    ieb.computeMod(altl_adc,csr,baq,beam_mask,percentOfBW);
	    ieb.setPri(1/prf);
	    ieb.setDutycycle(dutycycle);
	    ieb.computeWaveform();
//added 09/30/2011.  Same code from scat mode to check ESS 3rd limit.
     //       if (pulse == 0)
     //         {  // added by RW 070706 to support iapetus-049
                 // (buffer really 32767, but leave some margin 32000)
              // Zero pulses means compute pulse number to fill buffer or rtt.
              // the -1 leaves space for the 1 PRI between xmit and RW
              int pulse_plus_tro_rtt =
                (int)convert_to_double(2*prf*range_3dB_bore_(2,4)/speed_light)
                  - 1;
              int pulse_plus_tro_buf =
                (int)(32000/convert_to_double(altl_adc/prf));
              int pulse_plus_ess_limit = (1 + (70-100*dutycycle)/100) * (prf.getInUnits("Hz")-15.318)/13.221 * (0.9);
              if (pulse_plus_tro_rtt > pulse_plus_tro_buf)
                {  // rtt exceeds buffer, so use buffer limit with space for tro
                pulse = pulse_plus_tro_buf - tro;
                }
              else
                {  // buffer exceeds rtt, so use rtt
                pulse = pulse_plus_tro_rtt - tro;
                }
              if(pulse>pulse_plus_ess_limit)
                {
                  pulse = pulse_plus_ess_limit;

                }
       //       }
// add (above) 09/30/2011
	    ieb.setNumberOfPulses(pulse);
	   
	    //-------------------------------------------------------------
	    //Special method is needed when beam3's boresight points
	    // empty space, not targeted planet.
	    //In this case, use the previous Ieb's RWD,Tro,and CSF
	    //-------------------------------------------------------------
	    if(range_3dB_bore_(2,4) > Uvar(0,"km"))
	      {//beam3's boresight
		ieb.computeALT_Rwd_Tro_Csf(range_3dB_bore_,doppler_3dB_bore_);
	      }
	    else
	      {
		Ieb ieb_tmp= ieb_list_.back();
		ieb.setRwd(ieb_tmp.getRwd());
		ieb.setTroTime(ieb_tmp.getTro());
		ieb.setChirpStartFrequency(ieb_tmp.getChirpStartFrequency());
	      }
	    ieb.setMaxDataRate(data_rate);
	  } 
	else if(radar_mode=="scatterometer")
	  {  //scatt
	  if (beam_mask!=4) 
	    {
            cout<<"------------------------------------ "<<endl;
            cout<<"Scatterometer not using beam 3 only "<<endl;
            cout<<"------------------------------------ "<<endl;
	    }
	  if (baq!=5)
            {
            cout<<"-----------------------------"<<endl;
            cout<<"IebList.cpp::scatt mode not using 8-8 straight"<<endl;
            cout<<"-----------------------------"<<endl;
	    }
	  if (range_3dB_bore_(2,4) > Uvar(0,"km"))
            {  // beam3's boresight on the surface
	    ieb.computeMod(altl_adc,csr,baq,beam_mask,percentOfBW);
	    if (prf == 0)
              {  // added by RW 070706 to support iapetus-049
              // When division PRF is zero treat this scatterometer division
              // as an imaging scatt division which uses the PRF profile.
              prf=iebprofile.getScatPrf(epoch_time_ + epoch_relative_time);  
              }
	    else if (prf == -1)
	      { // added by RW 070706 to support scatterometry
		//cout<<"PRF PRofile Test:"<<endl;
		prf = computePRF(epoch_time_ + epoch_relative_time);
		cout<<"PRF "<<(epoch_relative_time)/60<<" "<<prf<<endl;
		
	      }
	    ieb.setPri(1/prf);
	    ieb.setDutycycle(iebprofile.limitDutyCycle(dutycycle,prf));
	    ieb.computeWaveform();
            if (pulse == 0)
              {  // added by RW 070706 to support iapetus-049
                 // (buffer really 32767, but leave some margin 32000)
              // Zero pulses means compute pulse number to fill buffer or rtt.
              // the -1 leaves space for the 1 PRI between xmit and RW
              int pulse_plus_tro_rtt =
                (int)convert_to_double(2*prf*range_3dB_bore_(2,4)/speed_light)
                  - 1;
              int pulse_plus_tro_buf =
                (int)(32000/convert_to_double(altl_adc/prf));
	      int pulse_plus_ess_limit = (1 + (70-100*dutycycle)/100) * (prf.getInUnits("Hz")-15.318)/13.221 * (0.9);
              if (pulse_plus_tro_rtt > pulse_plus_tro_buf)
                {  // rtt exceeds buffer, so use buffer limit with space for tro
                pulse = pulse_plus_tro_buf - tro;
                }
              else
                {  // buffer exceeds rtt, so use rtt
                pulse = pulse_plus_tro_rtt - tro;
                }
	      if(pulse>pulse_plus_ess_limit)
		{
		  pulse = pulse_plus_ess_limit;
		 
		}
              }
            ieb.setNumberOfPulses(pulse);
            if (csr !=1)
              {
	      if (set_csf_times_pri_integer_ && percentOfBW == 0)
                {  //special chirp selection
		Uvar special_chirp =
                  findToneChirpStartFrequency(ieb.getPri(), ieb.getAdc());
		ieb.computeToneALT_Rwd_Tro(range_3dB_bore_, special_chirp);
	        }
	      else
                {
		//ieb.computeALT_Rwd_Tro_Csf(range_3dB_bore_,doppler_3dB_bore_);
                ieb.computeALT_Rwd_Tro_Csf(range_3dB_bore_,
                  doppler_3dB_bore_,tro);
	        }
              }
            else
              {
              ieb.setRwd(0);
              ieb.setTroTime (0);
              ieb.setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz") +
                (altl_chirp_bandwidth*percentOfBW)/2.0 ) ;
              }
	    }
	  else
            {  // beam3's boresight off the disk
	    if (csr !=1)
              {
              Ieb ieb_tmp= ieb_list_.back();
              ieb.setRwd(ieb_tmp.getRwd());
              ieb.setTroTime(ieb_tmp.getTro());
              ieb.setChirpStartFrequency(ieb_tmp.getChirpStartFrequency());
	      } 
            else
              {
              ieb.computeMod(altl_adc,csr,baq,beam_mask,percentOfBW);
              ieb.setPri(1/prf);
              ieb.setDutycycle(dutycycle);
              ieb.computeWaveform();
              ieb.setNumberOfPulses(pulse);
              ieb.setRwd(0);
              ieb.setTroTime (0);
              ieb.setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz") +
              (altl_chirp_bandwidth*percentOfBW)/2.0 ) ;
	      }
	    }
	  ieb.setMaxDataRate(data_rate);
	  }
	else if(radar_mode=="altimeter")
	  {
	    if(beam_mask!=4) 
	      {
		cout<<"------------------------------------ "<<endl;
		cout<<"Altimeter not using beam 3 only "<<endl;
		cout<<"------------------------------------ "<<endl;
	      }
	    ieb.computeMod(alth_adc,csr,baq,beam_mask,percentOfBW);
   //gh 20110708 getprofile for alth mode 
         //   if (prf == -1)
         //     { // added by RW 070706 to support scatterometry
         //       cout<<"PRF PRofile Test:"<<endl;
         //       prf = computePRF(epoch_time_ + epoch_relative_time);
         //       cout<<"prf"<<prf<<endl;
         //     }
   //gh 20110708

	    ieb.setPri(1/prf);
	    ieb.setDutycycle(dutycycle);
	    ieb.computeWaveform();
	    ieb.setNumberOfPulses(pulse);

            if (csr==1) {
              ieb.setRwd(0);
              ieb.setTroTime (0);
              ieb.setChirpStartFrequency(Uvar(5375.0*1000.0,"Hz") + (alth_chirp_bandwidth*percentOfBW)/2.0 ) ;
            } else {
	      if(range_3dB_bore_(2,4) <= Uvar(0,"km"))
		ErrorMessage("IebList.cpp::No surface intercept for altitmeter mode:Pointing design error").throwMe();
	      ieb.computeALT_Rwd_Tro_Csf(range_3dB_bore_,doppler_3dB_bore_);
	    }
	    
	    ieb.setMaxDataRate(data_rate);
	  }
	else if(radar_mode=="sarl" || radar_mode=="sarh")
	  {
	    if(beam_mask!=31 && beam_mask!=4)// all Beams, but can also be Beam 3 only???? 
	      {
		cout<<"------------------------------------ "<<endl;
		cout<<"SAR not using all 5beams "<<endl;
		cout<<"------------------------------------ "<<endl;
	      }
	    if(radar_mode=="sarh")
	      {
		ieb.computeMod(sarh_adc,csr,baq,beam_mask,percentOfBW);
	      }
	    else
	      {
		ieb.computeMod(sarl_adc,csr,baq,beam_mask,percentOfBW);
	      }
	    if (csr == 1) prf=iebprofile.getDivisionPrf(div_name);
	    else{ 
	      //-----------------------------------
	      //PRF Setting for SAR
	      //----------------------------------

	      //-----------------------------------
	      //(Special case) When need to use a fixed PRI/PRF
	      //----------------------------------
	      if(prf > Uvar(0.0001,"Hz")){
		//Normal PRF values in config file should be 0
		if(fabs(epoch_relative_time)<Uvar(30*60,"s")){
		  //ErrorMessage("IebList.cpp: Trying to use a fixed PRF inside +-30 min w.r.t. the closest approach time").throwMe();
		  cout <<"Warning: IebList.cpp: Trying to use a fixed PRF inside +-30 min w.r.t. the closest approach time" << endl;
		  }
		ieb.setPri(1/prf);
	      }
	      //-----------------------
	      //(Nominal case) Use PRF profile
	      //--------------------------
	      else{
		prf=iebprofile.getSarPrf(epoch_relative_time);  
		ieb.setPri(1/prf);
	      }
	    }
	    ieb.setDutycycle(dutycycle);
	    ieb.computeWaveform();
	    
	    if (csr==1) {
              ieb.setNumberOfPulses(pulse);
	      ieb.setRwd(0);
 	      ieb.setTroTime (0);
	      
  	      if (radar_mode =="sarh") ieb.setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz") + (sarh_chirp_bandwidth*percentOfBW)/2.0 ) ;
  	      if (radar_mode =="sarl") ieb.setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz") + (sarl_chirp_bandwidth*percentOfBW)/2.0 );
	      
	    } else  ieb.computeSAR_Pul_Rwd_Tro_Csf(range_5dB_bore_,doppler_5dB_bore_);
	    
	    ieb.setMaxDataRate(data_rate);
	  }
	/////////////////////// Start Ping Pong Mode/////////////////////
	else if(radar_mode=="sar_ping_pong")
	  {
	    Uvar div_end_time = cfg_["div_end_time_" + div_name]; 
	    Uvar div_start_time = cfg_["div_start_time_" + div_name]; 
	    //double data_rate_factor;
	    //double noise_bit_factor;
	    
	    //---------------------
	    // Assigning values to:
	    // 1.data rate factor
	    // 2.noise bit factor
	    //---------------------
	    if (cfg_.keywordExists("Ping_Pong_Data_Rate_Factor") && cfg_.keywordExists("Ping_Pong_Data_Rate_Factor"))
	      {
		data_rate_factor = cfg_.getDouble("Ping_Pong_Data_Rate_Factor");
		noise_bit_factor = cfg_.getDouble("Ping_Pong_Noise_Bit_Factor"); 
	      }
	    else
	      {
		cout<<"Assigning default noise bit and data rate factor = 1"<<endl;
		data_rate_factor = 1;
		noise_bit_factor = 1;
	      }

	    
	    //-------------------------------
	    // Alternating Sar Hi and Sar Low
	    //-------------------------------
	    //while (epoch_relative_time < div_end_time)
	    
	    // {
	    ping_pong_counter = round_double(((epoch_relative_time-div_start_time)/time_step).getInUnits(""));
	   
	    
	    if(beam_mask!=31 && beam_mask!=4)// all Beams, but can also be Beam 3 only???? 
	      {
		cout<<"------------------------------------ "<<endl;
		cout<<"SAR not using all 5beams "<<endl;
		cout<<"------------------------------------ "<<endl;
	      }
	    cout<<"Epoch time: "<<epoch_relative_time<<endl;
	    cout<<"Ping pong counter: "<<ping_pong_counter<<endl;
	    if(ping_pong_counter%2 == 1)
	      {cout<<"Sar High"<<endl;
		ieb.computeMod(sarh_adc,csr,baq,beam_mask,percentOfBW);
	      }
	    else
	      {
		cout<<"Sar Low"<<endl;
		ieb.computeMod(sarl_adc,csr,baq,beam_mask,percentOfBW);
	      }
	    if (csr == 1) prf=iebprofile.getDivisionPrf(div_name);
	    else{ 
	      //-----------------------------------
	      //PRF Setting for SAR
	      //----------------------------------
	      
	      //-----------------------------------
	      //(Special case) When need to use a fixed PRI/PRF
	      //----------------------------------
	      if(prf > Uvar(0.0001,"Hz")){
		//Normal PRF values in config file should be 0
		if(fabs(epoch_relative_time)<Uvar(30*60,"s")){
		  //ErrorMessage("IebList.cpp: Trying to use a fixed PRF inside +-30 min w.r.t. the closest approach time").throwMe();
		  cout <<"Warning: IebList.cpp: Trying to use a fixed PRF inside +-30 min w.r.t. the closest approach time" << endl;
		}
		ieb.setPri(1/prf);
	      }
	      //-----------------------
	      //(Nominal case) Use PRF profile
	      //--------------------------
	      else{
		prf=iebprofile.getSarPrf(epoch_relative_time);  
		ieb.setPri(1/prf);
	      }
	    }
	    ieb.setDutycycle(dutycycle);
	    ieb.computeWaveform();
	    
	    if (csr==1) {
	      ieb.setNumberOfPulses(pulse);
	      ieb.setRwd(0);
	      ieb.setTroTime (0);
	      
	      if (ping_pong_counter%2 == 1)
		{
		  ieb.setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz") + (sarh_chirp_bandwidth*percentOfBW)/2.0 ) ;
		}
	      else
		{
		  ieb.setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz") + (sarl_chirp_bandwidth*percentOfBW)/2.0 );
		}
	      
	    } else  ieb.computeSAR_Pul_Rwd_Tro_Csf(range_5dB_bore_,doppler_5dB_bore_);
	    
	    if (ping_pong_counter%2 == 1)
	      {
		ieb.setMaxDataRate(data_rate);
	      }
	    else
	      {
		ieb.setMaxDataRate(data_rate*data_rate_factor);
	      }
	    /* //(7) set RIP and compute radiometer window count
	    //every ieb shares these two steps
	    ieb.setRip(rip);
	    ieb.computeRadandBpd(bursts_in_flight);
	    
	    //(8) set attenuators
	    
	    
	    //--------------------------------------
	    //attenuation setting
	    // when auto gain is enabled, attenuation setting will be done
	    // once when radar mode changes.  Otherwise, do it at every time step
	    //case: first element
	    //      mode change when auto gain is on
	    //      do it for normal operation
	    //------------------------------------
	    Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
	    if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
	    if(ping_pong_counter%2 == 1)
	      {
		ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_sarh_);
	      }
	    else 
	      {
		ieb.computeAttenuationSettings(noise_bit*noise_bit_factor, Pn, squared_deviation_of_system_noise_cal_sarl_);
	      } 
	    
	    //--------------------------------------------------
	    //(9) compute bii: bii can be reset after finishing the construction
	    // of entire ieb list
		//-------------------------------------
		//compute bii: int(ieb_time_step / bpd)
		//-----------------------------------------
	    unsigned int Nburst= (unsigned int) (time_step/ieb.getBpd()).getInUnits("");
	    if (Nburst==0)  Nburst = 1;//never set to 0 (continuous burst mode)
	    unsigned int Nbeam =ieb.getNumberOfBeams();
	    unsigned int Nbeam_remain = Nburst % Nbeam;
	    if (Nbeam_remain != 0) Nburst = (Nburst/Nbeam +1 )*Nbeam;
	    if(Nburst > 255) Nburst = 255;
	    ieb.setBii(Nburst);
	    ieb.encodeSlowFastfield();
	    ieb.encodePowerTncfield();
	    addIeb(ieb);
	    //cout<<"data points per pri "<<ieb_pr->get   
	    //(10) set next time for a new ieb   
	    epoch_relative_time +=time_step;*/
	    // }
	  }
	
	///////////////////// End Ping Pong Mode //////////////////////
	else
	  {
	    ErrorMessage("IebList.cpp::requested mode does not exist:"+radar_mode).throwMe();
	  }
	
	//if(ping_pong_counter == 0)
	// {
	//(7) set RIP and compute radiometer window count
	//every ieb shares these two steps
	ieb.setRip(rip);
	ieb.computeRadandBpd(bursts_in_flight);
	
	//(8) set attenuators
	cout<<"Set attenuators:"<<endl;
	if ( radar_mode !="radiometer")
	  {
	    
	    //ieb will set att for radiometer automatically
	    //--------------------------------------
	    //attenuation setting
	    // when auto gain is enabled, attenuation setting will be done
	    // once when radar mode changes.  Otherwise, do it at every time step
	    //case: first element
	    //      mode change when auto gain is on
	    //      do it for normal operation
	    //------------------------------------
	    
	    if(radar_mode=="scatterometer" || radar_mode=="scatterometer_compressed"){
	      Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
	      if(noise_bit!=4.0) cout<<"Warning: noise bit is not the expected setting of  4 "<<endl;
	      ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_altl_);
	    }
	    else if(radar_mode=="altimeter"){
	      Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
	      if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
	      ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_alth_);
	    }
	    else if(radar_mode=="sarl" ){
	    //else if(radar_mode=="sarl" || (radar_mode=="sar_ping_pong" && ping_pong_counter%2 == 0)){
	      Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
	      if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
	      ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_sarl_);
	    }
            else if(radar_mode=="sarh" ){
	    //else if(radar_mode=="sarh" ||(radar_mode=="sar_ping_pong" && ping_pong_counter%2 == 1) ){
	      Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
	      if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
	      ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_sarh_);
	    }
            else if( (radar_mode=="sar_ping_pong" && ping_pong_counter%2 == 0)){
              Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
              if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
              ieb.computeAttenuationSettings(noise_bit*noise_bit_factor, Pn, squared_deviation_of_system_noise_cal_sarl_);
            }
            else if((radar_mode=="sar_ping_pong" && ping_pong_counter%2 == 1) ){
              Uvar Pn = Tsys_ * boltzmann_constant*ieb.getRcv();
              if(noise_bit!=2.0) cout<<"Warning: noise bit is not the expected setting of  2 "<<endl;
              ieb.computeAttenuationSettings(noise_bit, Pn, squared_deviation_of_system_noise_cal_sarh_);
            }
	    else{
	      ErrorMessage("Invalid radar mode").throwMe();
	    }
	  }
	
	//--------------------------------------------------
	//(9) compute bii: bii can be reset after finishing the construction
	// of entire ieb list
	//-------------------------------------
	//compute bii: int(ieb_time_step / bpd)
	//-----------------------------------------
	unsigned int Nburst= (unsigned int) (time_step/ieb.getBpd()).getInUnits("");
	if (Nburst==0)  Nburst = 1;//never set to 0 (continuous burst mode)
	unsigned int Nbeam =ieb.getNumberOfBeams();
	unsigned int Nbeam_remain = Nburst % Nbeam;
	if (Nbeam_remain != 0) Nburst = (Nburst/Nbeam +1 )*Nbeam;
	if(Nburst > 255) Nburst = 255;
	ieb.setBii(Nburst);
	ieb.encodeSlowFastfield();
	ieb.encodePowerTncfield();
	addIeb(ieb);
	//cout<<"data points per pri "<<ieb_pr->get   
	//(10) set next time for a new ieb   
	epoch_relative_time +=time_step;
cout<<"Epoch time: "<<epoch_relative_time<<endl;
      }	
    // ping_pong_counter = 0;
    // }

	//-----
	//Do not delete
	//-------
	
	//Uvec Echo_power("echo power",5); 
	//
	//if (ieb.getAdc() == altl_adc)
	//  {// noise level is 20 dB higher than quantization noise level
	//ieb.computeAttenuationSettings(Pn, 20.0,
	//	   squared_deviation_of_system_noise_cal_altl_);
	//  }
	//else
	//  {//noise level is 3dB higher than quantization noise level
	//ieb.computeAttenuationSettings(Pn, 3.0,
	// squared_deviation_of_system_noise_cal_altl_);
	//  }
	///cout<<"attenuation setting "<<ieb_pr->getAt1()
	//  <<" "<<ieb_pr->getAt3()<<" "<<ieb_pr->getAt4()<<endl;
	//
	//----------------------------------------------
	//Do not delete
	//-------------------------------------------------
	//if (mode_change==1 || ieb_pr->getCsr()==0)
	//{
	//  for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
	//    {	      
	//    unsigned int bem = ieb_pr->getBem();
	//    unsigned short bem_mask;
	//    bitset(bem_mask,0,4,bem);
	//    unsigned int beam_status=bitget(bem_mask,i_beam,i_beam);
	//    if (beam_status == 1){
	//computeSarRadarEcho(i_beam,t,
	//		ieb_pr->getPri() * ieb_pr->getDutycycle(),
	//		Echo_power(i_beam));}
	//    else Echo_power(i_beam) = Pn;
	//    }
	//  ieb_pr->computeAttenuationSettings(Pn, 
	//  			     Echo_power,
	//     squared_deviation_of_system_noise_cal_altl_);
	// 
	//}
	//else
	//{
	//  ieb_pr->setAttenuation((ieb_pr-1)->getAt1(),
	//		   (ieb_pr-1)->getAt3(),	
	//		   (ieb_pr-1)->getAt4());
	//}
	//cout<<"attenuation setting "<<ieb_pr->getAt1()
	//  <<" "<<ieb_pr->getAt3()<<" "<<ieb_pr->getAt4()<<endl;
	
  }



//------------------------
//compute radar echo within beam3dB or pulserange
//----------------------------
void IebList::computeSarRadarEcho(const unsigned int& beam_number,
				  const Time& t, 
				  const Uvar& taup,
				  Uvar& echo)
  {
    StateVector sc_state;
    target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
    Uvar lambda = speed_light/carrier_frequency;
    echo = Uvar(0,"kg km km/(s s s)");
    //calculate the strength of return echo for each beam
    DirectionVector beam_boresight("boresight",Fvector_[beam_number],t,0,0,1); 
    TargetGeom tg_bore(t);
    tg_bore.setState(sc_state);
    tg_bore.setTarget(target_name_,target_frame_);
    tg_bore.setLookDirection(beam_boresight);
    Uvar beam_range = tg_bore.range();
    Uvar beam_inc_angle = tg_bore.incidenceAngle();
    Uvar radius = tg_bore.radius();
    Uvar azimuth_start,azimuth_end,elev_start,elev_end;
    azimuth_start = azimuth_3dB_(beam_number,0);
    azimuth_end = azimuth_3dB_(beam_number,1);
    if (azimuth_start > azimuth_end) ErrorMessage("Ieb-computeAtt: azimuth start is larger than azimuth end").throwMe();
    
    Uvar beam_3dB_extent = elevation_3dB_(beam_number,1)
      - elevation_3dB_(beam_number,0);
    if (beam_3dB_extent < Uvar(0,"rad")){ 
      ErrorMessage("Ieb-computeAtt: elevation start is larger than elevation end").throwMe();}
    
    double eff_elev_extent =  (Uvar("speed_light")*taup
			       /2.0/beam_range/tan(beam_inc_angle)).getInUnits("");
    elev_start = Uvar(-eff_elev_extent/2.0,"rad");
    elev_end = Uvar(eff_elev_extent/2.0,"rad");
    if ((elev_end - elev_start) > beam_3dB_extent)
      {
	elev_start = elevation_3dB_(beam_number,0);
	elev_end = elevation_3dB_(beam_number,1);
      }
    
    //cout<<"beam number elev start and end "<<beam_number<<" "
    //<<elev_start.getInUnits("deg")<<" "<<elev_end.getInUnits("deg")
    //<<endl;
    
    //calculate echo power for a pixel size of 10 x 10
    unsigned int Ngrid = 11;
    unsigned int Npixel = Ngrid-1;
    Array2D<PositionVector> Surface_position("surface position",Ngrid,Ngrid);
    Umat Surface_range("surface range",Ngrid,Ngrid);
    Dmat Beam_gain("beam gain",Ngrid,Ngrid);
    Umat Surface_incidenceAngle("surface incidence angle",Ngrid,Ngrid);
    
    TargetGeom tg_look(t);
    for (unsigned int i_azi = 0; i_azi < Ngrid;++i_azi)
    for(unsigned int i_elev = 0; i_elev < Ngrid;++i_elev)
      {
	Uvar azi = azimuth_start ;
	azi += (azimuth_end-azimuth_start)*double(i_azi)/double(Ngrid-1);
	Uvar elev = elev_start;
	elev += (elev_end-elev_start)*double(i_elev)/double(Ngrid-1);
	
	DirectionVector look("bore",Fvector_[beam_number],t,0,0,1);
	look.setAzimuthElevation(azi,elev);//set azimuth and elev
	
	tg_look.setState(sc_state);
	tg_look.setTarget(target_name_,target_frame_);
	tg_look.setLookDirection(look); 
	
	Beam_gain(i_azi,i_elev) = Bvector_[beam_number].getGainOneWay(look);
	Surface_range(i_azi,i_elev) = tg_look.range();
	Surface_position(i_azi,i_elev)=tg_look.surfaceIntercept();
	Surface_incidenceAngle(i_azi,i_elev)=tg_look.incidenceAngle();
	
	tg_look.reset(t);
      }//loop over azimuth and elevation

    //calculate return power for each pixel
    for (unsigned int i_azi = 0; i_azi < Npixel;++i_azi)
    for(unsigned int i_elev = 0; i_elev < Npixel;++i_elev)	
      {
	//gain
	double gain = (Beam_gain(i_azi,i_elev) 
		       + Beam_gain(i_azi+1,i_elev)
		       + Beam_gain(i_azi,i_elev+1) 
		       + Beam_gain(i_azi+1,i_elev+1))/4.0;
	Uvar angle = (Surface_incidenceAngle(i_azi,i_elev)
		      +Surface_incidenceAngle(i_azi+1,i_elev)
		      +Surface_incidenceAngle(i_azi,i_elev+1)
		      +Surface_incidenceAngle(i_azi+1,i_elev+1))/4.0;
	Uvar range_to_surface = ( Surface_range(i_azi,i_elev) 
				  +Surface_range(i_azi+1,i_elev) 
				  +Surface_range(i_azi,i_elev+1) 
				  +Surface_range(i_azi+1,i_elev+1))/4.0;
	Uvar area = surfacearea(Surface_position(i_azi,i_elev),
				Surface_position(i_azi+1,i_elev),
				Surface_position(i_azi,i_elev+1),
				radius);
	area += surfacearea(Surface_position(i_azi+1,i_elev+1),
			    Surface_position(i_azi+1,i_elev),
			    Surface_position(i_azi,i_elev+1),
			    radius);
	
	echo+=Pt_*gain*gain*lambda*lambda*area 
	  *muhleman_backscatter( muhleman_k1_, muhleman_k2_ ,angle)
	  /pow(4.0*pi,3)/pow(range_to_surface,4);	     
	} //loop over surface points	
    cout<<"beam number and echo "<< beam_number+1 <<" "<<echo<<endl;
  }

//------------------------------------------
//void computeUsableCrosstrackExtent(Uvar& extent,
//                                   Uvar& gap,
//                                   const Time& t,
//                                   const Time& time_offset,
//                                   const Frame& track_frame)
//  Basic assumption: assume beam1,3,5 send pulses at time t 
// and beam2 and 4 at t + time_offset
//------------------------------------------
void IebList::computeUsableCrosstrackExtentandGap(Uvar& extent,
					    Uvar& gap,
					    const Time& t, 
					    const Uvar& time_offset,
					    const Frame& track_frame)
  {
    //initialize variables
    extent = Uvar(0,"km");
    gap = Uvar(0,"km");
    //Build 5 ambiguityFNN class for 5 beams
    Uvar lambda = speed_light/carrier_frequency;
    amb.clear();
    amb.resize(5);
    nlooks.clear();
    nlooks.resize(5);
    ML_area.clear();
    ML_area.resize(5);
   
    Ieb ieb_t = getIeb(t);
    ieb_t.decodeSlowFastfield();
    Ieb ieb_t_offset=getIeb(t+time_offset);
    ieb_t_offset.decodeSlowFastfield();
    
   
    for (unsigned int i = 0; i < 5;++i)
      {
	unsigned int beam_num = i+ 1;
	Time burst_time= t;	
	Ieb ieb;
	if (beam_num == 1 || beam_num ==3 || beam_num==5)
	  {//beam 1 3 5
	    ieb = ieb_t;
	  }
	else
	  {//beam 2 and 4
	    ieb = ieb_t_offset;
	    burst_time += time_offset;//take into account time offest 
	  }
       
	//-----------------------------------------
	//calculate range and azimuth resolution
	//-------------------------------------------
	StateVector sc_state("sc_state");
	target_frame_.ephemeris(sc_state,"Cassini",burst_time,"NONE");
	DirectionVector boresight("boresight",Fvector_[i],burst_time,0,0,1);
	//make a correction to beam's max gain
	//Uvar azimuth=Bvector_[i].elevAndMaxGainToAzim(Uvar(0,"rad"));
	//boresight.setAzimuthElevation(azimuth,Uvar(0,"rad"));
	//set ground range azimuth resolutions
	FloatVector velocity = sc_state.velocity();
	velocity.representIn(Fvector_[i]);//beam frame
	DirectionVector velocity_dir = velocity;
	Uvar Vst_bore("along_track_speed");
	Uvar speed = velocity.magnitude();
	Vst_bore = velocity.magnitude();
	Vst_bore *= sqrt(1.0 - pow(dot(velocity_dir,boresight),2));
	
	TargetGeom tg(burst_time);
	tg.setState(sc_state);
	tg.setTarget(target_name_,target_frame_);
	tg.setLookDirection(boresight);
	Uvar range0 = tg.range();
	Uvar doppler0=tg.doppler(lambda);
	Uvar incidenceAngle = tg.incidenceAngle();
	int N_p = int(ieb.getPul());//number of transmitted pulses
	Uvar receive_window = double(N_p) *ieb.getPri();
	Uvar x_res = lambda * range0/(2.0 * Vst_bore * receive_window);//azi res
	Uvar r_res = speed_light/(2.0 * ieb.getBr());//range resolution
	Uvar rg_res = r_res/sin(incidenceAngle);//ground resolution  
	
	//----------------
	//pulse gate
	//---------------
	Uvar azimuth_range = range0 
	  * Bvector_[i].getAzimuthWidthOneWay().getInUnits("rad");
	Uvar elev_range = range0 
	  * Bvector_[i].getElevationWidthOneWay().getInUnits("rad");
	Uvar range_3db = elev_range * tan(incidenceAngle);
	Uvar tau_3db = 2.0 * range_3db/speed_light;//range gate
	Uvar pulse_gate =(ieb.getTaup()+tau_3db)*speed_light/2.0;
	//cout<<"pulse gate "<<pulse_gate<<endl;
	//cout<<"prf "<<1.0/ieb.getPri()<<endl;
	//cout<<"pbw "<<pbw_ratio_/ieb.getPri()<<endl;
	//cout<<"bore range "<< range0<<endl;
	amb[i].setTime(burst_time);
	amb[i].config(cfg_);//transfer grid size
	amb[i].setTarget(target_name_,target_frame_);
	amb[i].computeState();
	amb[i].setBeam(Bvector_[i]);
	amb[i].setProcessWindow(range0,doppler0,
				1.0/ieb.getPri(),pbw_ratio_/ieb.getPri(),
				pulse_gate);
	amb[i].setTrackFrame(track_frame);
	amb[i].calAmbGeometry();//obtain range, doppler, area, gain
	
	//cout<<"range "<< amb[i].range<<endl;
	
	//-----------------------
	//radar equation prefactor
	//-----------------------
	double Nr = (ieb.getTaup()).getInUnits("s") 
	  *(ieb.getBr()).getInUnits("Hz"); 
	double Ni = double(N_p)*Nr;
	Uvar X0 = Pt_ * Ni/pow(4.0*pi,3);	
	Uvar Pn = boltzmann_constant * Tsys_ * ieb.getRcv();
	//Uvar usable_area = amb[i].calUsableArea(X0,Pn,x_res,rg_res);
	//cout<<"usbel area"<<usable_area<<endl;
	//cout<<"noise "<<Pn<<endl;
	unsigned int Nbeam = ieb.getNumberOfBeams();
	Uvar SL_area=amb[i].calUsableArea(X0,Pn,x_res,rg_res);
        ML_area[i]
	  = amb[i].calMultilookUsableArea(X0,Pn,x_res,rg_res,double(Nbeam)*ieb.getBpd(),velocity,N_p);
	
	nlooks[i] = amb[i].getNumberofLooks();

	
      }//loop over all 5 beams
    
    //-------------------------------------------------------------------
    //once we calculated usable area, we can calculate usable crosstrack
    // extent of usable area's center track: strip
    // Data Structure
    //  thermal snr, amb ratio, beam_gain: square matrix 
    //  let's take the middle value: dop = dop_bore, 
    //                               range- inner most to outer most
    //-----------------------------------------------------------------
    unsigned int Nrange_bin = amb[0].getGridsize();
    unsigned int Ndop_bin= Nrange_bin;//square grid
    Umat cross_track_strip("cross track extent",5,Nrange_bin+1);
    Umat thermal_snr_strip("thermal snr strip",5,Nrange_bin);
    Umat amb_ratio_strip("amb strip",5,Nrange_bin);
    Umat beam_gain_strip("beam gain strip",5,Nrange_bin);
    Uvec min_cross_track_strip("min ",5);
    Uvec max_cross_track_strip("max",5);

    for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
      {
	for (unsigned int i = 0; i < Nrange_bin+1;++i)
	  { 
	    cross_track_strip(i_beam,i) =(amb[i_beam].crosstrack)(i,Ndop_bin/2);
	  }
	Uvar min, max;
	min_max(min,max,cross_track_strip);
	min_cross_track_strip(i_beam) = min;
	max_cross_track_strip(i_beam) = max;
      }
    for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
    for (unsigned int i = 0; i < Nrange_bin;++i)
      { 
	thermal_snr_strip(i_beam,i)=(amb[i_beam].thermal_snr_ML)(i,Ndop_bin/2);
	amb_ratio_strip(i_beam,i) = (amb[i_beam].amb_ratio_ML)(i,Ndop_bin/2);
	beam_gain_strip(i_beam,i)=(amb[i_beam].oneway_beam_gaindB)(i,Ndop_bin/2);
      }
    
    Uvar min_cross, max_cross;
    min_max(min_cross,max_cross,cross_track_strip);//min and max of crosstrack
    //cout<<"strip processor "<< min_cross<<" "<<max_cross<<endl;
    
    //---------------------------------------------------
    //The distance between max_cross and min_cross: max_cross - min_cross
    // number of grid between them: max_cross - min_cross + 1
    //----------------------------------------------------
    unsigned int Ncross_bin =(unsigned int)(max_cross- min_cross).getInUnits("km")+1;
    vector<int> good_area_counter;
    good_area_counter.clear();

    for (unsigned int i = 0; i < Ncross_bin;++i)
      {
	Uvar cross_track_value = min_cross 
	  + (max_cross-min_cross)*double(i)/double(Ncross_bin-1);
	
	//------------------------------------------------------
	//check the above value satisfies usable area conditions
	//  in terms of snr, amb,and beam gain for all 5beams!
	//-----------------------------------------------------
	list<Uvar> amb_collector, snr_collector,beam_collector;
	for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
	  {
	    if(cross_track_value >=min_cross_track_strip(i_beam)
	       && cross_track_value < max_cross_track_strip(i_beam))
		{
		  for (unsigned int j = 0; j < Nrange_bin;++j)
		  { 
		    if (cross_track_strip(i_beam,j) <= cross_track_strip(i_beam,j+1))
		    {//right looking geometry, crosstrack increases with j
		      if( cross_track_value >= cross_track_strip(i_beam,j) &&
			  cross_track_value < cross_track_strip(i_beam,j+1))
		      {
			amb_collector.push_back(amb_ratio_strip(i_beam,j));
			snr_collector.push_back(thermal_snr_strip(i_beam,j));
			beam_collector.push_back(beam_gain_strip(i_beam,j));
		      }
		    }
		    else 
		    {//left looking geometry,crosstrack decreases with j
		      if( cross_track_value <= cross_track_strip(i_beam,j) &&
			  cross_track_value > cross_track_strip(i_beam,j+1))
		      {
			amb_collector.push_back(amb_ratio_strip(i_beam,j));
			snr_collector.push_back(thermal_snr_strip(i_beam,j));
			beam_collector.push_back(beam_gain_strip(i_beam,j));
		      }
		    }
		  }//loop over beam and range bins
		}//if crosstrack is between min and max of each beam's strip
	  }//for each beam

	//snr is here defined as -(noise_equivalent_sigma0)
	if (amb_collector.size() == 0) 
	  {
	    extent += Uvar(0,"km");
	  }
	else if (amb_collector.size() == 1)
	  {
	    if (amb_collector.front() >= amb_ratio_ &&
		snr_collector.front() >= -noise_equiv_sigma0_ &&
		beam_collector.front() >=min_beam_gain_)
	      {
		extent+= Uvar(1,"km");//usable crosstrack
		good_area_counter.push_back(int(i));//good area
	      }
	  }
	else
	  {
	    //---------------------
	    //Here, we are taking max amb, snr, and beam gain of each
	    // pixel imaged by more than one beams
	    //----------------------
	    amb_collector.sort();
	    snr_collector.sort();
	    beam_collector.sort();
	    if ( amb_collector.back() >= amb_ratio_ &&
		 snr_collector.back() >= -noise_equiv_sigma0_ &&
		 beam_collector.back() >=min_beam_gain_) 
	      { 
		extent +=Uvar(1,"km");//usable crosstrack
		good_area_counter.push_back(int(i));//good area
	      }
	    else
	      {//bad area
		//cout<<"amb "<<amb_collector.back()<<endl;
		//cout<<"snr  "<<snr_collector.back()<<endl;
		//cout<<"beam "<<beam_collector.back()<<endl;
		//cout<<"bad area crosstrack value"<<cross_track_value<<endl;
	      }
	  }
      }//loop over next range bin
  
    // cout<<"track extent"<<extent<<endl;
    //find gap across good area
    if (good_area_counter.size() == 0 || good_area_counter.size() == 1)
      {
	gap = Uvar(0,"km");
      }
    else
      {
	for (vector<int>::iterator pr=good_area_counter.begin();
	     pr <(good_area_counter.end()-1);++pr)
	  {
	    int a = *(pr+1);
	    int b = *pr;
	    if ((a-b) != 1)
	      {//there is a gap, like 2,3,4,6 where 5 is missing 
		int c = a - b - 1;
		gap += Uvar(c ,"km");
	      }
	  }
      }
    //cout<<"gap "<<gap<<endl;    
  }


 //return ambiguity output

void 
IebList::getAmbiguity(vector<AmbiguityFNN>& amb_output){

  amb_output.clear();

  for(unsigned int i=0;i<5;++i){
    //cout<< amb[i].range<<endl;
    amb_output.push_back(amb[i]);
  }
  
}


//------------------------------------------
//void computeUsableCrosstrackExtent(Uvar& extent,
//                                   Uvar& gap,
//                                   const Time& t,
//                                   const Ieb& ieb135,
//                                   const Ieb& ieb24,
//                                   const Time& time_offset,
//                                   const Frame& track_frame)
//  Here, we do not update ieb while geometry is updated
// 
//------------------------------------------
void IebList::computeUsableCrosstrackExtentandGap(Uvar& extent,
						  Uvar& gap,
						  const Time& t, 
						  const Ieb& ieb135,
						  const Uvar& time_offset,
						  const Ieb& ieb24,
						  const Frame& track_frame)
  {
    //initialize variables
    extent = Uvar(0,"km");
    gap = Uvar(0,"km");
    //Build 5 ambiguityFNN class for 5 beams
    Uvar lambda = speed_light/carrier_frequency;
    amb.clear();
    amb.resize(5);
    nlooks.clear();
    nlooks.resize(5);
    ML_area.clear();
    ML_area.resize(5);
    for (unsigned int i = 0; i < 5;++i)
      {
	unsigned int beam_num = i+ 1;
	Time burst_time= t;	
	Ieb ieb;
	if (beam_num == 1 || beam_num ==3 || beam_num==5)
	  {//beam 1 3 5
	    ieb = ieb135;
	  }
	else
	  {//beam 2 and 4
	    ieb = ieb24;
	    burst_time += time_offset;
	  }
	
	//-----------------------------------------
	//calculate range and azimuth resolution
	//-------------------------------------------
	StateVector sc_state("sc_state");
	target_frame_.ephemeris(sc_state,"Cassini",burst_time,"NONE");
	DirectionVector boresight("boresight",Fvector_[i],burst_time,0,0,1);

	//make a correction to beam's max gain
	Uvar azimuth=Bvector_[i].elevAndMaxGainToAzim(Uvar(0,"rad"));
	boresight.setAzimuthElevation(azimuth,Uvar(0,"rad"));
	
	//set ground range azimuth resolutions
	FloatVector velocity = sc_state.velocity();
	velocity.representIn(Fvector_[i]);
	DirectionVector velocity_dir = velocity;
	Uvar Vst_bore("along_track_speed");
	Uvar speed = velocity.magnitude();
	Vst_bore = velocity.magnitude();
	Vst_bore *= sqrt(1.0 - pow(dot(velocity_dir,boresight),2));
	
	TargetGeom tg(burst_time);
	tg.setState(sc_state);
	tg.setTarget(target_name_,target_frame_);
	tg.setLookDirection(boresight);
	Uvar range0 = tg.range();
	Uvar doppler0=tg.doppler(lambda);
	Uvar incidenceAngle = tg.incidenceAngle();
	int N_p = int(ieb.getPul());
	Uvar receive_window = double(N_p) *ieb.getPri();
	Uvar x_res = lambda * range0/(2.0 * Vst_bore * receive_window);
	Uvar r_res = speed_light/(2.0 * ieb.getBr());
	Uvar rg_res = r_res/sin(incidenceAngle);//ground resolution  
	
	//----------------
	//pulse gate
	//---------------
	Uvar azimuth_range = range0 
	  * Bvector_[i].getAzimuthWidthOneWay().getInUnits("rad");
	Uvar elev_range = range0 
	  * Bvector_[i].getElevationWidthOneWay().getInUnits("rad");
	Uvar range_3db = elev_range * tan(incidenceAngle);
	Uvar tau_3db = 2.0 * range_3db/speed_light;
	Uvar pulse_gate =(ieb.getTaup()+tau_3db)*speed_light/2.0;
	//cout<<"pulse gate "<<pulse_gate<<endl;
	//cout<<"prf "<<1.0/ieb.getPri()<<endl;
	//cout<<"pbw "<<pbw_ratio_/ieb.getPri()<<endl;
	
	amb[i].setTime(burst_time);
	amb[i].config(cfg_);//transfer grid size
	amb[i].setTarget(target_name_,target_frame_);
	amb[i].setBeam(Bvector_[i]);
	amb[i].setProcessWindow(range0,doppler0,
				1.0/ieb.getPri(),pbw_ratio_/ieb.getPri(),
				pulse_gate);
	amb[i].setTrackFrame(track_frame);
	amb[i].calAmbGeometry();//obtain range, doppler, area, gain
	
	
	//-----------------------
	//radar equation prefactor
	//-----------------------
	double Nr = (ieb.getTaup()).getInUnits("s") 
	  *(ieb.getBr()).getInUnits("Hz"); 
	double Ni = double(N_p)*Nr;
	Uvar X0 = Pt_ * Ni/pow(4.0*pi,3);	
	Uvar Pn = boltzmann_constant * Tsys_ * ieb.getRcv();
	//Uvar usable_area = amb[i].calUsableArea(X0,Pn,x_res,rg_res);
	//cout<<"usbel area"<<usable_area<<endl;
	//cout<<"noise "<<Pn<<endl;
	unsigned int Nbeam = ieb.getNumberOfBeams();
	ML_area[i] 
	  = amb[i].calMultilookUsableArea(X0,Pn,x_res,rg_res,double(Nbeam)*ieb.getBpd(),velocity,N_p);
	
	nlooks[i] = amb[i].getNumberofLooks();
      //cout<<"beam number ml usable"<<i+1<<" "<<ML_usable_area<<endl;	
      //cout<<"looks "<<looks<<endl;
      }//loop over all 5 beams
    
    //-------------------------------------------------------------------
    //once we calculated usable area, we can calculate usable crosstrack
    // extent of usable of "center" of the track: strip
    // Data Structure
    //  thermal snr, amb ratio, beam_gain: square matrix 
    //  let's take the middle value: dop = dop_bore, 
    //                               range- inner most to outer most
    //-----------------------------------------------------------------
    unsigned int Nrange_bin = amb[0].getGridsize();
    unsigned int Ndop_bin= Nrange_bin;//square grid
    Umat cross_track_strip("cross track extent",5,Nrange_bin+1);
    Umat thermal_snr_strip("thermal snr strip",5,Nrange_bin);
    Umat amb_ratio_strip("amb strip",5,Nrange_bin);
    Umat beam_gain_strip("beam gain strip",5,Nrange_bin);
    
    for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
    for (unsigned int i = 0; i < Nrange_bin+1;++i)
      { 
	cross_track_strip(i_beam,i) =(amb[i_beam].crosstrack)(i,Ndop_bin/2);
      }
    for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
    for (unsigned int i = 0; i < Nrange_bin;++i)
      { 
      thermal_snr_strip(i_beam,i)=(amb[i_beam].thermal_snr_ML)(i,Ndop_bin/2);
      amb_ratio_strip(i_beam,i) = (amb[i_beam].amb_ratio_ML)(i,Ndop_bin/2);
      beam_gain_strip(i_beam,i)=(amb[i_beam].oneway_beam_gaindB)(i,Ndop_bin/2);
      }
    
    Uvar min_cross, max_cross;
    min_max(min_cross,max_cross,cross_track_strip);
    //cout<<"strip processor "<< min_cross<<" "<<max_cross<<endl;
    
    //---------------------------------------------------
    //The distance between max_cross and min_cross: max_cross - min_cross
    // number of grid between them: max_cross - min_cross + 1
    //----------------------------------------------------
    unsigned int Ncross_bin =(unsigned int)(max_cross- min_cross).getInUnits("km")+1;
    vector<int> good_area_counter;
    good_area_counter.clear();

    for (unsigned int i = 0; i < Ncross_bin;++i)
      {
	Uvar cross_track_value = min_cross 
	  + (max_cross-min_cross)*double(i)/double(Ncross_bin-1);
	
	//------------------------------------------------------
	//check the above value satisfies usable area conditions
	//  in terms of snr, amb,and beam gain
	//-----------------------------------------------------
	list<Uvar> amb_collector, snr_collector,beam_collector;
	for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
	for (unsigned int j = 0; j < Nrange_bin;++j)
	  { 
	   
	    if (cross_track_strip(i_beam,j) <= cross_track_strip(i_beam,j+1))
	      {//right looking geometry, crosstrack increases with j
		if( cross_track_value >= cross_track_strip(i_beam,j) &&
		    cross_track_value < cross_track_strip(i_beam,j+1))
		  {
		    amb_collector.push_back(amb_ratio_strip(i_beam,j));
		    snr_collector.push_back(thermal_snr_strip(i_beam,j));
		    beam_collector.push_back(beam_gain_strip(i_beam,j));
		  }
	      }
	    else 
	      {//left looking geometry,crosstrack decreases with j
		if( cross_track_value <= cross_track_strip(i_beam,j) &&
		    cross_track_value > cross_track_strip(i_beam,j+1))
		  {
		    amb_collector.push_back(amb_ratio_strip(i_beam,j));
		    snr_collector.push_back(thermal_snr_strip(i_beam,j));
		    beam_collector.push_back(beam_gain_strip(i_beam,j));
		  }
	      }
	  }//loop over beam and range bins
	

	//snr is here defined as -(noise_equivalent_sigma0)
	if (amb_collector.size() == 0) 
	  {
	    extent += Uvar(0,"km");
	  }
	else if (amb_collector.size() == 1)
	  {
	    if (amb_collector.front() >= amb_ratio_ &&
		snr_collector.front() >= -noise_equiv_sigma0_ &&
		beam_collector.front() >=min_beam_gain_)
	      {
		extent+= Uvar(1,"km");//usable crosstrack
		good_area_counter.push_back(int(i));//good area
	      }
	  }
	else
	  {
	    amb_collector.sort();
	    snr_collector.sort();
	    beam_collector.sort();
	    if ( amb_collector.back() >= amb_ratio_ &&
		 snr_collector.back() >= -noise_equiv_sigma0_ &&
		 beam_collector.back() >=min_beam_gain_) 
	      { 
		extent +=Uvar(1,"km");//usable crosstrack
		good_area_counter.push_back(int(i));//good area
	      }
	  }
      }//loop over next range bin
  
    // cout<<"track extent"<<extent<<endl;
    //find gap across good area
    if (good_area_counter.size() == 0 || good_area_counter.size() == 1)
      {
	gap = Uvar(0,"km");
      }
    else
      {
	for (vector<int>::iterator pr=good_area_counter.begin();
	     pr <(good_area_counter.end()-1);++pr)
	  {
	    int a = *(pr+1);
	    int b = *pr;
	    if ((a-b) != 1)
	      {//there is a gap, like 2,3,4,6 where 5 is missing 
		int c = a - b - 1;
		gap += Uvar(c ,"km");
	      }
	  }
      }
    //cout<<"gap "<<gap<<endl;    
  }


//------------------------------------
//compute Azimuth resolution
//------------------------------------
void IebList::computeAzimuthRangeResolution(Uvar& x_res,
					    Uvar& r_res,
					    Uvar& rg_res,
					    const Time& t)
  {
    //-----------------------------------------
    //calculate range and azimuth resolution
    //-------------------------------------------
    static bool zero_warning = false;
    Ieb ieb;
    ieb=getIeb(t);

    StateVector sc_state("sc_state");
    target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
    DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//beam 3
	
    //set ground range azimuth resolutions
    FloatVector velocity = sc_state.velocity();
    velocity.representIn(Fvector_[2]);//beam 3
    DirectionVector velocity_dir = velocity;
    Uvar Vst_bore("along_track_speed");
    Uvar speed = velocity.magnitude();
    Vst_bore = velocity.magnitude();
    Vst_bore *= sqrt(1.0 - pow(dot(velocity_dir,boresight),2));
	
    TargetGeom tg(t);
    tg.setState(sc_state);
    tg.setTarget(target_name_,target_frame_);
    tg.setLookDirection(boresight);
    Uvar range0 = tg.range();
    Uvar incidenceAngle = tg.incidenceAngle();

    int N_p = int(ieb.getPul());
    Uvar receive_window = double(N_p) *ieb.getPri();
    Uvar lambda = speed_light/carrier_frequency;
    if (receive_window == 0)
      {
      x_res = 0;
      if (!zero_warning)
        {
        cout << "Warning: can't compute azimuth resolution with zero rcv window"
          << endl;
        cout << "         Azimuth resolution set to zero if no rcv window"
          << endl;
        zero_warning = true;
        }
      }
    else
      {
      x_res = lambda* range0/(2.0 * Vst_bore * receive_window);
      }
    r_res = speed_light/(2.0 * ieb.getBr());
    if (incidenceAngle == 0)
      {
      rg_res = 0*r_res;
      }
    else
      {
      rg_res = r_res/sin(incidenceAngle);//ground resolution      
      }
  }

Uvar IebList:: findToneChirpStartFrequency(const Uvar& pri,const Uvar& adc)
  {

    Uvar csf=Uvar(0,"Hz");
    Uvar bore_doppler=doppler_3dB_bore_(2,4);//beam3 boresight

    //-------------------------
    //peak position window: between -30 and -95 KHz
    //--------------------------
    
    Uvar min_return_freq = min_frequency_of_returned_cw_echo_;
    Uvar max_return_freq = max_frequency_of_returned_cw_echo_;
    Uvar chirp_step_freq = Uvar(30e6,"Hz")/pow(2,16);
  
    unsigned int pri_number = (unsigned int) round_double( (adc*pri).getInUnits(""));
    pri_number /=2;
    pri_number *=2;
    //cout<<"pri number "<< pri_number<<endl;
   
    Uvar pri_choice = (1/adc)*pri_number;
    Uvar prf_choice= 1/pri_choice;

    Uvar tmp_csf;
   
    double tmp_csf_pri;
    unsigned int min_csf_number, max_csf_number;
    Uvar min_csf, max_csf;
    Uvar tmp_return;
    min_csf = slo_frequency - bore_doppler + min_return_freq;
    max_csf = slo_frequency - bore_doppler + max_return_freq;
    min_csf_number = (unsigned int) round_double( (min_csf /Uvar(30e6,"Hz")*pow(2,16)).getInUnits(""));
    max_csf_number = (unsigned int) round_double( (max_csf /Uvar(30e6,"Hz")*pow(2,16)).getInUnits(""));
    if(max_csf_number < min_csf_number) ErrorMessage("IebList::findToneChirpStartFrequency: max index is smaller than min index").throwMe();
    
    //cout<<"min and max csf number "<< min_csf_number<<" "<<max_csf_number<<endl;
    //cout<<"bore doppler "<< bore_doppler<<endl;
   

    unsigned int Ncsf_step = max_csf_number - min_csf_number+1;
    vector<Uvar> local_return_peak;
    vector<Uvar> local_remainder;

    local_return_peak.clear();
    local_remainder.clear();
    Dvec remainder("",Ncsf_step);
    Uvec return_peak_position("",Ncsf_step);
    for(unsigned int ii=0;ii<Ncsf_step;++ii){
      tmp_csf= Uvar(30e6,"Hz")/pow(2,16)*double(ii+min_csf_number);
      tmp_csf_pri =( tmp_csf *pri_choice).getInUnits("");
      remainder(ii)= tmp_csf_pri - int(tmp_csf_pri);
      if(remainder(ii)>0.5) remainder(ii)=1.0-remainder(ii);
      tmp_return= -slo_frequency + tmp_csf + bore_doppler;
      if( (tmp_return >Uvar(-101000,"Hz") && tmp_return<Uvar(-99000,"Hz")) ||
	  (tmp_return>Uvar(-51000,"Hz") && tmp_return<Uvar(-49000,"Hz")))
	remainder(ii)=0.5;
      return_peak_position(ii)= tmp_return;
      local_return_peak.push_back(tmp_return);
      local_remainder.push_back(remainder(ii));
    }

    //debug
    //Plot a;
    //a.addXY(local_return_peak,"KHz",local_remainder,"",line("none"),sym("circle","red",1));
    //a.setTitle("peak position vs remainder");
    //a.show("x");


    //find min remainder
    unsigned int min_index;
    remainder.min(min_index);
    //---------------------------------------------------------
    //find csf numbers whose return peak is close to -adc/4
    //cout<<"min remainder "<< remainder(min_index)<<endl;
    //if more than one minimum points, find one close to middle
    //------------------------------------------------------------
    double remainder_cutoff=remainder(min_index);
    unsigned int csf_number;
    
    vector<unsigned int> local_csf_number; 
    vector<Uvar> distance_from_middle;
    for(unsigned int ii=0;ii<Ncsf_step;++ii){
      if(remainder(ii)<= remainder_cutoff){
	local_csf_number.push_back(ii+min_csf_number);
	distance_from_middle.push_back(fabs(return_peak_position(ii)+adc/4.0));
      }
    }
    
    
    Uvar min_distance, max_distance;
    min_max(min_distance, max_distance, distance_from_middle);
    unsigned int NN=distance_from_middle.size();
    //cout<<"number of remainders below cutoff "<< NN<<endl;
    unsigned int MM;
    for(unsigned int ii=0;ii<NN;++ii){
      if(min_distance==distance_from_middle[ii]){
	MM=ii;
	break;
      }
    }
    csf_number =local_csf_number[MM] ;
    csf =Uvar(30e6,"Hz")/pow(2,16)*csf_number;
    tmp_return= -slo_frequency + csf + bore_doppler;
  
    //cout<<"min and max return "<< min_return_freq<<" "<<max_return_freq<<endl;
    //cout<<"doppler "<< doppler_3dB_bore_(2,4)<<endl;
    //cout<<"best value remainder "<< remainder(min_index)<<endl;
    //cout<<"optimum value "<< (csf*pri_choice).getInUnits("")-int( (csf*pri_choice).getInUnits("") )<<endl;
    //cout<<"pri number "<< pri_number<<endl;
    //cout<<"csf number "<< csf_number<<endl;
    //cout<<"return freq "<< tmp_return<<endl;
    //if(remainder(min_index)!=0.0) ErrorMessage("IebList:: could not find a chirp start frequency satisfying chir_start_freq*pri=N ").throwMe();
    return(csf); 
  }
//---------------------------------
//computeLooks(const Time& t) of beam 3
//----------------------------------
//unsigned int IebList::computeLooks(const Time& t)
//{  
//unsigned int looks = 0;
//Ieb ieb;
//ieb = getIeb(t);
//StateVector sc_state;
//target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
//FloatVector velocity=sc_state.velocity();
//Uvar lambda = speed_light/carrier_frequency;
//DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//use beam 3
//
//TargetGeom tg(t);
///tg.setState(sc_state);
//tg.setTarget(target_name_,target_frame_);
// tg.setLookDirection(boresight);
//Uvar doppler=tg.doppler(lambda);
//Uvar range = tg.range();
//double cos_theta = (doppler * lambda/
//	      (2.0*velocity.magnitude())).getInUnits("");
//unsigned int Nbeam = ieb.getNumberOfBeams();
//Uvar freq_bp = 2.0 * velocity.magnitude()*velocity.magnitude()
//  *double(Nbeam)*ieb.getBpd()/(lambda*range);//5 beam * bpd
//freq_bp *= (1.0-cos_theta*cos_theta);
//Uvar process_window = (1.0/ieb.getPri()) * 0.8;//80% of prf
//looks= (unsigned int) round_double((process_window/freq_bp).getInUnits(""));
//if (looks < 1) looks = 1;
//return(looks);
//}
//
//-------------------------------------
//compute looks for using a given ieb
//---------------------------------------
//unsigned int IebList::computeLooks(const Time& t, const Ieb& ieb)
//{
//Ieb ieb0 = ieb;
//unsigned int looks = 0;
//StateVector sc_state;
//target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
//FloatVector velocity=sc_state.velocity();
//Uvar lambda = speed_light/carrier_frequency;
//DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//use beam 3
//TargetGeom tg(t);
//tg.setState(sc_state);
//tg.setTarget(target_name_,target_frame_);
//tg.setLookDirection(boresight);
//Uvar doppler=tg.doppler(lambda);
//Uvar range = tg.range();
//double cos_theta = (doppler * lambda/
//	      (2.0*velocity.magnitude())).getInUnits("");
//unsigned int Nbeams = ieb0.getNumberOfBeams();
//Uvar freq_bp = 2.0 * velocity.magnitude()*velocity.magnitude()
//  *double(Nbeams)*ieb0.getBpd()/(lambda*range);//5 beam * bpd
//freq_bp *= (1.0-cos_theta*cos_theta);
//Uvar process_window = (1.0/ieb0.getPri()) * 0.8;//80% of prf
//looks= (unsigned int) round_double((process_window/freq_bp).getInUnits(""));
//if (looks < 1) looks = 1;
//return(looks);
//}

//---------------------------------
//computeBpd for a given number of looks
//----------------------------------
//Uvar IebList::computeBpd(const Time& t, const unsigned int& looks)
//{
//Ieb ieb;
//ieb=getIeb(t);
//Uvar process_window = (1.0/ieb.getPri()) * 0.8;//80% of prf 
//
//StateVector sc_state;
//target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
//FloatVector velocity=sc_state.velocity();
//Uvar lambda = speed_light/carrier_frequency;
//DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//use beam 3
//TargetGeom tg(t);
//tg.setState(sc_state);
//tg.setTarget(target_name_,target_frame_);
//tg.setLookDirection(boresight);
//Uvar doppler=tg.doppler(lambda);
//Uvar range = tg.range();
//double cos_theta = (doppler * lambda/
//	      (2.0*velocity.magnitude())).getInUnits("");
//
////freq_bp = process-window /looks
//unsigned int Nbeams = ieb.getNumberOfBeams();
//Uvar bpd = process_window/double(looks)*(lambda*range)
//  /(2.0* velocity.magnitude()*velocity.magnitude()
//    *double(Nbeams)*(1.0-cos_theta*cos_theta));
//bpd -= Uvar(1.1*mstos,"s");//to avoid possible problems associated with quantization error
//// Uvar bpd = Uvar(bpd_ms * mstos,"s");
//
///  Uvar freq_bp = 2.0 * velocity.magnitude()*velocity.magnitude()
//*double(ieb.getNumberOfBeams())*ieb.getBpd()/(lambda*range);//5 beam * bpd
//freq_bp *= (1.0-cos_theta*cos_theta);
//looks= (unsigned int) (process_window/freq_bp).getInUnits("");
//return(bpd);
//}

//-----------------------------
//readAllIebRecords(): read all the record in the file
//-----------------------------
void IebList::readAllIebRecords(const string& filename, const string& filetype) 
  {
  
  if (filetype=="wb"|| filetype=="w" )
    {
      ErrorMessage("Can't read from outputfile").throwMe();
    }
  cout<<"clearing internal ieb list before reading in a new sequence"<<endl;
  ieb_list_.clear();
  
  FileMgr file(filename,filetype);
  while(!file.eof())
    {
      double et_time_in_sec;
      file.read(et_time_in_sec);
      
      Ieb ieb(Uvar(et_time_in_sec,"s"));
      
      unsigned int iebMode;
      file.read(iebMode);
      ieb.setIebMod(iebMode);
      
      for(unsigned int i_slow=0; i_slow < Nslowfield;++i_slow){
	file.read(slowfield_[i_slow]);}
            
      for(unsigned int i_fast=0;i_fast<Nfastfield;++i_fast){
	file.read(fastfield_[i_fast]);}

      for(unsigned int i_power = 0; i_power <Npowerfield;++i_power){
	file.read(powerfield_[i_power]);}
      
      for(unsigned int i_tnc =0; i_tnc< Ntncfield;++i_tnc){
	file.read(tncfield_[i_tnc]);}

      ieb.loadSlowFastfield(slowfield_,fastfield_);  
      ieb.loadPowerfield(powerfield_);
      ieb.loadTncfield(tncfield_);
      ieb.decodeSlowFastfield();
      ieb.decodePowerTncfield();
      ieb_list_.push_back(ieb);          
    }
 
  //check whether time information is stored from past to future order
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
       ieb_pt < (ieb_list_.end()-1);++ieb_pt)
    {
      if ( (*ieb_pt).getTime() > (*(ieb_pt+1)).getTime()) 	
	ErrorMessage("IebList: time is not increasing with record number").throwMe();
    }
  
  cout<<"read all the records"<<endl;
  start_time_ = ieb_list_.begin()->getTime();
  end_time_ = (ieb_list_.end()-1)->getTime();
  file.close();
 
  }

//-----------------------------------------------------------
//writeRecord(): write time ordered fast/slow field into file
//-----------------------------------------------------------
void IebList::writeAllIebRecords(const string& filename, const string& filetype) 
  {
  if(ieb_list_.size()==0) return;//no record, do nothing
  if (filetype=="rb"|| filetype=="r" )
    {
      throw ErrorMessage("Can't write to inputfile");
    }
  
  FileMgr file(filename,filetype);
 
  //check time order
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin();
       ieb_pt < (ieb_list_.end()-1);++ieb_pt)
    {
     
      if ((*ieb_pt).getTime() > (*(ieb_pt+1)).getTime())
	{
	cout<<"time of one record before "<< ieb_pt->getTime().utc("ISOD")<<endl;
	cout<<"time of current record "<< (ieb_pt+1)->getTime().utc("ISOD")<<endl;
	ErrorMessage("IebList: time is not increasing with record number").throwMe();
	}
    }
 
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
       ieb_pt < ieb_list_.end() ;
       ++ieb_pt)
    {
      // gh: add some debug print statements
      //cout << "IebList(binwrite1)(time/TFI/Mod): " << ((ieb_pt->getTime()).et()).getInUnits("s") <<":(TFI:"<< ieb_pt->getFastTfi() <<"):"<< ieb_pt->getIebMod()<<"FIN"<< ieb_pt->getFin() <<endl;

      //check encoding complete
      //if(!ieb_pt->encodeSlowFastfieldComplete()) ieb_pt->encodeSlowFastfield();
      //if(!ieb_pt->encodePowerTncfieldComplete()) ieb_pt->encodePowerTncfield();
      // encode regardless to insure match 12/19/2003
      try{
	ieb_pt->encodeSlowFastfield();
	ieb_pt->encodePowerTncfield();
      }
      catch(ErrorMessage& e){
	cerr<<e.msg<<endl;
	cerr<<"encoding  failed "<<endl;
	cerr<<"this record will not be written into ieblist "<<endl;
	cerr<<"time information "<< ieb_pt->getTime().utc("ISOD")<<endl;
	continue;
      }
      try{
	ieb_pt->decodeSlowFastfield();
	ieb_pt->decodePowerTncfield();
      }
      catch(ErrorMessage& e){
	cerr<<e.msg<<endl;
	cerr<<"decocding failed "<<endl;
	cerr<<"this record will not be written into ieblist "<<endl;
	cerr<<"time information "<< ieb_pt->getTime().utc("ISOD")<<endl;
	continue;
      }
      
      //write time
      double time=((ieb_pt->getTime()).et()).getInUnits("s");
      file.write(time);
      

      //write IebMode
      unsigned int iebMod = ieb_pt->getIebMod();
      file.write(iebMod);

      slowfield_ = (*ieb_pt).getSlowfield();
      fastfield_ = (*ieb_pt).getFastfield();
      powerfield_ = (*ieb_pt).getPowerfield();
      tncfield_=(*ieb_pt).getTncfield();

      // gh: add some debug print statements
      //cout << "IebList(binwrite2)(time/TFI/Mod): " << ((ieb_pt->getTime()).et()).getInUnits("s") <<":(TFI:"<< ieb_pt->getFastTfi() <<"  MOD):"<< ieb_pt->getIebMod()<<"FIN"<< ieb_pt->getFin();


      for(unsigned int i_slow=0; i_slow < Nslowfield;++i_slow){
	file.write(slowfield_[i_slow]);
      }
      //cout << ".  SLOWfield0 = " << slowfield_[0] <<" :"<< slowfield_[1]  <<" :"<< slowfield_[2] <<" :HERE:->"<< ((slowfield_[2] & 0xff00)>>8) ;
	
      for(unsigned int i_fast=0;i_fast<Nfastfield;++i_fast){
	file.write(fastfield_[i_fast]);
      }
      //cout << "  FASTfield0 = " << fastfield_[0] <<" :"<< fastfield_[1]  <<" :HERE->"<< ((fastfield_[1] & 0x00ff)) << endl;

      for(unsigned int i_power = 0; i_power <Npowerfield;++i_power){
	file.write(powerfield_[i_power]);}
      
      for(unsigned int i_tnc =0; i_tnc< Ntncfield;++i_tnc){
	file.write(tncfield_[i_tnc]);}
    }
  file.close();
  } 

//-------------------------------------------
//{get} AllIebRecords: transfer internal ieb 
// records to target list
//--------------------------------------------
void IebList::getAllIebRecords(vector<Ieb>& target_ieb_list)
  {
    target_ieb_list.resize(ieb_list_.size());
    target_ieb_list = ieb_list_;
  }



//--------------------------- Not part of IebList---------------------


//------------------------
//GH:computeBII(): set BII, SIN, FIN
//------------------------
unsigned int IebList::computeBii()
  // return totla SIN,FIN count, by SIN * 1000 + FIN.
  // Therefore 200500. Equal 200 sins and 500 fins (max value.)
  {
    int SIN = 1;
    int FIN = 1 ;
    int fin_total = 0;
    int sin_total = 0;
    unsigned int new_bii = 0;
    unsigned short bem_mask ;
    unsigned int i_beam ;
    unsigned int beam_status ;
    unsigned int Nbeam ;
    Uvar end_inst_time = Uvar (0,"s");

    //cout << "in compute BII" << endl;

    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
       ieb_pt < (ieb_list_.end()-1);++ieb_pt)
    {
      
      if ( (*ieb_pt).getIebMod() == 4) {
	cout  << "ERROR: (IebList:computeBii): Found instruction with default IebMOD of 4 at " << (ieb_pt)->getFastTfi() <<":MOD:"<< (ieb_pt)->getIebMod() <<":Mode:"<< (ieb_pt)->getMod() << "BIIwas: " << (ieb_pt)->getBii() << "at FIN:" << (ieb_pt)->getFin() ;
      }

      //cout << "in computeBII at TFI: " <<  (ieb_pt)->getFastTfi() <<":MOD:"<< (ieb_pt)->getIebMod() <<":Mode:"<< (ieb_pt)->getMod() << "BIIwas: " << (ieb_pt)->getBii() << "at FIN:" << (ieb_pt)->getFin() ;

      // if power or TNC then
      //first set to old value, this way if power, tnc or Fast only, other value will be previous IEB
      // this help the binary file.  The ascii file treats these as don't cares
      //SIN/FIN does not increament, but set to old value for binary file
        (ieb_pt)->setFin(FIN) ;
	(ieb_pt)->setSin(SIN) ; 


      if ((ieb_pt)->getIebMod() == 2 || (ieb_pt)->getIebMod() == 3) {
        (ieb_pt)->setFin(FIN++) ;                                    // FIN alway increase (either Fast or Slow/Fast
	if ((*ieb_pt).getIebMod() == 3 ) (ieb_pt)->setSin(SIN++) ;   // Slow only increase on Slow/Fast

	// determine number of Beams
	bitset(bem_mask,0,4,ieb_pt->getBem());
	Nbeam = 0;
	for (i_beam = 0; i_beam < 5;++i_beam) {
	  beam_status=bitget(bem_mask,i_beam,i_beam);
	  if (beam_status == 1) Nbeam++;
	}
	if (Nbeam == 0) Nbeam = 1; // calibration mode may have beam = 0, treat is if beam = 1

	if (Nbeam < 0 || Nbeam > 5) ErrorMessage("IebList:computeBii:Beam counting error").throwMe();
	// give an additional second of pad prior to any power instruction
	if ((ieb_pt+1)->getIebMod() == 0) end_inst_time =  (ieb_pt+1)->getFastTfi() - Uvar(2,"s")  ;
	else end_inst_time =  (ieb_pt+1)->getFastTfi() ;

        new_bii =  (unsigned int) ( (unsigned int)((end_inst_time.getInUnits("s")-(ieb_pt)->getFastTfi().getInUnits("s"))*1000)  / ((unsigned int) ((ieb_pt)->getBpd().getInUnits("ms"))) )  ;

	// added 01/24/2006:  The -2 seecond pad above can cause the (unsigned int) substraction to 
	// go really large if end_time is less than Fast TFI.  And really the 2 second can be empty
	if (end_inst_time < (ieb_pt)->getFastTfi() ){
        	cout << "BII:  no time between these ILXes for more bursts." << endl;
        	new_bii = 0;
	}

	// adjust for number of beams
	new_bii = Nbeam*((int)(new_bii/Nbeam)) ;

	// case when BPD in greater that allocated time
	if (new_bii < 1 * Nbeam ) { // cannot get at least 1 cycle
	  cout << "INFO)IebList:computeBii-> Delete IEB.  Duration of IEB less than 1 burst(beam_cycle)!" << endl;
	  cout << "numbers>> " << new_bii <<")"<<(ieb_pt)->getFastTfi() << "::" << (ieb_pt+1)->getFastTfi() << ":B:" << Nbeam << ":Bpd:" << (ieb_pt)->getBpd() << "= "<<  ( (ieb_pt+1)->getFastTfi()- (ieb_pt)->getFastTfi() ) /  (ieb_pt)->getBpd() << " =" <<( (((ieb_pt+1)->getFastTfi().getInUnits("s")-(ieb_pt)->getFastTfi().getInUnits("s"))*1000)  /  ((ieb_pt)->getBpd().getInUnits("ms")) )  << endl ;
  
	  deleteIeb( *(ieb_pt) );
	} else if (new_bii > 255) {
	  cout << "Warning: BII going to  0: Bii was:" << new_bii <<" TFI:"<<(ieb_pt)->getFastTfi() << ":TFI+1:" << (ieb_pt+1)->getFastTfi() << ":Beams:" << Nbeam << ":Bpd:" << (ieb_pt)->getBpd() << "TFI-TFInext/Bpd "<<  ( (ieb_pt+1)->getFastTfi()- (ieb_pt)->getFastTfi() ) /  (ieb_pt)->getBpd() << " in sec:" <<  (unsigned int) ( ((ieb_pt+1)->getFastTfi().getInUnits("s")-(ieb_pt)->getFastTfi().getInUnits("s"))  /  (ieb_pt)->getBpd().getInUnits("s") ) << endl ;
	  new_bii = 0;
	  (*ieb_pt).setBii(new_bii);
	} else (*ieb_pt).setBii(new_bii);
	//cout << " BII now: " << (ieb_pt)->getBii() <<"FIN"<<(ieb_pt)->getFin()<< endl;
      }

      // this is the line that correct the sin/fin number from the user edit file.
      // exactly why I don't understand yet.... (gh: 12/18/2003)
      (ieb_pt)->encodeSlowFastfield();

      if (SIN == 256) {SIN = 0; sin_total = sin_total+256; } //  handle roll over
      if (FIN == 256) {FIN = 0; fin_total = fin_total+256; } // handle roll over

    }
  
      sin_total = sin_total + SIN -1 ;
      fin_total = fin_total + FIN -1 ;

      if (sin_total > 200) {cout << "ERROR: computeBII: Too many SINs: " << sin_total << endl ; }
      if (fin_total > 500) {cout << "ERROR: computeBII: Too many FINs: " << fin_total << endl ; }
      
      //cout << "1)---------SIN&FINS(bii)-------- " << sin_total << ":" << fin_total << endl;

      //cout << " compute BII.....DONE.." << endl;

     


      return(sin_total*1000+fin_total) ;
  }

//------------------------
//GH: compress slowfield.  If only different is waveform, merge waveform characteristics
//------------------------
unsigned int IebList::compressSlowfield(const int percentage)
  {
    float currentduty = 0.0;
    float nextduty = 0.0;
    unsigned int numbermerged = 0 ;
    int adc_multiplier = 1;
 
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
       ieb_pt < (ieb_list_.end()-1);++ieb_pt)
    {
      if ( (ieb_pt)->getIebMod() == 3 && (ieb_pt+1)->getIebMod() == 3) {

	cout << " in compressSlowfield:" << (ieb_pt)->getSlowTfi() ;
	
	if ( (ieb_pt)->getAdcEncodedValue() == 3 ) adc_multiplier = 10;
	     else adc_multiplier = 1; 

	currentduty = ( (ieb_pt)->getCsd().getInUnits("ns")+0.5 ) / 1000 * (ieb_pt)->getCsq() / (ieb_pt)->getPriEncodedValue() / (ieb_pt)->getAdc().getInUnits("MHz") * adc_multiplier ;

	nextduty = ( (ieb_pt+1)->getCsd().getInUnits("ns")+0.5 ) / 1000 * (ieb_pt+1)->getCsq() / (ieb_pt+1)->getPriEncodedValue() / (ieb_pt+1)->getAdc().getInUnits("MHz") * adc_multiplier ;

	//nextduty    = (ieb_pt+1)->getCsd().getInUnits("ns") * 133.333+0.5;
	//chirp_width = (int)((dptr->rpb.radar.sf_chirp_step_dur*133.333+0.5) / 1000 * dptr->rpb.radar.sf_chirp_step_qty) ;
	//XMTR = (float)chirp_width / dptr->rpb.radar.ff_pulse_rep_interval / adc_ten * ADC_mhz  ;

	
	  cout << " current DC = :" << currentduty <<" next DC = " << nextduty << endl << endl;

      } // mode IF

    } // vector IEBLIST for statement


    return(numbermerged) ;
  }

//------------------------
//write_wtk_Record(): write time ordered fast/slow field into file (in wtk format)
//------------------------
void IebList::writeAllwtkIebRecords(const string& filename, const int base)
  {
   // ---------------- header for wtk file ----------------

        char slowheader[] = "\ntype\tSTFI\tSTYP\tDTN\tSIN\tMOD\tCSR\tADC\tRCV\tTRO\tBAQ\tBEM\tAT1\tAT3\tAT4\tRIP\tRAD\tCSD\tCSQ\tCFS" ;

        char fastheader[] = "\tFTFI\tFTYP\tFIN\tBII\tPUL\tBPD\tRWD\tPRI\tCSF";

//        char pwrheader[] = "\tPTFI\tPTYP\tPMD";

//        char tncheader [] = "\tTTFI\tTTYP\tTNC\tTCA\tTCB\tTCC";

        char edsheader [] = "\tRECF\tDEL1\tDEL2\tDELB\tSCL1\tSCL2\tNCOF\tUPCN\tUPCA\tSIN\tFIN\tEDSM\tATSM\tATSA\tRATN";

 //       char engheader [] = "\tVALIDTYPE\tAUTO_RAD_ENABLE\tAUTO_GAIN_ENABLE\tCTU_TAIL_ENABLE\tENG_TAIL_ENABLE\tPICK_UP_MODE\tND_IP\tRL_IP\tAUTO_ND_MAX\tAUTO_RL_MAX\tPRODUCTION_RATE\tSDB_WORDS\tSABS_PRODUCED\tSAB_SIZE\tCUM_DATA_VOLUME";

//        char trajheader [] = "\tTARGET_RANGE\tTARGET_ANGLE";
        char tab [] = "\t";
        int tab_loop ;

// ------------------------- end of header ----------------------------


  string filetype = "w";

  ofstream file (filename.c_str());

// ----- write header on file
file << slowheader << fastheader ;
//<< pwrheader << tncheader ;
file << edsheader ;
//<< engheader << trajheader ;
file << endl ;

// Target Geom stuff
/*
 bool geometry_print = false;
 StateVector sc_state;
 Uvar lambda = speed_light/carrier_frequency;
*/

//StateVector sc;
//Uvar lambda = speed_light/carrier_frequency;  
    //Frame beam_frame ("CASSINI_RADAR_3", "Cassini");
    // Frame target_frame ("IAU_"+target_name_, target_name_); 

/*
  StateVector sc_state;
  target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
  FloatVector velocity=sc_state.velocity();
  Uvar lambda = speed_light/carrier_frequency;
  DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//use beam 3
  TargetGeom tg(t);
  tg.setState(sc_state);
  tg.setTarget(target_name_,target_frame_);
  tg.setLookDirection(boresight);
  Uvar doppler=tg.doppler(lambda);
  Uvar range = tg.range();
*/


 if (base == 16) file << std::hex ;       // write numbers as hex to wtk file

  //check time order
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin();
       ieb_pt < (ieb_list_.end()-1);++ieb_pt)
    {
        if ((*ieb_pt).getTime() > (*(ieb_pt+1)).getTime()) ErrorMessage("IebList(wtk): time is not ordered from past to future").throwMe();
    }

  

  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); 
       ieb_pt < ieb_list_.end() ;
       ++ieb_pt)
  {

    //  geometry_print = false; // default.  Slow or Slow/Fast should set to true
  int ilx_type = (*ieb_pt).getIebMod()  ;
  /*
  if (ilx_type == 2) //fast only
    cout << "wtkwriteFAST: " << (int)(*ieb_pt).getFastTfi().getInUnits("s") << ":FIN:" <<(int)(*ieb_pt).getFin()  << endl;
  if (ilx_type == 3) //S/F
    cout << "wtkwrite.S/F: " << (int)(*ieb_pt).getSlowTfi().getInUnits("s") <<":" << (int)(*ieb_pt).getFastTfi().getInUnits("s")<< ":SIN:" <<(int)(*ieb_pt).getSin()  << ":FIN:" <<(int)(*ieb_pt).getFin()<< endl;
  */

switch ( ilx_type ) {
case (0):				// -- POWER PARAMETERS --------
        //cout << "POWER ILX" << endl;
        file << "@" ;
        file << tab << (int)(*ieb_pt).getPowerTfi().getInUnits("s") ;
        file << tab << (*ieb_pt).getPowerTyp() ;
        file << tab << (*ieb_pt).getPmd() ;

        for (tab_loop = 0 ; tab_loop < (19-3) ; tab_loop++ ) { //slowfield tabs
	  if (base == 16) file << tab ;
	  else file << tab << 0;
        }

        for (tab_loop = 0 ; tab_loop < 9 ; tab_loop++ ) { //fastfield tabs
          if (base == 16) file << tab ;
	  else file << tab << 0;
 
        }
        break;

case (1):				// -- TNC PARAMETERS --------
  // cout << "TNC ILX" << endl;
        file << "^" ;
//      write TNC - hardcoded value for now.....
	// file << tab <<  (*ieb_pt).getTncTfi().getInUnits("s") << tab << "1" << tab<< "40" << tab<< "54" << tab<< "A9" << tab<< "54" ;
	 file << tab << (int)(*ieb_pt).getTncTfi().getInUnits("s")  << tab << (*ieb_pt).getTncTyp() << tab << (*ieb_pt).getTnc() ;
	 file << tab << (*ieb_pt).getTca() << tab << (*ieb_pt).getTcb() << tab << (*ieb_pt).getTcc() ;
          
        for (tab_loop =0; tab_loop < (19-6); tab_loop++){ //slowfield tabs
          if (base == 16) file << tab ;
          else file << tab << 0;
        }
        for (tab_loop =0; tab_loop < 9; tab_loop++){ //fastfield tabs
          if (base == 16) file << tab ;
          else file << tab << 0;
        }
        break;
case (2):				//-- FAST FIELD PARAMETERS --------
       // cout << "Fast Only ILX" << endl;
// -- SLOW FIELD PARAMETERS need to get 'tabbed' out --------
        file << "<"  ;
        for (tab_loop =0; tab_loop < 19; tab_loop++){ //slowfield tabs
          if (base == 16) file << tab ;
	  else file << tab << 0;
        }

        file << tab << (int)(*ieb_pt).getFastTfi().getInUnits("s") ;
        file << tab << (*ieb_pt).getFastTyp() ;
        file << tab << (*ieb_pt).getFin() ;
        file << tab << (*ieb_pt).getBii() ;
        file << tab << (*ieb_pt).getPul() ;
        file << tab << (*ieb_pt).getBpdEncodedValue() ;
        file << tab << (*ieb_pt).getRwdInPriUnits() ;
        file << tab << (*ieb_pt).getPriEncodedValue() ;
        file << tab << (*ieb_pt).getChirpStartFrequencyEncodedValue() ;

	//geometry_print = true ; // print geometery information at end of line



	//Time t = (*ieb_pt).getTime();
	//target_frame_.ephemeris(sc,"Cassini",(*ieb_pt).getTime(),"NONE");
    //DirectionVector boresight("boresight",Fvector_[2],(*ieb_pt).getTime(),0,0,1);//beam 3
 
    //DirectionVector boresight("NONE",beam_frame,(*ieb_pt).getTime() ,0,0,1); 
    //  //  TargetGeom tg ( (*ieb_pt).getTime() );
    // tg.setState(sc);
    //  tg.setTarget(target_name_,target_frame_);
    //  tg.setLookDirection(boresight);
    //cout <<"WTK:TARGETGEOM: Alt = " << tg.altitude() ;
    //cout << endl;

        break;

case (3):				// -- SLOW/FAST FIELD PARAMETERS --------
       // cout << "SlowFast ILX" << endl ;
// -- SLOW FIELD PARAMETERS --------
        file << ">" ;
        file << tab << (int)(*ieb_pt).getSlowTfi().getInUnits("s") ;
        file << tab << (*ieb_pt).getSlowTyp() ;
        file << tab << (*ieb_pt).getDtn() ;
        file << tab << (*ieb_pt).getSin() ;
        file << tab << (*ieb_pt).getMod() ;
        file << tab << (*ieb_pt).getCsr() ;
        file << tab << (*ieb_pt).getAdcEncodedValue() ;
        file << tab << (*ieb_pt).getRcvEncodedValue() ;
        file << tab << (*ieb_pt).getTroInPriUnits() ;
        file << tab << (*ieb_pt).getBaq() ;
        file << tab << (*ieb_pt).getBem() ;
        file << tab << ( ((*ieb_pt).getAt1N1dB()*0x80) + ((*ieb_pt).getAt1N2dB()*0x04) + ((*ieb_pt).getAt1N3dB()/4) ) ;
        file << tab << ( ((*ieb_pt).getAt3N1dB()*0x80) + ((*ieb_pt).getAt3N2dB()*0x04) + ((*ieb_pt).getAt3N3dB()/4) ) ;
        file << tab << ( ((*ieb_pt).getAt4N1dB()*0x80) + ((*ieb_pt).getAt4N2dB()*0x04) + ((*ieb_pt).getAt4N3dB()/4) ) ;
        file << tab << (*ieb_pt).getRipEncodedValue() ;
        file << tab << (*ieb_pt).getRad();
        file << tab << (*ieb_pt).getCsdEncodedValue() ;
        file << tab << (*ieb_pt).getCsq() ;
        file << tab << (*ieb_pt).getCfsEncodedValue();
// -- FAST FIELD PARAMETERS --------
        file << tab << (int)(*ieb_pt).getFastTfi().getInUnits("s") ;
        file << tab << (*ieb_pt).getFastTyp() ;
        file << tab << (*ieb_pt).getFin() ;
        file << tab << (*ieb_pt).getBii() ;
        file << tab << (*ieb_pt).getPul() ;
        file << tab << (*ieb_pt).getBpdEncodedValue() ;
        file << tab << (*ieb_pt).getRwdInPriUnits() ;
        file << tab << (*ieb_pt).getPriEncodedValue() ;
        file << tab << (*ieb_pt).getChirpStartFrequencyEncodedValue() ; 

	//geometry_print = true ; // print geometery information at end of line

        break;

default:
  cout << " ERROR with ILXTYPE (IebMod) > 3 at time :" << (*ieb_pt).getSlowTfi().getInUnits("s") << endl;

} // end of switch


// -- EDS (Echo Delay Simulator) PARAMETERS --------
 int tempint = 0;
 float tempcsf0 = 0;

 if ((*ieb_pt).getIebMod() > 1 ) { // only do for slow fast field items (not power and tnc)
 if ((*ieb_pt).getMod() == 1 || (*ieb_pt).getMod() == 9) tempint = 1; else tempint = 0;
 file << tab << tempint ;             // this is RECF

 tempint = 0;
 file << tab << tempint ;             // this is DEL1
 file << tab << tempint ;             // this is DEL2
 file << tab << tempint ;             // this is DELB
 tempint = 15;
 file << tab << tempint ;             // this is SCL1
 file << tab << tempint ;             // this is SCL2
 tempint = 0;
tempint = (*ieb_pt).getRcvEncodedValue() ;
 switch (tempint){
 case (0):
   tempcsf0 = 9884.375;// in khz
   break;
 case (1):
   tempcsf0 = 9537.5;
   break;
 case (2):
   tempcsf0 = 9075.0;
   break;
 case (3):
   tempcsf0 = 5375.0;
   break;
 }// end of NCOF switch
 tempint = int(((tempcsf0-(*ieb_pt).getChirpStartFrequency().getInUnits("KHz"))*1000+105000114)/1029.97-98305+0.5);
 file << tab << tempint ;             // this is NCOF


 tempint = 0;
 file << tab << tempint ;             // this is UPCN
 file << tab << tempint ;             // this is UPCA

 file << tab << (*ieb_pt).getSin() ;             // this is SIN
 file << tab << (*ieb_pt).getFin() ;             // this is FIN

 tempint = 3;
 file << tab << tempint ;             // this is EDSM

 tempint = 1;
 file << tab << tempint ;             // this is ATSM

 tempint = 0;
 file << tab << tempint ;             // this is ATSA
 file << tab << tempint ;             // this is RATN
 // end of EDS stuff
 /*
// Target Geometry printing
 if (geometry_print){  // should be true is Slow/Fast or Fast only case
   Time t = (*ieb_pt).getTime();
   target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
   DirectionVector boresight("boresight",Fvector_[2],t,0,0,1);//use beam 3
   TargetGeom tg(t);
   tg.setState(sc_state);
   tg.setTarget(target_name_,target_frame_);
   tg.setLookDirection(boresight);
   Uvar doppler=tg.doppler(lambda);
   Uvar range = tg.range();
   cout << "WTK:TARGET Altitude = " << tg.altitude();
   cout << " Range = " << tg.range();
   cout << " Doppler = " << tg.doppler(lambda);
   cout << " IncidenceAngle = " << tg.incidenceAngle().getInUnits("deg");
   cout << " Lat = " << tg.lat();
   cout << " Lon =  " << tg.lon();
   cout << endl ;

   
   file << std::dec ;       // return write these numbers in decimal no matter what!
   file << tab << tg.altitude().getInUnits("km");
   file << tab << tg.range().getInUnits("km");
   file << tab << tg.doppler(lambda).getInUnits("hz");
   file << tab << tg.incidenceAngle().getInUnits("deg");
   file << tab << tg.lat().getInUnits("deg");
   file << tab << tg.lon().getInUnits("deg");

   if (base == 16) file << std::hex ;       // write numbers as hex to wtk file (as base defines)
   

   geometry_print = false;
  }
 */

 // end of target geometry insert


 file << endl; // end of liine for all

    }   // end vector<Ieb> for loop
  } // end of if for EDS on slow fast file donly (not power & tnc)
 file << std::dec ;       // return write numbers to decimal from now on

  }


//-----------------------
//located PRF hopping time
//------------------------
void IebList::locatePrfHopping(Array1D<Time>& hopping_time)
  {
    //first check record length
    if(ieb_list_.size()==0) ErrorMessage("No ieb record has been stored").throwMe();
    vector<Ieb> sar_ieb;//store SAR iebs
    sar_ieb.clear();

    //collect SAR normal operation mode iebs
    for(unsigned int i=0; i< ieb_list_.size();++i){
      Ieb ieb = ieb_list_[i];
      ieb.decodeSlowFastfield();
      //check whether it is SAR mode
      if(!(ieb.getMod()==2 || ieb.getMod()==3 ||
	   ieb.getMod()==10 || ieb.getMod()==11)) continue;//SARmode
      if(!(ieb.getCsr()==0 || ieb.getCsr()==8)) continue;//Image mode
      if(ieb.getBem()!=31) continue;//all 5Beams used
      sar_ieb.push_back(ieb);
    }
    cout<<"out of "<<ieb_list_.size()<<", "<< sar_ieb.size()<<" were collected for SAR"<<endl;
    hopping_time.resize(1);
    if(sar_ieb.size()==0) return;//no SAR iebs

    //find PRF hopping
    // fin number change 
    // PRF difference > 10 %
    vector<Time> time_mark;
    for(unsigned int i=0;i<sar_ieb.size()-1;++i){
      Ieb ieb1 = sar_ieb[i];
      ieb1.decodeSlowFastfield();
      Ieb ieb2 = sar_ieb[i+1];
      ieb2.decodeSlowFastfield();
      if(ieb1.getFin()==ieb2.getFin()) continue;
      if(ieb1.getPri()==ieb2.getPri()) continue;
      Uvar prf1 =1.0/ieb1.getPri();
      Uvar prf2 = 1.0/ieb2.getPri();
      float diff =fabs( (float) ((prf2- prf1)/prf1).getInUnits("")*100.0);
      if(diff > 9.0){
	if(time_mark.size()==0){
	  time_mark.push_back(ieb1.getTime());
	  time_mark.push_back(ieb2.getTime());
	  //cout<<"prf1 and prf2 "<<prf1<<" "<<prf2<<endl;
	}
	else{
	  if( (ieb1.getTime() -time_mark.back()) >Uvar(10,"s"))
	    time_mark.push_back(ieb1.getTime());//only if they are differentW
	  time_mark.push_back(ieb2.getTime());
	  //cout<<"Prf1 and prf2 "<< prf1<<" "<<prf2<<endl;
	}
      }
    }
    if(time_mark.size()==0) return;//no PRF hopping
    
    //get rid of duplicate time
    cout<<"PRF hopping "<< time_mark.size()<<endl;
    hopping_time.resize(time_mark.size());
    for(unsigned int i=0;i<time_mark.size();++i){
      hopping_time(i)=time_mark[i];
      //cout<<"hopping time "<< hopping_time(i).utc("ISOD")<<endl;
    }
    return;
  }
  
//-----------------------
//located PRF hopping division
//------------------------
void IebList::locatePrfHoppingDivision(Array1D<Time>& prf_start_time,
				       Array1D<Time>& prf_end_time,
				       unsigned int& number_of_divisions)
  {
    Array1D<Time> prf_hopping_time("");
    locatePrfHopping(prf_hopping_time);
    number_of_divisions=0;
    if(prf_hopping_time.size()==1) return;
    //----------------------------------
    //compute how many groups of prf hopping
    //---------------------------------
    unsigned int hopping_group=0;
    hopping_group=1;
    for(unsigned int i=0; i< prf_hopping_time.size()-1;++i){
      if( (prf_hopping_time(i+1)-prf_hopping_time(i)) > Uvar(5,"s"))
	hopping_group++;
    }
    
    //cout<<"PRF hopping group "<< hopping_group<<endl;
    //-------------------------------------------------
    //compute beginning and end of prf hopping group
    //-----------------------------------------------  
    number_of_divisions=hopping_group;
    prf_start_time.resize(hopping_group);
    prf_end_time.resize(hopping_group);
    //located each group's start time
    prf_start_time(0)= prf_hopping_time(0);
    unsigned int count=1;
    for(unsigned int i=1; i< prf_hopping_time.size();++i){
      if( (prf_hopping_time(i) - prf_hopping_time(i-1)) > Uvar(10,"s")){
	prf_start_time(count)=prf_hopping_time(i);
	//cout<<count<<" "<<prf_start_time(count).utc("ISOD")<<endl;
	count++;
      }
    }
    //located each group's end time
    count=0;
    for(unsigned int i=0; i< prf_hopping_time.size()-1;++i){
      if( (prf_hopping_time(i+1) - prf_hopping_time(i)) > Uvar(10,"s")){
	prf_end_time(count)=prf_hopping_time(i);
	//cout<<count<<" "<<prf_end_time(count).utc("ISOD")<<endl;
	count++;
      }
    }
    prf_end_time(hopping_group-1)= prf_hopping_time(prf_hopping_time.size()-1);
    
    
    for(unsigned int i=0; i< hopping_group;++i)
      cout<<prf_start_time(i).utc("ISOD")<<" "
	  <<prf_end_time(i).utc("ISOD")<<endl;
  }

//--------------------------------------
//CV: Compute PRF for SCATTEROMETER division
//--------------------------------------
Uvar IebList::computePRF(const Uvar& t)
  {

    // Initialize
    unsigned int Nphi = 1000;
    Uvec dop_delta("",Nphi); 
    Uvec range_delta("",Nphi);
    Uvec g_azi("",Nphi);
    Uvec g_elev("",Nphi);

    // Geometry setup and dop_delta, range_dalta calcualtion
    computeDelta_Dop_Rng(dop_delta, g_azi,  range_delta, g_elev, t); //Compute Delta Doppler & Range

    // ---------------
    // PRF CALCULATION
    // ---------------

    // The ranges below are chosen with Cassini radar geometry in 
    // mind where we know that PRF's are going to be between these bounds.
    Uvar fp_max = Uvar(10000,"Hz");
    Uvar fp_min = Uvar(300,"Hz");

   
    // ------------------
    // BISECTIONAL SEARCH
    // ------------------
    Uvar fp_try; 
    Uvar fp_acc = Uvar(1,"Hz");
    Uvar ydiff;
    Uvar amb_ratio_dop;
    Uvar amb_ratio_range;

    unsigned int Nbisection = log(fabs(fp_max-fp_min)/fp_acc)/log(2)+1;
   
    for(unsigned int i=0; i< Nbisection;++i)
      {
	fp_try = 0.5*(fp_max+fp_min);
   
	gamb_ratio_dop_rng(amb_ratio_dop, amb_ratio_range, dop_delta, g_azi,  range_delta, g_elev, t,fp_try);

	ydiff = amb_ratio_dop - amb_ratio_range;
	if (ydiff < 0)
	 fp_min = fp_try;
	else
	  fp_max = fp_try;


      }
    //cout<<"DOP "<<(t-epoch_time_)/60<<" "<<10*log(amb_ratio_dop)/log(10)<<endl;
    //cout<<"RANGE "<<(t-epoch_time_)/60<<" "<<10*log(amb_ratio_range)/log(10)<<endl;
    
    if(fp_try>Uvar(1300,"Hz"))
      {  
	fp_try = Uvar(1300,"Hz");   
	gamb_ratio_dop_rng(amb_ratio_dop, amb_ratio_range, dop_delta, g_azi,  range_delta, g_elev, t,fp_try);
      }

    cout<<"DOP "<<(t-epoch_time_)/60<<" "<<10*log(amb_ratio_dop)/log(10)<<endl;
    cout<<"RANGE "<<(t-epoch_time_)/60<<" "<<10*log(amb_ratio_range)/log(10)<<endl;

    // Return PRF
    // ----------
    return(fp_try);
} 




//--------------------------------
//CV:  METHOD computeDelta_Dop_Rng
//--------------------------------
void IebList::computeDelta_Dop_Rng( Uvec& dop_delta, Uvec& g_azi, Uvec& range_delta, Uvec& g_elev, const Uvar& t)
{
  float g_doplevel_dB = -3;
  float g_doplevel = pow(10,g_doplevel_dB/10);
  float g_rangelevel_dB = -3;
  float g_rangelevel = pow(10,g_rangelevel_dB/10);
  Uvar r_prf;
  Uvar phi_width;
  unsigned int Nphi = dop_delta.size();
  Uvec phi("",Nphi);
  Uvec spherical_gain("",Nphi);
  Uvec spherical_gain_frame2("",Nphi);
  Uvec spherical_gain_frame3("",Nphi);
  Uvec gain_dB("",Nphi);
  Uvec gain_dB_frame2("",Nphi);
  Uvec fdop_iso_range("",Nphi);
  Uvec range_iso_dop("",Nphi);
  Uvec costheta_iso_range("",Nphi);
  DirectionVector cassini_x, cassini_y, cassini_z;

  time_t ticTimeStart;
  time_t ticTimeEnd;
  float tictoc;

  //ticTimeStart = clock();
  
  // ----------
  // Beam Frame
  // ----------
  Frame frameb3("CASSINI_RADAR_"+toStr(3),"Cassini");
  
  //---------------
  // Target Geomety
  // --------------
  StateVector sc_state;
  target_frame_.ephemeris(sc_state,"Cassini",t,"NONE");
  
  
  // Set look direction w.r.t beam3 frame
  DirectionVector look("look",frameb3,t,0,0,1); 
  
  
  TargetGeom tg0(t);
  tg0.setState(sc_state);
  tg0.setTarget(target_name_);
  tg0.setLookDirection(look); 
 

  Uvar radius=tg0.radius();
  Uvar incidenceAngle = tg0.incidenceAngle();
  Uvar alt = tg0.altitude();
  Uvar range = tg0.range();
  
  FloatVector vcas = sc_state.velocity();
  DirectionVector uvel = vcas; 
  Uvar speed = vcas.magnitude();
  
 
  Flyby flyby(cfg_);
  Frame track_frame = flyby.trackFrame();
  
  DirectionVector usurf = tg0.normal();
  DirectionVector i_z("i_z",tg0.nadir());
  PositionVector rsurf = tg0.radius()*usurf;
  DirectionVector i_y ("y",cross(i_z,rsurf));
  DirectionVector i_x ("x",cross(i_y,i_z));
  
  PositionVector rsc = sc_state.position();
  Frame internal_frame("internal_frame",i_x,i_y,i_z);
  DirectionVector uisorange(internal_frame,t,0,0,1);
  StateVector sat_state("sat_state");
  PositionVector rlook =  rsc -rsurf; 
  DirectionVector ulook = rsc -rsurf; 
  Uvar costhetav = dot(uvel,ulook);
  Uvar thetav = acos(costhetav);
 
  
  // ---------
  // S/C Frame
  // ---------
  Frame fsc("CASSINI_SC_COORD","Cassini");
  
  //---------------------
  //read in radar frequency and wavelength
  //--------------------
  Uvar lambda = speed_light/carrier_frequency;
  
  
  // ------------------------------------
  //  Doppler Ambiguity calculation setup
  // ------------------------------------
  Uvar theta = Uvar(acos((pow(radius,2)+pow(alt+radius,2)-pow(range,2))/(2*radius*(radius+alt))),"rad");
  
  Uvar radiusIsoRange = sin(theta)*radius;
  
  Uvar beamwidth_elev_max = Bvector_[2].getElevationWidthOneWay();
  
  
  
  //------Set phi width-------
  if (radiusIsoRange == Uvar(0,"km"))
    {
      phi_width = 2*pi;
    }
  else
    {
      phi_width = Uvar((beamwidth_elev_max *range/radiusIsoRange),"rad");
      
    }
  
  
  Uvar limit_phi = Uvar(360/10*pi/180,"rad");
  
  phi_width = MIN(phi_width,limit_phi);
  
  
  Uvec b3_azi("",Nphi);
  Uvec b3_elev("",Nphi);
  Uvec index("",Nphi);
  Uvec dopl_delta_1("",Nphi);
  Uvec spherical_gain_norm("",Nphi);
  Uvec spherical_gain_norm1("",Nphi);
  Uvec dop_delta_1("",Nphi);
  Uvec range_delta_1("",Nphi);

  Uvar amb_ratio_range; 
  PositionVector rsurf_iso_dop;
  Uvar max_gain = 0;
  
  
  // -------Set phi---------
  
  phi_width = MIN(phi_width,limit_phi);
  //cout<<"phi_width"<<phi_width<<endl;
  
  for(unsigned int i=0; i< Nphi;++i)
    {
      phi(i) = 12*phi_width/Nphi*(i+1-round(Nphi/2));
      index(i) = i;
      
      DirectionVector uisorange(internal_frame,t,0,0,1);
      uisorange.setSpherical(theta,phi(i));
      uisorange.representIn(usurf);
      PositionVector rlook_iso_range = radius*uisorange-rsc;
      DirectionVector ulook_iso_range(rlook_iso_range);
      
      costheta_iso_range(i) = dot(uvel,ulook_iso_range);
      fdop_iso_range(i) = -2.0 * costheta_iso_range(i) * speed / lambda;
      
      //ulook_iso_range.representIn(frameb3);

      
	// Gain from Frame 2
      spherical_gain_frame2(i)= Bvector_[2].getGainOneWay(ulook_iso_range);
      //spherical_gain_frame2(i) =0;
  
      gain_dB_frame2(i) = 10*log(spherical_gain_frame2(i));

      if(i==0)
	{
	  Uvar max_gain = spherical_gain_frame2(i);
	}
      if(spherical_gain_frame2(i)> max_gain)
	{
	  max_gain = spherical_gain_frame2(i);
	}
      
    }
  
  
  spherical_gain_norm = spherical_gain_frame2/max_gain;
  
 
  unsigned int j = 0;
  
  do
    {
      j++;
    } while (spherical_gain_norm(j)<g_doplevel);
  
  // Note that doppler shifts are always reported as positive.
  dop_delta_1 = fabs(fdop_iso_range - fdop_iso_range(j));
    
  //  Zeros out delta's in the neg phi dir to avoid a false match.
  for(unsigned int i=0; i< Nphi;++i)
    {	
      if(i<j)
	{
	  dop_delta(i) = 0;
	  
	}
      else
	{
	  dop_delta(i) = dop_delta_1(i);
	  
	}
    }
  
      
    
  // ---------------------------------
  // Range Ambiguity calculation setup
  // ---------------------------------
  
  // Radius of iso-doppler circle (related to radius of iso-range circle)
  Uvar radius_d = sqrt(pow(range,2) - pow(radius,2));
  PositionVector radii = default_target_radii;
  
  // For iso-doppler contours, set coordinate system 3 with z-axis parallel 
  // to s/c velocity vector (which is also body relative).
  DirectionVector z3_p = uvel;
  
  
  // Set x axis in plane of incidence and orthogonal to z, putting the boresight
  // at phi=0.
  DirectionVector y3_p("y3_p",cross(z3_p,ulook));
  DirectionVector x3_p("x3_p",cross(y3_p,z3_p));
  
  
  Frame internal_frame2("internal_frame2",x3_p,y3_p,z3_p);
  
   ticTimeStart = clock();
  
  
  // Trace iso-doppler circles
  
  if (radius_d == 0)
    phi_width = Uvar(2*pi,"rad");
  else
    phi_width = beamwidth_elev_max*range/radius_d;
  
  phi_width = MIN(phi_width,limit_phi);
  

  for(unsigned int i=0; i< Nphi;++i)
    {
      phi(i) = 50*phi_width/Nphi*(i+1-round(Nphi/2));
      DirectionVector uisodop(internal_frame2,t,0,0,1);
      uisodop.setSpherical(thetav,phi(i));
      uisodop.representIn(usurf);
      

      tg0.reset(t);
      tg0.setState(sc_state);
      tg0.setTarget(target_name_);
      tg0.setLookDirection(-uisodop);
      PositionVector rsurf_iso_dop=tg0.surfaceIntercept();
      
      //spice_surfpt(rsurf_iso_dop,spice_found,rsc, -uisodop,target_radii);

      PositionVector rlook_iso_dop = rsurf_iso_dop - rsc;
      DirectionVector ulook_iso_dop = rlook_iso_dop;
      range_iso_dop(i) = rlook_iso_dop.magnitude();
      //ulook_iso_dop.representIn(frameb3);
      
      spherical_gain_frame3(i)= Bvector_[2].getGainOneWay(ulook_iso_dop);
      //spherical_gain_frame3(i)= 0;
      if(i==0)
	{
	  Uvar max_gain = spherical_gain_frame3(i);
	}
      if(spherical_gain_frame3(i)> max_gain)
	{
	  max_gain = spherical_gain_frame3(i);
	}
    }
  Uvec spherical_gain_norm2("",Nphi);
  spherical_gain_norm2 = spherical_gain_frame3/max_gain;
  
  unsigned int k = 0;
  do
    {
      k++;
    } while (spherical_gain_norm2(k)<g_rangelevel);
  
  // Note that doppler shifts are always reported as positive.
  range_delta_1 = fabs(range_iso_dop - range_iso_dop(k));
  
  
  //  Zeros out delta's in the neg phi dir to avoid a false match.
  for(unsigned int i=0; i< Nphi;++i)
    {	
	if(i<k)
	  {
	    range_delta(i) = 0;
	    
	  }
	else
	  {
	    range_delta(i) = range_delta_1(i);
	    
	  }
    }

   
    g_azi = spherical_gain_norm;
    g_elev = spherical_gain_norm2;

}

//-----------------------
// CV: METHOD gamb_ratios 
//-----------------------
void IebList::gamb_ratio_dop_rng(Uvar& amb_ratio_dop, Uvar& amb_ratio_range,  const Uvec& dop_delta, const Uvec& g_azi, const Uvec& range_delta, const Uvec& g_elev, const Uvar& t, const Uvar& fp_try)
{
  int Nphi = 1000;
  float g_doplevel_dB = -3;
  float g_doplevel = pow(10,g_doplevel_dB/10);
  float g_rangelevel_dB = -3;
  float g_rangelevel = pow(10,g_rangelevel_dB/10);


 // gamb ratios azim (dop)
 int l = 0;
    
  do
    {
      l++;
      if (l==Nphi-1)
	{
	  l = Nphi-1;
	  break;
	}
    } while (dop_delta(l)<fp_try);
  
  int i_dop_ambig = l;
  if(g_azi(i_dop_ambig) == 0)
    amb_ratio_dop = 0.01;
  else if(i_dop_ambig == Nphi-1)
    amb_ratio_dop = 1000;
  else
    amb_ratio_dop = pow(g_doplevel/g_azi(i_dop_ambig),2);
  


 Uvar range_ambig_spacing = speed_light/2/fp_try;
  // gamb ratios elev (range)
  int m = 0; 
  
  do
    {
      m++;
      if (m==Nphi-1)
	{
	  m = Nphi-1;
	  break;
	  }
    } while (range_delta(m)<range_ambig_spacing);
  
  int i_range_ambig = m;

  if(g_elev(i_range_ambig) == 0)
    amb_ratio_range = 0.01;
  else if(i_range_ambig == Nphi-1)
    amb_ratio_range = 1000;
  else
    amb_ratio_range = pow(g_rangelevel / g_elev(i_range_ambig),2);

}

