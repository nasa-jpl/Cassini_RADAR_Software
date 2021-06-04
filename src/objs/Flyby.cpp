//----------------------------------------------------------------------------
// Flyby.cpp
//
// This file contains method definitions for the Flyby class.
// Flyby provides support for target flyby calculations.
// In particular, it contains named epoch times, and the ability to locate
// closest approach times.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char flybyspec_c[] =
  "@(#) $Id: Flyby.cpp,v 11.6 2012/09/27 20:01:58 richw Exp $";

#include <strings.h>
#include <string>
#include <map>
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Frame.h"
#include "Time.h"
#include "TargetGeom.h"
#include "TemplateUtils.h"
#include "Ckernel.h"
#include "Flyby.h"
#include "Plot.h"
#include "Utils.h"

using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;

//------------------------
// Class Flyby Methods
//------------------------

//------------------
// Automatic testing
//------------------

bool Flyby::selfTest()
  {
  return(true);
  }

//-------------
// Constructors
//-------------

//--------------------------------------------------------
// The only construction method is from a Config object.
//--------------------------------------------------------

Flyby::Flyby(Config& cfg)
  : cfg_filename_(cfg.filename()), percent_completed_(0),
    lowest_altitude_set_(false),track_frame_set_(false),
    computePolyfitCoeff_set_(false),local_offset_set_(false),
    altitude_polyfit_("altitude values for polyfit"),
    angle_polyfit_("angle values for polyfit"),
    angle_poly_coeff_time_("poly coeff time "),
    target_lat_lon_fit_("poly fit coeff for special target lat lon fit"),
    angle_offset_("offset angle"), tr1_("tr1"),tr2_("tr2"),tr3_("tr3"),tr4_("tr4"),tr12_mid_("mid time between tr1 and tr2"),tr34_mid_("mid time between tr3 and tr4"),tr_from_ideal_("tr from ideal"),tr_to_ideal_("transition toward ideal pointing")
 {

   config(cfg);
   Noffset_=0;
 }

//-------------
// Setup
//-------------

//--------------------------------------------------------------------------
// config(cfg)
//
// Read the config parameters in cfg which specify a time interval to deal
// with.  The results are returned as an epoch_time, and an epoch
// relative set of time steps from start_time to end_time stepping by
// interval.
//
// The following config keywords are processed by this function:
//
// epoch_selection (absolute,closest_flyby,designated_flyby,ckernel_valid)
// epoch_time (UTC time string)
// time_step
// start_time
// end_time
// epoch_accuracy
//--------------------------------------------------------------------------

void Flyby::config(Config& cfg)
  {
  epoch_selection_ = cfg.str("epoch_selection");
  epoch_accuracy_ = cfg["epoch_accuracy"];
  target_ = cfg.str("target");
  spacecraft_ = cfg.str("spacecraft");
  lowest_altitude_set_ = false;

  //--------------------------------------------------------------------
  // Determine ckernel start and end times for potential substitution
  // for the epoch time, or for the time interval boundaries.
  //--------------------------------------------------------------------

  list<string> ckernel_names;
  Frame::getCkernelList(ckernel_names);
  for (list<string>::const_iterator p=ckernel_names.begin();
       p != ckernel_names.end(); ++p)
    {
    Time ck_start_time,ck_end_time;
    Ckernel::validTimeRange(*p,spacecraft_,ck_start_time,ck_end_time);
    if (p == ckernel_names.begin())
      {
      ck_earliest_start_time_ = ck_start_time;
      ck_latest_end_time_ = ck_end_time;
      }
    if (ck_start_time < ck_earliest_start_time_)
      {
      ck_earliest_start_time_ = ck_start_time;
      }
    if (ck_end_time > ck_latest_end_time_)
      {
      ck_latest_end_time_ = ck_end_time;
      }
    }

  epoch_time_ = Uvar(0,"s");
  epoch_time_ = timeValueSubstitute(cfg,"epoch_time");

  if (epoch_selection_ == "closest_flyby")
    {
    //-------------------------------------------------------------------
    // User supplies an epoch_time close to the desired flyby, and epoch
    // relative start and end times and the interval.  The epoch of the
    // closest approach of the flyby closest to the input epoch_time is
    // located and used as the epoch_time.  In this case, epoch_time is
    // an input/output variable and is modified.
    //-------------------------------------------------------------------

    //-------------------------------------------------------------
    // Set epoch to the lowest altitude time.
    //-------------------------------------------------------------

    epoch_time_ = lowestAltitudeTime();
    }
  else if (epoch_selection_ != "absolute")
    {
    ErrorMessage e("Invalid epoch_selection parameter: " + epoch_selection_);
    e.throwMe();
    }

  //----------------------------------------------------------------
  // Get all of the time iterval specifications
  // Substitude epoch relative values for any absolute time strings
  //----------------------------------------------------------------

  interval_ = cfg["time_step"];
  start_time_ = timeValueSubstitute(cfg,"start_time");
  end_time_ = timeValueSubstitute(cfg,"end_time");

  cfg.suffixSet("start_time_",start_time_map_);
  cfg.suffixSet("end_time_",end_time_map_);
  cfg.suffixSet("time_step_",interval_map_);

  Config::strings str_start_time_map;
  cfg.suffixSet("start_time_",str_start_time_map);
  for (Config::strings::const_iterator pstart = str_start_time_map.begin();
       pstart != str_start_time_map.end(); ++pstart)
    {
    string suffix = pstart->first;
    Uvar value = timeValueSubstitute(cfg,"start_time_" + suffix);
    start_time_map_[suffix] = value;
    }
  
  Config::strings str_end_time_map;
  cfg.suffixSet("end_time_",str_end_time_map);
  for (Config::strings::const_iterator pend = str_end_time_map.begin();
       pend != str_end_time_map.end(); ++pend)
    {
    string suffix = pend->first;
    Uvar value = timeValueSubstitute(cfg,"end_time_" + suffix);
    end_time_map_[suffix] = value;
    }
 
  //--------------------------------------- 
  // Finish setting up the time intervals
  //--------------------------------------- 

  for (Config::numbers::const_iterator pstart = start_time_map_.begin();
       pstart != start_time_map_.end(); ++pstart)
    {
    string suffix = pstart->first;
    Config::numbers::const_iterator pend =
      end_time_map_.find(suffix);
    if (pend == end_time_map_.end())
      {
      ErrorMessage e("Missing end time for suffix: " + suffix);
      e.throwMe();
      }
    Config::numbers::const_iterator pinterval =
      interval_map_.find(suffix);
    if (pinterval == interval_map_.end())
      {
      ErrorMessage e("Missing time step for suffix: " + suffix);
      e.throwMe();
      }
    Uvar dur = pend->second - pstart->second;
    ntime_map_[suffix] = (unsigned int)(dur/pinterval->second).getInUnits("");
    resetTimeSteps(suffix);
    }
  
  Uvar duration = end_time_ - start_time_;
  ntime_ = (unsigned int)(duration / interval_).getInUnits("");

  resetTimeSteps();

  //------------------------------------------------------------------
  // The following config keywords are read for potential use by
  // the method computePolyfitCoeff().
  // SC orientation, azimuth rotation,sc_z_rotation, tracking option
  //-------------------------------------------------------------------

  sc_orientation_ =  cfg.str("sc_orientation");
  tracking_option_= cfg.str("tracking_option");
  azimuth_rotation_= cfg["set_azimuth_rotation"];
  sc_z_rotation_ = cfg["set_sc_z_rotation"];

  //------------------------------------------
  // Compute incidence angle profile coefficients
  // depending on incidence angle profile.
  // Needed by dlap_generate
  //------------------------------------------
  ivd_time_pad_ = cfg["ivd_time_pad"];
  angle_profile_option_=cfg.str("angle_profile_option");
  if (angle_profile_option_ == "constant_incidence_angle")
    {
    angle_polyfit_.resize(1);//single value
    angle_polyfit_(0)=cfg["constant_incidence_angle"];
    }
  else if (angle_profile_option_=="polynomial_fit")
    {
    unsigned int N = (unsigned int)
      cfg["number_of_altitude_incidenceAngle_variations"].getInUnits("");
    altitude_polyfit_.resize(N);
    angle_polyfit_.resize(N);

    unsigned int poly_power =
      (unsigned int) cfg["polynomial_power"].getInUnits("");
    poly_power += 1;
    angle_poly_coeff_time_.resize(poly_power);
    if (N < 2)
      {
      ErrorMessage
        ("For poly fit, number_of_incidenceAngle_variations > 1").throwMe();
      }
    for (unsigned int i = 0; i < N;++i)
      {
	if(cfg.keywordExists("altitude01")){//new format
	  altitude_polyfit_(i)= cfg["altitude"+toStr(i+1,2)];    
	  angle_polyfit_(i) = cfg["incidenceAngle"+toStr(i+1,2)];
	}
	else{//old format
	  altitude_polyfit_(i)= cfg["altitude"+toStr(i+1)];    
	  angle_polyfit_(i) = cfg["incidenceAngle"+toStr(i+1)];
	}
      }
    for (unsigned int i = 0; i < N-1;++i)
      {
      if(altitude_polyfit_(i) > altitude_polyfit_(i+1))
	{
	  ErrorMessage("Flybyb.cpp:: reorder altitude in increasing order").throwMe();
	}    
      }
    }
  else if (angle_profile_option_=="sech_function")
    {
    angle_at_ca_ = cfg["incidence_angle_at_lowest_altitude"];
    mid_alt_ =cfg["mid_altitude"];
    angle_at_mid_alt_=cfg["incidence_angle_at_mid_altitude"];
    residual_angle_=cfg["residual_incidence_angle_at_high_altitude"];
    }
 else
   {
   ErrorMessage("Flyby.cpp::cfg(): no valid angle profile option").throwMe();
   }

  //-------------------------------------------------
  //Special case when target_lat_lon option is "On"
  //-------------------------------------------------
  target_lat_lon_option_=cfg.str("target_lat_lon_option");
 
  if(strcasecmp(target_lat_lon_option_.c_str(),"on")==0){
    if(cfg.keywordExists("target_isodoppler_pointing")){
      if(strcasecmp(cfg.str("target_isodoppler_pointing").c_str(),"on")==0)
	target_lat_lon_option_="On";
      else {
	target_lat_lon_option_="Off";
      }
    }
  }
  else{
    target_lat_lon_option_="Off";
  }
  
  
  if (target_lat_lon_option_=="On")
    {
      cout<<"-------------------------"<<endl;
      cout<<"Isodoppler target lat lon option is selected"<<endl;
      cout<<"--------------------------"<<endl;
      off_centered_beam_time_lag_ = cfg["off_centered_beam_time_lag"];
      start_target_lat_lon_tracking_=cfg["start_target_lat_lon_tracking"];
      end_target_lat_lon_tracking_=cfg["end_target_lat_lon_tracking"];
      angle_at_closest_range_to_target_= cfg["angle_at_closest_range_to_target"];
      time_at_closest_range_to_target_=cfg["time_at_closest_range_to_target"];  
    }
  }

//-------------------------------------------------------------------------
// setTarget(target)
// setObserver(observer)
//
// Set a new target or observer which may be different from what was
// specified in the config file used to construct this Flyby object.
//-------------------------------------------------------------------------

void Flyby::setTarget(const string& target)
  {
  target_ = target;
  lowest_altitude_set_ = false;
  }

void Flyby::setObserver(const string& observer)
  {
  spacecraft_ = observer;
  lowest_altitude_set_ = false;
  }

//-------------------------------------------------------------------------
// setEpoch(epoch_selection,epoch_time,epoch_accuracy)
//
// Set a new epoch time which may be different from the original time
// set in the config file.  Epoch selection is a string set to either
// "absolute" (epoch_time is the epoch time to use), or "closest_flyby"
// (epoch_time is starting point in search for nearest closest approach
// to the target).  epoch_accuracy is the final accuracy of the search
// (if used).
//-------------------------------------------------------------------------

void Flyby::setEpoch(const string& epoch_selection, const Time& t,
  const Uvar& epoch_accuracy)
  {
  lowest_altitude_set_ = false;
  epoch_selection_ = epoch_selection;
  epoch_accuracy_ = epoch_accuracy;
  epoch_time_ = t;
  if (epoch_selection_ == "closest_flyby")
    {
    epoch_time_ = lowestAltitudeTime();
    }
  else if (epoch_selection_ != "absolute")
    {
    ErrorMessage e("Flyby::setEpoch: Invalid epoch selection: " +
        epoch_selection);
    e.throwMe();
    }
  }

//----------------
// Predicates
//----------------

bool Flyby::validSuffix(const string& suffix) const 
  {
  map<string,unsigned int>::const_iterator p = ntime_map_.find(suffix);
  return(p != ntime_map_.end());
  }


//-----------------------------------------------------------------------
// done()
// done(suffix)
//
// Returns true if the time steps have reached the end_time,
// false otherwise.
//-----------------------------------------------------------------------

bool Flyby::done() const
  {
  if (time_counter_ >= ntime_) return(true);
  return(false);
  }
  
bool Flyby::done(const string& suffix)
  {
  suffixCheck(suffix);
  if (time_counter_map_[suffix] >= ntime_map_[suffix]) return(true);
  return(false);
  }
  
//--------------------------
// Special Service Methods
//--------------------------

void Flyby::suffixCheck(const string& suffix) const 
  {
  if (!validSuffix(suffix))
    {
    ErrorMessage e("Time interval suffix (" + suffix +
      ") not present in " + cfg_filename_);
    e.throwMe();
    }
  }


//------------------------------------------------------------------
// start_time = startTime()
// end_time = endTime()
// interval = interval()
//
// Epoch relative start and end times.
// Time step interval.
//------------------------------------------------------------------

Uvar Flyby::startTime() const 
  {
  return(start_time_);
  }

Uvar Flyby::endTime() const 
  {
  return(end_time_);
  }

Uvar Flyby::interval() const 
  {
  return(interval_);
  }

Uvar Flyby::startTime(const string& suffix) const 
  {
  suffixCheck(suffix);
  return(start_time_map_[suffix]);
  }

Uvar Flyby::endTime(const string& suffix) const 
  {
  suffixCheck(suffix);
  return(end_time_map_[suffix]);
  }

Uvar Flyby::interval(const string& suffix) const 
  {
  suffixCheck(suffix);
  return(interval_map_[suffix]);
  }

//------------------------------------------------------------------
// epoch_time = epochTime()
// epoch_accuracy = epochAccuracy()
//
// Flyby epoch time (usually close to closest approach).
// Epoch accuracy (used by search for closest approach).
//------------------------------------------------------------------

Time Flyby::epochTime() const 
  {
  return(epoch_time_);
  }

Uvar Flyby::epochAccuracy() const 
  {
  return(epoch_accuracy_);
  }

//------------------------------------------------------------------
// N = numTimeSteps()
//
// Number of time steps specified for this flyby.
//------------------------------------------------------------------

unsigned int Flyby::numTimeSteps() const 
  {
  return(ntime_);
  }

unsigned int Flyby::numTimeSteps(const string& suffix) const 
  {
  suffixCheck(suffix);
  return(ntime_map_[suffix]);
  }

//------------------------------------------------------------------
// target_str = target()
//
// Name of target specified for this flyby.
//------------------------------------------------------------------

string Flyby::target() const 
  {
  return(target_);
  }

//------------------------------------------------------------------
// id_str = id()
//
// ID string for this flyby (eg., T7, E4 etc).
//------------------------------------------------------------------

string Flyby::id() const 
  {
  string id;
  return(id);
  }

//------------------------------------------------------------------
// a = lowestAltitude();
//
// Compute lowest altitude.
// Note that this method pulls spacecraft position from the SPICE
// ephemeris.
//------------------------------------------------------------------

Uvar Flyby::lowestAltitude()
  {
  if (lowest_altitude_set_)
    {
    return(lowest_altitude_);
    }

  //-------------------------------------------------------------
  // Locate nearest closest approach using golden section search.
  //-------------------------------------------------------------

  Uvar t_mid("t_mid");
  t_mid = epoch_time_ + Time(Uvar(0.001,"s"));
  Uvar t2("t2");
  Uvar final_time("final_time");
  TargetAltitude tga(spacecraft_,target_);
  Time t1 = epoch_time_;  // need a copy to avoid changing epoch_time_
  bracket_minimum(t1,t_mid,t2,tga);
  golden_section_search(t1,t_mid,t2,epoch_accuracy_,tga,
    final_time,lowest_altitude_);
  lowest_altitude_time_.setEt(final_time);
  lowest_altitude_set_ = true;

  return(lowest_altitude_);
  }

//------------------------------------------------------------------
// t = lowestAltitudeTime();
//
// Compute the time of lowest altitude.
// Note that this method pulls spacecraft position from the SPICE
// ephemeris.
//------------------------------------------------------------------

Time Flyby::lowestAltitudeTime()
  {
  if (!lowest_altitude_set_)
    {
    lowestAltitude();
    }

  return(lowest_altitude_time_);
  }

//------------------------------------------------------------------
// s = lowestAltitudeState();
//
// Return the spacecraft statevector at the point of closest approach.
// Note that this method pulls spacecraft position from the SPICE
// ephemeris.
//------------------------------------------------------------------

StateVector Flyby::lowestAltitudeState()
  {
  if (!lowest_altitude_set_)
    {
    lowestAltitude();
    }
  StateVector s;
  Frame fj2000("J2000",target_);
  fj2000.ephemeris(s,spacecraft_,lowest_altitude_time_,"NONE");
  return(s);
  }

//------------------------------------------------------------------
// ftrack = trackFrame();
//
// Return the track frame for the current flyby trajectory.
// The track frame puts the x-axis through the spacecraft nadir point
// at closest approach.  The y and z axes are then arranged so that the
// xy plane is tangential to the ground track at closest approach.
// Thus, the ground track will stay close to the xy plane.
//------------------------------------------------------------------

Frame Flyby::trackFrame()
  {
  if(!track_frame_set_) computeTrackFrame();
  return(track_frame_);
  }

void Flyby::computeTrackFrame()
  {
  if (!lowest_altitude_set_)
    {
    lowestAltitude();
    }
  if(track_frame_set_) 
    {
    ErrorMessage
      ("Flyby.cpp::computeTrackFrame(): trackFrmae is alreay set").throwMe();
    }
  TargetGeom tg(lowest_altitude_time_);
  StateVector s = lowestAltitudeState();
  tg.setState(s);
  tg.setTarget(target_);
  Frame target_frame("IAU_" + target_, target_);
  DirectionVector track_x("track_x",tg.nadir().representIn(target_frame));
  DirectionVector track_z("track_z",
    cross(s.velocity().representIn(target_frame),track_x));
  DirectionVector track_y("track_y",cross(track_z,track_x));
  track_frame_=Frame("track_frame", track_x, track_y, track_z);
  track_frame_set_ = true;
  }

//----------------------------------------------------------------
// LookDirection(): provide look direction(right or left looking)
//------------------------------------------------------------------

string Flyby::LookDirection()
  {
  StateVector s = lowestAltitudeState();
  DirectionVector s_pos = s.position();
  DirectionVector s_vel = s.velocity();
  DirectionVector s_cross= cross(s_pos,s_vel);

  Frame fbeam3("CASSINI_RADAR_3","Cassini");
  DirectionVector beam3("beam3 boresight direction",fbeam3,
    lowest_altitude_time_,0,0,1);
  
  double roll_direction = dot(beam3,s_cross);
  if (roll_direction < 0)
    {
    return("Right");
    }
  else
    {
    return("Left");
    }
  }

//-----------------------------------------
//timeOffset(const Time& time)
//  Calculate time offset required for beam2 to illuminate 
//  the same area beam3 illuminates earlier or later
//  depending on sc attitude
//  If beam2 or 4 does not have any surface intercept point,
//  Uvar(0,"s") will be returned. Such case can occur when the Radar
//  is in Radiometer or  ALTL mode.
// 
//-----------------------------------------

Uvar Flyby::timeOffset(const Time& t)
  {
  StateVector sc_state("sc_state");
  Frame ftitan("IAU_TITAN",target_);
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");
    
  //beam frame3 and boresight direction
  Frame fbeam3("CASSINI_RADAR_3","Cassini");
  DirectionVector look_bore_b3("bore",fbeam3,t,0,0,1);
  Uvar bore_along_b3,bore_cross_b3;
  TargetGeom tg_bore_b3(t);
  tg_bore_b3.setState(sc_state);
  tg_bore_b3.setTarget(target_,ftitan);
  tg_bore_b3.setTrackFrame(trackFrame());
  tg_bore_b3.setLookDirection(look_bore_b3);
  tg_bore_b3.interceptAlongCross(bore_along_b3,bore_cross_b3);
    
  //beam frame2 
  Frame  fbeam2("CASSINI_RADAR_2","Cassini");
  Frame  fbeam4("CASSINI_RADAR_4","Cassini");
  DirectionVector look_bore_b2("beam 2 bore",fbeam2,t,0,0,1);
  Uvar bore_along_b2,bore_cross_b2;
  TargetGeom tg_bore_b2(t);
  tg_bore_b2.setState(sc_state);
  tg_bore_b2.setTarget(target_,ftitan);
  tg_bore_b2.setTrackFrame(trackFrame());
  tg_bore_b2.setLookDirection(look_bore_b2);
  tg_bore_b2.interceptAlongCross(bore_along_b2,bore_cross_b2);
    
  //if beam 2 is behind beam 3, sign = 1
  //otherwise sign = -1
  int sign = 1;
  if (bore_along_b2 < bore_along_b3)
    {
    sign = 1;
    }
  else
    {
    sign = -1;
    }
   
  cout<<"b3 along cross "<<bore_along_b3<<" "<<bore_cross_b3<<endl;     
  cout<<"b2 along cross "<<bore_along_b2<<" "<<bore_cross_b2<<endl;
   
  unsigned int Ntime_lag_step = 120*4;
  Uvec b2_along("b2_along",Ntime_lag_step);
  Uvec b4_along("b4_along",Ntime_lag_step);
  Array1D<unsigned int> valid_b2("",Ntime_lag_step);
  Array1D<unsigned int> valid_b4("",Ntime_lag_step);

  valid_b2=0;
  valid_b4=0;

  Uvar crosstrack;
  Uvec time_lag("time_lag",Ntime_lag_step);
  for (unsigned int i_lag = 0; i_lag < Ntime_lag_step;++i_lag)
    {
    time_lag(i_lag) = Uvar(double(i_lag),"s")*double(sign)/4.0;
    Time t_lag = t + time_lag(i_lag);
    DirectionVector look_bore_b2("bore",fbeam2,t_lag,0,0,1);
    ftitan.ephemeris(sc_state,"Cassini",t_lag,"NONE");

    //beam2 boresight
    TargetGeom tg_b2_bore(t_lag);
    tg_b2_bore.setState(sc_state);
    tg_b2_bore.setTarget(target_,ftitan);
    tg_b2_bore.setTrackFrame(trackFrame());
    tg_b2_bore.setLookDirection(look_bore_b2);
    if(tg_b2_bore.foundSurfaceIntercept()){
      tg_b2_bore.interceptAlongCross(b2_along(i_lag),crosstrack);
    }
    else{
      // ErrorMessage("no surface intercept for beam 2").throwMe();
      valid_b2(i_lag)=1;
    }

    //beam4 boresight
    DirectionVector look_bore_b4("bore",fbeam4,t_lag,0,0,1);
    TargetGeom tg_b4_bore(t_lag);
    tg_b4_bore.setState(sc_state);
    tg_b4_bore.setTarget(target_,ftitan);
    tg_b4_bore.setTrackFrame(trackFrame());
    tg_b4_bore.setLookDirection(look_bore_b4);
    if(tg_b4_bore.foundSurfaceIntercept()){
      tg_b4_bore.interceptAlongCross(b4_along(i_lag),crosstrack);
    }
    else{
      // ErrorMessage("no surface intercept for beam 4").throwMe();
      valid_b4(i_lag)=1;
    }
    }
    
    //cout<<"sign "<< sign<<endl;
    unsigned int b2_lag_index,b4_lag_index;
  b2_lag_index = 0;
  b4_lag_index = 0;
  bool found_b2=false;
  bool found_b4=false;
  if (sign == 1)
    {
    for (unsigned int i_lag = 0; i_lag < Ntime_lag_step-1;++i_lag)
      {	
	if(valid_b2(i_lag)==1 || valid_b4(i_lag)==1) continue;
	if (b2_along(i_lag) <  bore_along_b3 
	    && b2_along(i_lag+1) > bore_along_b3)
	  {
	    b2_lag_index = i_lag;
	    found_b2=true;
	  }
	
	if (b4_along(i_lag) <  bore_along_b3 
	    && b4_along(i_lag+1) > bore_along_b3)
	  {
	    b4_lag_index = i_lag;
	    found_b4=true;
	  }
	if(found_b2 && found_b4) break;
      }
    }
  else 
    {
    for (unsigned int i_lag = 0; i_lag < Ntime_lag_step-1;++i_lag)
      {	
	if(valid_b2(i_lag)==1 || valid_b4(i_lag)==1) continue;
	if (b2_along(i_lag) >  bore_along_b3 
	    && b2_along(i_lag+1) < bore_along_b3)
	  {
	    b2_lag_index = i_lag;
	    found_b2=true;
	  }
	if (b4_along(i_lag) >  bore_along_b3 
	    && b4_along(i_lag+1) < bore_along_b3)
	  {
	    b4_lag_index = i_lag;
	    found_b4=true;
	  }
	if(found_b2 && found_b4) break;
      }
    }
  if (b2_lag_index==0 || b4_lag_index==0)
    {
      cout<<"0 s will be returned because it can not calculate time lag within "<< Ntime_lag_step/4/60<<"  min"<<endl;
      ErrorMessage("increase time step ").throwMe();
    }  
  
  cout<<"beam 3 bore "<< bore_along_b3<<endl;
  cout<<"index value "<< b2_lag_index<<" "<<b4_lag_index<<endl;
  cout<<"time lag 2 and 4 "<< time_lag(b2_lag_index)<<" "<<time_lag(b4_lag_index)<<endl;
  cout<<"beam 2 and 4 alongtrack "<<b2_along(b2_lag_index)<<" "<< b4_along(b4_lag_index)<<endl;

  return( (time_lag(b2_lag_index)+time_lag(b4_lag_index))/2.0 );   
  }

//-----------------------
// Spacecraft orientation
//-----------------------

string Flyby::sc_Orientation()
  {
  return(sc_orientation_);
  }

//---------------------------------
// Spacecraft tracking option
//---------------------------------

string Flyby::trackingOption()
  {
  return(tracking_option_);
  }

//----------------------------------
// Spacecraft azimuth rotation
//---------------------------------

Uvar Flyby::azimuthRotation()
  {
  return(azimuth_rotation_);
  }

//-------------------------------------------
// Space rotation w.r.t. z(beam3 boresight)
//--------------------------------------------

Uvar Flyby::sc_z_rotation()
  {
  return(sc_z_rotation_);
  }


//------------------------------------------------------
// Compute polynomial coeff for incidence angle vs time
//------------------------------------------------------

void Flyby::computePolyfitCoeff()
  {
  if(computePolyfitCoeff_set_) 
    {
    return;
    }
  if(angle_profile_option_ !="polynomial_fit")
    {
    cout<<"try to access polynomial fit when angle profile option does not ";
    cout<<" require to compute polynomial coefficients"<<endl;
    return;
    }
  cout<<"Converting input data of angle vs altitude into angle vs time"<<endl;
  cout<<"lowest altitude  "<<lowestAltitude()<<endl;
  Uvar min_altitude = lowestAltitude()-Uvar(10,"km");

  // Want to reduce incidence angle by 1/2 at 1000 km higher alt
  // than the last one
  Uvar max_alt1 = altitude_polyfit_(altitude_polyfit_.size()-1)+Uvar(1000,"km");
  Uvar angle_at_max_alt1 = angle_polyfit_(angle_polyfit_.size()-1)/2.0;

  vector<Uvar> alt_larger_than_min_altitude;
  vector<Uvar> angle_larger_than_min_altitude;
  alt_larger_than_min_altitude.clear();
  angle_larger_than_min_altitude.clear();
  for (unsigned int i = 0; i < altitude_polyfit_.size();++i)
    {
      if(altitude_polyfit_(i) > min_altitude)
	{
	  alt_larger_than_min_altitude.push_back(altitude_polyfit_(i));
	  angle_larger_than_min_altitude.push_back(angle_polyfit_(i));
	  //cout<<"alt and angle "<<altitude_polyfit_(i)<<" "<<angle_polyfit_(i)<<endl;
	}
    }
 
  cout<<"Number of points to be used for polyfit is "<<alt_larger_than_min_altitude.size()<<endl;
  cout<<"Number of input points from config file is "<< altitude_polyfit_.size()<<endl; 

  if(alt_larger_than_min_altitude.size() != altitude_polyfit_.size())
    {
      cout<<"Warning:: some input altitudes are higher than the lowest altitude "<<endl;
      cout<<"Warning:: input altitudes larger than lowest altitude will not be considered valid while computing alt vs time fit coefficient "<<endl;
      Uvar low_angle, low_alt;
      low_angle = angle_larger_than_min_altitude.front();
      low_alt = alt_larger_than_min_altitude.front();
      //assume we keep the low angle at the lowest altitude
      //without this "addition" of a single data, polyfit usually
      //produce a larger deviation at or near lowest altitude
      alt_larger_than_min_altitude.insert(alt_larger_than_min_altitude.begin(),lowestAltitude());
      angle_larger_than_min_altitude.insert(angle_larger_than_min_altitude.begin(),low_angle);
    }
 
  

  unsigned int Npoints = alt_larger_than_min_altitude.size();
  Uvec alt_tmp("alt_tmp",Npoints+1);
  Uvec angle_tmp("angle_tmp",Npoints+1);
  for(unsigned int i = 0; i<Npoints;++i)
    {
      alt_tmp(i)=alt_larger_than_min_altitude[i];
      angle_tmp(i)=angle_larger_than_min_altitude[i];
    }
  alt_tmp(Npoints) = max_alt1;
  angle_tmp(Npoints) = angle_at_max_alt1;
  
  
 
  //-----------------------------------
  //plot alt vs incidence angle
  //-----------------------------------
  Plot plot0;
  plot0.addXY(alt_tmp,"km",angle_tmp,"deg",line("none"),sym("circle","red",2));
  plot0.setXlabelandUnits("altitude");
  plot0.setYlabelandUnits("deg");
  //plot0.show("x");
     
  unsigned int N = alt_tmp.size();   
  Uvec time("time",2*N);
  Uvec angle("angle",2*N);
  Time epoch = lowestAltitudeTime();
  Uvar accuracy("accuray",0.1,"s");

  //find times for given set of altitude values
  //use golden section search 
  //find time for a given altitude before epoch
  for (unsigned int i = 0; i < N;++i)
    {
    angle(N-1-i) = angle_tmp(i);
    Time t1 = epoch - Uvar(30*60,"s");//start search 30 min before epoch
    Time t2 = epoch + Uvar(0.1,"s");//1 s after lowest altitude time
    Time t_mid = (t1 + t2)/2.0;
    TargetTimeAltitude tta(spacecraft_,target_,alt_tmp(i));
    Uvar final_time,final_altitude;
    golden_section_search(t1,t_mid,t2,accuracy,tta,final_time,final_altitude);
    time(N-1-i) = final_time-epoch;
    //cout<<"check result final time final alt "<<final_time-epoch<<" "<<final_altitude<<" "<<alt_tmp(i)<<endl;
    }

  //find time for a given altitude after epoch
  for (unsigned int i = 0; i < N;++i)
    {
    angle(i+N) = angle_tmp(i);
    Time t1 = epoch - Uvar(0.1,"s");
    Time t2 = epoch + Uvar(30*60,"s");
    Time t_mid = (t1 + t2)/2.0;
    TargetTimeAltitude tta(spacecraft_,target_,alt_tmp(i));
    Uvar final_time,final_altitude;
    golden_section_search(t1,t_mid,t2,accuracy,tta,final_time,final_altitude);
    time(i+N) = final_time-epoch;
    //cout<<"check result final time final alt "<<final_time-epoch<<" "<<final_altitude<<" "<<alt_tmp(i)<<endl;
    }

  Plot plot1;
  plot1.addXY(time,"min",angle,"deg",line("none"),sym("circle","red",2));
  plot1.setTitle("Data to be fit");
  //plot1.show("x");


  //fitting using double values
  // 1/angle(deg) vs time(min)
  //------------------------------------------
  //From experience, it is better fit the data with 6th power 
  // when lowest altitude > 1500 km
  // while with 8th power when lowest altitude < 1500 km
  //exception: if all input data were used, use 8th power
  //-------------------------------------------
 if(lowestAltitude() < Uvar(1500,"km")) angle_poly_coeff_time_.resize(9);
 else if( lowestAltitude() > Uvar(1500,"km")) angle_poly_coeff_time_.resize(7);
  else ErrorMessage("Flyby.cpp:lowest altitude no valid value").throwMe();

  Dvec time_in_min("time in min",2*N);
  Dvec angle_inv("1/angle",2*N);
  for (unsigned int i = 0; i < 2*N;++i)
    {
      time_in_min(i) = time(i).getInUnits("s")/60.0;
      angle_inv(i) = 1.0/(angle(i).getInUnits("rad")*radtodeg);
    }
  polyfit(angle_poly_coeff_time_,angle_inv,time_in_min);

  //cout<<"coeff "<<angle_poly_coeff_time_<<endl;
  //plot time vs angle-inv
  Plot plot_time_angle;
  plot_time_angle.addXY(time_in_min,angle_inv,line("solid","red",2),sym("none"));
  //plot_time_angle.show("x");


  Uvar start_time = startTime("ivd") - Uvar(15*60,"s");
  Uvar end_time = endTime("ivd") + Uvar(15*60,"s");
  Uvar time_step = interval("ivd");
  unsigned int Nmax = (unsigned int) ((end_time - start_time)/time_step).getInUnits("");
  Uvec angle_fit("angle fit",Nmax);
  Uvec time_fit("time fit",Nmax);
  for (unsigned int i = 0; i < Nmax;++i)
    {
    time_fit(i) = start_time+(end_time-start_time)*double(i)/double(Nmax-1);
    double t_in_min=time_fit(i).getInUnits("s")/60.0;
    double angle=0.0;
    for (unsigned int i_p = 0; i_p < angle_poly_coeff_time_.size();++i_p)
      {
      angle += angle_poly_coeff_time_(i_p)*pow(t_in_min,i_p);
      }
    angle_fit(i) = Uvar((1.0/angle)*degtorad,"rad");
    //cout<<"time and fit "<<time_fit(i)<<" "<<angle_fit(i).getInUnits("deg")<<endl;
    }
  Plot plot2;
  plot2.addXY(time,"min",angle,"deg",line("none"),sym("circle","red"));
  plot2.addXY(time_fit,"min",angle_fit,"deg",line("solid","red",2),sym("none"));
  plot2.setTitle("Data and fit");
  //plot2.show("x");

  //when speical target lat lon option is abled
  if(target_lat_lon_option_=="On")
    {
      //---------------------------------------
      //formular a + b*t + c*t^2+d*t^3+e*t^4
      // t: current time w.r.t. epoch(lowestaltitude)
      //t0:  time at the closest time to target w.r.t. epoch
      //start time and end time are values defined w.r.t. 
      //time_at_closest appraoch    
      //------------------------------------------

    double t0 = time_at_closest_range_to_target_.getInUnits("s")/60.0;//in min
    double t1 = start_target_lat_lon_tracking_.getInUnits("s")/60.0;
    double t2 = end_target_lat_lon_tracking_.getInUnits("s")/60.0;
    double t3 =  -fabs(off_centered_beam_time_lag_.getInUnits("s"))/60.0;
    double t4 =  fabs(off_centered_beam_time_lag_.getInUnits("s"))/60.0;

    //cout<<"time "<< t1<<" "<< t3<<" "<<t0 <<" "<<t4<<" "<<t2<<endl;

    unsigned int NN=7;
    target_lat_lon_fit_.resize(NN);
    Dmat A("matrix",NN,NN);
    Dvec angle_inv2("angle inv",NN);

    A(0,0) = 1; A(0,1) = 0; 
    A(0,2) = 0 ;A(0,3) = 0;A(0,4) = 0;
    A(0,5)=0; A(0,6)=0;

    A(1,0) = 1; A(1,1) = t1;
    A(1,2) = t1*t1;A(1,3)=t1*t1*t1;A(1,4)=t1*t1*t1*t1;
    A(1,5)=pow(t1,5); A(1,6)=pow(t1,6);

    A(2,0) = 1; A(2,1) = t2;
    A(2,2) = t2*t2;A(2,3)=t2*t2*t2;A(2,4)=t2*t2*t2*t2;
    A(2,5)=pow(t2,5); A(2,6)=pow(t2,6);

    A(3,0) = 0; A(3,1) = 1;A(3,2) = 2.0*t1;
    A(3,3)= 3.0*t1*t1;A(3,4) = 4.0*t1*t1*t1;
    A(3,5)=5.0*pow(t1,4);A(3,6)=6.0*pow(t1,5);

    
    A(4,0) = 0; A(4,1) = 1;A(4,2) = 2.0*t2;
    A(4,3)= 3.0*t2*t2;A(4,4) = 4.0*t2*t2*t2;
    A(4,5) = 5.0*pow(t2,4); A(4,6)= 6.0*pow(t2,5);

    A(5,0) = 1;       A(5,1) = t3; 
    A(5,2) = t3*t3;   A(5,3) = t3*t3*t3; A(5,4) = t3*t3*t3*t3;
    A(5,5)=pow(t3,5); A(5,6)=pow(t3,6);
    
    A(6,0) = 1;       A(6,1) = t4; 
    A(6,2) = t4*t4;   A(6,3) = t4*t4*t4; A(6,4) = t4*t4*t4*t4;
    A(6,5)=pow(t4,5); A(6,6)=pow(t4,6);

    A.inv(); 
    //--------------------------
    //matrix
    //----------------------------
    //for (unsigned int i = 0; i < 5;++i)
    //{
    //for (unsigned int j = 0; j < 5;++j)
    //{
    //  cout<<A(i,j)<<" ";
    //}
    //cout<<endl;
    //}
   
    //for (unsigned int i = 0; i < 5;++i)
    //{
    //for (unsigned int j = 0; j < 5;++j)
    //  {
    //    cout<<A(i,j)<<" ";
    //  }
    //cout<<endl;
    //}

    angle_inv2 = 0.0;//reset
    angle_inv2(0) = 1.0/angle_at_closest_range_to_target_.getInUnits("deg");
    //cout<<"inverse angle "<<angle_inv2(0)<<endl;

    for (unsigned int i_p = 0; i_p < angle_poly_coeff_time_.size();++i_p)
      {//1/angle at t1 and t2
	angle_inv2(1) += angle_poly_coeff_time_(i_p)*pow(t0+t1,i_p);
	angle_inv2(2) += angle_poly_coeff_time_(i_p)*pow(t0+t2,i_p);
      }
    for (unsigned int i_p = 1; i_p < angle_poly_coeff_time_.size();++i_p)
      {// derivative
	angle_inv2(3) += angle_poly_coeff_time_(i_p)*pow(t0+t1,i_p-1)*double(i_p);
	angle_inv2(4) += angle_poly_coeff_time_(i_p)*pow(t0+t2,i_p-1)*double(i_p);
      }

    //---------------------
    //hold the same angle
    //---------------------
    angle_inv2(5) = angle_inv2(0);
    angle_inv2(6) = angle_inv2(0);
    
    //for(unsigned int i=0;i<NN;++i){
    //cout<<"targeted angle "<< 1.0/angle_inv2(i)<<endl;
    //}

    target_lat_lon_fit_ = 0.0;
    for (unsigned int i = 0; i < NN;++i)
    for (unsigned int i_p = 0; i_p<NN;++i_p)
      {  
	target_lat_lon_fit_(i)+= A(i,i_p) * angle_inv2(i_p);
      }
    
    Uvar start_time = startTime("ivd") ;
    Uvar end_time = endTime("ivd") ;
    Uvar time_step = Uvar(1,"s");
    unsigned int Nmax = (unsigned int)
      ((end_time - start_time)/time_step).getInUnits("");

    Uvec angle_fit("angle fit",Nmax);
    Uvec time_fit("time fit",Nmax);

    for (unsigned int i = 0; i < Nmax;++i)
      {
      time_fit(i) = start_time+(end_time-start_time)*double(i)/double(Nmax-1);
      
      if(time_fit(i)<(start_target_lat_lon_tracking_+time_at_closest_range_to_target_)
	||time_fit(i)>(end_target_lat_lon_tracking_+time_at_closest_range_to_target_))
	{//normal fit
	double t_in_min=time_fit(i).getInUnits("s")/60.0;
	double angle=0.0;
	for (unsigned int i_p = 0; i_p < angle_poly_coeff_time_.size();++i_p)
	  {
	  angle += angle_poly_coeff_time_(i_p)*pow(t_in_min,i_p);
	  }
	angle_fit(i) = Uvar((1.0/angle)*degtorad,"rad");
	}
      else
	{//target lat_lon fit
	double t_in_min=time_fit(i).getInUnits("s")/60.0;
	double angle=0.0;
	for (unsigned int i_p = 0; i_p < NN;++i_p)
	  {
	  angle += target_lat_lon_fit_(i_p)*pow(t_in_min-t0,i_p);
	  }
	angle_fit(i) = Uvar((1.0/angle)*degtorad,"rad");
	//cout<<"time and angle "<<t_in_min<<" "<<angle_fit(i).getInUnits("deg")<<endl;
	}
      }
    
    //Plot plot3;
    //plot3.addXY(time_fit,"min",angle_fit,"deg",line("solid","red",2),sym("none"));
    //plot3.setTitle("Target lat lon");
    //plot3.show("x");
    }
  computePolyfitCoeff_set_ = true;
  }

//------------------------------------
// Incidence Angle at time t 
// Profile option was read in config method
//-----------------------------------
void Flyby::load_local_incidence_angle_offset_params(
	     const vector<Uvar> local_incidence_angle_offset,
	     const vector<Uvar> local_incidence_angle_offset_tr1,
	     const vector<Uvar> local_incidence_angle_offset_tr2,
	     const vector<Uvar> local_incidence_angle_offset_tr3,
	     const vector<Uvar> local_incidence_angle_offset_tr4)
{
  Noffset_=local_incidence_angle_offset.size();
  if(Noffset_ == 0) return;//do not bother
  local_offset_set_=true;
  angle_offset_.resize(Noffset_);
  tr1_.resize(Noffset_);
  tr2_.resize(Noffset_);
  tr3_.resize(Noffset_);
  tr4_.resize(Noffset_);
  tr12_mid_.resize(Noffset_);
  tr34_mid_.resize(Noffset_);
  
  for(unsigned int i=0;i<Noffset_;++i){
    angle_offset_(i)=local_incidence_angle_offset[i];
    tr1_(i)=local_incidence_angle_offset_tr1[i];
    tr2_(i)=local_incidence_angle_offset_tr2[i];
    tr3_(i)=local_incidence_angle_offset_tr3[i];
    tr4_(i)=local_incidence_angle_offset_tr4[i];
    
    tr12_mid_(i)= (tr1_(i) + tr2_(i))/2.0;
    tr34_mid_(i)= (tr3_(i) + tr4_(i))/2.0;
  }

  if(!computePolyfitCoeff_set_ ) ErrorMessage("Flyby.cpp::load_local_incidence_angle_offset_params: incidence angle polyfitcoeff is not done yet").throwMe();
  
  Uvar time_step=Uvar(2,"s");
  unsigned int NN=4;

  tr_from_ideal_.resize(Noffset_,4);
  tr_to_ideal_.resize(Noffset_,4);

  //transition off from nominal profile
  for (unsigned int i=0;i<Noffset_;++i){
    //let's compute polynomial fitting
   
    Dmat A("matrix",NN,NN);
    double t1= tr1_(i).getInUnits("min") - tr12_mid_(i).getInUnits("min");
    double t2= tr2_(i).getInUnits("min") - tr12_mid_(i).getInUnits("min");
      
    Dvec Y("angle value",NN);
    Y(0)= desiredIncidenceAngle(epoch_time_ + tr1_(i)).getInUnits("deg");
    Y(1)= desiredIncidenceAngle(epoch_time_ + tr2_(i)).getInUnits("deg")+angle_offset_(i).getInUnits("deg");

    //cout<<"targeted angle "<< Y(0)<<" "<<Y(1)<<endl;
    
    Y(2)=  desiredIncidenceAngle(epoch_time_ + tr1_(i) + time_step).getInUnits("deg");
    Y(2)-=  desiredIncidenceAngle(epoch_time_ + tr1_(i) - time_step).getInUnits("deg");
    Y(2) /= 2.0*time_step.getInUnits("min");

    Y(3) =  desiredIncidenceAngle(epoch_time_ + tr2_(i) + time_step).getInUnits("deg");
    Y(3) -=  desiredIncidenceAngle(epoch_time_ + tr2_(i) - time_step).getInUnits("deg");
    Y(3) /= 2.0*time_step.getInUnits("min");




    A(0,0) = 1;
    A(0,1) = t1;
    A(0,2) = t1*t1;
    A(0,3)=t1*t1*t1;
   
    
    A(1,0) = 1; 
    A(1,1) = t2;
    A(1,2) = t2*t2;
    A(1,3)=t2*t2*t2;
  

    A(2,0) = 0; 
    A(2,1) = 1;
    A(2,2) = 2.0*t1;
    A(2,3)= 3.0*t1*t1;
   
    
    A(3,0) = 0;
    A(3,1) = 1;
    A(3,2) = 2.0*t2;
    A(3,3)= 3.0*t2*t2;
   
  
    A.inv(); 

    for(unsigned int j=0;j<NN;++j){
      tr_from_ideal_(i,j)=0.0;

      for(unsigned int k=0;k<NN;++k)
	tr_from_ideal_(i,j) += A(j,k)*Y(k);
    }

  }  


  //transition to nomianl profile
  for (unsigned int i=0;i<Noffset_;++i){
    //let's compute polynomial fitting
    unsigned int NN=4;
    Dmat A("matrix",NN,NN);
    double t1= tr3_(i).getInUnits("min") - tr34_mid_(i).getInUnits("min");
    double t2= tr4_(i).getInUnits("min") - tr34_mid_(i).getInUnits("min");

    Dvec Y("angle value",NN);
    Y(0)= desiredIncidenceAngle(epoch_time_ + tr3_(i)).getInUnits("deg")+angle_offset_(i).getInUnits("deg");
    Y(1)= desiredIncidenceAngle(epoch_time_ + tr4_(i)).getInUnits("deg");

    //cout<<"targeted angle "<< Y(0)<<" "<<Y(1)<<endl;
    
    Y(2)=  desiredIncidenceAngle(epoch_time_ +tr3_(i) + time_step).getInUnits("deg");
    Y(2)-=  desiredIncidenceAngle(epoch_time_ + tr3_(i) - time_step).getInUnits("deg");
    Y(2) /= 2.0*time_step.getInUnits("min");

    Y(3) =  desiredIncidenceAngle(epoch_time_ + tr4_(i) + time_step).getInUnits("deg");
    Y(3) -=  desiredIncidenceAngle(epoch_time_ + tr4_(i) - time_step).getInUnits("deg");
    Y(3) /= 2.0*time_step.getInUnits("min");

    A(0,0) = 1;
    A(0,1) = t1;
    A(0,2) = t1*t1;
    A(0,3)=t1*t1*t1;
    
    A(1,0) = 1; 
    A(1,1) = t2;
    A(1,2) = t2*t2;
    A(1,3)=t2*t2*t2;
  
    A(2,0) = 0; 
    A(2,1) = 1;
    A(2,2) = 2.0*t1;
    A(2,3)= 3.0*t1*t1;
   
    A(3,0) = 0;
    A(3,1) = 1;
    A(3,2) = 2.0*t2;
    A(3,3)= 3.0*t2*t2;
   
    A.inv(); 

    for(unsigned int j=0;j<NN;++j){
      tr_to_ideal_(i,j)=0.0;

      for(unsigned int k=0;k<NN;++k)
	tr_to_ideal_(i,j) += A(j,k)*Y(k);
    }
  }  
  

  //debug plot
  /*  
  unsigned int index=0;
  unsigned int Nstep= (unsigned int) round_double( ( ( tr2_(index) - tr1_(index)) / time_step).getInUnits(""));
  Uvec angle("",Nstep),time("",Nstep);
  for(unsigned int i=0;i<Nstep;++i){
    time(i)= tr1_(index)+float(i)*time_step;
    double t= ( time(i) - tr12_mid_(index)).getInUnits("min");
    double theta=0.0;
    for(unsigned int j=0;j<NN;++j)
      theta += tr_from_ideal_(index,j)*pow(t,j);
    angle(i)=Uvar(theta*pi/180.0,"rad");
  }
  Plot A;
  A.addXY(time,"min",angle,"deg",line("none"),sym("circle","red",1));
  A.show("x");

  Nstep= (unsigned int) round_double( ( ( tr4_(index) - tr3_(index)) / time_step).getInUnits(""));
  for(unsigned int i=0;i<Nstep;++i){
    time(i)= tr3_(index)+float(i)*time_step;
    double t= ( time(i) - tr34_mid_(index)).getInUnits("min");
    double theta=0.0;
    for(unsigned int j=0;j<NN;++j)
      theta += tr_to_ideal_(index,j)*pow(t,j);
    angle(i)=Uvar(theta*pi/180.0,"rad");
  }
  Plot B;
  B.addXY(time,"min",angle,"deg",line("none"),sym("circle","red",1));
  B.show("x");
  */


}

Uvar Flyby::desiredIncidenceAngle(const Time& t)
  {
    if (angle_profile_option_ == "sech_function")
      {
      StateVector sc_state("sc_state");
      Frame ftitan("IAU_TITAN",target_);
      ftitan.ephemeris(sc_state,"Cassini",t,"NONE");
      TargetGeom tg(t);
      tg.setState(sc_state);
      tg.setTarget(target_,ftitan);
      Uvar alt = tg.altitude();
      Uvar  x1 =    angle_at_ca_- residual_angle_;
      x1 /= (angle_at_mid_alt_ -residual_angle_);//unitless
      if (x1 * x1 < Uvar(1.0," ")) 
	{
	  throw ErrorMessage("Angle profile option: no solution to meet request angle profile");
	}
      Uvar y = x1 + sqrt( (x1*x1)-1.0 );//unitless
      Uvar sech_function_coeff = (mid_alt_ - lowestAltitude())/log(y);
      Uvar desired_incidence_angle= residual_angle_;
      Uvar angle_mid_step = 2.0*(angle_at_ca_-residual_angle_);
      Uvar x =(alt - lowestAltitude())/sech_function_coeff;
      angle_mid_step /=(exp(x) + exp(-x));
      desired_incidence_angle += angle_mid_step;
      return(desired_incidence_angle);
      }    
    else if ( angle_profile_option_ =="constant_incidence_angle")
      {
	return(angle_polyfit_(0));
      }
    else if (angle_profile_option_ == "polynomial_fit")
      {
      if (!computePolyfitCoeff_set_) ErrorMessage("No coeff fitting was done").throwMe();

      Uvar t_wrt_epoch = t - lowestAltitudeTime();
      double t_in_min = t_wrt_epoch.getInUnits("s")/60.0;//min
      double angle_inv=0.0;
      if(t_wrt_epoch<=
	 (start_target_lat_lon_tracking_+time_at_closest_range_to_target_)
	 ||t_wrt_epoch>=
	 (end_target_lat_lon_tracking_+time_at_closest_range_to_target_))
	{//normal fit
	  for (unsigned int i_p = 0; i_p < angle_poly_coeff_time_.size();++i_p)
	  {
	  angle_inv += angle_poly_coeff_time_(i_p)*pow(t_in_min,i_p);
	  }
	}
      else
	{//target lat_lon fit
	  double t0=time_at_closest_range_to_target_.getInUnits("s")/60.0;//in min
	  for (unsigned int i_p = 0; i_p < target_lat_lon_fit_.size();++i_p)
	  {
	  angle_inv += target_lat_lon_fit_(i_p)*pow(t_in_min-t0,i_p);
	  }
	}
      return(Uvar((1.0/angle_inv)*degtorad,"rad"));
      }
    else
      {
      cout<<"Flyby::incidenceAngle(t): no valid option is chosen, so 20 deg incidence angle will be returned"<<endl;
      }
    return(Uvar(20*degtorad,"rad"));
  }

//------------------
//Return incidence angle with an addition of local bias
//-------------------
Uvar Flyby::desiredIncidenceAngle_with_local_bias(const Time& t){

  if(!local_offset_set_)
    return(desiredIncidenceAngle(t));
  
  Uvar time=t - epoch_time_;
  bool local_bias=false;
  unsigned int index=0;
  for(unsigned int i=0; i<Noffset_;++i){
    if(time>=tr1_(i) && time<=tr4_(i)){
      index=i;
      local_bias=true;
      break;
    }
  }

  if(!local_bias) return(desiredIncidenceAngle(t));

  Uvar value;
  unsigned int MM,NN;
  tr_from_ideal_.size(MM,NN);

  if(time>=tr1_(index) && time<=tr2_(index)){
    double t2= ( time - tr12_mid_(index)).getInUnits("min");
    double theta=0.0;
    for(unsigned int j=0;j<NN;++j)
      theta += tr_from_ideal_(index,j)*pow(t2,j);
    value=Uvar(theta*pi/180.0,"rad");
  }
  
  else if( time>tr2_(index) && time <tr3_(index)){
    value= desiredIncidenceAngle(t);
    value += angle_offset_(index);
  }
  else if(time>=tr3_(index) && time <=tr4_(index)){
    double t2= ( time - tr34_mid_(index)).getInUnits("min");
    double theta=0.0;
    for(unsigned int j=0;j<NN;++j)
      theta += tr_to_ideal_(index,j)*pow(t2,j);
    value=Uvar(theta*pi/180.0,"rad");
  }
  else{
    ErrorMessage("flyby.cpp::desiredincidenceangle_with_local_bias: time segmentation error").throwMe();
  }
  return(value);
}

//------------------------------------------------------------------
// resetTimeSteps()
//
// Set current time counter back to start time.
//------------------------------------------------------------------

void Flyby::resetTimeSteps()
  {
  time_counter_ = 0;
  cur_time_ = epoch_time_ + Time(start_time_);
  percent_completed_ = 0;
  }

void Flyby::resetTimeSteps(const string& suffix)
  {
  suffixCheck(suffix);
  time_counter_map_[suffix] = 0;
  cur_time_map_[suffix] = epoch_time_ +
    Time(start_time_map_[suffix]);
  percent_completed_map_[suffix] = 0;
  }

//------------------------------------------------------------------
// t = currentTime();
//
// Return current time step.  Also increment to next time step.
//------------------------------------------------------------------

Time Flyby::currentTime()
  {
  if (time_counter_ >= ntime_)
    {
    ErrorMessage e("Exceeded number of Flyby time steps");
    e.throwMe();
    }
  Time t = cur_time_;
  time_counter_++;
  cur_time_ += interval_;
  return(t);
  }

Time Flyby::currentTime(const string& suffix)
  {
  suffixCheck(suffix);
  if (time_counter_map_[suffix] >= ntime_map_[suffix])
    {
    ErrorMessage e("Exceeded number of Flyby time steps for suffix: " + suffix);
    e.throwMe();
    }
  Time t = cur_time_map_[suffix];
  time_counter_map_[suffix]++;
  cur_time_map_[suffix] += interval_map_[suffix];
  return(t);
  }

//-------------------------------------------------------------------------
// percentCompletedDisplay(stream,percent_step)
//
// Show the percentage of time steps completed at percent_step intervals
// on the supplied stream.
//-------------------------------------------------------------------------

void Flyby::percentCompletedDisplay(ostream& s, unsigned int percent_step)
  {
  if (time_counter_ % ((ntime_*percent_step)/100) == 0)
    {
    s << percent_completed_ << "%" << endl;
    percent_completed_ += percent_step;
    }
  }
  
void Flyby::percentCompletedDisplay(ostream& s, unsigned int percent_step,
  const string& suffix)
  {
  suffixCheck(suffix);
  if (time_counter_map_[suffix] % ((ntime_map_[suffix]*percent_step)/100) == 0)
    {
    s << percent_completed_map_[suffix] << "%" << endl;
    percent_completed_map_[suffix] += percent_step;
    }
  }

//------------------------------
// Internal support methods
//------------------------------

//--------------------------------------------------------------------
// value = timeValueSubstitute(cfg,keyword)
//
// Return a numeric value for the supplied keyword after making any
// defined string substitutions for the time value.  The recognized
// substitutions are defined in this function.  The value is always
// returned epoch relative.  This function assumes that epoch_time_,
// ck_earliest_start_time_, and ck_latest_end_time_ are all set before
// it is called.  Note that a valid UTC time string is a recognized
// value.
//--------------------------------------------------------------------

Uvar Flyby::timeValueSubstitute(Config& cfg, const string& keyword)
  {
  Uvar value;
  if (cfg.numberKeywordExists(keyword))
    {
    value = cfg[keyword];
    }
  else
    {
    string str = cfg.str(keyword);
    if (str == "ckernel_start")
      {
      value = ck_earliest_start_time_ - epoch_time_;
      }
    else if (str == "ckernel_end")
      {
      value = ck_latest_end_time_ - epoch_time_;
      }
    else if (isUtc(str))
      {
      Time t(str);
      value = t - epoch_time_;
      }
    else
      {
      ErrorMessage e("Invalid keyword, string value pair: (" +
        keyword + " " + str + ")");
      e.throwMe();
      }
    }
  return(value);
  }
//------------------------
// Class Orbit Methods
//------------------------

//------------------
// Automatic testing
//------------------

bool Orbit::selfTest()
  {
  return(true);
  }

//-------------
// Constructors
//-------------

//--------------------------------------------------------
// Orbit(target_name,state)
//
// Establish an orbit/trajectory for a spacecraft of insignificant
// mass and the named target body with the supplied StateVector
// providing one position velocity pair on the trajectory.
// Other points will be propagated from this one.
//--------------------------------------------------------

Orbit::Orbit(const string& target_name, const StateVector& state)
  {
  setTarget(target_name);
  setState(state);
  }

//-------------
// Setup
//-------------

//--------------------------------------------------------
// setTarget(target_name)
//
// Specifies a new target body name to compute a trajectory around.
// Calling this method will reset this Orbit object so that future
// calls for state vectors will require propagation.
//--------------------------------------------------------

void Orbit::setTarget(const string& target_name)
  {
  target_name_ = target_name;
  ftarget_I_ = Frame("J2000",target_name);
  if (states_.size() > 1)
    {  // Any propagated states for the prior target have to be discarded
    states_.clear();
    }
  // Fetch target gravitational constant from SPICE kernel
  spice_target_value(target_name,"GM",GM_);
  target_set_ = true;
  }

//--------------------------------------------------------
// setState(state)
//
// Specifies a defining point on a trajectory.
// Calling this method will reset this Orbit object so that future
// calls for state vectors will require propagation.
//--------------------------------------------------------

void Orbit::setState(const StateVector& state)
  {
  // Note that setTarget is only required once.  The previous target
  // will be used if setTarget is not called again.
  states_.clear();
  states_.push_back(state);  // first state is the defining point
  }

//-------------------
// Orbit Propagation
//-------------------

//------------------------------------------------------------------
// state = getState(t)
//
// Returns the spacecraft statevector at time t.
// Automatically performs any needed propagation and interpolation.
//------------------------------------------------------------------

StateVector Orbit::getState(const Time& t)
  {
  if (states_.size() == 0)
    {
    ErrorMessage e("Need to define orbit before requesting a StateVector");
    e.throwMe();
    }
  if (!target_set_)
    {
    ErrorMessage
      e("Need to define orbit target before requesting a StateVector");
    e.throwMe();
    }

  StateVector return_state("return_state");

  if (t < states_.front().time())
    {  // need to propagate backwards in time
    StateVector s_target;
    
    Time front_time = states_.front().time();
    unsigned int N =
      1 + (unsigned int)((front_time - t)/interval_).getInUnits("");

    PositionVector r(states_.front().position());
    FloatVector v(states_.front().velocity());
    for (unsigned int i=0; i < N; ++i)
      {
      // Two body separation at original front time
      ftarget_I_.ephemeris(s_target,target_name_,front_time,"NONE");
      PositionVector rr = s_target.position() - r;
      // Update position based on position and velocity at original front time
      r -= v*interval_;
      // Update velocity based on velocity and Newtons law of gravitation for
      // the two body separation at the original front time
      v -= rr*GM_*interval_/pow(rr.magnitude(),3);  // update velocity
      // Update times, build statevector, and add to list
      front_time -= interval_;
      r.setTime(front_time);
      v.setTime(front_time);
      StateVector new_state(r,v);
      states_.push_front(new_state);
      }

    // Interpolate first interval to requested time
    STATELIST::iterator p1 = states_.begin();
    STATELIST::iterator p2 = p1;
    p2++;
    return_state = *p2 - (p2->time()-t)/(p2->time() - p1->time())*(*p2 - *p1);
    }
  else if (t > states_.back().time())
    {  // need to propagate forwards in time
    StateVector s_target;

    Time back_time = states_.back().time();
    unsigned int N = 1 +
      (unsigned int)((t - back_time)/interval_).getInUnits("");

    PositionVector r(states_.back().position());
    FloatVector v(states_.back().velocity());
    for (unsigned int i=0; i < N; ++i)
      {
      // Two body separation at original back time
      ftarget_I_.ephemeris(s_target,target_name_,back_time,"NONE");
      PositionVector rr = s_target.position() - r;
      // Update position based on position and velocity at original back time
      r += v*interval_;
      // Update velocity based on velocity and Newtons law of gravitation for
      // the two body separation at the original back time
      v += rr*GM_*interval_/pow(rr.magnitude(),3);  // update velocity
      // Update times, build statevector, and add to list
      back_time += interval_;
      r.setTime(back_time);
      v.setTime(back_time);
      StateVector new_state(r,v);
      states_.push_back(new_state);
      }

    // Interpolate last interval to requested time
    STATELIST::iterator p2 = states_.end();
    STATELIST::iterator p1 = states_.end();
    p2--;  // backup to last member
    p1--;  // backup to last member
    p1--;  // backup to 2nd to last member
    return_state = *p2 - (p2->time()-t)*(*p2 - *p1)/(p2->time() - p1->time());
    }
  else
    {  // requested time is in range
    STATELIST::iterator p2;
    for (p2=states_.begin(); p2 != states_.end(); ++p2)
      {  // loop over states until we find the two that bracket the desired t
      if (t <= p2->time()) break;
      }
    STATELIST::iterator p1 = p2;
    p1--;
    // Interpolate between current and prior in list
    return_state = *p2 - (p2->time()-t)*(*p2 - *p1)/(p2->time() - p1->time());
    }

  return(return_state);
  }

//-----------------------
// Orbit Characteristics
//-----------------------

//----------------------------------------------------------------------------
// f = bPlaneFrame()
//
// Returns the B-plane frame.
// The B-plane is the plane perpendicular to the asymptotic velocity vector
// on approach, and which contains the target body center point.
// The B-plane frame has its origin at the point where the trajectory
// intercepts the B-plane.  The x-axis points away from the target center.
// The y-axis is orthogonal to the x-axis and lies in the B-plane.
// Positive y is chosen so that the resulting z-axis points along the
// direction of motion.
//----------------------------------------------------------------------------

Frame Orbit::bPlaneFrame()
  {
  PositionVector r_intercept = bPlaneIntercept().position();
  r_intercept.representIn(ftarget_I_);
  DirectionVector v_asymptote_dir = asymptoticApproachVelocity();
  v_asymptote_dir.representIn(ftarget_I_);
  DirectionVector x("x",r_intercept);
  DirectionVector y("y",cross(v_asymptote_dir,r_intercept));
  DirectionVector z("z",cross(x,y));
  return(Frame("bplane_frame",r_intercept,x,y,z));
  }

StateVector Orbit::bPlaneIntercept()
  {
  StateVector smax = lowestAltitudeState();
  return(smax);
  }

FloatVector Orbit::asymptoticApproachVelocity()
  {
  return(FloatVector("not done",ftarget_I_,lowest_altitude_time_,0,0,0));
  }

//------------------------------------------------------------------
// a = lowestAltitude();
//
// Compute lowest altitude.
// Note that this method pulls spacecraft position from *this Orbit
//------------------------------------------------------------------

Uvar Orbit::lowestAltitude()
  {
  if (lowest_altitude_set_)
    {
    return(lowest_altitude_);
    }

  //-------------------------------------------------------------
  // Locate nearest closest approach using golden section search.
  //-------------------------------------------------------------

  Time t1 = states_.front().time();
  Uvar t_mid("t_mid");
  t_mid = t1 + Time(Uvar(0.001,"s"));
  Uvar t2("t2");
  Uvar final_time("final_time");
  TargetAltitude tga(this,target_name_);
  bracket_minimum(t1,t_mid,t2,tga);
  golden_section_search(t1,t_mid,t2,interval_,tga,
    final_time,lowest_altitude_);
  lowest_altitude_time_.setEt(final_time);
  lowest_altitude_set_ = true;

  return(lowest_altitude_);
  }

//------------------------------------------------------------------
// t = lowestAltitudeTime();
//
// Compute the time of lowest altitude.
// Note that this method pulls spacecraft position from *this Orbit
//------------------------------------------------------------------

Uvar Orbit::lowestAltitudeTime()
  {
  if (!lowest_altitude_set_)
    {
    lowestAltitude();
    }

  return(lowest_altitude_time_);
  }

//------------------------------------------------------------------
// s = lowestAltitudeState();
//
// Return the spacecraft statevector at the point of closest approach.
// Note that this method pulls spacecraft position from *this Orbit
//------------------------------------------------------------------

StateVector Orbit::lowestAltitudeState()
  {
  if (!lowest_altitude_set_)
    {
    lowestAltitude();
    }
  return(getState(lowest_altitude_time_));
  }

//-----------------------
// Static data for Orbit
//-----------------------

const Uvar Orbit::interval_(1,"s");
