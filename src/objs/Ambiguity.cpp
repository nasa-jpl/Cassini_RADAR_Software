//==========================================================
// calculate  surface geomtric factors at time t and sc attitude
//==========================================================
//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Frame.h"
#include "TargetGeom.h"
#include "Units.h"
#include "Time.h"
#include "Constants.h"
#include "TemplateUtils.h"
#include "Ambiguity.h"
#include "Config.h"
#include "Array.h"
#include "Io.h"
#include "Units.h"
#include "Utils.h"

//-------------------------------------------------------------
//Calculate ambiguity at t for a beam 
// need to pass Beam information later
//-------------------------------------------------------------

using std::endl;
using std::cout;

Ambiguity::Ambiguity(const Time& t)
  :range("range"),
   doppler("dop"),
   oneway_beam_gain("oneway_gain"),
   incidenceangle("incidence_angle"),  
   area("area"),
   alongtrack("alongtrack"),
   crosstrack("crosstrack"),
   radar_geom_factor("radar_geom_factor"),
   radar_geom_factordB("radar_geom_factordB"),
   number_of_solution("number of solution"),
   oneway_beam_gain_mirror("oneway_gain"),
   incidenceangle_mirror("incidence_angle"),  
   area_mirror("area"),
   radar_geom_factor_mirror("radar_geom_factor_mirror"),
   number_of_solution_mirror("mirror number of solution"),
   good_bin("good_bin"),
   ML_good_bin("ML_good_bin"),
   amb_ratio_SL("amb_ratio_SL"),
   amb_ratio_ML("amb_ratio_ML"),
   thermal_snr_SL("thermal_snr_SL"),
   thermal_snr_ML("thermal_snr_ML"),
   time_set_(true),
   target_set_(false),
   beam_set_(false),
   ieb_set_(false),
   config_set_(false),
   mirror_geom_set_(false),
   beam_number_set_(false),
   beam_max_gain_set_(false),
   grid_number_set_(false),
   radar_geom_factor_set_(false),
   peak_radar_geom_factor_set_(false),
   usable_area_set_(false),
   ML_usable_area_set_(false),
   backscatt_model_set_(false),
   process_window_set_(false),
   track_frame_set_(false),
   t_(t),
   nlook_(1)
  {
  }

Ambiguity::Ambiguity()
  :range("range"),
   doppler("dop"),
   oneway_beam_gain("oneway_gain"),
   incidenceangle("incidence_angle"),
   area("area"),
   alongtrack("alongtrack"),
   crosstrack("crosstrack"),
   radar_geom_factor("radar_geom_factor"),
   radar_geom_factordB("radar_geom_factordB"),
   number_of_solution("number of solution"),
   oneway_beam_gain_mirror("oneway_gain"),
   incidenceangle_mirror("incidence_angle"),
   area_mirror("area"),
   radar_geom_factor_mirror("radar_geom_factor_mirror"),
   number_of_solution_mirror(" mirror number of solution"),
   good_bin("good_bin"),
   ML_good_bin("ML_good_bin"),
   amb_ratio_SL("amb_ratio_SL"),
   amb_ratio_ML("amb_ratio_ML"),
   thermal_snr_SL("thermal_snr_SL"),  
   thermal_snr_ML("thermal_snr_ML"),
   time_set_(false),
   target_set_(false),
   beam_set_(false),
   ieb_set_(false),
   config_set_(false),
   mirror_geom_set_(false),
   beam_number_set_(false),
   beam_max_gain_set_(false),
   grid_number_set_(false),
   radar_geom_factor_set_(false),
   peak_radar_geom_factor_set_(false),
   usable_area_set_(false),
   ML_usable_area_set_(false),
   backscatt_model_set_(false),
   process_window_set_(false),
   track_frame_set_(false),
   nlook_(1)
 
  { 
  }


//-----------------------
//config(): read in number of grid
//-----------------------
void Ambiguity::config(Config& cfg)
{
  Uvar x_tmp = cfg["Number_of_range_doppler_grid"];
  Ndat_ = (unsigned int) x_tmp.getValue();
  Ngrid_ = Ndat_ + 1;
  lambda_ = speed_light/carrier_frequency;
  max_prf_ = cfg["ambiguity_search_max_prf"];
  min_prf_ = cfg["ambiguity_search_min_prf"];
  roll_direction_ = cfg.str("look_direction");
 

  string mirror_ambiguity_option=cfg.str("calculate_mirror_ambiguity");
  if (mirror_ambiguity_option =="Yes" || mirror_ambiguity_option=="yes")
    {
    mirror_geom_set_ = true;
    }
  else
    {
    mirror_geom_set_ = false;
    }
  grid_number_set_ = true;
  config_set_ = true;
}

//-----------------------
//{get,set}Time(const Time& t)
//----------------------
Time  Ambiguity::getTime()
  {
  if(!time_set_)
    {
    ErrorMessage e("No time is set");
    e.throwMe();
    }
  return(t_);
  }

void  Ambiguity::setTime(const Time& t)
  {
  if(time_set_)
    {
    ErrorMessage e("time is already set");
    e.throwMe();
    }
  t_ = t;
  time_set_ = true;
  }


//--------------------------------
//setTarget(const string& target_name, const Frame& target_frame)
//--------------------------------
void Ambiguity::setTarget(const string& target_name,const Frame& target_frame)
  {
  if (!time_set_){ErrorMessage("Ambiguity::setTarget: no time is set").throwMe();}
  target_name_ = target_name;
  target_frame_ = target_frame;
  target_set_ = true;
  target_frame_.ephemeris(state_,"Cassini",t_,"NONE");
  }
 
//--------------------------------
//setBeam(const Beam& beam)
//--------------------------------
void Ambiguity::setBeam(const Beam& beam)
  {
  beam_ = beam;
  beam_num_ = beam_.getBeamNumber();
  beam_max_gain_ = beam_.getMaxGain();
  beam_set_ = true;
  beam_number_set_ = true;
  beam_max_gain_set_ = true;
  }
//---------------------------
//{get,set}Beamnumber()
//---------------------------
unsigned int Ambiguity::getBeamNumber()
  {
  if(!beam_number_set_)
    {
    ErrorMessage e("Ambiguity::getBeamnumber: No beam number  is set");
    e.throwMe();
    }
  return(beam_num_);
  }

void Ambiguity::setBeamNumber(const unsigned int& beam_num)
  {
  if(beam_number_set_)
    {
    ErrorMessage e("Ambiguity::getBeamnumber: Beam number is already set");
    e.throwMe();
    }
  beam_num_ = beam_num;
  beam_number_set_ = true;
  }

//--------------------------------
//{get,set}beam's max gain
//--------------------------------
double Ambiguity::getBeamMaxGain()
  {
  if (!beam_max_gain_set_)
    {
    ErrorMessage e("Ambiguity::getMaxGain: Beam max gain is not set");
    e.throwMe();
    }
  return(beam_max_gain_);
  }

void  Ambiguity::setBeamMaxGain(const double& beam_max_gain)
  {
  if (beam_max_gain_set_)
    {
    ErrorMessage e("Ambiguity::setMaxGain: Beam is already set");
    e.throwMe();
    }
  beam_max_gain_ = beam_max_gain;
  beam_max_gain_set_ = true;
  }

//--------------------------
//{set} track frame
//--------------------------
void Ambiguity::setTrackFrame(const Frame& ftrack)
  {
    if(track_frame_set_)
      {
	ErrorMessage("Track frame is set").throwMe();
      }
    ftrack_ = ftrack;
    track_frame_set_ = true;
  }

//-----------------------------
//{set,get} carrier wavelength
//-----------------------------
void Ambiguity::setWavelength(const Uvar& lambda)
  {
    lambda_ = lambda;
  }

Uvar Ambiguity::getWavelength()
  {
    if(!config_set_)
      {
	ErrorMessage("Ambiguity::getWavelength: no config parameters are set").throwMe();
      }
    return(lambda_);
  }

//------------------------------------------
//{get,set}Gridsize(): actual data size Ndat_
//--------------------------------------------
unsigned int Ambiguity::getGridsize()
  {
  if (!grid_number_set_)
    {
    ErrorMessage e("No grid numer is specified");
    e.throwMe();
    }
  return(Ndat_);
  }

void  Ambiguity::setGridsize(const unsigned int& Ndat)
  {
  if (grid_number_set_)
    {
    ErrorMessage e("Grid size is already set");
    e.throwMe();
    }
  Ndat_ = Ndat;
  range.resize(Ndat_);
  doppler.resize(Ndat_);
  oneway_beam_gain.resize(Ndat_,Ndat_);
  incidenceangle.resize(Ndat_,Ndat_);
 
  area.resize(Ndat_,Ndat_);
  alongtrack.resize(Ndat_+1,Ndat_+1);
  crosstrack.resize(Ndat_+1,Ndat_+1);
  radar_geom_factor.resize(Ndat_,Ndat_);
  radar_geom_factordB.resize(Ndat_,Ndat_);
  number_of_solution.resize(Ndat_,Ndat_);
  
  oneway_beam_gain_mirror.resize(Ndat_,Ndat_);
  incidenceangle_mirror.resize(Ndat_,Ndat_);
 
  area_mirror.resize(Ndat_,Ndat_);
  radar_geom_factor_mirror.resize(Ndat_,Ndat_);
  number_of_solution_mirror.resize(Ndat_,Ndat_);

  good_bin.resize(Ndat_,Ndat_);
  ML_good_bin.resize(Ndat_,Ndat_);
  thermal_snr_SL.resize(Ndat_,Ndat_);
  amb_ratio_SL.resize(Ndat_,Ndat_);
  thermal_snr_ML.resize(Ndat_,Ndat_);
  amb_ratio_ML.resize(Ndat_,Ndat_);

  //initialize result variables
  range = Uvar(0,"km");
  doppler = Uvar(0,"Hz");

  oneway_beam_gain = 0.0;

  incidenceangle = Uvar(0,"rad");  
  area =Uvar(0,"km km");
  alongtrack = Uvar(0,"km");
  crosstrack = Uvar(0,"km");

  radar_geom_factor = 0.0;
  radar_geom_factordB = -100.0;
  number_of_solution = 0;

  oneway_beam_gain_mirror = 0.0;

  incidenceangle_mirror = Uvar(0,"rad"); 
  area_mirror = Uvar(0,"km km");

  radar_geom_factor_mirror= 0.0;
  number_of_solution_mirror = 0;

  good_bin = 0;
  ML_good_bin=0;
  thermal_snr_SL = 0;
  amb_ratio_SL = 0;
  thermal_snr_ML = 0;
  amb_ratio_ML = 0;
  grid_number_set_ = true;
  }

//------------------------------------
//{set,get}mirror ambiguity setting
//-----------------------------------

void Ambiguity::setMirrorAmbiguity(const bool& mirror_set)
  {
    mirror_geom_set_= mirror_set;
  }
bool Ambiguity::getMirrorAmbiguity()
  {
    return(mirror_geom_set_);
  }

//-------------
//getNumberofLooks()
//-------------- 
unsigned int Ambiguity::getNumberofLooks()
{
  if(!ML_usable_area_set_)
    {
      ErrorMessage("Ambiguity:: Multilook usable area calculation has not been performed").throwMe();
    }
  return(nlook_);
}

//--------------------------
//{set} backscatt model and process window
// Will not check whether these parameters have been set already
//to allow users to reset at any time
//---------------------------
void Ambiguity::setBackscattCoefficient(const double& mulheman_k1,
			       const double& mulheman_k2)
{
  muhleman_k1_ = mulheman_k1;
  muhleman_k2_ = mulheman_k2;
  backscatt_model_set_ = true; 
}

//-----------------------------------
//{set} process window
//------------------------------------
void Ambiguity:: setProcessWindow(const Uvar& center_range,
				  const Uvar& center_doppler,
				  const Uvar& prf,
				  const Uvar& pbw, 
				  const Uvar& pulse_gate)
				
{
  center_range_=center_range;
  center_doppler_=center_doppler;
  prf_ = prf;
  pri_ = 1/prf;
  pbw_ = pbw;
  pulse_gate_ = pulse_gate;
 
  process_window_set_ = true;
}

//--------------------------------
//get crosstrack min and crosstrack max inside usable area
//-------------------------------
void Ambiguity::getCrosstrackExtent(Uvar& cross_min,Uvar& cross_max,Uvar& crosstrack_length)
  {
    if (!ML_usable_area_set_) ErrorMessage("Ambiguity::getCrossTrackExtent: no multilook usable area set has been performed").throwMe();
    cross_min =  ML_crosstrack_extent_min_;
    cross_max =  ML_crosstrack_extent_max_;
    crosstrack_length=ML_crosstrack_length_;
  }
//----------------------
//clear()
//---------------------
void Ambiguity::clear()
  {
  time_set_ = false;
  target_set_=false;
  beam_set_=false;  
  ieb_set_ = false;
  config_set_ = false;
  mirror_geom_set_ = false;
  beam_number_set_ = false;
  beam_max_gain_set_ = false;
  grid_number_set_=false;
  radar_geom_factor_set_ = false;
  peak_radar_geom_factor_set_=false;
  usable_area_set_= false;
  ML_usable_area_set_=false;
  backscatt_model_set_ = false;
  process_window_set_ = false;
 }

//------------------------------
// return  range/doppler index for a given input
//-----------------------------
unsigned int Ambiguity::RangeIndex(const Uvar& range_value)
  {
  if (!grid_number_set_)
    {
      ErrorMessage e("Ambiguity::getRangeIndex: No grid parameter is set");
      e.throwMe();
    }
 
  if (range_value< range(0) || range_value > range(Ndat_-1))
    {
      ErrorMessage e("Ambiguity::getRangeIndex:Range index is out of range");
      e.throwMe();
    }
  Uvar step=range(Ndat_/2)-range(Ndat_/2-1);
  unsigned int range_index = (unsigned int) 
    ((range_value - range(0))/(step)).getInUnits("");  
  return (range_index);
  }

unsigned int Ambiguity::DopplerIndex(const Uvar& doppler_value)
  {
    if (!grid_number_set_)
      {
	ErrorMessage("Ambiguity::getDopplerIndex: No grid parameter is set").throwMe();
      }
    if (doppler_value< doppler(0) || doppler_value > doppler(Ndat_-1)) 
      {
	ErrorMessage("Dop index is out of range").throwMe();
      }
    Uvar step=doppler(Ndat_/2)-doppler(Ndat_/2-1);
    unsigned int doppler_index = (unsigned int)
      ((doppler_value-doppler(0))/(step)).getInUnits("");
   
    return(doppler_index);    
  }

//-----------------------------
//get range /doppler interval of the grid
//-----------------------------
Uvar Ambiguity::RangeInterval()
  {
  if (!grid_number_set_)
      {
	ErrorMessage("Ambiguity::getDopplerIndex: No grid parameter is set").throwMe();
      }
  return(range(1) - range(0));
  }

Uvar Ambiguity::DopplerInterval()
  {
  if (!grid_number_set_)
      {
	ErrorMessage("Ambiguity::getDopplerIndex: No grid parameter is set").throwMe();
      }
  return(doppler(1) - doppler(0));
  }


//---------------------------
//calAmbiguity()
//---------------------------
void Ambiguity::calAmbGeometry()
  {
  if(!time_set_)
    {
    ErrorMessage("Ambiguity::calAmbiguity: no time is set").throwMe();
    }
  if(!target_set_)
    {
    ErrorMessage("Ambiguity::calAmbiguity: no target is set").throwMe();
    }
  if(!beam_set_)
    {
    ErrorMessage("Ambiguity::calAmbiguity: no beam is set").throwMe();
    }
  if(!config_set_)
    {
    ErrorMessage("Ambiguity::calAmbiguity: no config file is read in").throwMe();
    }
  if(!track_frame_set_)
    {
    ErrorMessage("Ambiguity::calAmbiguity:no track frame is set").throwMe();
    }
       
  //----------------------
  //Need three direction vector to find a correct solution
  //-----------------------
  DirectionVector dir_pos = state_.position();
  DirectionVector dir_vel = state_.velocity();
  DirectionVector dir_cross = cross(dir_vel,dir_pos); 


  //-----------------------------------------------
  //determine roll direction using beam3 boresight direction
  //------------------------------------------------
  Frame fbeam3("CASSINI_RADAR_3","Cassini");
  DirectionVector boresight3("boresight",fbeam3,t_,0,0,1);
  double roll_direction=dot(boresight3,dir_cross);
  //if (roll_direction > 0) cout<<"right looking "<<endl;
  //else cout<<"left looking "<<endl;
  

  //---------------------------------
  //obtain boresight range and doppler
  //--------------------------------
  unsigned int beam_num= beam_.getBeamNumber();
  Frame fbeam("CASSINI_RADAR_" + toStr(beam_num),"Cassini");
  DirectionVector boresight("boresight",fbeam,t_,0,0,1);
 
  TargetGeom tg(t_);
  tg.setState(state_);
  tg.setLookDirection(boresight);
  tg.setTarget(target_name_,target_frame_);
  range_bore= tg.range();
  doppler_bore=tg.doppler(lambda_);
  radius_ = tg.radius();
  altitude=tg.altitude();
  boresight_incidence_angle = tg.incidenceAngle();
 
  //-----------------------
  //Grid variables declaration/initiallization
  //----------------------
  Uvec range_grid("range",Ngrid_);
  Uvec doppler_grid("frequency",Ngrid_);
  
  Imat no_of_solution_grid("no_solution",Ngrid_,Ngrid_);
  Dmat oneway_beam_gain_grid("beam_gain",Ngrid_,Ngrid_); 
  Umat incidenceangle_grid("thetai",Ngrid_,Ngrid_); 
  Array2D<PositionVector> surface_intercept_grid("surface",Ngrid_,Ngrid_);

  Imat no_of_solution_mirror_grid("no_solution",Ngrid_,Ngrid_);
  Dmat oneway_beam_gain_mirror_grid("beam_gain",Ngrid_,Ngrid_); 
  Umat  incidenceangle_mirror_grid("thetai",Ngrid_,Ngrid_);
  Array2D<PositionVector> surface_intercept_mirror_grid("surface",Ngrid_,Ngrid_);
  
  //initialize all the grid variables
  range_grid=Uvar(0,"km");
  doppler_grid=Uvar(0,"Hz");

  no_of_solution_grid=0;
  oneway_beam_gain_grid = 0.0;
  incidenceangle_grid = Uvar(0,"rad");
 
  no_of_solution_mirror_grid=0;
  oneway_beam_gain_mirror_grid = 0.0;
  incidenceangle_mirror_grid = Uvar(0,"rad");
 
  //--------------------  
  //Pixel variables resize/initiallization
  // (alongtrack,crosstrack) will keep its grid size
  // 
  //---------------------
  range.resize(Ndat_);
  doppler.resize(Ndat_);

  oneway_beam_gain.resize(Ndat_,Ndat_);
  incidenceangle.resize(Ndat_,Ndat_); 
  area.resize(Ndat_,Ndat_);
  alongtrack.resize(Ndat_+1,Ndat_+1);
  crosstrack.resize(Ndat_+1,Ndat_+1);
 
  number_of_solution.resize(Ndat_,Ndat_);

  oneway_beam_gain_mirror.resize(Ndat_,Ndat_);
  incidenceangle_mirror.resize(Ndat_,Ndat_);
  area_mirror.resize(Ndat_,Ndat_);
  number_of_solution_mirror.resize(Ndat_,Ndat_);
  
  range = Uvar(0,"km");
  doppler = Uvar(0,"Hz");

  oneway_beam_gain = 0.0;
  incidenceangle = Uvar(0,"rad");
  area =Uvar(0,"km km");
  alongtrack = Uvar(0,"km");
  crosstrack = Uvar(0,"km");
  number_of_solution = 0;

  oneway_beam_gain_mirror = 0.0;
  incidenceangle_mirror = Uvar(0,"rad");
  area_mirror =  Uvar(0,"km km");
  number_of_solution_mirror = 0;

  //--------------------------------
  //Assign values: max_delta_doppler: max_prf
  //               max_delta_range : min_prf
  // based on the preliminary amb study
  // it is not likely that we are going use prf higher than 6
  // and lower than 1.5 
  //---------------------------------
  Uvar delta_doppler = 2.0 * max_prf_;
  Uvar delta_range = 2.0 * speed_light/2.0/min_prf_;

  for (unsigned int i_grid = 0; i_grid < Ngrid_;++i_grid)
    {
    range_grid(i_grid) = -delta_range + range_bore;
    range_grid(i_grid) += 2.0 * delta_range * double(i_grid)/double(Ngrid_ -1);
    doppler_grid(i_grid) = -delta_doppler+doppler_bore;
    doppler_grid(i_grid) += 2.0 * delta_doppler*double(i_grid)/double(Ngrid_ -1 );
    //cout<<"assigning "<< range_grid(i_grid)<<" "<<doppler_grid(i_grid)<<endl;
    }

 

  //------------------------------------------------------------
  //This setup is needed to speed up the coordinate 
  //transformation from target_body_frame to beam_frame or track_frame
  // This needs to be done only once as long as time is not changed
  //------------------------------------------------------------------
  double a,b,c;
  DirectionVector target_to_beam_x,target_to_beam_y,target_to_beam_z;
  DirectionVector target_to_track_x,target_to_track_y,target_to_track_z;
 
  fbeam.axialVectors(target_frame_,
		     t_,
		     target_to_beam_x,
		     target_to_beam_y,
		     target_to_beam_z);
  ftrack_.axialVectors(target_frame_,
		       t_,
		       target_to_track_x,
		       target_to_track_y,
		       target_to_track_z);
  DirectionVector lookInBeamFrame("beam frame direction",fbeam,t_,0,0,1);
  DirectionVector dirToSurface("track frame direction",target_frame_,t_,0,0,1);
  double radius_in_km = radius_.getInUnits("km");
 
  
  //-------------- spacecraft state vector------
  DirectionVector pos = state_.position();
  DirectionVector vel = state_.velocity();
  DirectionVector ulook("DirectionVector from sc to surface",target_frame_,t_,0,0,1);
  DirectionVector ulook_mirror("mirrorLook from sc to surfact",target_frame_,t_,0,0,1);
  //---------------dummy variables to be used for solving RangeDoppler
  double p1,p2,p3, v1,v2,v3;
  double sol;
  int sign=1;

  p1 = pos[DirectionVector::X];
  p2 = pos[DirectionVector::Y];
  p3 = pos[DirectionVector::Z];
  
  v1 = vel[DirectionVector::X];
  v2 = vel[DirectionVector::Y];
  v3 = vel[DirectionVector::Z];
  
 

  double u1,u2,u3;
  double A,B,C,D,E,F;
  double range_in_radius;
  double position_in_radius = (state_.position().magnitude()/radius_).getInUnits("");
  double compute_roll;
  //-------------------------------------------------------------------
  
  
  for (unsigned int i_range = 0; i_range<Ngrid_; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ngrid_; ++j_dop)
    {	
    //------------- Instead of relying on calling targetgeom, 
    //let's solve everyting here
    range_in_radius = (range_grid(i_range)/radius_).getInUnits("");
    A =  (1 - position_in_radius*position_in_radius 
	  - range_in_radius*range_in_radius)
      /(2.0 * range_in_radius * position_in_radius);
    B =  (lambda_ * doppler_grid(j_dop)/(2 * state_.velocity().magnitude())).getInUnits("");

    C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
    D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
    E = (A - p1 *C)/p3;
    F = -(p2 + p1*D)/p3;
    
    sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);
    if (sol >= 0)
      {
      no_of_solution_grid(i_range,j_dop) = 1;   	 
      u2 = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
      u1 = C + D*u2;
      u3 = E + F*u2;

      //set look vector from sc to surface in target frame
      ulook[DirectionVector::X] = u1;
      ulook[DirectionVector::Y] = u2;
      ulook[DirectionVector::Z] = u3;
 
      //--------------------------------
      //set surface intercept point
      //--------------------------------
      surface_intercept_grid(i_range,j_dop) = state_.position();
      surface_intercept_grid(i_range,j_dop) += range_grid(i_range)*ulook;
     
      //-------------------------------------------
      //Find directly the along/cross distances
      //by using axialVectors between ftrack and ftitan
      // surface_intercept_grid is defined in target frame
      // in order to calculate along/cross, we need to transform 
      //the surface vector
      // in ftrack frame
      //-------------------------------------------   
  
      dirToSurface = surface_intercept_grid(i_range,j_dop);//dir vector in target frame
      compute_roll=dir_cross[DirectionVector::X]*dirToSurface[DirectionVector::X]
	+dir_cross[DirectionVector::Y]*dirToSurface[DirectionVector::Y]
	+dir_cross[DirectionVector::Z]*dirToSurface[DirectionVector::Z];
      
      if (compute_roll*roll_direction > 0)
	{
	sign = 1;
	}
      else
	{
	sign=-1;
	u2 = ( -(C*D + E*F)+sqrt( sol))/(D*D + F*F + 1.0);
	u1 = C + D*u2;
	u3 = E + F*u2;
	
	//set look vector from sc to surface in target frame
	ulook[DirectionVector::X] = u1;
	ulook[DirectionVector::Y] = u2;
	ulook[DirectionVector::Z] = u3;
	
	//--------------------------------
	//set surface intercept point
	//--------------------------------
	surface_intercept_grid(i_range,j_dop) = state_.position();
	surface_intercept_grid(i_range,j_dop) += range_grid(i_range)*ulook;
	dirToSurface = surface_intercept_grid(i_range,j_dop);//dir vector in target frame
	}
      

      a = dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::X]
	+dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::X]
	+dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::X];
      
      b= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Y]
	+dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Y]
	+dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Y];
      
      c= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Z]
	+dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Z]
	+dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Z];
      
      //---------------------------------------------
      //theta = -acos(c)+pi/2;//change to latitude
      //phi = atan2(b,a);
      //-------------------------------------------------
      crosstrack(i_range,j_dop).setValue(radius_in_km * (pi/2.0-acos(c)));
      alongtrack(i_range,j_dop).setValue(-radius_in_km *atan2(b,a));
    
      //if (radius_in_km * (pi/2.0-acos(c)) == 0 ||
      //  radius_in_km*atan2(b,a) == 0) 
      //ErrorMessage("surface solution, but along/cross are 0").throwMe();
      //---------------------------------------------------
      //direction vector look is defined in target frame
      // need to transform into beam frame
      //reuse double variables a,b, and c
      //ulook is defined in target frame
      // so need a rotation to beam frame
      //---------------------------------------------------
      a = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	+ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	+ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
      
      b = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	+ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	+ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
      
      c = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	+ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	+ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];

      lookInBeamFrame[DirectionVector::X] = a;
      lookInBeamFrame[DirectionVector::Y] = b;
      lookInBeamFrame[DirectionVector::Z] = c;
      oneway_beam_gain_grid(i_range,j_dop)=
	beam_.getGainOneWay(lookInBeamFrame);

      //----------------------
      //set incidence angle
      //------------------------
      a = dot(-ulook,dirToSurface);//reuse a
      incidenceangle_grid(i_range,j_dop)=  Uvar(acos(a),"rad");


      //debug
      /*
      TargetGeom tg1(t_); 
      tg1.setTrackFrame(ftrack_);
      tg1.setState(state_);
      tg1.setTarget(target_name_,target_frame_);
      tg1.setRangeDopplerInTargetFrame(range_grid(i_range),
				       doppler_grid(j_dop),
				       lambda_);
      Uvar tmp_a,tmp_c;
      tg1.interceptAlongCross(tmp_a,tmp_c);
      cout<<"TG along cross "<<tmp_a<<" "<<tmp_c<<endl;
      cout<<" "<<alongtrack(i_range,j_dop)<<" "<<crosstrack(i_range,j_dop)<<endl;
      cout<<" "<<tg1.mirrorsurfaceIntercept()<<endl;
      */

      if(mirror_geom_set_)
	{
	no_of_solution_mirror_grid(i_range,j_dop) = 1;
	u2 = ( -(C*D + E*F) +double(sign)* sqrt( sol))
	  /(D*D + F*F + 1.0);
	u1 = C + D*u2;
	u3 = E + F*u2;
	ulook_mirror[DirectionVector::X] = u1;
	ulook_mirror[DirectionVector::Y] = u2;
	ulook_mirror[DirectionVector::Z] = u3;
	//-----------------------------
	//mirror surface intercept
	//----------------------------
	surface_intercept_mirror_grid(i_range,j_dop) = state_.position();
	surface_intercept_mirror_grid(i_range,j_dop)+= range_grid(i_range) * ulook_mirror;
	//cout<<"mirror "<<	surface_intercept_mirror_grid(i_range,j_dop)<<endl;
     
	//--------------------------
	//incidence angle
	//dirToSurface and ulook_mirror are  defined in target frame
	//-------------------------
	dirToSurface = surface_intercept_mirror_grid(i_range,j_dop);//reuse to_surface vector
	a = dot(-ulook_mirror,dirToSurface);
	incidenceangle_mirror_grid(i_range,j_dop) = Uvar(acos(a),"rad");
	
	//---------------------------
	//Again manual transformation 
	//-----------------------
	a = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	  +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	  +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	
	b = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	  +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	  +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	
	c = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	  +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	  +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];

	lookInBeamFrame[DirectionVector::X] = a;
	lookInBeamFrame[DirectionVector::Y] = b;
	lookInBeamFrame[DirectionVector::Z] = c;
	oneway_beam_gain_mirror_grid(i_range,j_dop)
	  =beam_.getGainOneWay(lookInBeamFrame);  
	}//only if mirror geom set is true
      }
    else
      {
	no_of_solution_grid(i_range,j_dop) = 0;//no surface intercept point
	no_of_solution_mirror_grid(i_range,j_dop) = 0;//no surface intercept point
      }
    }//Loop over i_range and j_dop
  
  //print current time
  //cout<<Time::getCurrentLocalTime()<<endl;

  //-----------------------------
  //Old Method that fully uses targetgeom method
  //-----------------------------
  //TargetGeom tg1(t_); 
  //tg1.setTrackFrame(ftrack_);
  //for (unsigned int i_range = 0; i_range<Ngrid_; ++i_range)
  //for (unsigned int j_dop =0; j_dop<Ngrid_; ++j_dop)
  //{	
  //tg1.setState(state_);
  //tg1.setTarget(target_name_,target_frame_);
  //tg1.setRangeDopplerInTargetFrame(range_grid(i_range),
  //			     doppler_grid(j_dop),
  //			     lambda_,
  //			     roll_direction_);
  //if (tg1.foundSurfaceIntercept()==false)
  //  {
  //   no_of_solution_grid(i_range,j_dop) = 0;//no surface intercept point
  //   no_of_solution_mirror_grid(i_range,j_dop) = 0;//no surface intercept point
  //  }
  //  else
  //{
  //no_of_solution_grid(i_range,j_dop) = 1;   
  //incidenceangle_grid(i_range,j_dop)= tg1.incidenceAngle(); 	
  //surface_intercept_grid(i_range,j_dop)=tg1.surfaceIntercept();
  //	tg1.interceptAlongCross(alongtrack(i_range,j_dop),
  //	crosstrack(i_range,j_dop));
  //
  //oneway_beam_gain_grid(i_range,j_dop)=
  //     	  beam_.getGainOneWay(tg1.lookDirection());
  //--------------------------
  //take care of mirror look only if mirror_geom_set_ == true
  //---------------------
  //if(mirror_geom_set_)
  //  {
  //  no_of_solution_mirror_grid(i_range,j_dop) = 1;
  //  incidenceangle_mirror_grid(i_range,j_dop)= tg1.mirrorincidenceAngle(); 
  //  surface_intercept_mirror_grid(i_range,j_dop)
  //    =tg1.mirrorsurfaceIntercept();
  //
  //  oneway_beam_gain_mirror_grid(i_range,j_dop)
  //    =beam_.getGainOneWay(tg1.mirrorlookDirection());  
  //	  }//only if mirror geom set is true
  //} 
  //tg1.reset(t_);
  //}//Loop over i_range and j_dop
  //
  //print current time
  //cout<<Time::getCurrentLocalTime()<<endl;



  //--------------------
  //pixel range and doppler
  //-----------------
  for (unsigned int i_range = 0; i_range<Ngrid_-1; ++i_range)
    {
      range(i_range)= (range_grid(i_range) + range_grid(i_range + 1))/2.0;
    }
  for (unsigned int j_dop =0; j_dop<Ngrid_-1; ++j_dop)
    {
      doppler(j_dop)=(doppler_grid(j_dop)+doppler_grid(j_dop+1))/2.0;
    }

  //------------------------------------------------
  //Finalize data set into Ndat_ X Ndat_ matrix form
  //construct Ndat_ x Ndat_ data set
  //-------------------------------------------------

  //------------------
  //For area calculation, i will use the fact that three position vectors
  // are defined in the same frame
  //------------------
  //use double a b c which are declared earlier as double dummy variables
  DirectionVector dir_s1,dir_s2,dir_s3;
  double side_a,side_b,side_c,cos_c;
 

  for (unsigned int i_range = 0; i_range<Ngrid_-1; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ngrid_-1; ++j_dop)
    {
    if (no_of_solution_grid(i_range,j_dop)
	*no_of_solution_grid(i_range+1,j_dop)
	*no_of_solution_grid(i_range,j_dop+1)
	*no_of_solution_grid(i_range+1,j_dop+1) == 1)
      {
      number_of_solution(i_range,j_dop) = 1;
      //------------------------------------
      //find Area element covered by (range,dop)
      //divide an area enclosed by two range and dop values
      // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
      //into two triangles and calclate their areas and add them
      //together,Well.  Instead, just calculate the area of one
      //triangle and muliply it by 2 to save time
      //-------------------------------------  
      //area(i_range,j_dop) = 2.0*
      //surfacearea(surface_intercept_grid(i_range,j_dop),
      //	    surface_intercept_grid(i_range+1,j_dop),
      //	    surface_intercept_grid(i_range,j_dop+1),
      //	    radius_);
      //--------------- Do not use ---- time consuming ------

      dir_s1 = surface_intercept_grid(i_range,j_dop);
      dir_s2 = surface_intercept_grid(i_range+1,j_dop);
      dir_s3 = surface_intercept_grid(i_range,j_dop+1);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_in_km* acos(a);
      side_b = radius_in_km* acos(b);
      side_c = radius_in_km* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //We calculate only half of the rectangle defined by
      // four corner values of position vectors and then multiply it by 2
      // to have the total area.
      area(i_range,j_dop).setValue(2.0*0.5*sin(acos(cos_c))*side_a*side_b);

      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain(i_range,j_dop) = (
		oneway_beam_gain_grid(i_range,j_dop)
		+oneway_beam_gain_grid(i_range+1,j_dop)
		+oneway_beam_gain_grid(i_range,j_dop+1)
		+oneway_beam_gain_grid(i_range+1,j_dop+1)
		)/4.0;
      
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle(i_range,j_dop) = (
		  incidenceangle_grid(i_range,j_dop) 
		  +incidenceangle_grid(i_range+1,j_dop)
		  +incidenceangle_grid(i_range,j_dop+1)
		  +incidenceangle_grid(i_range+1,j_dop+1)
		  )/4.0;
      }
    else
      {
	number_of_solution(i_range,j_dop) = 0;
      }//no solution, return zero
    }//Loop over range/dop values

  //print current time
  //cout<<Time::getCurrentLocalTime()<<endl;

  if (!mirror_geom_set_) {return;}
  //if no mirror option is chosen, return here
  //----------------------------------------------
  //take care of mirror amb  
  // if mirror_geom_set_ is false, return
  //----------------------------------------------
 
  for (unsigned int i_range = 0; i_range<Ngrid_-1; i_range++)
  for (unsigned int j_dop =0; j_dop<Ngrid_-1; j_dop++)
    {
    if (no_of_solution_mirror_grid(i_range,j_dop)
	*no_of_solution_mirror_grid(i_range+1,j_dop)
	*no_of_solution_mirror_grid(i_range,j_dop+1)
	*no_of_solution_mirror_grid(i_range+1,j_dop+1) == 1)
      {
      number_of_solution_mirror(i_range,j_dop) = 1;
      //------------------------------------
      //find Area element covered by (range,dop)
      //divide an area enclosed by two range and dop values
      // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
      //into two triangles and calclate their areas and add them
      //together.
      //-------------------------------------  
      //      area_mirror(i_range,j_dop) = 2.0*
      //surfacearea(surface_intercept_mirror_grid(i_range,j_dop),
      //	    surface_intercept_mirror_grid(i_range+1,j_dop),
      //	    surface_intercept_mirror_grid(i_range,j_dop+1),
      //	    radius_);
      //instead of calculating the areas of two triangles,
      // just do it once and multiply by 2.0      
      //area_mirror(i_range,j_dop) += 
      //surfacearea(surface_intercept_mirror_grid(i_range+1,j_dop+1),
      //	    surface_intercept_mirror_grid(i_range+1,j_dop),
      //	    surface_intercept_mirror_grid(i_range,j_dop+1),
      //	    radius_);
      //----------------- Do not use ----------------------
      dir_s1 = surface_intercept_mirror_grid(i_range,j_dop);
      dir_s2 = surface_intercept_mirror_grid(i_range+1,j_dop);
      dir_s3 = surface_intercept_mirror_grid(i_range,j_dop+1);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_in_km* acos(a);
      side_b = radius_in_km* acos(b);
      side_c = radius_in_km* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //multiply 2 for the other triangle
      //remember: we are calculating the area of a rectangle defined by
      // four corners
      // save output to area_mirror!
      area_mirror(i_range,j_dop).setValue(
			     2.0*0.5*sin(acos(cos_c))*side_a*side_b);

      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain_mirror(i_range,j_dop) = (
	     oneway_beam_gain_mirror_grid(i_range,j_dop)
	     +oneway_beam_gain_mirror_grid(i_range+1,j_dop)
	     +oneway_beam_gain_mirror_grid(i_range,j_dop+1)
	     +oneway_beam_gain_mirror_grid(i_range+1,j_dop+1)
	     )/4.0;
      
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle_mirror(i_range,j_dop) = (
	   incidenceangle_mirror_grid(i_range,j_dop) 
	   + incidenceangle_mirror_grid(i_range+1,j_dop)
	   +incidenceangle_mirror_grid(i_range,j_dop+1)
	   +incidenceangle_mirror_grid(i_range+1,j_dop+1)
	   )/4.0;    
      }
    else
      {
      number_of_solution_mirror(i_range,j_dop) = 0;
      }//no solution, return zero
    }//Loop over range/dop values
  }

//----------------------------------
//Calculate radar geometry
//-----------------------------------
void Ambiguity::calRadarGeometryFactor()
  {
  if (radar_geom_factor_set_)
    {
    ErrorMessage e("Ambiguity::calRadarGeometryFactor: Calculation is already done");
    e.throwMe();
    }
  if (!grid_number_set_)
    {
      ErrorMessage e("Ambiguity::calRadarGeometryFactor: No grid parameter is set");
      e.throwMe();
    }
  if(!backscatt_model_set_)
    {
      ErrorMessage("Ambiguity::calRadarGeometryFactor: Backscattering model coefficients are not set").throwMe();
    }
  for (unsigned int i = 0 ; i < Ndat_;++i)
  for(unsigned int j = 0; j < Ndat_;++j)
    {
      //cout<<"i j solution "<<i<<j<<number_of_solution(i,j)<<endl;
      if (number_of_solution(i,j) == 0)
	{
	  radar_geom_factor(i,j)=0;
	}
      else if (number_of_solution(i,j)== 1)
	{
	  Uvar X = oneway_beam_gain(i,j)
	    *oneway_beam_gain(i,j)*lambda_ *lambda_
	    * muhleman_backscatter(muhleman_k1_,
				   muhleman_k2_,
				   incidenceangle(i,j))
	    *area(i,j)/pow(range(i),4);
	  radar_geom_factor(i,j)=X.getInUnits("");
	}
      else
	{
	  ErrorMessage("Ambiguity::calRadarGeometryFactor: number of solution set is contaminated").throwMe();
	}
    }
  radar_geom_factor_set_ = true;
  
  //take care of mirror amb only if mirror_geom_set is true
  if (!mirror_geom_set_) {return;}
  for (unsigned int i = 0 ; i < Ndat_;++i)
  for(unsigned int j = 0; j < Ndat_;++j)
    {
      if (number_of_solution_mirror(i,j) == 0)
      	{
	  radar_geom_factor_mirror(i,j) =0.0;
	}
      else if (number_of_solution_mirror(i,j) ==1)
	{
	  Uvar X_mirror= oneway_beam_gain_mirror(i,j) 
	    * oneway_beam_gain_mirror(i,j) * lambda_ * lambda_
	    *muhleman_backscatter(muhleman_k1_,
				  muhleman_k2_,
				  incidenceangle_mirror(i,j))
	    * area_mirror(i,j)/ pow(range(i),4);
	  radar_geom_factor_mirror(i,j) =X_mirror.getInUnits("");
	}
      else
	{
	  ErrorMessage("Ambiguity::calRadarGeomFactor: number_of_solution_mirror is contaminated").throwMe();
	}
    }
  }

//---------------------------------------------------------
//Convert radar geom factor in dB with respect to the max value
// radar geom factor
//---------------------------------------------------------
 void Ambiguity::calRadarGeometryFactordB(const double& out_of_bounds)
  {
    double max = 0.0;
    for (unsigned int i = 0; i < Ndat_;++i)
    for (unsigned int j = 0; j < Ndat_;++j)
      {
	if(max < radar_geom_factor(i,j)) max = radar_geom_factor(i,j);
      }
   
    for (unsigned int i = 0; i < Ndat_;++i)
    for (unsigned int j = 0; j < Ndat_;++j)
      {
	if (number_of_solution(i,j) ==0)
	  {
	    radar_geom_factordB(i,j) = out_of_bounds;
	  }
	else if (number_of_solution(i,j) ==1)
	  {
	    radar_geom_factordB(i,j) = 10.0 * 
	      log(radar_geom_factor(i,j)/max)/log(10.0);
	    if (radar_geom_factordB(i,j) < out_of_bounds)
	      radar_geom_factordB(i,j) = out_of_bounds;
	  }
      }
  } 

//-------------------------------------------
//Calculate usable area for a given prf and pulse gate
//-------------------------------------------
Uvar Ambiguity::calUsableArea(const Uvar& X0, 
			      const Uvar& Pn,     
			      const Uvar& x_res,
			      const Uvar& rg_res,
			      const double& signal_to_ambdB, 
			      const double& noise_equivalent_sigma0dB, 
			      const double& min_gaindB)
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage e("Ambiguity::calUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
    if(!process_window_set_)
      {
	ErrorMessage("Ambiguity::calUsableArea: process window is not set").throwMe();
      }

    Uvar usable_area=Uvar(0.0,"km km");
    Uvar range_interval = range(1)- range(0);
    Uvar dop_interval = doppler(1) - doppler(0);
   
    //------------------------------------
    //collect data inside process window
    // to do that, let's find index values corresponding 
    // range_start, range_end, dop_start, dop_end
    //-----------------------------------
    Uvar range_start = center_range_-pulse_gate_/2.0;
    Uvar range_end = center_range_+ pulse_gate_/2.0;
    Uvar dop_start = center_doppler_ - pbw_/2.0;
    Uvar dop_end = center_doppler_ + pbw_/2.0;
    
    unsigned int i_range_start = RangeIndex(range_start);
    unsigned int i_range_end = RangeIndex(range_end);     
    unsigned int j_dop_start = DopplerIndex(dop_start);    
    unsigned int j_dop_end = DopplerIndex(dop_end);

     
    //---------------
    //examine data values falling inside 
    //range_start/range_end, dop_start/dop_end
    //-----------------
   
    good_bin = 0;//reset good bin indicator
    thermal_snr_SL = 0.0;
    amb_ratio_SL = 0.0;
    for (unsigned int i_range = i_range_start; i_range <i_range_end;++i_range)
    for (unsigned int j_dop =j_dop_start; j_dop < j_dop_end ;++j_dop)
      {
	//calculate amb ratio from four nearby amb patches
	//Left: -prf
	//Up: pri*c/2
	//Down: pri*c/2
	//Right: +prf
	unsigned int range_index_displacement = 
	  (unsigned int)(pri_*speed_light/2.0/range_interval).getInUnits("");
	unsigned int dop_index_displacement = 
	  (unsigned int)(prf_/dop_interval).getInUnits("");
	
	unsigned int range_Up = i_range + range_index_displacement;
	unsigned int range_Down=i_range - range_index_displacement;
	unsigned int dop_Left = j_dop + dop_index_displacement;
	unsigned int dop_Right = j_dop - dop_index_displacement;
	
	if (range_Up >= Ndat_ 
	    || range_Down>=Ndat_ 
	    || dop_Left >=Ndat_ 
	    || dop_Right >=Ndat_)  
	  {
	    ErrorMessage("Single look ambiguity:index is out of range").throwMe();
	  }
	
	double amb_radar_geom = radar_geom_factor(range_Up,j_dop)
	  +radar_geom_factor(range_Down,j_dop)
	  +radar_geom_factor(i_range,dop_Left)
	  +radar_geom_factor(i_range,dop_Right);
	if (mirror_geom_set_)
	  {//only if mirror_geom_set_ is true
	  amb_radar_geom = amb_radar_geom
	     +radar_geom_factor_mirror(range_Up,j_dop)
	     +radar_geom_factor_mirror(range_Down,j_dop)
	     +radar_geom_factor_mirror(i_range,dop_Left)
	     +radar_geom_factor_mirror(i_range,dop_Right)
	     +radar_geom_factor_mirror(i_range,j_dop);	
	  }
	//------------------------------------------
	//calculate amb ratio and noise equivalent sigma0
	//-----------------------------------------
	double ambratiodB = 10.0 * log(radar_geom_factor(i_range,j_dop)
				       /amb_radar_geom)/log(10.0);
	
	double area_ratio= (x_res * rg_res/area(i_range,j_dop)).getInUnits("");
	double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
	double signal_in_Watt = 
	  (X0
	   *radar_geom_factor(i_range,j_dop)
	   *area_ratio
	   /muhleman_backscatter(muhleman_k1_
				 ,muhleman_k2_
				 ,incidenceangle(i_range,j_dop))
	   ).getInUnits("kg km km/(s s s)");
	double nesigma0dB=10.0
	  * log(noise_in_Watt/signal_in_Watt)/log(10.0);
	double gaindB = 10.0 
	  * log(oneway_beam_gain(i_range,j_dop)/beam_max_gain_)/log(10.0);

	//------
	//Store thermal snr and amb ratio 
	//------
	thermal_snr_SL(i_range,j_dop) = -nesigma0dB;
	amb_ratio_SL(i_range,j_dop) = ambratiodB;
	//-------------
	//Pick usable area satisfying all the following conditions
	//-------------
	if (ambratiodB>= signal_to_ambdB &&
	    nesigma0dB<= noise_equivalent_sigma0dB &&
	    gaindB >= min_gaindB)
	  { 
	    usable_area+=area(i_range,j_dop);
	    good_bin(i_range,j_dop) = 1;
	  }
      }//loop over i_range, j_dop
    usable_area_set_ = true;
    return(usable_area);
  }

//---------------------------------------
//Calculate usable area for multilookings
//--------------------------------------
Uvar Ambiguity::calMultilookUsableArea(const Uvar& X0, 
				       const Uvar& Pn, 
				       const Uvar& x_res,
				       const Uvar& rg_res,
				       const Uvar& bpd,
				       const FloatVector& velocity,
				       const unsigned int& N_p,
				       const double& signal_to_ambdB, 	
				       const double& noise_equivalent_sigma0dB, 
				       const double& min_gaindB)
  {
  if (!radar_geom_factor_set_)
      {
	ErrorMessage e("Ambiguity::calMultilookUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
  if (!process_window_set_)
    {
      ErrorMessage("Ambiguity::calMultilookUsabeArea: process window is not set").throwMe();
    }
  Uvar usable_area=Uvar(0.0,"km km");
  Uvar range_interval = range(1)- range(0);
  Uvar dop_interval = doppler(1) - doppler(0);
  

  //------------------------------------
  //collect data inside process window
  // to do that, let's find index values corresponding 
  // range_start, range_end, dop_start, dop_end
  //-----------------------------------
  Uvar range_start = center_range_ - pulse_gate_/2.0;
  Uvar range_end = center_range_  + pulse_gate_/2.0;
  Uvar dop_start = center_doppler_ - pbw_/2.0;
  Uvar dop_end = center_doppler_ + pbw_/2.0;
    
  unsigned int i_range_start = RangeIndex(range_start);
  unsigned int i_range_end = RangeIndex(range_end);
   
  unsigned int j_dop_start =DopplerIndex(dop_start);
  unsigned int j_dop_end = DopplerIndex(dop_end);
  unsigned int j_dop_center = DopplerIndex(center_doppler_);
 
  //--------------------------------
  //Now calculate number of looks for data set
  //encloded by (i_range_start_ML,j_dop_start_ML)
  //and (i_range_end_ML, j_dop_end_ML)
  //--------------------------------
  double cos_theta = 
    (doppler_bore * lambda_ /(2.0 * velocity.magnitude())).getInUnits("");
  if (cos_theta > 1.0)
    {
      ErrorMessage("Ambiguity::cos theta value is larger than 1.0").throwMe();
    }

  Uvar ML_freq_bp = 2.0 * velocity.magnitude() 
                    * velocity.magnitude()*bpd/(lambda_*center_range_);
  ML_freq_bp *= (1.0-cos_theta * cos_theta);

  Uvar ML_freq_on = 1/(pri_ * double(N_p));
  //frequency resolution=1/(integration time(=pri * N_p))
 
  //cout<<"frequency shift "<<ML_freq_bp<<endl;
  //-------------- SH's simplied model produces the same result----------
  //double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  //Uvar ML_freq_bp = 
  //2.0 * velocity.magnitude()* Vst_bore * bpd*sin_theta/(lambda*range_bore);
  //note: Vst_bore = V * (1-cos(theta)^2)^0.5 --> same result!!!!!
  //for detailed information, see c:/ygim/radar_work/iso_doppler_profile
  //-------------------------------------------------------------------------

  //unsigned int ML_i_fon = 
  //(unsigned int) (ML_freq_on / dop_interval).getInUnits("");
 
  unsigned int ML_i_fbp = 
    (unsigned int) (ML_freq_bp / dop_interval).getInUnits("");

  if (ML_i_fbp ==0)
    {
      cout<<"you need longer burst period "<<endl;
      ML_i_fbp=1;
      //ErrorMessage("amb_param_scan:Need longer burst period to produce a shift in dop cell").throwMe();
      //infinite number of looks, unrealistic solution
    }

  vector<Uvar> doppler_container;
  doppler_container.clear();
  for (unsigned int i_range = i_range_start; i_range <i_range_end;++i_range)   
  for (unsigned int j_dop = j_dop_start; j_dop <j_dop_end;++j_dop)
    {
      if(10.0*log(oneway_beam_gain(i_range,j_dop)/beam_max_gain_)/log(10.0)>= min_gaindB)
	{
	  doppler_container.push_back(doppler(j_dop));
	}
    }
  
  if (doppler_container.size() == 0)
    {
      nlook_ = 1;
    }
  else if (doppler_container.size()==1)
    {
      nlook_ = 1;
    }
  else
    {
      sort(doppler_container.begin(),doppler_container.end());
      Uvar min_dop, max_dop;
      min_dop = doppler_container.front();
      max_dop = doppler_container.back();
      nlook_ = (unsigned int) round_double(((max_dop - min_dop)/ ML_freq_bp).getInUnits(""));
    }

  cout<<"nlook "<< nlook_<<endl;
  //------
  //setting up mask
  //----
  Array1D<unsigned int> mask("mask",Ndat_);
  for (unsigned int j_dop = 0; j_dop < Ndat_;++j_dop)
    {
      mask(j_dop) = j_dop % ML_i_fbp;	
    }  
  



  //---------------------------------
  //setup crosstrack variables to measure
  // ML usable area crosstrack extent(min and max crosstrack)
  // 
  //--------------------------------
  vector <Uvar> crosstrack_extent;
  crosstrack_extent.clear();
  Uvar crosstrack_length("crosstrack length ",0,"km");

  //----------------------------------
  //routine for calculating multilook
  //----------------------------------

  ML_good_bin = 0;//reset good bin indicator
  thermal_snr_ML = 0.0;
  amb_ratio_ML = 0.0;

  for (unsigned int i_range = i_range_start; i_range <i_range_end;++i_range)
    {     
      //-----------------------------
      //Here, I would like to calculate min signal_to_amb ration
      //after multilook
      //-----------------------------
      Uvec ML_signal_bin("ML signal bin", ML_i_fbp);
      Uvec ML_amb_bin("ML amb bin",ML_i_fbp); 
      Uvec ML_every_signal_bin("ML signal bin", ML_i_fbp);
      Uvec ML_every_amb_bin("ML amb bin",ML_i_fbp); 
      ML_signal_bin = 0.0;
      ML_amb_bin = 0.0;
      ML_every_signal_bin = 0.0;
      ML_every_amb_bin = 0.0;
      
   
      for (unsigned int j_dop =j_dop_start; j_dop < j_dop_end ;++j_dop)
	{
	  //calculate amb ratio from four nearby amb patches
	  //Left: -prf
	  //Up: pri*c/2
	  //Down: pri*c/2
	  //Right: +prf
	  unsigned int range_index_displacement = 
	    (unsigned int)(pri_*speed_light/2.0/range_interval).getInUnits("");
	  unsigned int dop_index_displacement = 
	    (unsigned int)(prf_/dop_interval).getInUnits("");
	  
	  unsigned int range_Up = i_range + range_index_displacement;
	  unsigned int range_Down=i_range - range_index_displacement;
	  unsigned int dop_Left = j_dop + dop_index_displacement;
	  unsigned int dop_Right = j_dop - dop_index_displacement;

	  double amb,signal;
	  amb=0.0;
	  signal=0.0;

	  if (range_Up >= Ndat_ 
	      || range_Down>=Ndat_ 
	      || dop_Left >=Ndat_ 
	      || dop_Right >=Ndat_)  
	    {
	      ErrorMessage("Single look ambiguity:index is out of range").throwMe();
	    }

	  signal = radar_geom_factor(i_range,j_dop);
	  amb = radar_geom_factor(range_Up,j_dop)
	    +radar_geom_factor(range_Down,j_dop)
	    +radar_geom_factor(i_range,dop_Left)
	    +radar_geom_factor(i_range,dop_Right)
	    +radar_geom_factor_mirror(range_Up,j_dop);
	  if (mirror_geom_set_)
	    {
	      amb += radar_geom_factor_mirror(range_Down,j_dop)
		+radar_geom_factor_mirror(i_range,dop_Left)
		+radar_geom_factor_mirror(i_range,dop_Right)
		+radar_geom_factor_mirror(i_range,j_dop);
	    }//only if mirror geom set is true



	  //for every bin whether there are inside or outside
	  // min beam gain
	  ML_every_amb_bin(mask(j_dop)) += amb;
	  ML_every_signal_bin(mask(j_dop)) += signal;
	  //special treatment for bins inside process window
	  if(10.0*log(oneway_beam_gain(i_range,j_dop)/beam_max_gain_)/log(10.0)>= min_gaindB)
	    {	  
	      ML_signal_bin(mask(j_dop)) +=signal;
	      ML_amb_bin(mask(j_dop))+=amb;
	    }//multilook only when oneway gain > min gain
	}//loop over j_dop

      //---------------------------------------
      //find multilook ambiguity ratio
      //---------------------------------------
      vector<Uvar> amb_inside_window;
      amb_inside_window.clear();
      for (unsigned int i_fbp = 0; i_fbp < ML_i_fbp;++i_fbp)
	{
	  if(ML_signal_bin(i_fbp)!= 0.0)
	    {
	    amb_inside_window.push_back(10.0 * log(ML_signal_bin(i_fbp)/ML_amb_bin(i_fbp))/log(10.0));
	    }
	}
     
      //------------------
      //find lowest ambiguity ratio while occurence_count > 1
      //------------------
      Uvar  min_amb_ratio_inside_min_gain ;      
      if(amb_inside_window.size()==0) min_amb_ratio_inside_min_gain = 0;
      else
	{
	  sort(amb_inside_window.begin(),amb_inside_window.end());
	  min_amb_ratio_inside_min_gain=amb_inside_window.front();
	}
      //------------------------
      //Assign amb values
      //----------------------------
      for (unsigned int j_dop =j_dop_start; j_dop < j_dop_end;++j_dop)
	{
	  double area_ratio= (x_res * rg_res/area(i_range,j_dop)).getInUnits("");
	  double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
	  double signal_in_Watt = 
	    (X0
	     *radar_geom_factor(i_range,j_dop)
	     *area_ratio
	     /muhleman_backscatter(muhleman_k1_,
				   muhleman_k2_,
				   incidenceangle(i_range,j_dop))
	     ).getInUnits("kg km km/(s s s)");
	  double nesigma0dB=10.0
	    * log(noise_in_Watt/signal_in_Watt/sqrt(double(nlook_)))/log(10.0);

	  //------
	  //Store thermal snr and amb ratio 
	  //------
	  thermal_snr_ML(i_range,j_dop) = -nesigma0dB;

	  //-----------------------------------
	  //Set amb ratio
	  //----------------------------------
 	  if(10.0*log(oneway_beam_gain(i_range,j_dop)/beam_max_gain_)/log(10.0)>= min_gaindB)
	    {
	      amb_ratio_ML(i_range,j_dop)= min_amb_ratio_inside_min_gain;
	    }
	  else
	    {
	      amb_ratio_ML(i_range,j_dop)= 10.0 
		* log(ML_every_signal_bin(mask(j_dop))/ML_every_amb_bin(mask(j_dop))/log(10.0));
	    }
	  //-------------
	  //Pick usable area satisfying all the following conditions
	  //-------------
	  if (amb_ratio_ML(i_range,j_dop).getInUnits("")>=signal_to_ambdB &&
	      thermal_snr_ML(i_range,j_dop) >=- noise_equivalent_sigma0dB &&
	     (10.0*log(oneway_beam_gain(i_range,j_dop)/beam_max_gain_)/log(10.0)) >= min_gaindB)
	    {
	      usable_area+=area(i_range,j_dop);
	      ML_good_bin(i_range,j_dop) = 1;
	      crosstrack_extent.push_back(crosstrack(i_range,j_dop));
	      crosstrack_extent.push_back(crosstrack(i_range+1,j_dop));
	      if(j_dop==j_dop_center)
		{
		  Uvar a=crosstrack(i_range+1,j_dop)-crosstrack(i_range,j_dop);
		  crosstrack_length +=Uvar(fabs(a.getInUnits("km")),"km");
		  //always take positve value
		}
	    }
	}//loop over j_dop
    }//loop over i_range

  //------------------------------
  //Let's extract crosstrack min and max
  //-------------------------------
  ML_crosstrack_extent_min_ = Uvar(0,"km");
  ML_crosstrack_extent_max_ = Uvar(0,"km");
  if (crosstrack_extent.size() > 0)
    {//only if ML usable area is not zero
      sort(crosstrack_extent.begin(),crosstrack_extent.end());//sorting
      ML_crosstrack_extent_min_ = crosstrack_extent.front();//first element
      ML_crosstrack_extent_max_ = crosstrack_extent.back();//last element
    }
  ML_crosstrack_length_=crosstrack_length;
  ML_usable_area_set_ = true;
  return(usable_area);
  }

//Static const member initiallization for AmbiguityFNN
// Ndat_ corresponds to pixel size
// Ngrid_ corresponds to grid size
// Therefore, pixel size should be 1 less than grid size

const unsigned int AmbiguityFNN::Ngrid_ = 41;
const unsigned int AmbiguityFNN::Ndat_ = 40;
const unsigned int AmbiguityFNN::Npatch_ = 5;
const unsigned int AmbiguityFNN::Ncenter_patch_ = 2;

AmbiguityFNN::AmbiguityFNN(const Time& t)
  :range("range"),
   doppler("dop"),
   oneway_beam_gaindB("oneway_gain"),
   incidenceangle("incidence_angle"),
   lat("latitude"),
   lon("longitude"),
   crosstrack("crosstrack"),
   alongtrack("alongtrack"),
   area("area"),
   radar_geom_factor("radar_geom_factor"),
   echo_power("echo power"),
   radar_geom_factor_mirror("radar_geom_factor mirror"),
   good_bin("good_bin"),
   ML_good_bin("ML_good_bin"),
   thermal_snr_SL("thermal_snr_SL"),
   amb_ratio_SL("amb_ratio_SL"),
   thermal_snr_ML("thermal_snr_ML"),
   amb_ratio_ML("amb_ratio_ML"),
   time_set_(true),
   target_set_(false),
   state_set_(false),
   beam_set_(false), 
   config_set_(false),
   mirror_geom_set_(false),
   beam_number_set_(false),
   beam_max_gain_set_(false),
   process_window_set_(false),
   radar_geom_factor_set_(false),
   peak_radar_geom_factor_set_(false),
   usable_area_set_(false),
   ML_usable_area_set_(false),
   track_frame_set_(false),
   t_(t),
   nlook_(1)
  {
  range.resize(Ndat_);
  doppler.resize(Ndat_);
  oneway_beam_gaindB.resize(Ndat_,Ndat_);
  incidenceangle.resize(Ndat_,Ndat_);
  lat.resize(Ngrid_,Ngrid_);
  lon.resize(Ngrid_,Ngrid_);
  crosstrack.resize(Ngrid_,Ngrid_);
  alongtrack.resize(Ngrid_,Ngrid_);
  area.resize(Ndat_,Ndat_);
  radar_geom_factor.resize(Ndat_,Ndat_,Npatch_);
  echo_power.resize(Ndat_,Ndat_);
  radar_geom_factor_mirror.resize(Ndat_,Ndat_,Npatch_);
  good_bin.resize(Ndat_,Ndat_);
  ML_good_bin.resize(Ndat_,Ndat_); 
  thermal_snr_SL.resize(Ndat_,Ndat_);
  amb_ratio_SL.resize(Ndat_,Ndat_);
  thermal_snr_ML.resize(Ndat_,Ndat_);
  amb_ratio_ML.resize(Ndat_,Ndat_);

  //initialize result variables
  range = Uvar(0,"km");
  doppler = Uvar(0,"Hz");
  oneway_beam_gaindB = -100.0;
  incidenceangle = Uvar(0,"rad");
  lat = Uvar(0,"rad");
  lon = Uvar(0,"rad");
  crosstrack =  Uvar(0,"km");
  alongtrack =  Uvar(0,"km");
  area = Uvar(0,"km km");
  radar_geom_factor = 0.0;
  echo_power = 0.0;
  radar_geom_factor_mirror = 0.0;
  good_bin = 0;
  ML_good_bin=0;
  thermal_snr_SL = 0;
  amb_ratio_SL = 0;  
  thermal_snr_ML = 0;
  amb_ratio_ML = 0;  
  prf_=Uvar(0,"Hz");
  pri_=Uvar(0,"s");
  pbw_=Uvar(0,"Hz");
  pulsegate_=Uvar(0,"km");
  }


AmbiguityFNN::AmbiguityFNN()
  :range("range"),
   doppler("dop"),
   oneway_beam_gaindB("oneway_gain"),
   incidenceangle("incidence_angle"),
   lat("latitude"),
   lon("longitude"),
   crosstrack("crosstrack"),
   alongtrack("alongtrack"),
   area("area"),
   radar_geom_factor("radar_geom_factor"),
   echo_power("echo_power"),
   radar_geom_factor_mirror("radar_geom_factor mirror"),
   good_bin("good_bin"),
   ML_good_bin("ML_good_bin"),
   thermal_snr_SL("thermal_snr_SL"),
   amb_ratio_SL("amb_ratio_SL"),
   thermal_snr_ML("thermal_snr_ML"),
   amb_ratio_ML("amb_ratio_ML"),
   time_set_(false),
   target_set_(false),
   state_set_(false),
   beam_set_(false), 
   config_set_(false),
   mirror_geom_set_(false),
   beam_number_set_(false),
   beam_max_gain_set_(false),
   process_window_set_(false),
   radar_geom_factor_set_(false),
   peak_radar_geom_factor_set_(false),
   usable_area_set_(false),
   ML_usable_area_set_(false),
   track_frame_set_(false),
   nlook_(1)
  {
  //resize with Ndat_ and Ngrid_ values
  range.resize(Ndat_);
  doppler.resize(Ndat_);
  oneway_beam_gaindB.resize(Ndat_,Ndat_);
  incidenceangle.resize(Ndat_,Ndat_);
  lat.resize(Ngrid_,Ngrid_);
  lon.resize(Ngrid_,Ngrid_);
  crosstrack.resize(Ngrid_,Ngrid_);
  alongtrack.resize(Ngrid_,Ngrid_);
  area.resize(Ndat_,Ndat_);
  radar_geom_factor.resize(Ndat_,Ndat_,Npatch_);
  echo_power.resize(Ndat_,Ndat_);
  radar_geom_factor_mirror.resize(Ndat_,Ndat_,Npatch_); 
  good_bin.resize(Ndat_,Ndat_);
  ML_good_bin.resize(Ndat_,Ndat_); 
  thermal_snr_SL.resize(Ndat_,Ndat_);
  amb_ratio_SL.resize(Ndat_,Ndat_);  
  thermal_snr_ML.resize(Ndat_,Ndat_);
  amb_ratio_ML.resize(Ndat_,Ndat_);

  //initialize result variables
  range = Uvar(0,"km");
  doppler = Uvar(0,"Hz");
  oneway_beam_gaindB = -100.0;
  incidenceangle = Uvar(0,"rad");
  lat = Uvar(0,"rad");
  lon = Uvar(0,"rad");
  crosstrack =  Uvar(0,"km");
  alongtrack =  Uvar(0,"km");
  area = Uvar(0,"km km");
  radar_geom_factor = 0.0;
  echo_power = 0.0;
  radar_geom_factor_mirror = 0.0;
  good_bin = 0;
  ML_good_bin=0;
  thermal_snr_SL = 0;
  amb_ratio_SL = 0;
  thermal_snr_ML = 0;
  amb_ratio_ML = 0;
  prf_=Uvar(0,"Hz");
  pri_=Uvar(0,"s");
  pbw_=Uvar(0,"Hz");
  pulsegate_=Uvar(0,"km");
  //cout<<"Warning:: process window's center is located at the boresight "<<endl;
  }



//-----------------------
//config(): read in lambda and mirror ambiguity 
//-----------------------
void AmbiguityFNN::config(Config& cfg)
{  
  if(config_set_) ErrorMessage("config is set").throwMe();
  lambda_ = speed_light/carrier_frequency;
  string mirror_option=cfg.str("calculate_mirror_ambiguity");
  if (mirror_option =="Yes" || mirror_option=="yes")
    {
    mirror_geom_set_ = true;
    }
  else
    {
    mirror_geom_set_ = false;
    }

  muhleman_k1_= cfg["Muhleman_backscatt_model_k1"].getInUnits("");
  muhleman_k2_= cfg["Muhleman_backscatt_model_k2"].getInUnits("");
  signal_to_ambdB_ =cfg["signal_to_amb_ratiodB"].getInUnits(""); 
  noise_equivalent_sigma0dB_ = cfg["noise_equivalent_sigma0dB"].getInUnits("");
  min_gaindB_ = cfg["min_oneway_gaindB_wrt_peak"].getInUnits("");
  roll_direction_=cfg.str("look_direction");
  config_set_ = true;
 
 
}


//--------------------------------
//setTarget(const string& target_name)
//--------------------------------
void AmbiguityFNN::setTarget(const string& target_name,const Frame& target_frame)
  {
    if(target_set_) ErrorMessage("target is set").throwMe();
    target_name_ = target_name;
    target_frame_ = target_frame;
    target_set_ = true;
  }
//-----------------------------
//setState
//-------------------------------
void AmbiguityFNN::setState(const StateVector& sc_state){
  if(time_set_){
    if(t_ != sc_state.time()) 
      ErrorMessage("time and statetime do not match").throwMe();
    else
      state_ = sc_state;
  }
  else{
    t_ = sc_state.time();
    time_set_ = true;
    state_ = sc_state;
  }
  state_set_ = true;
}

//------------------------------
//compute state
//-------------------------------
void AmbiguityFNN::computeState()
{
  if(!time_set_) ErrorMessage("time is not set").throwMe();
  if(!target_set_) ErrorMessage("target and target frame are not set").throwMe();
  target_frame_.ephemeris(state_,"Cassini",t_,"NONE");
  state_set_= true;
}

//--------------------------------
//setBeam(const Beam& beam)
//--------------------------------
void AmbiguityFNN::setBeam(const Beam& beam)
  {
  beam_ = beam;
  beam_num_ = beam_.getBeamNumber();
  beam_max_gain_ = beam_.getMaxGain();
  beam_set_ = true;
  beam_number_set_ = true;
  beam_max_gain_set_ = true;
  }

//----------------------------------
//setProcesswindow
//---------------------------------
 void AmbiguityFNN::setProcessWindow(const Uvar& central_range,
				     const Uvar& central_doppler,
				     const Uvar& prf, const Uvar& pbw,
				     const Uvar& pulsegate)
  { 
    center_range_ = central_range;
    center_doppler_ = central_doppler;
    prf_ = prf;
    pri_ = 1/prf_;
    pbw_ = pbw;
    pulsegate_ = pulsegate;
    process_window_set_ = true;
  }

//--------------------------------
//{set} track frame to calculate along/cross in a fixed frame
// named as "track frame"
//-------------------------------
void AmbiguityFNN::setTrackFrame(const Frame& ftrack)
  {
    ftrack_ = ftrack;
    track_frame_set_ = true;
  }

//-----------------------
//{get,set}Time(const Time& t)
//----------------------
Time  AmbiguityFNN::getTime() const
  {
  if(!time_set_)
    {
    ErrorMessage e("No time is set");
    e.throwMe();
    }
  return(t_);
  }

void AmbiguityFNN::setTime(const Time& t)
  {
    if(time_set_)
      {
	ErrorMessage("Time has been set").throwMe();
      }
    t_ = t;
    time_set_ = true;
  }
//---------------------------
//{get,set}Beamnumber()
//---------------------------
unsigned int AmbiguityFNN::getBeamNumber() const
  {
  if(!beam_number_set_)
    {
    ErrorMessage e("AmbiguityFNN::getBeamnumber: No beam number  is set");
    e.throwMe();
    }
  return(beam_num_);
  }

void  AmbiguityFNN::setBeamNumber(const unsigned int& beam_num) 
  {
  if(beam_number_set_)
    {
    ErrorMessage e("AmbiguityFNN::getBeamnumber:  beam number  is set");
    e.throwMe();
    }
  beam_num_ = beam_num;
  beam_number_set_ = true;
  }

//--------------------------------
//{get}beam's max gain
//--------------------------------
double AmbiguityFNN::getBeamMaxGain() const
  {
  if (!beam_max_gain_set_)
    {
    ErrorMessage e("AmbiguityFNN::getMaxGain: Beam max gain is not set");
    e.throwMe();
    }
  return(beam_max_gain_);
  }

//-----------------------------
//{get}Gridsize(): actual data size Ndat_
//-----------------------------
unsigned int AmbiguityFNN::getGridsize() const
  {  
  return(Ndat_);
  }

//-------------
//{get,set}NumberofLooks()
//-------------- 
unsigned int AmbiguityFNN::getNumberofLooks() const
{
  return(nlook_);
}

void AmbiguityFNN::setNumberofLooks(const unsigned int& nlook)
{
  nlook_ = nlook;
}

//------------------------------
// range/doppler index
//-----------------------------
unsigned int AmbiguityFNN::RangeIndex(const Uvar& range_value)
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
  if (range_value< range(0) || range_value > range(Ndat_-1))
    {
      ErrorMessage e("AmbiguityFNN::getRangeIndex:Range index is out of range");
      e.throwMe();
    }
  unsigned int range_index = (unsigned int) 
    ((range_value - range(0))/(range(1)- range(0))).getInUnits("");  
  return (range_index);
  }

unsigned int AmbiguityFNN::DopplerIndex(const Uvar& doppler_value)
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
    if (doppler_value< doppler(0) || doppler_value > doppler(Ndat_-1)) 
      {
	ErrorMessage("Dop index is out of range").throwMe();
      }
    unsigned int doppler_index = (unsigned int)
      ((doppler_value-doppler(0))/(doppler(1) - doppler(0))).getInUnits("");
   
    return(doppler_index);    
  }

//-----------------------------
//get range /doppler interval of the grid
//-----------------------------
Uvar AmbiguityFNN::RangeInterval()
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
  return(range(1) - range(0));
  }

Uvar AmbiguityFNN::DopplerInterval()
  {  
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
    return(doppler(1) - doppler(0));
  }

//------------------------------
//{get} center patch number
//----------------------------
unsigned int AmbiguityFNN::getCenterPatchNumber() const
{
  return(Ncenter_patch_);
}

//-----------------------------
//get crosstrack extent
//-----------------------------
void AmbiguityFNN:: getCrosstrackExtent(Uvar& cross_min,
					Uvar& cross_max, 
					Uvar& cross_length)

  {
    if(!ML_usable_area_set_)
      ErrorMessage("AmbiguityFNN:No ML usable area set").throwMe();
    cross_min= ML_crosstrack_extent_min_;
    cross_max= ML_crosstrack_extent_max_;
    cross_length=ML_crosstrack_length_;
  }
//---------------------------
//calAmbiguityFNN()
//---------------------------
void AmbiguityFNN::calAmbGeometry()
  {
  if(!time_set_)
    {
    ErrorMessage("AmbiguityFNN::calAmbiguityFNN: no time is set").throwMe();
    }
  if(!target_set_)
    {
    ErrorMessage("AmbiguityFNN::calAmbiguityFNN: no target is set").throwMe();
    }
  if(!state_set_)
    {
    ErrorMessage("AmbiguityFNN::calAmbiguityFNN: no state is set").throwMe();
    }

  if(!beam_set_)
    {
    ErrorMessage("AmbiguityFNN::calAmbiguityFNN: no beam is set").throwMe();
    }
  if(!config_set_)
    {
    ErrorMessage("AmbiguityFNN::calAmbiguityFNN: no config file is read in").throwMe();
    }
  if(!process_window_set_)
    {
      ErrorMessage("AmbiguityFNN:calAmbiguityFNN: no process window has been set up").throwMe();
    }
  if(!track_frame_set_)
    {
      ErrorMessage("AmbiguityFNN.cpp:calAmbiguityFNN: no track frame set").throwMe();
    }
  
  //----------------------
  //Need three direction vector to find a correct solution
  //-----------------------
  DirectionVector dir_pos = state_.position();
  DirectionVector dir_vel = state_.velocity();
  DirectionVector dir_cross = cross(dir_vel,dir_pos); 


  //-----------------------------------------------
  //determine roll direction using beam3 boresight direction
  //------------------------------------------------
  Frame fbeam3("CASSINI_RADAR_3","Cassini");
  DirectionVector boresight3("boresight",fbeam3,t_,0,0,1);
  double roll_direction=dot(boresight3,dir_cross);
  //if (roll_direction > 0) cout<<"right looking "<<endl;
  //else cout<<"left looking "<<endl;
  

  //---------------------------------
  //obtain boresight range and doppler
  //--------------------------------
  unsigned int beam_num= beam_.getBeamNumber();
  Frame fbeam("CASSINI_RADAR_" + toStr(beam_num),"Cassini");
  DirectionVector boresight("boresight",fbeam,t_,0,0,1);
  
  TargetGeom tg(t_);
  tg.setState(state_);
  tg.setLookDirection(boresight);
  tg.setTarget(target_name_,target_frame_);
  radius_ = tg.radius();
  Uvar altitude = tg.altitude();//altitude

  //-----------------------
  //grid variables(local variables)
  //----------------------
  Umat range_grid("range",Ngrid_,Npatch_);
  Umat doppler_grid("frequency",Ngrid_,Npatch_);

  Imat no_of_solution_grid("no_solution",Ngrid_*Ngrid_,Npatch_);
  Dmat oneway_beam_gain_grid("beam_gain",Ngrid_*Ngrid_,Npatch_);
  Umat incidenceangle_grid("thetai",Ngrid_*Ngrid_,Npatch_);
  Array2D<PositionVector> surface_intercept_grid("surface",Ngrid_*Ngrid_,Npatch_);

  Imat no_of_solution_mirror_grid("no_solution",Ngrid_*Ngrid_,Npatch_);
  Dmat oneway_beam_gain_mirror_grid("beam_gain",Ngrid_*Ngrid_,Npatch_);
  Umat  incidenceangle_mirror_grid("thetai",Ngrid_*Ngrid_,Npatch_);
  Array2D<PositionVector> surface_intercept_mirror_grid("surface",Ngrid_*Ngrid_,Npatch_);
  
  
  //initialize all the grid variables
 
  range_grid=Uvar(0,"km");
  doppler_grid=Uvar(0,"Hz");

  no_of_solution_grid=0;
  oneway_beam_gain_grid = 0.0; 
  incidenceangle_grid = Uvar(0,"rad");
 
  no_of_solution_mirror_grid=0;
  oneway_beam_gain_mirror_grid = 0.0; 
  incidenceangle_mirror_grid = Uvar(0,"rad");
 
  //-----------------------------
  //Pixel variables
  //-----------------------------
  Umat range_pixel("range",Ndat_,Npatch_);
  Umat doppler_pixel("frequency",Ndat_,Npatch_);

  D3D oneway_beam_gain_pixel("beam_gain",Ndat_,Ndat_,Npatch_);
  U3D area_pixel("area pixel",Ndat_,Ndat_,Npatch_);
  U3D incidenceangle_pixel("thetai",Ndat_,Ndat_,Npatch_);
  
  D3D oneway_beam_gain_mirror_pixel("beam_gain",Ndat_,Ndat_,Npatch_);
  U3D area_mirror_pixel("area mirror pixel",Ndat_,Ndat_,Npatch_);
  U3D incidenceangle_mirror_pixel("thetai",Ndat_,Ndat_,Npatch_);
  
  range_pixel=Uvar(0,"km");
  doppler_pixel=Uvar(0,"Hz");

  oneway_beam_gain_pixel = 0.0; 
  area_pixel = Uvar(0.0,"km km");
  incidenceangle_pixel = Uvar(0,"rad");

  oneway_beam_gain_mirror_pixel = 0.0;
  area_mirror_pixel = Uvar(0,"km km");
  incidenceangle_mirror_pixel = Uvar(0,"rad");
  



 

  //--------------------------------
  //Assign values: 
  // process window  - pbw/2 to pbw/2
  //                 -pulsegate/2 to pulsegate/2
  // nearby amb patches are separated by prf in the doppler axis and c*pri/2
  // in the range axis
  //---------------------------------
  Uvar delta_range = pulsegate_/2.0;
  Uvar delta_doppler = pbw_/2.0;
  for (unsigned int i_patch = 0 ; i_patch < Npatch_;++i_patch)
    {
      Uvar range_offset,doppler_offset; 
      switch(i_patch)
	{
	case 0:
	  range_offset = Uvar(0,"km");	
	  doppler_offset = -prf_;
	  break;
	case 1:
	  range_offset = speed_light * pri_/2.0;
	  doppler_offset = Uvar(0,"Hz");	
	  break;
	case Ncenter_patch_:
	  range_offset = Uvar(0,"km");
	  doppler_offset = Uvar(0,"Hz");
	  break;
	case 3:
	  range_offset = -speed_light * pri_/2.0;
	  doppler_offset = Uvar(0,"Hz");	
	  break;
	case 4:
	  range_offset = Uvar(0,"km");	
	  doppler_offset = prf_;
	  break;
	default://should not happen
	  ErrorMessage("Wrong patch assignment").throwMe();
	}
      
      for (unsigned int i_grid = 0; i_grid < Ngrid_;++i_grid)
	{
	  range_grid(i_grid,i_patch) = -delta_range 
	    + center_range_ + range_offset;

	  range_grid(i_grid,i_patch) += 2.0 * delta_range 
	    * double(i_grid)/double(Ngrid_ -1);

	  doppler_grid(i_grid,i_patch) = -delta_doppler
	    +center_doppler_+doppler_offset;

	  doppler_grid(i_grid,i_patch) += 2.0 * delta_doppler
	    *double(i_grid)/double(Ngrid_ -1 );
	}

      if(i_patch==Ncenter_patch_){
	for(unsigned int i_range=0;i_range<Ndat_;++i_range)
	  range(i_range) = (range_grid(i_range,i_patch) + range_grid(i_range+1,i_patch))/2.0;
	
	for(unsigned int j_dop=0;j_dop<Ndat_;++j_dop)
	  doppler(j_dop) = (doppler_grid(j_dop,i_patch) + doppler_grid(j_dop+1,i_patch))/2.0;
      }
    }
  
 


  //------------------------------------------------------------
  //This setup is needed to speed up the coordinate 
  //transformation from target_body_frame to beam_frame or track_frame
  // This needs to be done only once as long as time is not changed
  //------------------------------------------------------------------
  double a,b,c;
  DirectionVector target_to_beam_x,target_to_beam_y,target_to_beam_z;
  DirectionVector target_to_track_x,target_to_track_y,target_to_track_z;
 
  fbeam.axialVectors(target_frame_,
		     t_,
		     target_to_beam_x,
		     target_to_beam_y,
		     target_to_beam_z);
  ftrack_.axialVectors(target_frame_,
		       t_,
		       target_to_track_x,
		       target_to_track_y,
		       target_to_track_z);
  DirectionVector lookInBeamFrame("beam frame direction",fbeam,t_,0,0,1);
  DirectionVector dirToSurface("track frame direction",target_frame_,t_,0,0,1);
  double radius_in_km = radius_.getInUnits("km");
 
  
  //-------------- spacecraft state vector------
  DirectionVector pos = state_.position();
  DirectionVector vel = state_.velocity();
  DirectionVector ulook("DirectionVector from sc to surface",target_frame_,t_,0,0,1);
  DirectionVector ulook_mirror("mirrorLook from sc to surfact",target_frame_,t_,0,0,1);
  //---------------dummy variables to be used for solving RangeDoppler
  double p1,p2,p3, v1,v2,v3;
  double sol;
  int sign=1;

  p1 = pos[DirectionVector::X];
  p2 = pos[DirectionVector::Y];
  p3 = pos[DirectionVector::Z];
  
  v1 = vel[DirectionVector::X];
  v2 = vel[DirectionVector::Y];
  v3 = vel[DirectionVector::Z];
  
  int roll_direction_indicator = 0;
  if (roll_direction_ =="Right" || roll_direction_=="right")
    {
      roll_direction_indicator = -1;
    }
  else if (roll_direction_=="Left" || roll_direction_=="left")
    {
      roll_direction_indicator = 1;
    }
  else
    {
      ErrorMessage("TargetGeom::invalid roll direction").throwMe();
    }

  double u1,u2,u3;
  double A,B,C,D,E,F;
  double range_in_radius;
  double position_in_radius = (state_.position().magnitude()/radius_).getInUnits("");
  double compute_roll;
  

  //recycle Targetgeom obj tg0
  for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Ngrid_; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ngrid_; ++j_dop)
    {
    unsigned int bin_index = i_range * Ngrid_ + j_dop;	

    //----------------------------------------------------
    //  Instead of relying on calling targetgeom, 
    //let's solve everyting here
    //------------------------------------------------------

    range_in_radius = (range_grid(i_range,i_patch)/radius_).getInUnits("");
    A =  (1.0 - position_in_radius*position_in_radius 
	  - range_in_radius*range_in_radius)
      /(2.0 * range_in_radius * position_in_radius);
    B =  ((lambda_ * doppler_grid(j_dop,i_patch))
	  /(2.0 * state_.velocity().magnitude())).getInUnits("");

    C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
    D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
    E = (A - p1 *C)/p3;
    F = -(p2 + p1*D)/p3;
    
    sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);

    if (sol >= 0)
      {
	no_of_solution_grid(bin_index,i_patch) = 1;   	 
	u2 = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
	u1 = C + D*u2;
	u3 = E + F*u2;
	
	//set look vector from sc to surface in target frame
	ulook[DirectionVector::X] = u1;
	ulook[DirectionVector::Y] = u2;
	ulook[DirectionVector::Z] = u3;
	
	//--------------------------------
	//set surface intercept point
	//--------------------------------
	surface_intercept_grid(bin_index,i_patch) = state_.position();
	surface_intercept_grid(bin_index,i_patch) += range_grid(i_range,i_patch)*ulook;
	dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	
	compute_roll=dir_cross[DirectionVector::X]*dirToSurface[DirectionVector::X]
	  +dir_cross[DirectionVector::Y]*dirToSurface[DirectionVector::Y]
	  +dir_cross[DirectionVector::Z]*dirToSurface[DirectionVector::Z];
	if (compute_roll*roll_direction > 0)
	  {
	    sign = 1;
	  }
	else
	  {
	    //cout<<"alt solution "<<endl;
	    sign=-1;
	    u2 = ( -(C*D + E*F)+sqrt( sol))/(D*D + F*F + 1.0);
	    u1 = C + D*u2;
	    u3 = E + F*u2;
	    
	    //set look vector from sc to surface in target frame
	    ulook[DirectionVector::X] = u1;
	    ulook[DirectionVector::Y] = u2;
	    ulook[DirectionVector::Z] = u3;
	    
	    //--------------------------------
	    //set surface intercept point
	    //--------------------------------
	    surface_intercept_grid(bin_index,i_patch) = state_.position();
	    surface_intercept_grid(bin_index,i_patch) += range_grid(i_range,i_patch)*ulook;
	    dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	  }
	
	//---------------------------------
	//check the result with ordinary method
	// do not activate this part unless necessary
	// it will burn your cpu time very quickly
	//--------------------------------
	//TargetGeom tg0(t_);
	//tg0.setState(state_);
	//tg0.setTarget(target_name_,target_frame_);
	//tg0.setLookDirection(ulook);
	//cout<<"range diff "<< range_grid(i_range,i_patch) - tg0.range()<<endl;
	//cout<<"dop diff "<<doppler_grid(j_dop,i_patch) - tg0.doppler(lambda_)<<endl;
	//--------------------------------
	//set surface intercept point
	//--------------------------------
	surface_intercept_grid(bin_index,i_patch) = state_.position();
	surface_intercept_grid(bin_index,i_patch) += range_grid(i_range,i_patch)*ulook;
	dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	
	//----------------------
	//set incidence angle
	//------------------------
	a = dot(-ulook,dirToSurface);//reuse a
	incidenceangle_grid(bin_index,i_patch)=  Uvar(acos(a),"rad");
	
	if(i_patch ==Ncenter_patch_)
	  {
	    //calculate crosstrack and alongtrack for locations
	    // inside the process window
	    //-------------------------------------------
	    //Find directly the along/cross distances
	    //by using axialVectors between ftrack and ftitan
	    // surface_intercept_grid is defined in target frame
	    // in order to calculate along/cross, we need to transform 
	    //the surface vector
	    // in ftrack frame
	    //-------------------------------------------    
	    a=dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::X]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::X]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::X];
	    
	    b= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Y];
	    
	    c= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Z];
	    
	    //---------------------------------------------
	    //theta = -acos(c)+pi/2;//change to latitude
	    //phi = atan2(b,a);
	    //-------------------------------------------------
	    crosstrack(i_range,j_dop).setValue(radius_in_km * (pi/2.0-acos(c)));
	    alongtrack(i_range,j_dop).setValue(-radius_in_km *atan2(b,a));
	  }
     
	//---------------------------------------------------
	//direction vector look is defined in target frame
	// need to transform into beam frame
	//reuse double variables a,b, and c
	//ulook is defined in target frame
	// so need a rotation to beam frame
	//---------------------------------------------------
	a = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	
	b = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	
	c = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];
	
	lookInBeamFrame[DirectionVector::X] = a;
	lookInBeamFrame[DirectionVector::Y] = b;
	lookInBeamFrame[DirectionVector::Z] = c;
	oneway_beam_gain_grid(bin_index,i_patch)=
	  beam_.getGainOneWay(lookInBeamFrame);
	
	//debugging
	//
	//if(beam_num ==5 &&
	// j_dop == Ndat_/2 &&
	// i_patch ==Ncenter_patch_)
	// 
	//{
	//  cout<<"check beam number"<< beam_.getBeamNumber()<<endl;
	//  cout<<"direction vector "<<lookInBeamFrame<<endl;
	//  cout<<"i range jdop "<<i_range<<" "<<j_dop<<endl;
	//  cout<<"gain db "<<10.0*log(oneway_beam_gain_grid(bin_index,i_patch)/beam_.getMaxGain())/log(10.0)  <<endl;
	//}


	if(mirror_geom_set_)
	  {
	    no_of_solution_mirror_grid(bin_index,i_patch) = 1;
	    u2 = ( -(C*D + E*F) +double(sign)* sqrt( sol))
	      /(D*D + F*F + 1.0);
	    u1 = C + D*u2;
	    u3 = E + F*u2;
	    ulook_mirror[DirectionVector::X] = u1;
	    ulook_mirror[DirectionVector::Y] = u2;
	    ulook_mirror[DirectionVector::Z] = u3;
	    //-----------------------------
	    //mirror surface intercept
	    //----------------------------
	    surface_intercept_mirror_grid(bin_index,i_patch) = state_.position();
	    surface_intercept_mirror_grid(bin_index,i_patch)+= range_grid(i_range,i_patch) * ulook_mirror;
	    
	    //--------------------------
	    //incidence angle
	    //dirToSurface and ulook_mirror are  defined in target frame
	    //-------------------------
	    dirToSurface = surface_intercept_mirror_grid(bin_index,i_patch);//reuse to_surface vector
	    a = dot(-ulook_mirror,dirToSurface);
	    incidenceangle_mirror_grid(bin_index,i_patch) = Uvar(acos(a),"rad");
	    
	    //---------------------------
	    //Again manual transformation 
	    //-----------------------
	    a = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	    
	    b = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	    
	    c = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];
	    
	    lookInBeamFrame[DirectionVector::X] = a;
	    lookInBeamFrame[DirectionVector::Y] = b;
	    lookInBeamFrame[DirectionVector::Z] = c;
	    oneway_beam_gain_mirror_grid(bin_index,i_patch)
	      =beam_.getGainOneWay(lookInBeamFrame);  
	  }//only if mirror geom set is true
      }//if there is a solution
    else
      {
	no_of_solution_grid(bin_index,i_patch) = 0;//no surface intercept point
	no_of_solution_mirror_grid(bin_index,i_patch) = 0;   
      }//no solution   
    }//loop over patch, range, doppler
    
 
  
  //------------------------------------------------
  //Finalize data set into Ndat_ X Ndat_ matrix form
  //construct Ndat_ x Ndat_ data set
  //-------------------------------------------------

  DirectionVector dir_s1,dir_s2,dir_s3;
  double side_a,side_b,side_c,cos_c;

  I3D no_of_solution_pixel("no_solution",Ndat_,Ndat_,Npatch_);
  for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Ngrid_-1; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ngrid_-1; ++j_dop)
    {
    unsigned int bin_index = i_range * Ngrid_ + j_dop;	
    unsigned int bin_Up,  bin_Left,bin_Diag;
    bin_Up = (i_range+1)*Ngrid_ + j_dop;	
    bin_Left = i_range*Ngrid_ + (j_dop+1);
    bin_Diag = (i_range+1)*Ngrid_ + j_dop+1;
    if (no_of_solution_grid(bin_index,i_patch)
	*no_of_solution_grid(bin_Up,i_patch)
	*no_of_solution_grid(bin_Left,i_patch)
	*no_of_solution_grid(bin_Diag,i_patch) == 1)
      {
      no_of_solution_pixel(i_range,j_dop,i_patch) = 1;

      dir_s1 = surface_intercept_grid(bin_index,i_patch);
      dir_s2 = surface_intercept_grid(bin_Up,i_patch);
      dir_s3 = surface_intercept_grid(bin_Left,i_patch);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_in_km* acos(a);
      side_b = radius_in_km* acos(b);
      side_c = radius_in_km* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //We calculate only half of the rectangle defined by
      // four corner values of position vectors and then multiply it by 2
      // to have the total area.
      area_pixel(i_range,j_dop,i_patch).setValue(2.0*0.5*sin(acos(cos_c))*side_a*side_b);
     
      //--------------------
      //cell range and doppler
      //-----------------
      range_pixel(i_range,i_patch)= (range_grid(i_range,i_patch) 
			       + range_grid(i_range + 1,i_patch))/2.0;
      doppler_pixel(j_dop,i_patch)=(doppler_grid(j_dop,i_patch)
			      +doppler_grid(j_dop+1,i_patch))/2.0;
      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain_pixel(i_range,j_dop,i_patch) = 
	(oneway_beam_gain_grid(bin_index,i_patch)
	 +oneway_beam_gain_grid(bin_Up,i_patch)
	 +oneway_beam_gain_grid(bin_Left,i_patch)
	 +oneway_beam_gain_grid(bin_Diag,i_patch))
	/4.0;
      
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle_pixel(i_range,j_dop,i_patch) = 
	(incidenceangle_grid(bin_index,i_patch) 
	 + incidenceangle_grid(bin_Up,i_patch)
	 +incidenceangle_grid(bin_Left,i_patch)
	 +incidenceangle_grid(bin_Diag,i_patch))/4.0;
      
      if (i_patch == Ncenter_patch_)
	{//keep center patch(targeted area) values
	  area(i_range,j_dop)= area_pixel(i_range,j_dop,i_patch);
	  oneway_beam_gaindB(i_range,j_dop)=
	    10.0*
	    log(oneway_beam_gain_pixel(i_range,j_dop,i_patch)/beam_max_gain_)
	    /log(10.0) ;
	  incidenceangle(i_range,j_dop)=
	    incidenceangle_pixel(i_range,j_dop,i_patch); 
	}
      }
    else
      {
	no_of_solution_pixel(i_range,j_dop,i_patch) = 0;
      }//no solution, return zero
    }//Loop over range/dop values


  for (unsigned int i_patch = 0; i_patch< Npatch_;++i_patch)
  for (unsigned int i = 0 ; i < Ndat_;++i)
  for(unsigned int j = 0; j < Ndat_;++j)
    {
      if (no_of_solution_pixel(i,j,i_patch) == 0)
	{
	  radar_geom_factor(i,j,i_patch) = 0.0;
	}
      else
	{
	  Uvar X = oneway_beam_gain_pixel(i,j,i_patch)
	    *oneway_beam_gain_pixel(i,j,i_patch)*lambda_*lambda_
	    * muhleman_backscatter(muhleman_k1_,
				   muhleman_k2_,
				   incidenceangle_pixel(i,j,i_patch))
	    *area_pixel(i,j,i_patch)/pow(range_pixel(i,i_patch),4);
	  radar_geom_factor(i,j,i_patch)=X.getInUnits("");		
	}     
    }

  radar_geom_factor_set_= true;
  //----------------------------------------------
  //take care of mirror amb  
  // if mirror_geom_set_ is false, return
  //----------------------------------------------
  if (!mirror_geom_set_) {return;}
  I3D no_of_solution_mirror_pixel("no_solution",Ndat_,Ndat_,Npatch_);
  for (unsigned int i_patch = 0; i_patch <Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Ngrid_-1; i_range++)
  for (unsigned int j_dop =0; j_dop<Ngrid_-1; j_dop++)
    {
    unsigned int bin_index = i_range * Ngrid_ + j_dop;	
    unsigned int bin_Up,  bin_Left,bin_Diag;
    bin_Up = (i_range+1)*Ngrid_ + j_dop;	
    bin_Left = i_range*Ngrid_ + (j_dop+1);
    bin_Diag = (i_range+1)*Ngrid_ + j_dop+1;
    if (no_of_solution_mirror_grid(bin_index,i_patch)
	*no_of_solution_mirror_grid(bin_Up,i_patch)
	*no_of_solution_mirror_grid(bin_Left,i_patch)
	*no_of_solution_mirror_grid(bin_Diag,i_patch) == 1)
      {
      no_of_solution_mirror_pixel(i_range,j_dop,i_patch) = 1;
      //------------------------------------
      //find Area element covered by (range,dop)
      //divide an area enclosed by two range and dop values
      // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
      //into two triangles and calclate their areas and add them
      //together.
      //-------------------------------------  
	    
      dir_s1 = surface_intercept_mirror_grid(bin_index,i_patch);
      dir_s2 = surface_intercept_mirror_grid(bin_Up,i_patch);
      dir_s3 = surface_intercept_mirror_grid(bin_Left,i_patch);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_in_km* acos(a);
      side_b = radius_in_km* acos(b);
      side_c = radius_in_km* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //multiply 2 for the other triangle
      //remember: we are calculating the area of a rectangle defined by
      // four corners
      // save output to area_mirror!
      area_mirror_pixel(i_range,j_dop,i_patch).setValue(
			     2.0*0.5*sin(acos(cos_c))*side_a*side_b);
      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain_mirror_pixel(i_range,j_dop,i_patch) = 
	(oneway_beam_gain_mirror_grid(bin_index,i_patch)
	 +oneway_beam_gain_mirror_grid(bin_Up,i_patch)
	 +oneway_beam_gain_mirror_grid(bin_Left,i_patch)
	 +oneway_beam_gain_mirror_grid(bin_Diag,i_patch))
	/4.0;
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle_mirror_pixel(i_range,j_dop,i_patch) = 
      (incidenceangle_mirror_grid(bin_index,i_patch) 
       + incidenceangle_mirror_grid(bin_Up,i_patch)
       +incidenceangle_mirror_grid(bin_Left,i_patch)
       +incidenceangle_mirror_grid(bin_Diag,i_patch))/4.0;      
      }
    else
      {
	no_of_solution_mirror_pixel(i_range,j_dop,i_patch) =0;
      }//no solution, return zero
    }//Loop over range/dop values


  for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
  for (unsigned int i = 0 ; i < Ndat_;++i)
  for (unsigned int j = 0; j < Ndat_;++j)
    {
      if (no_of_solution_mirror_pixel(i,j,i_patch) == 0)
	{
	 radar_geom_factor_mirror(i,j,i_patch) =0.0;
	} 
      else
	{
	  Uvar X_mirror= oneway_beam_gain_mirror_pixel(i,j,i_patch) 
	    * oneway_beam_gain_mirror_pixel(i,j,i_patch) * lambda_ * lambda_
	    *muhleman_backscatter(muhleman_k1_
				  ,muhleman_k2_,
				  incidenceangle_mirror_pixel(i,j,i_patch))
	    * area_mirror_pixel(i,j,i_patch)/ pow(range_pixel(i,i_patch),4);
	  radar_geom_factor_mirror(i,j,i_patch) =X_mirror.getInUnits("");
	}
    }
  }

//-----------------------------------
//MaxRadarGeomFactor():return max value of radar geom factor
//---------------------------------
double AmbiguityFNN::MaxRadarGeomFactor()
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("AmbiguityFNN.cpp:MaxRadarGeomFactor:radar geom factor is not set").throwMe(); 
      }

    if (peak_radar_geom_factor_set_) {return (peak_radar_geom_factor_);}
    peak_radar_geom_factor_ = radar_geom_factor(0,0,Ncenter_patch_);

    for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
    for (unsigned int i = 0; i < Ndat_;++i)
    for (unsigned int j = 0; j < Ndat_;++j)
      {
	if (peak_radar_geom_factor_< radar_geom_factor(i,j,i_patch))
	  peak_radar_geom_factor_ = radar_geom_factor(i,j,i_patch);
	if (peak_radar_geom_factor_< radar_geom_factor_mirror(i,j,i_patch))
	  peak_radar_geom_factor_ = radar_geom_factor_mirror(i,j,i_patch);
      }
    return(peak_radar_geom_factor_);
  }


//-------------------------------------------
//Calculate usable area for a given prf and pulse gate
//-------------------------------------------
Uvar AmbiguityFNN::calUsableArea(const Uvar& X0, 
				 const Uvar& Pn,     
				 const Uvar& x_res,
				 const Uvar& rg_res)
				
{
    if (!radar_geom_factor_set_)
      {
	ErrorMessage e("AmbiguityFNN::calUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
   
    usable_area_SL = Uvar(0,"km km");
    //---------------
    //examine data values falling inside 
    //range_start/range_end, dop_start/dop_end
    //-----------------  
    good_bin = 0;//reset good bin indicator
    for (unsigned int i_range = 0 ; i_range <Ndat_;++i_range)
    for (unsigned int j_dop =0 ; j_dop < Ndat_ ;++j_dop)
      {
	double amb_radar_geom =0.0;
	for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
	  {
	    amb_radar_geom += radar_geom_factor(i_range,j_dop,i_patch);
	    if (mirror_geom_set_)
	      {
		amb_radar_geom += radar_geom_factor_mirror(i_range,j_dop,i_patch);
	      }
	  }
	amb_radar_geom -= radar_geom_factor(i_range,j_dop,Ncenter_patch_);
	
	//------------------------------------------
	//calculate amb ratio and noise equivalent sigma0
	//-----------------------------------------
	double ambratiodB = 10.0 
	  * log(radar_geom_factor(i_range,j_dop,Ncenter_patch_)
		/amb_radar_geom)/log(10.0);
	Uvar area_ratio= x_res * rg_res/area(i_range,j_dop);
	double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
	echo_power(i_range,j_dop) = X0
	  *radar_geom_factor(i_range,j_dop,Ncenter_patch_)
	  *area_ratio.getInUnits("")
	  /muhleman_backscatter(muhleman_k1_,
				muhleman_k2_,
				incidenceangle(i_range,j_dop));
	double signal_in_Watt = 
	  echo_power(i_range,j_dop).getInUnits("kg km km/(s s s)");
	double nesigma0dB=10.0
	  * log(noise_in_Watt/signal_in_Watt)/log(10.0);

	//------
	//Store thermal snr and amb ratio 
	//------
	thermal_snr_SL(i_range,j_dop) = -nesigma0dB;
	amb_ratio_SL(i_range,j_dop) = ambratiodB;
	
	//------------------------------------
	//IN case there is no solution (surface intercept point)
	//-----------------------------------
	if (radar_geom_factor(i_range,j_dop,Ncenter_patch_)==0.0)
	  {
	    amb_ratio_SL(i_range,j_dop) = -100.0;
	    thermal_snr_SL(i_range,j_dop) = -100.0;
	  }
 
	//-------------
	//Pick usable area satisfying all the following conditions
	//-------------
	if (amb_ratio_SL(i_range,j_dop).getInUnits("") >=signal_to_ambdB_ &&
	    thermal_snr_SL(i_range,j_dop) >= -noise_equivalent_sigma0dB_ &&
	    oneway_beam_gaindB(i_range,j_dop).getInUnits("") >= min_gaindB_)
	  { 
	    usable_area_SL+=area(i_range,j_dop);
	    good_bin(i_range,j_dop) = 1;
	  }
      }//loop over i_range, j_dop
    usable_area_set_ = true;
    //cout<<"single look area "<< usable_area_SL;
    return(usable_area_SL);
  }


//---------------------------------------
//Calculate usable area for multilookings
//--------------------------------------
Uvar AmbiguityFNN::calMultilookUsableArea(const Uvar& X0, 
					  const Uvar& Pn, 
					  const Uvar& x_res,
					  const Uvar& rg_res,
					  const Uvar& bpd,
					  const FloatVector& velocity,
					  const unsigned int& N_p)
  {
  if (!radar_geom_factor_set_)
      {
	ErrorMessage e("AmbiguityFNN::calMultilookUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
  
  
  usable_area_ML=Uvar(0.0,"km km");
  Uvar dop_interval = doppler(1) - doppler(0);
  unsigned int j_dop_center=DopplerIndex(center_doppler_);
  
  
  //--------------------------------
  //Now calculate number of looks for data set
  //encloded by (i_range_start_ML,j_dop_start_ML)
  //and (i_range_end_ML, j_dop_end_ML)
  //--------------------------------
  double cos_theta = 
    (center_doppler_ * lambda_ /(2.0 * velocity.magnitude())).getInUnits("");
  if (cos_theta > 1.0)
    {
      ErrorMessage e("AmbiguityFNN::cos theta value is larger than 1.0");
      e.throwMe();
    }
  Uvar ML_freq_bp = 2.0 * velocity.magnitude() 
                   * velocity.magnitude()*bpd/(lambda_*center_range_);
  ML_freq_bp *= (1.0-cos_theta * cos_theta);

  Uvar ML_freq_on = 1/(pri_ * double(N_p));

  //frequency resolution=1/(integration time(=pri * N_p))
 
  //cout<<"frequency shift "<<ML_freq_bp<<endl;
  //-------------- SH's simplied model produces the same result----------
  //double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  //Uvar ML_freq_bp = 
  //2.0 * velocity.magnitude()* Vst_bore * bpd*sin_theta/(lambda*range_bore);
  //note: Vst_bore = V * (1-cos(theta)^2)^0.5 --> same result!!!!!
  //for detailed information, see c:/ygim/radar_work/iso_doppler_profile
  //-------------------------------------------------------------------------

  //unsigned int ML_i_fon = 
  //(unsigned int) (ML_freq_on / dop_interval).getInUnits("");
  
  unsigned int ML_i_fbp = 
  (unsigned int) (ML_freq_bp / dop_interval).getInUnits("");
  
  
  if (ML_i_fbp ==0)
    {
      ErrorMessage("amb_param_scan:Need longer burst period to produce a shift in dop cell").throwMe();//infinite number of looks, unrealistic solution
      ML_i_fbp=1;
    }

  cout<<"ML_i_fbp "<< ML_i_fbp<<endl;
  //else if(ML_i_fbp > Ndat_)
  //{
  //  nlook_ = 1;
  //  //doppler frequency change after each burst is greater 
  //  //than process bandwidth
  //}
  //else
  //{
  //  nlook_=(unsigned int) Ndat_/ML_i_fbp;
  //}

  //calculate number of looks
  //usable process window > min_gain
  vector<Uvar> doppler_container;
  doppler_container.clear();
  for (unsigned int i_range = 0; i_range <Ndat_;++i_range)   
  for (unsigned int j_dop = 0; j_dop <Ndat_;++j_dop)
    {
      if(oneway_beam_gaindB(i_range,j_dop).getInUnits("") >= min_gaindB_)
	{
	  doppler_container.push_back(doppler(j_dop));
	}
    }

  if (doppler_container.size() == 0)
    {
      nlook_ = 1;
    }
  else if (doppler_container.size()==1)
    {
      nlook_ = 1;
    }
  else
    {
      sort(doppler_container.begin(),doppler_container.end());
      Uvar min_dop, max_dop;
      min_dop = doppler_container.front();
      max_dop = doppler_container.back();
      nlook_ = (unsigned int) round_double(((max_dop - min_dop)/ ML_freq_bp).getInUnits(""));
    }
  cout<<"nlook "<< nlook_<<endl; 
  //------
  //setting up mask
  //----
  Array1D<unsigned int> mask("mask",Ndat_);
  for (unsigned int j_dop = 0; j_dop < Ndat_;++j_dop)
    {
      mask(j_dop) = j_dop % ML_i_fbp;	
    }  
  

  ML_good_bin = 0;//reset good bin indicator
  thermal_snr_ML= 0.0;
  amb_ratio_ML = 0.0;
  vector<Uvar> crosstrack_extent;
  crosstrack_extent.clear();
  Uvar crosstrack_length("crosstrack length",0,"km");

  //------------------------------------------
  //Assume we process 3dB or 5dB points inside
  //------------------------------------------
  for (unsigned int i_range = 0; i_range <Ndat_;++i_range)
    {
      //-----------------------------
      //Here, I would like to calculate min signal_to_amb ration
      //after multilook
      //-----------------------------
      Uvec ML_signal_bin("ML signal bin", ML_i_fbp);
      Uvec ML_amb_bin("ML amb bin",ML_i_fbp); 
      Uvec ML_every_signal_bin("ML signal bin", ML_i_fbp);
      Uvec ML_every_amb_bin("ML amb bin",ML_i_fbp); 
      ML_signal_bin = 0.0;
      ML_amb_bin = 0.0;
      ML_every_signal_bin = 0.0;
      ML_every_amb_bin = 0.0;
      Dvec amb_ratio_inside_window("amb ratio inside window",Ndat_);     
      amb_ratio_inside_window= 0.0;
      for (unsigned int j_dop = 0; j_dop <Ndat_;++j_dop)
	{
	  double amb_radar_geom=0.0;
	  double signal =0.0;
	  for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
	    {
	      if (i_patch != Ncenter_patch_) 
		{
		  amb_radar_geom += radar_geom_factor(i_range,j_dop,i_patch);
		}
	      else
		{
		  signal = radar_geom_factor(i_range,j_dop,Ncenter_patch_);
		}
	      if (mirror_geom_set_)
		{
		  amb_radar_geom += radar_geom_factor_mirror(i_range,j_dop,i_patch);
		}
	    }

	  ML_every_signal_bin(mask(j_dop)) += signal;
	  ML_every_amb_bin(mask(j_dop))+= amb_radar_geom;

	  //if gain is larger than min gain
	  //need special treatment
	  if(oneway_beam_gaindB(i_range,j_dop) >= min_gaindB_)
	    {
	      ML_signal_bin(mask(j_dop))+= signal;
	      ML_amb_bin(mask(j_dop))+=amb_radar_geom;
	    }
	}//loop over j_dop

      //---------------------------------------
      //find multilook ambiguity ratio inside min beam gain
      //---------------------------------------
      vector<Uvar> amb_inside_window;
      amb_inside_window.clear();
      for (unsigned int i_fbp = 0; i_fbp < ML_i_fbp;++i_fbp)
	{
	  if(ML_signal_bin(i_fbp)!= 0.0)
	    {
	    amb_inside_window.push_back(10.0 * log(ML_signal_bin(i_fbp)/ML_amb_bin(i_fbp))/log(10.0));
	    }
	}

      //------------------
      //find lowest ambiguity ratio inside min beam gain
      //------------------
      Uvar  min_amb_ratio_inside_min_gain ;      
      if(amb_inside_window.size()==0) min_amb_ratio_inside_min_gain = 0;
      else
	{
	  sort(amb_inside_window.begin(),amb_inside_window.end());
	  min_amb_ratio_inside_min_gain=amb_inside_window.front();
	}
      
      
     

      //------------------------
      //Assign amb values
      //----------------------------
      //find min amb ratio inside min beam gain
      for (unsigned int j_dop =0; j_dop < Ndat_;++j_dop)
       {
	 Uvar area_ratio= x_res * rg_res/area(i_range,j_dop);
	 double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
	 double signal_in_Watt = 
	   (X0
	    *radar_geom_factor(i_range,j_dop,Ncenter_patch_)
	    *area_ratio.getInUnits("")
	    /muhleman_backscatter(muhleman_k1_,
				  muhleman_k2_,
				  incidenceangle(i_range,j_dop))
	    ).getInUnits("kg km km/(s s s)");
	 double nesigma0dB=10.0
	   * log(noise_in_Watt/signal_in_Watt/sqrt(double(nlook_)))/log(10.0);
	
	 //------
	 //Store thermal snr and amb ratio 
	 //------
	 thermal_snr_ML(i_range,j_dop) = -nesigma0dB;	 
	 //-----------------------------------
	 //Set amb ratio
	 //----------------------------------
	 if(oneway_beam_gaindB(i_range,j_dop) >= min_gaindB_)
	   {
	     amb_ratio_ML(i_range,j_dop)= min_amb_ratio_inside_min_gain;
	   }
	 else
	   {
	     amb_ratio_ML(i_range,j_dop)= 10.0 
	       * log(ML_every_signal_bin(mask(j_dop))/ML_every_amb_bin(mask(j_dop))/log(10.0));
	   }

 
	 //-------------
	 //Pick usable area satisfying all the following conditions
	 //-------------
	 if (amb_ratio_ML(i_range,j_dop).getInUnits("") >=signal_to_ambdB_ &&
	     thermal_snr_ML(i_range,j_dop) >= -noise_equivalent_sigma0dB_ &&
	     oneway_beam_gaindB(i_range,j_dop).getInUnits("") >= min_gaindB_)
	   { 
	     usable_area_ML+=area(i_range,j_dop);
	     ML_good_bin(i_range,j_dop) = 1;
	     crosstrack_extent.push_back(crosstrack(i_range,j_dop));
	     crosstrack_extent.push_back(crosstrack(i_range+1,j_dop));
	     if(j_dop==j_dop_center)
	       {
		 Uvar a=crosstrack(i_range+1,j_dop)-crosstrack(i_range,j_dop);
		 crosstrack_length +=Uvar(fabs(a.getInUnits("km")),"km");
		 //always take positve value	     
	       }
	   }
       }//loop over j_dop
    }//loop over i_range
  

  //------------------------------
  //Let's extract crosstrack min and max
  //-------------------------------
  ML_crosstrack_extent_min_ = Uvar(0,"km");
  ML_crosstrack_extent_max_ = Uvar(0,"km");
  ML_crosstrack_length_=crosstrack_length;
  if (crosstrack_extent.size() > 0)
    {//only if ML usable area is not zero
      sort(crosstrack_extent.begin(),crosstrack_extent.end());//sorting
      ML_crosstrack_extent_min_ = crosstrack_extent.front();//first element
      ML_crosstrack_extent_max_ = crosstrack_extent.back();//last element
    }
  ML_usable_area_set_ = true;
  //cout<<"ML look area "<< usable_area_ML;
  return(usable_area_ML);
  }


//--------------------------------
//return mirror_geom_set_
//---------------------------------
bool AmbiguityFNN::mirrorGeomSet()
  {
  return(mirror_geom_set_);
  }


//---------------------
//clear method
//---------------------
void AmbiguityFNN::clear()
{
  time_set_=false;
  state_set_=false;
  beam_set_=false;
  beam_number_set_=false;
  beam_max_gain_set_=false;
  process_window_set_=false;

  radar_geom_factor_set_=false;
  peak_radar_geom_factor_set_=false;
  usable_area_set_=false;
  ML_usable_area_set_=false;
}





const unsigned int SARAmb::Npatch_ = 5;
const unsigned int SARAmb::Ncenter_patch_ = 2;
/*
SARAmb::SARAmb(const Time& t)
  :surface_intercept_grid("surface"),
   surface_intercept_mirror_grid("surface"),
   time_set_(true),
   target_set_(false),
   state_set_(false),
   beam_set_(false), 
   config_set_(false),
   mirror_geom_set_(false),
   beam_number_set_(false),
   beam_max_gain_set_(false),
   process_window_set_(false),
   radar_geom_factor_set_(false),
   peak_radar_geom_factor_set_(false),
   usable_area_set_(false),
   ML_usable_area_set_(false),
   track_frame_set_(false),
   t_(t),
   nlook_(1)
  {
    beam_.clear();//empty the container
    beam_max_gain_.clear();
  }
*/


SARAmb::SARAmb()
  : surface_intercept_grid("surface"),
    surface_intercept_mirror_grid("surface"),
    time_set_(false),
    target_set_(false),
    state_set_(false),
    beam_set_(false), 
    config_set_(false),
    mirror_geom_set_(false),
    beam_number_set_(false),
    beam_max_gain_set_(false),
    process_window_set_(false),
    radar_geom_factor_set_(false),
    peak_radar_geom_factor_set_(false),
    usable_area_set_(false),
    ML_usable_area_set_(false),
    track_frame_set_(false),
    grid_set_(false),
    nlook_(1)
  {
   
    //-------------------
    //Beam parameters
    //-------------------
    beam_.clear();
    beam_max_gain_.clear();
  }

//---------------------
//destruction
//----------------------
SARAmb::~SARAmb()
  {
  if(grid_set_) deleteArrays();
  }

//-----------------------
//config(): read in lambda and mirror ambiguity 
//-----------------------
void SARAmb::config(Config& cfg)
{  
  if(config_set_) ErrorMessage("config is set").throwMe();
  lambda_ = (speed_light/carrier_frequency).getInUnits("km");
  string mirror_option=cfg.str("calculate_mirror_ambiguity");
  if (mirror_option =="Yes" || mirror_option=="yes")
    {
    mirror_geom_set_ = true;
    }
  else
    {
    mirror_geom_set_ = false;
    }

  muhleman_k1_= cfg["Muhleman_backscatt_model_k1"].getInUnits("");
  muhleman_k2_= cfg["Muhleman_backscatt_model_k2"].getInUnits("");
  signal_to_ambdB_ =cfg["signal_to_amb_ratiodB"].getInUnits(""); 
  noise_equivalent_sigma0dB_ = cfg["noise_equivalent_sigma0dB"].getInUnits("");
  min_gaindB_ = cfg["min_oneway_gaindB_wrt_peak"].getInUnits("");
  roll_direction_=cfg.str("look_direction");

  //set up beams
  for (unsigned int i = 0; i < 5;i++)
    {
      Beam beam(i+1,cfg,true);//use beam pattern file
      beam_max_gain_.push_back(beam.getMaxGain());
      beam_.push_back(beam);
    }

  int n = cfg.getInt("number_of_range_bins");
  if(n<1) ErrorMessage("Invalid number of range bins "+toStr(n));
  unsigned int Nrange = (unsigned int) n;
  
  n = cfg.getInt("number_of_dop_bins");
  if(n<1) ErrorMessage("Invalid number of doppler bins "+toStr(n));
  unsigned int Ndoppler=(unsigned int) n;
  
  createGrid(Nrange, Ndoppler);//only after config_set
  
  beam_max_gain_set_ = true;
  beam_set_ = true;
  config_set_ = true;
 
}

//------------------------
//crate grid
//-------------------------
void SARAmb::createGrid(const unsigned int& Nrange, const unsigned int& Ndoppler)
{
  Nrange_ = Nrange;
  Ndoppler_ = Ndoppler;
  Nrange_grid_ = Nrange_ + 1;
  Ndoppler_grid_ = Ndoppler_+1;
  createArrays();
}
//--------------------------------
//setTarget(const string& target_name)
//--------------------------------
void SARAmb::setTarget(const string& target_name,const Frame& target_frame)
  {
    if(target_set_) ErrorMessage("target is set").throwMe();
    target_name_ = target_name;
    target_frame_ = target_frame;
    TargetGeom tg;
    tg.setTarget(target_name);
    radius_ = tg.radius().getInUnits("km");
    target_set_ = true;
  }
//-----------------------------
//setState
//-------------------------------
void SARAmb::setState(const StateVector& sc_state){
  if(time_set_){
    if(t_ != sc_state.time()) 
      ErrorMessage("time and statetime do not match").throwMe();
    else
      state_ = sc_state;
  }
  else{
    t_ = sc_state.time();
    time_set_ = true;
    state_ = sc_state;
  }
  state_set_ = true;
}

//------------------------------
//compute state
//-------------------------------
void SARAmb::computeState()
{
  if(!time_set_) ErrorMessage("time is not set").throwMe();
  if(!target_set_) ErrorMessage("target and target frame are not set").throwMe();
  target_frame_.ephemeris(state_,"Cassini",t_,"NONE");
  state_set_= true;
}

//--------------------------------
//setBeamNumber(const Beam& beam)
//--------------------------------
void SARAmb::setBeamNumber(const unsigned int& beam_number)
  {
    if(beam_number < 1 || beam_number > 5)
      ErrorMessage("Invalid beam number").throwMe();
    beam_id_ = beam_number -1;
    beam_number_set_ = true;
  }
//---------------------------
//{get,set}Beamnumber()
//---------------------------
unsigned int SARAmb::getBeamNumber() const
  {
  if(!beam_number_set_)
    {
    ErrorMessage e("SARAmb::getBeamnumber: No beam number  is set");
    e.throwMe();
    }
  return(beam_id_+1);
  }

//--------------------------------
//{get}beam's max gain
//--------------------------------
double SARAmb::getBeamMaxGain() const
  {
    if(!beam_set_) ErrorMessage("Beam is not set").throwMe();
    if(!beam_max_gain_set_) ErrorMessage("max beam gain is not set").throwMe();
    if (!beam_number_set_)   ErrorMessage("beam is not set").throwMe();
    return(beam_max_gain_[beam_id_]);
  }


//----------------------------------
//setProcesswindow
//---------------------------------
 void SARAmb::setProcessWindow(const Uvar& central_range,
			       const Uvar& central_doppler,
			       const Uvar& lower_range_wrt_central_range,
			       const Uvar& upper_range_wrt_central_range,
			       const Uvar lower_doppler_wrt_central_doppler,
			       const Uvar upper_doppler_wrt_central_doppler,
			       const Uvar& prf)
				     
  { 
    if( (lower_range_wrt_central_range > upper_range_wrt_central_range) ||
	(lower_doppler_wrt_central_doppler > upper_doppler_wrt_central_doppler))
      {
	ErrorMessage("Lower sides of process window is larger than upper size of process window").throwMe();
      }						
    center_range_ = central_range.getInUnits("km");
    center_doppler_ = central_doppler.getInUnits("Hz");
    prf_ = prf.getInUnits("Hz");
    pri_ = 1/prf_;
    lower_range_ =lower_range_wrt_central_range.getInUnits("km");
    upper_range_= upper_range_wrt_central_range.getInUnits("km");
    lower_doppler_=  lower_doppler_wrt_central_doppler.getInUnits("Hz");
    upper_doppler_=  upper_doppler_wrt_central_doppler.getInUnits("Hz");
    process_window_set_ = true;
  }

//--------------------------------
//{set} track frame to calculate along/cross in a fixed frame
// named as "track frame"
//-------------------------------
void SARAmb::setTrackFrame(const Frame& ftrack)
  {
    ftrack_ = ftrack;
    track_frame_set_ = true;
  }

//-----------------------
//{get,set}Time(const Time& t)
//----------------------
Time  SARAmb::getTime() const
  {
  if(!time_set_)
    {
    ErrorMessage e("No time is set");
    e.throwMe();
    }
  return(t_);
  }

void SARAmb::setTime(const Time& t)
  {
    if(time_set_)
      {
	ErrorMessage("Time has been set").throwMe();
      }
    t_ = t;
    time_set_ = true;
  }




//-----------------------------
//{get}Gridsize(): actual data size Ndat_
//-----------------------------
unsigned int SARAmb::getRangeGridsize() const
{
   if(!grid_set_) ErrorMessage("no grid set").throwMe();
   return(Nrange_);
}
unsigned int SARAmb::getDopplerGridsize() const
{
   if(!grid_set_) ErrorMessage("no grid set").throwMe();
   return(Ndoppler_);
}
//-------------
//{get,set}NumberofLooks()
//-------------- 
unsigned int SARAmb::getNumberofLooks() const
{
  if(!ML_usable_area_set_) ErrorMessage("No multilook computation done").throwMe();
  return(nlook_);
}

void SARAmb::setNumberofLooks(const unsigned int& nlook)
{
  nlook_ = nlook;
}

//------------------------------
// range/doppler index
//----------------------------
unsigned int SARAmb::RangeIndex(const Uvar& range)
  {
    return(RangeIndex(range.getInUnits("km")));
  }

unsigned int SARAmb::RangeIndex(const double& range_in_km)
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
    if (range_in_km< range[0] || range_in_km > range[Nrange_-1])
    {
      ErrorMessage e("SARAmb::getRangeIndex:Range index is out of range");
      e.throwMe();
    }
  double step= range[Nrange_/2]-range[Nrange_/2-1];
  unsigned int range_index = (unsigned int) round_double(
    ((range_in_km - range[0])/(step))); 
  return (range_index);
  }

unsigned int SARAmb::DopplerIndex(const Uvar& doppler)
{
  return(DopplerIndex(doppler.getInUnits("Hz")));
}

unsigned int SARAmb::DopplerIndex(const double& doppler_in_Hz)
  {
   
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
    if (doppler_in_Hz< doppler[0] || doppler_in_Hz > doppler[Ndoppler_-1]) 
      {
	ErrorMessage("Dop index is out of range").throwMe();
      }
   
    double step= doppler[Ndoppler_/2]-doppler[Ndoppler_/2-1];
    unsigned int doppler_index = (unsigned int) round_double(
      ((doppler_in_Hz-doppler[0])/(step)));
    
    return(doppler_index);    
  }

//-----------------------------
//get range /doppler interval of the grid
//-----------------------------
Uvar SARAmb::RangeInterval()
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
  return(Uvar(range[1] - range[0],"km"));
  }

Uvar SARAmb::DopplerInterval()
  {  
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("Process window has not been set up").throwMe();
      }
    return(Uvar(doppler[1] - doppler[0],"Hz"));
  }

//------------------------------
//{get} center patch number
//----------------------------
unsigned int SARAmb::getCenterPatchNumber() const
{
  return(Ncenter_patch_);
}

//-----------------------------
//get crosstrack extent
//-----------------------------
void SARAmb:: getCrosstrackExtent(Uvar& cross_min,
				  Uvar& cross_max, 
				  Uvar& cross_length)

  {
    if(!ML_usable_area_set_)
      ErrorMessage("SARAmb:No ML usable area set").throwMe();
    cross_min= Uvar(ML_crosstrack_extent_min_,"km");
    cross_max= Uvar(ML_crosstrack_extent_max_,"km");
    cross_length=Uvar(ML_crosstrack_length_,"km");
  }
//---------------------------
//calSARAmb()
//---------------------------
void SARAmb::calAmbGeometry()
  {
  if(!time_set_)
    {
    ErrorMessage("SARAmb::calSARAmb: no time is set").throwMe();
    }
  if(!target_set_)
    {
    ErrorMessage("SARAmb::calSARAmb: no target is set").throwMe();
    }
  if(!state_set_)
    {
    ErrorMessage("SARAmb::calSARAmb: no state is set").throwMe();
    }

  if(!beam_set_)
    {
    ErrorMessage("SARAmb::calSARAmb: no beam is set").throwMe();
    }
  if(!config_set_)
    {
    ErrorMessage("SARAmb::calSARAmb: no config file is read in").throwMe();
    }
  if(!process_window_set_)
    {
      ErrorMessage("SARAmb:calSARAmb: no process window has been set up").throwMe();
    }
  if(!track_frame_set_)
    {
      ErrorMessage("SARAmb.cpp:calSARAmb: no track frame set").throwMe();
    }
  if(!grid_set_)
    {
      ErrorMessage("SARAmb.cpp:calSARAmbNo grid is set").throwMe();
    }  
  //----------------------
  //Need three direction vector to find a correct solution
  //-----------------------
  DirectionVector dir_pos = state_.position();
  DirectionVector dir_vel = state_.velocity();
  DirectionVector dir_cross = cross(dir_vel,dir_pos); 


  //-----------------------------------------------
  //determine roll direction using beam3 boresight direction
  //------------------------------------------------
  //Frame fbeam3("CASSINI_RADAR_3","Cassini");
  //DirectionVector boresight3("boresight",fbeam3,t_,0,0,1);
  //if (roll_direction > 0) cout<<"right looking "<<endl;
  //else cout<<"left looking "<<endl;
  

  //---------------------------------
  //obtain boresight range and doppler
  //--------------------------------
  unsigned int beam_num= getBeamNumber();
  Frame fbeam("CASSINI_RADAR_" + toStr(beam_num),"Cassini");
  DirectionVector boresight("boresight",fbeam,t_,0,0,1);
  double roll_direction=dot(boresight,dir_cross);

  //TargetGeom tg(t_);
  //tg.setState(state_);
  //tg.setLookDirection(boresight);
  //tg.setTarget(target_name_,target_frame_);

  
 
 
  //--------------------------------
  //Assign values: 
  //Expand edges of the process window by one pixel
  //---------------------------------
  double delta_r = (upper_range_ - lower_range_)/double(Nrange_grid_);
  upper_range_ += delta_r;
  lower_range_ -= delta_r;
  double delta_f = (upper_doppler_ - lower_doppler_)/double(Ndoppler_grid_);
  upper_doppler_ += delta_f;
  lower_doppler_ -=delta_f;
  double d_speed_light = speed_light.getInUnits("km/s");
   for (unsigned int i_patch = 0 ; i_patch < Npatch_;++i_patch)
    {
      double range_offset,doppler_offset; 
      switch(i_patch)
	{
	case 0:
	  range_offset = 0.0;	
	  doppler_offset = -prf_;
	  break;
	case 1:
	  range_offset = d_speed_light * pri_/2.0;
	  doppler_offset = 0.0;	
	  break;
	case Ncenter_patch_:
	  range_offset = 0;
	  doppler_offset =0;
	  break;
	case 3:
	  range_offset = -d_speed_light * pri_/2.0;
	  doppler_offset = 0;	
	  break;
	case 4:
	  range_offset = 0;	
	  doppler_offset = prf_;
	  break;
	default://should not happen
	  ErrorMessage("Wrong patch assignment").throwMe();
	}

      //  double delta_range = pulsegate_/2.0;
      //double  delta_doppler = pbw_/2.0;
     
      for (unsigned int i_grid = 0; i_grid < Nrange_grid_;++i_grid)
	{
	  range_grid[i_grid][i_patch] =  center_range_ + range_offset;
	  range_grid[i_grid][i_patch] += lower_range_
	  +(upper_range_ - lower_range_) * double(i_grid)/double(Nrange_grid_ -1);
	}
      for(unsigned int i_grid=0;i_grid<Ndoppler_grid_;++i_grid)
	{
	  doppler_grid[i_grid][i_patch] =center_doppler_+doppler_offset;
	  doppler_grid[i_grid][i_patch] += lower_doppler_ +
	    +(upper_doppler_- lower_doppler_)*double(i_grid)/double(Ndoppler_grid_ -1 );
	}
    }
  
 


  //------------------------------------------------------------
  //This setup is needed to speed up the coordinate 
  //transformation from target_body_frame to beam_frame or track_frame
  // This needs to be done only once as long as time is not changed
  //------------------------------------------------------------------
  double a,b,c;
  DirectionVector target_to_beam_x,target_to_beam_y,target_to_beam_z;
  DirectionVector target_to_track_x,target_to_track_y,target_to_track_z;
 
  fbeam.axialVectors(target_frame_,
		     t_,
		     target_to_beam_x,
		     target_to_beam_y,
		     target_to_beam_z);
  ftrack_.axialVectors(target_frame_,
		       t_,
		       target_to_track_x,
		       target_to_track_y,
		       target_to_track_z);
  DirectionVector lookInBeamFrame("beam frame direction",fbeam,t_,0,0,1);
  DirectionVector dirToSurface("track frame direction",target_frame_,t_,0,0,1);
 
 
  
  //-------------- spacecraft state vector------
  DirectionVector pos = state_.position();
  DirectionVector vel = state_.velocity();
  DirectionVector ulook("DirectionVector from sc to surface",target_frame_,t_,0,0,1);
  DirectionVector ulook_mirror("mirrorLook from sc to surfact",target_frame_,t_,0,0,1);
  //---------------dummy variables to be used for solving RangeDoppler
  double p1,p2,p3, v1,v2,v3;
  double sol;
  int sign=1;

  p1 = pos[DirectionVector::X];
  p2 = pos[DirectionVector::Y];
  p3 = pos[DirectionVector::Z];
  
  v1 = vel[DirectionVector::X];
  v2 = vel[DirectionVector::Y];
  v3 = vel[DirectionVector::Z];
  
  int roll_direction_indicator = 0;
  if (roll_direction_ =="Right" || roll_direction_=="right")
    {
      roll_direction_indicator = -1;
    }
  else if (roll_direction_=="Left" || roll_direction_=="left")
    {
      roll_direction_indicator = 1;
    }
  else
    {
      ErrorMessage("TargetGeom::invalid roll direction").throwMe();
    }

  double u1,u2,u3;
  double A,B,C,D,E,F;
  double range_in_radius;
  double position_in_radius = state_.position().magnitude().getInUnits("km")/radius_;
  double compute_roll;
  

  //recycle Targetgeom obj tg0
  double flatform_speed = state_.velocity().magnitude().getInUnits("km/s");
  for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Nrange_grid_; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ndoppler_grid_; ++j_dop)
    {
    unsigned int bin_index = i_range * Ndoppler_grid_ + j_dop;	

    //----------------------------------------------------
    //  Instead of relying on calling targetgeom, 
    //let's solve everyting here
    //------------------------------------------------------

    range_in_radius = range_grid[i_range][i_patch]/radius_;
    A =  (1.0 - position_in_radius*position_in_radius 
	  - range_in_radius*range_in_radius)
      /(2.0 * range_in_radius * position_in_radius);
    B =  (lambda_ * doppler_grid[j_dop][i_patch])
	  /(2.0 * flatform_speed);

    C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
    D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
    E = (A - p1 *C)/p3;
    F = -(p2 + p1*D)/p3;
    
    sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);

    if (sol >= 0)
      {
	no_of_solution_grid[bin_index][i_patch] = 1;   	 
	u2 = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
	u1 = C + D*u2;
	u3 = E + F*u2;
	
	//set look vector from sc to surface in target frame
	ulook[DirectionVector::X] = u1;
	ulook[DirectionVector::Y] = u2;
	ulook[DirectionVector::Z] = u3;
	
	//--------------------------------
	//set surface intercept point
	//--------------------------------
	surface_intercept_grid(bin_index,i_patch) = state_.position();
	surface_intercept_grid(bin_index,i_patch) += Uvar(range_grid[i_range][i_patch],"km")*ulook;
	dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	
	compute_roll=dir_cross[DirectionVector::X]*dirToSurface[DirectionVector::X]
	  +dir_cross[DirectionVector::Y]*dirToSurface[DirectionVector::Y]
	  +dir_cross[DirectionVector::Z]*dirToSurface[DirectionVector::Z];
	if (compute_roll*roll_direction > 0)
	  {
	    sign = 1;
	  }
	else
	  {
	    //cout<<"alt solution "<<endl;
	    sign=-1;
	    u2 = ( -(C*D + E*F)+sqrt( sol))/(D*D + F*F + 1.0);
	    u1 = C + D*u2;
	    u3 = E + F*u2;
	    
	    //set look vector from sc to surface in target frame
	    ulook[DirectionVector::X] = u1;
	    ulook[DirectionVector::Y] = u2;
	    ulook[DirectionVector::Z] = u3;
	    
	    //--------------------------------
	    //set surface intercept point
	    //--------------------------------
	    surface_intercept_grid(bin_index,i_patch) = state_.position();
	    surface_intercept_grid(bin_index,i_patch) += Uvar(range_grid[i_range][i_patch],"km")*ulook;
	    dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	  }
	
	//---------------------------------
	//check the result with ordinary method
	// do not activate this part unless necessary
	// it will burn your cpu time very quickly
	//--------------------------------
	//TargetGeom tg0(t_);
	//tg0.setState(state_);
	//tg0.setTarget(target_name_,target_frame_);
	//tg0.setLookDirection(ulook);
	//cout<<"range diff "<< range_grid(i_range,i_patch) - tg0.range()<<endl;
	//cout<<"dop diff "<<doppler_grid(j_dop,i_patch) - tg0.doppler(lambda_)<<endl;
	//--------------------------------
	//set surface intercept point
	//--------------------------------
	surface_intercept_grid(bin_index,i_patch) = state_.position();
	surface_intercept_grid(bin_index,i_patch) += Uvar(range_grid[i_range][i_patch],"km")*ulook;
	dirToSurface = surface_intercept_grid(bin_index,i_patch);//dir vector in target frame
	
	//----------------------
	//set incidence angle
	//------------------------
	a = dot(-ulook,dirToSurface);//reuse a
	incidenceangle_grid[bin_index][i_patch]=  acos(a);
	
	if(i_patch ==Ncenter_patch_)
	  {
	    //calculate crosstrack and alongtrack for locations
	    // inside the process window
	    //-------------------------------------------
	    //Find directly the along/cross distances
	    //by using axialVectors between ftrack and ftitan
	    // surface_intercept_grid is defined in target frame
	    // in order to calculate along/cross, we need to transform 
	    //the surface vector
	    // in ftrack frame
	    //-------------------------------------------    
	    a=dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::X]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::X]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::X];
	    
	    b= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Y];
	    
	    c= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Z];
	    
	    //---------------------------------------------
	    //theta = -acos(c)+pi/2;//change to latitude
	    //phi = atan2(b,a);
	    //-------------------------------------------------
	    crosstrack[i_range][j_dop]=radius_ * (pi/2.0-acos(c));
	    alongtrack[i_range][j_dop]=-radius_ *atan2(b,a);
	  }
     
	//---------------------------------------------------
	//direction vector look is defined in target frame
	// need to transform into beam frame
	//reuse double variables a,b, and c
	//ulook is defined in target frame
	// so need a rotation to beam frame
	//---------------------------------------------------
	a = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	
	b = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	
	c = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	  +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	  +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];
	
	lookInBeamFrame[DirectionVector::X] = a;
	lookInBeamFrame[DirectionVector::Y] = b;
	lookInBeamFrame[DirectionVector::Z] = c;
	oneway_beam_gain_grid[bin_index][i_patch]=
	  beam_[beam_num-1].getGainOneWay(lookInBeamFrame);
	
	//debugging
	//
	//if(beam_num ==5 &&
	// j_dop == Ndoppler_/2 &&
	// i_patch ==Ncenter_patch_)
	// 
	//{
	//  cout<<"check beam number"<< beam_.getBeamNumber()<<endl;
	//  cout<<"direction vector "<<lookInBeamFrame<<endl;
	//  cout<<"i range jdop "<<i_range<<" "<<j_dop<<endl;
	//  cout<<"gain db "<<10.0*log(oneway_beam_gain_grid(bin_index,i_patch)/beam_.getMaxGain())/log(10.0)  <<endl;
	//}


	if(mirror_geom_set_)
	  {
	    no_of_solution_mirror_grid[bin_index][i_patch] = 1;
	    u2 = ( -(C*D + E*F) +double(sign)* sqrt( sol))
	      /(D*D + F*F + 1.0);
	    u1 = C + D*u2;
	    u3 = E + F*u2;
	    ulook_mirror[DirectionVector::X] = u1;
	    ulook_mirror[DirectionVector::Y] = u2;
	    ulook_mirror[DirectionVector::Z] = u3;
	    //-----------------------------
	    //mirror surface intercept
	    //----------------------------
	    surface_intercept_mirror_grid(bin_index,i_patch) = state_.position();
	    surface_intercept_mirror_grid(bin_index,i_patch)+= Uvar(range_grid[i_range][i_patch],"km") * ulook_mirror;
	    
	    //--------------------------
	    //incidence angle
	    //dirToSurface and ulook_mirror are  defined in target frame
	    //-------------------------
	    dirToSurface = surface_intercept_mirror_grid(bin_index,i_patch);//reuse to_surface vector
	    a = dot(-ulook_mirror,dirToSurface);
	    incidenceangle_mirror_grid[bin_index][i_patch] = acos(a);
	    
	    //---------------------------
	    //Again manual transformation 
	    //-----------------------
	    a = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	    
	    b = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	    
	    c = ulook_mirror[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	      +ulook_mirror[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	      +ulook_mirror[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];
	    
	    lookInBeamFrame[DirectionVector::X] = a;
	    lookInBeamFrame[DirectionVector::Y] = b;
	    lookInBeamFrame[DirectionVector::Z] = c;
	    oneway_beam_gain_mirror_grid[bin_index][i_patch]
	      =beam_[beam_num-1].getGainOneWay(lookInBeamFrame);  
	  }//only if mirror geom set is true
      }//if there is a solution
    else
      {
	no_of_solution_grid[bin_index][i_patch] = 0;//no surface intercept point
	no_of_solution_mirror_grid[bin_index][i_patch] = 0;   
      }//no solution   
    }//loop over patch, range, doppler
    
 
  
  //------------------------------------------------
  //Finalize data set into Ndat_ X Ndat_ matrix form
  //construct Ndat_ x Ndat_ data set
  //-------------------------------------------------

  DirectionVector dir_s1,dir_s2,dir_s3;
  double side_a,side_b,side_c,cos_c;


  for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Nrange_grid_-1; ++i_range)
  for (unsigned int j_dop =0; j_dop<Ndoppler_grid_-1; ++j_dop)
    {
    unsigned int bin_index = i_range * Ndoppler_grid_ + j_dop;	
    unsigned int bin_Up,  bin_Left,bin_Diag;
    bin_Up = (i_range+1)*Ndoppler_grid_ + j_dop;	
    bin_Left = i_range*Ndoppler_grid_ + (j_dop+1);
    bin_Diag = (i_range+1)*Ndoppler_grid_ + j_dop+1;


    //--------------------
    //cell range and doppler (set whether or not this range,dop is on ground)
    //-----------------
    range_pixel[i_range][i_patch]= (range_grid[i_range][i_patch] 
				    + range_grid[i_range + 1][i_patch])/2.0;
    doppler_pixel[j_dop][i_patch]=(doppler_grid[j_dop][i_patch]
				   +doppler_grid[j_dop+1][i_patch])/2.0;


    if(i_patch== Ncenter_patch_){
	  range[i_range] = range_pixel[i_range][i_patch];
	  doppler[j_dop] = doppler_pixel[j_dop][i_patch];
    }

    if (no_of_solution_grid[bin_index][i_patch]
	*no_of_solution_grid[bin_Up][i_patch]
	*no_of_solution_grid[bin_Left][i_patch]
	*no_of_solution_grid[bin_Diag][i_patch] == 1)
      {
      no_of_solution_pixel[i_range][j_dop][i_patch] = 1;

      dir_s1 = surface_intercept_grid(bin_index,i_patch);
      dir_s2 = surface_intercept_grid(bin_Up,i_patch);
      dir_s3 = surface_intercept_grid(bin_Left,i_patch);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_* acos(a);
      side_b = radius_* acos(b);
      side_c = radius_* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //We calculate only half of the rectangle defined by
      // four corner values of position vectors and then multiply it by 2
      // to have the total area.
      area_pixel[i_range][j_dop][i_patch]=2.0*0.5*sin(acos(cos_c))*side_a*side_b;
     

      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain_pixel[i_range][j_dop][i_patch] = 
	(oneway_beam_gain_grid[bin_index][i_patch]
	 +oneway_beam_gain_grid[bin_Up][i_patch]
	 +oneway_beam_gain_grid[bin_Left][i_patch]
	 +oneway_beam_gain_grid[bin_Diag][i_patch])
	/4.0;
      
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle_pixel[i_range][j_dop][i_patch] = 
	(incidenceangle_grid[bin_index][i_patch] 
	 + incidenceangle_grid[bin_Up][i_patch]
	 +incidenceangle_grid[bin_Left][i_patch]
	 +incidenceangle_grid[bin_Diag][i_patch])/4.0;
      
      if (i_patch == Ncenter_patch_)
	{//keep center patch(targeted area) values
	  area[i_range][j_dop]= area_pixel[i_range][j_dop][i_patch];
	  oneway_beam_gaindB[i_range][j_dop]=
	    10.0*
	    log(oneway_beam_gain_pixel[i_range][j_dop][i_patch]/beam_max_gain_[beam_num-1])
	    /log(10.0) ;
	  incidenceangle[i_range][j_dop]=
	    incidenceangle_pixel[i_range][j_dop][i_patch]; 
	}
      }
    else
      {
	no_of_solution_pixel[i_range][j_dop][i_patch] = 0;
      }//no solution, return zero
    }//Loop over range/dop values


  for (unsigned int i_patch = 0; i_patch< Npatch_;++i_patch)
  for (unsigned int i = 0 ; i < Nrange_;++i)
  for(unsigned int j = 0; j < Ndoppler_;++j)
    {
      if (no_of_solution_pixel[i][j][i_patch] == 0)
	{
	  radar_geom_factor[i][j][i_patch] = 0.0;
	}
      else
	{
	  radar_geom_factor[i][j][i_patch] = oneway_beam_gain_pixel[i][j][i_patch]
	    *oneway_beam_gain_pixel[i][j][i_patch]*lambda_*lambda_
	    * muhleman_backscatter(muhleman_k1_,
				   muhleman_k2_,
				   Uvar(incidenceangle_pixel[i][j][i_patch],"rad"))
	    *area_pixel[i][j][i_patch]/pow(range_pixel[i][i_patch],4);
	}     
    }

  radar_geom_factor_set_= true;
  //----------------------------------------------
  //take care of mirror amb  
  // if mirror_geom_set_ is false, return
  //----------------------------------------------
  if (!mirror_geom_set_) {return;}
  
  for (unsigned int i_patch = 0; i_patch <Npatch_;++i_patch)
  for (unsigned int i_range = 0; i_range<Nrange_grid_-1; i_range++)
  for (unsigned int j_dop =0; j_dop<Ndoppler_grid_-1; j_dop++)
    {
    unsigned int bin_index = i_range * Ndoppler_grid_ + j_dop;	
    unsigned int bin_Up,  bin_Left,bin_Diag;
    bin_Up = (i_range+1)*Ndoppler_grid_ + j_dop;	
    bin_Left = i_range*Ndoppler_grid_ + (j_dop+1);
    bin_Diag = (i_range+1)*Ndoppler_grid_ + j_dop+1;
    if (no_of_solution_mirror_grid[bin_index][i_patch]
	*no_of_solution_mirror_grid[bin_Up][i_patch]
	*no_of_solution_mirror_grid[bin_Left][i_patch]
	*no_of_solution_mirror_grid[bin_Diag][i_patch] == 1)
      {
      no_of_solution_mirror_pixel[i_range][j_dop][i_patch] = 1;
      //------------------------------------
      //find Area element covered by (range,dop)
      //divide an area enclosed by two range and dop values
      // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
      //into two triangles and calclate their areas and add them
      //together.
      //-------------------------------------  
	    
      dir_s1 = surface_intercept_mirror_grid(bin_index,i_patch);
      dir_s2 = surface_intercept_mirror_grid(bin_Up,i_patch);
      dir_s3 = surface_intercept_mirror_grid(bin_Left,i_patch);

      a = dir_s1[DirectionVector::X] * dir_s2[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s2[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s2[DirectionVector::Z];
      
      b = dir_s1[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s1[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s1[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      c = dir_s2[DirectionVector::X] * dir_s3[DirectionVector::X]
	+dir_s2[DirectionVector::Y] * dir_s3[DirectionVector::Y]
	+dir_s2[DirectionVector::Z] * dir_s3[DirectionVector::Z];
      
      side_a = radius_* acos(a);
      side_b = radius_* acos(b);
      side_c = radius_* acos(c);

      cos_c = (side_a * side_a + side_b*side_b - side_c*side_c)
	       /(2.0*side_a*side_b);
      //multiply 2 for the other triangle
      //remember: we are calculating the area of a rectangle defined by
      // four corners
      // save output to area_mirror!
      area_mirror_pixel[i_range][j_dop][i_patch]=
			     2.0*0.5*sin(acos(cos_c))*side_a*side_b;
      //-----------------------------
      //store beam gain and sigm0
      //----------------------------
      oneway_beam_gain_mirror_pixel[i_range][j_dop][i_patch] = 
	(oneway_beam_gain_mirror_grid[bin_index][i_patch]
	 +oneway_beam_gain_mirror_grid[bin_Up][i_patch]
	 +oneway_beam_gain_mirror_grid[bin_Left][i_patch]
	 +oneway_beam_gain_mirror_grid[bin_Diag][i_patch])
	/4.0;
      //-----------------------------------
      //incidenceangle
      //----------------------------------
      incidenceangle_mirror_pixel[i_range][j_dop][i_patch] = 
      (incidenceangle_mirror_grid[bin_index][i_patch] 
       + incidenceangle_mirror_grid[bin_Up][i_patch]
       +incidenceangle_mirror_grid[bin_Left][i_patch]
       +incidenceangle_mirror_grid[bin_Diag][i_patch])/4.0;      
      }
    else
      {
	no_of_solution_mirror_pixel[i_range][j_dop][i_patch] =0;
      }//no solution, return zero
    }//Loop over range/dop values


  for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
  for (unsigned int i = 0 ; i < Nrange_;++i)
  for (unsigned int j = 0; j < Ndoppler_;++j)
    {
      if (no_of_solution_mirror_pixel[i][j][i_patch] == 0)
	{
	 radar_geom_factor_mirror[i][j][i_patch] =0.0;
	} 
      else
	{
	  radar_geom_factor_mirror[i][j][i_patch]= oneway_beam_gain_mirror_pixel[i][j][i_patch] 
	    * oneway_beam_gain_mirror_pixel[i][j][i_patch] * lambda_ * lambda_
	    *muhleman_backscatter(muhleman_k1_
				  ,muhleman_k2_,
				  Uvar(incidenceangle_mirror_pixel[i][j][i_patch],"rad"))
	    * area_mirror_pixel[i][j][i_patch]/ pow(range_pixel[i][i_patch],4);
	}
    }


  }

//-----------------------------------
//MaxRadarGeomFactor():return max value of radar geom factor
//---------------------------------
double SARAmb::MaxRadarGeomFactor()
  {
    if (!radar_geom_factor_set_)
      {
	ErrorMessage("SARAmb.cpp:MaxRadarGeomFactor:radar geom factor is not set").throwMe(); 
      }

    if (peak_radar_geom_factor_set_) {return (peak_radar_geom_factor_);}
    peak_radar_geom_factor_ = radar_geom_factor[0][0][Ncenter_patch_];

    for (unsigned int i_patch=0; i_patch < Npatch_;++i_patch)
    for (unsigned int i = 0; i < Nrange_;++i)
    for (unsigned int j = 0; j < Ndoppler_;++j)
      {
	if (peak_radar_geom_factor_< radar_geom_factor[i][j][i_patch])
	  peak_radar_geom_factor_ = radar_geom_factor[i][j][i_patch];
	if (peak_radar_geom_factor_< radar_geom_factor_mirror[i][j][i_patch])
	  peak_radar_geom_factor_ = radar_geom_factor_mirror[i][j][i_patch];
      }
    return(peak_radar_geom_factor_);
  }


//-------------------------------------------
//Calculate usable area for a given prf and pulse gate
//-------------------------------------------
Uvar SARAmb::calUsableArea(const Uvar& X0, 
				 const Uvar& Pn,     
				 const Uvar& x_res,
				 const Uvar& rg_res)
				
{
    if (!radar_geom_factor_set_)
      {
	ErrorMessage e("SARAmb::calUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
 
    for(unsigned int i=0;i<Nrange_;++i)
      for(unsigned int j=0;j<Ndoppler_;++j){
	good_bin[i][j] = 0;//reset good bin indicator
	thermal_snr_SL[i][j]= 0.0;
	amb_ratio_SL[i][j] = 0.0;
      }
   
    double d_X0 = X0.getInUnits("kg km km/(s s s)");
    double d_area_res = (x_res*rg_res).getInUnits("km km");
    double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
    usable_area_SL = Uvar(0,"km km");
    //---------------
    //examine data values falling inside 
    //range_start/range_end, dop_start/dop_end
    //-----------------  

   
    for (unsigned int i_range = 0 ; i_range <Nrange_;++i_range)
    for (unsigned int j_dop =0 ; j_dop < Ndoppler_ ;++j_dop)
      {
	double amb_radar_geom =0.0;
	for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
	  {
	    amb_radar_geom += radar_geom_factor[i_range][j_dop][i_patch];
	    if (mirror_geom_set_)
	      {
		amb_radar_geom += radar_geom_factor_mirror[i_range][j_dop][i_patch];
	      }
	  }
	amb_radar_geom -= radar_geom_factor[i_range][j_dop][Ncenter_patch_];
	
	//------------------------------------------
	//calculate amb ratio and noise equivalent sigma0
	//-----------------------------------------
	double ambratiodB = 10.0 
	  * log(radar_geom_factor[i_range][j_dop][Ncenter_patch_]
		/amb_radar_geom)/log(10.0);
	double  area_ratio= d_area_res/area[i_range][j_dop];
	echo_power[i_range][j_dop] = d_X0
	  *radar_geom_factor[i_range][j_dop][Ncenter_patch_]
	  *area_ratio
	  /muhleman_backscatter(muhleman_k1_,
				muhleman_k2_,
				Uvar(incidenceangle[i_range][j_dop],"rad"));
	double signal_in_Watt = echo_power[i_range][j_dop];
	double nesigma0dB=10.0
	  * log(noise_in_Watt/signal_in_Watt)/log(10.0);

	//------
	//Store thermal snr and amb ratio 
	//------
	thermal_snr_SL[i_range][j_dop] = -nesigma0dB;
	amb_ratio_SL[i_range][j_dop] = ambratiodB;
	
	//------------------------------------
	//IN case there is no solution (surface intercept point)
	//-----------------------------------
	if (radar_geom_factor[i_range][j_dop][Ncenter_patch_]==0.0)
	  {
	    amb_ratio_SL[i_range][j_dop] = -100.0;
	    thermal_snr_SL[i_range][j_dop] = -100.0;
	  }
 
	//-------------
	//Pick usable area satisfying all the following conditions
	//-------------
	if (amb_ratio_SL[i_range][j_dop] >=signal_to_ambdB_ &&
	    thermal_snr_SL[i_range][j_dop] >= -noise_equivalent_sigma0dB_ &&
	    oneway_beam_gaindB[i_range][j_dop] >= min_gaindB_)
	  { 
	    usable_area_SL+=Uvar(area[i_range][j_dop],"km km");
	    good_bin[i_range][j_dop] = 1;
	  }
      }//loop over i_range, j_dop
    usable_area_set_ = true;
    //cout<<"single look area "<< usable_area_SL;
    return(usable_area_SL);
  }


//---------------------------------------
//Calculate usable area for multilookings
//--------------------------------------
Uvar SARAmb::calMultilookUsableArea(const Uvar& X0, 
					  const Uvar& Pn, 
					  const Uvar& x_res,
					  const Uvar& rg_res,
					  const Uvar& bpd,
					  const FloatVector& velocity,
					  const unsigned int& N_p)
  {
  if (!radar_geom_factor_set_)
      {
	ErrorMessage e("SARAmb::calMultilookUsableArea: radar geom factor should be calculated first");
	e.throwMe();
      }
  
  
  double d_X0 = X0.getInUnits("kg km km/(s s s)");
  double d_area_res = (x_res * rg_res).getInUnits("km km");
  double noise_in_Watt = Pn.getInUnits("kg km km/(s s s)");
  usable_area_ML=Uvar(0.0,"km km");
  double dop_interval = doppler[1] - doppler[0];
  unsigned int j_dop_center=DopplerIndex(center_doppler_);
  
  
  //--------------------------------
  //Now calculate number of looks for data set
  //encloded by (i_range_start_ML,j_dop_start_ML)
  //and (i_range_end_ML, j_dop_end_ML)
  //--------------------------------
  double speed = velocity.magnitude().getInUnits("km/s");
  double d_bpd = bpd.getInUnits("s");
  double cos_theta = (center_doppler_ * lambda_ /(2.0 * speed));

  if (cos_theta > 1.0)
    {
      ErrorMessage e("SARAmb::cos theta value is larger than 1.0");
      e.throwMe();
    }
  double ML_freq_bp = 2.0 * speed
                   * speed * d_bpd/(lambda_*center_range_);
  ML_freq_bp *= (1.0-cos_theta * cos_theta);

  //  double  ML_freq_on = 1/(pri_ * double(N_p));

  //frequency resolution=1/(integration time(=pri * N_p))
 
  //cout<<"frequency shift "<<ML_freq_bp<<endl;
  //-------------- SH's simplied model produces the same result----------
  //double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  //Uvar ML_freq_bp = 
  //2.0 * velocity.magnitude()* Vst_bore * bpd*sin_theta/(lambda*range_bore);
  //note: Vst_bore = V * (1-cos(theta)^2)^0.5 --> same result!!!!!
  //for detailed information, see c:/ygim/radar_work/iso_doppler_profile
  //-------------------------------------------------------------------------

  //unsigned int ML_i_fon = 
  //(unsigned int) (ML_freq_on / dop_interval).getInUnits("");
  
  unsigned int ML_i_fbp = 
  (unsigned int) round_double( ML_freq_bp / dop_interval);
  
  
  if (ML_i_fbp ==0)
    {
      ErrorMessage("amb_param_scan:Need longer burst period to produce a shift in dop cell").throwMe();//infinite number of looks, unrealistic solution
      ML_i_fbp=1;
    }

  //else if(ML_i_fbp > Ndoppler_)
  //{
  //  nlook_ = 1;
  //  //doppler frequency change after each burst is greater 
  //  //than process bandwidth
  //}
  //else
  //{
  //  nlook_=(unsigned int) Ndoppler_/ML_i_fbp;
  //}

  //calculate number of looks
  //usable process window > min_gain
  vector<double> doppler_container;
  doppler_container.clear();
  for (unsigned int i_range = 0; i_range <Nrange_;++i_range)   
  for (unsigned int j_dop = 0; j_dop <Ndoppler_;++j_dop)
    {
      if(oneway_beam_gaindB[i_range][j_dop] >= min_gaindB_)
	{
	  doppler_container.push_back(doppler[j_dop]);
	}
    }

  if (doppler_container.size() == 0)
    {
      nlook_ = 1;
    }
  else if (doppler_container.size()==1)
    {
      nlook_ = 1;
    }
  else
    {
      sort(doppler_container.begin(),doppler_container.end());
      double min_dop, max_dop;
      min_dop = doppler_container.front();
      max_dop = doppler_container.back();
      nlook_ = (unsigned int) round_double((max_dop - min_dop)/ ML_freq_bp);
    }
 
  //------
  //setting up mask
  //----
  unsigned int* mask= new unsigned int[Ndoppler_];
  for (unsigned int j_dop = 0; j_dop < Ndoppler_;++j_dop)
    {
      mask[j_dop] = j_dop % ML_i_fbp;	
    }  
  
  for(unsigned int i=0;i<Nrange_;++i)
    for(unsigned int j=0;j<Ndoppler_;++j){
      ML_good_bin[i][j] = 0;//reset good bin indicator
      thermal_snr_ML[i][j]= 0.0;
      amb_ratio_ML[i][j] = 0.0;
    }
  vector<double> crosstrack_extent;
  crosstrack_extent.clear();
  double  crosstrack_length;

  //------------------------------------------
  //Assume we process 3dB or 5dB points inside
  //------------------------------------------
  for (unsigned int i_range = 0; i_range <Nrange_;++i_range)
    {
      //-----------------------------
      //Here, I would like to calculate min signal_to_amb ration
      //after multilook
      //-----------------------------
      double* ML_signal_bin=new double[ML_i_fbp];
      double* ML_amb_bin=new double[ML_i_fbp]; 
      double* ML_every_signal_bin=new double[ML_i_fbp];
      double* ML_every_amb_bin= new double[ML_i_fbp]; 
       

      for(unsigned int ii=0;ii<ML_i_fbp;++ii){
	ML_signal_bin[ii]=0.0;
	ML_amb_bin[ii]=0.0;
	ML_every_signal_bin[ii]=0.0;
	ML_every_amb_bin[ii]=0.0;
      }     
      for (unsigned int j_dop = 0; j_dop <Ndoppler_;++j_dop){
	double amb_radar_geom=0.0;
	double signal =0.0;
	for (unsigned int i_patch = 0; i_patch < Npatch_;++i_patch)
	  {
	    if (i_patch != Ncenter_patch_) 
	      {
		amb_radar_geom += radar_geom_factor[i_range][j_dop][i_patch];
	      }
	    else
	      {
		signal = radar_geom_factor[i_range][j_dop][Ncenter_patch_];
	      }
	    if (mirror_geom_set_)
	      {
		amb_radar_geom += radar_geom_factor_mirror[i_range][j_dop][i_patch];
	      }
	  }
	
	ML_every_signal_bin[mask[j_dop]] += signal;
	ML_every_amb_bin[mask[j_dop]]+= amb_radar_geom;
	
	//if gain is larger than min gain
	//need special treatment
	if(oneway_beam_gaindB[i_range][j_dop] >= min_gaindB_)
	  {
	    ML_signal_bin[mask[j_dop]]+= signal;
	    ML_amb_bin[mask[j_dop]]+=amb_radar_geom;
	  }
      }//loop over j_dop
      
      //---------------------------------------
      //find multilook ambiguity ratio inside min beam gain
      //---------------------------------------
      vector<double> amb_inside_window;
      amb_inside_window.clear();
      for (unsigned int i_fbp = 0; i_fbp < ML_i_fbp;++i_fbp)
	{
	  if(ML_signal_bin[i_fbp]!= 0.0)
	    {
	    amb_inside_window.push_back(10.0 * log(ML_signal_bin[i_fbp]/ML_amb_bin[i_fbp])/log(10.0));
	    }
	}

      //------------------
      //find lowest ambiguity ratio inside min beam gain
      //------------------
      double  min_amb_ratio_inside_min_gain ;      
      if(amb_inside_window.size()==0) min_amb_ratio_inside_min_gain = 0;
      else
	{
	  sort(amb_inside_window.begin(),amb_inside_window.end());
	  min_amb_ratio_inside_min_gain=amb_inside_window.front();
	}

     

      //------------------------
      //Assign amb values
      //----------------------------
      //find min amb ratio inside min beam gain
      for (unsigned int j_dop =0; j_dop < Ndoppler_;++j_dop)
       {
	 double area_ratio=d_area_res/area[i_range][j_dop];
	 double signal_in_Watt = d_X0
	   *radar_geom_factor[i_range][j_dop][Ncenter_patch_]
	   *area_ratio
	   /muhleman_backscatter(muhleman_k1_,
				 muhleman_k2_,
				 Uvar(incidenceangle[i_range][j_dop],"rad"));
	   
	   double nesigma0dB=10.0
	   * log(noise_in_Watt/signal_in_Watt/sqrt(double(nlook_)))/log(10.0);
	
	 //------
	 //Store thermal snr and amb ratio 
	 //------
	 thermal_snr_ML[i_range][j_dop] = -nesigma0dB;	 
	 //-----------------------------------
	 //Set amb ratio
	 //----------------------------------
	 if(oneway_beam_gaindB[i_range][j_dop] >= min_gaindB_)
	   {
	     amb_ratio_ML[i_range][j_dop]= min_amb_ratio_inside_min_gain;
	   }
	 else
	   {
	     amb_ratio_ML[i_range][j_dop]= 10.0 
	       * log(ML_every_signal_bin[mask[j_dop]]/ML_every_amb_bin[mask[j_dop]]/log(10.0));
	   }

 
	 //-------------
	 //Pick usable area satisfying all the following conditions
	 //-------------
	 if (amb_ratio_ML[i_range][j_dop]>=signal_to_ambdB_ &&
	     thermal_snr_ML[i_range][j_dop] >= -noise_equivalent_sigma0dB_ &&
	     oneway_beam_gaindB[i_range][j_dop]>= min_gaindB_)
	   { 
	     usable_area_ML+=Uvar(area[i_range][j_dop],"km km");
	     ML_good_bin[i_range][j_dop] = 1;
	     crosstrack_extent.push_back(crosstrack[i_range][j_dop]);
	     crosstrack_extent.push_back(crosstrack[i_range+1][j_dop]);
	     if(j_dop==j_dop_center)
	       {
		 double a=crosstrack[i_range+1][j_dop]-crosstrack[i_range][j_dop];
		 crosstrack_length +=fabs(a);
		 //always take positve value	     
	       }
	   }
       }//loop over j_dop
      delete[] ML_signal_bin;
      delete[] ML_amb_bin;
      delete[] ML_every_signal_bin;
      delete[] ML_every_amb_bin;
      //delete[]  amb_ratio_inside_window;
    }//loop over i_range
  

  //------------------------------
  //Let's extract crosstrack min and max
  //-------------------------------
  ML_crosstrack_extent_min_ = 0;
  ML_crosstrack_extent_max_ =0 ;
  ML_crosstrack_length_=crosstrack_length;
  if (crosstrack_extent.size() > 0)
    {//only if ML usable area is not zero
      sort(crosstrack_extent.begin(),crosstrack_extent.end());//sorting
      ML_crosstrack_extent_min_ = crosstrack_extent.front();//first element
      ML_crosstrack_extent_max_ = crosstrack_extent.back();//last element
    }
  ML_usable_area_set_ = true;
  //cout<<"ML look area "<< usable_area_ML;


  delete[] mask;

  return(usable_area_ML);

  }


//--------------------------------
//return mirror_geom_set_
//---------------------------------
bool SARAmb::mirrorGeomSet()
  {
  return(mirror_geom_set_);
  }

//--------------------
//is good bin (works for single look)
//---------------------
bool SARAmb::isGoodBin(const Uvar& range_input, const Uvar& dop_input,
		      const  bool& use_multilook_ambiguity)
  {
    double r = range_input.getInUnits("km");
    double f= dop_input.getInUnits("Hz");
    return(isGoodBin(r,f,use_multilook_ambiguity));
  }
 
bool SARAmb::isGoodBin(const double& range_in_km, const double& doppler_in_Hz,
		       const bool& use_ML)
{
  if((range_in_km<range[0] || range_in_km>range[Nrange_-1])|| 
     (doppler_in_Hz<doppler[0] || doppler_in_Hz>doppler[Ndoppler_-1])){
    return(false);
    }
  unsigned int range_index = RangeIndex(range_in_km);
  unsigned int doppler_index=DopplerIndex(doppler_in_Hz);
  bool good = false;
  if(use_ML){
    if( ML_good_bin[range_index][doppler_index] == 1) good = true;
  }
  else if( good_bin[range_index][doppler_index] == 1) good = true;
  return(good);
}

double SARAmb::getAmbigRatio(const double& range_in_km, const double& doppler_in_hz, const bool& use_ML){
  // assume points outside of grid have an ambig to main ratio of 100
  // this is primarily a way to falg such points as VERY BAD
  if((range_in_km<range[0] || range_in_km>range[Nrange_-1])|| 
     (doppler_in_hz<doppler[0] || doppler_in_hz>doppler[Ndoppler_-1])){
    return(100.0);
  }

  unsigned int range_index = RangeIndex(range_in_km);
  unsigned int doppler_index=DopplerIndex(doppler_in_hz);
  if(use_ML){
    return(pow(10,-0.1*amb_ratio_ML[range_index][doppler_index]));
  }
  else{
    return(pow(10,-0.1*amb_ratio_SL[range_index][doppler_index]));
  }
}
void SARAmb::pointReport(ofstream& afs, float range_in_km,
			  float doppler_in_hz, bool use_ML)
{ 
  afs << "at doppler=" << doppler_in_hz << " Hz and range="
			 << range_in_km  << " km." << endl;

  if((range_in_km<range[0] || range_in_km>range[Nrange_-1])|| 
     (doppler_in_hz<doppler[0] || doppler_in_hz>doppler[Ndoppler_-1])){
    afs << "Outside of SARAmb grid" << endl;
    return;
  }

  afs << "SARAmb parameters are:" << endl;
  unsigned int range_index = RangeIndex(range_in_km);
  unsigned int doppler_index=DopplerIndex(doppler_in_hz);
  afs << "Grid size (range,doppler)=("<<Nrange_<<","<<Ndoppler_<<")" << endl;
  afs << "range_index= " << range_index << " doppler_index= " << doppler_index
      << endl;

  afs << "oneway_beam_gaindB=" 
      << oneway_beam_gaindB[range_index][doppler_index]
      << " dB Usable is >= " << min_gaindB_ << endl;
  if(use_ML){
    afs << "ML_good_bin=" << ML_good_bin[range_index][doppler_index] << endl;

    afs << "amb_ratio_ML= " << amb_ratio_ML[range_index][doppler_index]
	<< " dB Usable is >= " << signal_to_ambdB_ << endl;
    afs << "thermal_snr_ML= " << thermal_snr_ML[range_index][doppler_index]
	<< " dB Usable is >= " << -noise_equivalent_sigma0dB_ << endl;
    
  }
  else{
    afs << "SL_good_bin=" << good_bin[range_index][doppler_index] << endl;

    afs << "amb_ratio_SL= " << amb_ratio_SL[range_index][doppler_index]
	<< " dB Usable is >= " << signal_to_ambdB_ << endl;
    afs << "thermal_snr_SL= " << thermal_snr_SL[range_index][doppler_index]
	<< " dB Usable is >= " << -noise_equivalent_sigma0dB_ << endl;
    
  }

    afs << "SARAmb thinks the beam number is " << beam_id_+1 << endl;

}
//---------------------
//clear method
//---------------------
void SARAmb::clear()
{
  time_set_=false;
  state_set_=false;
  beam_number_set_=false;
  process_window_set_=false;

  radar_geom_factor_set_=false;
  peak_radar_geom_factor_set_=false;
  usable_area_set_=false;
  ML_usable_area_set_=false;
}



void SARAmb::createArrays()
 {
  
  
   //-----------------------------------------
   //range and doppler grid  variables
   //---------------------------------------------
   surface_intercept_grid.resize(Nrange_grid_*Ndoppler_grid_,Npatch_); 
   surface_intercept_mirror_grid.resize(Nrange_grid_*Ndoppler_grid_,Npatch_);   
   range= new double[Nrange_];
   doppler= new double[Ndoppler_];
   oneway_beam_gaindB= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   incidenceangle= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   
   crosstrack=  (double**)make_array(sizeof(double),2,Nrange_grid_,Ndoppler_grid_);
   alongtrack=  (double**)make_array(sizeof(double),2,Nrange_grid_,Ndoppler_grid_);
   area= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   radar_geom_factor= (double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   echo_power= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   radar_geom_factor_mirror= (double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   good_bin=(int**)make_array(sizeof(int),2,Nrange_,Ndoppler_);
   ML_good_bin= (int**)make_array(sizeof(int),2,Nrange_,Ndoppler_);
   thermal_snr_SL= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   amb_ratio_SL= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   thermal_snr_ML= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   amb_ratio_ML= (double**)make_array(sizeof(double),2,Nrange_,Ndoppler_);
   
   
   
   range_grid=(double**)make_array(sizeof(double),2,Nrange_grid_,Npatch_);
   doppler_grid=(double**)make_array(sizeof(double),2,Ndoppler_grid_,Npatch_);
   no_of_solution_grid=(int**)make_array(sizeof(int),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   oneway_beam_gain_grid=(double**)make_array(sizeof(double),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   incidenceangle_grid=(double**)make_array(sizeof(double),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   no_of_solution_mirror_grid=(int**)make_array(sizeof(int),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   oneway_beam_gain_mirror_grid=(double**)make_array(sizeof(double),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   incidenceangle_mirror_grid=(double**)make_array(sizeof(double),2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   
   
   //-----------------------------
   //Pixel variables
   //-----------------------------
   range_pixel=(double**) make_array(sizeof(double),2,Nrange_,Npatch_);
   doppler_pixel=(double**) make_array(sizeof(double),2,Ndoppler_,Npatch_);
   oneway_beam_gain_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   area_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   incidenceangle_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   oneway_beam_gain_mirror_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   area_mirror_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);
   incidenceangle_mirror_pixel=(double***)make_array(sizeof(double),3,Nrange_,Ndoppler_,Npatch_);      
   no_of_solution_pixel=(int***)make_array(sizeof(int),3,Nrange_,Ndoppler_,Npatch_);
   no_of_solution_mirror_pixel=(int***)make_array(sizeof(int),3,Nrange_,Ndoppler_,Npatch_);   
   grid_set_ = true;
 }


void SARAmb::deleteArrays()
 {
   delete[] range;
   delete[] doppler;
   
   free_array((void**)oneway_beam_gaindB,2,Nrange_,Ndoppler_);
   free_array((void**) incidenceangle,2,Nrange_,Ndoppler_);
   free_array((void**)crosstrack,2,Nrange_grid_,Ndoppler_grid_);
   free_array((void**)alongtrack,2,Nrange_grid_,Ndoppler_grid_);
   free_array((void**)area,2,Nrange_,Ndoppler_);
   free_array((void***)radar_geom_factor,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void**)echo_power,2,Nrange_,Ndoppler_);
   free_array((void***)radar_geom_factor_mirror,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void**)good_bin,2,Nrange_,Ndoppler_);
   free_array((void**)ML_good_bin,2,Nrange_,Ndoppler_);
   free_array((void**)thermal_snr_SL,2,Nrange_,Ndoppler_);
   free_array((void**)amb_ratio_SL,2,Nrange_,Ndoppler_);
   free_array((void**)thermal_snr_ML,2,Nrange_,Ndoppler_);
   free_array((void**)amb_ratio_ML,2,Nrange_,Ndoppler_);
   
   //-----------------------
   //grid variables(local variables)
   //----------------------
   free_array((void**) range_grid,2,Nrange_grid_,Npatch_);
   free_array((void**) doppler_grid,2,Ndoppler_grid_,Npatch_);
   free_array((void**) no_of_solution_grid ,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   free_array((void**) oneway_beam_gain_grid,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   free_array((void**) incidenceangle_grid,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   free_array((void**) no_of_solution_mirror_grid,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   free_array((void**) oneway_beam_gain_mirror_grid,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   free_array((void**)  incidenceangle_mirror_grid,2,Nrange_grid_*Ndoppler_grid_,Npatch_);
   
   //-----------------------------
   //Pixel variables
   //-----------------------------
   free_array((void**) range_pixel,2,Nrange_,Npatch_);
   free_array((void**) doppler_pixel,2,Ndoppler_,Npatch_);
   free_array((void***) oneway_beam_gain_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) area_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) incidenceangle_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) oneway_beam_gain_mirror_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) area_mirror_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) incidenceangle_mirror_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) no_of_solution_pixel,3,Nrange_,Ndoppler_,Npatch_);
   free_array((void***) no_of_solution_mirror_pixel,3,Nrange_,Ndoppler_,Npatch_);
 }


///////////////////////// DO NOT USE ANY MORE /////////////////////

//-------------------------------------------------
//Earlier version of ambiguity calculation tools
//All of the functions are implemented in Ambiguity Class Object
//---------------------------------------------------



void cal_ambiguity(const string& target_name,
		   Beam& beam,
		   const Frame& ftitan,
		   const StateVector& sc_state,
		   const Frame& track_frame,
		   const Uvar& prf,
		   const Uvar& pulse_gate,
		   const Uvar& pbw,
		   const Uvar& lambda,
		   const Uvar& X0,
		   const Uvar& Pn,
		   const Uvar& x_res,
		   const Uvar& rg_res,
		   const double& Ni,
		   const double&  muhleman_k1,
		   const double&  muhleman_k2,
		   const unsigned int& Nrange_patch,
		   const unsigned int& Ndop_patch,
		   const unsigned int& Nrange_bin,
		   const unsigned int& Ndop_bin,
		   const Uvar& range_offset,
		   const Uvar& frequency_offset,
		   ofstream& outputfile) 
                  throw(Unit::UnitError,ErrorMessage)
  { 
    
    typedef Array2D<DirectionVector> Dirmat;    
    Uvar usable_area = Uvar(0);
    Time t = sc_state.time();
    unsigned int i_range_center = int (Nrange_patch/2);
    unsigned int j_dop_center = int(Ndop_patch/2);
        

    //get pri
    Uvar pri = 1.0/prf;

    //get max beam gain
    //double max_gain = beam.getMaxGain();
    Uvar azimuth_width = beam.getAzimuthWidthOneWay();
    Uvar elevation_width=beam.getElevationWidthOneWay();
    //------------------------------------------
    //First use Frame method to get range and dop
    //------------------------------------------
    unsigned int beamnum = beam.getBeamNumber();
    Frame fbeam("CASSINI_RADAR_" + toStr(beamnum),"Cassini");
    DirectionVector boresight("boresight",fbeam,t,0,0,1);
    Uvar range_F, dop_F, thetai_F,altitude_F;
    TargetGeom tg(t);
    tg.setState(sc_state);
    tg.setLookDirection(boresight);
    tg.setTarget(target_name,ftitan);
    Uvar radius = tg.radius();
    range_F = tg.range();
    dop_F = tg.doppler(lambda); 
    thetai_F = tg.incidenceAngle();
    altitude_F = tg.altitude();
    
    //-----------------------------------------
    //Then, use non-frame method(vector equations)
    // to calculate boresight range and dop
    //------------------------------------------
    
    PositionVector sc_pos = sc_state.position();
    sc_pos.representIn(ftitan);
    DirectionVector pos_dir = sc_pos;
    pos_dir.representIn(fbeam);   
    
    FloatVector velocity = sc_state.velocity();
    velocity.representIn(fbeam);
    DirectionVector velocity_direction = velocity;
   
    //---------
    //doppler
    //---------
    double v_b = dot(velocity_direction,boresight);
    Uvar dop0 = 2.0* v_b * velocity.magnitude()/lambda;
   
    //---------
    //range
    //---------
    double p_b = dot(pos_dir,boresight);
    Uvar in_sqrt = sc_pos.magnitude() *sc_pos.magnitude() * p_b * p_b;
    in_sqrt += radius * radius - sc_pos.magnitude()*sc_pos.magnitude();
    Uvar range0 = Uvar(0.0);
    if(in_sqrt.getValue() > 0.0) 
      {
	range0 = -sc_pos.magnitude() * p_b - sqrt(in_sqrt);
      }
    else
      {
	throw ErrorMessage("No intercepting point in the boresight direction");
	
      }  
   
  
    //-----------------------------------------
    //compare frame-method with non frame method
    //------------------------------------------
    double diff_range = range0.getInUnits("km") - range_F.getInUnits("km");
    double diff_dop   = dop0.getInUnits("Hz") - dop_F.getInUnits("Hz");
    
    if (fabs(diff_range)>0.001)
      {
	throw ErrorMessage("Frame method does not produce the same range");
      }
    if (fabs(diff_dop) > 0.001)
      {
	throw ErrorMessage("Frame method does not produce the same dop");
      }
    

    //----------------------------------------------------------
    //Determine range and dop centers for each patch
    //----------------------------------------------------------      
    Uvec range_center("range_center",Nrange_patch);
    Uvec dop_center("dop_center",Ndop_patch);
    Umat range_bin("range_bin",Nrange_patch,Nrange_bin);
    Umat dop_bin("dop_bin",Ndop_patch,Ndop_bin);
    range_center = Uvar(0.0);
    dop_center = Uvar(0.0);
    range_bin = Uvar(0.0);
    dop_bin = Uvar(0.0);
   

    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      range_center(i) = range0;
      range_center(i) += (double(i) -double(i_range_center))*speed_light*pri/2.0;
      range_center(i) -= range_offset;
      for (unsigned int i_range=0; i_range< Nrange_bin; i_range++)
      { 
	range_bin(i,i_range) = range_center(i) -0.5*pulse_gate;
	range_bin(i,i_range)+= double(i_range)/double(Nrange_bin-1)*pulse_gate;
      }
    }

    for (unsigned int j = 0; j < Ndop_patch; j++)
    {
      dop_center(j) = dop0 + (double(j)- double(j_dop_center))*prf;
      dop_center(j) -= frequency_offset;
      for (unsigned int j_dop=0; j_dop<Ndop_bin; j_dop++)
      {
	dop_bin(j,j_dop) = dop_center(j) -0.5 * pbw ;
	dop_bin(j,j_dop) += double(j_dop)/double(Ndop_bin-1)*pbw;
      }
    }

    //--------------------------------------------------------
    //Display range and dop settings on the screen
    //-------------------------------------------------------  
    /*
    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      cout<<"range center "<<range_center(i)<<endl ;
      if (i == i_range_center)
      {
	for (unsigned int i_range=0; i_range< Nrange_bin; i_range++)
	{ 
	  //cout<<"range bin "<<range_bin(i,i_range) <<endl;
	}
      }
    }

    for (unsigned int i = 0; i < Ndop_patch; i++)
    {
      cout<<"dop center "<<dop_center(i)<<endl;
      if (i == j_dop_center)
      {
	for (unsigned int i_dop=0; i_dop<Ndop_bin; i_dop++)
	{
	  //cout<<"dop_bin "<<dop_bin(i,i_dop)<<endl;
	}
      }
    }
    */



    Imat no_of_solution("no_solution",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat look_azimuth("look_azimuth_angle",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat look_elevation("look_elevation_angle",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  lat("bin_patch_lat",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  lon("bin_patch_lon",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  range("range",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  doppler("doppler",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  thetai("thetai",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  alongtrack("alongtrack",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat crosstrack("crosstrack",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat bin_power_main("bin_power_main",Nrange_bin,Ndop_bin);
    Umat bin_power_amb("bin_power_amb",Nrange_bin,Ndop_bin);
    Umat normalized_bin_power("normalized_bin_power",Nrange_bin,Ndop_bin);
    Dmat  bin_backscatter("bin_backscatter",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dmat bin_antenna("bin_antenna",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dmat center_patch_antenna("center_bin_antenna",Nrange_bin,Ndop_bin);
    Umat patch_sum("patch_strength",Nrange_bin,Ndop_bin);
    Umat total_noise("total_nose_summed",Nrange_bin,Ndop_bin);
    Umat center_area("center_area",Nrange_bin,Ndop_bin);
    Umat area("area_element",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dirmat surface_intercept("surface_intercept",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);

   

    //Initialize variables
    no_of_solution= 0;
    look_azimuth=Uvar(0.0);
    look_elevation= Uvar(0.0);
    lat= Uvar(0.0);
    lon= Uvar(0.0);
    range= Uvar(0.0);
    doppler= Uvar(0.0);
    thetai= Uvar(0.0);
    alongtrack=Uvar(0.0);
    crosstrack=Uvar(0.0);
    bin_power_main= Uvar(0.0);
    bin_power_amb= Uvar(0.0);
    normalized_bin_power= Uvar(0.0);
    bin_backscatter= 0.0;
    bin_antenna= 0.0;
    center_patch_antenna= 0.0;
    patch_sum= Uvar(0.0);
    total_noise= Pn; // Pn = k * Tsys * Bn
    center_area= Uvar(0.0); 
    area= Uvar(0.0);

    //------------------------------------------------------------------
    //Find  intercepting points corresponding range and dop values
    //no_interception: store how many range and dop bins have a solution
    //------------------------------------------------------------------
  
    TargetGeom tg1(t);
    tg1.setTrackFrame(track_frame);
    for (unsigned int i = 0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
	{
	unsigned int patch_index = i * Ndop_patch+j;   
	for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
	for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
	  {
	    unsigned int bin_index = i_range * Ndop_bin+j_dop;	
	    tg1.setState(sc_state);
	    tg1.setTarget(target_name,ftitan);
	    tg1.setRangeDopplerInBeam3Frame(range_bin(i,i_range), dop_bin(j,j_dop),lambda);
	    DirectionVector look = tg1.lookDirection();
	    if (tg1.foundSurfaceIntercept()==false)
	      {
	      no_of_solution(patch_index,bin_index) = 0;//no surface intercept point
	      }
	    else
	      {
	      no_of_solution(patch_index,bin_index) = 1;   
	      DirectionVector new_look= look.representIn(fbeam);
	      new_look.getAzimuthElevation(look_azimuth(patch_index,bin_index),
				     look_elevation(patch_index,bin_index));
	      lat(patch_index,bin_index)= tg1.lat();
	      lon(patch_index,bin_index) = tg1.lon();
	      range(patch_index,bin_index)= tg1.range();
	      doppler(patch_index,bin_index) =tg1.doppler(lambda);
	      thetai(patch_index,bin_index)= tg1.incidenceAngle(); 
	      tg1.interceptAlongCross(alongtrack(patch_index,bin_index),
				      crosstrack(patch_index,bin_index));
	      surface_intercept(patch_index,bin_index)=tg1.surfaceIntercept();
	      //DirectionVector mirror_look
	      } 
	    tg1.reset(t);
	  }//Loop over i_range and j_dop
	}//Loop over range_patch and dop_patch

    //    cout<<"Finish finding all surface intercept points"<<endl;
    //    cout<<"Start calculating area elements and radar return signal "<<endl;
    //------------------------------------------	
    //End of finding all range and dop solutions
    //if there is a surface intercept point, no_of_s
    //----------------------------------------------
    //****************************************
    //    Main Ambiguity Calculation
    //there is a solution for the given  range and dop value
    //range_bin(i,i_range) and dop_bin(j,j_dop)
    //let's try to get radar return signal strenth here
    //*******************************
       
    for (unsigned int i =0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
      {
      unsigned int patch_index = i * Ndop_patch+j;   
      for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
      for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
	{
	unsigned int bin_index = i_range * Ndop_bin+j_dop;	
	if(i_range < Nrange_bin -1 
	   && j_dop <  Ndop_bin -1 
	   && no_of_solution(patch_index,bin_index)==1)
	  { 
	  unsigned int bin_Up,  bin_Left,bin_Diag;
	  bin_Up = (i_range+1)*Ndop_bin+j_dop;	
	  bin_Left = i_range*Ndop_bin+(j_dop+1);
	  bin_Diag = (i_range+1)*Ndop_bin+j_dop+1;
	 
	  //------------------------------------
	  //find Area element covered by (range,dop)
	  //divide an area enclosed by two range and dop values
	  // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
	  //into two triangles and calclate their areas and add them
	  //together
	  //-------------------------------------
	  if (no_of_solution(patch_index,bin_Up) * 
	      no_of_solution(patch_index,bin_Diag)*
	      no_of_solution(patch_index,bin_Left) == 1)
	    {
	      //Here I will calculate surface area, beam , backscatter
	      //radar return.  Beware that I do not include 
	      //i_range = Nrang_bin-1 j_dop = Ndop_bin-1
	      // their contributions are counted when i_range= Nrange_bin-2
	      // and j_dop = Ndop_bin -2
	    Uvar side_a = radius * 
	      fabs((surface_intercept(patch_index,bin_index)
		    .angle(surface_intercept(patch_index,bin_Left))
		    .getInUnits("rad")));
	    Uvar side_b = radius * 
	      fabs((surface_intercept(patch_index,bin_index)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    Uvar side_c = radius * 
	      fabs((surface_intercept(patch_index,bin_Left)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    Uvar cos_c = (side_a*side_a + side_b*side_b - side_c*side_c)
	      /(2.0 * side_a*side_b);
	    Uvar angle_c = Uvar(acos(cos_c.getValue()),"rad");
	    
	    area(patch_index,bin_index) = 
	      0.5 * sin(angle_c) * side_a * side_b;
	    
	    Uvar side_d = radius * 
	      fabs((surface_intercept(patch_index,bin_Left)
		    .angle(surface_intercept(patch_index,bin_Diag))
		    .getInUnits("rad")));
	    Uvar side_e = radius * 
	      fabs((surface_intercept(patch_index,bin_Diag)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    cos_c = (side_d*side_d + side_e*side_e - side_c*side_c)
	      /(2.0 * side_d*side_e);
	    angle_c = Uvar(acos(cos_c.getValue()),"rad");
	    area(patch_index,bin_index) += 
	      0.5 * sin(angle_c) * side_d * side_e;	    
	    //--------
	    //averaged beam gain over bin_index, bin_Up,bin_Left,bin_Diag
	    //---------
	    
	    double beam_gain = beam.bilinear(
				      look_azimuth(patch_index,bin_index),
				    look_elevation(patch_index,bin_index));
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Up),
				    look_elevation(patch_index,bin_Up));		
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Left),
			            look_elevation(patch_index,bin_Left));
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Diag),
				    look_elevation(patch_index,bin_Diag));
	    beam_gain /= 4.0;
	        
	    //--------------------
	    //averaged back scattering coeff
	    //--------------------
	    double sigma0_elem = muhleman_k1 
	      * cos(thetai(patch_index,bin_index))
	      / pow(sin(thetai(patch_index,bin_index))
		    +muhleman_k2*cos(thetai(patch_index,bin_index))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Up))
	      / pow(sin(thetai(patch_index,bin_Up))
		    +muhleman_k2*cos(thetai(patch_index,bin_Up))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Left))
	      / pow(sin(thetai(patch_index,bin_Left))
		    +muhleman_k2*cos(thetai(patch_index,bin_Left))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Diag))
	      / pow(sin(thetai(patch_index,bin_Diag))
		    +muhleman_k2*cos(thetai(patch_index,bin_Diag))
		    ,3.0);
	    sigma0_elem /= 4.0;
	   
	    
	    //-----------------------------
	    //store beam gain and sigm0
	    //----------------------------
	    bin_antenna(patch_index,bin_index) = beam_gain;
	    bin_backscatter(patch_index,bin_index) = sigma0_elem;
	    //---------------------------
	    //average range
	    //---------------------------
	    Uvar avg_range =(range_bin(i,i_range) +range_bin(i,i_range+1))/2.0;    
	    		
	    //---------------------------------------
	    //find X = G^2 * lambda^2 * sigma0 / R^4
	    //--------------------------------------
	    Uvar X = beam_gain * beam_gain * lambda * lambda* sigma0_elem
	      /pow(avg_range,4);  	
	  
	    if (i == i_range_center && j == j_dop_center)
	      {
	      bin_power_main(i_range,j_dop)=X*X0*area(patch_index,bin_index);  
	      normalized_bin_power(i_range,j_dop)= bin_power_main(i_range,j_dop);
	      normalized_bin_power(i_range,j_dop) /=sigma0_elem;  
	      patch_sum(i,j) += X*X0*area(patch_index,bin_index);  
	      center_patch_antenna(i_range,j_dop)= beam_gain;
	      center_area(i_range,j_dop)=area(patch_index,bin_index); 
	      
	      //Debugging routine
	      /*
	      cout<<"beam number range and dop "<<beamnum<<" "
		  <<i_range<<" "<<j_dop<<endl;
	      cout<<look_azimuth(patch_index,bin_index).getInUnits("deg")<<" "
	       <<look_elevation(patch_index,bin_index).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Up).getInUnits("deg")<<
	       " "<<look_elevation(patch_index,bin_Up).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Left).getInUnits("deg")<<" "
		<<look_elevation(patch_index,bin_Left).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Diag).getInUnits("deg")<<" "
	       <<look_elevation(patch_index,bin_Diag).getInUnits("deg")<<endl;
	      cout<<"gain area sigma0 "<<10.0*log(beam_gain/max_gain)/log(10.0)
	       <<" " <<area(patch_index,bin_index)<<" "<<sigma0_elem<<endl<<endl;
	      */
	      // End of Debugging routine
	      }
	    else
	      {
	      bin_power_amb(i_range,j_dop) += X * X0 
		* area(patch_index,bin_index);
	      patch_sum(i,j) +=X * X0 
		* area(patch_index,bin_index);
	      }
	    }//if there is a solution at Up, Left, 
	     //Diag and (i_range != Nrange_bin -1) and (j_dop !=Ndob_bin-1)
	  }//if there is a solution   
	}//Loop over range/dop bin
      }//Loop over range/dop patch
      


   
      
    //cout<<"Starting counting usable area "<<endl;
   
    //-----------------------------------------------------------------
    //Determine whether each bin is good for imaging
    // criteria
    // (1) amb_ratio > 14 dB
    // (2) thermal noise equivalent sigma 0 < -10 dB
    // (3) total snr > thermal snr + 3db
    //-----------------------------------------------------------------
    Imat is_good_bin("is_this_a_good_bin",Nrange_bin,Ndop_bin);   
    Imat is_good_bin_amb("amb_ratio_requirements",Nrange_bin,Ndop_bin);   
    Imat is_good_bin_nesigma0("noise_equivalent_sigma0",Nrange_bin,Ndop_bin);  
    Imat is_good_bin_totalSNR("total_snr_to_theraml_snr",Nrange_bin,Ndop_bin);
    Imat is_good_bin_pn_to_amb("ration_of_pn_to_amb",Nrange_bin,Ndop_bin);
    is_good_bin = 0;
    is_good_bin_amb = 0;
    is_good_bin_nesigma0=0;
    is_good_bin_totalSNR = 0;
    is_good_bin_pn_to_amb = 0;

    for (unsigned int i_range = 0; i_range < Nrange_bin-1; i_range++)
    for (unsigned int j_dop = 0 ; j_dop < Ndop_bin-1; j_dop++)
    { 
      double bin_ratio;
      double bin_ratiodB;
      unsigned int patch_index = i_range_center * Ndop_patch+j_dop_center;
      unsigned int bin_index = i_range * Ndop_bin+j_dop;	

      if (no_of_solution(patch_index,bin_index) == 1)
	{
	  if (bin_power_main(i_range,j_dop)==Uvar(0.0))
	  {
	  bin_ratio = 1e-100;
	  }
	else
	  {
	  bin_ratio = 
	    bin_power_main(i_range,j_dop).getInUnits("kg km km/(s s s)");
	  bin_ratio /= 
	    bin_power_amb(i_range,j_dop).getInUnits("kg km km/(s s s)");	
	  }
	}
      else if (no_of_solution(patch_index,bin_index) == 0)
	{
	bin_ratio = 1e-100;	
	}
      else
	{
	throw ErrorMessage("no_solution data set has been contaminated");
	}
      bin_ratiodB = 10.0*  log(bin_ratio)/log(10.0);

      Uvar power_sigma0_normalized =  normalized_bin_power(i_range,j_dop) ;  
      power_sigma0_normalized *= x_res * rg_res /area(patch_index,bin_index);

      double ne0 = Pn.getInUnits("kg km km/(s s s)");
      ne0 /= power_sigma0_normalized.getInUnits("kg km km/(s s s)");
      double ne0_dB = 10.0 * log(ne0)/log(10.0);
      
      Uvar main_power_normalized= bin_power_main(i_range,j_dop);
      main_power_normalized *= x_res * rg_res /area(patch_index,bin_index);     
      double thermal_snr = main_power_normalized.getInUnits("kg km km/(s s s)")/
      	Pn.getInUnits("kg km km/(s s s)");
      double thermal_snrdB = 10.0 * log(thermal_snr)/log(10.0);

      Uvar amb_power_normalized = bin_power_amb(i_range,j_dop);
      amb_power_normalized *= x_res * rg_res /area(patch_index,bin_index);

      double total_snr = main_power_normalized.getInUnits("kg km km/(s s s)");
      total_snr /=(amb_power_normalized.getInUnits("kg km km/(s s s)")
		   +Pn.getInUnits("kg km km/(s s s)"));
      double total_snrdB = 10.0 * log(total_snr)/log(10.0);
      double pn_to_amb = Pn.getInUnits("kg km km/(s s s)")
                /amb_power_normalized.getInUnits("kg km km/(s s s)");

      if (bin_ratiodB >= 14.0) is_good_bin_amb(i_range,j_dop) = 1;
      if (ne0_dB <= -10.0)    is_good_bin_nesigma0(i_range,j_dop)=1;
      if ((total_snrdB + 3.0)>= thermal_snrdB) 
               is_good_bin_totalSNR(i_range,j_dop) = 1;
      if (pn_to_amb >= 1.0) is_good_bin_pn_to_amb(i_range,j_dop) = 1;

      if (is_good_bin_amb(i_range,j_dop)==1
	  && is_good_bin_nesigma0(i_range,j_dop)==1
	  )
      {
	usable_area += area(patch_index,bin_index);
	is_good_bin(i_range,j_dop) = 1;
      }
    }

    
    //---------------------------------------------------------
    //Generate range and dop range for each patch for outputfile
    //---------------------------------------------------------
    Uvec range_lower("range_lower",Nrange_patch);
    Uvec range_upper("range_upper",Nrange_patch);
    Uvec dop_lower("dop_lower",Ndop_patch);
    Uvec dop_upper("dop_upper",Ndop_patch);

 
    for (unsigned int i = 0 ; i <Nrange_patch; i++)
      {
      range_lower(i) = range_center(i) -0.5*pulse_gate;
      range_upper(i)=range_center(i)+0.5*pulse_gate;
      }
    for (unsigned int j = 0; j < Ndop_patch; j++)
      {      
      dop_lower(j) = dop_center(j) -0.5 * pbw ;
      dop_upper(j) = dop_center(j) + 0.5*pbw;      
      }

    //---------------------------------------------------------
    //Find average incident angle, antenna_gain, backscatter
    // azi
    //---------------------------------------------------------
    Umat avg_incidence("avg_incidence",Nrange_patch,Ndop_patch);
    Dmat avg_backscatter("avg_backsctter",Nrange_patch,Ndop_patch);
    Dmat avg_antenna("avg_antenna_gain",Nrange_patch,Ndop_patch);
    Umat avg_azi("avg_azi_angle",Nrange_patch,Ndop_patch);
    Umat avg_elev("avg_elev",Nrange_patch,Ndop_patch);
    Umat avg_area("avg_area",Nrange_patch,Ndop_patch);

    avg_incidence = Uvar(0.0);
    avg_backscatter = 0.0;
    avg_antenna = 0.0;
    avg_azi = Uvar(0.0);
    avg_elev = Uvar(0.0);
    avg_area = Uvar(0.0);

    for (unsigned int i =0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
      {
      unsigned int patch_index = i*Ndop_patch+j;
      unsigned int bin_count = 0;
      for (unsigned int i_range = 0; i_range<Nrange_bin-1; i_range++)
      for (unsigned int j_dop =0; j_dop<Ndop_bin-1; j_dop++)
	{
        unsigned int bin_index = i_range*Ndop_bin+j_dop;	
	if (no_of_solution(patch_index,bin_index) == 1)
	  {
	  avg_incidence(i,j) += thetai(patch_index,bin_index);
	  avg_backscatter(i,j) += bin_backscatter(patch_index,bin_index);
	  avg_antenna(i,j) +=bin_antenna(patch_index,bin_index);
	  avg_azi(i,j) += look_azimuth(patch_index,bin_index);
	  avg_elev(i,j) += look_elevation(patch_index,bin_index);
	  avg_area(i,j) += area(patch_index,bin_index);
	  bin_count ++;
	  }
	}//Loop over range/dop_bin
      if (bin_count != 0)
	{
	avg_incidence(i,j) /= double(bin_count);
	avg_backscatter(i,j) /=double(bin_count);
	avg_antenna(i,j)/=double(bin_count);
	avg_azi(i,j) /= double(bin_count);
	avg_elev(i,j) /=double(bin_count);
	avg_area(i,j) /= double(bin_count);
	}
      }//Loop over range/dop_patch
  


    //------------------------------------------
    //Output data file format
    //Let's write down the output to the file
    //-------------------------------------------   
    outputfile<<usable_area.getInUnits("km km")<<endl;

    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      outputfile<<range_lower(i).getInUnits("km")<<endl; 
      outputfile<<range_upper(i).getInUnits("km")<<endl;
    }

    for (unsigned int j = 0; j < Ndop_patch; j++)
    {      
      outputfile<<dop_lower(j).getInUnits("Hz")<<endl;
      outputfile<<dop_upper(j).getInUnits("Hz")<<endl;  
    }

    for (unsigned int i = 0; i<Nrange_patch;i++)
    for(unsigned int j = 0; j <Ndop_patch;j++)
      {
      outputfile<<avg_incidence(i,j).getInUnits("deg")<<endl;
      outputfile<<avg_backscatter(i,j)<<endl;
      outputfile<<avg_antenna(i,j)<<endl;
      outputfile<<avg_azi(i,j).getInUnits("deg")<<endl; 
      outputfile<<avg_elev(i,j).getInUnits("deg")<<endl;
      outputfile<<patch_sum(i,j).getInUnits("kg km km/(s s s)")<<endl;
      outputfile<<avg_area(i,j).getInUnits("km km")<<endl;
      }
    

    for (unsigned int i =0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
    for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
    for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
      {
      unsigned int patch_index = i*Ndop_patch+j;
      unsigned int bin_index = i_range*Ndop_bin+j_dop;	
      outputfile<<lat(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<lon(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<thetai(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<range(patch_index,bin_index).getInUnits("km")<<endl;
      outputfile<<doppler(patch_index,bin_index).getInUnits("Hz")<<endl;
      outputfile<<alongtrack(patch_index,bin_index).getInUnits("km")<<endl;
      outputfile<<crosstrack(patch_index,bin_index).getInUnits("km")<<endl;
      }
        
    for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
    for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
    {
      outputfile<<bin_power_main(i_range,j_dop).getInUnits("kg km km/(s s s)")<<endl;   
      outputfile<<bin_power_amb(i_range,j_dop).getInUnits("kg km km/(s s s)")<<endl;
      outputfile<<normalized_bin_power(i_range,j_dop).getInUnits("kg km km/(s s s)")<<endl;
      outputfile<<center_area(i_range,j_dop).getInUnits("km km")<<endl;
      outputfile<<center_patch_antenna(i_range,j_dop)<<endl;
      outputfile<<is_good_bin_amb(i_range,j_dop) <<endl;
      outputfile<<is_good_bin_nesigma0(i_range,j_dop)<<endl;
      outputfile<<is_good_bin_totalSNR(i_range,j_dop) <<endl;
      outputfile<<is_good_bin_pn_to_amb(i_range,j_dop)<<endl;
      outputfile<<is_good_bin(i_range,j_dop)<<endl;
    }
  }//END of cal_amb


//------------------------------------------------
//Calculate amb table for BS
//-----------------------------------------------
void cal_ambiguity_geometry_table(const string& target_name,
				  Beam& beam,
				  const Frame& ftitan,
				  const StateVector& sc_state,
				  const Uvar& lambda,
				  const double&  muhleman_k1,
				  const double&  muhleman_k2,
				  const Uvec& dop_freq_values,
				  const Uvec& range_values,
				  Dmat& cell_one_way_gain,
				  Dmat& cell_radar_geom_factor,
				  Umat& cell_cross_track,
				  Umat& cell_along_track,
				  Umat& cell_incidenceangle,
				  Umat& cell_area,
				  Umat& cell_lat,
				  Umat& cell_lon,
				  Dmat& cell_sigma0,
				  Uvec& range_axis,
				  Uvec& dop_axis)
  throw(Unit::UnitError,ErrorMessage)
  {
  typedef Array1D<DirectionVector> Dirvec;    
  
  Time t = sc_state.time();
      


  //get max beam gain
  //double max_gain = beam.getMaxGain();
  Uvar azimuth_width = beam.getAzimuthWidthOneWay();
  Uvar elevation_width=beam.getElevationWidthOneWay();
  
  //------------------------------------------
  //First use Frame method to get range and dop
  //------------------------------------------
  unsigned int beamnum = beam.getBeamNumber();
  Frame fbeam("CASSINI_RADAR_" + toStr(beamnum),"Cassini");
  DirectionVector boresight("boresight",fbeam,t,0,0,1);
  Uvar range_F, dop_F, thetai_F,altitude_F;
  TargetGeom tg(t);
  tg.setState(sc_state);
  tg.setLookDirection(boresight);
  tg.setTarget(target_name,ftitan);
  Uvar radius = tg.radius();
  range_F = tg.range();
  dop_F = tg.doppler(lambda); 
  thetai_F = tg.incidenceAngle();
  altitude_F = tg.altitude();
  
  //-----------------------------------------
  //Then, use non-frame method(vector equations)
  // to calculate boresight range and dop
  //------------------------------------------  
  PositionVector sc_pos = sc_state.position();
  sc_pos.representIn(ftitan);
  DirectionVector pos_dir = sc_pos;
  pos_dir.representIn(fbeam);   
  
  FloatVector velocity = sc_state.velocity();
  velocity.representIn(fbeam);
  DirectionVector velocity_direction = velocity;
  
  //---------
  //doppler
  //---------
  double v_b = dot(velocity_direction,boresight);
  Uvar dop0 = 2.0* v_b * velocity.magnitude()/lambda;
  
  //---------
  //range
  //---------
  double p_b = dot(pos_dir,boresight);
  Uvar in_sqrt = sc_pos.magnitude() *sc_pos.magnitude() * p_b * p_b;
  in_sqrt += radius * radius - sc_pos.magnitude()*sc_pos.magnitude();
  Uvar range0 = Uvar(0.0);
  if(in_sqrt.getValue() > 0.0) 
    {
      range0 = -sc_pos.magnitude() * p_b - sqrt(in_sqrt);
    }
  else
    {
      throw ErrorMessage("No intercepting point in the boresight direction");
      
    }  
  
  //-----------------------------------------
  //compare frame-method with non frame method
  //------------------------------------------
  double diff_range = range0.getInUnits("km") - range_F.getInUnits("km");
  double diff_dop   = dop0.getInUnits("Hz") - dop_F.getInUnits("Hz");  
  if (fabs(diff_range)>0.001)
    {
      throw ErrorMessage("Frame method does not produce the same range");
    }
  if (fabs(diff_dop) > 0.001)
    {
      throw ErrorMessage("Frame method does not produce the same dop");
    }
  
  
  
 
  unsigned int Ngrid = dop_freq_values.size();
  if (Ngrid != range_values.size()) throw ErrorMessage("Need square matrix");

 
	

    Ivec no_of_solution("no_solution",Ngrid*Ngrid);
    Uvec look_azimuth("look_azimuth_angle",Ngrid*Ngrid);
    Uvec look_elevation("look_elevation_angle",Ngrid*Ngrid);
    Uvec  lat("bin_patch_lat",Ngrid*Ngrid);
    Uvec  lon("bin_patch_lon",Ngrid*Ngrid);
    Uvec  range("range",Ngrid*Ngrid);
    Uvec  doppler("doppler",Ngrid*Ngrid);
    Uvec  thetai("thetai",Ngrid*Ngrid);
    Uvec  alongtrack("alongtrack",Ngrid*Ngrid);
    Uvec crosstrack("crosstrack",Ngrid*Ngrid);
    Uvec area("area_element",Ngrid*Ngrid);
    Dirvec surface_intercept("surface_intercept",Ngrid*Ngrid);
    Dvec beam_gain("beam_gain",Ngrid*Ngrid);
    Dvec sigma0_elem("sigma0_elem",Ngrid*Ngrid);
    
    //Initialize variables
    no_of_solution= 0;
    look_azimuth=Uvar(0.0);
    look_elevation= Uvar(0.0);
    lat= Uvar(0.0);
    lon= Uvar(0.0);
    range= Uvar(0.0);
    doppler= Uvar(0.0);
    thetai= Uvar(0.0);
    alongtrack=Uvar(0.0);
    crosstrack=Uvar(0.0);
    area= Uvar(0.0);
    beam_gain= 0.0;
    sigma0_elem= 0.0;

   

    //------------------------------------------------------------------
    //Find  intercepting points corresponding range and dop values
    //no_interception: store how many range and dop bins have a solution
    //------------------------------------------------------------------  
    TargetGeom tg1(t); 
    for (unsigned int i_range = 0; i_range<Ngrid; ++i_range)
    for (unsigned int j_dop =0; j_dop<Ngrid; ++j_dop)
      {
      unsigned int bin_index = i_range * Ngrid + j_dop;	
      tg1.setState(sc_state);
      tg1.setTarget(target_name,ftitan);
      tg1.setRangeDopplerInBeam3Frame(range0+range_values(i_range), dop0+dop_freq_values(j_dop),lambda);
      DirectionVector look = tg1.lookDirection();
  
      if (tg1.foundSurfaceIntercept()==false)
	{
	no_of_solution(bin_index) = 0;//no surface intercept point
	lat(bin_index)= Uvar(0,"rad");
	lon(bin_index) =  Uvar(0,"rad");
	range(bin_index)= range0+range_values(i_range);
	doppler(bin_index) = dop0+dop_freq_values(j_dop);
	thetai(bin_index)= Uvar(0,"rad"); 
	alongtrack(bin_index)=Uvar(0,"km");
	crosstrack(bin_index)=Uvar(0,"km");
	surface_intercept(bin_index)=DirectionVector();
	beam_gain(bin_index) = 0.0;
	sigma0_elem(bin_index) = 0.0;
	}
      else
	{
	no_of_solution(bin_index) = 1;   
	DirectionVector new_look= look.representIn(fbeam);
	new_look.getAzimuthElevation(look_azimuth(bin_index),
				     look_elevation(bin_index));
	lat(bin_index)= tg1.lat();
	lon(bin_index) = tg1.lon();
	range(bin_index)= tg1.range();
	doppler(bin_index) =tg1.doppler(lambda);
	thetai(bin_index)= tg1.incidenceAngle(); 
	tg1.interceptAlongCross(alongtrack(bin_index),
				crosstrack(bin_index));
	surface_intercept(bin_index)=tg1.surfaceIntercept();
	beam_gain(bin_index) =  beam.bilinear(look_azimuth(bin_index),
					      look_elevation(bin_index));
	sigma0_elem(bin_index)= muhleman_k1 
	  * cos(thetai(bin_index))
	  / pow(sin(thetai(bin_index))
		+muhleman_k2*cos(thetai(bin_index))
		,3.0);
	} 
      tg1.reset(t);
      }//Loop over i_range and j_dop
 

    //    cout<<"Finish finding all surface intercept points"<<endl;
    //    cout<<"Start calculating area elements and radar return signal "<<endl;
    //------------------------------------------	
    //End of finding all range and dop solutions
    //if there is a surface intercept point, no_of_s
    //----------------------------------------------
    //****************************************
    //    Main Ambiguity Calculation
    //there is a solution for the given  range and dop value
    //range_bin(i,i_range) and dop_bin(j,j_dop)
    //let's try to get radar return signal strenth here
    //*******************************
    for (unsigned int i_range = 0; i_range<Ngrid-1; ++i_range)
    for (unsigned int j_dop =0; j_dop<Ngrid-1; ++j_dop)
      {
      unsigned int bin_index = i_range * Ngrid + j_dop;	
      unsigned int bin_Up,  bin_Left,bin_Diag;
      bin_Up = (i_range+1)*Ngrid+j_dop;	
      bin_Left = i_range*Ngrid+(j_dop+1);
      bin_Diag = (i_range+1)*Ngrid+j_dop+1;
      if (no_of_solution(bin_index)
	  *no_of_solution(bin_Up)
	  *no_of_solution(bin_Left)
	  *no_of_solution(bin_Diag) == 1)
	{
	//------------------------------------
	//find Area element covered by (range,dop)
	//divide an area enclosed by two range and dop values
        // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
        //into two triangles and calclate their areas and add them
        //together
        //-------------------------------------           
	Uvar side_a = radius * 
	  fabs((surface_intercept(bin_index)
		.angle(surface_intercept(bin_Left))
		.getInUnits("rad")));
	Uvar side_b = radius * 
	  fabs((surface_intercept(bin_index)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	Uvar side_c = radius * 
	  fabs((surface_intercept(bin_Left)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	Uvar cos_c = (side_a*side_a + side_b*side_b - side_c*side_c)
	  /(2.0 * side_a*side_b);
	Uvar angle_c = Uvar(acos(cos_c.getValue()),"rad");
	
	area(bin_index) = 0.5 * sin(angle_c) * side_a * side_b;
	
	Uvar side_d = radius * 
	  fabs((surface_intercept(bin_Left)
		.angle(surface_intercept(bin_Diag))
		.getInUnits("rad")));
	Uvar side_e = radius * 
	  fabs((surface_intercept(bin_Diag)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	cos_c = (side_d*side_d + side_e*side_e - side_c*side_c)
	  /(2.0 * side_d*side_e);
	angle_c = Uvar(acos(cos_c.getValue()),"rad");
	area(bin_index) += 0.5 * sin(angle_c) * side_d * side_e;	    
	
	
	//-----------------------------
	//Store results into Ngrid-1 X Ngrid-1 matrix form
	//-----------------------------    
	
	//----------------------------
	//cell area
	//---------------------------
	cell_area(i_range,j_dop) = area(bin_index);
	
	//--------------------
	//cell range and doppler
	//-----------------
	range_axis(i_range) = 
	  (range(bin_index)+range(bin_Up)
	   +range(bin_Left)+range(bin_Diag))/4.0;
	

	dop_axis(j_dop) = 
	  (doppler(bin_index)+doppler(bin_Up)
	   +doppler(bin_Left)+doppler(bin_Diag))/4.0;

	//-----------------------------
	//store beam gain and sigm0
	//----------------------------
	cell_one_way_gain(i_range,j_dop) = 
	  (beam_gain(bin_index)+beam_gain(bin_Up)
	   +beam_gain(bin_Left)+beam_gain(bin_Diag))/4.0;
	
	
	cell_sigma0(i_range,j_dop) = 
	  (sigma0_elem(bin_index) + sigma0_elem(bin_Up)
	   +sigma0_elem(bin_Left) + sigma0_elem(bin_Diag))/4.0;
	
	//-----------------------------------
	//incidenceangle
	//----------------------------------
	cell_incidenceangle(i_range,j_dop) = 
	  (thetai(bin_index) + thetai(bin_Up)
	   +thetai(bin_Left)+thetai(bin_Diag))/4.0;
	
	
	//------------------------------------------------------
	//crosstrack and alongtrack distances
	//------------------------------------------------------ 
	cell_cross_track(i_range,j_dop) =
	  ( crosstrack(bin_index)+crosstrack(bin_Up)
	    +crosstrack(bin_Left)+crosstrack(bin_Diag))/4.0;
	
	cell_along_track(i_range,j_dop) = 
	  (alongtrack(bin_index)+alongtrack(bin_Up)
	   +alongtrack(bin_Left)+alongtrack(bin_Diag))/4.0;
	
	//------------------------
	//lat and lon values
	//--------------------------
	cell_lat(i_range,j_dop)=
	  (lat(bin_index) + lat(bin_Up)
	   + lat(bin_Left)+ lat(bin_Diag))/4.0;
	
	cell_lon(i_range,j_dop)=
	  (lon(bin_index) + lon(bin_Up)
	   + lon(bin_Left)+ lon(bin_Diag))/4.0;
	
	//---------------------------------------
	//find X = G^2 * lambda^2 * sigma0 / R^4
	//--------------------------------------
	
	Uvar X = 
	  cell_one_way_gain(i_range,j_dop) * cell_one_way_gain(i_range,j_dop)
	  * lambda * lambda
	  *cell_sigma0(i_range,j_dop)
	  /pow(range_axis(i_range),4);  	     
	
	cell_radar_geom_factor(i_range,j_dop)= 
	  (X*cell_area(i_range,j_dop)).getValue();  
	}
      else
	{
	cell_area(i_range,j_dop) = Uvar(0,"km km");
	range_axis(i_range) = 
	  (range(bin_index)+range(bin_Up)
	   +range(bin_Left)+range(bin_Diag))/4.0;
	
	dop_axis(j_dop) = 
	  (doppler(bin_index)+doppler(bin_Up)
	   +doppler(bin_Left)+doppler(bin_Diag))/4.0;	  
	cell_one_way_gain(i_range,j_dop) = 0.0;
	cell_sigma0(i_range,j_dop) = 0.0;
	cell_incidenceangle(i_range,j_dop) = Uvar(0,"deg");
	cell_cross_track(i_range,j_dop) = Uvar(0,"km");
	cell_along_track(i_range,j_dop)=Uvar(0,"km");
	cell_lat(i_range,j_dop) = Uvar(0,"deg");
	cell_lon(i_range,j_dop) = Uvar(0,"deg");
	cell_radar_geom_factor(i_range,j_dop) = 0.0;
	}//no solution, return zero
      }//Loop over range/dop values   

    //------------------------
    //display range and dop axis on the screen
    //-----------------------
    //for (unsigned int i = 0; i < Ngrid-1;++i)
    //{
    //cout<<"range and dop "<< range_axis(i)<<" "<<dop_axis(i)<<endl;
    //}

  }




//-------------------------------------------------
// calculate sinc pattern beam gain 
//-------------------------------------------------
double sincgain(const Uvar& azim, 
		const Uvar& elev, 
		const Uvar& azimuth_width,
		const Uvar& elevation_width,
		const double& max_gain) throw(Unit::UnitError,ErrorMessage)
  {
    Uvar au,bv;
    double sau,sbv;
    
    au = 0.88 * azim;
    au /= azimuth_width;
    bv = 0.88* elev;
    bv /= elevation_width;
    sau = sinc(au.getValue()); 
    sbv = sinc(bv.getValue());
    return( sau * sau * sbv * sbv *max_gain); 
  }












