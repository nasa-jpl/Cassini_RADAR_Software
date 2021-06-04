//----------------------------------------------------------------------------
// TargetGeom.cpp
//
// This file contains method definitions for the TargetGeom class.
// These classes use the NAIF Spice toolkit to perform most calculations.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_targetgeom_c[] =
  "@(#) $Id: TargetGeom.cpp,v 11.10 2017/03/13 22:22:25 richw Exp $";

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <iomanip>
#include "TargetGeom.h"
#include "Frame.h"
#include "Constants.h"
#include "TemplateUtils.h"

using std::cout;
using std::endl;
using std::cerr;
using std::setprecision;

//------------------------------
// Methods for class TargetGeom
//------------------------------

//--------------
// Constructors
//--------------

//-----------------------------------------------------------------------------
// TargetGeom()
//
// Default constructor.  Time needs to be set before this object
// can be used.
//-----------------------------------------------------------------------------

TargetGeom::TargetGeom()
  : state_set_(false),
    target_set_(false),
    ulook_set_(false),
    mirrorlook_set_(false),
    lat_set_(false), lon_set_(false), unormal_set_(false),
    surface_intercept_set_(false),mirror_surface_intercept_set_(false), 
    range_set_(false), rlook_set_(false),
    thetai_set_(false),mirror_thetai_set_(false),polang_set_(false),
    doppler_set_(false), rnadir_set_(false),
    along_cross_set_(false),
    ulook_call_(false), lat_call_(false), lon_call_(false),
    unormal_call_(false), surface_intercept_call_(false),
    range_call_(false), rlook_call_(false),
    thetai_call_(false), doppler_call_(false), along_cross_call_(false),
    spice_found_(false),
    time_set_(false),
    track_frame_set_(false),
    ring_target_(false),
    t_(Uvar(0,"s")),
    state_("TargetGeom::state_"),
    ulook_("TargetGeom::look_"),
    surface_intercept_("TargetGeom::surface_intercept_"),
    target_frame_("J2000","EARTH"),
    target_radii_("TargetGeom::target_radii_"),
    unormal_("TargetGeom::unormal_"),
    rlook_("TargetGeom::rlook_"),
    rnadir_("TargetGeom::rnadir_")
  {
  }

//-----------------------------------------------------------------------------
// TargetGeom(const Time& t)
//
// This constructor sets up a TargetGeom object at the input time.
//-----------------------------------------------------------------------------

TargetGeom::TargetGeom(const Time& t)
  : state_set_(false),
    target_set_(false),
    ulook_set_(false),
    mirrorlook_set_(false),
    lat_set_(false), lon_set_(false), unormal_set_(false),
    surface_intercept_set_(false), mirror_surface_intercept_set_(false),
    range_set_(false), rlook_set_(false),
    thetai_set_(false), mirror_thetai_set_(false),polang_set_(false),
    doppler_set_(false), rnadir_set_(false),
    along_cross_set_(false),
    ulook_call_(false), lat_call_(false), lon_call_(false),
    unormal_call_(false), surface_intercept_call_(false),
    range_call_(false), rlook_call_(false),
    thetai_call_(false), polang_call_(false),
    doppler_call_(false), along_cross_call_(false),
    spice_found_(false),
    time_set_(true),
    track_frame_set_(false),
    ring_target_(false),
    t_(t),
    state_("TargetGeom::state_"),
    ulook_("TargetGeom::look_"),
    surface_intercept_("TargetGeom::surface_intercept_"),
    target_frame_("J2000","EARTH"),
    target_radii_("TargetGeom::target_radii_"),
    unormal_("TargetGeom::unormal_"),
    rlook_("TargetGeom::rlook_"),
    rnadir_("TargetGeom::rnadir_")
  {
  }

//--------------
// Self Test
//--------------

bool TargetGeom::selfTest() // No exceptions
  {
/*** not needed if Frame::selfTest runs first ....
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/000728_SK_SOI_T45_82.bsp");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas00049.tsc");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/naif0007.tls");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas_radar_v11.ti");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas_v31.tf");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/pck00006.tpc");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/981005_PLTEPH-DE405S.bsp");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/naif/sat083.bsp");
  Frame::spiceLoad("/u/svejk0/richw/dat/cassini/tour_sim/naif/pdt_c_kernel_p8e.bc");
**/

try
  {
  Time t("2005-301T04:30:37.840");
  cout << "t=" << t << endl;
  Frame ftitan("IAU_TITAN","Titan");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");
  cout << "sc_state=" << sc_state << endl;

  Frame fbeam("CASSINI_RADAR_3","Cassini");
  DirectionVector boresight("boresight",fbeam,t,0,0,1);
  TargetGeom tg(t);
  tg.setState(sc_state);
  tg.setLookDirection(boresight);
  tg.setTarget("Titan",ftitan);
  Uvar alongtrack("alongtrack");
  Uvar crosstrack("crosstrack");
  tg.interceptAlongCross(alongtrack,crosstrack);
  cout << alongtrack << " " << crosstrack << endl;
  }
catch(const ErrorMessage& e)
  {
  std::cout << e.msg << std::endl;
  }
catch(...)
  {
  std::cout << "Exception caught" << std::endl;
  }
  
  return(true);
  }

//--------------
// Reset
//--------------

void TargetGeom::reset(const Time& t) // No exceptions
  {
  state_set_ = false;
  target_set_ = false;
  ulook_set_ = false;
  mirrorlook_set_ = false;
  lat_set_ = false;
  lon_set_ = false;
  unormal_set_ = false;
  surface_intercept_set_ = false;
  mirror_surface_intercept_set_ =false;
  range_set_ = false;
  rlook_set_ = false;
  thetai_set_ = false;
  mirror_thetai_set_ = false;
  polang_set_ = false;
  doppler_set_ = false;
  rnadir_set_ = false;
  along_cross_set_ = false;
  track_frame_set_ = false;
  ring_target_ = false;

  ulook_call_ = false;
  lat_call_ = false;
  lon_call_ = false;
  unormal_call_ = false;
  surface_intercept_call_ = false;
  range_call_ = false;
  rlook_call_ = false;
  thetai_call_ = false;
  polang_call_ = false;
  doppler_call_ = false;
  along_cross_call_ = false;

  spice_found_ = false;
  
  t_ = t;
  time_set_ = true;
  }

//--------------
// Predicates
//--------------

bool TargetGeom::foundSurfaceIntercept() 
  {
    if(surface_intercept_set_)
      {
	return(spice_found_);
      }
    else
      {
	surfaceIntercept();
	return(spice_found_);
      }
  }

//--------------
// Other methods
//--------------

//----------------------------------------------------------------------
// setState(s)
//
// Set the current state of the observer.
// If the time of s is different from the time of this TargetGeom,
// then an exception is thrown.
//----------------------------------------------------------------------

void TargetGeom::setState(const StateVector& s)
  {
  if (state_set_)
    {
    ErrorMessage e("TargetGeom::setState: State already set");
    e.throwMe();
    }

  if (!s.timeMatch(t_))
    {
    ErrorMessage e("TargetGeom time doesn't match observer state vector: " +
      s.name());
    e.throwMe();
    }

  state_ = s;
  state_set_ = true;
  }

//----------------------------------------------------------------------
// setState(observer_name)
//
// Sets the current state of the observer using the observer name and
// the current time of this TargetGeom.
//----------------------------------------------------------------------

void TargetGeom::setState(const string& observer_name)
  {
  if (state_set_)
    {
    ErrorMessage e("TargetGeom::setState: State already set");
    e.throwMe();
    }
  
  if(target_set_){
    Frame f("J2000",target_name_);
    f.ephemeris(state_,observer_name,t_,"NONE");
    state_set_ = true;
  }
  else{
    // This uses Saturn as a default body if available otherwise
    // uses the sun
    try{
      Frame f("J2000","Saturn");
      f.ephemeris(state_,observer_name,t_,"NONE");
      state_set_ = true;
    }
    catch(GeomError e){
      cerr << "Warning: Saturn ephemeris is unavailable using the Sun as a default body" << endl;
      Frame f("J2000","Sun");
      f.ephemeris(state_,observer_name,t_,"NONE");
      state_set_ = true;  
    }
  }
  }
 


//----------------------------------------------------------------------
// setLatLon(lat,lon)
//
// Set the current look direction to point to (lat,lon) on the target's
// surface at the current time.
//----------------------------------------------------------------------

void TargetGeom::setLatLon(const Uvar& lat, const Uvar& lon){
  if (ulook_set_)
    {
    ErrorMessage
      e("TargetGeom::setLatLon: Look direction already set");
    e.throwMe();
    }
   if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setLatLon: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: Time is not set");
      e.throwMe();
    }
  if (ring_target_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: Can't use with ring target");
      e.throwMe();
    }
  // set lat and lon
  lat_=lat;
  lon_=lon;
  lat_set_=true;
  lon_set_=true;

  // set normal direction to surface point 
  unormal_=DirectionVector("",target_frame_,t_,0,0,1);
  unormal_.setPlanetodetic(lat_,lon_);
  unormal_set_=true;
  
  // set surface intercept
  surface_intercept_=target_radii_*unormal_;
  surface_intercept_set_=true;
  spice_found_= true;
 
  // position vector in target frame
  PositionVector pos=state_.position().representIn(target_frame_);

 // compute look direction
  rlook_=surface_intercept_-pos;
  rlook_set_=true;
  ulook_=rlook_;
  ulook_set_=true;

  // check to see if spacecraft is above or below tangential plane
  // if below then (lat,lon) is over the horizon from the spacecraft
  // and then spice_found_ = false
 
  if(dot(pos,unormal_)<surface_intercept_.magnitude()){
    //TargetGeomError e(TargetGeomError::beyond_horizon);
    //e.throwMe();
    spice_found_ = false;
    surface_intercept_ = PositionVector();
    surface_intercept_set_ = true;
    unormal_ = DirectionVector();
    unormal_set_ = true;
    rlook_ = PositionVector();
    rlook_set_ = true;
    ulook_ = DirectionVector();
    ulook_set_ = true;
  }
  
  
  
}



//----------------------------------------------------------------------
// setLatLonHeight(lat,lon,ht)
//
// Set the current look direction to point to (lat,lon,ht) above the target's
// surface at the current time.
//----------------------------------------------------------------------

void TargetGeom::setLatLonHeight(const Uvar& lat, const Uvar& lon, const Uvar& ht){
  if (ulook_set_)
    {
    ErrorMessage
      e("TargetGeom::setLatLon: Look direction already set");
    e.throwMe();
    }
   if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setLatLon: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: Time is not set");
      e.throwMe();
    }
  if (ring_target_)
    {
      ErrorMessage e(
	 "TargetGeom::setLatLon: Can't use with ring target");
      e.throwMe();
    }
  // set lat and lon
  lat_=lat;
  lon_=lon;
  lat_set_=true;
  lon_set_=true;

  // set normal direction to surface point 
  unormal_=DirectionVector("",target_frame_,t_,0,0,1);
  unormal_.setPlanetodetic(lat_,lon_);
  unormal_set_=true;
  
  // set surface intercept
  surface_intercept_=(target_radii_+ht)*unormal_;
  surface_intercept_set_=true;
  spice_found_= true;
 
  // position vector in target frame
  PositionVector pos=state_.position().representIn(target_frame_);

 // compute look direction
  rlook_=surface_intercept_-pos;
  rlook_set_=true;
  ulook_=rlook_;
  ulook_set_=true;

  // check to see if spacecraft is above or below tangential plane
  // if below then (lat,lon) is over the horizon from the spacecraft
  // and then spice_found_ = false
 
  if(dot(pos,unormal_)<surface_intercept_.magnitude()){
    //TargetGeomError e(TargetGeomError::beyond_horizon);
    //e.throwMe();
    spice_found_ = false;
    surface_intercept_ = PositionVector();
    surface_intercept_set_ = true;
    unormal_ = DirectionVector();
    unormal_set_ = true;
    rlook_ = PositionVector();
    rlook_set_ = true;
    ulook_ = DirectionVector();
    ulook_set_ = true;
  }  
  
  
}

//----------------------------------------------------------------------
// setLookDirection(ulook)
//
// Set the current look direction.
// If the time of ulook is different from the time of this TargetGeom,
// then an exception is thrown.
//----------------------------------------------------------------------

void TargetGeom::setLookDirection(const DirectionVector& ulook)
  
  {
  //cout << "setLookDirection:" << endl;
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setLookDirection: Look direction already set");
    e.throwMe();
    }

  if (!ulook.timeMatch(t_))
    {
    ErrorMessage e("TargetGeom time doesn't match look direction vector: " +
      ulook.name());
    e.throwMe();
    }

  if (!target_set_)
    {
    ErrorMessage e("TargetGeom Need target set before look direction: " +
      ulook.name());
    e.throwMe();
    }

  ulook_ = ulook;
  ulook_.representIn(target_frame_);
  ulook_set_ = true;
  //cout << "ulook in = " << ulook << endl;
  //cout << "ulook out = " << ulook_ << endl;
  }


//----------------------------------------------------------------------
// setTarget()
//
// Set the current target name and related target data using default
// target information which was read form the config file in
// Frame::config
//----------------------------------------------------------------------

void TargetGeom::setTarget(){
  if (target_set_)
    {
      ErrorMessage e(
		     "TargetGeom::setTarget: Target already set");
      e.throwMe();
    }
  
  // Set ring target if config option for ring sweep is on
  ring_target_ = default_ring_target;

  // target_name_ is not set

  // Get spice id for this body
  target_id_=default_target_spice_id;

  // Setup target frame to use with this target body
  target_frame_ = 
    Frame(default_target_frame_spice_id,default_target_spice_id);

  // Get radii of this body
  target_radii_=default_target_radii;
  target_radii_.setTime(t_);

  target_set_ = true;
}

//----------------------------------------------------------------------
// setTarget(target_name)
//
// Set the current target name and related target data.
// If the name is unknown to the Spice system,
// then an exception is thrown.
// This method uses the IAU body fixed frame for the target body to
// express target quantities (eg., lat lon).
//----------------------------------------------------------------------

void TargetGeom::setTarget(const string& target_name)
  
  {
  setTarget(target_name, Frame("IAU_" + target_name, target_name));
  }

//----------------------------------------------------------------------
// setTarget(target_name,target_frame)
//
// Set the current target name and related target data.
// If the name is unknown to the Spice system,
// then an exception is thrown.
// This method uses a user supplied frame to express target quantities
// in (eg., lat, lon).
//----------------------------------------------------------------------

void TargetGeom::setTarget(const string& target_name, const Frame& target_frame)
  
  {
  if (target_set_)
    {
    ErrorMessage e(
      "TargetGeom::setTarget: Target already set");
    e.throwMe();
    }

  target_name_ = target_name;

  // Get spice id for this body
  spice_target_id(target_name,target_id_);

  // Setup target frame to use with this target body
  target_frame_ = target_frame;

  // Get radii of this body
  spice_target_radii(target_name,target_frame_,t_,target_id_,target_radii_);
  target_set_ = true;

  // Set ring target if config option for ring sweep is on
  ring_target_ = default_ring_target;
  }

//----------------------------------------------------------------------
// setAzimuthIncidence(azimuth,incidence,side,x,y,z)
//
// Set a desired attitude by specifying the incidence angle, the delta
// azimuth angle, and the side (left or right relative to velocity).
// The desired attitude is returned in 3 DirectionVectors x,y,z, that
// give the 3 coordinate axes of the spacecraft such that -z is pointed
// at a point on the surface with the desired delta azimuth and incidence angle.
// This method also sets the look direction to the -z axis direction so
// that subsequent use of this TargetGeom object will generate data for
// that same point.  Thus, this method is used in place of setLookDirection
// to fully specify a look geometry.
// Delta azimuth is the y-rotation away from the target-center pointing
// attitude.
//
// Note that a frame is not returned because this method is not intended
// to replace the Spice frame CASSINI_SC_COORD.  Rather, it is intended
// to support calculation of data to put into a ckernel to represent that
// frame.
// 
// This method just sets a look direction.  The user still needs to check
// whether there is a surface intercept point separately.
//----------------------------------------------------------------------

void TargetGeom::setAzimuthIncidence(const Uvar& azimuth,
				     const Uvar& thetai, 
				     const string& side,
				     DirectionVector& x, 
				     DirectionVector& y, 
				     DirectionVector& z)
  {
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setAzimuthIncidence: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setAzimuthIncidence: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setAzimuthIncidence: Time is not set");
      e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage
      e("TargetGeom::setSquintIncidence: Look direction already set");
    e.throwMe();
    }

  //---------------------------------------------------------------------
  // Start by pointing -z towards nadir and y orthogonal to the velocity
  // direction, putting x close the the velocity direction.
  // Note that when building DirectionVectors from PositionVectors, the
  // result depends on the origin of the PositionVector (PositionVectors
  // don't float), so care must be taken to ensure that the PositionVector
  // is in the correct frame.  In this case, state_ is forced into
  // the target frame (even though it is probably already there)
  // to make this issue obvious.
  //---------------------------------------------------------------------

  // enforce the use of ftitan
  DirectionVector sc_z("sc_z",state_.position().representIn(target_frame_));
  DirectionVector sc_y("sc_y",
    cross(state_.position().representIn(target_frame_),
          state_.velocity().representIn(target_frame_)));
  DirectionVector sc_x("sc_x",cross(sc_y,sc_z));

  Frame fvelocity("velocity_frame",
    state_.position().representIn(target_frame_),
    sc_x, sc_y, sc_z);

  //------------------------------------------------------------------
  // Apply desired delta azimuth angle by a y-rotation, slewing the
  // boresights behind/ahead of the nadir point.
  //------------------------------------------------------------------

  Uvar zero_deg("zero_deg",0.0,"rad");
  Rotation azimuth_rot(zero_deg, azimuth, zero_deg, 0,1,2);
  Frame fazimuth("delta_azimuth_frame", fvelocity, azimuth_rot);

  //------------------------------------------------------------------
  // Find the x-rotation that will give the desired incidence angle
  // using an iterative search.  Sense of the rotation is determined
  // by the choice of left or right looking (with respect to the
  // velocity direction).
  //------------------------------------------------------------------

  int sign = 0;
  if (side == "right")
    {
    sign = -1;
    }
  else if (side == "left")
    {
    sign = 1;
    }
  else
    {
    ErrorMessage e(
      "TargetGeom::setSquintIncidence: Invalid side = " + side);
    e.throwMe();
    }

  //--------------------------------------------------------------------
  // Bisection search for x-rotation that gives desired incidence angle.
  //--------------------------------------------------------------------

  Uvar roll_accuracy("roll_accuracy");
  roll_accuracy = coerce_base(0.01,"deg");
  Uvar roll_min("roll_min",0.0,"rad");
  Uvar roll_max("roll_max");
  roll_max = limbRoll();
  Uvar roll("roll");
  // Compute number of iterations to achieve roll_accuracy above
  unsigned int Nbisection =
    (unsigned int)(log((roll_max - roll_min)/roll_accuracy)/log(2.0) + 1);

  for (unsigned int i=1; i <= Nbisection; ++i)
    {
    roll = 0.5*(roll_max + roll_min);
    // Note that sign of roll is handled separately.
    Rotation inc_rot(sign*roll, zero_deg, zero_deg, 0,1,2);
    Frame fsc("spacecraft_frame", fazimuth, inc_rot);
    DirectionVector boresight("boresight",fsc,t_,0,0,-1);
    // Need to use a new TargetGeom object for each look direction because
    // TargetGeom doesn't allow the look direction to be reset.
    TargetGeom tg(t_);
    tg.setState(state_);
    tg.setTarget(target_name_,target_frame_);
    tg.setLookDirection(boresight);
    if (tg.incidenceAngle() < thetai)
      {
      roll_min = roll;
      }
    else
      {
      roll_max = roll;
      }
    }

  //--------------------------------------------------------------
  // Use the final rotated frame to specify the coordinate axes.
  //--------------------------------------------------------------

  Rotation inc_rot(sign*roll, zero_deg, zero_deg, 0,1,2);
  Frame fsc("spacecraft_frame", fazimuth, inc_rot);
  x = DirectionVector("x",fsc,t_,1,0,0);
  y = DirectionVector("y",fsc,t_,0,1,0);
  z = DirectionVector("z",fsc,t_,0,0,1);

  //--------------------------------------------------
  // Set the look direction to the -z axis direction
  //--------------------------------------------------

  ulook_ = -z;
  ulook_set_ = true;//just set look direction
  ulook_call_ = false;
  }


//----------------------------------------------------------------------
//  setIsodopplerIncidence(const Uvar& thetai, const string& side,
//   const Uvar& set_sc_z_rotation
//    DirectionVector& x, DirectionVector&y,DirectionVector& z) 
//    
//
// Set a desired attitude to produce a constant doppler shift 
//  along the beam elevation axis for a given incidence angle, altitude
//  and a target with sphere geometry.  
// This method also sets the look direction to the -z axis direction so
//  that subsequent use of this TargetGeom object will generate data for
//  that same point.  Thus, this method is used in place of setLookDirection
//  to fully specify a look geometry.
//
// Note that a frame is not returned because this method is not intended
// to replace the Spice frame CASSINI_SC_COORD.  Rather, it is intended
// to support calculation of data to put into a ckernel to represent that
// frame.
//
//
// This method just sets a look direction.  The user still needs to check
// whether there is a surface intercept point separately.
//----------------------------------------------------------------------


void TargetGeom::setIsodopplerIncidence(const Uvar& thetai, 
					const string& side,
					const Uvar& set_sc_z_rotation,
					DirectionVector& x, 
					DirectionVector& y,
					DirectionVector& z) 
  {
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setIsodopplerIncidence: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setIsodopplerIncidence: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setIsodopplerIncidence: Time is not set");
      e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setIsodopplerIncidence: Look direction already set");
    e.throwMe();
    }
  
  if (thetai >= Uvar(pi/2.0,"rad"))
    {
    ErrorMessage e(
     "TargetGeom::setIsodopplerIncidence: Incidence angle larger than 90 deg");
    e.throwMe();
    }

  //---------------------------------------------------------------------
  // Start by pointing -z towards nadir and y orthogonal to the velocity
  // direction, putting x close the the velocity direction.
  // Note that when building DirectionVectors from PositionVectors, the
  // result depends on the origin of the PositionVector (PositionVectors
  // don't float), so care must be taken to ensure that the PositionVector
  // is in the correct frame.  In this case, state_ is forced into
  // the target frame (even though it is probably already there)
  // to make this issue obvious.
  //---------------------------------------------------------------------

 
  DirectionVector sc_z("sc_z",state_.position().representIn(target_frame_));
  DirectionVector sc_y("sc_y",
		  cross(state_.position().representIn(target_frame_),
		        state_.velocity().representIn(target_frame_)));
  DirectionVector sc_x("sc_x",cross(sc_y,sc_z));

  Frame ftcn("track_crosstrack_nadir_coordinate",
    state_.position().representIn(target_frame_),
    sc_x, sc_y, sc_z);

  //----------------------------------------------------------------
  //Find a target radius and altitude
  //----------------------------------------------------------------  
  TargetGeom tg(t_);
  tg.setState(state_);
  tg.setTarget(target_name_,target_frame_);
  Uvar radius = tg.radius();
  Uvar altitude = tg.altitude();
 
  //------------------------------------------------------------------
  //Find a theta_r (angle from -N-axis) angle for a given incidence angle
  //theta_c is an angle correpsonding to a crosstrack distance from nadir
  //--------------------------------------------------------------------
  Uvar A = 1.0 + altitude/radius;
  Uvar theta_r = asin(sin(thetai)/A);
  Uvar theta_c = thetai - theta_r;
  

  //----------------------------------------------------
  //Calculate range from sc to surface intercept point
  //---------------------------------------------------
  //Uvar range = (altitude + radius* (1.0 - cos(theta_c)))/cos(theta_r);

  //------------------------------------------------------
  //sc velocity component in ftcn frame
  //-----------------------------------------------------
  FloatVector sc_vel = state_.velocity().representIn(target_frame_);
  sc_vel.representIn(ftcn);
  DirectionVector vel_dir = sc_vel;
  Uvar theta_p = asin(Uvar(vel_dir[DirectionVector::Z]));

  //------------------------------------------------------------------
  // Find a location to produce isodoppler (constant doppler)
  // along the range direction in terms of theta_a(azimuth angle in ftcn frame)
  // detailed calculations are documented in power point format
  //(ask Y. Gim)
  //-------------------------------------------------------------------
  
  double x_top =sin(thetai)*sin(theta_p)*cos(theta_r) 
              - sin(theta_c)*sin(theta_p);
  double x_bottom = sin(thetai)*cos(theta_p)*sin(theta_r) 
                  - cos(theta_c)*cos(theta_p);
  
  Uvar theta_a = Uvar(pi/2,"rad");//default value
  if(x_bottom==0.0) 
    {
      ulook_ = DirectionVector();
      ulook_call_ = false;
      spice_found_ = false;
      cout<<"No solution for the given incidence angle"<<endl;
      return;
    }
  else
    {
      if (fabs(x_top/x_bottom) <= 1.0)
	{ 
	  theta_a = acos(Uvar(x_top/x_bottom));
	}    
      else
	{
	  ulook_ = DirectionVector();
	  ulook_call_ = false;
	  spice_found_ = false;
	  cout<<"No solution for the given incidence angle"<<endl;
	  return;
	}
    }
 
  
  //----------------------------------------
  //Right looking: positive theta_a : right side w.r.t. velocity
  //Left looking: negative theta_a : left side w.r.t. velocity
  //---------------------------------------


  //-----------------------------------------------------------------
  //azimuth angle is defined from the T-axi(major velocity direction)
  //, which is negative rotation w.r.t. N-axis(-nadir direction).
  //
  //Right looking corresponds to a rotation by  (90 - theta_a) around N-axis
  //followed by  roll_rotation around T_axis
  //
  //Left looking corresponds to a rotation by -(90-theta_a) around N-axis
  // followed by roll_rot around T_axis, 
  //----------------------------------------------------------------

  //------------------------------------------------------------------
  // Sense of the rotation is determined
  // by the choice of left or right looking (with respect to the
  // velocity direction).
  //------------------------------------------------------------------

  

  Uvar zero_deg("zero_deg",0.0,"rad");
  Uvar azim_rot = Uvar(pi/2,"rad")- theta_a;//no change for right look
  Uvar roll_rot = theta_r;// no change for left look
  if (side == "right")
    {
      //right looking
      //azimuth angle rotation -- around "positive"  N-axis
      //look angle rotation --  around "negative"  T-axis
      roll_rot = -roll_rot;
    }
  else if (side == "left")
    {
      //left looking   
      //azimuth angle rotation -- around "negative" N-axis
      //look angle rotation-- around "positive" T-axis
      azim_rot = -azim_rot;   
    }
  else
    {
      ErrorMessage("look direction should be either right or left").throwMe();
      
    }
 
  Rotation azimuth_rot(zero_deg,zero_deg,azim_rot,0,1,2);
  Frame fazimuth("azimuth_rot",ftcn,azimuth_rot);
  Rotation look_rot(roll_rot,zero_deg,zero_deg,0,1,2);
  Frame fsc("space_craft",fazimuth,look_rot);
  Rotation sc_z_rotation(zero_deg,zero_deg,set_sc_z_rotation,0,1,2);
  Frame fsc_z_rotation("space craft with z_rotation",fsc,sc_z_rotation);

  x = DirectionVector("x",fsc_z_rotation,t_,1,0,0);
  y = DirectionVector("y",fsc_z_rotation,t_,0,1,0);
  z = DirectionVector("z",fsc_z_rotation,t_,0,0,1);
  
  //--------------------------------------------------
  // Set the look direction to the -z axis direction
  //--------------------------------------------------    
  ulook_ = -z;
  ulook_set_ = true;//just set look direction, does not guarantee surface interception
  ulook_call_ =false;
  }


//----------------------
// setPitchBiasedIncidenceAngle()
//
// Set a desired attitude in such a way that we can increase
// surface coverage in the alongtrack direction.  This method
// initially set an isodoppler point and roll the spacecraft
// about its y-axis away from ideal pointing.
//
// This pointing, correctly applied, is expected to produce 
// more alongtrack coverage with reduced looks at high altitude
// and about same alongtrack coverage near CA, if pitch angles decreases
// to zero at CA
//
// Testing stage: Feb. 28, 2005
//
//----------------------------------

void
TargetGeom:: setPitchBiasedIncidence(const Uvar& thetai, 
			       const string& side,
			       const Uvar& set_sc_z_rotation,
			       const Uvar& pitch,
			       DirectionVector& x, 
			       DirectionVector& y,
			       DirectionVector& z) 
{
  if (!target_set_)
    {
      ErrorMessage e(
		     "TargetGeom::setPitchBiasedIncidence: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
		     "TargetGeom::setPitchBiasedIncidence: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
		     "TargetGeom::setPitchBiasedIncidence: Time is not set");
      e.throwMe();
    }
  if (ulook_set_)
    {
      ErrorMessage e(
		     "TargetGeom::setPitchBiasedIncidence: Look direction already set");
      e.throwMe();
    }
  
  if (thetai >= Uvar(pi/2.0,"rad"))
    {
      ErrorMessage e(
		     "TargetGeom::setPitchBiasedIncidence: Incidence angle larger than 90 deg");
      e.throwMe();
    }
  
  //---------------------------------------------------------------------
  // Start by pointing -z towards nadir and y orthogonal to the velocity
  // direction, putting x close the the velocity direction.
  // Note that when building DirectionVectors from PositionVectors, the
  // result depends on the origin of the PositionVector (PositionVectors
  // don't float), so care must be taken to ensure that the PositionVector
  // is in the correct frame.  In this case, state_ is forced into
  // the target frame (even though it is probably already there)
  // to make this issue obvious.
  //---------------------------------------------------------------------
  
  
  DirectionVector sc_z("sc_z",state_.position().representIn(target_frame_));
  DirectionVector sc_y("sc_y",
		  cross(state_.position().representIn(target_frame_),
		        state_.velocity().representIn(target_frame_)));
  DirectionVector sc_x("sc_x",cross(sc_y,sc_z));

  Frame ftcn("track_crosstrack_nadir_coordinate",
    state_.position().representIn(target_frame_),
    sc_x, sc_y, sc_z);

  //----------------------------------------------------------------
  //Find a target radius and altitude
  //----------------------------------------------------------------  
  TargetGeom tg(t_);
  tg.setState(state_);
  tg.setTarget(target_name_,target_frame_);
  Uvar radius = tg.radius();
  Uvar altitude = tg.altitude();
 
  //------------------------------------------------------------------
  //Find a theta_r (angle from -N-axis) angle for a given incidence angle
  //theta_c is an angle correpsonding to a crosstrack distance from nadir
  //--------------------------------------------------------------------
  Uvar A = 1.0 + altitude/radius;
  Uvar theta_r = asin(sin(thetai)/A);
  Uvar theta_c = thetai - theta_r;
  

  //----------------------------------------------------
  //Calculate range from sc to surface intercept point
  //---------------------------------------------------
  //Uvar range = (altitude + radius* (1.0 - cos(theta_c)))/cos(theta_r);

  //------------------------------------------------------
  //sc velocity component in ftcn frame
  //-----------------------------------------------------
  FloatVector sc_vel = state_.velocity().representIn(target_frame_);
  sc_vel.representIn(ftcn);
  DirectionVector vel_dir = sc_vel;
  Uvar theta_p = asin(Uvar(vel_dir[DirectionVector::Z]));

  //------------------------------------------------------------------
  // Find a location to produce isodoppler (constant doppler)
  // along the range direction in terms of theta_a(azimuth angle in ftcn frame)
  // detailed calculations are documented in power point format
  //(ask Y. Gim)
  //-------------------------------------------------------------------
  
  double x_top =sin(thetai)*sin(theta_p)*cos(theta_r) 
              - sin(theta_c)*sin(theta_p);
  double x_bottom = sin(thetai)*cos(theta_p)*sin(theta_r) 
                  - cos(theta_c)*cos(theta_p);
  
  Uvar theta_a = Uvar(pi/2,"rad");//default value
  if(x_bottom==0.0) 
    {
      ulook_ = DirectionVector();
      ulook_call_ = false;
      spice_found_ = false;
      cout<<"No solution for the given incidence angle"<<endl;
      return;
    }
  else
    {
      if (fabs(x_top/x_bottom) <= 1.0)
	{ 
	  theta_a = acos(Uvar(x_top/x_bottom));
	}    
      else
	{
	  ulook_ = DirectionVector();
	  ulook_call_ = false;
	  spice_found_ = false;
	  cout<<"No solution for the given incidence angle"<<endl;
	  return;
	}
    }
 
  
  //----------------------------------------
  //Right looking: positive theta_a : right side w.r.t. velocity
  //Left looking: negative theta_a : left side w.r.t. velocity
  //---------------------------------------


  //-----------------------------------------------------------------
  //azimuth angle is defined from the T-axi(major velocity direction)
  //, which is negative rotation w.r.t. N-axis(-nadir direction).
  //
  //Right looking corresponds to a rotation by  (90 - theta_a) around N-axis
  //followed by  roll_rotation around T_axis
  //
  //Left looking corresponds to a rotation by -(90-theta_a) around N-axis
  // followed by roll_rot around T_axis, 
  //----------------------------------------------------------------

  //------------------------------------------------------------------
  // Sense of the rotation is determined
  // by the choice of left or right looking (with respect to the
  // velocity direction).
  //------------------------------------------------------------------

  

  Uvar zero_deg("zero_deg",0.0,"rad");
  Uvar azim_rot = Uvar(pi/2,"rad")- theta_a;//no change for right look
  Uvar roll_rot = theta_r;// no change for left look
  if (side == "right")
    {
      //right looking
      //azimuth angle rotation -- around "positive"  N-axis
      //look angle rotation --  around "negative"  T-axis
      roll_rot = -roll_rot;
    }
  else if (side == "left")
    {
      //left looking   
      //azimuth angle rotation -- around "negative" N-axis
      //look angle rotation-- around "positive" T-axis
      azim_rot = -azim_rot;   
    }
  else
    {
      ErrorMessage("look direction should be either right or left").throwMe();
      
    }
 
  Rotation azimuth_rot(zero_deg,zero_deg,azim_rot,0,1,2);
  Frame fazimuth("azimuth_rot",ftcn,azimuth_rot);

  Rotation pitch_rot(zero_deg,pitch, zero_deg,0,1,2);
  Frame fpitch("pitch rotation", fazimuth,pitch_rot);
  
  Rotation look_rot(roll_rot,zero_deg,zero_deg,0,1,2);
  Frame fsc("space_craft",fpitch,look_rot);
  
  Rotation sc_z_rotation(zero_deg,zero_deg,set_sc_z_rotation,0,1,2);
  Frame fsc_z_rotation("space craft with z_rotation",fsc,sc_z_rotation);
  
  x = DirectionVector("x",fsc_z_rotation,t_,1,0,0);
  y = DirectionVector("y",fsc_z_rotation,t_,0,1,0);
  z = DirectionVector("z",fsc_z_rotation,t_,0,0,1);
  
  //--------------------------------------------------
  // Set the look direction to the -z axis direction
  //--------------------------------------------------    
  ulook_ = -z;
  ulook_set_ = true;//just set look direction, does not guarantee surface interception
  ulook_call_ =false;
}






//----------------------------------------------------------------------
//  setContiguousSurfaceCoverage(const Uvar& thetai, const string& side,
//    DirectionVector& x, DirectionVector&y,DirectionVector& z)     
//
// Set a desired attitude so that beam projection on the surface
// is perpendicular to the motion of the spacecraft.
// This method will ensure that groundswath of one beam 3db is contiguous
// across beam
// Caution: This steering will cause non-ideal doppler frequency profile 
// along the beam elevation direction unlike setIsodopplerincidence method 
//
// This method also sets the look direction to the -z axis direction so
//  that subsequent use of this TargetGeom object will generate data for
//  that same point.  Thus, this method is used in place of setLookDirection
//  to fully specify a look geometry.
//
// Note that a frame is not returned because this method is not intended
// to replace the Spice frame CASSINI_SC_COORD.  Rather, it is intended
// to support calculation of data to put into a ckernel to represent that
// frame.
//----------------------------------------------------------------------
void TargetGeom::setContiguousSurfaceCoverage(const Uvar& thetai,
					      const string& side,
					      DirectionVector& x, 
					      DirectionVector& y,
					      DirectionVector& z) 
  {
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setContiguousSurfaceCoverage: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setContiguousSurfaceCoverage: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setContiguousSurfaceCoverage: Time is not set");
      e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setContiguousSurfaceCoverage: Look direction already set");
    e.throwMe();
    }
  
  if (thetai >= Uvar(pi/2.0,"rad"))
    {
    ErrorMessage e(
     "TargetGeom::setContiguousSurfaceCoverage: Incidence angle larger than 90 deg");
    e.throwMe();
    }

  //---------------------------------------------------------------------
  // Start by pointing -z towards nadir and y orthogonal to the velocity
  // direction, putting x close the the velocity direction.
  // Note that when building DirectionVectors from PositionVectors, the
  // result depends on the origin of the PositionVector (PositionVectors
  // don't float), so care must be taken to ensure that the PositionVector
  // is in the correct frame.  In this case, state_ is forced into
  // the target frame (even though it is probably already there)
  // to make this issue obvious.
  //--------------------------------------------------------------------- 
  DirectionVector sc_z("sc_z",state_.position().representIn(target_frame_));
  DirectionVector sc_y("sc_y",
		  cross(state_.position().representIn(target_frame_),
		        state_.velocity().representIn(target_frame_)));
  DirectionVector sc_x("sc_x",cross(sc_y,sc_z));

  Frame ftcn("track_crosstrack_nadir_coordinate",
    state_.position().representIn(target_frame_),
    sc_x, sc_y, sc_z);

  //----------------------------------------------------------------
  //Find a target radius and altitude
  //----------------------------------------------------------------  
  TargetGeom tg(t_);
  tg.setState(state_);
  tg.setTarget(target_name_,target_frame_);
  Uvar radius = tg.radius();
  Uvar altitude = tg.altitude();

  //------------------------------------------------------------------
  //Find a theta_r (angle from -N-axis) angle for a given incidence angle
  //theta_c is an angle correpsonding to a crosstrack distance from nadir
  //--------------------------------------------------------------------
  Uvar A = 1.0 + altitude/radius;
  Uvar theta_r = asin(sin(thetai)/A);
  Uvar theta_c = thetai - theta_r;
  
  //------------------------------------------------------
  //sc velocity component in ftcn frame
  //-----------------------------------------------------
  FloatVector sc_vel = state_.velocity().representIn(target_frame_);
  sc_vel.representIn(ftcn);
  DirectionVector vel_dir = sc_vel;
  Uvar theta_p = asin(Uvar(vel_dir[DirectionVector::Z]));

  //------------------------------------------------------------------
  // Find a location to produce isodoppler (constant doppler)
  // along the range direction in terms of theta_a(azimuth angle in ftcn frame)
  //-------------------------------------------------------------------
  
  double x_top =sin(theta_p) * sin(theta_c);
  double x_bottom = cos(theta_c)*cos(theta_p);
  
  Uvar theta_a = Uvar(pi/2,"rad");//default value
  if(x_bottom == 0.0)
    {
      ulook_ = DirectionVector();
      spice_found_ = false;
      ulook_call_ = false;
      cout<<"no solution "<<endl;
      return;
    }
  else
    {
      if (fabs(x_top/x_bottom) <= 1.0)
	{ 
	  theta_a = acos(Uvar(x_top/x_bottom));
	}    
      else
	{
	  ulook_ = DirectionVector();
	  spice_found_ =false;
	  ulook_call_ = false;
	  cout<<"no solution "<<endl;
	  return;
	}
    }
 
  
  //----------------------------------------
  //Right looking: positive theta_a 
  //Left looking: negative theta_a
  //---------------------------------------


  //-----------------------------------------------------------------
  //azimuth angle is defined from the T-axi(major velocity direction)
  //, which is negative rotation w.r.t. N-axis(-nadir direction).
  //Right looking corresponds to a rotation by  -(90 - theta_a) around N-axis
  //Left looking corresponds to a rotation by +(90-theta_a) around N-axis, 
  //----------------------------------------------------------------

  //------------------------------------------------------------------
  // Sense of the rotation is determined
  // by the choice of left or right looking (with respect to the
  // velocity direction).
  //------------------------------------------------------------------

  Uvar zero_deg("zero_deg",0.0,"rad");
  Uvar azim_rot = Uvar(pi/2,"rad")- theta_a;
  Uvar roll_rot = theta_r;
  if (side == "right")
    {
      //right looking
      //azimuth angle rotation -- around "positive"  N-axis
      //look angle rotation --  around "negative"  T-axis
      roll_rot = -roll_rot;
    }
  else if (side == "left")
    {
      //left looking   
      //azimuth angle rotation -- around "negative" N-axis
      //look angle rotation-- around "positive" T-axis
      azim_rot = -azim_rot;   
    }
  else
    {
      ErrorMessage e( "TargetGeom::setContiguousSurfaceCoverage: Invalid side = " 
		      + side);
      e.throwMe();
    }
 
  Rotation azimuth_rot(zero_deg,zero_deg,azim_rot,0,1,2);
  Frame fazimuth("azimuth_rot",ftcn,azimuth_rot);
  Rotation look_rot(roll_rot,zero_deg,zero_deg,0,1,2);
  Frame fsc("space_craft",fazimuth,look_rot);

  x = DirectionVector("x",fsc,t_,1,0,0);
  y = DirectionVector("y",fsc,t_,0,1,0);
  z = DirectionVector("z",fsc,t_,0,0,1);
  
  //--------------------------------------------------
  // Set the look direction to the -z axis direction
  //--------------------------------------------------    
  ulook_ = -z;
  ulook_set_ = true;//just set look direciton
  ulook_call_ = false;
  }

//----------------------------------------------------------------------
// setRangeDopplerInBeam3Frame(range,doppler,lambda)
//  This method will return a look vector defined in "beam3" corresponding
//   to given range and doppler values.  It is assumed that the target 
//   has a sphere geometry.
//----------------------------------------------------------------------

void TargetGeom::setRangeDopplerInBeam3Frame(const Uvar& range, 
					     const Uvar& doppler,
					     const Uvar& lambda)
  {
  if (range_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeDopplerInBeam3Frame: Range already set");
    e.throwMe();
    }
  if (doppler_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeDopplerInBeam3Frame: Doppler already set");
    e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Look direction already set");
    e.throwMe();
    }
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setRangeDopplerInBeam3Frame: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInBeam3Frame: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInBeam3Frame: Time is not set");
      e.throwMe();
    }
  
  range_ = range;
  range_set_ = true;
  doppler_ = doppler;
  doppler_set_ = true;
  
  double p1,p2, p3, v1,v2,v3;
  DirectionVector dir_pos,dir_vel;
  int sign = 1;
 
  //______________________________
  //Determine direction vectors of
  //sc's position and velocity
  // in beam frame (use beam 3)
  //______________________________    
  Frame fbeam("CASSINI_RADAR_3","Cassini");  
  DirectionVector beam3_bore("bore",fbeam,t_,0,0,1);
  //enforce the use of target_frame  
 
  dir_pos = state_.position();     
  dir_vel = state_.velocity();
  DirectionVector dir_cross = cross(dir_pos,dir_vel);
  double roll_direction = dot(dir_cross,beam3_bore);
  

  dir_pos.representIn(fbeam); 
  dir_vel.representIn(fbeam);

  p1 = dir_pos[DirectionVector::X];
  p2 = dir_pos[DirectionVector::Y];
  p3 = dir_pos[DirectionVector::Z];
       
  v1 = dir_vel[DirectionVector::X];
  v2 = dir_vel[DirectionVector::Y];
  v3 = dir_vel[DirectionVector::Z];

  //_______________________________________
  //Determine roll direction from direction
  //of sc_state in beam frame
  //_______________________________________
  
  if (roll_direction < 0.0) 
    {
    //roll_direction is "right";
    sign = -1;
    }
  else  
    {
    //roll_direction is "left";
    sign = 1;
    }
  
  //__________________________________
  //determine unitless dist1 and freq1
  //__________________________________
  Uvar  dist1 =(
       target_radii_.magnitude()/sqrt(3.0)*target_radii_.magnitude()/sqrt(3.0) 
       - state_.position().magnitude()*state_.position().magnitude()
       -range*range)
                   /(2.0 * range *state_.position().magnitude());
  Uvar freq1 =  doppler * lambda / (2.0 * state_.velocity().magnitude());	 	  

  //_______________________________________
  //get magnitude of each values
  //_______________________________________
  double dist = dist1.getInUnits("");
  double freq = freq1.getInUnits("");
  
  //-------------------------------------------------------------------
  //Equations to find a look direction (b1,b2,b3)
  //freq = vx * sin(azi) + vy*cos(azi)*sin(elev) + vz*cos(azi)*cos(elev)
  //dist=  px * sin(azi) + py*cos(azi)*sin(elev) + pz*cos(azi)*cos(elev)
  //--------------------------------------------------------------------
  //simplified solution
  // v1 * b1 + v2 * b2 + v3 * b3 = v_b = freq
  // p1 * b1 + p2 * b2 + p3 * b3 = p_b = dist
  // b1^2 + b2^2 + b3^2 = 1.0
  //usually v1 and p3 are close to +-1
  // b1 = (freq - v2 * b2 - v3*b3)/v1
  // p1 * (freq - v2 * b2 - v3*b3)/v1 + p2 * b2 + p3 * b3 = p_b = dist
  // b2* (p2 - p1 * v2/v1) + b3*(p3 -p1 *v3/v1)= dist - p1*freq/v1
  //*************** Danger: p1 is very small number
  // A * b2 + B * b3 = C
  //************* C is virtually independent of freq because p1 is small
  // b3 = (C - A *b2)/B = D + E*b2
  // b1 = F * b2 + G 
  // v1 * b1 + v2 * b2 + v3 *(D + E*b2) = freq
  // v1* b1 = -(v2 + v3*E) b2 +freq - v3 *D
  // b1 = -(v2/v1 + v3*E/v1) + (freq - v3*D)/v1
  // b1 = F * b2 + G 
  //************* F is very small 
  // b1^2 + b2^2 + b3^2 = 1
  // (F*b2 + G)^2 + b2^2 + (D+E*b2)^2 = 1.0
  //(F^2 + 1 + E^2)* b2^2 + 2*b2 *(FG+DE) + G^2+D^2-1 =0
  //H*b^2 + 2I *b2 + J = 0
  // b2 =[ -I +(-)sqrt(I*I - H*J)]/H	
  
  double  D,E,F,G,H,I,J;
  double sol,ans;
  double b1, b2, b3;
  //A = p2 - p1 * v2/v1;
  //B = p3 -p1 *v3/v1 ;
  //C= dist - p1*freq/v1;

  D =(dist - p1*freq/v1)/(p3 -p1 *v3/v1) ;//(C/B)
  E = -(p2 - p1 * v2/v1)/(p3 -p1 *v3/v1);  // -(A/B)
  F = -v2/v1 - v3*E/v1;
  G = freq/v1 - v3*D/v1;    
  H = F*F + 1.0 + E*E;
  I = F*G + D*E;
  J = G*G + D*D - 1.0;    
  sol = I*I - H * J;
  

  if (sol >= 0.0)
    {
    spice_found_ = true;
    ans = (-I + double(sign)* sqrt(sol))/H;   
    b2 = ans;
    b1 = F*b2 + G;
    b3 = D + E*b2;    
    ulook_ = DirectionVector("TargetGeom::look_",fbeam,t_,b1,b2,b3); 
    ulook_set_ = true;
    surfaceIntercept();
    //As a matter of fact, there is another solution that produces 
    // the same range and doppler shift.  Its location, with respect to
    // the look direction,  is located
    // at the opposite side.  This is what we call "Mirror Site."
    //For ambiguity calculation, it is specially important when the incidence
    // angle is small and therefore its pointing is very close to nadir.
    //In this case, original look direction and mirror direction is very close
    // and therby ambiguity is as strong as signal
   
    ans = (-I - double(sign)* sqrt(sol))/H;   
    b2 = ans;
    b1 = F*b2 + G;
    b3 = D + E*b2;    
    mirrorlook_ = DirectionVector("TargetGeom::look_",fbeam,t_,b1,b2,b3); 
    mirrorlook_set_ = true;
    }
  else
    {     
    spice_found_ = false;
    ulook_ = DirectionVector();//no solution
    }
  ulook_call_ = false;
  }



//----------------------------------------------------------------------
// setRangeDopplerInTargetFrame(range,doppler,lambda,roll_direction)
//  This method will find surface intercept points for a given set of 
// (range,doppler) and spacecraft roll direction
//  roll direction information can be obtained from flyby object: 
//   string s = flyby.LookDirection, s = "Right" or "Left"
//  look direction (from spacecrft to surface point) will be defined 
// in the target frame
//----------------------------------------------------------------------
 void TargetGeom::setRangeDopplerInTargetFrame(const Uvar& range, 
					       const Uvar& doppler, 
					       const Uvar& lambda)
  {
  if (range_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeDopplerInTargetFrame: Range already set");
    e.throwMe();
    }
  if (doppler_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeDopplerInTargetFrame: Doppler already set");
    e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Look direction already set");
    e.throwMe();
    }
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setRangeDopplerInTargetFrame: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: Time is not set");
      e.throwMe();
    }
  //  Do not Delete
  // solve in the body fixed frame
  //basic assumption: Target is a sphere
  // Surface_intercept_vector = sc_state.position() + Range*unit_look_vector
  // Doppler = 2 *( sc_state.velocity * unit_look_vector)/lambda
  // sc_state.position() = sc_state.position().magnitude()*(p1,p2,p3)
  // sc_state.velocity() = sc_state.velocity().magnitude()*(v1,v2,v3)


    DirectionVector pos = state_.position();
    DirectionVector vel = state_.velocity();
    DirectionVector dir_cross = cross(vel,pos);
    
    //______________________________
    //Determine direction vectors of
    //beam3 boresight
    //______________________________    
    Frame fbeam("CASSINI_RADAR_3","Cassini");  
    DirectionVector beam3_bore("bore",fbeam,t_,0,0,1);
    double roll_direction=dot(beam3_bore,dir_cross);
    //cout<<"roll direction +(right) -(left) "<<roll_direction<<endl;

    double p1,p2,p3, v1,v2,v3;
    double u1,u2,u3;
   
    p1 = pos[DirectionVector::X];
    p2 = pos[DirectionVector::Y];
    p3 = pos[DirectionVector::Z];

    v1 = vel[DirectionVector::X];
    v2 = vel[DirectionVector::Y];
    v3 = vel[DirectionVector::Z];
    
   

    double range_in_radius = (range
		     /(target_radii_.magnitude()/sqrt(3.0))).getInUnits("");
    double position_in_radius = (state_.position().magnitude()
		     /(target_radii_.magnitude()/sqrt(3.0))).getInUnits("");

    double A =  (1 - position_in_radius*position_in_radius 
		 - range_in_radius*range_in_radius)
      /(2.0 * range_in_radius * position_in_radius);
    double B =  (lambda * doppler/(2 * state_.velocity().magnitude())).getInUnits("");

    double C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
    double D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
    double E = (A - p1 *C)/p3;
    double F = -(p2 + p1*D)/p3;
    
    double sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);
    if (sol > 0)
      {
      spice_found_ = true;
      u2 = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
      u1 = C + D*u2;
      u3 = E + F*u2;
    
      ulook_ = DirectionVector("TargetGeom::look",target_frame_,t_,u1,u2,u3);  
      surface_intercept_ = state_.position();
      surface_intercept_ += range* ulook_;
      DirectionVector to_surface = surface_intercept_;
      double look_confirm=dot(dir_cross,to_surface);
      if(look_confirm*roll_direction >0.0) 
	{
	  double a = dot(-ulook_,to_surface);
	  thetai_ = Uvar(acos(a),"rad");
	  u2 = ( -(C*D + E*F) +  sqrt( sol))/(D*D + F*F + 1.0);
	  u1 = C + D*u2;
	  u3 = E + F*u2;
	  mirrorlook_ = DirectionVector("TargetGeom::mirrorlook",target_frame_,t_,u1,u2,u3);
	  mirror_surface_intercept_ = state_.position();
	  mirror_surface_intercept_ += range*mirrorlook_; 
	  to_surface = mirror_surface_intercept_;//reuse to_surface vector
	  a = dot(-mirrorlook_,to_surface);
	  mirror_thetai_ = Uvar(acos(a),"rad");	
	}
      else
	{
	  //cout<<"alt sol "<<endl;
	  u2 = ( -(C*D + E*F) +sqrt( sol))/(D*D + F*F + 1.0);
	  u1 = C + D*u2;
	  u3 = E + F*u2;
	  ulook_ = DirectionVector("TargetGeom::look",target_frame_,t_,u1,u2,u3);  
	  surface_intercept_ = state_.position();
	  surface_intercept_ += range* ulook_;
	  DirectionVector to_surface = surface_intercept_;
	  double a = dot(-ulook_,to_surface);
	  thetai_ = Uvar(acos(a),"rad");

	  u2 = ( -(C*D + E*F) - sqrt( sol))/(D*D + F*F + 1.0);
	  u1 = C + D*u2;
	  u3 = E + F*u2;
	  mirrorlook_ = DirectionVector("TargetGeom::mirrorlook",target_frame_,t_,u1,u2,u3);
	  mirror_surface_intercept_ = state_.position();
	  mirror_surface_intercept_ += range*mirrorlook_; 
	  to_surface = mirror_surface_intercept_;//reuse to_surface vector
	  a = dot(-mirrorlook_,to_surface);
	  mirror_thetai_ = Uvar(acos(a),"rad");	
	}
      ulook_set_ = true;
      ulook_call_ = false;    
      surface_intercept_set_ = true;
      surface_intercept_call_ = false;    
      thetai_set_ = true;
      thetai_call_ = false;
      mirrorlook_set_ = true;
      mirror_surface_intercept_set_ = true;   
      mirror_thetai_set_ = true;
      }
    else
      {     
	spice_found_ = false;
	ulook_ = DirectionVector();//no solution
      }
    ulook_call_ = false;//do not have to call again
  }


//----------------------------------------------------------------------
// setRangeAzimuth(range,azimuth,fbeam)
//  This method will find elevation angle for a given set of 
// (range,azimuth,fbeam).
// look direction (from spacecrft to surface point) will be defined 
// in the target frame
//----------------------------------------------------------------------
 void TargetGeom::setRangeAzimuth(const Uvar& range, 
				  const Uvar& azi0, 
				  const Frame& fbeam)
  {
  if (range_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Range already set");
    e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Look direction already set");
    e.throwMe();
    }
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setRangeDopplerInTargetFrame: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: Time is not set");
      e.throwMe();
    }
  
  Frame ftitan("IAU_TITAN","Titan");

  //---- Beam axis in titan body coordinates ----
  double M_T2n[3][3];
  DirectionVector aa, bb, cc;
  ftitan.axialVectors(fbeam,t_,aa,bb,cc);
  M_T2n[0][0] = aa[DirectionVector::X];
  M_T2n[0][1] = aa[DirectionVector::Y];
  M_T2n[0][2] = aa[DirectionVector::Z];

  M_T2n[1][0] = bb[DirectionVector::X];
  M_T2n[1][1] = bb[DirectionVector::Y];
  M_T2n[1][2] = bb[DirectionVector::Z];

  M_T2n[2][0] = cc[DirectionVector::X];
  M_T2n[2][1] = cc[DirectionVector::Y];
  M_T2n[2][2] = cc[DirectionVector::Z];
  
  // --- Nadir pointing frame in titan body coordinates ---
  double px = M_T2n[2][0];
  double py = M_T2n[2][1];
  double pz = M_T2n[2][2];
  
  DirectionVector P0 = -nadir();

  // define beam direction
  DirectionVector P_n = -nadir();
  P_n.representIn(fbeam);
  Uvar nadir_azi,nadir_elev;
  P_n.getAzimuthElevation(nadir_azi,nadir_elev);
  int beam_direction = 1;
  if(nadir_elev>Uvar(0.,"rad")) beam_direction = -1;
  //cout<<"beam direction "<<beam_direction<<endl;
 
  double cx = P0[DirectionVector::X];
  double cy = P0[DirectionVector::Y];
  double cz = P0[DirectionVector::Z];
  
  double b0 = sqrt(pow(cy*pz-py*cz,2.)+pow(cz*px-pz*cx,2.)+pow(cx*py-px*cy,2.));
  double bx = (cy*pz-py*cz)/b0;
  double by = (cz*px-pz*cx)/b0;
  double bz = (cx*py-px*cy)/b0;
  
  double ax = by*cz-cy*bz;
  double ay = bz*cx-cz*bx;
  double az = bx*cy-cx*by;
  
  // Matrix product M_T2n * Trans(T_T2r)
  double m11 = M_T2n[0][0]*ax +  M_T2n[0][1]*ay +  M_T2n[0][2]*az;
  double m12 = M_T2n[0][0]*bx +  M_T2n[0][1]*by +  M_T2n[0][2]*bz;
  double m13 = M_T2n[0][0]*cx +  M_T2n[0][1]*cy +  M_T2n[0][2]*cz;
  double m21 = M_T2n[1][0]*ax +  M_T2n[1][1]*ay +  M_T2n[1][2]*az;
  double m22 = M_T2n[1][0]*bx +  M_T2n[1][1]*by +  M_T2n[1][2]*bz;
  double m23 = M_T2n[1][0]*cx +  M_T2n[1][1]*cy +  M_T2n[1][2]*cz;
  double m31 = M_T2n[2][0]*ax +  M_T2n[2][1]*ay +  M_T2n[2][2]*az;
  double m32 = M_T2n[2][0]*bx +  M_T2n[2][1]*by +  M_T2n[2][2]*bz;
  double m33 = M_T2n[2][0]*cx +  M_T2n[2][1]*cy +  M_T2n[2][2]*cz;  
    
  //---- Calculate xi ----
  Uvar Xsc=state_.position()[PositionVector::X];
  Uvar Ysc=state_.position()[PositionVector::Y];
  Uvar Zsc=state_.position()[PositionVector::Z];
  Uvar R0 = sqrt(pow(Xsc,2.)+pow(Ysc,2.)+pow(Zsc,2.));
  Uvar RT = radius();

  if ((range<=R0-RT)||(range>=sqrt(R0*R0-RT*RT)))
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Range is either too large or too small");
    e.throwMe();
    }
  
  double xi = ((pow(R0,2.) + pow(range,2.) - pow(RT,2.)) / (2.*range*R0)).getInUnits("");
  //cout << "   xi = "<<xi<<endl;
  
  //----- Calculate elevation ------
  double a = (m11-m31*tan(azi0))*sqrt(1-xi*xi);
  double b = (m12-m32*tan(azi0))*sqrt(1-xi*xi);
  double c = (m13-m33*tan(azi0))*xi;
  
  
  double sol = a*a + b*b -c*c;
  if (sol >= 0)
    {
      spice_found_ = true;
      double sin_phi = (double(beam_direction)*a*sqrt(sol) - b*c)/(a*a + b*b);
      double cos_phi = -(b*sin_phi+c)/a;
      //cout<<"sin_phi: "<<sin_phi<<" cos_phi: "<<cos_phi<<endl;
      double sin_ele0 = (m21*sqrt(1-xi*xi)*cos_phi + m22*sqrt(1-xi*xi)*sin_phi + m23*xi);
      Uvar ele0 = Uvar(asin(sin_ele0),"rad");
      //cout <<"   azi0= "<<azi0<<"  ele0= "<<ele0<<endl; 
      
      //set look vector from sc to surface in target frame
      ulook_ = DirectionVector("TargetGeom::look",fbeam,t_,0,0,1);
      ulook_.setAzimuthElevation(azi0,ele0);
    }
  else
    {     
      spice_found_ = false;
      ulook_ = DirectionVector();//no solution
    }
  ulook_set_ = true;
  ulook_call_ = false;
  }


//----------------------------------------------------------------------
// setRangeElevation(range,elevation,fbeam)
//  This method will find elevation angle for a given set of 
// (range,azimuth,fbeam).
// look direction (from spacecrft to surface point) will be defined 
// in the target frame
//----------------------------------------------------------------------
 void TargetGeom::setRangeElevation(const Uvar& range, 
				  const Uvar& ele0, 
				  const Frame& fbeam)
  {
  if (range_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Range already set");
    e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Look direction already set");
    e.throwMe();
    }
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setRangeDopplerInTargetFrame: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRangeDopplerInTargetFrame: Time is not set");
      e.throwMe();
    }
  
  Frame ftitan("IAU_TITAN","Titan");
  //---- Beam axis in titan body coordinates ----
  double M_T2n[3][3];
  DirectionVector aa, bb, cc;
  ftitan.axialVectors(fbeam,t_,aa,bb,cc);
  M_T2n[0][0] = aa[DirectionVector::X];
  M_T2n[0][1] = aa[DirectionVector::Y];
  M_T2n[0][2] = aa[DirectionVector::Z];

  M_T2n[1][0] = bb[DirectionVector::X];
  M_T2n[1][1] = bb[DirectionVector::Y];
  M_T2n[1][2] = bb[DirectionVector::Z];

  M_T2n[2][0] = cc[DirectionVector::X];
  M_T2n[2][1] = cc[DirectionVector::Y];
  M_T2n[2][2] = cc[DirectionVector::Z];
  
  // --- Nadir pointing frame in titan body coordinates ---
  double px = M_T2n[2][0];
  double py = M_T2n[2][1];
  double pz = M_T2n[2][2];
  
  DirectionVector P0 = -nadir();

  // define beam direction
  DirectionVector P_n = -nadir();
  P_n.representIn(fbeam);
  Uvar nadir_azi,nadir_elev;
  P_n.getAzimuthElevation(nadir_azi,nadir_elev);
  int beam_direction = 1;
  if(nadir_elev>Uvar(0.,"rad")) beam_direction = -1;
  //cout<<"beam direction "<<beam_direction<<endl;
 
  double cx = P0[DirectionVector::X];
  double cy = P0[DirectionVector::Y];
  double cz = P0[DirectionVector::Z];
  
  double b0 = sqrt(pow(cy*pz-py*cz,2.)+pow(cz*px-pz*cx,2.)+pow(cx*py-px*cy,2.));
  double bx = (cy*pz-py*cz)/b0;
  double by = (cz*px-pz*cx)/b0;
  double bz = (cx*py-px*cy)/b0;
  
  double ax = by*cz-cy*bz;
  double ay = bz*cx-cz*bx;
  double az = bx*cy-cx*by;
  
  // Matrix product M_T2n * Trans(T_T2r)
  double m11 = M_T2n[0][0]*ax +  M_T2n[0][1]*ay +  M_T2n[0][2]*az;
  double m12 = M_T2n[0][0]*bx +  M_T2n[0][1]*by +  M_T2n[0][2]*bz;
  double m13 = M_T2n[0][0]*cx +  M_T2n[0][1]*cy +  M_T2n[0][2]*cz;
  double m21 = M_T2n[1][0]*ax +  M_T2n[1][1]*ay +  M_T2n[1][2]*az;
  double m22 = M_T2n[1][0]*bx +  M_T2n[1][1]*by +  M_T2n[1][2]*bz;
  double m23 = M_T2n[1][0]*cx +  M_T2n[1][1]*cy +  M_T2n[1][2]*cz;
  //double m31 = M_T2n[2][0]*ax +  M_T2n[2][1]*ay +  M_T2n[2][2]*az;
  //double m32 = M_T2n[2][0]*bx +  M_T2n[2][1]*by +  M_T2n[2][2]*bz;
  //double m33 = M_T2n[2][0]*cx +  M_T2n[2][1]*cy +  M_T2n[2][2]*cz;  
    
  //---- Calculate xi ----
  Uvar Xsc=state_.position()[PositionVector::X];
  Uvar Ysc=state_.position()[PositionVector::Y];
  Uvar Zsc=state_.position()[PositionVector::Z];
  Uvar R0 = sqrt(pow(Xsc,2.)+pow(Ysc,2.)+pow(Zsc,2.));
  Uvar RT = radius();

  if ((range<=R0-RT)||(range>=sqrt(R0*R0-RT*RT)))
    {
    ErrorMessage e(
      "TargetGeom::setRangeAzimuth: Range is either too large or too small");
    e.throwMe();
    }
  
  double xi = ((pow(R0,2.) + pow(range,2.) - pow(RT,2.)) / (2.*range*R0)).getInUnits("");
  //cout << "   xi = "<<xi<<endl;

  
  //----- Calculate azimuth -----
  double a = m21*sqrt(1-xi*xi);
  double b = m22*sqrt(1-xi*xi);
  double c = m23*xi - sin(ele0);

  
  double sol = a*a + b*b -c*c;
  if (sol >= 0)
    {
      spice_found_ = true;
      double sin_phi = (double(beam_direction)*a*sqrt(sol) - b*c)/(a*a + b*b);
      double cos_phi = -(b*sin_phi+c)/a;
      
      double sin_azi0 = (m11*sqrt(1-xi*xi)*cos_phi + m12*sqrt(1-xi*xi)*sin_phi + m13*xi);
      Uvar azi0 = Uvar(asin(sin_azi0/cos(ele0)),"rad");

      //cout <<"   azi0= "<<azi0<<"  ele0= "<<ele0<<endl; 
      
      //set look vector from sc to surface in target frame
      ulook_ = DirectionVector("TargetGeom::look",fbeam,t_,0,0,1);
      ulook_.setAzimuthElevation(azi0,ele0);
    }
  else
    {     
      spice_found_ = false;
      ulook_ = DirectionVector();//no solution
    }
  ulook_set_ = true;
  ulook_call_ = false;
  }




//----------------------------------------------------------------------
// setRangeAzimuth(range,azimuth,i_beam)
//  This method will find elevation angle for a given set of 
// (range,azimuth,fbeam).
// look direction (from spacecrft to surface point) will be defined 
// in the target frame
//----------------------------------------------------------------------
 void TargetGeom::setRangeAzimuth(const Uvar& range, 
				  const Uvar& azi0, 
				  const unsigned int& beam_number)
  {
    Frame fbeam(beam_frame_spice_id[beam_number-1],cassini_spice_id);
    setRangeAzimuth(range, azi0, fbeam);
  }


//----------------------------------------------------------------------
// setRangeElevation(range,elevation,i_beam)
//  This method will find azimuth angle for a given set of 
// (range,azimuth,fbeam).
// look direction (from spacecrft to surface point) will be defined 
// in the target frame
//----------------------------------------------------------------------
 void TargetGeom::setRangeElevation(const Uvar& range, 
				  const Uvar& ele0, 
				  const unsigned int& beam_number)
  {
    Frame fbeam(beam_frame_spice_id[beam_number-1],cassini_spice_id); 
    setRangeElevation(range,ele0,fbeam);
  }



//-----------------------------------
//void setAlongCross(const Uvar& alongtrack, const Uvar& crosstrack)
//
// It is assumed that the target has a sphere geometry
//-----------------------------------
void TargetGeom::setAlongCross(const Uvar& alongtrack, const Uvar& crosstrack)
  
{
  if(!target_set_)
    {ErrorMessage("TargetGeom::setAlongCross:Target is not set").throwMe();}
  if(!state_set_)
    {ErrorMessage("TargetGeom::setAlongCross:State is not set").throwMe();}
  if(!time_set_)
    {ErrorMessage("TargetGeom::setAlongCross:Time is not set").throwMe();}
  if(along_cross_set_)
    {
      ErrorMessage e("TargetGeom::setAlongCross:along_cross is already set");
      e.throwMe();
    }
  if(!rnadir_set_)
    {nadir();}
  if(!track_frame_set_)
    {
    ErrorMessage e("TargetGeom::setAlongCross: track frame is not set");
    e.throwMe();
    }
  along_ = alongtrack;
  cross_ = crosstrack;
  along_cross_set_ = true;


  ulook_ =DirectionVector("u",track_frame_,t_,0,0,1);
  Uvar along_in_angle = Uvar((-alongtrack/target_radii_[PositionVector::X]).getInUnits(""),"rad");
  Uvar cross_in_angle = Uvar((crosstrack/target_radii_[PositionVector::Z]).getInUnits(""),"rad");
  if ( fabs(along_in_angle.getInUnits("deg")) > 90 || fabs(cross_in_angle.getInUnits("deg"))> 90)
   {
   ErrorMessage e("Alongtrack/Crosstrack position is out of sight from spacecraft");
   e.throwMe();
   }
  //debug
  //cout<<"inside target geom "<< cross_in_angle<<" "<<along_in_angle<<endl;
  ulook_.setPlanetodetic(cross_in_angle,along_in_angle);
  ulook_.representIn(target_frame_);
 
  surface_intercept_ = radius() * ulook_;
  spice_found_ =true;
  surface_intercept_set_ = true; 
  ulook_set_ = true;
  ulook_call_ = false;
  }

//----------------------------------------------------------------------
//  setRingRadiusSweepClosest(t0,t1,t2,t3,time_pad,
//    radius0,radius1,radius2,radius3,
//    x_dir,y_dir,z_dir);
//
// Set a desired attitude to point the HGA (-z dir) towards the desired
// ring radius which is closest to the spacecraft at a given time.
// The desired radius is specified by a piecewise linear fit to the
// supplied time,radius points.
//   t0,radius0 - start of first segment
//   t1,radius1 - end of first segment, start of second segment
//   t2,radius2 - end of second segment
//   t3,radius3 - end of third segment
//   t4,radius4 - end of fourth segment
//   time_pad - extra time on ends padded with inertial extensions
//
// The spacecraft X-axis is kept parallel to the ring plane. 
//----------------------------------------------------------------------


void TargetGeom::setRingRadiusSweepClosest
  (
  const Uvar& t_epoch, 
  const Uvar& t0, 
  const Uvar& t1, 
  const Uvar& t2, 
  const Uvar& t3, 
  const Uvar& time_pad, 
  const Uvar& radius0, 
  const Uvar& radius1, 
  const Uvar& radius2, 
  const Uvar& radius3, 
  DirectionVector& x, 
  DirectionVector& y,
  DirectionVector& z
  ) 
  {
  bool in_pad = false;
  if (!target_set_)
    {
      ErrorMessage e(
        "TargetGeom::setRingRadiusSweepClosest: Target is not set");
      e.throwMe();
    }
  if (!state_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRingRadiusSweepClosest: State is not set");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e(
	 "TargetGeom::setRingRadiusSweepClosest: Time is not set");
      e.throwMe();
    }
  if (ulook_set_)
    {
    ErrorMessage e(
      "TargetGeom::setRingRadiusSweepClosest: Look direction already set");
    e.throwMe();
    }
 
  Frame fj2000("J2000",target_name_);
  Uvar radius = radius0;
  Uvar trel = t_ - t_epoch;
  if (trel < t0-time_pad)
    {
    ErrorMessage e(
      "TargetGeom::setRingRadiusSweepClosest: Time too early");
    e.throwMe();
    }
  else if (trel < t0)
    {  // In the starting time pad - use vectors for t0
    //cout << "in starting time pad" << endl;
    //cout << "trel = " << trel << endl;
    //cout << "t0 = " << t0 << endl;
    // Need to keep recursive call from ending up here again.
    Time new_t = t0 + Uvar(1,"ms") + t_epoch;
    //cout << "new_t - t_epoch = " << new_t - t_epoch << endl;
    TargetGeom tg(new_t);
    tg.setTarget(target_name_);
    // Make s/c position in target_frame_ for convenient use later
    StateVector sc_state("sc_state");
    target_frame_.ephemeris(sc_state,"Cassini",new_t,"NONE");
    tg.setState(sc_state);
    tg.setRingRadiusSweepClosest(t_epoch,t0,t1,t2,t3,time_pad,
      radius0,radius1,radius2,radius3,x,y,z);
    // Rebrand the result vectors with the current time
    x.setTime(t_);
    y.setTime(t_);
    z.setTime(t_);
    in_pad = true;
    }
  else if (trel < t1)
    {
    if (t1 == t0 && radius1 == radius0)
      {
      radius = radius0;
      }
    else
      {
      Uvar slope = (radius1 - radius0)/(t1 - t0);
      radius = slope*(trel - t0) + radius0;
      }
    }
  else if (trel < t2)
    {
    Uvar slope = (radius2 - radius1)/(t2 - t1);
    radius = slope*(trel - t1) + radius1;
    }
  else if (trel <= t3)
    {
    if (t2 == t3 && radius2 == radius3)
      {
      radius = radius2;
      }
    else
      {
      Uvar slope = (radius3 - radius2)/(t3 - t2);
      radius = slope*(trel - t2) + radius2;
      }
    }
  else if (trel <= t3+time_pad)
    {  // In the ending time pad - use the vectors at t3.
    //cout << "in ending time pad" << endl;
    // Need to keep recursive call from ending up here again.
    Time new_t = t3 - Uvar(1,"ms") + t_epoch;
    TargetGeom tg(new_t);
    tg.setTarget(target_name_);
    // Make s/c position in target_frame_ for convenient use later
    StateVector sc_state("sc_state");
    target_frame_.ephemeris(sc_state,"Cassini",new_t,"NONE");
    tg.setState(sc_state);
    tg.setRingRadiusSweepClosest(t_epoch,t0,t1,t2,t3,time_pad,
      radius0,radius1,radius2,radius3,x,y,z);
    // Rebrand the result vectors with the current time
    x.setTime(t_);
    y.setTime(t_);
    z.setTime(t_);
    in_pad = true;
    }
  else if (trel > t3+time_pad)
    {
    //cout << "trel = " << trel << endl;
    //cout << "t3 = " << t3 << endl;
    ErrorMessage e(
      "TargetGeom::setRingRadiusSweepClosest: Time too late");
    e.throwMe();
    }
  
/*
  cout << "t_epoch = " << t_epoch << endl;
  cout << "t_ = " << t_ << endl;
  cout << "trel = " << trel << endl;
  cout << "t0 = " << t0 << endl;
  cout << "t1 = " << t1 << endl;
  cout << "t2 = " << t2 << endl;
  cout << "t3 = " << t3 << endl;
  cout << "radius0 = " << radius0 << endl;
  cout << "radius1 = " << radius1 << endl;
  cout << "radius2 = " << radius2 << endl;
  cout << "radius3 = " << radius3 << endl;
  cout << "radius = " << radius << endl;
*/
  if (in_pad)
    {  // x,y,z already set above
    // Finish setting up this TargetGeom object
    ulook_ = -z;
    }
  else
    {
    surface_intercept_ = state_.position();
    surface_intercept_[PositionVector::Z] = 0.0;
    DirectionVector utarget("utarget",surface_intercept_);
    surface_intercept_ = radius*utarget;
    PositionVector rlook_target = surface_intercept_ - state_.position();
    DirectionVector ulook_target(rlook_target);

    DirectionVector zb_p(ulook_target);
    // Keep s/c x-axis parallel to ring plane throughout
    DirectionVector zp("zp",target_frame_,t_,0,0,1);
    DirectionVector xb_p("xb_p",cross(zb_p,zp));
    DirectionVector yb_p("yb_p",cross(zb_p,xb_p));
/*
    cout << "target_frame_ name = " << target_frame_.name() << endl;
    cout << "state.position = " << state_.position() << endl;
    cout << "state.position (target_frame) = " << state_.position().representIn(target_frame_) << endl;
    cout << "utarget = " << utarget << endl;
    cout << "rlook_target = " << rlook_target << endl;
    cout << "ulook_target = " << ulook_target << endl;
    cout << "xb_p = " << xb_p << endl;
    cout << "yb_p = " << yb_p << endl;
    cout << "zb_p = " << zb_p << endl;
*/
  
    x = xb_p;
    y = -yb_p;
    z = -zb_p;

    // Put in inertial coordinates using the current time
    x.representIn(fj2000);
    y.representIn(fj2000);
    z.representIn(fj2000);

    ulook_ = -ulook_target;
    }
  
  //--------------------------------------------------
  // Set the look direction to the -z axis direction
  //--------------------------------------------------    

  ulook_set_ = true;
  ulook_call_ =false;
  ring_target_ = true;
  //spice_found_ = true;
  //surface_intercept_set_ = true; 

/*
  cout << "rtarget = " << surface_intercept_ << endl;
  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  cout << "z = " << z << endl;
*/
  }


//-------------------------------------------------------------------
//setTrackFrame(const Frame& ftrack)
// frame = getTrackFrame()
//
// Set/get the frame used to define along and cross track coordinates.
//-------------------------------------------------------------------

void TargetGeom::setTrackFrame(const Frame& ftrack)
  {
  if (track_frame_set_)
    {
    ErrorMessage("TargetGeom::setTrackFrame: already set").throwMe();
    }
  track_frame_ = ftrack;
  track_frame_set_ = true;
  }

Frame TargetGeom::getTrackFrame()
  {
  if (!track_frame_set_)
    {
    ErrorMessage("TargetGeom::getTrackFrame: track frame not set").throwMe();
    }
  return(track_frame_);
  }

//---------------------------------------------------------------------
// lat()
//
// Compute planetodetic latitude of surface intercept point.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

Uvar TargetGeom::lat()
  {
  if (lat_set_)
    {
    return(lat_);
    }

  if (lat_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute lat");
    e.throwMe();
    }
  lat_call_ = true;

  if (!unormal_set_)
    {
    normal();
    }

  if (foundSurfaceIntercept())
    {
    if (ring_target_)
      {  // Latitude is always 0 in ring plane
      DirectionVector ovec(surface_intercept_);
      ovec.getPlanetodetic(lat_,lon_);
      }
    else
      {
      unormal_.getPlanetodetic(lat_,lon_);
      }
    }
  else
    {
    lat_ = Uvar(0,"rad");
    lon_ = Uvar(0,"rad");
    }
  lat_set_ = true;
  lon_set_ = true;
  lat_call_ = false;
  return(lat_);  
  }

//---------------------------------------------------------------------
// lon()
//
// Compute planetodetic longitude of surface intercept point.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

Uvar TargetGeom::lon()
  {
  if (lon_set_)
    {
    return(lon_);
    }

  if (lon_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute lon");
    e.throwMe();
    }
  lon_call_ = true;

  if (!unormal_set_)
    {
    normal();
    }

  if (foundSurfaceIntercept())
    {
    if (ring_target_)
      {  // Latitude is always 0 in ring plane
      DirectionVector ovec(surface_intercept_);
      ovec.getPlanetodetic(lat_,lon_);
      }
    else
      {
      unormal_.getPlanetodetic(lat_,lon_);
      }
    }
  else
    {
    lat_ = Uvar(0,"rad");
    lon_ = Uvar(0,"rad");
    }
  lat_set_ = true;
  lon_set_ = true;
  lon_call_ = false;
  return(lon_);  
  }

//---------------------------------------------------------------------
// normal()
//
// Compute normal direction vector at the surface intercept point.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

DirectionVector TargetGeom::normal()
  {
  //cout << "normal:" << endl;
  if (unormal_set_)
    {
    return(unormal_);
    }

  if (unormal_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute normal");
    e.throwMe();
    }
  unormal_call_ = true;

  if (!target_set_)
    {
    ErrorMessage e("TargetGeom: Need target set to compute normal");
    e.throwMe();
    }

  if (!surface_intercept_set_)
    {
    surfaceIntercept();
    }

  if (!target_set_)
    {
    ErrorMessage e("TargetGeom: Not enough data to compute normal");
    e.throwMe();
    }

  if (foundSurfaceIntercept())
    {
    if (ring_target_)
      {  // Ring plane normal sign set according to which side the s/c is on
      if (state_.position()[PositionVector::Z] > 0.0)
        {
        unormal_ = DirectionVector(target_frame_,t_,0,0,1);
        }
      else
        {
        unormal_ = DirectionVector(target_frame_,t_,0,0,-1);
        }
      }
    else
      {
      PositionVector rnormal(surface_intercept_);
      unormal_=rnormal / target_radii_;
      }
    }
  else
    {
    unormal_ = DirectionVector();
    }
  unormal_set_ = true;
  unormal_call_ = false;
  return(unormal_);
  }
    
//---------------------------------------------------------------------
// surfaceIntercept()
//
// Compute the surface intercept point for the current look geometry.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

PositionVector TargetGeom::surfaceIntercept()
  {
  //cout << "surfaceIntercept:" << endl;
  if (surface_intercept_set_)
    {
    return(surface_intercept_);
    }

  if (surface_intercept_call_)
    {  // Stop infinite recursion
      ErrorMessage e("TargetGeom: Not enough data to compute surface_intercept");
      e.throwMe();
    }
  surface_intercept_call_ = true;

  if (!ulook_set_)
    {
    lookDirection();
    }

  if (!(state_set_ && target_set_))
    {
    ErrorMessage e("TargetGeom: Not enough data to compute surface intercept");
    e.throwMe();
    }

  // Call spice interface routine to do the work
  if (ring_target_)
    {
    cout << "surface_intercept on ring target" << endl;
    cout << "state_.position = " << state_.position() << endl;
    cout << "ulook = " << ulook_ << endl;
    spice_ringplane_intercept(surface_intercept_,spice_found_,
      state_.position(),ulook_,target_radii_);
    cout << "surface_intercept = " << surface_intercept_ << endl;
    cout << "spice_found = " << spice_found_ << endl;
    }
  else
    {
    spice_surfpt(surface_intercept_,spice_found_,state_.position(),
      ulook_,target_radii_);
    }

  surface_intercept_set_ = true;
  surface_intercept_call_ = false;
  return(surface_intercept_);
  }

//---------------------------------------------------------------------
// mirrorsurfaceIntercept()
//
// Return mirrorsurfaceIntercept() calculated in setRangeDopplerInTargetFrame
//---------------------------------------------------------------------
PositionVector TargetGeom::mirrorsurfaceIntercept()
  {
    if (!mirror_surface_intercept_set_)
      {
	ErrorMessage("TargetGeom::mirrorsurfaceIntercept: no mirror surface intercept is set").throwMe();
      }
    return(mirror_surface_intercept_);
  }    

//---------------------------------------------------------------------
// surfaceVelocity()
//
// Compute the velocity of the surface intercept point along the surface.
// ie., the component of the velocity orthogonal to the look direction.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

FloatVector TargetGeom::surfaceVelocity()
  {
  if (!(state_set_ && target_set_))
    {
    ErrorMessage e("TargetGeom: Not enough data to compute surface velocity");
    e.throwMe();
    }

  if (!unormal_set_)
    {
    normal();
    }

  FloatVector v;
  if (foundSurfaceIntercept())
    {
    // Form surface normal velocity compoenent by projecting s/c velocity
    // on surface normal direction
    v = unormal_ * dot(unormal_,state_.velocity());
    // The surface tangential velocity is the rest of the velocity vector after
    // the normal component is subtracted out.
    v = state_.velocity() - v;
    }
  else
    {  // no surface intercept, so just return 0
    v = FloatVector();
    }

  return(v);
  }

//---------------------------------------------------------------------
// range()
//
// Compute the range from the observer to the surface intercept point
// for the current look geometry.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//---------------------------------------------------------------------

Uvar TargetGeom::range()
  {
  if (range_set_)
    {
    return(range_);
    }

  if (range_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute range");
    e.throwMe();
    }
  range_call_ = true;

  if (!rlook_set_)
    {
    lookVector();
    }

  if (foundSurfaceIntercept())
    {
    range_ = rlook_.magnitude();
    }
  else
    {
    range_ = Uvar(0,"km");
    }
  range_set_ = true;
  range_call_ = false;
  return(range_);
  }

//-------------------------------------------------------------------------
// lookVector()
//
// Compute the look vector from the surface intercept point to the
// observer for the current look geometry.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//-------------------------------------------------------------------------

PositionVector TargetGeom::lookVector()
  {
  //cout << "lookVector:" << endl;
  if (rlook_set_)
    {
    return(rlook_);
    }

  if (rlook_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute look vector");
    e.throwMe();
    }
  rlook_call_ = true;

  if (!surface_intercept_set_)
    {
    surfaceIntercept();
    }

  if (!state_set_)
    {
    ErrorMessage e("TargetGeom: Not enough data to compute look vector");
    e.throwMe();
    }

  if (foundSurfaceIntercept())
    {
    rlook_ = surface_intercept_ - state_.position();
    }
  else
    {
    rlook_ = PositionVector();
    }
  rlook_set_ = true;
  rlook_call_ = false;
  return(rlook_);
  }

//-------------------------------------------------------------------------
// lookDirection()
//
// Compute the look direction from the surface intercept point to the
// observer for the current look geometry.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//-------------------------------------------------------------------------

DirectionVector TargetGeom::lookDirection()
  {
  //cout << "lookDirection:" << endl;
  if (ulook_set_)
    {
    return(ulook_);
    }

  if (ulook_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute look direction");
    e.throwMe();
    }
  ulook_call_ = true;

  if (!rlook_set_)
    {
    lookVector();
    }

  if (foundSurfaceIntercept())
    {
    ulook_ = rlook_;
    }
  else
    {
    ulook_ = DirectionVector();
    }
  ulook_set_ = true;
  ulook_call_ = false;
  return(ulook_);
  }

//-------------------------------------------------------------------------
// mirrorlookDirection()
//
// Return mirrorlookDirection only if setRangeDoppler sets
// the mirrorlook direction
//-------------------------------------------------------------------------

DirectionVector TargetGeom::mirrorlookDirection()
  {
  if (!mirrorlook_set_ )
    {
    ErrorMessage e("Mirror look is not determined by TargetGeom::setRangeDoppler method");
    e.throwMe();
    }  
  return(mirrorlook_);
  }

//-------------------------------------------------------------------------
// incidenceAngle()
//
// Compute the incidence angle which is the angle between the normal at
// the surface intercept point and the look vector.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//-------------------------------------------------------------------------

Uvar TargetGeom::incidenceAngle()
  {
  //cout << "incidenceAngle:" << endl;
  if (thetai_set_)
    {
    return(thetai_);
    }

  if (thetai_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute thetai");
    e.throwMe();
    }
  thetai_call_ = true;

  if (!ulook_set_)
    {
    lookDirection();
    }

  if (!unormal_set_)
    {
    normal();
    }

  if (foundSurfaceIntercept())
    {
    thetai_ = Uvar(acos(dot(-ulook_,unormal_)),"rad");
    }
  else
    {
    thetai_ = Uvar(0,"rad");
    }
  thetai_set_ = true;
  thetai_call_ = false;
/*
  cout << "incidenceAngle" << endl;
  cout << "  state = " << state_ << endl;
  cout << "  surface_intercept = " << surface_intercept_ << endl;
  cout << "  rlook = " << rlook_ << endl;
  cout << "  ulook = " << ulook_ << endl;
  cout << "  unormal = " << unormal_ << endl;
  cout << "  thetai = " << thetai_.getInUnits("deg") << endl;
*/
  return(thetai_);
  }

//-------------------------------------------------------------------------
// mirrorincidenceAngle()
//
//Return mirror incidence angle if calculated in setRangeDopplerInTargetFrame
//-------------------------------------------------------------------------

Uvar TargetGeom::mirrorincidenceAngle()
  {
    if(!mirror_thetai_set_)
      {
	ErrorMessage("TargetGeom::mirrorincidenceAngle:no mirror incidence angel is set").throwMe();
      }
    return(mirror_thetai_);
  }

//-------------------------------------------------------------------------
// polarizationAngle()
//
// Compute the polarization angle which is the angle between the E-field
// (also the +x axis of the spacecraft) and the plane of incidence.
// Thus, pure H-pol has a polarization angle of 0 deg, and pure V-pol has
// a polarization angle of 90 deg.
// The plane of incidence is the plane containing the normal at
// the surface intercept point and the look vector.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//-------------------------------------------------------------------------

Uvar TargetGeom::polarizationAngle()
  {
  if (polang_set_)
    {
    return(polang_);
    }

  if (polang_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute polarization angle");
    e.throwMe();
    }
  polang_call_ = true;

  if (!ulook_set_)
    {
    lookDirection();
    }

  if (!unormal_set_)
    {
    normal();
    }

  if (foundSurfaceIntercept())
    {
    DirectionVector dir_h = cross(ulook_, unormal_);
    DirectionVector dir_v = cross(dir_h, ulook_);
    // Note that efield is in direction of s/c x-axis
    DirectionVector efield("efield",Frame("CASSINI_SC_COORD","Cassini"),
      t_,1,0,0);
    polang_ = Uvar(atan2(dot(dir_v,efield),dot(dir_h,efield)),"rad");
    }
  else
    {
    polang_ = Uvar(0,"rad");
    }
  polang_set_ = true;
  polang_call_ = false;
  return(polang_);
  }

//-------------------------------------------------------------------------
// doppler()
//
// Compute the two-way doppler shift between the observer and
// the surface intercept point.
// This is part of the staged set of look geometry calculations which
// are performed as needed.
//-------------------------------------------------------------------------

Uvar TargetGeom::doppler(const Uvar& lambda)
  {
  if (doppler_set_)
    {
    return(doppler_);
    }

  if (doppler_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute doppler");
    e.throwMe();
    }
  doppler_call_ = true;

  if (!ulook_set_)
    {
    lookDirection();
    }

  if (!state_set_)
    {
     ErrorMessage e("TargetGeom: Not enough data to compute doppler");
     e.throwMe();
    }

  //-------------------------------------------------------------------------
  // Technically, we don't need to compute the surface intercept because
  // the doppler doesn't depend on it (because we are working in the
  // target body-fixed frame).  However, without a surface intercept, there
  // wouldn't be any echo with a doppler to care about.  Thus, it seems
  // to be useful to have doppler fail (return 0) when no surface intercept
  // is found.
  //-------------------------------------------------------------------------

  if (!surface_intercept_set_)
    {
    surfaceIntercept();
    }

  if (!unormal_set_)
    {
    normal();
    }

  if (foundSurfaceIntercept())
    {
    if (ring_target_)
      {
      // Assuming zero eccentriciy, and central mass dominance the orbital
      // speed is:  v = sqrt(GM/r) at radius r for central mass M
      // G is gravitational constant.
      cout << "TargetGeom::doppler" << endl;
      cout << "  state_ = " << state_ << endl;
      cout << "  ulook_ = " << ulook_ << endl;
      cout << "  ulook_.frame = " << ulook_.frameName() << endl;
      cout << "  surface_intercept_ = " << surface_intercept_ << endl;
      cout << "  surface_intercept_.frame = " << surface_intercept_.frameName() << endl;
      cout << "  unormal_ = " << unormal_ << endl;
      cout << "  unormal.frame = " << unormal_.frameName() << endl;
      Frame fj2000("J2000","Saturn");
      Uvar gm;
      spice_target_value("Saturn", "GM", gm);
      Uvar radius = surface_intercept_.magnitude();
      Uvar speed = sqrt(gm/radius);
      // Make ring particle velocity vector in J2000 frame.
      DirectionVector upos(surface_intercept_);
      cout << "  upos = " << upos << endl;
      cout << "  upos.frame = " << upos.frameName() << endl;
      FloatVector target_vel = speed *
        cross(unormal_.representIn(fj2000),upos.representIn(fj2000));
      // Doppler comes from relative motion in the look direction.
      // All figured in J2000 frame.
      FloatVector vrel = state_.velocity().representIn(fj2000) - target_vel;
      Uvar vrange_bore = dot(vrel,ulook_.representIn(fj2000));

      // accounting for special relativity
      Uvar carrier_freq = speed_light/lambda;
      doppler_ = 2.0*vrange_bore*carrier_freq/(speed_light-vrange_bore);
      cout << setprecision(16) << "  t_ = " << t_ << endl;
      cout << "  state_ = " << state_ << endl;
      cout << "  ulook_ = " << ulook_ << endl;
      cout << "  ulook_.frame = " << ulook_.frameName() << endl;
      cout << "  surface_intercept_ = " << surface_intercept_ << endl;
      cout << "  surface_intercept_.frame = " << surface_intercept_.frameName() << endl;
      cout << "  rlook_ = " << rlook_ << endl;
      cout << "  rnadir_ = " << rnadir_ << endl;
      cout << "  range_ = " << range_ << endl;
      cout << "  thetai_ = " << thetai_ << endl;
      cout << "  radius = " << radius << endl;
      cout << "  speed = " << speed << endl;
      cout << "  upos = " << upos << endl;
      cout << "  upos.frame = " << upos.frameName() << endl;
      cout << "  unormal_ = " << unormal_ << endl;
      cout << "  unormal.frame = " << unormal_.frameName() << endl;
      cout << "  vrel = " << vrel << endl;
      cout << "  vrange_bore = " << vrange_bore << endl;
      cout << "  doppler = " << doppler_ << endl;
      }
    else
      {  // rigid body, use body fixed frame
      Uvar vrange_bore = dot(state_.velocity(),ulook_);

      // minus sign is already in ulook_ (because it is directed towards target)
      // doppler_ = 2.0*vrange_bore/lambda; // classical version

      // accounting for special relativity
      Uvar carrier_freq = speed_light/lambda;
      doppler_ = 2.0*vrange_bore*carrier_freq/(speed_light-vrange_bore);
      }
    }
  else
    {
    doppler_ = Uvar(0,"Hz");
    }

  doppler_set_ = true;
  doppler_call_ = false;
  return(doppler_);
  }

//-------------------------------------------------------------------------
// nadir()
//
// Compute the subspacecraft point on the target body using spice routine.
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

PositionVector TargetGeom::nadir()
  {
  if (rnadir_set_)
    {
    return(rnadir_);
    }

  if (!state_set_ || !target_set_)
    {
    ErrorMessage e("TargetGeom: Not enough data to compute nadir point");
    e.throwMe();
    }

  // Call spice interface routine to do the work
  spice_nearpt(rnadir_,altitude_,state_.position(),target_radii_);

  rnadir_set_ = true;

  return(rnadir_);
  }

//-------------------------------------------------------------------------
// altitude()
//
// Compute the altitude of the spacecraft above the subspacecraft point
// on the target body.
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

Uvar TargetGeom::altitude()
  {
  if (rnadir_set_)
    {
    return(altitude_);
    }
  else
    {
    nadir();
    return(altitude_);
    }
  }

//-------------------------------------------------------------------------
// nadirLatLon()
//
// Compute the lat,lon of the nadir intercept point.
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

void TargetGeom::nadirLatLon(Uvar& lat, Uvar& lon)
  {
  if (!rnadir_set_)
    {
    nadir();
    }
  DirectionVector unormal = rnadir_ / target_radii_;
  unormal.getPlanetodetic(lat,lon);
  }

//-------------------------------------------------------------------------
// nadirAlongCross()
//
// Compute the alongtrack and crosstrack coordinates of the nadir
// intercept point.
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

void TargetGeom::nadirAlongCross(Uvar& along, Uvar& cross)
  {
  if (!rnadir_set_)
    {
    nadir();
    }
  DirectionVector unormal = rnadir_.representIn(track_frame_) / target_radii_;
  unormal.getPlanetodetic(cross,along);
  // We assume spherical target, but this could be generalized to handle
  // an ellipsoidal target....
  cross = cross.getInUnits("rad")*target_radii_[PositionVector::Z];
  // Minus to put positive along track in velocity direction
  along = -along.getInUnits("rad")*target_radii_[PositionVector::X];
  }

//-------------------------------------------------------------------------
// nadirDoppler()
//
// Compute the two-way doppler shift between the observer and
// the nadir surface intercept point.
//-------------------------------------------------------------------------

Uvar TargetGeom::nadirDoppler(const Uvar& lambda)
  {
  if (!state_set_)
    {
    ErrorMessage e("TargetGeom: Not enough data to compute nadir doppler");
    e.throwMe();
    }

  if (!rnadir_set_)
    {
    nadir();
    }

  DirectionVector nadir_look = rnadir_ - state_.position();
  Uvar vrange_nadir = dot(state_.velocity(),nadir_look);
  // minus sign is already in ulook_ (because it is directed towards target)
  return(2.0*vrange_nadir/lambda);
  }

//-------------------------------------------------------------------------
// limbRoll()
//
// Compute the x-rotation needed to put the look vector tangential to the
// surface (ie., a limb grazing look vector).
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

Uvar TargetGeom::limbRoll()
  {
  if (!rnadir_set_)
    {
    nadir();
    }

  //------------------------------------------------------
  // This solution is only valid for a spherical target!
  //------------------------------------------------------
 
  return(asin(target_radii_[PositionVector::X] /
	       (target_radii_[PositionVector::X] + altitude_)
             )
        );
  }

//-------------------------------------------------------------------------
// limbRange()
//
// Compute the range to the limb points.
// This is a special purpose geometry routine.
//-------------------------------------------------------------------------

Uvar TargetGeom::limbRange()
  
  {
  if (!rnadir_set_)
    {
    nadir();
    }

  //------------------------------------------------------
  // This solution is only valid for a spherical target!
  //------------------------------------------------------

  return(sqrt(pow(target_radii_[PositionVector::X] + altitude_, 2) -
    pow(target_radii_[PositionVector::X],2)));
  }

//-------------------------------------------------------------------------
// interceptAlongCross()
//
// Compute the along track and cross track position of the surface intercept
// point in a fixed "swathframe."
// The reason why we decide to use swathframe is to get a time_independent 
// along/cross track values for a flyby.  
// This is part of the staged set of look geometry calculations which
// are performed as needed.
// This routine assumes a spherical target body.
//-------------------------------------------------------------------------

void TargetGeom::interceptAlongCross(Uvar& alongtrack, Uvar& crosstrack)
  
  {
  if (along_cross_set_)
    {
    alongtrack = along_;
    crosstrack = cross_;
    }

  if (along_cross_call_)
    {  // Stop infinite recursion
    ErrorMessage e("TargetGeom: Not enough data to compute along track");
    e.throwMe();
    }

  if(!track_frame_set_)
    {// no swath frame to determine along/cross
    ErrorMessage e("TargetGeom::interceptAlongCross: no track frame is set");
    e.throwMe();
    }

  along_cross_call_ = true;

  if (!surface_intercept_set_)
    {
    surfaceIntercept();
    }

  if (foundSurfaceIntercept())
    {   
    DirectionVector si = surface_intercept_;  // change local copy only
    DirectionVector u("u",si.representIn(track_frame_));
    u.getPlanetodetic(cross_,along_);
    // We assume spherical target, but this could be generalized to handle
    // an ellipsoidal target....
    cross_ = cross_.getInUnits("rad")*target_radii_[PositionVector::Z];
    // Minus to put positive along track in velocity direction
    along_ = -along_.getInUnits("rad")*target_radii_[PositionVector::X];
    }
  else
    {
    along_ = Uvar(0,"km");
    cross_ = Uvar(0,"km");
    }
  along_cross_set_ = true;
  along_cross_call_ = false;
  alongtrack = along_;
  crosstrack = cross_;
  }


//--------------------------------------
// radius() - return radius of the target
//-------------------------------------

Uvar TargetGeom::radius()  
  {
  if (!target_set_)
    {
    ErrorMessage e("TargetGeom::radius(): no target set");
    e.throwMe();
    }
  return(target_radii_.magnitude()/sqrt(3.0));    
  }



//------------------------------------------
//compute doppler spread
//-------------------------------------------
 Uvar TargetGeom::dopplerSpreadBeam3(const Uvar& lambda,const Uvar& beam_width)
  {
  if (!target_set_)
    {
      ErrorMessage e("TargetGeom::radius(): no target set");
      e.throwMe();
    }
  
  if (!state_set_)
    {
      ErrorMessage e("TargetGeom: Not enough data to compute look vector");
      e.throwMe();
    }
  if (!time_set_)
    {
      ErrorMessage e("TargetGeom::setRangeDopplerInBeam3Frame: Time is not set");
      e.throwMe();
    }
  
  Uvar return_value(0,"Hz");
  Uvar diameter =2.0*radius();
  Uvar beam_size= state_.position().magnitude()*beam_width.getInUnits("rad");
  
  //Frame fbeam("CASSINI_RADAR_3","CASSINI");
  Uvar speed=state_.velocity().magnitude();
  state_.position().representIn(target_frame_);
  state_.velocity().representIn(target_frame_);
  DirectionVector dir_pos= state_.position();
  DirectionVector dir_vel=state_.velocity();
  Uvar angle_vel_center = dir_vel.angle(-dir_pos);

  if(beam_size< diameter){
    Uvar delta_theta = beam_width;
    Uvar angle1 = (2.0*angle_vel_center + delta_theta)/2.0;
    Uvar angle2 = (2.0*angle_vel_center - delta_theta)/2.0;
    return_value= 2.0*speed/lambda*fabs((cos(angle1)-cos(angle2)));
  }
  else{//when beam is larger than disk: beam size > diameter
    Uvar delta_theta = Uvar( (diameter/state_.position().magnitude()).getInUnits(""),"rad");
    Uvar angle1 = (2.0*angle_vel_center + delta_theta)/2.0;
    Uvar angle2 = (2.0*angle_vel_center - delta_theta)/2.0;
    return_value = 2.0*speed/lambda*fabs((cos(angle1) - cos(angle2)));
  }
  return(return_value);
  }

//---------------------------------------------------------------
// radialDistance()
//   Return distance in xy plane from origin to the target point
//---------------------------------------------------------------

Uvar TargetGeom::radialDistance()  
  {
  if (!target_set_)
    {
    ErrorMessage e("TargetGeom::radius(): no target set");
    e.throwMe();
    }
  PositionVector rsurf(surface_intercept_);
  rsurf[PositionVector::Z] = 0.0;
  return(rsurf.magnitude());    
  }

//---------------------------------------------------------------
// isoRangeRadiusAngle()()
//   Return angle between iso-range and iso-radius lines at the
//   target point in the xy plane of the body fixed frame.
//---------------------------------------------------------------

Uvar TargetGeom::isoRangeRadiusAngle()  
  {
  if (!target_set_)
    {
    ErrorMessage e("TargetGeom::radius(): no target set");
    e.throwMe();
    }
  if (foundSurfaceIntercept())
    {
    DirectionVector unormal_xy("unormal_xy",target_frame_,t_,0,0,1);
    DirectionVector u_iso_range = cross(unormal_xy,ulook_);
    DirectionVector usurf(surface_intercept_);
    DirectionVector u_iso_radius = cross(unormal_xy,usurf);
    return(u_iso_range.angle(u_iso_radius));
    }
  else
    {  // No surface intercept so return 0.
    return(Uvar(0.0,"rad"));
    }
  }


//-------------------------------------
// TargetGeomError method definitions 
//-------------------------------------

TargetGeomError::TargetGeomError(TargetGeomError::errorE err_type) 
  // No exceptions
  : error_type(err_type)
  {
  if (error_type == unspecified)
    msg = "Unspecified Target Geometry Error";
  else if (error_type == beyond_horizon)
    msg = "Target Geometry Error: Target location beyond horizon.";
  }

TargetGeomError::TargetGeomError(const string& emsg, 
  TargetGeomError::errorE err_type) 
  //No exceptions
  : error_type(err_type)
  {
  msg = emsg;
  }

void TargetGeomError::throwMe() 
  {
  throw *this;
  }

