//----------------------------------------------------------------------------
// CassiniSim.cpp
//
// This file contains method definitions for the CassiniSim class.
// CassiniSim provides spacecraft related data and calculations.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_cassini_sim_c[] =
  "@(#) $Id: CassiniSim.cpp,v 11.7 2016/09/20 17:12:28 richw Exp $";

#include "CassiniSim.h"
#include "Constants.h"
#include "Array.h"
#include "DebugInfo.h"
#include<iostream>

using std::cout;
using std::cerr;
using std::endl;

//---------------
// Constructors
//---------------
CassiniSim::CassiniSim()
: time_set_(false), sc_state_("sc_state"), cur_beam_(-1)
  {
    configured_=false;
  }

CassiniSim::CassiniSim(Config& cfg)
  : time_set_(false), sc_state_("sc_state"), cur_beam_(-1)
  {
  config(cfg);
  configured_=true;
  }
  
CassiniSim::CassiniSim(unsigned int beamnum, Config& cfg)
  : time_set_(false), sc_state_("sc_state"), cur_beam_(-1)
  { 
  config(cfg);
  setBeam(beamnum);
  }

//---------------------------------------------------------------------------- 
// config(cfg)
//
// Reads in spacecraft related config parameters including beam related info.
//---------------------------------------------------------------------------- 

void CassiniSim::config(Config& cfg)
  {
  // Compute lambda
  Uvar frequency("frequency");
  frequency = cfg["carrier_frequency"];
  lambda_ = Uvar("speed_light")/frequency;

  spacecraft_ = "Cassini";
  target_name_ = default_target_name;
  if (target_name_ != "none")
    {
    ftarget_ = Frame(default_target_frame_spice_id, default_target_spice_id);
    }
  else
    {  // if not target, use J2000 for ephemeris calls
    ftarget_ = Frame(j2000_frame_spice_id, sun_spice_id);
    }

  // Setup beams
  bool read_beam_pattern;
  string beam_pattern_source = cfg.str("beam_pattern_source");
  if (beam_pattern_source == "file")
    {
    read_beam_pattern = true;
    }
  else
    {
    read_beam_pattern = false;
    }

  beam_.resize(5);
  for (unsigned int i=1; i <= 5; ++i)
    {
    beam_[i-1] = Beam(i,cfg,read_beam_pattern);
    }
  configured_=true;
  }

//------------------
// Get/Set methods
//------------------

void CassiniSim::getDopplerProfile(const Uvar& lambda_chirp, 
				   Uvar& slope, Uvar& doppler_boresight,
				   Uvar& range_boresight, Uvar& rate_slope,
				   Uvar& doppler_rate_boresight)
{
  // construct objects once if possible
  static Uvar range2, doprate2, doppler2, nextdop_bs, nextdop2;
  static PositionVector pos1, pos2;
  static TargetGeom tg(t_);
  static Time t2;
  static Uvar dt(1,"s");

  ///---------------------------------------------------------
  // Check to make sure CassiniSim has been properly initialized
  //---------------------------------------------------------
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::getDopplerProfile: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::getDopplerProfile: Error, Can't do this with target = none");
    e.throwMe();
    }

  //-----------------------------------------------------------------------//
  // Step 1: Compute two ranges (boresight and 0.2 degrees off boresight)  //
  //-----------------------------------------------------------------------//


  // compute boresight range
  range_boresight = tg_.range();


  // compute range for 0.2 deg difference in elevation
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  Uvar azim=0;
  Uvar elev=Uvar(0.2*degtorad,"rad");
  Frame beam_frame=Frame(beam_frame_spice_id[cur_beam_],cassini_spice_id);
  DirectionVector newlook(beam_frame,t_,0,0,0);
  newlook.setAzimuthElevation(azim,elev);
  tg.setLookDirection(newlook);
  range2= tg.range();
 
  //-----------------------------------------------------------------------//
  // Step 2: Compute doppler centroid and surface position vector          //
  // for each range                                                        //
  //-----------------------------------------------------------------------//

  getDopplerCentroid(lambda_chirp,range_boresight, pos1, doppler_boresight);
  getDopplerCentroid(lambda_chirp,range2,pos2,doppler2);

  //-----------------------------------------------------------------------//
  // Step 3: Compute slope of doppler vs. range                            //
  //-----------------------------------------------------------------------//
  slope=(doppler2-doppler_boresight)/(range2-range_boresight);

  
  //-----------------------------------------------------------------------//
  // Step 3: Compute doppler at ranges and time t_ + 1 second              //
  //-----------------------------------------------------------------------//

  t2=t_+dt;
  getFixedTargetDopplerAtTime(lambda_chirp,pos1,t2,nextdop_bs);
  getFixedTargetDopplerAtTime(lambda_chirp,pos2,t2,nextdop2);

  //-----------------------------------------------------------------------//
  // Step 4: Numerically estimate dop_rates at each range                  //
  //-----------------------------------------------------------------------//

  doppler_rate_boresight=(nextdop_bs-doppler_boresight)/dt;
  doprate2=(nextdop2-doppler2)/dt;

  //-----------------------------------------------------------------------//
  // Step 5: Compute slope of doppler_rate vs. range                       //
  //-----------------------------------------------------------------------// 

  rate_slope=(doprate2-doppler_rate_boresight)/(range2-range_boresight);

}

void 
CassiniSim::getDopplerCentroid(const Uvar& lambda_chirp,
			       const Uvar& range, PositionVector& pos,
			       Uvar& fdc)
{
  // Construct intermediate objects once
  static TargetGeom tg;
  static DirectionVector newlook;
  static Uvar azim, elev;

  //------------------------------------------------
  // compute doppler centroid for a given range
  //-------------------------------------------------

  // first iteration get elevation for range and azim=0
  azim=0;
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setRangeAzimuth(range,azim,cur_beam_+1);
  newlook=tg.lookDirection();
  newlook.getAzimuthElevation(azim,elev);

  // second iteration get azimuth at maximum gain and elevation
  azim=beam_[cur_beam_].maxGainLineToAzim(elev);

  // compute doppler
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  newlook.setAzimuthElevation(azim,elev);
  tg.setLookDirection(newlook);
  fdc=tg.doppler(lambda_chirp);

  // compute position
  pos=tg.surfaceIntercept();
  pos.representIn(ftarget_);      
}

//-----------------------------------------------------------------
// Routine computes doppler at some arbitrary time for a fixed point
// on the target.
// Input PositionVector must be in target body fixed frame
// If not routine gives wrong answer.
// User beware!
//-----------------------------------------------------------------
void CassiniSim::getFixedTargetDopplerAtTime(const Uvar& lambda_chirp, 
				 const PositionVector& surface_pos,
				 const Time& t, Uvar& fdop){

  // construct intermediate objects once
  static TargetGeom tg;
  static PositionVector same_pos_different_time, sc_pos, look_vector;
  static DirectionVector look_direction;
  static StateVector scstate;


  
  // change time of position vector
  same_pos_different_time=surface_pos;
  same_pos_different_time.setTime(t);
  
  // compute state vector at time t
  double et;
  t.getEt(et);
  ftarget_.ephemeris(scstate,cassini_spice_id,et);
  
  // compute look direction to fixed target at time t
  sc_pos=scstate.position();
  look_vector=same_pos_different_time-sc_pos;
  look_direction=look_vector;

  tg.reset(t);
  tg.setTarget();
  tg.setState(scstate);
  tg.setLookDirection(look_direction);
  
  // get doppler
  fdop=tg.doppler(lambda_chirp);  
}				    

void CassiniSim::setMidPulseImpactTime(
     const Time& pulse_transmit_time,
     const Uvar& lat, 
     const Uvar& lon){

  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  setBeamTimeLatLon(cur_beam_+1,pulse_transmit_time,lat,lon);

  // iterate once to get reasonably accurate time
  Uvar range = tg_.range();
  Uvar delay = range/speed_light;
  setBeamTimeLatLon(cur_beam_+1,pulse_transmit_time+delay,lat,lon);
  
}

//----------------------------------------------//
// setBeamTimeLatLon(beamnum,t,lat,lon)
// set beam time and look direction with latlon
// maintains correct order and resets targetgeom
// correctly                                    
//-----------------------------------------------//

void CassiniSim::setBeamTimeLatLon(unsigned int beam_no, const Time& t, 
				   const Uvar& lat, const Uvar& lon)
{
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  setTime(t);
  // set Beam Number but don't make TargetGeom point to boresight
  // which is why setBeam() is not used.
  cur_beam_ = beam_no - 1;
  if(cur_beam_<0 || cur_beam_> 4){
    ErrorMessage e("CassiniSim: Bad beam number in setBeamTimeLatLon");
    e.throwMe();
  }
  tg_.setLatLon(lat,lon); 
}



//----------------------------------------------//
// setBeamTimeLatLonHeight(beamnum,t,lat,lon,ht)
// set beam time and look direction with latlon
// maintains correct order and resets targetgeom
// correctly                                    
//-----------------------------------------------//

void CassiniSim::setBeamTimeLatLonHeight(unsigned int beam_no, const Time& t, 
				   const Uvar& lat, const Uvar& lon,
				   const Uvar& ht)
{
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  setTime(t);
  // set Beam Number but don't make TargetGeom point to boresight
  // which is why setBeam() is not used.
  cur_beam_ = beam_no - 1;
  if(cur_beam_<0 || cur_beam_> 4){
    ErrorMessage e("CassiniSim: Bad beam number in setBeamTimeLatLon");
    e.throwMe();
  }
  tg_.setLatLonHeight(lat,lon,ht); 
}


//------------------------------------------------------------------//
// setBeamTimeAlongCross(beamnum,t,along,cross,ftrack)
// Set beam time and look direction with along,cross coordinates.
// Maintains correct order and resets targetgeom correctly.
// Track frame supplied is used to define along and cross-track coordinates
//-------------------------------------------------------------------//

void CassiniSim::setBeamTimeAlongCross(unsigned int beam_no, const Time& t, 
  const Uvar& along, const Uvar& cross, const Frame& ftrack)
  {
  if (!configured_)
    {
    ErrorMessage e("CassiniSim not configured");
    e.throwMe();
    }
  setTime(t);
  // set Beam Number but don't make TargetGeom point to boresight
  // which is why setBeam() is not used.
  cur_beam_ = beam_no - 1;
  if (cur_beam_ < 0 || cur_beam_ > 4)
    {
    ErrorMessage e("CassiniSim: Bad beam number in setBeamTimeLatLon");
    e.throwMe();
    }
  setTrackFrame(ftrack);
  tg_.setAlongCross(along,cross); 
  }

//-----------------------------------------------------------------------
// setBeamTimeAzimuthElevation(beamnum,t,azim,elev)
//
// Set Beam and time, and set look direction to the indicated offset
// in beam coordinates.  (Not the usual boresight setting).
// If beamnum and t are not included, then use the current beam and time.
// 
//-----------------------------------------------------------------------

void CassiniSim::setBeamTimeAzimuthElevation(unsigned int beam_no,
  const Time& t, const Uvar& azim, const Uvar& elev)
{
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  setTime(t);
  // set Beam Number but don't make TargetGeom point to boresight
  // which is why setBeam() is not used.
  cur_beam_ = beam_no - 1;
  if(cur_beam_<0 || cur_beam_> 4){
    ErrorMessage e("CassiniSim: Bad beam number in setBeamTimeLatLon");
    e.throwMe();
  }

  Frame bf(beam_frame_spice_id[cur_beam_],cassini_spice_id);

  // Set look vector
  DirectionVector look(bf,t,0,0,1);
  look.setAzimuthElevation(azim,elev);
  tg_.setLookDirection(look);
}

//----------------------------------------------------------------------
// setTime(t), reorient(t)
//
// Set the current time.  This time establishes all the SPICE related
// geometry.
//----------------------------------------------------------------------

void CassiniSim::setTime(const Time& t)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  t_=t;
  tg_.reset(t_);
  double et;
  t_.getEt(et);
  ftarget_.ephemeris(sc_state_,cassini_spice_id,et);
  tg_.setState(sc_state_);
  if (target_name_ != "none")
    {
    tg_.setTarget();
    }
  for (unsigned int i=1; i <= 5; ++i)
    {
    beam_[i-1].setBoresight(t_);
    }
  time_set_ = true;
  }

void CassiniSim::reorient(const Time& t)
  {
  if (!configured_)
    {
    ErrorMessage e("CassiniSim not configured");
    e.throwMe();
    }
  setTime(t);
  }

//-------------------------------------------------------------------
// setTrackFrame(ftrack)
//
// Set the track frame for the internal TargetGeom object.
// This method is called after setting the beam and time if
// the user wants to call any track-related methods.
//-------------------------------------------------------------------

void CassiniSim::setTrackFrame(const Frame& ftrack)
  {
  if (!configured_)
    {
    ErrorMessage e("CassiniSim not configured");
    e.throwMe();
    }
  tg_.setTrackFrame(ftrack);
  }

//---------------------------------------------------------------
// setTimeMismatch  set a mismatch between beam (ckernel) and
// ephemeris time
//---------------------------------------------------------------

void CassiniSim::setTimeMismatch(Uvar terr)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  for (unsigned int i=1; i <= 5; ++i)
    {
    beam_[i-1].setBoresight(t_,terr);
    }
  }

void CassiniSim::setBeam(unsigned int beamnum)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (beamnum < 1 || beamnum > 5)
    {
    ErrorMessage e("CassiniSim::setBeam: Error, beamnum=" + toStr(beamnum)
      + " out of range");
    e.throwMe();
    }

  // Set up current Beam index
  cur_beam_ = beamnum - 1;

  // Set beam related stuff
  DirectionVector boresight = beam_[cur_beam_].getBoresight();
  tg_.setLookDirection(boresight);
  }

//-------------------------------------------------------------------------
// setBeamTime(beamnum,t)
//
// Sets both current beam and current time.
// This method ensures that TargetGeom is handled appropriately unlike
// setBeam and setTime which can produce errors if setBeam is called more
// than once after setTime. 
//-------------------------------------------------------------------------

void CassiniSim::setBeamTime(unsigned int beamnum, const Time& t)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  setTime(t);
  setBeam(beamnum);
  }

TargetGeom CassiniSim::getTargetGeom(Uvar azim_frac, Uvar elev_frac)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1)
    {
    ErrorMessage
      e("CassiniSim::getTargetGeom: Error, No beam specified for CassiniSim");
    e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::getTargetGeom: Error, Time not set in CassiniSim");
    e.throwMe();
    }

  TargetGeom tg(t_);
  tg.setState(sc_state_);
  DirectionVector v = beam_[cur_beam_].beamFractionToDirection(azim_frac,elev_frac);
  tg.setLookDirection(v);
  tg.setTarget();
  return(tg);
  }

Uvar CassiniSim::getDoppler(Uvar azim_frac, Uvar elev_frac)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  TargetGeom tg = getTargetGeom(azim_frac,elev_frac);
  return(tg.doppler(lambda_));
  }

Uvar CassiniSim::getDoppler(Uvar lambda_chirp){
  return(tg_.doppler(lambda_chirp));
}
//-----------------------------------------------------------------------
// d = getDirectionTo(target)
//
// Return a direction vector pointing towards the indicated target or
// special direction at the current time.
// Special names:
//    SC_PX - positive x axis
//    SC_NX - negative x axis
//    SC_PY,SC_NY,SC_PZ,SC_NZ
//    SC_VEL - direction of velocity vector
//    SC_NVEL - negative of direction of velocity vector
//-----------------------------------------------------------------------

DirectionVector CassiniSim::getDirectionTo(const string& target)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  Frame fsc(spacecraft_ + "_SC_COORD", spacecraft_);

  if (target == "SC_PX")
    {
    return(DirectionVector("SC_PX",fsc,t_,1,0,0));
    }
  else if (target == "SC_NX")
    {
    return(DirectionVector("SC_NX",fsc,t_,-1,0,0));
    }
  else if (target == "SC_PY")
    {
    return(DirectionVector("SC_PY",fsc,t_,0,1,0));
    }
  else if (target == "SC_NY")
    {
    return(DirectionVector("SC_NY",fsc,t_,0,-1,0));
    }
  else if (target == "SC_PZ")
    {
    return(DirectionVector("SC_PZ",fsc,t_,0,0,1));
    }
  else if (target == "SC_NZ")
    {
    return(DirectionVector("SC_NZ",fsc,t_,0,0,-1));
    }
  else if (target == "SC_VEL")
    {
    return(DirectionVector("SC_VEL",sc_state_.velocity()));
    }
  else if (target == "SC_NVEL")
    {
    return(DirectionVector("SC_NVEL",-sc_state_.velocity()));
    }

  //------------------------------------------------------------
  // Not a special direction, so look for a SPICE body.
  // Light time correction is applied here, so we are assuming
  // passive observation of the body.
  //------------------------------------------------------------

  StateVector s;
  fsc.ephemeris(s,target,t_,"NONE");
  return(DirectionVector(target,s.position().representIn(fsc)));
  }


//--------------------------------
// Boresight related calculations
//--------------------------------

//-----------
// Boresight ground location in km
//----------

void CassiniSim::boresightPositionInKm(double pos[3]){
 if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightLat: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightLat: Error, Can't do this with target = none");
    e.throwMe();
    }
  PositionVector p=tg_.surfaceIntercept();
  pos[0]=get_in_base_units(p[PositionVector::X]);
  pos[1]=get_in_base_units(p[PositionVector::Y]);
  pos[2]=get_in_base_units(p[PositionVector::Z]);
}
//----------
// Latitude
//----------


Uvar CassiniSim::boresightLat()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightLat: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightLat: Error, Can't do this with target = none");
    e.throwMe();
    }

  return(tg_.lat());
  }

//------------
// Longitude
//------------

Uvar CassiniSim::boresightLon()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightLon: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightLon: Error, Can't do this with target = none");
    e.throwMe();
    }

  return(tg_.lon());
  }

//------------
// Slant range
//------------

Uvar CassiniSim::boresightRange()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightRange: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightRange: Error, Can't do this with target = none");
    e.throwMe();
    }

  return(tg_.range());
  }

//----------------
// Incidence Angle
//----------------

Uvar CassiniSim::boresightIncidenceAngle()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightIncidenceAngle: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightIncidenceAngle: Error, Not with target = none");
    e.throwMe();
    }

  return(tg_.incidenceAngle());
  }

//----------------
// Polarization Angle
//----------------

Uvar CassiniSim::boresightPolarizationAngle()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightPolarizationAngle: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightPolarizationAngle: Error, Not with target = none");
    e.throwMe();
    }

  return(tg_.polarizationAngle());
  }

//----------
// Doppler
//----------

Uvar CassiniSim::boresightDoppler()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightDoppler: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::boresightDoppler: Error, Not with target = none");
    e.throwMe();
    }

  return(tg_.doppler(lambda_));
  }

//---------------------------------------------------------
// Doppler Spread
// This method uses TargetGeom::dopplerSpreadBeam3 which
// assumes that the boresight is centered on the target.
// Currently only beam 3 is supported.
//---------------------------------------------------------

Uvar CassiniSim::beamDopplerSpread()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::beamDopplerSpread: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::beamDopplerSpread: Error, Not with target = none");
    e.throwMe();
    }

  if (cur_beam_ != 2)
    {
    ErrorMessage
      e("CassiniSim::beamDopplerSpread: Error, Only Beam 3 allowed");
    e.throwMe();
    }

  return(tg_.dopplerSpreadBeam3(lambda_,
    beam_[cur_beam_].getAzimuthWidthTwoWay()));
  }

//----------------------------
// Slew Rate in inertial space
//----------------------------

Uvar CassiniSim::boresightInertialSlewRate()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightSlewRate: Error, Beam and time not set");
    e.throwMe();
    }

  Uvar delta_t(1,"s");
  Frame fj2000("J2000",spacecraft_);
  DirectionVector b1 = beam_[cur_beam_].getBoresight(t_);
  DirectionVector b2 = beam_[cur_beam_].getBoresight(t_ + delta_t);
  b1.representIn(fj2000);
  b2.representIn(fj2000);
  return(b1.angle(b2)/delta_t);
  }

//--------------------------------
// Boresight speed on the surface
//--------------------------------

Uvar CassiniSim::boresightSpeed()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::boresightSlewRate: Error, Beam and time not set");
    e.throwMe();
    }

  return(tg_.surfaceVelocity().magnitude());
  }



//--------------------------------
// Area computations
//--------------------------------

Uvar CassiniSim::areaFootprintBoundingRectangle(
       const Time& t, int b)
{
  DebugInfo dbg("CassiniSim::areaFootprintBoundingRectangle");
  if(dbg.level){
    dbg.file << "CassiniSim::areaFootprintBoundingRectangle:"
	     << "  time:=" << t.utc("ISOD")
             << "  beam =" << b << endl;
  }
  Uvar target_radius=(default_target_radii.magnitude())/sqrt(3.0);

  // get widths in azimuth and elevation
  Uvar az=beam_[b-1].getAzimuthWidthTwoWay();
  Uvar el=beam_[b-1].getElevationWidthTwoWay();
  
  if(dbg.level){
    dbg.file << "  Two Way Beam Widths in degrees (azim,elev)=" << "("
	     << az.getInUnits("deg") << "," << el.getInUnits("deg")
	     << ")" << endl;      
  }
  // set up "corner points" for rectangle in beam frame
  PositionVector p[4];
  int p_idx=0;
  for(int a=0;a<2;a++){
    Uvar az2=a*az-az/2;
    for(int e=0;e<2;e++){
      Uvar el2=e*el-el/2;
      setBeamTimeAzimuthElevation(b,t,az2,el2);
      p[p_idx]=tg_.surfaceIntercept();
      if(dbg.level){
	dbg.file << "  Corner point " << p_idx+1 << " in TBF: "
		 << p[p_idx] << endl;
      }
      p_idx++;

    }
  }
  Uvar dummy1,dummy2,dummy3;
  Uvar area=getSphericalTriangleArea(p[0],p[1],p[2],target_radius,
				   dummy1,dummy2,dummy3);
  if(dbg.level){
    dbg.file << "  First triangle area =" << area;
  }
  area+=getSphericalTriangleArea(p[1],p[2],p[3],target_radius,
				   dummy1,dummy2,dummy3);
  if(dbg.level){
    dbg.file << "  Total area =" << area << endl << endl;
  }
  return(area);
} 

//---------------
// Nadir methods
//---------------

//------------------------------------------------
// a = altitude()
//
// Returns altitude above surface of target body.
//------------------------------------------------

Uvar CassiniSim::altitude()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::altitude: Error, time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::altitude: Error, Can't do this with target = none");
    e.throwMe();
    }

  return(tg_.altitude());
  }

//-------------------------------
// Nadir Latitude and Longitude
//-------------------------------

void CassiniSim::nadirLatLon(Uvar& lat, Uvar& lon)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::nadirLatLon: Error, time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::nadirLatLon: Error, Can't do this with target = none");
    e.throwMe();
    }

  tg_.nadirLatLon(lat,lon);
  }

//-----------------------------------------------
// Nadir Along track and Cross track coordinates 
//-----------------------------------------------

void CassiniSim::nadirAlongCross(Uvar& along, Uvar& cross)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::nadirAlongCross: Error, time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::nadirAlongCross: Error, Can't do this with target = none");
    e.throwMe();
    }

  tg_.nadirAlongCross(along,cross);
  }

//----------------
// Nadir Doppler
//----------------

Uvar CassiniSim::nadirDoppler()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::nadirDoppler: Error, time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::nadirDoppler: Error, Not with target = none");
    e.throwMe();
    }

  return(tg_.nadirDoppler(lambda_));
  }

//---------------
// Other methods
//---------------

Uvar CassiniSim::lambda()
  {
  if(!configured_)
    {
    ErrorMessage e("CassiniSim not configured");
    e.throwMe();
    }
  return(lambda_);
  }

//-------------------------------------------------------------------------
// a = subtendedAngle(target1,target2)
//
// Returns angle subtended between the direction vectors pointing from
// the spacecraft to each of the designated targets.
// Recognized target strings include SPICE recognized bodies, and
// special strings identifying spacecraft body axes and velocity vector.
// Special names:
//    SC_PX - positive x axis
//    SC_NX - negative x axis
//    SC_PY,SC_NY,SC_PZ,SC_NZ
//    SC_VEL - direction of velocity vector
//-------------------------------------------------------------------------

Uvar CassiniSim::subtendedAngle(const string& target1, const string& target2)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::subtendedAngle: Error, time not set");
    e.throwMe();
    }

  DirectionVector d1 = getDirectionTo(target1);
  DirectionVector d2 = getDirectionTo(target2);
  return(d1.angle(d2));
  }

//-----------------------------------------------------------------------
// getRADEC(target,ra,dec)
//
// Compute the RA and DEC of the specified target at the current time.
//-----------------------------------------------------------------------

void CassiniSim::getRADEC(const string& target, Uvar& ra, Uvar& dec)
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  Frame fj2000("J2000",spacecraft_);
  DirectionVector d = getDirectionTo(target);
  d.representIn(fj2000);
  d.getRADEC(ra,dec);
  }

//------------------------------------------------------------------------
// crossTrack2WayBeamWidth()
//
// Compute width of beam (2-way) perpendicular to the boresight
// track direction on the target surface.
// The track frame must be set first.
//------------------------------------------------------------------------

Uvar CassiniSim::crossTrack2WayBeamWidth()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::crossTrackBeamWidth: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::crossTrackBeamWidth: Error, Can't do this with no target");
    e.throwMe();
    }

  // Construct intermediate objects once
  static TargetGeom tg;
  static Uvar azim, elev;
  static Uvar along1, cross1;
  static Uvar along2, cross2;

  // Get along,cross track coordinates of upper elevation point
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  Frame beam_frame(beam_frame_spice_id[cur_beam_],cassini_spice_id);
  DirectionVector ulook(beam_frame,t_,0,0,0);
  azim = Uvar(0,"rad");
  elev = beam_[cur_beam_].getElevationWidthTwoWay();
  ulook.setAzimuthElevation(azim,0.5*elev);
  tg.setLookDirection(ulook);
  tg.interceptAlongCross(along1,cross1);

  // Get along,cross track coordinates of lower elevation point
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  ulook.setAzimuthElevation(azim,-0.5*elev);
  tg.setLookDirection(ulook);
  tg.interceptAlongCross(along2,cross2);

  return(fabs(cross2 - cross1));
  }

//------------------------------------------------------------------------
// crossTrack1Way5BeamWidth()
//
// Compute width of all 5 beams (1-way 3 dB) perpendicular to the boresight
// track direction on the target surface.
// The track frame must be set first.
//------------------------------------------------------------------------

Uvar CassiniSim::crossTrack1Way5BeamWidth()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::crossTrack1Way5BeamWidth: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::crossTrack1Way5BeamWidth: Error, Can't do this with no target");
    e.throwMe();
    }

  // Construct intermediate objects once
  static TargetGeom tg;
  static Uvar azim, elev;
  static Uvar along1a, along1b, cross1a, cross1b;
  static Uvar along5a, along5b, cross5a, cross5b;

  Frame beam_frame1(beam_frame_spice_id[0],cassini_spice_id);
  Frame beam_frame5(beam_frame_spice_id[4],cassini_spice_id);

  // Get along,cross track coordinates of beam 1 elevation points
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  DirectionVector ulook1(beam_frame1,t_,0,0,0);
  azim = Uvar(0,"rad");
  elev = beam_[0].getElevationWidthOneWay();
  ulook1.setAzimuthElevation(azim,0.5*elev);
  tg.setLookDirection(ulook1);
  tg.interceptAlongCross(along1a,cross1a);

  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  elev = beam_[0].getElevationWidthOneWay();
  ulook1.setAzimuthElevation(azim,-0.5*elev);
  tg.setLookDirection(ulook1);
  tg.interceptAlongCross(along1b,cross1b);

  // Get along,cross track coordinates of beam 5 elevation points
  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  DirectionVector ulook5(beam_frame5,t_,0,0,0);
  elev = beam_[4].getElevationWidthOneWay();
  ulook5.setAzimuthElevation(azim,0.5*elev);
  tg.setLookDirection(ulook5);
  tg.interceptAlongCross(along5a,cross5a);

  tg.reset(t_);
  tg.setTarget();
  tg.setState(sc_state_);
  tg.setTrackFrame(tg_.getTrackFrame());
  elev = beam_[4].getElevationWidthOneWay();
  ulook5.setAzimuthElevation(azim,-0.5*elev);
  tg.setLookDirection(ulook5);
  tg.interceptAlongCross(along5b,cross5b);

  // Return the maximum difference which is the cross track extent of all 5.
  if (cross1a > cross5a)
    {
    if (cross1a > cross1b)
      {
      return(cross1a - cross5b);
      }
    else
      {
      return(cross1b - cross5a);
      }
    }
  else
    {
    if (cross5a > cross5b)
      {
      return(cross5a - cross1b);
      }
    else
      {
      return(cross5b - cross1a);
      }
    }
  }

//------------------------------------------------------------------------
// alongTrackSpeed
//
// Compute speed of boresight point in the along track direction on
// the target surface.
//------------------------------------------------------------------------

Uvar CassiniSim::alongTrackSpeed()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (cur_beam_ == -1 || !time_set_)
    {
    ErrorMessage
      e("CassiniSim::alongTrackSpeed: Error, Beam and time not set");
    e.throwMe();
    }
  if (target_name_ == "none")
    {
    ErrorMessage
      e("CassiniSim::alongTrackSpeed: Error, Can't do this with no target");
    e.throwMe();
    }

  // Construct intermediate objects once
  static TargetGeom tg;
  static Uvar dt(0.1,"s");
  static StateVector dstate;
  static Uvar along1, cross1;
  static Uvar along2, cross2;
  static Uvar a_speed;

  // Get current boresight along,cross track coordinates
  tg_.interceptAlongCross(along1,cross1);

  // Get boresight point a small delta time away from now
  Time t(t_ + dt);
  tg.reset(t);
  double et;
  t.getEt(et);
  ftarget_.ephemeris(dstate,cassini_spice_id,et);
  tg.setState(dstate);
  tg.setTarget();
  tg.setTrackFrame(tg_.getTrackFrame());
  Frame beam_frame(beam_frame_spice_id[cur_beam_],cassini_spice_id);
  DirectionVector ulook(beam_frame,t,0,0,1);
  tg.setLookDirection(ulook);
  tg.interceptAlongCross(along2,cross2);

  // Along track speed of boresight point
  return(fabs(along2 - along1)/dt);
  }

//-------------------------------------------------------------------------
// a = radialDistance()
//
// Returns radialDistance for current TargetGeom.
//-------------------------------------------------------------------------

Uvar CassiniSim::radialDistance()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::radialDistance: Error, time not set");
    e.throwMe();
    }
  return(tg_.radialDistance());
  }

//-------------------------------------------------------------------------
// a = isoRangeRadiusAngle()
//
// Returns isoRangeRadiusAngle for current TargetGeom.
//-------------------------------------------------------------------------

Uvar CassiniSim::isoRangeRadiusAngle()
  {
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  if (!time_set_)
    {
    ErrorMessage
      e("CassiniSim::radialDistance: Error, time not set");
    e.throwMe();
    }
  return(tg_.isoRangeRadiusAngle());
  }

void CassiniSim::adjustSpacecraftAltitude(Uvar offset)
{
  // do nothing
}

void CassiniSim::adjustSpacecraftCrossTrackPosition(Uvar offset)
{
  // do nothing
}



