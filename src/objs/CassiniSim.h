//============================================================================
// CassiniSim.h
//
// This file contains the CassiniSim class declaration.
// The CassiniSim class provides spacecraft related calculations and
// data.  It is based fundamentally on time like most of the RAS classes.
// This means that calculations are performed at a particular time with the
// SPICE system specifying the positions,velocities of various objects at
// that time.
// In the future, an array of times may be supported.
//=============================================================================

#ifndef CASSINI_SIM_H
#define CASSINI_SIM_H

#include "Frame.h"
#include "Beam.h"
#include "Units.h"
#include "Config.h"
#include "Time.h"
#include "TargetGeom.h"
#include "Array.h"
#include "SARFunctions.h"
#include <vector>

class CassiniSim
  {
  public:

  //---------------
  // Construction
  //---------------

  CassiniSim();
  CassiniSim(Config& cfg);
  CassiniSim(unsigned int beamnum, Config& cfg);

  //---------------
  // Configuration
  //---------------

  void config(Config& cfg);

  //------------------
  // Get/Set methods
  //------------------

  // reorient reorients frames to projected position at time t
  void setTime(const Time& t);  // same as reorient
  void reorient(const Time& t);

  void setTrackFrame(const Frame& ftrack);

  // setTimeMismatch  set a mismatch between beam (ckernel) and
  // ephemeris time
  void setTimeMismatch(Uvar terr);
  void setBeam(unsigned int beamnum);
  void setBeamTime(unsigned int beamnum, const Time& t);
  void setMidPulseImpactTime( const Time& pulse_transmit_time,
			      const Uvar& lat,  const Uvar& lon);
  void setBeamTimeLatLon(unsigned int beam_no, 
			 const Time& t, const Uvar& lat, const Uvar& lon);
  void setBeamTimeLatLonHeight(unsigned int beam_no, 
			       const Time& t, const Uvar& lat, const Uvar& lon,                                const Uvar& ht);
  void setBeamTimeAzimuthElevation(unsigned int beam_no, 
			 const Time& t, const Uvar& azim, const Uvar& elev);
  void setBeamTimeAlongCross(unsigned int beam_no,
    const Time& t, const Uvar& along, const Uvar& cross, const Frame& ftrack);
  // getTargetGeom returns a target geom object given an azimuth and elevation
  // fraction
  TargetGeom getTargetGeom(Uvar azim_frac, Uvar elev_frac);
  void getDopplerProfile(const Uvar& lambda_chirp,
			 Uvar& slope, Uvar& doppler_boresight,
			 Uvar& range_boresight,Uvar& rate_slope,
			 Uvar& doppler_rate_boresight);

  void getDopplerCentroid(const Uvar& lambda_chirp,
			       const Uvar& range, PositionVector& pos,
			  Uvar& fdc);

  void getFixedTargetDopplerAtTime(const Uvar& lambda_chirp, 
				   const PositionVector& surface_pos, 
				   const Time& t, Uvar& fdop);

  //---------------------------------
  // Boresight related calculations
  //---------------------------------

  void boresightPositionInKm(double pos[3]);
  Uvar boresightLat();
  Uvar boresightLon();
  Uvar boresightRange();
  Uvar boresightIncidenceAngle();
  Uvar boresightPolarizationAngle();
  Uvar boresightDoppler();
  Uvar beamDopplerSpread();
  Uvar boresightInertialSlewRate();
  Uvar boresightSpeed();
  

  //--------------------------------
  // Area computations
  //--------------------------------

  Uvar areaFootprintBoundingRectangle(const Time& t, int b);

  //---------------------------------
  // Nadir related calculations
  //---------------------------------

  Uvar altitude();
  void nadirLatLon(Uvar& lat, Uvar& lon);
  void nadirAlongCross(Uvar& along, Uvar& cross);
  Uvar nadirDoppler();
  
  //---------------
  // Other methods
  //---------------

  Uvar lambda();
  Uvar getDoppler(Uvar lambda_chirp);
  Uvar getDoppler(Uvar azim_frac, Uvar elev_frac);
  DirectionVector getDirectionTo(const string& target);
  DirectionVector getDirectionTo(const Uvar& lat, const Uvar& lon);
  Uvar subtendedAngle(const string& target1, const string& target2);
  void getRADEC(const string& target, Uvar& ra, Uvar& dec);
  Uvar crossTrack2WayBeamWidth();
  Uvar crossTrack1Way5BeamWidth();
  Uvar alongTrackSpeed();
  Uvar radialDistance();
  Uvar isoRangeRadiusAngle();

  // adjustSpacecraftAltitude 
  void adjustSpacecraftAltitude(Uvar offset);

  // adjustSpacecraftCrossTrackPosition
  void adjustSpacecraftCrossTrackPosition(Uvar offset);

  protected:

  //-------------------------- 
  // Internal representation
  //-------------------------- 

  string spacecraft_;
  string target_name_;
  Frame ftarget_;
  Time t_;
  bool time_set_;
  bool configured_;
  StateVector sc_state_;
  vector<Beam> beam_;
  int cur_beam_;
  Uvar lambda_;
  TargetGeom tg_;
};

#endif



