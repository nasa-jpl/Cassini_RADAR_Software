//==============================================================================
// TargetGeom.h
//
// This file contains the class TargetGeom and associated functions that
// provide geomatric calculation support.
// This header comment summarizes the interface.
// For details about a specific function, look at the declarations in this file
// or the function leader comment in the .cpp file.
//
// The NAIF SPICE toolkit forms the foundation of the automatic frame handling
// system, therefore the appropriate SPICE kernel files need to be loaded
// before any frame conversions can be handled.  A static member function
// called Frame::spiceLoad is provided to load kernel files.
//
// Interface summary:
//
// class TargetGeom;
//   Geometric calculation helper.  This class computes many of the needed
//   geometric quantities when sufficient input data is provided.
//
//
// **************************************************************************
//
// TargetGeom Methods:
//
// Construction:
//
//   TargetGeom(t);    Setup for geometry calculations at time t
//
// Other methods:
//
//
//
//   
//==============================================================================

#ifndef TargetGeom_H
#define TargetGeom_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Frame.h"

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//----------------------
// Forward declarations
//----------------------

using std::string;

//-----------------------------
// Class TargetGeom declaration
//-----------------------------

class TargetGeom
  {
  public:

  enum indexE {X, Y, Z};

  //--------------
  // Construction
  //--------------

  TargetGeom();
  TargetGeom(const Time& t);

  //---------
  // Testing
  //---------

  static bool selfTest() ; // No exceptions

  //-----------
  // Reset
  //-----------

  void reset(const Time& t) ; // No exceptions

  //-----------
  // Predicates
  //-----------

  bool foundSurfaceIntercept() ; // No exceptions

  //-------------------------------------------------
  // Primary inputs - all of these must be specified
  //-------------------------------------------------

  void setState(const StateVector& state) ;
  void setState(const string& observer_name) ;

  void setTarget(); // uses defaults set by Frame::config
  void setTarget(const string& target_name) ;
  void setTarget(const string& target_name, const Frame& target_frame)
    ;

  //---------------------------------------------------
  // Secondary inputs - one of these must be specified
  //---------------------------------------------------

  void setLatLon(const Uvar& lat, const Uvar& lon) ;
  void setLatLonHeight(const Uvar& lat, const Uvar& lon, const Uvar& ht) ;
  void setLookDirection(const DirectionVector& look) ;
  void setAzimuthIncidence(const Uvar& azimuth, 
			   const Uvar& thetai,
			   const string& side,
			   DirectionVector& x, 
			   DirectionVector& y, 
			   DirectionVector& z);
  void setIsodopplerIncidence(const Uvar& thetai, 
			      const string& side,
			      const Uvar& set_sc_z_rotation,
			      DirectionVector& x, 
			      DirectionVector& y,
			      DirectionVector& z) ;

  void setPitchBiasedIncidence(const Uvar& thetai, 
			       const string& side,
			       const Uvar& set_sc_z_rotation,
			       const Uvar& pitch,
			       DirectionVector& x, 
			       DirectionVector& y,
			       DirectionVector& z) ;

  void setContiguousSurfaceCoverage(const Uvar& thetai, 
				    const string& side,
				    DirectionVector& x, 
				    DirectionVector& y,
				    DirectionVector& z) ;

  void setRingRadiusSweepClosest
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
    );
 
  void setRangeDopplerInBeam3Frame(const Uvar& range, const Uvar& doppler, 
    const Uvar& lambda);
  void setRangeDopplerInTargetFrame(const Uvar& range, const Uvar& doppler, 
		 const Uvar& lambda);
  void setRangeAzimuth(const Uvar& range, const Uvar& azi0, const Frame& fbeam);
  void setRangeElevation(const Uvar& range, const Uvar& ele0, const Frame& fbeam);
  void setRangeAzimuth(const Uvar& range, const Uvar& azi0, const unsigned int& beam_number);
  void setRangeElevation(const Uvar& range, const Uvar& ele0, const unsigned int& beam_number);
  void setAlongCross(const Uvar& alongtrack, const Uvar& crosstrack);

  //---------------------------------------------------------------
  //Special method to define alongtrack/crosstrack in a fixed frame
  //To calculate/access alongtrack/crosstrack, ftrack must be
  //set in advance. 
  //--------------------------------------------------------------

  void setTrackFrame(const Frame& ftrack);

  //----------------
  // Output methods
  //----------------

  Frame	          getTrackFrame();
  Uvar            lat() ;
  Uvar            lon() ;
  DirectionVector normal() ;
  PositionVector  surfaceIntercept() ;
  PositionVector  mirrorsurfaceIntercept();
  FloatVector     surfaceVelocity();
  Uvar            range() ;
  PositionVector  lookVector() ;
  DirectionVector lookDirection() ;
  DirectionVector mirrorlookDirection();
  Uvar            incidenceAngle() ;
  Uvar            mirrorincidenceAngle();
  Uvar            polarizationAngle() ;
  Uvar            doppler(const Uvar& lambda);
  PositionVector  nadir() ;
  StateVector     state() {return state_;}
  Uvar            altitude() ;
  void	          nadirLatLon(Uvar& lat, Uvar& lon);
  void	          nadirAlongCross(Uvar& along, Uvar& cross);
  Uvar            nadirDoppler(const Uvar& lambda);
  Uvar            limbRoll();
  Uvar            limbRange();
  void            interceptAlongCross(Uvar& alongtrack, Uvar& crosstrack);
  Uvar            radius();
  Uvar            dopplerSpreadBeam3(const Uvar& lambda,
				     const Uvar& beam_width);
  Uvar radialDistance();
  Uvar isoRangeRadiusAngle();
  
  private:

  //-------------------------
  // Internal representation
  //-------------------------

  // Flags indicating precomputed or set data available
  bool state_set_;
  bool target_set_;
  bool ulook_set_;
  bool mirrorlook_set_;
  bool lat_set_;
  bool lon_set_;
  bool unormal_set_;
  bool surface_intercept_set_;
  bool mirror_surface_intercept_set_;
  bool range_set_;
  bool rlook_set_;
  bool thetai_set_;
  bool mirror_thetai_set_;
  bool polang_set_;
  bool doppler_set_;
  bool rnadir_set_;
  bool along_cross_set_;

  // Flags indicating computation call in progress.
  // These are used to prevent infinite recursion.
  bool ulook_call_;
  bool lat_call_;
  bool lon_call_;
  bool unormal_call_;
  bool surface_intercept_call_;
  bool range_call_;
  bool rlook_call_;
  bool thetai_call_;
  bool polang_call_;
  bool doppler_call_;
  bool along_cross_call_;

  bool spice_found_;
  bool time_set_;
  bool track_frame_set_;
  bool ring_target_;

  Time t_;
  StateVector state_;
  DirectionVector ulook_;
  DirectionVector mirrorlook_;
  PositionVector surface_intercept_;
  PositionVector mirror_surface_intercept_;
  Frame target_frame_;
  Frame track_frame_;
  PositionVector target_radii_;
  DirectionVector unormal_;
  PositionVector rlook_;
  PositionVector rnadir_;
  string target_name_;
  SpiceInt target_id_;
  Uvar lat_;
  Uvar lon_;
  Uvar range_;
  Uvar thetai_;
  Uvar mirror_thetai_;
  Uvar polang_;
  Uvar doppler_;
  Uvar altitude_;
  Uvar along_;
  Uvar cross_;
  
  };

class TargetGeomError : public ErrorMessage
{
 public:
  enum errorE {unspecified, beyond_horizon};
  TargetGeomError(const string& emsg, errorE err_type = unspecified); 
  TargetGeomError(errorE err_type = unspecified); 
  void throwMe();
  string msg;
  errorE error_type;
};
#endif






