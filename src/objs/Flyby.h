//============================================================================
// Flyby.h
//
// This file contains the declaration for classes Flyby, Orbit, and ObjectData.
// Flyby provides support for target flyby calculations.
// In particular, it contains named epoch times, and the ability to locate
// closest approach times.
// Orbit provides orbit propagation and parameterization services.
// ObjectData provides fixed data on various planetary objects.  It pulls
// most of its data from the SPICE system.
//============================================================================

#ifndef Flyby_H
#define Flyby_H

#include <string>
#include <map>
#include "Units.h"
#include "Error.h"
#include "Array.h"
#include "TargetGeom.h"

//-------------------------------
// Flyby
// Timing support for flyby's
//-------------------------------

class Flyby
  {
  public:

  //------------------
  // Automatic testing
  //------------------

  static bool selfTest();

  //--------------
  // Constructors
  //--------------

  Flyby(Config& cfg);

  //--------------
  // Setup
  //--------------

  void config(Config& cfg);
  void setTarget(const string& target);
  void setObserver(const string& observer);
  void setEpoch(const string& epoch_selection, const Time& t,
    const Uvar& accuracy);

  //--------------
  // Predicates
  //--------------

  bool validSuffix(const string& suffix) const;
  bool done() const;
  bool done(const string& suffix);

  //------------------------
  // Special service methods
  //------------------------
  void suffixCheck(const string& suffix) const;
  Uvar startTime() const;
  Uvar endTime() const;
  Uvar interval() const;
  Uvar startTime(const string& suffix) const;
  Uvar endTime(const string& suffix) const;
  Uvar interval(const string& suffix) const;
  Time epochTime() const;
  Uvar epochAccuracy() const;
  unsigned int numTimeSteps() const;
  unsigned int numTimeSteps(const string& suffix) const;
  string target() const;
  string id() const;
  Uvar lowestAltitude();
  Time lowestAltitudeTime();
  StateVector lowestAltitudeState();
  string LookDirection();
  Uvar timeOffset(const Time& t);
  string sc_Orientation(); 
  string trackingOption();
  Uvar azimuthRotation();
  Uvar sc_z_rotation();

 

  //--------------------------
  // Ask flyby to compute polyfit coeff for angle profile and target_lat_lon
  //---------------------------
  void computePolyfitCoeff();
  void load_local_incidence_angle_offset_params(const vector<Uvar> local_incidence_angle_offset,
						const vector<Uvar> local_incidence_angle_offset_tr1,
						const vector<Uvar> local_incidence_angle_offset_tr2,
						const vector<Uvar> local_incidence_angle_offset_tr3,
						const vector<Uvar> local_incidence_angle_offset_tr4);
 

  //------------------------------------
  // Ask incidence angle at time t
  // for polyfit, computePolyCoeff() should be performed
  //-------------------------------------
  Uvar desiredIncidenceAngle(const Time& t);
  Uvar desiredIncidenceAngle_with_local_bias(const Time& t);
  //---------------------------
  // trackFrame related method: Needed to be done just once for a flyby
  //---------------------------
  Frame trackFrame();//return computed track frame
  void computeTrackFrame();//compute track frame if not computed before

  //---------------------
  // Time stepping
  //---------------------
  void resetTimeSteps();
  Time currentTime();
  void percentCompletedDisplay(ostream& s, unsigned int percent_step);
  void resetTimeSteps(const string& suffix);
  Time currentTime(const string& suffix);
  void percentCompletedDisplay(ostream& s, unsigned int percent_step,
    const string& suffix);

  

  private:

  //-------------------------
  // Internal support methods
  //-------------------------
  Uvar timeValueSubstitute(Config& cfg, const string& keyword);
  
  //------------------------
  // Internal Representation
  //------------------------
  string cfg_filename_;
  unsigned int time_counter_;
  map<string,unsigned int> time_counter_map_;
  Time cur_time_;
  map<string,Time> cur_time_map_;
  unsigned int percent_completed_;
  map<string,unsigned int> percent_completed_map_;
  Uvar start_time_;
  Uvar end_time_;
  Time epoch_time_;
  Uvar interval_;
  Uvar ivd_time_pad_;
  string target_;
  string spacecraft_;
  string sc_orientation_;
  string tracking_option_;
  Uvar azimuth_rotation_;
  Uvar sc_z_rotation_;
  bool lowest_altitude_set_;
  bool track_frame_set_;
  bool computePolyfitCoeff_set_;
  bool local_offset_set_;
  Uvar lowest_altitude_;
  Time lowest_altitude_time_;
  unsigned int ntime_;
  mutable map<string,unsigned int> ntime_map_;
  string epoch_selection_;
  Uvar epoch_accuracy_;
  Time ck_earliest_start_time_;
  Time ck_latest_end_time_;
  Frame track_frame_;
  unsigned int Noffset_;
  mutable Config::numbers start_time_map_;
  mutable Config::numbers end_time_map_;
  mutable Config::numbers interval_map_;
  
 
  //parameters needed for  polynomial fit of incidence angle
  string angle_profile_option_;
  Uvec altitude_polyfit_;
  Uvec angle_polyfit_;
  Dvec angle_poly_coeff_time_;
  Uvar maximum_polyfit_altitude_;
  Uvar angle_at_maximum_polyfit_altitude_;
  
 
  //parameters needed by sech-function type incidence angle fit
  Uvar angle_at_ca_;
  Uvar angle_at_mid_alt_;
  Uvar mid_alt_;
  Uvar residual_angle_;

  //parameter needed by special lat-lon tracking option
  string target_lat_lon_option_;
  Uvar off_centered_beam_time_lag_;
  Uvar start_target_lat_lon_tracking_;
  Uvar end_target_lat_lon_tracking_;
  Uvar angle_at_closest_range_to_target_;
  Uvar time_at_closest_range_to_target_;
  Dvec target_lat_lon_fit_;
  Uvec angle_offset_;
  Uvec tr1_,tr2_,tr3_,tr4_;
  Uvec tr12_mid_,tr34_mid_;
  Dmat tr_from_ideal_;
  Dmat tr_to_ideal_;
  };

//--------------------------------------------------
// Orbit
// Trajectory propagation and parameter calculation
//--------------------------------------------------

class Orbit
  {
  public:

  typedef list<StateVector> STATELIST;

  //------------------
  // Automatic testing
  //------------------

  static bool selfTest();

  //--------------
  // Constructors
  //--------------

  Orbit(const string& target_name, const StateVector& state);

  //--------------
  // Setup
  //--------------

  void setTarget(const string& target_name);
  void setState(const StateVector& state);

  //---------------------
  // Orbit Propagation
  //---------------------

  StateVector getState(const Time& t);

  //----------------------
  // Orbit Characteristics
  //----------------------

  Frame bPlaneFrame();
  StateVector bPlaneIntercept();
  FloatVector asymptoticApproachVelocity();
  Uvar lowestAltitude();
  Uvar lowestAltitudeTime();
  StateVector lowestAltitudeState();

  private:

  //------------------------
  // Internal Representation
  //------------------------

  bool target_set_;
  bool lowest_altitude_set_;
  Uvar lowest_altitude_;
  Time lowest_altitude_time_;
  static const Uvar interval_;
  string target_name_;
  Frame ftarget_I_;
  Uvar GM_;
  STATELIST states_;

  };

//-----------------------------------------------
// Helper function object used by Flyby and Orbit
//-----------------------------------------------

class TargetAltitude
  {
  public:

  // Construct using observer name and pull state from SPICE ephemeris
  TargetAltitude(const string& observer_name, const string& target_name)
    : observer_name_(observer_name), porbit_(NULL), target_name_(target_name)
    { }

  // Construct from user supplied Orbit (ptr)
  TargetAltitude(Orbit* porbit, const string& target_name)
    : observer_name_(""), porbit_(porbit), target_name_(target_name)
    { }

  // Altitude at time t
  Uvar operator()(const Uvar& t)
    {
    TargetGeom tg(t);
    tg.setTarget(target_name_);
    if (porbit_ != NULL)
      {
      tg.setState(porbit_->getState(t));
      }
    else
      {
      tg.setState(observer_name_);
      }
    return(tg.altitude());
    }

  private:

  string observer_name_;
  Orbit* porbit_;
  string target_name_;
  };



class TargetTimeAltitude
  {
  //return time when the sc is at the target altitude
  // since titan flyby is hyperbolic, it is the users responsibility to set 
  // proper time window in order to have a single solution for a given altitude
  //

  public:

  // Construct using observer name and pull state from SPICE ephemeris
  TargetTimeAltitude(const string& observer_name, const string& target_name,const Uvar& target_altitude)
    : observer_name_(observer_name),  target_name_(target_name),target_altitude_(target_altitude)
    { 
    ftitan_ = Frame("IAU_TITAN",target_name_);
    }

  

  // Altitude difference between target altitude, and altitude at time t
  Uvar operator()(const Uvar& t)
    {
    ftitan_.ephemeris(sc_state_,observer_name_,t,"NONE");
    TargetGeom tg(t);
    tg.setState(sc_state_);
    tg.setTarget(target_name_,ftitan_);
    Uvar diff = tg.altitude() - target_altitude_;
    if (diff < Uvar(0,"km")) diff = -diff;
    return(diff);
    }

  private:
  string observer_name_;
  Frame ftitan_;
  StateVector sc_state_;
  string target_name_;
  Uvar target_altitude_;
  };




class TargetLatLon
  {
  public:

  // Construct using observer name and pull state from SPICE ephemeris
    TargetLatLon(const string& observer_name, const string& target_name,
		 const Time& t, const Uvar& lat, const Uvar& lon, 
		 const string& look)
    :observer_name_(observer_name), target_name_(target_name),
    t_(t),lat_(lat),lon_(lon),look_(look)
    { 
      ftitan_ = Frame("IAU_TITAN",target_name_);
      ftitan_.ephemeris(sc_state_,observer_name_,t_,"NONE");
      dir_to_surface_=DirectionVector("dir to target",ftitan_,t_,0,0,1);
      dir_to_surface_.setPlanetodetic(lat_,lon_);
    }

  
 
  // Surface distance for a given incidence angle
  Uvar operator()(const Uvar& incidenceAngle)
    {
    TargetGeom tg(t_);
    tg.setState(sc_state_);
    tg.setTarget(target_name_,ftitan_);
    DirectionVector x_dir,y_dir,z_dir;
    tg.setIsodopplerIncidence(incidenceAngle,
			      look_,
			      Uvar(0,"rad"),
			      x_dir,
			      y_dir,
			      z_dir);
    PositionVector p = tg.surfaceIntercept();
    surface_= dir_to_surface_*tg.radius();
    return( tg.radius()
	    *fabs(p.angle(surface_).getInUnits("rad")));
    }

  private:
  string observer_name_;
  string target_name_;
  Time t_;
  Uvar lat_;
  Uvar lon_;
  string look_;
  StateVector sc_state_;
  Frame ftitan_;
  PositionVector surface_;
  DirectionVector dir_to_surface_;
  };



class MinimizeGap
  {
  public:

    //start with alt change
    MinimizeGap(const Dvec alt,
		const Dvec angle,
		const Dvec& alt_polyfit, 
		const double& fuzzy)
      :alt_polyfit_("alt_polyfit"),alt_("alt"),angle_("angle")
      {
	N_ = alt.size();
	alt_.resize(N_);
	angle_.resize(N_);
	alt_polyfit_.resize(alt_polyfit.size());
	
	alt_ = alt;
	angle_ = angle;
	alt_polyfit_ = alt_polyfit;
	
	fuzzy_ = fuzzy;
      }

    double operator() (const Dvec& ptry)
      {
	//from trajectory
	// alt = alt_polypolyfit(i) * time_in_min^(i)

	double angle_Error = 0.0;
	for (unsigned int i = 0; i < N_;++i)
	  {
	    double sum = 0.0;
	    for (unsigned int j = 0; j < ptry.size();++j)
	      {
		sum += ptry(j) * pow(alt_(i),j);
	      }
	    angle_Error += (angle_(i) - sum)*(angle_(i) - sum);
	  }
	
	//derivative
	double der_Error=0;
	for (unsigned int i = 0; i < N_;++i)
	  {
	    double dtheta_dh = 0.0;
	    double dh_dt=0;
	    for (unsigned int j = 1; j < ptry.size();++j)
	      {
		dtheta_dh += ptry(j) *double(j)* pow(alt_(i),j-1);
	      }
	    for (unsigned int j = 1; j <alt_polyfit_.size();++j)
	      {
		dh_dt += alt_polyfit_(j) * double(j) * pow(double(i),j-1);
	      } 
	    der_Error += dtheta_dh *dtheta_dh* dh_dt*dh_dt;
	  }
	return(angle_Error+fuzzy_*der_Error);
      }
    
  private:
    Dvec alt_polyfit_;
    Dvec alt_;
    Dvec angle_;
    unsigned int N_;
    double fuzzy_;
  };

#endif
