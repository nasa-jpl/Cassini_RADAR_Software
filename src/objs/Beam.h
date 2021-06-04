//-----------------------------------------------------------------------------
// Beam.h
//
// This file contains the Beam class declaration.
// The Beam class provides an interface to obtain an antenna beam pattern 
// from a file or to  generate a beam pattern based on beamwiths of azi 
// and elev angles.
// The Beam class also implements low level error handling facilities
// described in the comments in Beam.cpp.
//-----------------------------------------------------------------------------

#ifndef Beam_H
#define Beam_H

#include <string>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Config.h"

//----------------------------
//Forward declaration
//----------------------------

//class Beam;

using std::string;

//-------------------------
// Class Beam declaration
//-------------------------

class Beam
  {
  public:
 
  enum BeamExtremaType {CENTER, TOP, BOTTOM, LEFT, RIGHT};

  //======================================================
  // Interface routines that don't need setBoresight first
  //======================================================

  //--------------
  // Construction
  //-------------- 

  // Dummy constructor
  Beam(); 
  // Constructor from Config file
  Beam(unsigned int beamnum, Config& cfg);
  // Constructor from Config file with user over-ride
  Beam(unsigned int beamnum, Config& cfg, bool read_pattern);

  // configure Beam object from Config object
  void config(Config& cfg, unsigned int beamnum);
  // configure Beam object from Config object with user over-ride
  void config(Config& cfg, unsigned int beamnum, bool read_pattern);

  //----------------
  // Get/Set methods
  //----------------

  // Antenna gain methods
  double getGainOneWay(DirectionVector u);
  // return peak one-way gain
  double getMaxGain();

  // Boresight
  DirectionVector getBoresight(const Time& t) const;

  // Pattern parameters
  unsigned int getBeamNumber() const throw(ErrorMessage); 
  unsigned int getNazim(){return(Nazim_);}
  unsigned int getNelev(){return(Nelev_);}
  Uvar getAzimStep(){return(azim_step_antenna_);}
  Uvar getElevStep(){return(elev_step_antenna_);}
  Uvar getAzimStart(){return(azim_start_antenna_);}
  Uvar getElevStart(){return(elev_start_antenna_);}

  // return beam widths
  Uvar getAzimuthWidthOneWay();
  Uvar getElevationWidthOneWay();
  Uvar getAzimuthWidthTwoWay();
  Uvar getElevationWidthTwoWay();

  // beam frame
  const Frame&  getFrame(Time& t) const throw(ErrorMessage);

  //---------------------------
  // Other supporting functions
  //---------------------------

  void estimateBeamWidthsFromPattern() throw (ErrorMessage);
  // load a beam pattern from an external file
  void loadBeamPatternFromFile();
  // compute an artificial beam pattern using sinc function
  void computeSincBeamPattern();

 
  
  double getAzimuthIndex(const Uvar& azim) const; 
  double getElevationIndex(const Uvar& elev) const; 
 
  //integrate beam pattern
  void integrateBeamPattern(double& solid_angle, double& directivity_dB, double& beam_integral);
  

  double bilinear(const Uvar& azim, const Uvar& elev);
  double bilinear(double azim_in_rad, double elev_in_rad);

  //======================================================================
  // Special Interface routines:
  // These methods require setBoresight to be called first.  This sets a
  // time inside the Beam object.  The time
  // is then used to support some special operations where two times
  // are involved.
  //======================================================================

  // Constructor from Config operator and beam number and time
  // Does same as previous plus sets boresight in the approariate 
  // direction for time t.
  Beam(unsigned int beamnum, Config& cfg, bool read_pattern, Time t);
 
  //------------------------------------------------------------------------
  // Set boresight direction vector appropriately for time t.
  // If terr is assigned it simulates the effect on the boresight of
  // a time mismatch error so that the beam frame location is what it
  // should be at time t but the boresight look direction (in the inertial
  // frame is what it would be at time t+terr).
  // Returns the boresight as it was set.
  //------------------------------------------------------------------------

  DirectionVector setBoresight(const Time& t, const Uvar& terr) 
    throw(ErrorMessage);
  DirectionVector setBoresight(const Time& t) 
    throw(ErrorMessage);

  // Returns current  boresight
  DirectionVector getBoresight() const;

  const Frame&  getNominalFrame() const throw(ErrorMessage);

  // Returns direction vector for extrema locations
  // For example:
  // etype=TOP yields direction vector at positive elevation zero azimuth
  //           3 dB two-way point.
  DirectionVector getExtrema(BeamExtremaType etype);

  //---------------------------
  // Other supporting functions
  //---------------------------

  // Get direction vector from two-way beam width fractions
  DirectionVector beamFractionToDirection(Uvar azim_frac, Uvar elev_frac);

  // Get Two-Way beam fractions from direction vector
  Uvec directionToBeamFraction(DirectionVector v);
  
  // Determine location from gain
  Uvec azimAndGainToElev(const Uvar& azim, double gain_in_dB);
  Uvec elevAndGainToAzim(const Uvar& elev, double gain_in_dB); 

  // For a given elevation angle, locate azimuth angle for max gain
  // using geometry
  Uvar elevAndMaxGainToAzim(const Uvar& elev);

  // For a given elevation angle, locate azimuth angle for max gain
  // using pre-computed max gain line
  Uvar maxGainLineToAzim(const Uvar& elev);

  // Calculate max gain line in elev and azimuth
  void computeMaxGainLine();

  //compute best fit ellipse for a gain contour
  void computeBestFitEllipse(Uvec& azimuth, Uvec& elevation, const double& gain_db);
 
  //Get one-way (azimuth,elevation) points for a specific gain 
  void get1wayBeamGain(Uvec& azimuth, Uvec& elevation, const double& gain_db);


 
  // Not yet implemented
  unsigned int findex(const Uvar& target_value,const Uvec& coarse_array);


 
  private:

  //constant class data -- Hardcodes size of beam pattern array
  // unfortunately I can't get const or static to work for Uvar
  Uvar azim_step_antenna_;
  Uvar elev_step_antenna_;
  Uvar azim_start_antenna_;
  Uvar elev_start_antenna_;

  // double version of start and step information
  double azim_step_antenna_in_rad_;
  double elev_step_antenna_in_rad_;
  double azim_start_antenna_in_rad_;
  double elev_start_antenna_in_rad_;

  Uvar max_gain_line_azim_offset_;
  Uvar max_gain_line_slope_;
  const static unsigned int Nazim_;
  const static unsigned int Nelev_;
  // constant value used to compute beam pattern from sinc function
  const static double beam_broadening_factor_;  
  

  // Booleans for determining to what degree the beam object has been
  // configured
  bool compute_beam_pattern_;
  bool read_beam_pattern_;
  bool beam_parameter_set_;
  bool beam_widths_computed_;
  bool beam_num_set_;
  bool boresight_set_;
  bool pattern_loaded_;

  // location vector and frame for electrical boresight of beam
  DirectionVector boresight_;
  Frame nominal_frame_; // nominal beam frame
  Frame frame_; // actual beam frame (static relative to J2000)
  Time t_; // time for which frame_ is valid in special interface

  // parameters needed to create artificial (sinc) beam pattern
  // (beam_number_ is also used to determine which pattern file to read)
  // Reading in a beam pattern sets the beam width variables accordingly
  // widths are (re)calculated after gain pattern is computed or read in
  unsigned int beam_number_;
  Uvar beam_width_azim_oneway_;  
  Uvar beam_width_elev_oneway_;
  Uvar beam_width_azim_twoway_;  
  Uvar beam_width_elev_twoway_;

  //----------------------------------------
  // BEAM relative one-way beam pattern array 
  // maximum value is 1.
  //-----------------------------------------
  Dmat beam_gain_; 
  double beam_max_gain_; // peak one-way gain
  double beam_min_gain_; // minimum gain

  // beam pattern file name
  string pattern_filename_;

  string spacecraft_;
  
  };


#endif






