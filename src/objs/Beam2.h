#ifndef BEAM2_H
#define BEAM2_H
#include"Frame.h"
#include"Config.h"
#include"Units.h"
#include"Time.h"

typedef UnitVar<double> Uvar;
typedef Array1D<Uvar> Uvec;

class Beam2{

 public:

  enum BeamExtremaType {CENTER, TOP, BOTTOM, LEFT, RIGHT};
  // Dummy constructor
  Beam2();

  // Constructor using beamnumber and Config Object
  Beam2(unsigned int beamnum, const Config& cfg);

  // Constructor using beamnumber, Config object and time
  Beam2(unsigned int beamnum, const Config& cfg, Time t);

  // Set BeamNumber and assign widths and frame accordingly
  void setBeamNumber(unsigned int beamnum, const Config& cfg);

  // Set Boresight from time
  DirectionVector setBoresight(Time t, Uvar terr=Uvar(0,"s"));

  // Return Frame
  const Frame&  getFrame(){ return(frame_);};

  // Returen boresight
  DirectionVector getBoresight(){ return(boresight_);};
  // Get direction vector for an extrema point
  DirectionVector getExtrema(BeamExtremaType etype) const throw(ErrorMessage);

  DirectionVector beamFractionToDirection(Uvar azim_frac, Uvar elev_frac) const;

  Uvec directionToBeamFraction(DirectionVector v) const;

 private:

  DirectionVector boresight_;
  Uvar beam_width_azim_;
  Uvar beam_width_elev_;
  Frame frame_;
  unsigned int beam_number_;
};






#endif
