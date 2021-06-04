#include"Beam2.h"

Beam2::Beam2()
 : frame_("CASSINI_RADAR_3","Cassini")
{
  return;
}
Beam2::Beam2(unsigned int beamnum, const Config& cfg)
 : frame_("CASSINI_RADAR_3","Cassini")
{
  setBeamNumber(beamnum,cfg);
}
Beam2::Beam2(unsigned int beamnum, const Config& cfg, Time t)
 : frame_("CASSINI_RADAR_3","Cassini")
{
  setBeamNumber(beamnum,cfg);
  setBoresight(t);
}
void Beam2::setBeamNumber(unsigned int beamnum, const Config& cfg) 
{
  char c[20];
  
  // Set BeamNumber
  beam_number_=beamnum;

  // Assign widths accordingly
  sprintf(c,"beam_azi_%d",beam_number_);
  beam_width_azim_=cfg[c];
  sprintf(c,"beam_elev_%d",beam_number_);
  beam_width_elev_=cfg[c];

  // Assign frame accordingly
  frame_=Frame("CASSINI_RADAR_" + toStr(beam_number_),"Cassini");
}

DirectionVector Beam2::setBoresight(Time t, Uvar terr) 
{ 
 boresight_=DirectionVector("boresight_",frame_,t+terr,0,0,1);
 Frame inertial_frame("J2000","Earth");
 boresight_.representIn(inertial_frame);
 boresight_.setTime(t);
 boresight_.representIn(frame_);
 return(boresight_);
}

// Returns vectors for various points in the gain pattern
// CENTER = boresight
// TOP = + 1/2 3dB (2-way) elevation width
// BOTTOM = - 1/2 3dB (2-way) elevation width
// LEFT = - 1/2 3dB (2-way) azimuth width
// RIGHT = + 1/2 3dB (2-way) azimuth width
DirectionVector Beam2::getExtrema(BeamExtremaType etype) const throw(ErrorMessage)
{
  DirectionVector tmp=boresight_;
  Uvar azi_step=0.5*beam_width_azim_/sqrt(2);
  Uvar elev_step=0.5*beam_width_elev_/sqrt(2);
  Uvar zero = 0*azi_step;
  switch (etype){
  case CENTER:
    break;
  case TOP:
    tmp.setAzimuthElevation(zero,elev_step);
    break;
  case BOTTOM:
    tmp.setAzimuthElevation(zero,-elev_step);
    break;
  case LEFT:
    tmp.setAzimuthElevation(-azi_step,zero);
    break;
  case RIGHT:
    tmp.setAzimuthElevation(azi_step,zero);
    break;
  default:
    throw ErrorMessage("Beam2::getExtrema: Bad input enum_type value");
  }
  return(tmp);
}

// Returns DirectionVector for fractions of a beam width
DirectionVector Beam2::beamFractionToDirection(Uvar azim_frac, Uvar elev_frac) const
{
  DirectionVector tmp=boresight_;
  Uvar azim, elev;
  boresight_.getAzimuthElevation(azim,elev);
  azim+=azim_frac*beam_width_azim_/sqrt(2);
  elev+=elev_frac*beam_width_elev_/sqrt(2);
  tmp.setAzimuthElevation(azim,elev);
  return(tmp);
}

// Returns DirectionVector for fractions of a beam width
Uvec Beam2::directionToBeamFraction(DirectionVector v) const
{
  Uvec beam_frac("beam_frac",2);
  Uvar azim, elev;
  v.getAzimuthElevation(azim,elev);
  boresight_.getAzimuthElevation(beam_frac(0),beam_frac(1));
  beam_frac(0)=(azim-beam_frac(0))/(beam_width_azim_/sqrt(2));
  beam_frac(1)=(elev-beam_frac(1))/(beam_width_elev_/sqrt(2));
  return(beam_frac);
}
