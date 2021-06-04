//=============================================================================
// Constants.h
//
// This file contains physical constants in an un-named namespace.
// They will appear as global UnitVar's in any file that includes Constants.h
// 
//=============================================================================

#ifndef Constants_H
#define Constants_H

#include "Units.h"
#include "Frame.h"

namespace 
  {
  Uvar speed_light(2.997924580e5,"km/s");
  double speed_light_in_kps=2.997924580e5;
  Uvar boltzmann_constant(1.3806503e-29,"kg km km/(s s K)");
  Uvar gravitational_constant(6.673e-20,"km km km/(s s kg)");
  Uvar free_space_permittivity(8.854187818e-3,"uF/km");
  Uvar carrier_frequency(1.378e10,"Hz");
  Uvar slo_frequency(1e7,"Hz");
  Uvar sarh_adc(2e6,"Hz");
  Uvar sarl_adc(1e6,"Hz");
  Uvar alth_adc(1e7,"Hz");
  Uvar altl_adc(2.5e5,"Hz");
  Uvar sarh_rcv_bandwidth(9.35e5,"Hz");
  Uvar sarl_rcv_bandwidth(4.68e5,"Hz");
  Uvar alth_rcv_bandwidth(4.675e6,"Hz");
  Uvar altl_rcv_bandwidth(1.17e5,"Hz");
  Uvar sarh_chirp_bandwidth(8.5e5,"Hz"); 
  Uvar sarl_chirp_bandwidth(4.25e5,"Hz");
  Uvar alth_chirp_bandwidth(4.25e6,"Hz");
  Uvar altl_chirp_bandwidth(1.06e5,"Hz");
  string km_str="km";
  string km_per_s_str="km/s";
  string cassini_str="Cassini";
  string seconds_str="s";

  }

const double pi=M_PI;

// Unit  conversion constants for use with UnitVar no_unit_support mode
// Should use coerce_base instead unless speed is really an issue.
const double WtoMW = 1e-6;
const double mAtoA = 1e-3;
const double MHztoHz = 1e6;
const double mstos = 1e-3;
const double nstos = 1e-9;
const double ustos = 1e-6;
const double VtoMV = 1e-6;
const double degtorad=pi/180.0;
const double radtodeg= 180.0/pi;


// Global variables set by Frame::config
extern  SpiceInt cassini_spice_id;
extern  SpiceInt titan_spice_id;
extern  SpiceInt default_target_spice_id;
extern  SpiceInt sun_spice_id;
extern  SpiceInt earth_spice_id;
extern  SpiceInt iau_titan_frame_spice_id;
extern  SpiceInt j2000_frame_spice_id;
extern  SpiceInt beam_frame_spice_id[5];
extern  SpiceInt default_target_frame_spice_id;
extern  PositionVector default_target_radii;
extern  string default_target_name;
extern  bool default_ring_target;



#endif 



