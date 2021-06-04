#ifndef RADFUNCTIONS_H
#define RADFUNCTIONS_H

#include <string>
#include <vector>
#include "Units.h"
#include "Error.h"
#include "Array.h"


//calibrate system temperature using YFACTOR means when sc is looking cold sky
// i_beam_rec: system temperature for record number i_beam_rec
// Y_factor: Y = PN1/PN2 = (Ts + T1)/(Ts+T2)
// avg_time: take a mean value of system temperature averaged over this time
// Tcoldsky: cold sky temperature
// beam_rel_time: time tag
// beam_ncnt_rl: RL count
// beam_ncnt_radio: beam radiometer count
// beam_rlotmp: RL temperature sensor reading
// is_beam_coldsky: boolian vector (true: cold sky pointing, false: on target)
//-----------------------------------------------------------------------
//this method does not work well because RL input path and Radiometer input
// path are different.  Radiometer input path has an antenna and front end
// waveguid, which add noises due to losses
//------------------------------------------------------------------------
Uvar  CalibrateTsysandYFactor(const unsigned int& i_beam_rec, 
			      double& Y_factor,
			      const Uvar& avg_time,
			      const Uvar& Tcoldsky,
			      const vector<Uvar>& beam_rel_time,
			      const vector<Uvar>& beam_ncnt_rl, 
			      const vector<Uvar>& beam_ncnt_radio,
			      const vector<Uvar>& beam_rlotmp, 
			      const vector<bool>& is_beam_coldsky);


//----------------------------------------
//This method computes system temperature while assuming system 
// gain is constant
// Not a good method in most of time.
// We found out that radiometer output heavily depends on system gain, not 
// system temperature
// Could be useful when gain calibration is not avail
//----------------------------------------
Uvar  CalibrateTsysUsingConstGain(const unsigned int& i_beam_rec,
				  const Uvar& gain_sys,
				  const Uvar& avg_time, 
				  const Uvar& Tcoldsky,
				  const vector<Uvar>& beam_rel_time,
				  const vector<Uvar>& beam_ncnt_radio,
				  const vector<bool>& is_beam_coldsky) ;


//--------------------------------------------------------------------
//Compute system gain by assuming system noise temperature is constant
//--------------------------------------------------------------------
Uvar  CalibrateGsysUsingConstTsys(const unsigned int& i_beam_rec,
				  const Uvar& Tsys,
				  const Uvar& avg_time, 
				  const Uvar& Tcoldsky,
				  const vector<Uvar>& beam_rel_time,
				  const vector<Uvar>& beam_ncnt_radio,
				  const vector<bool>& is_beam_coldsky) ;


//-----------------------------------
//compute average radiometer count and RL count
//-------------------------------------
void  computeAvg_radio_rl_rlotmp(const unsigned int& i_beam_rec,
				 const Uvar& avg_time, 
				 const vector<Uvar>& beam_rel_time,
				 const vector<Uvar>& beam_ncnt_radio,
				 const vector<Uvar>& beam_ncnt_rl,
				 const vector<Uvar>& beam_rlotmp,
				 const vector<bool>& is_beam_coldsky,
				 Uvar& avg_cnt_radio,
				 Uvar& avg_cnt_rl,
				 Uvar& avg_rlotmp) ;

//----------------------------------
//compute system gain by using cold load (RL)
//----------------------------------
void compute_GT_using_RL(const unsigned int& i_beam_rec,
			 const double& frontend_loss_figure ,
			 const Uvar& avg_time,
			 const Uvar& Tcoldsky,
			 const vector<Uvar>& beam_rel_time,
			 const vector<Uvar>& beam_ncnt_radio,
			 const vector<Uvar>& beam_ncnt_rl,
			 const vector<Uvar>& beam_rlotmp,
			 const vector<Uvar>& beam_wgbtmp,
			 const vector<bool>& is_beam_coldsky,
			 Uvar& system_gain,
			 Uvar& receiver_temp);

//-------------------------
//Compute time-averaged RL counts
//---------------------------
void compute_time_averaged_rl(const Uvar& avg_time,
			      const vector<Uvar>& beam_ncnt_rl,
			      const vector<Uvar>& beam_rel_time,
			      vector<Uvar>&  avg_beam_ncnt_rl);

#endif
