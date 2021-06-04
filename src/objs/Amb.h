//=============================================================================
// Amb.h
//
// This file contains a subroutine to calculate ambiguity
//
//*************************************
//void cal_ambiguity()
//    This program calculate ambiguity performance  of main processing window 
// defined by pulse_bandwidth(pbw = 0.8 * prf) and pulse gate at time t
// set by sc_state.  It computes radar return power from  range and dop patches
// (Nrange_patch x Ndop_patch of them).  Each patch has Nrange_bin 
// x Ndop_bin bins.  
//    The program calculates radar return from each bin based on 
// backscattering model, range, area, and antenna beam gain.  Each bin
// in indexed according to range(i) and doppler frequency(j).  For bins outside
// the main process window, bins with same index(i,j)from different range and 
// dop patches are added up to produce ambiguity power(i,j) and then compared
// with radar return from bin(i,j) inside the main process window.  In addition
// to ambiguity power ratio, the program also calculates noise equivalent
// sigma0 and thernal SNR.
//   As an output, the program write down  the following information 
//into a file:
//   min and max of range/dop inside each patch
//   average incidenceangle, antenna gain, azi, elev, 
//      radar return from each patch, area of each patch
//   lat, lon, incidence angle, range, doppler, along/cross track of each bin
//   radar return power of each bin inside the main process window
//   ambiguity power of each bin from range/dop patches excluding the main
//           process window
//   noise equivalent sigmal0 of each bin
//   good bin indictor
//     (1) good amb ratio?
//     (2) good noise equivalent sigma0?
//     (3) good total SNR?
//     (4) good Pn/Pamb ratio?
//     (5) satisfying all the above conditions?
//***************************************

//---------------------------------------------
//void cal_ambiguity_search_ieb_par()
//  This program uses exactly the same methods as cal_ambiguity but
// stores only small part of information. This routine will be used
// by search_ieb_par, which is expected to produce about 1000 files
// containint amb calculation results for each different altitude,
// incidenceangle, prf, and dutycycle.
//----------------------------------------------

#ifndef Amb_H
#define Amb_H
#include <stdio.h>
#include <string>
#include <strstream>
#include <iostream>
#include <fstream>
#include "Beam.h"
#include "Error.h"
#include "Units.h"
#include "Frame.h"
#include "Array.h"
#include "Time.h"
#include "Config.h"

using std::string;
using std::ostrstream;


void cal_ambiguity(const string& target_name,
		   Beam& beam,
		   const Frame& ftitan,
		   const StateVector& sc_state,
		   const Uvar& prf,
		   const Uvar& pulse_gate,
		   const Uvar& pbw,
		   const Uvar& lambda,
		   const Uvar& X0,
		   const Uvar& Pn,
		   const Uvar& x_res,
		   const Uvar& rg_res,
		   const double& Ni,
		   const double&  muhleman_k1,
		   const double&  muhleman_k2,
		   const unsigned int& Nrange_patch,
		   const unsigned int& Ndop_patch,
		   const unsigned int& Nrange_bin,
		   const unsigned int& Ndop_bin,
		   const Uvar& range_offset,
		   const Uvar& frequency_offset,
		   ofstream& outputfile )
     throw(Unit::UnitError,ErrorMessage);


void cal_ambiguity_geometry_table(const string& target_name,
				   Beam& beam,
				  const Frame& ftitan,
				  const StateVector& sc_state,
				  const Uvar& lambda,
				  const double&  muhleman_k1,
				  const double&  muhleman_k2,
				  const Uvec& dop_freq_values,
				  const Uvec& range_values,
				  Dmat& cell_one_way_gain,
				  Dmat& cell_radar_geom_factor,
				  Umat& cell_cross_track,
				  Umat& cell_along_track,
				  Umat& cell_incidenceangle,
				  Umat& cell_area,
				  Umat& cell_lat,
				  Umat& cell_lon,
				  Dmat& cell_sigma0,
				  Uvec& range_axis,
				  Uvec& dop_axis)
     throw(Unit::UnitError,ErrorMessage);




 double sincgain(const Uvar& azim, 
		 const Uvar& elev, 
		 const Uvar& azimuth_width,
		 const Uvar& elevation_width,
		 const double& max_gain) throw(Unit::UnitError,ErrorMessage);
#endif 






