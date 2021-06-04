//=============================================================================
// AmbiguityData.h
//
// This file contains the AmbiguityData class declarations.
// The AmbiguityData class provides a convenient package to hold vectors
// of ambiguity data.
//
// Interface summary:
//
// AmbiguityData methods:
//
//   Construction:
//
//   AmbiguityData() //use default constructor
//--------------
//I/O
//-------------
//  void readAmbiguityData(ifstream& amb_data_file) throw(ErrorMessage);
//             read ambiguitydata containing detailed information about
//             range and dop patches
//  void displayAmbiguitySummary(ofstream& outputsummary) throw(ErrorMessage); 
//             make an output file of ambiguity calculation results

//----------------
//get method to access interval variables
//---------------
//  unsigned int getNrangepatch() throw(ErrorMessage);
//                return # of range patches
//  unsigned int getNdoppatch() throw (ErrorMessage);
//                return # of doppler patches
//  unsigned int getNrangebin() throw(ErrorMessage);
//                return # of range bins
//  unsigned int getNdopbin() throw (ErrorMessage);
//                return # of dop bins
//  string getEpochtime() throw(ErrorMessage);
//                return epoch time
//  Uvar getSarimagingtime() throw(ErrorMessage);
//                return time when ambiguity calculation is done, w.r.t. epoch
//--------------
//Supporting functions
//--------------
//  void BeamSwathLonLat(const unsigned int& i_beam, Uvec& lon_track, 
//		       Uvec& lat_track) throw(ErrorMessage);
//              return ground swath of all 5 beams in terms of 
//              planet lat and lon coordinates
//  void BeamSwathCrossAlong(const unsigned int& i_beam, Uvec& cross_track,
//			   Uvec& along_track) throw(ErrorMessage);
//              return ground swatch of all 5 beams in terms of 
//              crosstrack/alongtrack distance
//  void CrossAlongPatches(const unsigned int& i_beam, 
//			   const unsigned int& i_patch, 
//			   Uvec& cross_track,
//			   Uvec& along_track) throw(ErrorMessage);
//              return range and dop patches of beam i_beam in terms of
//              crosstrack/alongtrack distance
//  void SignalvsCrosstrack(const unsigned int& i_beam,
//			  const unsigned int& alongtrack_index,
//			  Uvec& return_power,
//			  Uvec& cross_track) throw(ErrorMessage); 
//               return radar  signal vs crosstrack distnace
//  void AmbiguityvsCrosstrack(const unsigned int& i_beam,
//			  const unsigned int& alongtrack_index,
//			  Uvec& return_power,
//			  Uvec& cross_track) throw(ErrorMessage);
//               return ambiguity power vs crosstrack disntace
//  void AmbratiodBvsCrosstrack(const unsigned int& i_beam,
//			  const unsigned int& alongtrack_index,
//			  Uvec& ratio_dB,
//			  Uvec& cross_track) throw(ErrorMessage);
//                return P_main_process_window/Pamb_patches vs crosstrack
//  void AntennadBvsCrosstrack(const unsigned int& i_beam,
//			  const unsigned int& alongtrack_index,
//			  Uvec& antenna_dB,
//			  Uvec& cross_track) throw(ErrorMessage);
//                return antenna gain vs crosstrack
//  void SignalvsAlongtrack(const unsigned int& i_beam,
//			  const unsigned int& crosstrack_index,
//			  Uvec& return_power,
//			  Uvec& along_track) throw(ErrorMessage);
//                 return radar signal vs alongtrack distnace
//  void AmbiguityvsAlongtrack(const unsigned int& i_beam,
//			  const unsigned int& crosstrack_index,
//			  Uvec& return_power,
//			  Uvec& along_track) throw(ErrorMessage);
//                 return ambiguity power vs alongtrack distance
//  void AmbratiodBvsAlongtrack(const unsigned int& i_beam,
//			  const unsigned int& crosstrack_index,
//			  Uvec& ratio_dB,
//			  Uvec& along_track) throw(ErrorMessage);
//                 return radar_signal/amb_power vs alongtrack distance
//  void AntennadBvsAlongtrack(const unsigned int& i_beam,
//			  const unsigned int& crosstrack_index,
//			  Uvec& antenna_dB,
//			  Uvec& along_track) throw(ErrorMessage);
//                 return antenna gain vs alongtrack distance
//  void Goodbins(const unsigned int& i_beam,Uvec& cross_track,
//		Uvec& along_track) throw(ErrorMessage);
//                  return locations of good bin in terms of alongtrack
//                  crosstrack distance
//  void Goodbins_Pn_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
//		   Uvec& along_track) throw(ErrorMessage);
//                  return locations of bins whose Pn/Pamb is larger than 1.0
//  void Goodbins_Signal_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
//		   Uvec& along_track) throw(ErrorMessage);
//                  return locations of bins whose P_signal/P_amb is 
//                 larger than 14 dB
//  void Goodbins_Nesigma0(const unsigned int& i_beam, Uvec& cross_track,
//		   Uvec& along_track) throw(ErrorMessage);
//                  return locations of bins whose noise_equivalent_sigma0
//                  is smaller than -10 db
//  void Processwindowedges(const unsigned int& i_beam,Uvec& beam_edges_cross,
//		  Uvec& beam_edges_along) throw(ErrorMessage);
//                 return locations of four edges of the main process window 
//  void Is_good_bin_Alongtrack(const unsigned int& i_beam,
//		      const unsigned int& i_indicator,
//		      const unsigned int& crosstrack_index, 
//		      Uvec& bin_indicator,
//		      Uvec& along_track) throw(ErrorMessage);
//                 return locations of good bins vs alongtrack
//  void Is_good_bin_Crosstrack(const unsigned int& i_beam,
//		      const unsigned int& i_indicator,
//		      const unsigned int& alongtrack_index, 
//		      Uvec& bin_indicator,
//		      Uvec& cross_track) throw(ErrorMessage);
//                 return locations of good bins vs crosstrack
//  Uvar total_usable_area() throw(ErrorMessage);
//                  return usable area of all 5 beams combined, excluding
//                 overlapping areas
//  Uvar usable_area_wo_Pn_Pamb_requirement() throw (ErrorMessage)
//                 return usable area without the requirement of Pn being
//                 larger than Pamb
//
//
//=============================================================================

#ifndef AMBIGUITYDATA_H
#define AMBIGUITYDATA_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "Array.h"
#include "Units.h"
#include "Error.h"


using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

class AmbiguityData
  {
  public:



  //--------------
  // Constructors
  //--------------

  
  AmbiguityData() throw(ErrorMessage);

  

  //--------------
  //I/O
  //-------------
  void readAmbiguityData(ifstream& amb_data_file) throw(ErrorMessage);
  void displayAmbiguitySummary(ofstream& outputsummary) throw(ErrorMessage);  
  //----------------
  //get method to access interval variables
  //---------------
  unsigned int getNrangepatch() throw(ErrorMessage);
  unsigned int getNdoppatch() throw (ErrorMessage);
  unsigned int getNrangebin() throw(ErrorMessage);
  unsigned int getNdopbin() throw (ErrorMessage);
  string getEpochtime() throw(ErrorMessage);
  Uvar getSarimagingtime() throw(ErrorMessage);

  //--------------
  //Supporting functions
  //--------------
  void BeamSwathLonLat(const unsigned int& i_beam, Uvec& lon_track, 
		       Uvec& lat_track) throw(ErrorMessage);
  void BeamSwathCrossAlong(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage);
  void CrossAlongPatches(const unsigned int& i_beam, 
			   const unsigned int& i_patch, 
			   Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage);
  void SignalvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& return_power,
			  Uvec& cross_track) throw(ErrorMessage); 
  void AmbiguityvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& return_power,
			  Uvec& cross_track) throw(ErrorMessage);
  void AmbratiodBvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& ratio_dB,
			  Uvec& cross_track) throw(ErrorMessage);
  void ThermalSNRdBvsCrosstrack(const unsigned int& i_beam,
				const unsigned int& alongtrack_index,
				Uvec& ratio_dB,
				Uvec& cross_track) throw (ErrorMessage);
  void AntennadBvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& antenna_dB,
			  Uvec& cross_track) throw(ErrorMessage);
  void SignalvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& return_power,
			  Uvec& along_track) throw(ErrorMessage);
  void AmbiguityvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& return_power,
			  Uvec& along_track) throw(ErrorMessage);
  void AmbratiodBvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& ratio_dB,
			  Uvec& along_track) throw(ErrorMessage);
  void ThermalSNRdBvsAlongtrack(const unsigned int& i_beam,
				const unsigned int& crosstrack_index,
				Uvec& ratio_dB,
				Uvec& along_track) throw(ErrorMessage);
  void AntennadBvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& antenna_dB,
			  Uvec& along_track) throw(ErrorMessage);
  void Goodbins(const unsigned int& i_beam,Uvec& cross_track,
		Uvec& along_track) throw(ErrorMessage);
  void Goodbins_Pn_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage);

  void Goodbins_Signal_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage);

  void Goodbins_Nesigma0(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage);

  void Processwindowedges(const unsigned int& i_beam,Uvec& beam_edges_cross,
			  Uvec& beam_edges_along) throw(ErrorMessage);
   
  void Is_good_bin_Alongtrack(const unsigned int& i_beam,
			      const unsigned int& i_indicator,
			      const unsigned int& crosstrack_index, 
			      Uvec& bin_indicator,
			      Uvec& along_track) throw(ErrorMessage);

  void Is_good_bin_Crosstrack(const unsigned int& i_beam,
			      const unsigned int& i_indicator,
			      const unsigned int& alongtrack_index, 
			      Uvec& bin_indicator,
			      Uvec& cross_track) throw(ErrorMessage);

  Uvar total_usable_area() throw(ErrorMessage);
  Uvar usable_area_wo_Pn_Pamb_requirement() throw (ErrorMessage);

  //--------------------
  // Public data vectors
  //--------------------

  Dvec max_beam_gain;
  Uvar lambda;
  Uvar frequency,tau_p,Pt,Tsys;
  Uvar BR, Bn, fadc, radius, Pn;    
  Uvar prf, min_bpd, pbw;
  double  pulse_bandwidth_ratio;
  Uvec tau_3db, pulse_gate;
  Uvar  receive_window;
  Uvar bpd;
  unsigned int N_p, Ns_rw;
  double duty_cycle;
  
  Uvar rtt;
  double muhleman_k1;
  double  muhleman_k2;
  Uvec range_offset,frequency_offset;
  Uvec range_3db, X0;
  Uvar x_res,rg_res;
 
 
  //dummy Imat and Dmat to make I3D and D3D work
  //sounds silly, but this is the way to make things work now
  Imat dummImat;
  Dmat dummDmat;
  
  //usable area of each beam
  Uvec usable_area;
  
  //range and dop boundaries of each range/dop patches for each beam
  Umat range_lower,range_upper;
  Umat dop_lower,dop_upper;
  
  //range/dop center of each beam's processing window
  Uvec range_center;
  Uvec dop_center;

  //area of process window for each beam
  Uvec area_process_window;
  
  U3D avg_incidence;
  D3D avg_backscatter;
  D3D avg_antenna;
  U3D avg_azi,avg_elev;
  U3D avg_area;
  U3D lat, lon, range;
  U3D alongtrack,crosstrack;
  U3D doppler, thetai; 
  U3D center_area;
  D3D bin_ratiodB;
  U3D normalized_bin_power;
  D3D bin_antenna;
  D3D bin_antennadB;
  I3D is_good_bin;
  I3D is_good_bin_amb;
  I3D is_good_bin_nesigma0;
  I3D is_good_bin_totalSNR;	
  I3D is_good_bin_pn_to_amb;
  
  U3D patch_sum;    
  U3D bin_power_main;
  U3D bin_power_amb;
  U3D Return_signal;
  U3D total_noise;

  U3D center_patch_along;
  U3D center_patch_cross;

  
 
  
  
  private:
  
  //-----------------
  //internal variables
  //-----------------
  bool read_data_;
  unsigned int Nrange_patch_;
  unsigned int Ndop_patch_;
  unsigned int i_range_center_;
  unsigned int j_dop_center_;
  unsigned int i_patch_center_;
  unsigned int Nrange_bin_;
  unsigned int Ndop_bin_;
  string epoch_time_;
  Uvar sar_imaging_time_;
  
  
  
    
  };

#endif





