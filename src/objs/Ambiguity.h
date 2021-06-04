//=============================================================================
// Ambiguity.h
//
// This file contains the Ambiguity Class declaration.
// This Ambiguity class calculates variables depending on target geometry
// at a given time t.
//
// I want to make this program flexible enought so that, once geometry related 
// parameters are calculated, they can be used again for ambiguity calculation
// with different prf, pbw, and pulsegate values.  In order to do that,
// this program will set up a large range/doppler grid system and
// calculate process-window (and backscattering model) 
// independent parameters such as gain, range , doppler, area, lat, lon
// Once process window and backscattering model  have been set,
// calRadarGeomtryFactor and calUsableArea will be used to construc radar echo
// and to determine usabel area pixel by pixel
//
//
// There are some cases in which the user want to calculate ambiguity and 
// thermal snr quickly with fixed process window variables (prf,pbw,pulsegate).
// For such case, I also provde a new calss called AmbiguityFNN(ambiguity 
// using four nearest neighbor).  This AmbiguityFNN will consider only four 
// nearby amb patches separated by prf, and c*pri/2 from the process window.
// After calculating radar echos from four patches, the results will be used
// to compare them with echos from the process window in order to calculate
// signal to amb ratio.
// 
//
//
//
//--------------
//  Ambiguity(const Time& t);//set sar image time and target
//  Ambiguity(); meanted to be used for array of ambiguity, but does not work 
//             because Ambiguity itself uses many array formated data inside
//
//
//  void config(Config& cfg): read beam parameters from config file 
//
//  Time getTime() // get internal time
//  void setTime(const Time& t) 
//               set internal time when ambiguity is calculated
//  void setTarget(const string& target_name, const Frame& target_frame) 
//               set target
//  void setBeam(const Beam& beam) 
//               set the beam to calculate parameters
//  unsigned int getBeamNumber()
//  void setBeamNumber(const unsigned int& beam_num);
//  double getBeamMaxGain();
//  void setBeamMaxGain(const double& beam_max_gain);
//
//  Uvar getAltitude();
//  void setAltitude(const Uvar& altitude);
//      Ambiguity calculate keeps information such as time, altitude
//      which do not depend on spacecraft attitude
//  Uvar getIncidenceAngle();//get boresight incidence angle
//  void setIncidenceAngle(const Uvar& incidenceangle);
//      set boresight incidence angle
//
//  Uvar getMinPrf()
//  Uvar getMaxPrf()
//      Those two will set up range , doppler grid
//      range = range_bore + - ( 2 * c / min_pri)
//      doppler = doppler_bore + - ( 2 * max_prf)
//
//  unsigned int getGridsize()
//  void setGridsize(const unsigned int& Ndat)
//
//  unsigned int getNumberofLooks();
//
//  void clear();//set all the flags false
//
//--------------
//Supporting functions
//--------------
//  unsigned int RangeIndex(const Uvar& range_value);
//  unsigned int DopplerIndex(const Uvar& doppler_value);
//      for a given range and doppler value, they will return index values
//  Uvar RangeInterval()
//  Uvar DopplerInterval()
//     range and doppler interval for uniformly setup (range, doppler) grid
//     range of range / Ngrid, range of doppler / Ngrid
//  void calAmbGeometry()
//     main calculation
//  void calRadarGeometryFactor(const double& muhleman_k1,
//		      const double& muhleman_k2,
//		      const Uvar& wavelength);
//     using backscattering model, it will calculate radar return power
//  double MaxRadarGeomFactor()
//  void calRadarGeometryFactordB(const double& out_of_bounds);
//     convert to dB scale using its max value
//  Uvar calUsableArea(const Uvar& prf,
//	     const Uvar& pbw,
//	     const Uvar& range_bore,
//	     const Uvar& doppler_bore,
//	     const Uvar& pulse_gate,
//	     const Uvar& X0, 
//	     const Uvar& Pn, 
//	     const Uvar& x_res,
//	     const Uvar& rg_res,
//	     const double& muhleman_k1,
//	     const double& muhleman_k2,
//	     const double& signal_to_ambdB, 
//	     const double& noise_equivalent_sigma0dB, 
//	     const double& min_gaindB);
//     calculate usable area using command parameters

//  void calCrosstrackExtent(Uvar& goodbin_crosstrack_start, 
//		   Uvar& goodbin_crosstrack_end); 
//    calculate extent of cross track of usable area
//  Uvar calMultilookUsableArea(const Uvar& prf,
//		      const Uvar& pbw,
//		      const Uvar& range_bore,
//		      const Uvar& doppler_bore,
//		      const Uvar& pulse_gate,
//		      const Uvar& X0, 
//		      const Uvar& Pn, 
//		      const Uvar& x_res,
//		      const Uvar& rg_res,
//		      const double& muhleman_k1,
//		      const double& muhleman_k2,
//		      const Uvar& bpd,
//		      const FloatVector& velocity,
//		      const unsigned int& N_p,
//		      const Uvar& lambda,
//		      const double& signal_to_ambdB, 
//		      const double& noise_equivalent_sigma0dB, 
//		      const double& min_gaindB);
//        calculate usable area for multi-looking
//  void calMLCrosstrackExtent(Uvar& goodbin_crosstrack_start, 
//		     Uvar& goodbin_crosstrack_end);
//         calculate extent of crosstrack of usable area after multilooking    
//  bool mirrorGeomSet();
//      calculate mirror ambiguity if it is asked to so
//      flag will be read in from config file  
//----------------------------------------------------------------------------

#ifndef Ambiguity_H
#define Ambiguity_H
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "Beam.h"
#include "Error.h"
#include "Units.h"
#include "Frame.h"
#include "Array.h"
#include "SimpleArray.h"
#include "Time.h"
#include "Config.h"
#include "Ieb.h"
#include "Io.h"
#include "Config.h"
 
using std::string;
using std::ostringstream;
using std::ofstream;

class Ambiguity
  {
  public:
  //--------------
  // Constructors
  //--------------
  Ambiguity(const Time& t);//set sar image time and target
  Ambiguity();
 
  //--------------------
  //setconfig
  //-------------------
  void config( Config& cfg);

  //---------------------
  //{get/set} method to acess internal variables
  //--------------------
  //time information
  Time getTime();
  void setTime(const Time& t);
  
  //target information
  void setTarget(const string& target_name, const Frame& target_frame);

  //beam
  void setBeam(const Beam& beam);
  unsigned int getBeamNumber();
  void setBeamNumber(const unsigned int& beam_num);
  double getBeamMaxGain();
  void setBeamMaxGain(const double& beam_max_gain);

  //set track frame
  void setTrackFrame(const Frame& ftrack);

  // carrier wavelength 
  void setWavelength(const Uvar& lambda);
  Uvar getWavelength();

  //grid size
  unsigned int getGridsize();
  void setGridsize(const unsigned int& Ndat);

  //mirror ambiguity setting
  void setMirrorAmbiguity(const bool& mirror_set);
  bool getMirrorAmbiguity();

  //number of looks
  unsigned int getNumberofLooks();

  //set backscattering model
  void setBackscattCoefficient(const double& muhleman_k1,
			       const double& muhleman_k2);
  //set process window
  void setProcessWindow(const Uvar& center_range,
			const Uvar& center_doppler,
			const Uvar& prf,
			const Uvar& pbw, 
			const Uvar& pulse_gate);
  //get crosstrack extent min..max
  void getCrosstrackExtent(Uvar& cross_min,Uvar& cross_max, Uvar& crosstrack_length);

  //reset all the flags flase
  void clear();
  

  //--------------
  //Supporting functions for computation
  //--------------
  
  //main job: calculate gain, area, range
  void calAmbGeometry();

  //request interval of grid in range/doppler
  Uvar RangeInterval();
  Uvar DopplerInterval();

  //request index number for a given value
  unsigned int RangeIndex(const Uvar& range_value);
  unsigned int DopplerIndex(const Uvar& doppler_value);

  //calculate G^2 area lambda^2 sigma0 / R^4
  void calRadarGeometryFactor();
 
  //convert to dB scale using the max value
  void calRadarGeometryFactordB(const double& out_of_bounds);

  //calculate usable area using process parameters
  Uvar calUsableArea(const Uvar& X0, 
		     const Uvar& Pn, 
		     const Uvar& x_res,
		     const Uvar& rg_res,
		     const double& signal_to_ambdB, 
		     const double& noise_equivalent_sigma0dB, 
		     const double& min_gaindB);
  //Extent of usable area crosstrack
 
  //calculate multilook usable area using process parameters
  Uvar calMultilookUsableArea(const Uvar& X0, 
			      const Uvar& Pn, 
			      const Uvar& x_res,
			      const Uvar& rg_res,
			      const Uvar& bpd,
			      const FloatVector& velocity,
			      const unsigned int& N_p,
			      const double& signal_to_ambdB, 
			      const double& noise_equivalent_sigma0dB, 
			      const double& min_gaindB);
  

  //check whether mirror ambiguity calculation is requested
  
  
  //--------------------
  // Public data vectors
  //--------------------
  //boresight information
  Uvar range_bore;
  Uvar doppler_bore;
  Uvar altitude;
  Uvar boresight_incidence_angle;

  //pixel 
  Uvec range;
  Uvec doppler;

  Dmat oneway_beam_gain;
  Umat incidenceangle;
  Umat area;
  Umat alongtrack;
  Umat crosstrack;
  Dmat radar_geom_factor;
  Umat radar_geom_factordB;
  Imat number_of_solution;

  Dmat oneway_beam_gain_mirror;
  Umat incidenceangle_mirror; 
  Umat area_mirror;
  Dmat radar_geom_factor_mirror;
  Imat number_of_solution_mirror;
 
  Imat good_bin;
  Imat ML_good_bin;

  Umat amb_ratio_SL;
  Umat amb_ratio_ML;

  Umat thermal_snr_SL;
  Umat thermal_snr_ML;
 

  private:  
  //-----------------
  //internal variables
  //-----------------
  bool time_set_;
  bool target_set_;
  bool beam_set_;
  bool ieb_set_;
  bool config_set_;
  bool mirror_geom_set_;
  bool beam_number_set_;
  bool beam_max_gain_set_;
  bool grid_number_set_;
  bool radar_geom_factor_set_;
  bool peak_radar_geom_factor_set_;
  bool usable_area_set_;
  bool ML_usable_area_set_;
  bool backscatt_model_set_;
  bool process_window_set_;
  bool track_frame_set_;

  Time t_;
  unsigned int nlook_;

  Frame target_frame_;
  Frame ftrack_;
  string target_name_;
  Uvar radius_;

  StateVector state_;

  unsigned int Ngrid_;
  unsigned int Ndat_;

  Uvar lambda_;
  double peak_radar_geom_factor_;

  Beam beam_;
  unsigned int beam_num_;
  double beam_max_gain_;


  Uvar max_prf_;
  Uvar min_prf_;

  Uvar prf_;
  Uvar pri_;
  Uvar pbw_;
  Uvar pulse_gate_; 
  Uvar center_range_;
  Uvar center_doppler_;
  double muhleman_k1_;
  double muhleman_k2_;
  Uvar ML_crosstrack_extent_min_;
  Uvar ML_crosstrack_extent_max_;
  Uvar ML_crosstrack_length_;
  string roll_direction_;
  };



//----------------------------------------
//Unlike aforementioned class, Ambiguity, AmbiguityFNN
// assumes that we know what kind of process parameters going
// to be used.  Therefore, it is not necessary to set a huge grid system
// of many different possible combinations of process parameters.
// As a result, this class will calculate ambiguity and thermal snr
// by setting up one targeted area and four nearby neighboring(FNN)
// ambiguity patches to calculate ambiguity ratio.  Each neighboring
// patch is separated from the targeted area by prf and c * pri /2.
//
//----------------------------------------
class AmbiguityFNN
  {
  public:
  //--------------
  // Constructors
  //--------------
  AmbiguityFNN(const Time& t);//set sar image time and target
  AmbiguityFNN();
  //--------------------
  //setconfig
  //-------------------
  void config(Config& cfg);

  //------------------------------------------------
  //Set methods to calculate Radar Geometry Factor and Along/Cross track
  //------------------------------------------------
  void setTarget(const string& target_name, const Frame& target_frame);
  void setState(const StateVector& sc_state);
  void computeState();
  void setBeam(const Beam& beam);
  void setProcessWindow(const Uvar& central_range, const Uvar& central_doppler,
			const Uvar& prf, const Uvar& pbw, 
			const Uvar& pulsegate);
  void setTrackFrame(const Frame& ftrack);

  //------------------
  //Get method to access some private vairables
  //------------------
  //time
  Time getTime() const;
  void setTime(const Time& t);

  //beam number
  unsigned int getBeamNumber() const;
  double getBeamMaxGain() const;
  void setBeamNumber(const unsigned int& beam_num);
  
  //grid size inside patch
  unsigned int getGridsize() const;
  unsigned int RangeIndex(const Uvar& range_value);//return range index 
  unsigned int DopplerIndex(const Uvar& doppler_value);//return doppler index
  Uvar RangeInterval();//return range interval
  Uvar DopplerInterval();//return doppler interval
 
  //number of looks
  unsigned int getNumberofLooks() const;
  void setNumberofLooks(const unsigned int& nlook);
  
  //center patch(targeted area) number (unsigned int) total_patch_number/2
  unsigned int getCenterPatchNumber() const;
  
  //get crosstrack extent min..max
  void getCrosstrackExtent(Uvar& cross_min,Uvar& cross_max, Uvar& cross_length);
  
  //--------------
  //Supporting functions
  //
  //--------------
  //main engine: do most of the job
  void calAmbGeometry();

  //double ask max value of radar geom factor
  double MaxRadarGeomFactor(); 
 
  //calculate usable area with radar power, noise power , ground resolution
  Uvar calUsableArea(const Uvar& X0, 
		     const Uvar& Pn, 
		     const Uvar& x_res,
		     const Uvar& rg_res);

  //multilook
  Uvar calMultilookUsableArea(const Uvar& X0, 
			      const Uvar& Pn, 
			      const Uvar& x_res,
			      const Uvar& rg_res,
			      const Uvar& bpd,
			      const FloatVector& velocity,
			      const unsigned int& N_p);

  //Is mirror ambiguity requested?  
  bool mirrorGeomSet();


  //-----------------------
  //clear  method: 
  //  set most of the flags false, except
  // config_set, target_set_, ftrack_set_
  // mirror_geom_set_
  //------------------------  
  void clear();


  //--------------------
  // Public data vectors of AmbiguityFNN
  //--------------------
 
  Uvec range;
  Uvec doppler;
  Umat oneway_beam_gaindB;
  Umat incidenceangle;
  Umat lat;
  Umat lon;
  Umat crosstrack;
  Umat alongtrack;
  Umat area;
  D3D radar_geom_factor;
  Umat echo_power;

  //Umat oneway_beam_gain_mirror;
  //Umat incidenceangle_mirror;
  //Umat lat_mirror;
  //Umat lon_mirror;
  //Umat crosstrack_mirror;
  //Umat alongtrack_mirror;
  //Umat area_mirror;
  D3D radar_geom_factor_mirror;
   
  Imat good_bin;
  Imat ML_good_bin;

  Umat thermal_snr_SL;
  Umat amb_ratio_SL;
  Umat thermal_snr_ML;
  Umat amb_ratio_ML;

  Uvar usable_area_SL ;
  Uvar usable_area_ML;
  private:  
  //-----------------
  //internal variables of AmbiguityFNN
  //-----------------
  const static unsigned int Ngrid_;
  const static unsigned int Ndat_;
  const static unsigned int Npatch_;
  const static unsigned int Ncenter_patch_;

  bool time_set_;
  bool target_set_;
  bool state_set_;
  bool beam_set_;
  bool config_set_;
  bool mirror_geom_set_;
  bool beam_number_set_;
  bool beam_max_gain_set_;
  bool process_window_set_;

  bool radar_geom_factor_set_;
  bool peak_radar_geom_factor_set_;
  bool usable_area_set_;
  bool ML_usable_area_set_;
  bool track_frame_set_;

  Time t_;
  Frame ftrack_;
  unsigned int nlook_;

  Frame target_frame_;
  string target_name_;
  Uvar radius_;

  StateVector state_;

  Uvar prf_, pri_,pbw_,pulsegate_;
  Uvar lambda_;
  //Uvar range_bore_,doppler_bore_,thetai_bore_;
  Uvar center_range_,center_doppler_;
  Uvar ML_crosstrack_extent_min_, ML_crosstrack_extent_max_;
  Uvar ML_crosstrack_length_;
  double muhleman_k1_, muhleman_k2_; 
  double peak_radar_geom_factor_;

  Beam beam_;
  unsigned int beam_num_;
  double beam_max_gain_;
  double signal_to_ambdB_ ;
  double noise_equivalent_sigma0dB_; 
  double min_gaindB_; 
  string roll_direction_;
  };




class SARAmb
  {
  public:
  //--------------
  // Constructors
  //--------------
    SARAmb();
  //------
  //destructor
  //--------
  ~SARAmb();
  //--------------------
  //set config
  //-------------------
  void config(Config& cfg);
  void createGrid(const unsigned int& Nrange, const unsigned int& Ndoppler);
  //------------------------------------------------
  //Set methods to calculate Radar Geometry Factor and Along/Cross track
  //------------------------------------------------
  void setTarget(const string& target_name, const Frame& target_frame);
  void setState(const StateVector& sc_state);
  void computeState();
  void setProcessWindow(const Uvar& central_range,
			       const Uvar& central_doppler,
			       const Uvar& lower_range_wrt_central_range,
			       const Uvar& upper_range_wrt_central_range,
			       const Uvar lower_doppler_wrt_central_doppler,
			       const Uvar upper_doppler_wrt_central_doppler,
				const Uvar& prf);

  void setTrackFrame(const Frame& ftrack);

  //------------------
  //Get method to access some private vairables
  //------------------
  //time
  Time getTime() const;
  void setTime(const Time& t);

  //beam number
  void setBeamNumber(const unsigned int& beam_num);
  unsigned int getBeamNumber() const;
  double getBeamMaxGain() const;
  
  //grid size inside patch
  unsigned int getRangeGridsize() const;
  unsigned int getDopplerGridsize() const;
  
   unsigned int RangeIndex(const Uvar& range);//return range index 
  unsigned int DopplerIndex(const Uvar& doppler);//return doppler index
  unsigned int RangeIndex(const double& range_in_km);//return range index 
  unsigned int DopplerIndex(const double& doppler_in_Hz);//return doppler index
  Uvar RangeInterval();//return range interval
  Uvar DopplerInterval();//return doppler interval
 
  //number of looks
  unsigned int getNumberofLooks() const;
  void setNumberofLooks(const unsigned int& nlook);
  
  //center patch(targeted area) number (unsigned int) total_patch_number/2
  unsigned int getCenterPatchNumber() const;
  
  //get crosstrack extent min..max
  void getCrosstrackExtent(Uvar& cross_min,Uvar& cross_max, Uvar& cross_length);
  
  //--------------
  //Supporting functions
  //
  //--------------
  //main engine: do most of the job
  void calAmbGeometry();

  //double ask max value of radar geom factor
  double MaxRadarGeomFactor(); 
 
  //calculate usable area with radar power, noise power , ground resolution
  Uvar calUsableArea(const Uvar& X0, 
		     const Uvar& Pn, 
		     const Uvar& x_res,
		     const Uvar& rg_res);

  //multilook
  Uvar calMultilookUsableArea(const Uvar& X0, 
			      const Uvar& Pn, 
			      const Uvar& x_res,
			      const Uvar& rg_res,
			      const Uvar& bpd,
			      const FloatVector& velocity,
			      const unsigned int& N_p);

  //Is mirror ambiguity requested?  
  bool mirrorGeomSet();

  //Is good bin
  bool isGoodBin(const Uvar& range_input, const Uvar& dop_input, const bool& use_multilook_ambiguity);
  bool isGoodBin(const double& range_in_km, const double& doppler_in_Hz, const bool& use_multilook_ambiguity);

  // Return linear scale ambiguity ratio (ambig to main lobe ratio)
  double getAmbigRatio(const double& range_in_km, const double& doppler_in_Hz, const bool& use_multilook_ambiguity);


  // report information for a single point in doppler and range to a file"
  void pointReport(ofstream& afs, float range_in_km,float doppler_in_hz,
		   bool use_ML);

  //-----------------------
  //clear  method: 
  //  set most of the flags false, except
  // config_set, target_set_, ftrack_set_
  // mirror_geom_set_
  //------------------------  
  void clear();
  // set lambda to be called immediately after clear() if desired
  // avoid using lambda from config file
  // Useful for accounting for doppler due to chirp frequency
  void setLambda(double lambda_in_km){lambda_=lambda_in_km;}

  //--------------------
  // Public data vectors of AmbiguityFNN
  //-------------------- 
  double* range;
  double* doppler;
  double** oneway_beam_gaindB;
  double** incidenceangle;
 
  double** crosstrack;
  double** alongtrack;
  double** area;
  double*** radar_geom_factor;
  double** echo_power;

 
  double*** radar_geom_factor_mirror;
   
  int** good_bin;
  int** ML_good_bin;

  double** thermal_snr_SL;
  double** amb_ratio_SL;
  double** thermal_snr_ML;
  double** amb_ratio_ML;

  Uvar usable_area_SL ;
  Uvar  usable_area_ML;


  //------------------------------------------------
  //temporary container: declared once to save time
  //-----------------------------------------------

  double** range_grid;
  double** doppler_grid;

  int** no_of_solution_grid;
  double** oneway_beam_gain_grid;
  double** incidenceangle_grid;
  Array2D<PositionVector> surface_intercept_grid;

  int** no_of_solution_mirror_grid;
  double** oneway_beam_gain_mirror_grid;
  double**  incidenceangle_mirror_grid;;
  Array2D<PositionVector> surface_intercept_mirror_grid;
  

  double** range_pixel;
  double** doppler_pixel;

  double*** oneway_beam_gain_pixel;
  double*** area_pixel;
  double*** incidenceangle_pixel;
  
  double*** oneway_beam_gain_mirror_pixel;
  double*** area_mirror_pixel;
  double*** incidenceangle_mirror_pixel;
  
  int***  no_of_solution_pixel;
  int*** no_of_solution_mirror_pixel;
  private:  
  //-----------------
  //internal variables of AmbiguityFNN
  //-----------------
  const static unsigned int Npatch_;
  const static unsigned int Ncenter_patch_;

  bool time_set_;
  bool target_set_;
  bool state_set_;
  bool beam_set_;
  bool config_set_;
  bool mirror_geom_set_;
  bool beam_number_set_;
  bool beam_max_gain_set_;
  bool process_window_set_;

  bool radar_geom_factor_set_;
  bool peak_radar_geom_factor_set_;
  bool usable_area_set_;
  bool ML_usable_area_set_;
  bool track_frame_set_;
  bool grid_set_;

  Time t_;
  Frame ftrack_;
  unsigned int nlook_;

  Frame target_frame_;
  string target_name_;
  double  radius_;

  StateVector state_;

  double prf_, pri_;
  
  double lambda_;
  double center_range_, center_doppler_;
  double  lower_doppler_, upper_doppler_;
  double  lower_range_, upper_range_;
  double ML_crosstrack_extent_min_, ML_crosstrack_extent_max_;
  double ML_crosstrack_length_;
  double muhleman_k1_, muhleman_k2_; 
  double peak_radar_geom_factor_;

  vector<Beam>  beam_;
  vector<double> beam_max_gain_;

  unsigned int beam_id_;
  double signal_to_ambdB_ ;
  double noise_equivalent_sigma0dB_; 
  double min_gaindB_; 
  unsigned int Nrange_;
  unsigned int Ndoppler_;
  unsigned int Nrange_grid_;
  unsigned int Ndoppler_grid_;
  string roll_direction_;
  void createArrays();
  void deleteArrays();
  };





//-------------------------------------------------------
//Earlier versions of ambiguity calculation tools: not associated 
//with Ambiguity calss
// This tells the history of ambiguity calculation evolution.
// In the beginning, it started from subroutine and then was upgrated
// into Class 
//-------------------------------------------------------
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
//
//void cal_ambiguity_search_ieb_par()
//  This program uses exactly the same methods as cal_ambiguity but
// stores only small part of information. This routine will be used
// by search_ieb_par, which is expected to produce about 1000 files
// containint amb calculation results for each different altitude,
// incidenceangle, prf, and dutycycle.
//----------------------------------------------

void cal_ambiguity(const string& target_name,
		   Beam& beam,
		   const Frame& ftitan,
		   const StateVector& sc_state,
		   const Frame& track_frame,
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



