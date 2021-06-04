#ifndef SARPROCPARAMS_H
#define SARPROCPARAMS_H

static const char rcs_id_sarprocparams_h[] =
  "@(#) $Id: SARProcParams.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include"Config.h"
#include"Units.h"
#include"BurstData.h"
#include"Projections.h"
#include "Flyby.h"
#include "Ambiguity.h"
#include "RasList.h"

class L1I;
class SARProcParams{

 public:

  enum regionTypeE {FULL, BURSTNO, LATLON, LINENO, 
		    ABSTIME, CLATIME, EPOCHTIME, TRIGGERTIME};
  enum signalTypeE {NOMINAL, THEOR_NOISE, DATA_NOISE};
  bool simBaqForNoise;

  SARProcParams(Config& cfg,int max_sab_counter);
  ~SARProcParams();
 
  void topoInit(float oneway_gain_cutoff, float prf_percent);
  const Uvar& Xcal(const Time& time, int beam_num, int radar_mode);
 const Uvar& powerToDataNumber(int radar_mode);
 const Uvar& getTsys(int radar_mode);
 
  void replaceSignalWithNoise(int radar_mode, L1I& b, float noise_var_est) const;
  void estimateNoiseSubtractionParameters(L1I& l1i, RasList* raslist);
  Uvar NoiseVarianceToSystemTemperature(float noise_var, L1I& b) const;
  // determine whether or not a burst is whether the desired processing region
  bool inRegion(const OblCylProj&, const BurstData&);
  void setDopplerProfile(const Time& time, int beam_num,
			 Uvar& slope, Uvar& doppler_boresight, 
			 Uvar& range_boresight,
			 Uvar& rate_slope,
			 Uvar& doppler_rate_boresight,
			 int sab_counter=-1) const;
  // Determine whether or not to skip aburst in processing
  bool skipBurst(const OblCylProj&, const BurstData&); // perform checks
  bool skipBurst(int sab_counter); // simplified check performed later during processing
  void reportSkippedBursts(ofstream& afs);

  //Compute Useable area mask
  bool isUseable(double fdop_in_Hz, double range_in_km); 

  // return Ambiguity Ratio (not in dB)
  double getAmbigRatio(double fdop_in_Hz, double range_in_km); 

  //Initialize Useable Area mask calculation
  void initUseableCalculation(const BurstData& data,
			      const StateVector& sc_state,
			      const Uvar& bore_range,
			      const Uvar& bore_doppler,
			      double lambda_in_km,
			      double look_repeat_time_in_s);
  float getGain2Thresh(const Time& time, int beam_num, int sab_counter);
  StateVector getClosestApproachState();
  string LookDirection();

  // Processor Algorithm Options
  bool single_beam_mode;
  int beam_number; // only used if in single beam mode
  bool process_SAR_mode_only;
  bool process_SARH_mode_only;
  bool process_SARL_mode_only;
  bool process_scat_mode_only;
  bool process_alt_mode_only;
  bool strip_processor_on;
  bool azimuth_range_orthogonalization_on;
  bool azimuth_deramp_processing_on;
  bool calibration_on;
  bool gain_cutoff_enabled;
  bool azimuth_width_cutoff_enabled;
  float azimuth_width_percent;
  bool full_useable_calc_enabled;
  bool use_multilook_ambiguity;
  bool set_doppler_centroid; // sets doppler (and range) from config file
  bool linear_dopran_approx;
  bool force_bidr_region;
  bool remove_doppler_phase_ramp;
  bool enable_data_dependent_cal;
  bool set_burst_range_dop;
  signalTypeE signal_type;
  double theor_noise_var_dn;
  FILE* data_noise_out_fp;
  int forced_first_line;
  int forced_last_line;
  Uvar azimuth_res_coeff;
  Uvar range_res_coeff;
  double maximal_gain_corr;

  Uvar transmit_power, Xcorr_constant,  carrier_frequency, lambda;
 // SARAmb object used for full up useability calculation if desired
  SARAmb amb_sar;
  int num_sabs;
  char* quality_override; // 0=BAD, 1=GOOD, -1=NO OVERRIDE;
  float* doppler_offset;  // array for burst by burst doppler offsets if
                          // desired
  float* range_offset;    // array for burst by burst range offsets if desired
  string cout_line;  // holder for standard info output for each burst by
                     // sar_proc;


  // noise subtraction values
  float thermal_noise_offset;
  float quant_noise_offset;
  float quant_noise_scale;
  bool noise_subtraction_on;
  float sarproc_gain;
 private:

  void readQualityOverrideFile(string file);
  void readBurstRangeAndDopplerFile(string file);
  // parameters for Xcal computation
  Uvar Xcal_;       // calibration (non-geometry) portion of Xfactor
  Uvar Xconstant_;  // constant part of Xcal
  Uvar Xbeam_[5];   // beam dependent correction to Xcal
  Uvar Xmode_[16];  // radar mode dependent correction to Xcal
  Uvar Tsys_[16];  // radar mode dependent correction to Xcal
  Uvar Xtimevar_;   // time varying correction to Xcal

  int Xtimevar_poly_order_; // order of polynomial describing Xtimevar
                           // 0= constant, 1=linear , etc.
  Uvec Xtimevar_coeff_; // coefficients of polynomial describing Xtimevar
                       // as a function of time from closest approach
                       // (May want to change to time from turn on ....)

  // Doppler centroid and range as a function of time polynomials
  // used when SET_DOPPLER_CENTROID config parameter is non-zero
  // Only used for debugging .....
  // Should not be used for normal ops ... 
  int fdopc_poly_order_;
  Umat fdopc_coeff_;
  Uvar fdopc_std;
  Uvar fdopc_mean;
  Uvar time_std;
  Uvar time_mean;
  
  int rangec_poly_order_;
  Umat rangec_coeff_;
  Uvar rangec_mean;
  Uvar rangec_std;

  // flyby object used for time from closest approach computations
  Flyby flyby_;

  
  Dmat gain_cutoff_coeff_;
  int gain_cutoff_poly_order_;

  regionTypeE region_type_;

  Time trigger_time_;
  int burstno_end, burstno_start;
  Time abstime_start, abstime_end;
  Uvar reltime_start, reltime_end;
  Uvar lat_south, lat_north, elon_east, elon_west;

  // Burst skip statistics
  int num_bursts_processed;
  int num_bursts_skipped;
  int min_sab_counter_processed;
  int max_sab_counter_processed;
  int num_bursts_outside_region;
  int num_bursts_forced_bad;
  int num_bursts_forced_good;
  int num_bursts_wrong_beam;
  int num_bursts_bad_ckernel;
  int num_bursts_bad_ephemeris;
  int num_bursts_off_limb;
  int num_bursts_bad_downlink;
  int num_bursts_wrong_rmode;
  int num_calibration_bursts;
  int num_bursts_no_pulses;
  

};

#endif
