//============================================================================
// PointTargetSim.h
//
// This file contains the PointTargetSim class declaration.
// The PointTargetSim class implements a point target simulation
// It is based fundamentally on time like most of the RAS classes.
// This means that calculations are performed at a particular time with the
// SPICE system specifying the positions,velocities of various objects at
// that time.
// In the future, an array of times may be supported.
//=============================================================================

#ifndef POINT_TARGET_SIM_H
#define POINT_TARGET_SIM_H

#include"CassiniSim.h"
#include"L1B.h"
#include"Profiling.h"
#include"Baq.h"
#include"Ieb.h"
#include"IebList.h"
#include"Flyby.h"
#include"MFile.h"
#include"CAGTable.h"

class PointTargetSim
: private CassiniSim
{
 public:
  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  enum ParamSourceE { CONFIG, L1BFILE, IEB };
  PointTargetSim(Config& cfg, ParamSourceE p);
  ~PointTargetSim();
  void run(Config& cfg); // needs config to rewrite L1B

 private:
  //-------------------------- 
  // Internal representation
  //-------------------------- 

  bool target_info_in_diag_file;
  bool target_beyond_horizon_;
  bool ideal_tracking_on_;
  bool compute_echo_amplitude_;
  bool simulate_backscatter_;
  bool simulate_area_;
  bool area_is_relative_;
  bool ideal_gain_setting_;
  bool use_config_bpd_;
  bool simulate_noise_;
  bool simulate_quantization_;
  bool use_random_phase_shift_;
  bool use_random_location_shift_;
  bool simulate_baq_;
  bool only_write_targeted_bursts;

  int sample_skip;
  Uvar max_target_angle;

  // Backscatter model coefficients

  double muhleman_k1_;
  double muhleman_k2_;
  Time start_time;
  Time burst_start_time;
  Time end_time;
  L1B* l1b_template;
  L1B* l1b;
  IebList* ieb_file;
  Uvar ieb_delta_trigger_time;
  unsigned int beam_order[6];
  unsigned int beam_order_idx;
  string target;
  unsigned int num_targets;
  Uvec target_lat;
  Uvec target_lon;
  Uvec target_height; // height above surface
  Uvec target_ampl;
  int  target_type_id;
  Uvar grid_azim_center;
  Uvar grid_azim_res;
  Uvar grid_elev_center;
  Uvar grid_elev_res;
  Uvec target_azim;
  Uvec target_elev;
  Ivec target_beam;
  Uvec target_alt;  // spacecraft altitude for which target in centered
  Ivec target_inbound;
  Ivec target_special; // flag to set targets near the boresight as special
                       // for debugging purposes
  Uvar updown_shift;
  Uvar filter_start_freq;
  Uvar filter_start_freq_offset;
  bool filter_centered;
  Uvar transmit_power;
  Uvar Tsys;
  Uvar chirp_rate;
  Uvar lambda_chirp;
  ParamSourceE radar_param_source;
  Baq  baq;
  Flyby flyby;
  MFile diag_file;
  float start_time_margin_in_s;
  float end_time_margin_in_s;
  float start_freq_margin_in_Hz;
  float end_freq_margin_in_Hz;
  float min_rtt_in_s;
  float max_rtt_in_s;
  float min_dop_in_Hz;
  float max_dop_in_Hz;

  unsigned int burst_no;
  unsigned int data_take_id;  

  // double-valued params for use in speeding up processing
  double rc_bw_in_Hz;
  double adc_in_Hz;
  double filter_start_freq_in_Hz;
  double chirp_length_in_s;
  double chirp_rate_in_Hz_per_s;
  double carrier_freq_in_Hz;
  double rwd_in_s;
  double fast_csf_in_Hz;
  double updown_shift_in_Hz;
  double pri_in_s;


  // intermediate variables
  Uvar fdop;
  Uvar delay;
  Uvar delay_previous_iteration;
  Uvar delay_difference;
  Uvar scale_uvar;
  Uvar power_receive;
  Uvar sar_cal_coeff;
  Uvar receive_range;
  Uvar transmit_range;
  double fdop_in_Hz;
  double delay_in_s;
  double scale;
  Uvar azimuth;
  Uvar elevation;
  Time dummy_time;
  DirectionVector look_direction;
  Uvar dXdA;
  Uvar noise_var;
  Ieb ieb;
  CAGTable cag_table_;

  // Debugging setup
  bool check_this_burst;
  bool perform_single_burst_check;
  int check_burst_idx;
  Uvar special_debug_pulse_duty;
  int special_debug_attr_gain_bias_db;
  int special_tro_in_pri;
  void config(Config& cfg);
  // put source controlled parameters in L1B file 
  void getRadarParams(Config& cfg);   
  void getRadarParamsFromIeb();   
  void updateRadarParams();
  void computeDoubleParams();
  void altBAEToLatLon(const int target_idx);
  void timeBAEToLatLon(const Time&, const int target_idx);
  bool checkTargetInRange(const int target_idx);
  void setChirpRate();
  void getDelayDopplerAndScale(double t_receive, const int target_idx);

  float addEcho(const int target_idx);
  void addNoise();
  void convertToByte();
  void simulateBaq();
  void chirpRate(Uvar& retval);
  void  idealDopplerRangeTrack();
  double getPulseTimeAndNumber(double  t_return_in_s, const Dvec& delay_in_s,
			     unsigned int& pulsenum, int first_sample); 
};
#endif









