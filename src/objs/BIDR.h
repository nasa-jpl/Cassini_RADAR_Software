#ifndef BIDR_H
#define BIDR_H

static const char rcs_id_bidr_h[] =
  "@(#) $Id: BIDR.h,v 11.10 2013/10/15 23:02:19 bstiles Exp $";

enum BIDRModeE {NORMAL, STEREO, TOPO, TOPO_PLUS, INVALID_MODE};
enum BIDRTypeE {INC, BEAMMASK, S0_UNCORR, S0_CORR, LAT, LON, START_BURST_NUM, END_BURST_NUM, NUM_LOOKS, S0_CORR_DB, DOPPLER, RANGE, NOMINAL_BURST_NUM, S0NSQT_UNCORR, S0_STD, S0_NOISE_EQUIV, INVALID};

#define BIDRTYPECOUNT 17
#define ISIS_NULL "16#FF7FFFFB#"
#define BAD_VALUE_HEX 0xff7ffffb
#define BITS_PER_BYTE 8
#define MAX_BIDR_LINES 50000

#include "PDSLabel.h"
#include "L1I.h"
#include "SARProcParams.h"
#include "Time.h"
#include "Io.h"
#include "LatLonGrid.h"
#include "Projections.h"

enum corrTypeE {VENUS, DBLINE041104, RKIRK01TA, HAGHAG, ENCELADUS_WYE1, RHEA_WYE1, DIONE_WYE1, INVALID_CORR};
  float incAngleCorrect(float s0_uncorr, float inc, float slat, float slon, corrTypeE inc_correction_model);
float theorBackscatter(float inc,corrTypeE inc_correction_model);

class BIDRLabel : public PDSLabel
{
 private:
  int dummy_param;
};

class BIDR
{
 public:
  BIDR(const string& prefix,
       Config& cfg, L1B& l1b, SARProcParams& spp, const string& mode, BIDRModeE proc_mode); 
  ~BIDR();
  void stripProcessAndWrite(L1I& l1i, SARProcParams& spp, bool single_burst_mode=false);
  void topoProcessAndWrite(L1I& l1i, SARProcParams& spp, bool single_burst_mode=false, float trial_height=0);
  void smoothOverlap(int** col, int** width);
  void setOverlap(int** col, int** width);
  void flush(const SARProcParams& spp);
  void outputLine(const SARProcParams& spp);
  void outputBeamDeltas(const SARProcParams& spp);
  OblCylProj& getProjection(){return(proj);};
  void writeHeaders();
  void rewriteHeaders();
  bool longitudeLineComplete(int i, const SARProcParams& spp);
  void  multipleBeamMerge(const SARProcParams& spp, int i, int j, 
			  unsigned char& beam_mask, float& s0_uncorr, 
                          float& inc, int& burstno1,int& burstno2, 
			  int& numlooks, float& X, float& s0ns, float& s0ne,
			  float& doppler, float& range);
  void  updateTopoPlusArrays(const SARProcParams& spp, int i, int j, int pn, int abs_i, int col); 
  int checksum(const BIDRTypeE bidr_type);
  static double positiveWestLon(double d) { return (360.0 - d); };
  static float positiveWestLon(float d) { return ((float) (360.0 - d)); };



  // topo specific public buffers
  float** dsigma0;
  float** topo_inc1;
  float** topo_inc2;
  float** topo_lat;
  float** topo_lon;
  float** topo_wnl;
  float** mid_ds0;
  float** slope_ds0;
  float** rms_ds0;
  float** power1;
  float** power2;
  float** noisefloor1;
  float** noisefloor2;
  float** ambrat1;
  float** ambrat2;
  float** sumsigma0;
  float** dsigma0std;
  float** dsigma0bias;
  int** overlap_col;
  int** overlap_width;
  float*** topoplus_s01;
  float*** topoplus_s02;
  float*** topoplus_s0ne1;
  float*** topoplus_s0ne2;
  float*** topoplus_s0baqscale1;
  float*** topoplus_s0baqscale2;
  float*** topoplus_s0amb1;
  float*** topoplus_s0amb2;
#define MAX_TOPO_OVERLAP_WIDTH 700

  // common public buffers
  float* area_pixel;
 private:

  // private methods
  void initLongitudeLine(int i);
  void allocateArrays(const SARProcParams& spp);
  void updateTopoArrays(int first_i);
  void openFiles(const string& prefix, const string& mode,int enable_noise_sub=0);
  void closeFiles();
  void incrementBurstsSinceUpdate();
  void writeHeader(FILE* f, BIDRTypeE e);
  void rewriteHeader(FILE* f, BIDRTypeE e);
  string constructLabel(BIDRTypeE e);
  double mapScale();
  string productID(const BIDRTypeE bidr_type);
  string filePrefix(const string filename);

  // sar stereo processing specific methods
  void initProcMode();
  corrTypeE inc_correction_model;
  // fields
  BIDRModeE   procMode;
  BIDRLabel   label;
  OblCylProj  proj; // Contains parameter values, routines for getting
                    // latitude and longitude from range, doppler, and time
  LatLonGrid grid; // includes size of Grid starting lat and lon and steps 
              // method for going from lat/lon to indices

  float*** inc_bfr;
  float*** s0_uncorr_bfr; 
  float*** s0_nsub_uncorr_bfr; 
  float*** s0ne_bfr; 
  float*** X_bfr;
  float*** max_X_bfr;
  float*** dop_bfr;
  float*** range_bfr;
  float*** resfactor;
  float*** noise_equiv_s0;
  float*** Xambig_bfr;
  int* bursts_since_update;
  bool** beams_found;
  int*** start_burst_num;
  int*** end_burst_num;
  int*** nom_burst_num;
  int*** num_looks;
  
  double** slon_bfr;
  double** slat_bfr;
  
  // The incidence angle corrected sigma0 is computed just prior to output
  // And Beam mask/selection/weighting is performed just prior to output
  // And geographic latitude longitude values computed just prior to output
  // 1-D arrays are needed to enable efficient I/O
  char* beammask_io_bfr;
  float* s0_corr_io_bfr;
  float* lat_io_bfr;
  float* lon_io_bfr;
  float* inc_io_bfr;
  float* X_io_bfr;
  float* std_io_bfr;
  float* s0ne_io_bfr;
  float* dop_io_bfr;
  float* range_io_bfr;
  float* s0_uncorr_io_bfr; 
  int* start_burst_num_io_bfr;
  int* end_burst_num_io_bfr;
  int* nom_burst_num_io_bfr;
  unsigned char* num_looks_io_bfr;

 

  // file pointers FILE* is used instead of FileMgr in order to
  // speed up I/O and allow flexibility in which files to open and write

  FILE* inc_file;
  FILE* beammask_file;
  FILE* s0_uncorr_file;
  FILE* s0_corr_file;
  FILE* lat_file;
  FILE* lon_file;
  FILE* start_burst_num_file;
  FILE* end_burst_num_file;
  FILE* num_looks_file;
  FILE* X_file;
  FILE* std_file;
  FILE* s0ne_file;
  FILE* beam_deltas_file;
  FILE* dop_file;
  FILE* range_file;
  FILE* nom_burst_num_file;
  FILE* topoplus_file;
  // options should probably go to SARProcParams
  bool use_oblique_cyl; // True:= outputs on Oblique Cyl Grid
                        // False:= output on standard lat/lon grid
  bool ignore_gaps_mode;
  bool allow_decreasing_lon_mode;
  bool replace_X_mode;
  bool feather_beams;
  bool overlap_set;
  string replace_X_param;
  bool dominant_beam3;
  bool mark_nadir_amb;
  bool output_beam_deltas;
  float replace_X_value;
  int pixelsperdegree;
  int num_lines_of_output;
  int num_looks_required;
  int max_looks_allowed;
  int num_looks_to_skip;
  bool enable_max_looks;
  int num_beams;
  int Xcheckdir;
  Time sclk_start_;
  Time sclk_stop_;
  int sclk_stop_raw_;
  
  // PDS label values
  static string pds_version_id_;
  static string record_type_;

  static string data_set_id_;
  static string data_set_name_;
  static string producer_institution_name_;
  static string producer_id_;
  string producer_full_name_;
  int data_take_id_;
  int product_version_;
  static string instrument_host_name_;
  static string instrument_host_id_;
  static string instrument_name_;
  static string instrument_id_;
  string target_name_;
  string source_product_id_;
  string mission_phase_name_;
  static string mission_name_;
  string software_version_;

  static string sample_type_[];
  bool output_enable_[BIDRTYPECOUNT];
  static int sample_bits_[];
  int *checksum_; // don't use float; more precision needed
  static double scaling_factor_;
  static double offset_;
  static string missing_constant_[];
  static string note_[];

  string data_set_map_projection_catalog_; // name of the PDS volume catalog
                                           // describing the map projection,
                                           // not the name of the map projection
  static string map_projection_type_;
  static string first_standard_parallel_;
  static string second_standard_parallel_;
  static string positive_longitude_direction_;
  static double proj_center_latitude_;
  static double proj_center_longitude_;
  static int line_first_pixel_;
  static int sample_first_pixel_;
  float maximum_latitude_;  // degrees
  float minimum_latitude_;  // degrees
  float maximum_longitude_; // degrees, positive-east (not -west)
  float minimum_longitude_; // degrees, positive-east (not -west)
  float first_longitude_;   // degrees, positive-east (not -west)

  static double map_projection_rotation_;
  string look_direction_;
  static string coordinate_system_type_;
  static string coordinate_system_name_;


  // bad value float
  float BAD_VALUE;

  // INCIDENCE ANGLE CORRECTION PARAMETERS
  
  float venus_k1, venus_k2, venus_s045;
  float dbline_offset, dbline_slope, dbline_s045;
  float invsin_coeff,invsin_s045;
  

  bool auto_overwrite_bidr_; // true => overwrite existing BIDR output files
                             // false => if BIDR output files exist, ask user
			     // to confirm overwrite

  bool check_pds_string_lengths_;
  bool use_config_data_take_number_; // true => read data take number from
                                     // config file
                                     // false => read number from LBDR file
  string flyby_id_pass_;
  int segment_id_;
};

#endif
