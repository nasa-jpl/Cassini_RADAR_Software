//=============================================================================
// L1I.h
//
// This file contains the LIB class declarations.  The L1I class services
// Cassini RADAR L1I files.  In particular, it defines the file format,
// and provides methods to read and write the header and records.
//
// An L1I object is designed to support either reading or writing; the
// constructor sets the type at the beginning.
//
// Interface summary:
//
//-----------------------------------------------------------------------

#ifndef L1I_H
#define L1I_H

// FILE SIZE DEFINITIONS
#define L1I_HEADER_LENGTH 40
#define L1I_RECORD_LENGTH 3145740
#define L1I_DOPPLER_WIDTH 256
#define L1I_RANGE_WIDTH 4096
#define GENERIC_ARRAY_WIDTH (32*1024)
#define FFT_INTERP_SIZE 4096
#include <string>
#include "Units.h"
#include "Error.h"
#include "SimpleArray.h"
#include "L1B.h"
#include "BurstData.h"
#include "GPHData.h"
#include "Frame.h"
#include "CassiniSim.h"
#include "MFile.h"
#include "SARProcParams.h"
#include "SincTable.h"
#include "DebugInfo.h"
#include "RangeDopplerFile.h"
#include "Baq.h"
#include "RasList.h"
#include "TopoMap.h"

class BIDRFile;

class L1I : public BurstData
{

 public: 
  friend class L1B;
  enum errorE {unspecified, read_error, write_error, internal_error, 
	       unknown_parameter, bad_l1b, not_configured, invalid_time};
  class L1IError;

  //--------------
  // Constructors
  //--------------

  // HACK
  void computeFakeWindows();
  L1I(const std::string& filename, const std::string& mode);
  void allocateArrays();
  ~L1I(); 

  //-----
  // I/O
  //-----
  void readRecord();
  // routine for writing record
  void writeRecord();

  //  header handling routines
  void writeHeader(BurstData& l1b);
  void writeHeader();
  void writeHeader(int num_records);
  void rewriteHeader();
  void readHeader();
  void outputRangeDoppler(RangeDopplerFile& rdf);
  // Parameter access methods
  double getSpecialParam(const char* name);
  void enableSpecialParam(const char* name);
  void disableSpecialParam(const char* name);


  // Read information from processor config file
  void config(Config& cfg);

  //---------------------------------------------------------
  // Primary Burst Processor methods called from main routine
  //---------------------------------------------------------
  void zeroArrays();
  void computeBurstParams(const SARProcParams& spp, ofstream& afs, RasList* raslist=NULL);
  float estimateBaqInputRatio();
  float estimateBaqOutputRatio(Ivec& Thresh);
  float estimateBaqOutputRatio();
  float estimateBaqOutputRatio(const L1B& b); // modifies L1I, copies echo
  float estimateBaqOutputRatio(RasList* raslist);
  float numBlankPulsesForBAQ();
  float fracHighBit(Ivec& Thresh);
  float FHBToBAQBias(float fhb, float numblank);
  void rangeCompressHensley(const SARProcParams& spp, ofstream& afs);
  void azimuthCompress(const SARProcParams& spp, ofstream& afs);
  void computeUsableExtent(const SARProcParams& spp);
  void calibrate(SARProcParams& spp, ofstream& afs,float height_in_km=0); // spp not const because Xcal is updated
  int HISARcalibrate(SARProcParams& spp, double radii[3], double freq_shift, int ambigmode, double tol); // spp not const because Xcal is updated
  float correlateRawWithX();
  int getGoodSigma0s(float* val,int N);
  float getSARProcGain();
  void fftinterp();
  void removeDopplerPhaseRamp(const SARProcParams& spp, ofstream& afs); // spp not const because Xcal is updated


  //---------------------------------------------------
  // Routine for Computing SAR Ancilary data to go in LBDR
  //---------------------------------------------------

  void computeSARAncillaryData(const SARProcParams& spp, OblCylProj& proj);

  //---------------------------------------------------------
  // Doppler/Range to Latitude/Longitude conversion routines
  // called by Main routine or BIDR object methods
  //---------------------------------------------------------
  void  initLatLonTransformation(const SARProcParams& spp, OblCylProj& proj);
  void  getOblCylLatLon(OblCylProj& proj,double range, double doppler,double * lat, double* lon);
  void  initAzimElevTransformation();
  void computeLatLonBounds(const OblCylProj& proj);
  void reportBoundary(const OblCylProj& proj);
  bool goodSideOfNadir(double slon, double slat);
  bool goodSideOfNadir(double pos[3],double eang);
  void latLonToDopplerAndRange(double slon, 
				   double slat, double& fdop_in_Hz,
				   double& range_in_km);

  bool dopplerRangeHeightToLatLon(double fdop_in_Hz,double range_in_km, double height_in_km,double& slat, double& slon);
 
  void BAQEncodeAndDecode(RasList* raslist);
  float minimumBAQThresholdVariance(RasList* raslist);

  //--------------------------------------------------------
  // Methods for producing output for anomalous data
  //--------------------------------------------------------
  void setBadSARData();
  void fixRerouteChirpBurst(); // 10/24/2010

  //------------------------------------------------
  // Incidence Angle Computation
  //------------------------------------------------
  float getIncidenceAngle(float& slat, float& slon);


  //------------------------------------------------
  // Routines for interpolating in doppler/range space
  //------------------------------------------------
  bool  initInterpolate(double fdop_in_Hz, double range_in_km);
  float sigma0Interpolate();
  float sigma0Interpolate(double fdop_in_Hz, double range_in_km);
  bool nadirAmb(double fdop_in_Hz, double range_in_km);
  float XInterpolate();
  float XAmbigInterpolate();
  void initGeometryCalculations();

  //------------------------
  // routines for getting the true unadulterated range and Doppler
  // associated with a pixel (Range and Doppler centroid offsets removed)
  // used to produce range and doppler backplanes for stereo processing
  //------------------------
  double rawRange(double fdop_in_Hz, double range_in_km);
  double rawDoppler(double fdop_in_Hz);
  double resolutionArea(SARProcParams& spp);
  // Public booleans and enum
  bool output_scatterometer_info;
  bool quick_look_dopran;
  bool orthorectify;
  bool replace_X_mode;
  bool nadir_mode;
  enum replace_X_typeE {NONE, RANGE_RES, AZIMUTH_RES, LOOKVEL_ANGLE,
  RES_ELEMENT_ANGLE, AMBIG_RATIO, ORTHO_LOCAL_RAD, ORTHO_HEIGHT } replace_X_param;


  // Public Burst specific data fields


  float* radar_data;

  // Objects for noise subtraction
  Baq baq;
  

  // Pulse segmentation parameters
  unsigned int num_samples_per_window;
  unsigned int window_step;
  unsigned int num_range_bins;
  int start_range_bin;
  int end_range_bin;
  Uvar window_start_delay;
  unsigned int first_pulse_idx;
  unsigned int window_start_idx;
 
  // FFT sizes
  int Nfft_r; // length of FFT used in range compression
  int Nfft_az; // length of FFT used in azimuth compression

  // Quantities used in normalization of compressed image
  // Normalization Is Done By multiplying each detected pixel by
  // mean_square_inputs/effective_duty_cycle/image_sum
  float image_sum; // sum of detected compressed image
                              // pixels 
  float effective_duty_cycle; // Duty cycle of received pulses excludes 
                              // pulse spreading, includes window overlap
  float mean_square_inputs; // Mean square value of uncompressed 
  float sum_square_dechirp;
  float sum_square_fulldechirp;
  float rangeCompressGain;
  float azimuthCompressGain;
  float beam_azioff_rad[5];
  int num_signal_samples;

  float overlap_ratio;
  // pulse-segmented input image ( includes overlap between windows )

  
  // Bounding boxes
  // lonmin,lonmax,latmin,latmax
  double lonlat_bounds[4];
  double dopran_bounds[4];

  // Weighted sinc interpolation filters
  SincTable sinc_table_range;
  SincTable sinc_table_doppler;
  TopoMap ortho_tmap;
  Uvar nominal_delay;  // delay to boresight halfway between centers of Tx/Rx
  Uvar ortho_height_tolerance; // height tolerance for orthorectification

  bool use_HISAR_calibrate;
  bool exclude_mirror;

  protected: 

  // private doppler conversion/computation routines
  double absDoppler(int range_i, int doppler_j);
  double relDoppler(int j);

  // private routines called by I/o methods
  void resetRecord();
  void readSARData();
  void writeSARData();

  // private routines used by computeBurstParameters
  void computeWindows();
  void computeDoubleParams();
  void computeDopplerAndRange(const SARProcParams& spp);
  double getFastDopplerCentroid(double look[3]);
  double getFastDopplerRate(double look[3],double range_in_km, double fdop);
  int getFastLookVector(double range_in_km,float azim_rad,double look[3]);

  // private routine called by rangeCompressHensley
  void computeMatchedFilterHensley();
  
  // calibration sub-routines
  float computeGain2(double doppler, double range, const DebugInfo& dbg, 
		     float& azim_deg, float& elev_deg,double height_in_km=0);
  float getPixelAreaOnSphere(double doppler, double range, double height_in_km=0);

  float getBinAreaOnTriaxialBody(double radii[3],double azim0, double azim_res, double elev0, double elev_res);

  // ASCII dump of pixel by pixel values
  void pixelDump(OblCylProj& proj);

  // Boolean flags
  bool data_complete; // burst processing is complete

  //----------------------------------------------
  // private objects used in burst data processing
  //----------------------------------------------
  CassiniSim cassini_sim; // Object for handling beam and target geometry
                          // during a Cassini data take
 
  vector<Beam> beam_;       // Beam objects-- trying to get rid of dependency on
                       // CassiniSim

  StateVector scstate_; // State Vector object -- likewise
  StateVector scstate_tmp1_;
  StateVector scstate_tmp2_;
  FloatVector scacc_vec_; // s/c acceleration

  // Doppler and range arrays
  double* fdop_centroid;  // Hz
  double* fdop_rate; // Hz/s
  double* range;     // km

  // doppler and range step in compressed image
  double doppler_step_;
  double range_step_;

  // parameters used for estimating pixel size and resolution
  double doppler_ground_width_;
  double range_ground_width_;
  // Array used in azimuthCompress
  complex<float>* azimuth_ref_func;

  // Arrays used in rangeCompress
  complex<float>* matched_filter;
  complex<float>* fft_matched_filter;


  // generic complex<float> arrays
  complex<float>* tmp_complex_array;
  complex<float>* tmp_complex_array2;

  // output arrays
  complex<float>** dechirped;  // output of rangeCompress
  complex<float>** raw_sar_data; // output of azimuthCompress
  float** Xfactor; // output of calibrate
  float** Xambig; // output of HISARcalibrate
  float** ests0; // output of HISARcalibrate
  float** ests0ambig; // output of HISARcalibrate
  int** num_hs_bins;  // intermed values for HISARcalibrate
  int** num_hs_bins_amb;// intermed values for HISARcalibrate
  int HISAR_ambig_mode;
  float signal_to_amb_thresh;
  float min_oneway_gain;
  Uvar nadir_exclude_ang;
  complex<float>** sigma0; // output of calibrate
  complex<float>** interp_sigma0; // output of fftinterp
  bool fft_interpd_;
  double* fftinterp_fdop_centroid;
  float  fftinterp_range_step;
  float  fftinterp_doppler_step;
  float  fftinterp_range_min;
  float  fftinterp_reldop_min;
  int fftinterp_num_doppler_bins;
  int fftinterp_num_range_bins;
  int  fft_range_idx_;
  int  fft_doppler_idx_;
  float** area_buf;
  float** g2_buf;
  float** dop_buf;
  float** azim_buf;
  float** elev_buf;
  double nadir_range_in_km;
  double nadir_fdop_in_Hz;
  double range_window_in_km;
  // intermediate data file pointers
  FILE* sigma0fp; // sigma0
  FILE* is0fp; // interpolated sigma0
  FILE* Xfp; // Xfactor
  FILE* Xambfp; // Xambig
  FILE* dcfp; // dechirped data
  FILE* rfp; // raw_sar_data
  FILE* mffp; // matched filter (reference function for range compression)
  FILE* g2fp; // gain squared data
  FILE* areafp; // area data
  FILE* rangefp; // 1-D range data
  FILE* fdopfp; // 2-D doppler data
  FILE* estrdfp; // ASCII estimated doppler/range offset file
  FILE* azimfp; // 2-D antenna azimuth in degrees files
  FILE* elevfp; // 2-D antenna elevation in degrees files
  FILE* pixdfp; // ASCII pixel by pixel dump file
  FILE* rspecfp; // Range spectrum dump file
  double pixd_noise_factor; // estimate of uncalibrated detected SAR data for noise only
  
  // Skip parameters for estimating range/doppler
  int estrd_iskip, estrd_nskip, estrd_dopwid, estrd_ranwid;
  // Parameters specifying doppler and

  Uvar nominal_range; // range to boresight halfway between centers of Tx/Rx
  Uvar nominal_doppler_centroid;
  Uvar nominal_doppler_rate;
  Uvar nominal_pulse_spread;
  Uvar rc_time_shift;
  Uvar ac_dop_shift;
  double rc_range_shift_in_km,ac_dop_shift_in_Hz;
  Uvar ddopdrange;
  Uvar ddopratedrange;
  int nominal_range_idx;

  // Burst specific values
  Time ground_impact_time;// midpoint between center of transmit gate
                          // and center of receive gate
  Uvar lambda_chirp; // wavelength of center of transmitted spectrum

  // Instrument Parameters -- should probably move to SARProcParams
  Uvar updown_shift;
  Uvar filter_start_freq;
  Uvar filter_start_freq_offset;
  bool filter_centered;
  Uvar transmit_power;

  // SAR processing parameters and options -- should definitely move to
  // SARProcParams
  double range_ref_window_param;
  double azimuth_ref_window_param;
  double pulse_spread_in_pri;
  int DC_null_bins;
  bool use_upper_band;
  bool use_linear_chirp;
  bool use_trivial_pulse_segmentation_;
  bool disable_warnings_;


  //--------------------------------------
  // Individual Burst Checking parameters
  //--------------------------------------
  bool perform_single_burst_check;  // flag individual burst checking on
  bool check_this_burst; // flag burst checking on for current burst
  int check_burst_idx; // burst number to check
  int burst_no; // current burst number
  string check_burst_filename; // filename to output burst check info
  MFile cbf; // Object which handles writes to burst check file


  // Intermediate variable for antenna gain and pixel area calculations
  double rotmat_[3][3]; // Rotation matrix from Target Body Fixed to 
                        // Beam Frame at ground impact time
                        // this is the inverse of the rotation stored
                        // in the intermediate geometry file
  double rotmat_rev_[3][3]; // Beam frame to target body fixed;


  double scpos_[3]; // s/c position in TBF (km)
  double scvel_[3]; // s/c velocity in TBF (km/s)
  double scacc_[3]; // s/c velocity in TBF (km/s/s)

  double target_radius_; // radius of spherical (assumed) target (km)
  double position_in_radius_; // magnitude of scpos_ in units of target radii
  double scspeed_; // s/c speed in km/s
  
  double boresight_in_tbf_[3]; // beam boresight direction in Target Body Fixed
  double borepos_[3]; // ground location of boresight

  // double versions of Uvars
  double rc_bw_in_Hz;
  double adc_in_Hz;
  double chirp_length_in_s;
  double chirp_rate_in_Hz_per_s;
  double default_target_radii_in_km[3];
  // Correction term for error in range compression due to doppler
  // range(true)=range(SAR)+RDCC*(fdop-nominal_doppler_centroid)
  double range_doppler_correction_coeff; 


  double rwd_in_s;
  double fast_csf_in_Hz;
  double updown_shift_in_Hz;
  double pri_in_s;
  double prf_in_Hz;
  double slow_cfs_in_Hz;
  double csd_in_s;
  double fdop_max;
  double fdop_min;
  double range_min;
  double range_max;  
  double lambda_chirp_in_km;
  double lambda_special;  // effective lambda that account for special relativity when used in fast look vector computation
  double chirp_freq_in_Hz;  
  double lambda_in_km;
  double ground_impact_time_in_s; 
  double nominal_doppler_centroid_in_Hz;

  // interpolation parameters
  int current_range_idx_;
  int current_doppler_idx_;
  int use_sinc_interpolation_;  // 0 => use nearest neighbor
                                // else use sinc interpolation
  int use_fft_interpolation_;

  // bounds on indices of samples to be used in sinc interpolation calculations
  int min_range_idx_;
  int max_range_idx_;
  int min_doppler_idx_;
  int max_doppler_idx_;
  // buffers for interpolation in range and Doppler
  complex<float> *interp_range_buf;
  complex<float> *interp_doppler_buf;


  // data and routines for realap SARTopo processing
 public:
  float ambigfrac;
  float** s0frombidr;
  float means0;
  void getBIDRSigma0s(BIDRFile& bf);
  void realap_sartopo_calibrate(SARProcParams& spp, ofstream& afs, float h);  
  void sartopo_calibrate(SARProcParams& spp, ofstream& afs, float h, bool init_mode=false);  
  float getPowerRMSError(SARProcParams& spp,float h,bool dbg);
  void allocateRealAppSartopoArrays();
};

//------------------------------------------------------
// Error handling class for L1B (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class L1I::L1IError : public BurstDataError
  {
  public:

  // Constructors

    L1IError(errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified L1I Error";
    else if (error_type == read_error)
      msg = "L1I read error";
    else if (error_type == write_error)
      msg = "L1I write error";
    else if (error_type == internal_error)
      msg = "L1I internal error";
    else if (error_type == unknown_parameter)
      msg = "L1I unknown parameter name";
    else if (error_type == bad_l1b)
      msg = "L1I attempted to convert bad L1B";
    else if (error_type == not_configured)
      msg = "L1I not yet configured";
    else if (error_type == invalid_time)
      msg = "L1I time error";
    }

    L1IError(const std::string& emsg, errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    msg = emsg;
    }

    void throwMe() // No exceptions
      {
	throw *this;
      }
  // Public type flag
  L1I::errorE error_type;
  };



#endif


