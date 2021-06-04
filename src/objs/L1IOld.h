//=============================================================================
// L1IOld.h
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

#ifndef L1I_OLD_H
#define L1IOld_H

// FILE SIZE DEFINITIONS
#define L1IOld_HEADER_LENGTH 40
#define L1IOld_RECORD_LENGTH 2008944
#define L1IOld_DOPPLER_WIDTH 100
#define L1IOld_RANGE_WIDTH 1000

#include <string>
#include "Units.h"
#include "Error.h"
#include "Array.h"
#include "L1B.h"
#include "BurstData.h"
#include "GPHData.h"
#include "Frame.h"
#include "CassiniSim.h"
#include "MFile.h"

class L1IOld : public BurstData
{

 public: 
  enum errorE {unspecified, read_error, write_error, internal_error, 
	       unknown_parameter, bad_l1b, not_configured, invalid_time};
  class L1IOldError;

  //--------------
  // Constructors
  //--------------

  
  L1IOld(const std::string& filename, const std::string& mode);

  ~L1IOld(); 

  //-----
  // I/O
  //-----

  void config(Config& cfg);

  // methods for assigning data to record
  void convertL1B(L1B& l1b);
  void readRecord();


  // routine for writing record
  void writeRecord();

  //  header handling routines
  void writeHeader(BurstData& l1b);
  void writeHeader();
  void writeHeader(int num_records);
  void rewriteHeader();
  void readHeader();

  // Parameter access methods
  double L1IOld::getSpecialParam(const char* name);
  void L1IOld::enableSpecialParam(const char* name);
  void L1IOld::disableSpecialParam(const char* name);

  // Public L1IOld specific data fields
  float rms_radar_data;
  unsigned long int num_samples_per_window;
  unsigned long int window_step;
  unsigned long int num_range_bins;
  int Nfft;
  unsigned long int num_pulses_received;
  unsigned long int first_pulse_idx;
  Fmat raw_sar_data;
  Umat Xfactor;
  Umat sigma0;
  Umat range;
  Umat doppler;
  Uvar window_start_delay;
  unsigned long int window_start_idx;
  CFvec matched_filter_time;
  complex<float>* matched_filter_hensley;
  complex<float>* fft_matched_filter_hensley;
  complex<float>* tmp_complex_array;
  complex<float>* tmp_complex_array2;
  fdata doppler_filter_freq;
  Uvar fraction_energy_in_window;

  protected: 

  void resetRecord();
  void readSARData();
  void writeSARData();
  void computeWindows();
  void rangeCompressHensley();
  void computeMatchedFilter();
  void computeMatchedFilterHensley();
  void computeDoubleParams();
  void computeDopplerAndRange();
  void calibrate();


  bool estimate_range_for_windowing_;
  bool estimate_doppler_centroid_;
  bool use_match_filter_weights_;
  bool use_doppler_filter_weights_;
  bool data_complete_;
  unsigned long int Nradar_data;
  fdata  radar_data; 
  CassiniSim cassini_sim;
  Uvar nominal_doppler_centroid;
  Uvar nominal_delay;
  Uvar nominal_pulse_spread;
  Uvar updown_shift;
  Uvar filter_start_freq;
  Uvar filter_start_freq_offset;
  Uvar transmit_power;
  Uvar antenna_impedance;
  double range_ref_window_param;
  int DC_null_bins;
  bool use_hensley_range_compress;
  bool use_upper_band;
  bool use_linear_chirp;
  bool perform_single_burst_check;
  bool check_this_burst;
  int check_burst_idx;
  int burst_no;
  string check_burst_filename;
  MFile cbf;
  Uvar system_loss;
  CFmat dechirped;

  // double versions of Uvars
  double rc_bw_in_Hz;
  double adc_in_Hz;
  double chirp_length_in_s;
  double rwd_in_s;
  double fast_csf_in_Hz;
  double updown_shift_in_Hz;
  double pri_in_s;
  double slow_cfs_in_Hz;
  double csd_in_s;
};

//------------------------------------------------------
// Error handling class for L1B (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class L1IOld::L1IOldError : public BurstDataError
  {
  public:

  // Constructors

    L1IOldError(errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified L1IOld Error";
    else if (error_type == read_error)
      msg = "L1IOld read error";
    else if (error_type == write_error)
      msg = "L1IOld write error";
    else if (error_type == internal_error)
      msg = "L1IOld internal error";
    else if (error_type == unknown_parameter)
      msg = "L1IOld unknown parameter name";
    else if (error_type == bad_l1b)
      msg = "L1IOld attempted to convert bad L1B";
    else if (error_type == not_configured)
      msg = "L1IOld not yet configured";
    else if (error_type == invalid_time)
      msg = "L1IOld time error";
    }

    L1IOldError(const std::string& emsg, errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    msg = emsg;
    }

    void throwMe() // No exceptions
      {
	throw *this;
      }
  // Public type flag
  L1IOld::errorE error_type;
  };

#endif


