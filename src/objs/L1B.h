//=============================================================================
// L1B.h
//
// This file contains the L1B class declarations.  The L1B class services
// Cassini RADAR L1B files.  In particular, it defines the file format,
// and provides methods to read and write the header and records.
//
// An L1B object is designed to support either reading or writing; the
// constructor sets the type at the beginning.
//
// Interface summary:
//
//-----------------------------------------------------------------------

#ifndef L1B_H
#define L1B_H

static const char rcs_id_l1b_h[] =
  "@(#) $Id: L1B.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

// FILE SIZE DEFINITIONS
#define L1B_HEADER_LENGTH 40
#define L1B_PASSIVE_RECORD_LENGTH 1272

#include <string>
#include <map>
#include "Units.h"
#include "Sab.h"
#include "Error.h"
#include "Array.h"
#include "BurstData.h"
#include "GPHData.h"
#include "Frame.h"

// Other constants
#define AUTO_GAIN_MASK 0xfff7

enum BODPTypeE { SBDR, LBDR, ABDR, INVALID_BODP };
#define BODPTYPECOUNT 4

class L1I;
class L1B : public BurstData
{ 
  public:

  enum errorE {unspecified, read_error, write_error, internal_error, 
	       unknown_parameter, incomplete_sab, empty_sab, undecoded_sab,
	       not_configured, file_not_read, geom_not_computed,
	       invalid_time};

  class L1BError; 

  //--------------
  // Constructors
  //--------------

  
  L1B(const std::string& filename, const std::string& mode,
      const std::string& filetype, bool abdr_flag=false);

  ~L1B();

  //-----
  // I/O
  //-----

  void config(Config& cfg);

  // methods for assigning data to record   
  void loadSab(const Sab& sab);
  void setNum_bursts_in_flight(const unsigned int& num_bursts);
  void locate(const Umat& azim_1way3dB_ellipse_fit, 
	      const Umat& elev_1way3dB_ellipse_fit,
	      const Umat& azim_2way3dB_ellipse_fit, 
	      const Umat& elev_2way3dB_ellipse_fit);
  void copy(L1B& l1b);
  void copyPassive(L1B& l1b);
  void setDataTakeID(int id){data_take_id_=id; BurstData::data_take_id_=id;} 
  // Only used for special cases of corrupt data_take
  // for example the stuff_enceladus.cpp program

  // method for reading and writing record (even creating)
  void createRecord(const Time& t0);
  void readRecord();
  void writeRecord(int abs_record_no=-1);

  // header handling routines
  void writeHeader();
  void rewriteHeader();
  void readHeader();

  // routine for outputting SAR ancillary data from L1I object to L1B file
  void outputSARAncillaryData(const L1I& l1i); 

  // read one record to a L1I or (BurstData object, an array and int)
  void copyTo(L1I& l1i);
  void copyTo(BurstData& bd, float* r, unsigned int& N);
  void copyTo(BurstData& bd, fdata& r, unsigned int& N);
  void checkNradar_data();
  
  // set values needed for PDS labels
  void recordStartTime(Time& t);
  void recordStopTime(Time& t);
  void recordRadarMode(unsigned int csr, Uvar adc);

  // PDS label and PDS value handling routines
  void readPdsValuesFromLabel();
  void copyPdsValues(const L1B& dat);
  void copyPdsEndValues(const L1B& dat);
  void combineRecordCounts(const L1B& dat1, const L1B& dat2);

  // various aceess method  to data record
 
  void  mapRecords();//must be done before callling the following functions
  bool foundSABRecord(const unsigned int& sab_counter); //access file poistion via sab counter
  void findSABCounter_at_T(const Time& t, unsigned int& sab_counter,
		     bool& foundRecord);//access file position via time_from_ca
  unsigned int& operator[] (const unsigned int& n);//access sab counter via record number: 0 to Nrec-1
  void getEcho(const unsigned int& sab_counter, 
	       const unsigned int& Nradar_data, 
	       fdata& echo,
	       bool& valid);
  void getEcho_fast(const unsigned int& sab_counter, 
	       const unsigned int& Nradar_data, 
	       float* echo,
	       bool& valid);
  string productID();

  // public active mode data
  fdata radar_data;

  protected:

  //------------------
  // protected methods
  //------------------

  void resetRecord();

  string constructLabel(BODPTypeE bodp_type);
  string createProductID(BODPTypeE bodp_type);

  // parameter access

  double getSpecialParam(const char* name);
  void enableSpecialParam(const char* name);
  void disableSpecialParam(const char* name);  

  //-------------------
  // Internal variables
  //-------------------

  bool abdr_flag_;
  bool data_complete_;
  bool  map_complete_;
  unsigned  int Nradar_data_max_;

  //---------------------
  //record access: sab counter vs file position
  //               time_from_ca vs file position
  //               record id vs sab counter
  //---------------------
  map<unsigned int, unsigned int> record_counter_map_;
  map<unsigned int, unsigned int>::const_iterator record_counter_map_ptr_;

  map<unsigned int, int> sab_counter_map_;
  map<unsigned int, int>::const_iterator sab_counter_map_ptr_;

  map<Time, unsigned int> time_map_;
  map<Time, unsigned int>::const_iterator time_map_ptr_;;
 

  //---------------------------------------
  // PDS label values and related variables
  //---------------------------------------

  Time sclk_start_;
  Time sclk_stop_;

  static int RADIOMETER_MODE;
  static int SCATTEROMETER_MODE;
  static int ALTIMETER_MODE;
  static int SAR_MODE;

  static string bodp_type_name_[];

  static string pds_version_id_;
  static string record_type_;

  static string data_set_version_id_;
  static string data_set_id_[];
  static string data_set_name_[];
  static string producer_institution_name_;
  static string producer_id_;
  string producer_full_name_;
  int data_take_id_;
  int product_version_;
  int pds_radar_mode_;
  static string instrument_host_name_;
  static string instrument_host_id_;
  static string instrument_name_;
  static string instrument_id_;
  static string mission_name_;
  string software_version_;
  static string description_prefix_[];
  string flyby_id_;
  Time pds_trigger_time_;
  static string processing_history_text_[];

  static string interchange_format_;
  static string table_description_[];

  static int table_column_count_[];

  bool check_pds_string_lengths_;
  bool config_set_;

  string product_id_;
};




//------------------------------------------------------
// Error handling class for L1B (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class L1B::L1BError : public BurstDataError
  {
  public:

  // Constructors

    L1BError(errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified L1B Error";
    else if (error_type == read_error)
      msg = "L1B read error";
    else if (error_type == write_error)
      msg = "L1B write error";
    else if (error_type == internal_error)
      msg = "L1B internal error";
    else if (error_type == unknown_parameter)
      msg = "L1B unknown parameter name";
    else if (error_type == incomplete_sab)
      msg = "L1B attempted to load incomplete SAB";
    else if (error_type == empty_sab)
      msg = "L1B attempted to load empty SAB";
    else if (error_type == undecoded_sab)
      msg = "L1B attempted to load undecoded SAB";
    else if (error_type == not_configured)
      msg = "L1B not yet configured";
    else if (error_type == file_not_read)
      msg = "L1B has not yet read data file";
   else if (error_type == geom_not_computed)
      msg = "L1B has not yet computed geometry info";
   else if (error_type == invalid_time)
      msg = "L1B time error";
    }

    L1BError(const std::string& emsg, errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    msg = emsg;
    }

    void throwMe() // No exceptions
      {
	throw *this;
      }
  // Public type flag
  L1B::errorE error_type;
  };

#endif









