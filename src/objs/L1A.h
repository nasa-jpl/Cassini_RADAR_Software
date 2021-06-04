//=============================================================================
// L1A.h
//
// This file contains the L1A class declarations.  The L1A class services
// Cassini RADAR L1A files.  In particular, it defines the file format,
// and provides methods to read and write the header and records.
//
// An L1A object is designed to support either reading or writing; the
// constructor sets the type at the beginning.
//
// Interface summary:
//
// L1A methods:
//
//   Construction:
//
//   L1A(filename,"x"); Construct for reading (x = rb) or writing (x = wb)
//
//   I/O:
//
//   loadSab(sab);  Copies relevant data from sab to this L1A object
//   writeHeader(int num_records=0); Write header to binary file
//   readHeader();  Read header from binary file
//   writeRecord(); Write current record to binary file
//   readRecord();  Read record from binary file
//   show_passive_data_on_screen();   
//
//   Predicates:
//
//   eof();     queries status of associated file
//   isParam(); checks for valid parameter name
//
//   Parameter Access:
//
//   getParam(); Fetch the specified parameter by name
//
//   Other methods:
//
//   recordCount();  Return the number of L1A records in the associated file.
//
// L1A public Sab data:
//
// The data fields in the record are directly accessible.
// Some are UnitVar's which are coerced to particular units whenever they
// are read or written.
// The method getParam() provides another way to access the data fields.
//
// L1A Error handling:
//   L1A::L1AError(err_type);  Construct an L1A error handling object
//=============================================================================

#ifndef L1A_H
#define L1A_H

// FILE SIZE DEFINITIONS
#define L1A_HEADER_LENGTH 40
#define L1A_PASSIVE_RECORD_LENGTH 568

#include <string>
#include "Units.h"
#include "Sab.h"
#include "Error.h"
#include "Array.h"
#include "BurstData.h"

class L1A 
: public BurstData
  {
  public:

  enum errorE {unspecified, read_error, write_error, internal_error, 
	       unknown_parameter,incomplete_sab, empty_sab, undecoded_sab};

  class L1AError; 
  friend class L1B; // allows L1B to read L1A file directly

  //--------------
  // Constructors
  //--------------

  L1A(const std::string& filename, const std::string& mode,
      const std::string& filetype);

  ~L1A(); // destructor rewrites header if mode_ is "wb"

  //-----
  // I/O
  //-----

  void loadSab(const Sab& sab);
  void readRecord();
  void writeRecord(int abs_record_no=-1);
  void writeHeader();
  void rewriteHeader();
  void readHeader();
  void show_passive_data_on_screen();

  // parameter access

  double getSpecialParam(const char* name);
  void enableSpecialParam(const char* name);
  void disableSpecialParam(const char* name);
  
  // special Parameter calculation routines
  Uvar normalizeCntRL();
  Uvar normalizeCntND();
  Uvar normalizeCntRadio();

  // Public radar data Fields
  fdata  decoded_data;
  unsigned int Ndecoded;


  private:

  //-------------------
  // Internal variables
  //-------------------

  bool radar_data_loaded_;
  unsigned  int Ndat_max_;


  };

//------------------------------------------------------
// Error handling class for L1A (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class L1A::L1AError : public BurstDataError
  {
  public:

  // Constructors

    L1AError(errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified L1A Error";
    else if (error_type == incomplete_sab)
      msg = "L1A tried to load incomplete Sab";
    else if (error_type == empty_sab)
      msg = "L1A tried to load empty Sab";
    else if (error_type == undecoded_sab)
      msg = "L1A tried to load undecoded Sab";
    else if (error_type == read_error)
      msg = "L1A read error";
    else if (error_type == write_error)
      msg = "L1A write error";
    else if (error_type == internal_error)
      msg = "L1A internal error";
    else if (error_type == unknown_parameter)
      msg = "L1A unknown parameter name";
    }

    L1AError(const std::string& emsg, errorE err_type = unspecified) // No exceptions
    : error_type(err_type)
    {
    msg = emsg;
    }

    void throwMe() 
      {
	throw *this;
      }
  // Public type flag
  L1A::errorE error_type;
  };

#endif









