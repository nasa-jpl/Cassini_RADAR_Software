//==============================================================================
// Units.h
//
// This file contains the interface classes and functions for automatic
// unit handling.  The public interface consists of the templatized class
// UnitVar and some associated functions.  This header comment summarizes
// the interface.  For details about a specific function, look at the
// declarations in this file.
//
// The setup function specifyUnitConversions(filename) must be called before
// any unit conversions can be handled.  The filename refers to an ascii file
// containing the various unit conversions to be used.
//
// Interface summary:
//
// template<class T_numtype> class UnitVar:
//   Allows units to be attached to a numeric type <T_numtype>.
//   These units are automatically merged/converted with other units
//   when arithmetic involving UnitVar's is performed.  Unit mis-matches
//   cause exceptions to be thrown.  Units are usually assigned by the
//   constructor, but they can be added later.  Unit conversions can be
//   forced when desired.
//
// UnitVar methods:
//
//   Construction:
//
//   UnitVar();          Construct with no units and value of 0.
//   UnitVar(name);      Construct with no units and value of 0 and given name.
//   UnitVar(d);         Construct with no units and value of d.
//   UnitVar(d,ustr)     Construct with units ustr and value d.
//   UnitVar(name,d,ustr)Construct with units ustr and value d and given name.
//   UnitVar(uvar_str);  Construct from string representation uvar_str.
//   UnitVar(file);      Construct from binary file
//
//   Testing:
//   selfTest(unit_conversion_filename);  Test UnitVar functionality
//
//   Setup:
//   specifyUnitConversions(filename); Read external table of unit conversions.
//
//   Unary operators: +=, -=, *=, /=, -
//
//   Unary Predicates: hasUnits()
//   Binary Predicates: ==, !=, <, <=, >, >=
//
//   I/O:
//   void write(file,usize); Write binary file with 3 options for units.
//   void read(file,usize);  Read binary file with 3 options for units.
//   void read(file,usize,ustr);  Read binary file and coerce units to ustr.
//   void readFloat(file,usize);  Read binary file with 3 options for units.
//   void readFloat(file,usize,ustr);  Read binary file and coerce units.
//   void writeFloat(file,usize);  Write binary file with 3 options for units.
//   void writeFloat(file,usize,ustr);  Write binary file and coerce units.
//
//   Mathematical methods that modify *this and return *this:
//
//   pow(b); raise *this to b power handling units appropriately
//   sqrt(); compute sqrt(*this) handling units appropriately
//   exp(); compute exp(*this) handling units appropriately
//   log(); compute log(*this) handling units appropriately
//   abs(); compute absolute value.
//   negate(); negate the magnitude
//   asin(); compute asin(*this) handling units appropriately
//   acos(); compute acos(*this) handling units appropriately
//   atan(); compute atan(*this) handling units appropriately
//   
//   Mathematical methods that do not modify *this:
//
//   s = sin(); compute sin(*this) handling units appropriately
//   s = cos(); compute sin(*this) handling units appropriately
//   s = tan(); compute sin(*this) handling units appropriately
//
//   Other methods:
//   T_numtype coerce(ustr); Convert units to match ustr, also return the value.
//   T_numtype getInUnits(ustr); Return value after matching units to ustr.
//   T_numtype getValue();  Just return the value ignoring any units.
//   string getUnits();  Return string representation of *this's units.
//
// UnitVar Error Handling: Unit::UnitError
//
// Other supporting functions:
//
// Binary operators for UnitVar: +, -, *, /
// ASCII stream I/O operators: <<, >>
// Other overloaded template functions: pow(a,b)
//==============================================================================

#ifndef Units_H
#define Units_H

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_h[] =
  "@(#) $Id: TUnits.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <string>
#include <iostream>
#include <strstream>
#include <list>
#include <vector>
#include <functional>
//#include <algorithm>
#include <map>
//#include <pair.h>
//#include <utility.h>
#include <math.h>
#include "Error.h"
#include "Io.h"

//-------------------------------------------------------
// Supporting declarations (not part of public interface)
//-------------------------------------------------------

#include "Units_impl1.h"

//---------------------------------------------------------
// Class UnitVar (Primary user interface for unit handling)
//---------------------------------------------------------

template<class T_numtype>
class UnitVar
  {
  public:

  //--------------
  // Construction
  //--------------

  UnitVar();
  UnitVar(const T_numtype& d);
  UnitVar(const string& namestr);
  UnitVar(const T_numtype& d, const std::string& ustr) throw(Unit::UnitError);
  UnitVar(const string& name, const T_numtype& d, const std::string& ustr)
    throw(Unit::UnitError);
  UnitVar(const FileMgr& file, int usize = -1)
    throw(Unit::UnitError, FileMgr::IoError);
  UnitVar(const UnitVar<T_numtype>& uv) throw(); // copy constructor
  ~UnitVar();

  //---------
  // Testing
  //---------

  static bool selfTest(const std::string& filename);

  //------
  // Setup
  //------

  static void specifyUnitConversions(const std::string& filename)
    throw(Unit::UnitError);

  //----------------
  // Operators
  //----------------

  UnitVar<T_numtype> operator-() const throw();

  UnitVar<T_numtype>& operator+=(const UnitVar<T_numtype>& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator-=(const UnitVar<T_numtype>& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator*=(const UnitVar<T_numtype>& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator/=(const UnitVar<T_numtype>& arg)
    throw(Unit::UnitError);

  UnitVar<T_numtype>& operator+=(const T_numtype& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator-=(const T_numtype& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator*=(const T_numtype& arg)
    throw(Unit::UnitError);
  UnitVar<T_numtype>& operator/=(const T_numtype& arg)
    throw(Unit::UnitError);

  // Assignment
  UnitVar<T_numtype>& operator=(const UnitVar<T_numtype>& uv) throw();

  //------------
  // Predicates
  //------------

  bool operator==(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool operator!=(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool operator<(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool operator<=(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool operator>(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool operator>=(const UnitVar<T_numtype>& b) const throw(Unit::UnitError);
  bool hasUnits() const throw();

  //-------------------------------------------------------------
  // I/O
  // Also see the appropriate constructor above.
  //-------------------------------------------------------------

  //-------------------------------------------------------------------------
  // write(file,usize), read(file,usize)
  //
  // Write a UnitVar to a file and read a UnitVar from a file.
  // The numeric component is written/read first,
  // followed by the unit string.
  //
  // file = FileMgr object to write to
  // usize = the number of characters to write/read for the unit string.
  //   If usize is 0, then no units are written/read.
  //   If usize is -1, then the actual size of the unit string comes before
  //      the string itself.  The length of the unit string is written as a
  //      unsigned char, followed by the string itself (no trailing 0).
  //      If the unit string exceeds 255 chars, then an exception is thrown.
  //      (this is the default value if usize is not explicitely passed)
  //   If usize is positive, then precisely usize characters are written
  //      after the size field, and no trailing 0.  If the unit string is
  //      larger than usize chars, then an exception is thrown.
  //      Excess characters are padded with spaces.
  //-------------------------------------------------------------------------

  void write(const FileMgr& file, int usize = -1) const
    throw(Unit::UnitError, FileMgr::IoError);
  void read(const FileMgr& file, int usize = -1)
    throw(Unit::UnitError, FileMgr::IoError);

  //--------------------------------------------------------------------------
  // readFloat(file,usize)
  //
  // This method reads a float from the file regardless of the underlying
  // type of UnitVar.  If the UnitVar is an int, then a conversion from float
  // to int will occur.
  //--------------------------------------------------------------------------

  void UnitVar<T_numtype>::readFloat(const FileMgr& file, int usize)
    throw(Unit::UnitError, FileMgr::IoError);
  void UnitVar<T_numtype>::writeFloat(const FileMgr& file, int usize)
    throw(Unit::UnitError, FileMgr::IoError);

  // For the read option below, the UnitVar is coerced to have the units
  // contained in ustr after being read from the file.
  void read(const FileMgr& file, int usize, const std::string& ustr)
    throw(Unit::UnitError, FileMgr::IoError);
  // See readFloat above
  void UnitVar<T_numtype>::readFloat(const FileMgr& file, int usize,
    const std::string& ustr)
    throw(Unit::UnitError, FileMgr::IoError);
  void UnitVar<T_numtype>::writeFloat(const FileMgr& file, int usize,
    const std::string& ustr)
    throw(Unit::UnitError, FileMgr::IoError);

  friend std::ostream& operator<< <>
    (std::ostream& s, const UnitVar<T_numtype>& uvar);
  friend std::istream& operator>> <>
    (std::istream& s, UnitVar<T_numtype>& uvar);

  // These two are for bcc32 which can't handle the << >> operators correctly.
  std::ostream& streamOut(std::ostream& s) throw();
  std::istream& streamIn(std::istream& s) throw();

  // Read an entire UnitVar from a string.
  void stringIn(const string& uvar_str) throw(Unit::UnitError);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  // pow
  UnitVar<T_numtype>& pow(const UnitVar<T_numtype>& b) throw(Unit::UnitError);
  UnitVar<T_numtype>& pow(const T_numtype& b) throw(Unit::UnitError);
  UnitVar<T_numtype>& pow(const int& n) throw();

  // sqrt
  UnitVar<T_numtype>& sqrt() throw(Unit::UnitError);

  // exp and log
  UnitVar<T_numtype>& exp() throw(Unit::UnitError);
  UnitVar<T_numtype>& log() throw(Unit::UnitError);

  // abs and negate
  UnitVar<T_numtype>& abs() throw();
  UnitVar<T_numtype>& negate() throw();

  // Trigonometry
  T_numtype sin() const throw(Unit::UnitError);
  T_numtype cos() const throw(Unit::UnitError);
  T_numtype tan() const throw(Unit::UnitError);
  UnitVar<T_numtype>& asin() throw(Unit::UnitError);
  UnitVar<T_numtype>& acos() throw(Unit::UnitError);
  UnitVar<T_numtype>& atan() throw(Unit::UnitError);

  //---------------
  // Other methods
  //---------------

  T_numtype coerce(const std::string& ustr) throw(Unit::UnitError);
  T_numtype getInUnits(const std::string& ustr) const throw(Unit::UnitError);
  T_numtype getValue() const throw();
  string getUnits() const throw();
  string name() const throw();

  protected:

  // Internal representation
  string name_;
  T_numtype a_;
  Units u_;

  };

//------------------------------------------------------
// Error handling class for UnitVar (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class Unit::UnitError : public ErrorMessage
  {
  public:

  // Constructors

  UnitError(errorE err_type = unspecified) throw()
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified UnitVar Error";
    else if (error_type == read_error)
      msg = "UnitVar read error";
    else if (error_type == write_error)
      msg = "UnitVar write error";
    else if (error_type == internal_error)
      msg = "UnitVar internal error";
    else if (error_type == mismatch)
      msg = "UnitVar unit mismatch error";
    else if (error_type == formation)
      msg = "UnitVar unit formation error";
    else if (error_type == reduction_error)
      msg = "UnitVar unit reduction error";
    }

  UnitError(const std::string& emsg, errorE err_type = unspecified) throw()
    : error_type(err_type)
    {
    msg = emsg;
    }

  // Public type flag
  errorE error_type;
  };

//------------------------------
// Binary Operators for UnitVar.
//------------------------------

/**
// +
template<class T_numtype>
UnitVar<T_numtype>
  operator+(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError);

// +
template<class T_numtype, class T_numtype2>
UnitVar<T_numtype>
  operator+(const UnitVar<T_numtype>& a, const T_numtype2& b)
  throw(Unit::UnitError);

// +
template<class T_numtype, class T_numtype2>
UnitVar<T_numtype>
  operator+(const T_numtype2& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError);

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const UnitVar<T_numtype>& a, const float& b)
  throw(Unit::UnitError);

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const float& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const UnitVar<T_numtype>& a, const double& b)
  throw(Unit::UnitError);

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const double& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError);

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const float& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const UnitVar<T_numtype>& a, const float& b)
  throw(Unit::UnitError);

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const double& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const UnitVar<T_numtype>& a, const double& b)
  throw(Unit::UnitError);

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError);

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const UnitVar<T_numtype>& a, const float& b)
  throw(Unit::UnitError);

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const float& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const UnitVar<T_numtype>& a, const double& b)
  throw(Unit::UnitError);

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const double& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);
**/

// >>
template<class T_numtype>
std::istream& operator>>(std::istream& s, UnitVar<T_numtype>& uvar);

// <<
template<class T_numtype>
std::ostream& operator<<(std::ostream& s, const UnitVar<T_numtype>& uvar);

//------------------------------------
// Other overloaded template functions
//------------------------------------

// pow
template<class T_numtype, class T_exptype>
UnitVar<T_numtype>
  pow(const UnitVar<T_numtype>& a, const T_exptype& b)
  throw(Unit::UnitError);

// sqrt
template<class T_numtype>
UnitVar<T_numtype> sqrt(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// exp
template<class T_numtype>
T_numtype exp(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// log
template<class T_numtype>
T_numtype log(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// abs
template<class T_numtype>
UnitVar<T_numtype> abs(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// sin
template<class T_numtype>
T_numtype sin(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// cos
template<class T_numtype>
T_numtype cos(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// tan
template<class T_numtype>
T_numtype tan(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// asin
template<class T_numtype>
UnitVar<T_numtype> asin(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// acos
template<class T_numtype>
UnitVar<T_numtype> acos(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// atan
template<class T_numtype>
UnitVar<T_numtype> atan(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError);

// The remaining declarations and required definitions are placed in a separate
// file.  They are not considered part of the normal interface.
#include "Units_impl2.h"

typedef UnitVar<double> Uvar;

#endif
