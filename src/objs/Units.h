//==============================================================================
// Units.h
//
// This file contains the interface classes and functions for automatic
// unit handling.  The public interface consists of the class
// UnitVar and some associated functions.  This header comment summarizes
// the interface.  For details about a specific function, look at the
// declarations in this file.
//
// The unit tables are held internally and managed by classes UnitTable
// and UnitDerivedTable.  The user should not need access to these, however,
// the conversion factors and recognized unit terms can be seen in the
// initializers for static data for these two classes in Units.cpp.
//
// Interface summary:
//
// class UnitVar:
//   Allows units to be attached to a double precision float.
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
//   fabs(); compute absolute value.
//   negate(); negate the magnitude
//   asin(); compute asin(*this) handling units appropriately
//   acos(); compute acos(*this) handling units appropriately
//   atan(); compute atan(*this) handling units appropriately
//   atan2(); compute atan(*this) handling units appropriately
//   
//   Mathematical methods that do not modify *this:
//
//   s = sin(); compute sin(*this) handling units appropriately
//   s = cos(); compute sin(*this) handling units appropriately
//   s = tan(); compute sin(*this) handling units appropriately
//
//   Other methods:
//   double coerce(ustr); Convert units to match ustr, also return the value.
//   double getInUnits(ustr); Return value after matching units to ustr.
//   double getValue();  Just return the value ignoring any units.
//   double setValue(const double& d);  Set the value ignoring any units.
//   string getUnits();  Return string representation of *this's units.
//
// UnitVar Error Handling: Unit::UnitError
//
// Other supporting functions:
//
// Binary operators for UnitVar: +, -, *, /
// ASCII stream I/O operators: <<, >>
// Other overloaded functions: pow(a,b), exp, log, trig, invtrig
//==============================================================================

#ifndef Units_H
#define Units_H

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_h[] =
  "@(#) $Id: Units.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <string>
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <functional>
#include <map>
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

class UnitVar
  {
  //----------------------------------------------------------------
  // Static Interface - Class services as opposed to object services
  //----------------------------------------------------------------

  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  //------------------------
  // General class behavior
  //------------------------

  static void setMode(const string& mode_str);

  //---------
  // Testing
  //---------

  static bool selfTest(const string& filename);

  private:

  //-----------------------
  // Static implementation
  //-----------------------

  static bool automatic_;
  static bool enforce_base_units_;
  static bool no_unit_support_;


  //-------------------------
  // Normal object interface
  //-------------------------

  public:

  // Friend classes that need access to UnitVar's representation
  friend class Units;
  friend class Unit::UnitError;

  //--------------
  // Construction
  //--------------

  UnitVar();
  UnitVar(const double& d);
  UnitVar(const string& namestr);
  UnitVar(const double& d, const Units& u);
  UnitVar(const double& d, const string& ustr);
  UnitVar(const string& name, const double& d, const std::string& ustr);
  UnitVar(const FileMgr& file, int usize = -1);
  UnitVar(const UnitVar& uv); // copy constructor, no exceptions
  UnitVar(const UnitVar& uv, const string& ustr); // coercion constructor
  ~UnitVar();

  //------------------------------------------------------------
  // Conversion functions that need access to UnitVar internals
  //------------------------------------------------------------

  friend double convert_to_double(double d);
  friend double convert_to_double(const UnitVar& u);
  friend int convert_to_int(const UnitVar& u);
  friend double get_in_units(const UnitVar& uv, const string& ustr);
  
  //----------------
  // Operators
  //----------------

  UnitVar operator-() const ; // no exceptions

  UnitVar& operator+=(const UnitVar& arg);
  UnitVar& operator-=(const UnitVar& arg);
  UnitVar& operator*=(const UnitVar& arg);
  UnitVar& operator/=(const UnitVar& arg);
  UnitVar& operator*=(const double& arg);
  UnitVar& operator/=(const double& arg);

  // Assignment
  UnitVar& operator=(const UnitVar& uv); // no exceptions
  UnitVar& operator=(const double& b);

  //------------
  // Predicates
  //------------

  bool operator==(const UnitVar& b) const;
  bool operator!=(const UnitVar& b) const;
  bool operator<(const UnitVar& b) const;
  bool operator<=(const UnitVar& b) const;
  bool operator>(const UnitVar& b) const;
  bool operator>=(const UnitVar& b) const;

  bool hasUnits() const ; // no exceptions
  bool hasNonSuppressedUnits() const ; // no exceptions
  double matchConversion(const UnitVar& b) const;
  double matchConversionAndResetUnits(const UnitVar& b);

  //-------------------------------------------------------------
  // I/O
  // Also see the appropriate constructor above.
  //-------------------------------------------------------------

  void read(const FileMgr& file, int usize = -1);
  void write(const FileMgr& file, int usize = -1) const;

  void read(const FileMgr& file, int usize, const string& ustr);
  void write(const FileMgr& file, int usize, const string& ustr);


  void readFloat(const FileMgr& file, int usize);
  void writeFloat(const FileMgr& file, int usize);

  void readFloat(const FileMgr& file, int usize,
		 const std::string& ustr);
  void writeFloat(const FileMgr& file, int usize,
		  const std::string& ustr);

  friend std::ostream& operator<<
    (std::ostream& s, const UnitVar& uvar);
  friend std::istream& operator>>
    (std::istream& s, UnitVar& uvar);

  // These two are for bcc32 which can't handle the << >> operators correctly.
  std::ostream& streamOut(std::ostream& s) ; // no exceptions
  std::istream& streamIn(std::istream& s) ; // no exceptions
  void show() const;

  // Read an entire UnitVar from a string.
  void stringIn(const string& uvar_str) ;
  // Return strings
  string stringOut();
  string stringOut(const std::string& ustr);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  // pow
  UnitVar& pow(const UnitVar& b) ;
  UnitVar& pow(const double& b) ;
  UnitVar& pow(const int& n) ; // no exceptions

  // sqrt
  UnitVar& sqrt() ;

  // exp and log
  UnitVar& exp() ;
  UnitVar& log() ;

  // abs and negate
  UnitVar& fabs() ; // no exceptions
  UnitVar& negate() ; // no exceptions

  // Trigonometry
  double sin() const ;
  double cos() const ;
  double tan() const ;
  UnitVar& asin() ;
  UnitVar& acos() ;
  UnitVar& atan() ;
  UnitVar& atan2(const UnitVar& b) ;

  //---------------
  // Other methods
  //---------------

  double coerce(const string& ustr);
  double coerce(const UnitVar& uv);
  double coerceBase();
  double getInUnits(const string& ustr) const;
  double getInUnits(const UnitVar& uv) const;
  double getValue() const; // no exceptions
  void setValue(const double& d); // no exceptions
  string getUnits() const; // no exceptions
  void setkm(double d);
  double km() const;

  string name() const; // no exceptions
  UnitVar& floor(const double& a);
  UnitVar& floor(const UnitVar& a);
  UnitVar& ceiling(const double& a);
  UnitVar& ceiling(const UnitVar& a);

  private:

  // Internal methods
  double coerce(const Units& u);

  // Internal representation
  string* name_;
  double a_;
  Units* u_;
  static Units* dummy_units();
  void setUnitsPointer();
  };


//------------------------------------------------------
// Error handling class for UnitVar (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class Unit::UnitError : public ErrorMessage
  {
  public:

  typedef std::ostringstream OSTRINGSTREAM;

  // Constructors

  // should not throw any exceptions
  UnitError(errorE err_type = unspecified) 
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

  // should throw any exceptions 
  UnitError(const std::string& emsg, errorE err_type = unspecified) 
    : error_type(err_type)
    {
    msg = emsg;
    }

  UnitError(const UnitVar& a, const UnitVar& b)
    : error_type(mismatch)
    {
    msg = "UnitVar Mismatch: (name=" + a.name() + ",units=" +
      a.u_->toString() + "), (name=" + b.name() + ",units=" + 
      b.u_->toString() + ")";
    }

  void throwMe()
    {
    throw *this;
    }

  // Public type flag
  errorE error_type;
  };

//------------------------------
// Binary Operators for UnitVar.
//------------------------------

// +
UnitVar operator+(UnitVar a, const UnitVar& b);

// -
UnitVar operator-(UnitVar a, const UnitVar& b);

// *
UnitVar operator*(UnitVar a, const UnitVar& b);

// /
UnitVar operator/(UnitVar a, const UnitVar& b);

// >>
std::istream& operator>>(std::istream& s, UnitVar& uvar);

// <<
std::ostream& operator<<(std::ostream& s, const UnitVar& uvar);

//------------------------------------
// Other overloaded functions
//------------------------------------

// pow
UnitVar pow(UnitVar a, const double& b);
double pow(const double& a, const UnitVar& b);

// sqrt
UnitVar sqrt(UnitVar a);

// exp
double exp(UnitVar a);

// log
double log(UnitVar a);

// fabs
UnitVar fabs(UnitVar a);

// sin
double sin(UnitVar a);

// cos
double cos(UnitVar a);

// tan
double tan(UnitVar a);

// asin
UnitVar asin(UnitVar a);

// acos
UnitVar acos(UnitVar a);

// atan
UnitVar atan(UnitVar a);

// atan2
UnitVar atan2(UnitVar a, const UnitVar& b);

// floor
UnitVar floor(UnitVar a, const double& b);

// floor
UnitVar floor(UnitVar a, const UnitVar& b);

// ceiling
UnitVar ceiling(UnitVar a, const double& b);

// ceiling
UnitVar ceiling(UnitVar a, const UnitVar& b);

// overloaded function to convert UnitVar (or double) to a double
// This is useful for templates in which it is unknown whether a variable
// is a double or a UnitVar, so UnitVar::getValue() cannot be used

double convert_to_double(double d);
double convert_to_double(const UnitVar& u);
int convert_to_int(const UnitVar& u);

UnitVar coerce_base(const double& d, const string& ustr);
UnitVar coerce(const double& d, const string& ustr);
UnitVar coerce(UnitVar uv, const string& ustr);
double get_in_base_units(UnitVar uv);
double get_in_base_units(const double& d);
double get_in_units(const UnitVar& uv, const string& ustr);

typedef UnitVar Uvar;


#endif
