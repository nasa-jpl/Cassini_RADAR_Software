//----------------------------------------------------------------------------
// Units.cpp
//
// This file contains method definitions for the classes supporting UnitVar
// and supporting function definitions.  Template definitions need to be
// included in the .h file, so they are put in Units_impl2.h.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_c[] =
  "@(#) $Id: Units.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <iostream>
#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include "Units.h"
#include "Io.h"
#include "Constants.h"
#include "Utils.h"
#include <cmath>

using std::list;
using std::map;
using std::string;
using std::istream;
using std::ostream;
using std::ifstream;
using std::min;
using std::cerr;
using std::cout;
using std::endl;
using std::streampos;

//----------------------
// Class UnitVar Methods
//----------------------

//----------------------
// Static initialization
//----------------------

bool UnitVar::automatic_ = true;
bool UnitVar::enforce_base_units_ = false;
bool UnitVar::no_unit_support_ = false;




//-------------------------
// Static interface methods
//-------------------------

void UnitVar::setMode(const string& mode_str)
  {
  if (mode_str == "automatic")
    {
    automatic_ = true;
    enforce_base_units_ = false;
    no_unit_support_ = false;
    }
  else if (mode_str == "enforce_base_units")
    {
    automatic_ = false;
    enforce_base_units_ = true;
    no_unit_support_ = false;
    }
  else if (mode_str == "no_unit_support")
    {
    automatic_ = false;
    enforce_base_units_ = false;
    no_unit_support_ = true;
    }
  else
    {
    cerr << "UnitVar::setMode: Invalid mode string: " + mode_str << endl;
    exit(-1);
    }
  }

//------------------------------------------------------------------------
// Initialize Units pointer u_
//--------------------------------------------------------------------
void UnitVar::setUnitsPointer()
{
  if(no_unit_support_) u_=dummy_units();
  else{
    u_= new Units;
  }
}
//---------------------------------------------------------------
// Static routine called to intialize pointer in no_unit_support case
//--------------------------------------------------------------
Units* UnitVar::dummy_units()
{
  static Units* p = new Units();
  return(p);
}

//-----------------------
// Normal object methods
//-----------------------

//--------------
// Construction
//--------------


UnitVar::UnitVar()
  : name_(NULL), a_(0) { 
  setUnitsPointer();
}

//-------------------------------------------------------------------------
// UnitVar(name)
//
// Construct from a name string.  If the name is a recognized constant,
// then the constant value (see below) is set along with the name.
// If the name is unknown, then only the name is set, while the UnitVar
// is initialized by default to 0.0 without units.
//-------------------------------------------------------------------------

UnitVar::UnitVar(const string& namestr)
  : name_(NULL)
  {
  name_=new string;
  *name_=namestr;
  setUnitsPointer();
  if (namestr == "speed_light")
    {
    stringIn("2.997924580e5 km/s");
    }
  else if (namestr == "boltzmann_constant")
    {
    stringIn("1.3806503e-23 J/K");
    }
  else if (namestr == "gravitational_constant")
    {
    stringIn("6.673e-20 km km km/(s s kg)");
    }
  else if (namestr == "free_space_permittivity")
    {
    stringIn("8.854187818e-3 uF/km");
    }
  }

// Conversion constructors
UnitVar::UnitVar(const double& d)
  : name_(NULL), a_(d) {
  setUnitsPointer();
}

UnitVar::UnitVar(const double& d, const Units& u)
  : name_(NULL), a_(d)
  {
  setUnitsPointer();
  if (!no_unit_support_)
    {
    *u_ = u;
    }
  }

UnitVar::UnitVar(const double& d, const string& ustr)
  : name_(NULL), a_(d)
  {
  setUnitsPointer();
  if (!no_unit_support_)
    {
    *u_ = Units(ustr);
    if (enforce_base_units_ && *u_ != u_->baseUnits())
      {
      Unit::UnitError e("UnitVar(" + toStr(d) + "," + ustr +
        "): Not in base units", Unit::mismatch);
      e.throwMe();
      }
    }
  }

UnitVar::UnitVar(const string& namestr, const double& d, const string& ustr)
  : name_(NULL), a_(d)
  {
  name_=new string;
  *name_=namestr; 
  setUnitsPointer();
  if (!no_unit_support_)
    {
    *u_ = Units(ustr);
    if (enforce_base_units_ && *u_ != u_->baseUnits())
      {
      Unit::UnitError e(namestr + "(" + toStr(d) + "," + ustr +
        "): Not in base units", Unit::mismatch);
      e.throwMe();
      }
    }
  }

// copy constructor
// should not throw any exceptions
UnitVar::UnitVar(const UnitVar& uv) 
  : name_(NULL), a_(uv.a_)
  {
  setUnitsPointer();
  if (!no_unit_support_) *u_ = *(uv.u_);
  }

//------------------------------------------------------------------------
// UnitVar(uv,ustr)
//
// Coercion constructor.  Builds a copy of uv coerced to have the units
// specifed by ustr.  Just like coerce except that a UnitVar is returned
// instead of a double.
// In enforce_base_unit mode, ustr must be base units.  In this case,
// the only call that would have any effect would be if uv has no units
// in which case the specified units are attached.
//------------------------------------------------------------------------

UnitVar::UnitVar(const UnitVar& uv, const string& ustr) 
  : name_(NULL), a_(uv.a_)
  {
  setUnitsPointer();
  if (!no_unit_support_) *u_ = *(uv.u_);
  coerce(ustr);
  }

// destructor
UnitVar::~UnitVar()
  { 
    if(name_) delete name_;
    if(!no_unit_support_) delete u_;
  }

//---------------------------
// Construct from binary file
//---------------------------

UnitVar::UnitVar(const FileMgr& file, int usize)
  : name_(NULL)
  {
  // Order must match write below.
  setUnitsPointer();
  file.read(a_);
  u_->read(file,usize);
  }

//----------------------------------------------------------
// selfTest
//
// Exercize the UnitVar type and verify proper operation.
// Actual work done in Units::selfTest.
//----------------------------------------------------------

bool UnitVar::selfTest(const string& filename)
  {
  return(Units::selfTest(filename));
  }

//-----------------
// Member Operators
//-----------------

// unary -
// should not throw exceptions
UnitVar UnitVar::operator-() const
  {
  UnitVar val = *this;
  val.a_ = -val.a_;
  return(val);
  }

// +=
UnitVar& UnitVar::operator+=(const UnitVar& arg)
  {
  if (no_unit_support_)
    {
    a_ += arg.a_;
    return(*this);
    }
  double c = matchConversionAndResetUnits(arg);
  a_ += c * arg.a_;
  return(*this);
  }


// -=
UnitVar& UnitVar::operator-=(const UnitVar& arg)
  { 
  if (no_unit_support_)
    {
    a_ -= arg.a_;
    return(*this);
    }
  double c = matchConversionAndResetUnits(arg);
  a_ -= c * arg.a_;
  return(*this);
  }

// *=
UnitVar& UnitVar::operator*=(const double& arg)
  {
  a_ *= arg;
  return(*this);
  }

// *=
UnitVar& UnitVar::operator*=(const UnitVar& arg)
  { 
  if (no_unit_support_)
    {
    a_ *= arg.a_;
    return(*this);
    }
  // Merge arg's units with *this's units.
  double c = u_->factorMerge(*(arg.u_),Units::product);
  a_ *= c * arg.a_;
  return(*this);
  }

// /=
UnitVar& UnitVar::operator/=(const double& arg)
  {
  a_ /= arg;
  return(*this);
  }

// /=
UnitVar& UnitVar::operator/=(const UnitVar& arg)
  {
  if (no_unit_support_)
    {
    a_ /= arg.a_;
    return(*this);
    }
  // Merge arg's units with *this's units.
  double c = u_->factorMerge(*(arg.u_),Units::quotient);
  a_ *= c / arg.a_;
  return(*this);
  }

// Assignment
// should not throw any exceptions
UnitVar& UnitVar::operator=(const UnitVar& arg) 
  {
  if (no_unit_support_)
    {
    a_ = arg.a_;
    return(*this);
    }
  // Self assignment is harmless.
  // Assignment leaves *this's name unchanged, but copies everything else.
  a_ = arg.a_;
  *u_ = *(arg.u_);
  return(*this);
  }

// Assignment
UnitVar& UnitVar::operator=(const double& arg)
  {
  // does not reassign name.
  if (no_unit_support_)
    {
    a_ = arg;
    return(*this);
    }
  a_=arg;
  *u_= Units("");
  return(*this);
  }

//------------
// Predicates
//------------

// ==
bool UnitVar::operator==(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ == b.a_);
    }
  double c = matchConversion(b);
  return(a_ == c * b.a_);
  }

// !=
bool UnitVar::operator!=(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ != b.a_);
    }
  return(!(UnitVar::operator==(b)));
  }

// <
bool UnitVar::operator<(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ < b.a_);
    }
  double c = matchConversion(b);
  return(a_ < c * b.a_);
  }

// >
bool UnitVar::operator>(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ > b.a_);
    }
  double c = matchConversion(b);
  return(a_ > c * b.a_);
  }

// >=
bool UnitVar::operator>=(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ >= b.a_);
    }
  return(!(UnitVar::operator<(b)));
  }

// <=
bool UnitVar::operator<=(const UnitVar& b) const
  {
  if (no_unit_support_)
    {
    return(a_ <= b.a_);
    }
  return(!(UnitVar::operator>(b)));
  }


// hasUnits
bool UnitVar::hasUnits() const 
  {
  if (no_unit_support_)
    {
    Unit::UnitError e("UnitVar::hasUnits: N/A when unit support is disabled");
    e.throwMe();
    }
  return(u_->hasUnits());
  }

//--------------------------------------------------------------------------
// hasNonSuppressedUnits()
//
// Returns true if this UnitVar has any unit terms besides the suppressed
// unit term.
//--------------------------------------------------------------------------

bool UnitVar::hasNonSuppressedUnits() const 
  {
  if (no_unit_support_)
    {
    Unit::UnitError
      e("UnitVar::hasNonSuppressedUnits: N/A when unit support is disabled");
    e.throwMe();
    }
  return(u_->hasNonSuppressedUnits());
  }

//------------------------------------------------------------------
// UnitVar matchConversion routine
// Calls Units::matchConversion and implements special case in which
// one operand (either this or b) is a zero with no units
//-------------------------------------------------------------------

 
double UnitVar::matchConversion(const UnitVar& b) const
  {
  double conversion_factor=1.0;
  bool this_has_units = hasUnits();
  bool b_has_units = b.hasUnits();
  // Nominal case:: both operands have units
  if (this_has_units && b_has_units)
    {
    try
      {
      return(u_->matchConversion(*(b.u_)));
      }
    catch(Unit::UnitError& e)
      {
      Unit::UnitError e2(*this,b);
      e2.throwMe();
      }
    }
  else if (!this_has_units && !b_has_units)
    {
    // trivial case: Neither operand has units
    conversion_factor=1.0;
    }
  else if (!this_has_units && a_==0)
    {
    // special case 1: this is the ZERO UnitVar a_=0, no units
    conversion_factor=1.0;
    }
  else if (!b_has_units && b.a_==0)
    {
    // special case 2: b is the ZERO UnitVar b.a_=0, no units
    conversion_factor=1.0;
    }
  else if (u_->hasNonSuppressedUnits() || b.u_->hasNonSuppressedUnits())
    {
    // erroneous case: One operand has at least one non-suppressed unit,
    // the other doesn't, neither is ZERO.
    Unit::UnitError e2(*this,b);
    e2.throwMe();
    }

  return(conversion_factor);
  }

//----------------------------------------------------------------
// Non-const version of matchConversion
// The units of *this are  set to the units of b if and only if
// *this is zero with no units. This routine is used for += and -= , etc
// but not for < > >= etc.
//------------------------------------------------------------------

double UnitVar::matchConversionAndResetUnits(const UnitVar& b) 
  {
  double c = matchConversion(b);
  if(!hasUnits() && a_==0)
    {
    *u_ = *(b.u_);   
    }
  return(c);
  }

//-----
// I/O
//-----

//-------------------------------------------------------------------------
// write(file,usize), read(file,usize)
//
// Write a UnitVar to a file.  The numeric component is written first,
// followed by the unit string.
//
// file = FileMgr object to write to
// usize = the number of characters to write for the unit string.
//   If usize is 0, then no units are written.
//   If usize is -1, then the actual size of the unit string is written.
//      The length of the unit string is written as a unsigned char, followed
//      by the string itself (no trailing 0).  If the unit string exceeds 255
//      chars, then an exception is thrown.
//      (this is the default value if usize is not explicitely passed)
//   If usize is positive, then precisely usize characters are written
//      after the size field, and no trailing 0.  If the unit string is
//      larger than usize chars, then an exception is thrown.
//      Excess characters are padded with spaces.
//-------------------------------------------------------------------------

void UnitVar::write(const FileMgr& file, int usize) const
  {
  if (UnitVar::no_unit_support_ && usize != 0)
    {
    Unit::UnitError e("Can't auto write units without unit support: " +
      name() + ".write(" + file.name() + "," + toStr(usize) + ")",
      Unit::mismatch);
    e.throwMe();
    }
  else if (UnitVar::enforce_base_units_)
    {
    if (*u_ != u_->baseUnits())
      {
      string ustr;
      u_->toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".write(" + file.name() + "," + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }

  file.write(a_);
  u_->write(file,usize);
  }

void UnitVar::read(const FileMgr& file, int usize)
  {
  file.read(a_);
  u_->read(file,usize);

  if (!UnitVar::automatic_ && usize != 0)
    {  // check to see if units just read are base units
    if (*u_ != u_->baseUnits())
      {
      string ustr;
      u_->toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".read(" + file.name() + "," + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }

  }

//--------------------------------------------------------------------------
// readFloat(file,usize), writeFloat(file,usize)
//
// These methods read/write a float from/to the file regardless of the
// underlying type of UnitVar.  If the UnitVar is an int,
// then a conversion from/to float to/from int will occur.
//--------------------------------------------------------------------------

void UnitVar::readFloat(const FileMgr& file, int usize)
  {
  float a;
  file.read(a);
  a_ = (double)a;    // convert float to underlying UnitVar type
  u_->read(file,usize);

  if (!UnitVar::automatic_ && usize != 0)
    {  // check to see if units just read are base units
    if (*u_ != u_->baseUnits())
      {
      string ustr;
      u_->toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".readFloat(" + file.name() + "," + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }

  }

void UnitVar::writeFloat(const FileMgr& file, int usize)
  {
  if (UnitVar::no_unit_support_ && usize != 0)
    {
    Unit::UnitError e("Can't auto write units without unit support: " +
      name() + ".writeFloat(" + file.name() + "," + toStr(usize) + ")",
      Unit::mismatch);
    e.throwMe();
    }
  else if (UnitVar::enforce_base_units_)
    {
    if (*u_ != u_->baseUnits())
      {
      string ustr;
      u_->toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".writeFloat(" + file.name() + "," + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }

  float a = (float)a_;  // convert underlying UnitVar type to float
  file.write(a);
  u_->write(file,usize);
  }

//--------------------------------------------------------------------------
// read(file,usize,ustr),write(file,usize,ustr)
//
// These read/write options will attempt to coerce the UnitVar to have the
// units in ustr (after being read from the file, or before being written
// to the file).
//--------------------------------------------------------------------------

void UnitVar::read(const FileMgr& file, int usize,
  const std::string& ustr)
  {
  file.read(a_);
  u_->read(file,usize);
  coerce(ustr);
  }

void UnitVar::write(const FileMgr& file, int usize, const string& ustr)
  {
  // Convert units first
  double a = coerce(ustr);
  file.write(a);
  u_->write(file,usize);
  }

// See readFloat above
void UnitVar::readFloat(const FileMgr& file, int usize,
  const std::string& ustr)
  {
  float a;
  file.read(a);
  a_ = a;
  u_->read(file,usize);
  coerce(ustr);
  }

// See writeFloat above
void UnitVar::writeFloat(const FileMgr& file, int usize,
  const std::string& ustr)
  {
  // Convert units and convert underlying UnitVar type to float
  float a = (float)coerce(ustr);
  file.write(a);
  u_->write(file,usize);
  }

// streamIn
std::istream& UnitVar::streamIn(std::istream& s)
  {
  double a;
  if (!(s >> a)) return s;
  Units u;
  // Units don't have to be present, hence no return here.
  if (!(s >> u)) s.clear();
  a_ = a;
  *u_ = u;
  return s;
  }

// streamOut
// should not throw any exceptions
std::ostream& UnitVar::streamOut(std::ostream& s)
  {
  s << a_ << " " << *u_;
  return s;
  }

//------------------------------------------------
// stringIn(uvar_str)
//
// Build from string (including value and units).
//------------------------------------------------

void UnitVar::stringIn(const std::string& uvar_str)
  {
  ISTRINGSTREAM uvar_str_stream(uvar_str.c_str());
  if (!(uvar_str_stream >> *this))
    {
    Unit::UnitError e("Bad unit number in string", Unit::formation);
      e.throwMe();
    }
  }

//-------------------------------------------------------------------
// str = stringOut()
// str = stringOut(ustr)
//
// Return string representation after coercing to indicated unit str
// if present.
//-------------------------------------------------------------------

string UnitVar::stringOut()
  {
  OSTRINGSTREAM uvar_str_stream;
  uvar_str_stream << *this;
  return(uvar_str_stream.str());
  }

string UnitVar::stringOut(const std::string& ustr)
  {
  OSTRINGSTREAM uvar_str_stream;
  coerce(ustr);
  uvar_str_stream << *this;
  return(uvar_str_stream.str());
  }

void UnitVar::show() const
  {
  cout << *this << endl;
  }

//--------------------------
// Other arithmetic methods
//--------------------------

// pow

UnitVar& UnitVar::pow(const UnitVar& b)
  {
  if (no_unit_support_)
    {
    a_ =  std::pow(a_,b.a_);
    return(*this);
    }

  if (b.hasUnits())
    {
    Unit::UnitError e("Exponent can't have units in " +
      name() + ".pow(" + b.name() + ")", Unit::mismatch);
    e.throwMe();
    }
  return(pow(b.a_));  // call number pow()
  }

UnitVar& UnitVar::pow(const double& b)
  {
  if (no_unit_support_)
    {
    a_ =  std::pow(a_,b);
    return(*this);
    }

  int n = (int)b;
  if (b - n == 0.0)
    {  // number^integer
    return(pow(n));  // call integer pow()
    }
  else if (hasUnits())
    {
    OSTRINGSTREAM os;
    os << "Can't raise unit to non-integer power in "
       << name() << ".pow(" << b << ")";
    Unit::UnitError e(toStr(os), Unit::mismatch);
    e.throwMe();
    }
  else
    {  // Apply regular pow function (hopefully defined for double)
    a_ = std::pow(a_,b);
    }
  return(*this);
  }

UnitVar& UnitVar::pow(const int& n)
  {  // number^integer
  if (no_unit_support_)
    {
    a_ =  std::pow(a_,n);
    return(*this);
    }

  a_ = std::pow(a_,n);
  if (n == 0)
    {  // A power of 0 means a unitless result.
    u_->clear();
    }
  else
    {  // Replicate the units n times.
    Units u(*u_);
    for (int i=1; i < std::abs(n); ++i) u_->combine(u);
    }
  // If the power is negative then need to invert units
  if (n < 0) u_->invert();
  return(*this);
  }

// sqrt

UnitVar& UnitVar::sqrt()
  {
  if (no_unit_support_)
    {
    a_ =  std::sqrt(a_);
    return(*this);
    }

  //------------------------------------------------------------------------
  // Need to force base units (even in automatic mode) so that Units::reduce
  // will not have to figure out combinations like "km m" which would
  // take a lot of work.
  //------------------------------------------------------------------------

  coerceBase();

  try
    {
    u_->reduce(2);
    }
  catch(Unit::UnitError& e)
    {
    if (name() != "")
      {
	Unit::UnitError e2("[sqrt(" + name() + ")]" + e.msg, e.error_type);
	e2.throwMe();
      }
    else
      {
	e.throwMe();
      }
    }
  a_ = std::sqrt(a_);
  return(*this);
  }

// exp

UnitVar& UnitVar::exp()
  {
  if (no_unit_support_)
    {
    a_ =  std::exp(a_);
    return(*this);
    }

  if (hasUnits())
    {
      Unit::UnitError e("exp method needs unitless input in " +
      name() + ".exp()", Unit::mismatch);
      e.throwMe();
    }
  a_ = std::exp(a_);
  return(*this);
  }

// log

UnitVar& UnitVar::log()
  {
  if (no_unit_support_)
    {
    a_ =  std::log(a_);
    return(*this);
    }

  if (hasUnits())
    {
      Unit::UnitError e("log method needs unitless input in " +
		      name() + ".log()", Unit::mismatch);
      e.throwMe();
    }
  a_ = std::log(a_);
  return(*this);
  }

// fabs
// Should not throw exceptions
UnitVar& UnitVar::fabs()
  {
  a_ = std::fabs(a_);
  return(*this);
  }

// negate
// Should not throw exceptions
UnitVar& UnitVar::negate()
  {
  a_ = -a_;
  return(*this);
  }

// sin

double UnitVar::sin() const
  {
  if (no_unit_support_)
    {
    return(std::sin(a_));
    }

  if (!hasUnits())
    {
    if (a_==0) return(std::sin(0.0));
    else Unit::UnitError("sin(" + name() + ") needs units (rad)",
				 Unit::mismatch).throwMe();
    }
  return(std::sin(getInUnits("rad")));  // call number sin()
  }

// cos

double UnitVar::cos() const
  {
  if (no_unit_support_)
    {
    return(std::cos(a_));
    }

  if (!hasUnits())
    {
    if (a_==0) return(std::cos(0.0));
    else Unit::UnitError("cos(" + name() + ") needs units (rad)",
				 Unit::mismatch).throwMe();
    }
  return(std::cos(getInUnits("rad")));  // call number cos()
  }

// tan

double UnitVar::tan() const
  {
  if (no_unit_support_)
    {
    return(std::tan(a_));
    }

  if (!hasUnits())
    {
    if (a_==0) return(std::tan(0.0));
    else Unit::UnitError("tan(" + name() + ") needs units (rad)",
				 Unit::mismatch).throwMe();
    }
  return(std::tan(getInUnits("rad")));  // call number tan()
  }

// asin

UnitVar& UnitVar::asin()
  {
  if (no_unit_support_)
    {
    a_ = std::asin(a_);
    return(*this);
    }

  if (hasUnits())
    {
    Unit::UnitError e("asin(" + name() + ") needs unitless input",
      Unit::mismatch);
    e.throwMe();
    }
  a_ = std::asin(a_);
  coerce("rad");
  return(*this);
  }

// acos

UnitVar& UnitVar::acos()
  {
  if (no_unit_support_)
    {
    a_ = std::acos(a_);
    return(*this);
    }

  if (hasUnits())
    {
      Unit::UnitError e("acos(" + name() + ") needs unitless input",
		      Unit::mismatch);
      e.throwMe();
    }
  a_ = std::acos(a_);
  coerce("rad");
  return(*this);
  }

// atan

UnitVar& UnitVar::atan()
  {
  if (no_unit_support_)
    {
    a_ = std::atan(a_);
    return(*this);
    }

  if (hasUnits())
    {
      Unit::UnitError e("atan(" + name() + ") needs unitless input",
		      Unit::mismatch);
      e.throwMe();
    }
  a_ = std::atan(a_);
  coerce("rad");
  return(*this);
  }

// atan2

UnitVar& UnitVar::atan2(const UnitVar& a)
  {
  if (no_unit_support_)
    {
    a_ = std::atan2(a_,a.a_);
    return(*this);
    }

  if (hasUnits() && a.hasUnits())
    {
    // why doesn't this work?? coerce(a.u_);  // match units, throw error if mismatch
    Unit::UnitError e("atan2 needs unitless inputs", Unit::mismatch);
    e.throwMe();
    }
  else if (hasUnits())
    {  // both need to have units if one does
    Unit::UnitError e("atan2 needs matching inputs", Unit::mismatch);
    e.throwMe();
    }
  else if (a.hasUnits())
    {  // both need to have units if one does
    Unit::UnitError e("atan2 needs matching inputs", Unit::mismatch);
    e.throwMe();
    }
  a_ = std::atan2(a_,a.a_);
  coerce("rad");
  return(*this);
  }

//---------------
// Other methods
//---------------

//----------------------------------------------------------------------------
// coerce(ustr)
// coerce(uv)
// coerce(u) - private method (u is const Units&)
//
// Automatic mode:
// If the UnitVar has units, then they are forced to convert to the units
// specified by ustr (or by the units of uv).
// An exception is thrown if conversion is impossible because of a mismatch.
// If the UnitVar has no units, then it is given the specified units.
// The value of the UnitVar (double) is returned after any needed
// conversion has been applied.
// Enforce_base_units mode:
// Check to see if *this is in base units.  If not, throw an error exception.
// Also check the argument for base units.  If not, throw an error exception.
// Otherwise, return the value.
// No_unit_support mode:
// Just return the value.  Assume that the units are in base units.
//----------------------------------------------------------------------------

double UnitVar::coerce(const string& ustr)
  {
  if (no_unit_support_) return(a_);

  Units u(ustr);
  if (enforce_base_units_)
    {
    if (u.hasUnits() && u != u.baseUnits())
      {
      Unit::UnitError e("Need base units in: [" + name() +
        ".coerce(" + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    else if (u_->hasUnits() && *u_ != u_->baseUnits())
      {
      Unit::UnitError e("Target needs to be in base units in: [" + name() +
        ".coerce(" + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    *u_ = u;  // match result's units with those specifed by ustr
    }
  else
    {  // automatic mode - so do the conversion
    if (u_->hasUnits())
      {
      // Match units specified by ustr to *this's units (if possible).
      double c;
      try
        {
        c = u.matchConversion(*u_);
        }
      catch(Unit::UnitError& e)
        {
        Unit::UnitError e2("[" + name() + ".coerce(" + ustr + ")]" + e.msg,
          e.error_type);
        e2.throwMe();
        }
      a_ *= c;
      }
    *u_ = u;  // match result's units with those specifed by ustr
    }
  return(a_);
  }

double UnitVar::coerce(const UnitVar& uv)
  {
  if (no_unit_support_) return(a_);

  if (enforce_base_units_)
    {
    if (uv.u_->hasUnits() && (*(uv.u_) != uv.u_->baseUnits()))
      {
      string ustr;
      uv.u_->toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".coerce(" + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    else if ((uv.u_->hasUnits() && !u_->hasUnits()) ||
             (!uv.u_->hasUnits() && u_->hasUnits()))
      {
      string ustr;
      uv.u_->toString(ustr);
      Unit::UnitError e("No units is not a base unit: [" + name() +
        ".coerce(" + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }
  else if (u_->hasUnits())
    {
    // Match units specified by ustr to *this's units (if possible).
    double c;
    try
      {
      c = uv.u_->matchConversion(*u_);
      }
    catch(Unit::UnitError& e)
      {
      string ustr;
      uv.u_->toString(ustr);
      Unit::UnitError e2("[" + name() + ".coerce(UnitVar(" +
        ustr + "))]" + e.msg, e.error_type);
      e2.throwMe();
      }
    a_ *= c;
    }
  *u_ = *(uv.u_);  // match result's units with those specifed by uv
  return(a_);
  }

double UnitVar::coerce(const Units& u)
  {
  if (enforce_base_units_)
    {
    if (u.hasUnits() && u != u.baseUnits())
      {
      string ustr;
      u.toString(ustr);
      Unit::UnitError e("Need base units in: [" + name() +
        ".coerce<private>(" + ustr + ")]", Unit::mismatch);
      e.throwMe();
      }
    }

  if (u_->hasUnits())
    {
    // Match units specified by ustr to *this's units (if possible).
    double c;
    try
      {
      c = u_->matchConversion(u);
      }
    catch(Unit::UnitError& e)
      {
      string ustr;
      u.toString(ustr);
      Unit::UnitError e2("[" + name() + ".coerce<private>(" + ustr +
        ")]" + e.msg, e.error_type);
      e2.throwMe();
      }
    a_ /= c;
    }
  *u_ = u;  // match result's units with those specifed by ustr
  return(a_);
  }

//-------------------------------------------------------------
// d = coerceBase()
//
// Force the units of *this to be in base units.
// Returns the numerical value after forcing base units.
//-------------------------------------------------------------

double UnitVar::coerceBase()
  {
  
    // if(no_unit_support_) return(a_);

  Units u2 = u_->baseUnits();
  double c;
  try
    {
    c = u_->matchConversion(u2);
    }
  catch(Unit::UnitError& e)
    {
    Unit::UnitError e2(*name_ + ".coerceBase():" + e.msg, e.error_type);
    e2.throwMe();
    }
  a_ /= c;
  *u_ = u2;
  return(a_);
  }

//----------------------------------------------------------------------------
// getInUnits(ustr)
// getInUnits(uv)
//
// The value of this UnitVar is returned after any needed
// conversion has been applied to express it in the units specified by ustr
// (or the units of uv).
// An exception is thrown if conversion is impossible because of a mismatch.
//----------------------------------------------------------------------------

double UnitVar::getInUnits(const string& ustr) const
  {
  // special case for a_==0;
  if(a_==0)
    {
    if(no_unit_support_) return(0);
    if(!hasUnits()) return(0);
    }
  Units u(ustr);
  // Match units specified by ustr to *this's units (if possible).
  double c;
  try
    {
    c = u.matchConversion(*u_);
    }
  catch(Unit::UnitError& e)
    {
    Unit::UnitError e2("[" + name() + ".getInUnits(" + ustr + ")]" + e.msg,
      e.error_type);
    e2.throwMe();
    }
  return(c*a_);
  }

double UnitVar::getInUnits(const UnitVar& uv) const
  {
  // Match units specified by uv to *this's units (if possible).
  double c;
  try
    {
    c = uv.u_->matchConversion(*u_);
    }
  catch(Unit::UnitError& e)
    {
    string ustr;
    uv.u_->toString(ustr);
    Unit::UnitError e2("[" + name() + ".getInUnits(UnitVar(" + ustr +
      "))]" + e.msg, e.error_type);
    e2.throwMe();
    }
  return(c*a_);
  }

double UnitVar::km() const
{
  if(no_unit_support_) return(a_);
  else if(enforce_base_units_){
    if(a_!=0 && *u_!=Units(km_str)){
      ErrorMessage e("UnitVar::km Units are not km.");
      e.throwMe();
    }
  }
  return(getInUnits(km_str));
}

void UnitVar::setkm(double d){
  a_=d;
  if(!no_unit_support_){
    *u_=Units("km");
  }
}
//----------------------------------------------------------------------------
// getValue()
//
// This method returns the value regardless of any units present.
// It is intended for situations where the value is desired even though
// the unit category is not known and coerce or getInUnits can't be used.
// When possible, it is better to use coerce() or getInUnits().
// This inspector method cannot fail.
//----------------------------------------------------------------------------

// should not throw exceptions
double UnitVar::getValue() const 
  {
  return(a_);
  }

//----------------------------------------------------------------------------
// setValue()
//
// This method sets the value regardless of any units present.
// It is intended for situations where the value needs to be set even though
// the units are not known and coerce can't be used.
// When possible, it is better to use coerce().
//----------------------------------------------------------------------------

// should not throw exceptions
void UnitVar::setValue(const double& d)
  {
  a_ = d;
  }

//----------------------------------------------------------------------------
// getUnits()
//
// This method returns a string representation of the units for this UnitVar.
// If no units are present, an empty string is returned.
// This inspector method cannot fail.
//----------------------------------------------------------------------------

string UnitVar::getUnits() const 
  {
  if (no_unit_support_ || enforce_base_units_)
    {
    Unit::UnitError e("UnitVar::getUnits: Can't do this without unit support");
    e.throwMe();
    }
  string ustr;
  u_->toString(ustr);
  return(ustr);
  }

//----------------------------------------------------------------------------
// name()
//
// This method returns the name in a string.
// This inspector method cannot fail.
//----------------------------------------------------------------------------

//should not throw any exceptions
string UnitVar::name() const 
  {
  if (name_) return(*name_);
  else return("");
  }

//------------------------------------------------------------------------
// floor(a)
//
// Sets this UnitVar to the value of a if *this is less than a.
// a can be a UnitVar, in which case normal unit conversions apply,
// or it can be the value 0 passed as a double.  The value zero will then
// acquire the units of *this.
// A non-zero double will cause an exception if *this has units.
//------------------------------------------------------------------------

UnitVar& UnitVar::floor(const double& a)
  {
  if (*this < a)
    {
    a_ = a;
    }
  return(*this);
  }

UnitVar& UnitVar::floor(const UnitVar& a)
  {
  if (*this < a)
    {
    *this = a;
    }
  return(*this);
  }

//------------------------------------------------------------------------
// ceiling(a)
//
// Sets this UnitVar to the value of a if *this is greater than a.
// a can be a UnitVar, in which case normal unit conversions apply,
// or it can be the value 0 passed as a double.  The value 0 will then
// acquire the units of *this.
// A non-zero double will cause an exception if *this has units.
//------------------------------------------------------------------------

UnitVar& UnitVar::ceiling(const double& a)
  {
  if (*this > a)
    {
    a_ = a;
    }
  return(*this);
  }

UnitVar& UnitVar::ceiling(const UnitVar& a)
  {
  if (*this > a)
    {
    *this = a;
    }
  return(*this);
  }

//--------------------//
// Supporting Classes //
//--------------------//

// Check for equal Unit's based soley on their categories.
class CategoryEqual : public std::unary_function<Unit,bool>
  {
  public:
  explicit CategoryEqual(const Unit& cc) : c(cc) { }
  bool operator()(const Unit& u) const { return(u.categoryMatch(c)); }
  protected:
  Unit c;
  };  

//------------------------------
// Supporting functions
//------------------------------

// Binary operators for UnitVar.  Those that need access to the
// representation are friends of UnitVar.

// +
UnitVar operator+(UnitVar a, const UnitVar& b)
  {
  return(a += b);
  }

// -
UnitVar operator-(UnitVar a, const UnitVar& b)
  {
  return(a -= b);
  }

// *
UnitVar operator*(UnitVar a, const UnitVar& b)
  {
  return(a *= b);
  }

// /
UnitVar operator/(UnitVar a, const UnitVar& b)
  {
  return(a /= b);
  }

//-------------------------------------------------------------------
// >>
// In enforce_base_units mode, non base units are converted to base
// units right after reading.
//-------------------------------------------------------------------

std::istream& operator>>(std::istream& s, UnitVar& uvar)
  {
  double a;
  if (!(s >> a)) return s;
  Units u;
  // Units don't have to be present, hence no return here.
  if (!(s >> u)) s.clear();
  uvar.a_ = a;
  *(uvar.u_) = u;
  if (!UnitVar::automatic_)
    {
    uvar.coerceBase();
    }
  return s;
  }

// <<
std::ostream& operator<<(std::ostream& s, const UnitVar& uvar)
  {
  s << uvar.a_ << " " << *(uvar.u_);
  return s;
  }

//------------------------------------
// Other overloaded functions
//------------------------------------

// pow
UnitVar pow(UnitVar a, const double& b)
  {
  return(a.pow(b));
  }

// pow
double pow(const double& a, const UnitVar& b)
  {
  if (b.hasUnits())
    {
    Unit::UnitError e("Exponent can't have units in pow(base,uvar)",
      Unit::mismatch);
    e.throwMe();
    }
  return(pow(a,b.getValue()));
  }

// sqrt
UnitVar sqrt(UnitVar a)
  {
  return(a.sqrt());
  }

// exp
double exp(UnitVar a)
  {
  return(a.exp().getValue());
  }

// log
double log(UnitVar a)
  {
  return(a.log().getValue());
  }

// fabs
UnitVar fabs(UnitVar a)
  {
  return(a.fabs());
  }

// sin
double sin(UnitVar a)
  {
  return(a.sin());
  }

// cos
double cos(UnitVar a)
  {
  return(a.cos());
  }

// tan
double tan(UnitVar a)
  {
  return(a.tan());
  }

// asin
UnitVar asin(UnitVar a)
  {
  return(a.asin());
  }

// acos
UnitVar acos(UnitVar a)
  {
  return(a.acos());
  }

// atan
UnitVar atan(UnitVar a)
  {
  return(a.atan());
  }

// atan2
UnitVar atan2(UnitVar a, UnitVar b)
  {
  return(a.atan2(b));
  }

// floor
UnitVar floor(UnitVar a, const double& b)
  {
  return(a.floor(b));
  }

// floor
UnitVar floor(UnitVar a, const UnitVar& b)
  {
  return(a.floor(b));
  }

// ceiling
UnitVar ceiling(UnitVar a, const double& b)
  {
  return(a.ceiling(b));
  }

// ceiling
UnitVar ceiling(UnitVar a, const UnitVar& b)
  {
  return(a.ceiling(b));
  }

//---------------------------------------------------------------------------
// convert_to_double(d), convert_to_double(u)
//
// Overloaded function to convert UnitVar (or double) to a double
// This is useful for templates in which it is unknown whether a variable
// is a double or a UnitVar, so UnitVar::getValue() cannot be used
// convert_to_double
//---------------------------------------------------------------------------

double convert_to_double(double d)
  {
  return(d);
  }

double convert_to_double(const UnitVar& u)
  {
  if (!UnitVar::no_unit_support_ && u.hasNonSuppressedUnits())
    {
    Unit::UnitError e("convert_to_double requires unitless input",
		      Unit::mismatch);
    e.throwMe();
    }
  return(u.getValue());
  }

int convert_to_int(const UnitVar& u)
  {
  // hasUnits instead of hasNonSuppressedUnits (see above) because it seems
  // unlikely that integers and "rad" will be mixed together!
  if (!UnitVar::no_unit_support_ && u.hasUnits())
    {
    Unit::UnitError e("convert_to_int requires unitless input",
		      Unit::mismatch);
    e.throwMe();
    }
  double d=u.getValue();
  int i=(int) d;
  if (i!=d)
    {
    Unit::UnitError e("convert_to_int found non-integer value",
		      Unit::unspecified);
    e.throwMe();
    }
  return(i);
  }

//---------------------------------------------------------------------
// d = coerce_base(d,ustr)
//
// Build a UnitVar with the input value and units, converting to base
// units as needed and then return the value of the converted result.
//---------------------------------------------------------------------

Uvar coerce_base(const double& d, const string& ustr)
  {
  Units uu(ustr);
  Units u2 = uu.baseUnits();
  double c;
  try
    {
    c = uu.matchConversion(u2);
    }
  catch(Unit::UnitError& e)
    {
    Unit::UnitError e2("coerce_base():" + e.msg, e.error_type);
    e2.throwMe();
    }
  Uvar u(d/c,u2);
  return(u);
  }

//---------------------------------------------------------------------
// uv = coerce(d,ustr)
// uv = coerce(uv1,ustr)
//
// Build a UnitVar with the input value and units.
// Unlike coerce_base above, this function does not convert to base
// units.  If base unit mode is set, and ustr is not a base unit, then
// an error exception will be thrown.  In no_unit_support mode, the ustr
// argument is ignored.  When the 1st argument is also a UnitVar, then
// it is converted to the units in ustr (is possible).  If the input UnitVar
// has no units, then the ustr units are attached.
// Note that this function returns a UnitVar.  If you want the converted
// double result, use get_in_units() below instead.
//---------------------------------------------------------------------

Uvar coerce(const double& d, const string& ustr)
  {
  return(UnitVar(d,ustr));
  }

Uvar coerce(Uvar uv, const string& ustr)
  {
  uv.coerce(ustr);
  return(uv);
  }

//---------------------------------------------------------------------
// d = get_in_base_units(u)
// d = get_in_base_units(d)
//
// Convert the input UnitVar to base units and return the value.
// If input is a double, just return the value.
//---------------------------------------------------------------------

double get_in_base_units(UnitVar uv)
  {
    // Right now this only works for no unit support need to fix
  return(uv.getValue());
  }

double get_in_base_units(const double& d)
  {
  return(d);
  }

//---------------------------------------------------------------------
// d = get_in_units(u,"ustr")
//
// Convert the UnitVar u to the units in ustr and return the value.
//---------------------------------------------------------------------

double get_in_units(const UnitVar& uv, const string& ustr)
  {
  return(uv.getInUnits(ustr));
  }

//-------------------
// Class Unit Methods
//-------------------

//-----------------
// Construction
//-----------------

Unit::Unit() { }
Unit::Unit(UTI utp) : utp_(utp) { }

// Predicates

//-------------------------------------------------------------------------
// categoryMatch
//
// Category matches lie in the same unit category (eg., length) but are not
// necessarily exact matches (eg., m and km).
//-------------------------------------------------------------------------

// should not throws any exceptions
bool Unit::categoryMatch(const Unit& u) const 
  {
  return(utp_->second.category_index == u.utp_->second.category_index);
  }

//---------------------------------
// == and != require exact matches
//---------------------------------

// should not throws any exceptions
bool Unit::operator==(const Unit& u) const 
  {
  return(utp_ == u.utp_);
  }

// should not throws any exceptions
bool Unit::operator!=(const Unit& u) const 
  {
  return(utp_ != u.utp_);
  }

//--------------------------------------------------------------------
// isUnit(ustr)
// Static member function returns true if ustr is a valid string
// representation of a single unit term.
//--------------------------------------------------------------------

// should not throw any exceptions
bool Unit::isUnit(const string& ustr) 
  {
  if (ustr == "1") return(true);  // 1 is a special case

  Unit::UTI utp = Unit::unit_table().find(ustr);
  if (utp != Unit::unit_table().end())
    {  // found a matching unit name in the unit table
    return(true);
    }

  //--------------------------------------------------------------------
  // Didn't find the term in the unit table, so check the derived unit
  // table.
  //--------------------------------------------------------------------

  Unit::UDI udp = Unit::unit_derived().find(ustr);
  if (udp != Unit::unit_derived().end())
    {  // found a matching unit name in the derived unit table
    return(true);
    }

  //--------------------------------------------------------------------
  // Didn't find the term in the derived unit table either, so it isn't
  // a known unit term.
  //--------------------------------------------------------------------

  return(false);
  }

//--------------------
// Class Units Methods
//--------------------

// Construction/Destruction

Units::Units() { }

//-------------------------------------------------------------------------
// Construct from string
// This constructor converts a string into a Units object.
// An exception is thrown if the unit string is unrecognized, or formatted
// incorrectly.
//-------------------------------------------------------------------------

Units::Units(const string& ustr) 
  {
  string ustr1 = ustr + " ";
  ISTRINGSTREAM ist(ustr1.c_str());
  streamUnitsIn(ist);
  }

//----------------------------------------------------------------------------
// Construct from file
// This constructor reads and converts a character string from the indicated
// file.  By default, the number of characters is stored first as an
// unsigned char, followed by the characters themselves.
// Like above, an exception is thrown for bad unit strings.
// See read().
//----------------------------------------------------------------------------

Units::Units(const FileMgr& file, int usize)
  {
  read(file,usize);
  }

//---------
// Testing
//---------

//----------------------------------------------------------
// selfTest
//
// Exercize the UnitVar type and verify proper operation.
//----------------------------------------------------------

bool Units::selfTest(const string& filename)
  {
//  Unit::specifyUnitConversions(filename);
  UnitVar::setMode("automatic");

  UnitVar cc("speed_light");
  string s = "2.997924580e5 km/s";  // check reading of exponential notation
  ISTRINGSTREAM ist1(s.c_str());
  ist1 >> cc;
  if (cc != UnitVar(2.997924580e5,"km/s")) return(false);
  // check builtin constant
  if (cc != UnitVar("speed_light")) return(false);

  s = "1.3806503e-23 J/K";  // check reading of exponential notation
  ISTRINGSTREAM ist2(s.c_str());
  UnitVar bc;
  ist2 >> bc;
  if (bc != UnitVar(1.3806503e-23,"J/K")) return(false);
  // check builtin constant
  if (bc != UnitVar("boltzmann_constant")) return(false);

  UnitVar frequency("frequency");
  frequency = UnitVar(14,"GHz");
  if (!frequency.hasUnits()) return(false);
  if (std::fabs(frequency.getInUnits("Hz") - 14e9) > 1e-5) return(false);
  UnitVar vel("vel",4,"m/s");
  UnitVar mass("mass",3,"kg");
  UnitVar ke("ke");
  ke = mass*pow(vel,2);
  UnitVar kematch(3,"J");
  ke += UnitVar(3,"J");  // check for proper matching of derived units
  UnitVar a;
  if (a.coerce("") != 0) return(false);
  if (a.hasUnits()) return(false);
  UnitVar b(3.5);
  if (b.coerce("km") != 3.5) return(false);
  if (!b.hasUnits()) return(false);
  // Can't mix different UnitVar types yet.
  // Need to figure out how to generalize the template definitions, or supply
  // specializations.
  //  UnitVar<float> c(7.25,"m");
  UnitVar c(7.25,"km");
  if (fabs(c.coerce("cm") - 725000) > 1e-8) return(false);
  if (c.getInUnits("mm") != 7250000) return(false);
  if (!(b+c == UnitVar(10.75,"km"))) return(false);
  UnitVar t("t");
  t.stringIn("240 s");
  if (t.coerce("min") != 4) return(false);
  a = UnitVar(0.5);
  b = asin(a);
  if (fabs(b.getInUnits("deg") - 30.0) > 1e-8) return(false);
  bool errflag = false;
  try
    {
    a = c + t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  // Operator tests
  float g = 5.2;
  UnitVar g1 = g*t;
  UnitVar g2 = t*g;
  g1 *= g;
  g1 *= t;
  g1 = g/t;
  g2 = t/g;
  g1 /= g;
  g1 /= t;

  try
    {
    g1 = g+t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  try
    {
    g1 = g-t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  UnitVar theta(30,"deg");
  double d = sin(theta);
  if (std::fabs(d - 0.5) > 1e-10) return(false);
  UnitVar theta2 = asin(UnitVar(d));
  if (fabs(theta2 - theta).getValue() > 1e-10) return(false);

  theta = UnitVar(60,"deg");
  d = cos(theta);
  if (std::fabs(d - 0.5) > 1e-10) return(false);
  theta2 = acos(UnitVar(d));
  if (fabs(theta2 - theta).getValue() > 1e-10) return(false);

  d = tan(theta);
  if (std::fabs(d - sqrt(3.0)) > 1e-10) return(false);
  theta2 = atan(UnitVar(d));
  if (fabs(theta2 - theta).getValue() > 1e-10) return(false);

  double aa = exp(UnitVar(3.5));
  double bb = log(aa);
  if (fabs(bb - 3.5) > 1e-10) return(false);

  string ustr = "4 m m/(s s)";
  ISTRINGSTREAM ist3(ustr.c_str());
  ist3 >> a;
  if (a != UnitVar(4,"m m/(s s)")) return(false);

  b = sqrt(a);
  c = UnitVar(2,"m / s");
  if (b != c) return(false);
  a = UnitVar(9,"km m / s min");
  b = sqrt(a);
  if (fabs(b.getInUnits("m/min") - sqrt(9*1000*60)) > 1e-6) return(false);

  a = UnitVar(3.5,"km"); 
  b = UnitVar(3600,"m");
  if (a > b) return(false);
  if (a >= b) return(false);
  if (a == b) return(false);
  if (!(a != b)) return(false);

  // check suppressed unit
  a = UnitVar(2,"rad");
  if (a > 2.001) return(false);

  // Now check the enforce_base_units mode
  UnitVar::setMode("enforce_base_units");
  UnitVar x1(8,"km");
  UnitVar t1(16,"s");
  UnitVar v1 = x1/t1;
  UnitVar v2(0.5,"km/s");
  if (v1 != v2) return(false);
  double d1 = x1.getInUnits("m");
  if (d1 != 8000) return(false);
  s = "15000000 m m/(s min)";
  ISTRINGSTREAM ist4(s.c_str());
  ist4 >> a;
  if (a.getValue() != 0.25) return(false);

  s = "1.3806503e-23 J/K";  // check reading of exponential notation
  ISTRINGSTREAM ist5(s.c_str());
  ist5 >> bc;
  if (fabs(bc - UnitVar(1.3806503e-29,"kg km km/(s s K)"))/bc > 1e-10)
    {
    return(false);
    }
  // check builtin constant
  if (fabs(bc - UnitVar("boltzmann_constant"))/bc > 1e-10) return(false);
  // check external constant
  if (fabs(bc - boltzmann_constant)/bc > 1e-10) return(false);

  // Now check the no_unit_support mode
  UnitVar::setMode("no_unit_support");
  UnitVar t2(90,"s");
  double d2 = t2.getInUnits("min");
  if (d2 != 1.5) return(false);

  UnitVar::setMode("automatic");
  
  return(true); 
  }

// Predicates

//--------------------------------------------------------------
// hasUnits
// Returns true if there are any numerator or denominator units.
// Should not throw any exceptions.
//--------------------------------------------------------------

bool Units::hasUnits() const 
  {
  return(units_top_.size() != 0 || units_bot_.size() != 0);
  }

//----------------------------------------------------------------------------
// hasNonSuppressedUnits()
//
// Returns true if this has any unit terms besides the suppressed unit.
// Return false if units are empty or only have instances of the suppressed
// unit.
//----------------------------------------------------------------------------

bool Units::hasNonSuppressedUnits() const
  {
  for (list<Unit>::const_iterator p = units_top_.begin();
       p != units_top_.end(); ++p)
    {
    if (*p != Unit::suppressed_unit()) return(true);
    }
  for (list<Unit>::const_iterator p = units_bot_.begin();
       p != units_bot_.end(); ++p)
    {
    if (*p != Unit::suppressed_unit()) return(true);
    }
  return(false);
  }

bool Units::operator!=(const Units& u) const 
  {
  return(units_top_ != u.units_top_ || units_bot_ != u.units_bot_);
  }

// I/O

//-------------------------------------------------------------------------
// write(file,usize)
//
// Convert this Units object to a unit string, and write the string to the
// indicated file.
//
// file = FileMgr object to write to
// usize = the number of characters to write for the unit string.
//   If usize is 0, then no units are written.
//   If usize is -1, then the actual size of the unit string is written.
//      The length of the unit string is written as a unsigned char, followed
//      by the string itself (no trailing 0).  If the unit string exceeds 255
//      chars, then an exception is thrown.
//      (this is the default value if usize is not explicitely passed)
//   If usize is positive, then precisely usize characters are written
//      after the size field, and no trailing 0.  If the unit string is
//      larger than usize chars, then an exception is thrown.
//      Excess characters are padded with spaces.
//-------------------------------------------------------------------------

void Units::write(const FileMgr& file, int usize) const
  {
  if (usize == 0) return;

  string ustr;
  toString(ustr);
  if (ustr.length() > 255)
    {
      Unit::UnitError e("More than 255 characters in units to write to file " +
      file.name());
      e.throwMe();
    }
  unsigned char len = ustr.length();

  if (usize == -1)
    {
    file.write(len);
    file.write(ustr);
    }
  else if (len > usize)
    {
    Unit::UnitError e("More than " + toStr(usize) +
      " characters in units to write to file " + file.name());
    e.throwMe();
    }
  else if (len == usize)
    {
    file.write(len);
    file.write(ustr);
    }
  else
    {
    string s(usize-len,' ');
    file.write(len);
    file.write(ustr);
    file.write(s);
    }

  }

//-------------------------------------------------------------------------
// read(file,usize)
//
// Read a unit string from the file.  usize determines how the string is
// stored just like write() above.
//-------------------------------------------------------------------------

void Units::read(const FileMgr& file, int usize)
  {
  if (usize == 0) return;

  if (usize == -1)
    {
    unsigned char len;
    file.read(len);
    string ustr(' ',len);
    file.read(ustr);
    ISTRINGSTREAM ist(ustr.c_str());
    streamUnitsIn(ist);
    }
  else if (usize > 0)
    {
    string ustr(' ',usize);
    file.read(ustr);
    ISTRINGSTREAM ist(ustr.c_str());
    streamUnitsIn(ist);
    }
  else
    {
      Unit::UnitError e("Bad usize parameter passed to Units::read");
      e.throwMe();
    }

  }

void Units::show() const
  {
  cout << *this << endl;
  }

//------------------------
// Unit conversion support
//------------------------

//-------------------------------------------------------------------------
// matchConversion(u2)
//
// In automatic mode:
// Compute the unit conversion factor which needs to multiply the argument
// units (u2) to make them match *this's units (u_).  If all unit factors
// can't be matched, then a UnitMismatch exception is thrown.
// In enforce_base_units_ mode:
// Check to see if u2 is in base units.  If not, throw an exception.
// Also compute and return the conversion as described above.
// In no_unit_support mode:
// Assume that u2 matches *this in base units, and compute the conversion
// to *this.
//-------------------------------------------------------------------------

double Units::matchConversion(const Units& u2) const 
  {
  Units u2use;
  if (UnitVar::no_unit_support_)
    {
    // Assume *this matches argument in base units, but don't use u2.
    u2use = baseUnits();
    }
  else if (UnitVar::enforce_base_units_)
    {
    // Make sure u2 is in base units.
    if (u2 != u2.baseUnits())
      {
      Unit::UnitError e2(Unit::mismatch);
      e2.throwMe();
      }
    u2use = u2;
    }
  else
    {
    u2use = u2;
    }

  // For each u_ top/bot unit, find a category match in u2 and accumulate
  // the unit conversions (top/bot separately).
    
  double conv_factor_top = 1.0;

  list<Unit> unit_match_list = u2use.units_top_;  // make copy
  list<Unit>::iterator u2p;
  for (list<Unit>::const_iterator u1p=units_top_.begin();
       u1p != units_top_.end(); ++u1p)
    {
    if (factorMatch(*u1p,unit_match_list,conv_factor_top,u2p))
      {
      unit_match_list.erase(u2p);
      }
    else if (*u1p != Unit::suppressed_unit())
      {
      // No match found, so this is a unit mismatch because terms require
      // a category match for every unit factor.
      Unit::UnitError e(Unit::mismatch);
      e.throwMe();
      }
    }

  if (unit_match_list.size() != 0)
    {
    removeSuppressedUnit(unit_match_list);
    // If there are any unmatched units remaining in the match list, then we
    // have a mismatch because every unit in u2 should be matched in *this.
    if (unit_match_list.size() != 0)
      {
      Unit::UnitError e(Unit::mismatch);
      e.throwMe();
      }
    }

  // Again for the denominators

  unit_match_list = u2use.units_bot_;  // make copy
  double conv_factor_bot = 1.0;
  for (list<Unit>::const_iterator u1p=units_bot_.begin();
       u1p != units_bot_.end(); ++u1p)
    {
    if (factorMatch(*u1p,unit_match_list,conv_factor_bot,u2p))
      {
      unit_match_list.erase(u2p);
      }
    else if (*u1p != Unit::suppressed_unit())
      {
      // No match found, so this is a unit mismatch because terms require
      // a category match for every unit factor.
      Unit::UnitError e(Unit::mismatch);
      e.throwMe();
      }
    }

  if (unit_match_list.size() != 0)
    {
    removeSuppressedUnit(unit_match_list);
    if (unit_match_list.size() != 0)
      {
      Unit::UnitError e(Unit::mismatch);
      e.throwMe();
      }
    }
  
  return(conv_factor_top/conv_factor_bot);
  }

//--------------------------------------------------------------------------
// factorMerge
//
// Merge the unit terms of u2 into *this's units.
// The terms are multiplied together if mtype is product, and divided
// together if mtype is quotient.
// The return value is the accumulated conversion factor that was applied
// to cancel unit terms from the same categories.
// Should not throw any exceptions.
//--------------------------------------------------------------------------

double Units::factorMerge(const Units& u2, mergeE mtype)
  {
  double conv_factor_top = 1.0;
  double conv_factor_bot = 1.0;

  // Make copies which are stripped of matches below.
  list<Unit> u2_copy_top;
  list<Unit> u2_copy_bot;

  if (mtype == product)
    {
    u2_copy_top = u2.units_top_;
    u2_copy_bot = u2.units_bot_;
    }
  else
    {
    u2_copy_top = u2.units_bot_;
    u2_copy_bot = u2.units_top_;
    }

  // Process u1 (ie., *this) numerator.

  list<Unit>::iterator u1p = units_top_.begin();
  list<Unit>::iterator u2p;
  while (u1p != units_top_.end())
    {
    if (factorMatch(*u1p,u2_copy_bot,conv_factor_bot,u2p))
      {  // found a match, so cancel these units
      u1p = units_top_.erase(u1p);
      u2p = u2_copy_bot.erase(u2p);
      }
    else
      {
      ++u1p;
      }
    }

  // Process u1 (ie., *this) denominator.

  u1p = units_bot_.begin();
  while (u1p != units_bot_.end())
    {
    if (factorMatch(*u1p,u2_copy_top,conv_factor_top,u2p))
      {  // found a match, so cancel these units
      u1p = units_bot_.erase(u1p);
      u2p = u2_copy_top.erase(u2p);
      }
    else
      {
      ++u1p;
      }
    }
  
  // Merge the remaining units in u2 into the units of u1.
  units_top_.splice(units_top_.end(),u2_copy_top);
  units_bot_.splice(units_bot_.end(),u2_copy_bot);

  return(conv_factor_top/conv_factor_bot);
  }

//-------------------------------------------------------------------------
// bu = baseUnits()
//
// Return a unit object containing the base units for *this.
// The base unit for each category is the first unit name for that category
// in the unit table.
//-------------------------------------------------------------------------

Units Units::baseUnits() const
  {
  Units base_units;
  for (list<Unit>::const_iterator p=units_top_.begin();
       p != units_top_.end(); ++p)
    {
    base_units.units_top_.push_back(p->utp_->second.base_unit);
    }

  for (list<Unit>::const_iterator p=units_bot_.begin();
       p != units_bot_.end(); ++p)
    {
    base_units.units_bot_.push_back(p->utp_->second.base_unit);
    }

  return(base_units);
  }

//---------------
// Other services
//---------------

//-------------------
// clear()
// Removes any units
//-------------------

// Should not throw any exceptions
void Units::clear() 
  {
  units_top_.clear();
  units_bot_.clear();
  }

//------------------------------------------------------
// invert()
// Swaps the numerator units with the denominator units.
//------------------------------------------------------

// should not throw any exceptions
void Units::invert() 
  {
  std::swap(units_top_,units_bot_);
  }

//-------------------------------------------------------------------
// combine()
// Merges the units (num and den) of supplied Units with *this units.
//-------------------------------------------------------------------

// should not throw any exceptions
void Units::combine(Units ucopy) 
  {
  units_top_.splice(units_top_.end(),ucopy.units_top_);
  units_bot_.splice(units_bot_.end(),ucopy.units_bot_);
  }

//---------------------------------------------------------
// toString(ustr)
// s = toString()
// Convert this Units object into a string representation.
// (the opposite of toUnits)
//---------------------------------------------------------

void Units::toString(string& ustr) const 
  {
  string s;
  if (units_top_.size() == 0 && units_bot_.size() == 0)
    {
    return;
    }
  if (units_top_.size() == 0)
    {
    s += "1";
    }
  else if (units_top_.size() == 1)
    {
    s += (units_top_.front()).utp_->first;
    }
  else
    {
    list<Unit>::const_iterator p = units_top_.begin();
    for (unsigned int i=0; i < units_top_.size()-1; ++i)
      {
      s += p->utp_->first + " ";
      ++p;
      }
    s += p->utp_->first;
    }

  if (units_bot_.size() == 0)
    {
    ustr = s;
    return;
    }
  else if (units_bot_.size() == 1)
    {
    s += "/" +  (units_bot_.front()).utp_->first;
    }
  else
    {
    s += "/(";
    list<Unit>::const_iterator p = units_bot_.begin();
    for (unsigned int i=0; i < units_bot_.size()-1; ++i)
      {
      s += p->utp_->first + " ";
      ++p;
      }
    s += p->utp_->first + ")";
    }

  ustr = s;
  }

string Units::toString() const 
  {
  string ustr;
  toString(ustr);
  return(ustr);
  }

/****
This stream version doesn't work even though it should.  There appears to be
some errors in the strstream library implementation on cygwin/egcs.  Also,
the str() method of a stringstream doesn't return a string like it should.
Instead, it returns a char*.

void Units::toString(string& ustr) const
  {
  OSTRINGSTREAM s;
  s.flush();
  cout << "toString:";
//  cout << s.str() << " ";
  if (units_top_.size() == 0 && units_bot_.size() == 0)
    {
    return;
    }
  if (units_top_.size() == 0)
    {
    s << "1";
    }
  else if (units_top_.size() == 1)
    {
    s << units_top_.front();
    s.flush();
    cout << (units_top_.front()).utp_->first;
//    cout << s.str() << " ";
    }
  else
    {
    list<Unit>::const_iterator p = units_top_.begin();
    for (unsigned int i=0; i < units_top_.size()-1; ++i)
      {
      s << *p << " ";
      ++p;
      }
    s << *p;
    }

  if (units_bot_.size() == 0)
    {
    ustr = toStr(s);
    return;
    }
  else if (units_bot_.size() == 1)
    {
    s << "/" << units_bot_.front();
    s.flush();
    cout << (units_bot_.front()).utp_->first;
//    cout << s.str() << " ";
    }
  else
    {
    s << "/(";
    list<Unit>::const_iterator p = units_bot_.begin();
    for (unsigned int i=0; i < units_bot_.size()-1; ++i)
      {
      s << *p << " ";
      ++p;
      }
    s << *p << ")";
    }

//  cout << s.str() << " ";
  if (toStr(s)) ustr = toStr(s);
  cout << ustr << " " << s.str() << ":toString" << endl;
  cout << ustr.size() << endl;
  exit(1);
  }
*****/

//---------------------------------------------------------
// toUnits()
// Convert the supplied string into *this Units.
// (the opposite of toString)
// Throws an exception if the string is not valid units.
//---------------------------------------------------------

void Units::toUnits(const string& ustr) 
  {
  string ustr1 = ustr + " ";
  ISTRINGSTREAM ist(ustr1.c_str());
  streamUnitsIn(ist);
  }

//-------------------------------------------------
// reduce(n)
//
// Reduce the units in *this by the n'th power.
// Thus, "m m/(s s)" is reduced to "m/s" if n = 2.
// Empty units mean no action.
// n = 0 or 1 means no action.
//-------------------------------------------------

void Units::reduce(unsigned int n) 
  {
  if (!hasUnits()) return;
  reduceList(units_top_, n);
  reduceList(units_bot_, n);
  }

//----------------------------------------------------------------------
// removeSuppressedUnit()
//
// Remove any instances of the suppressed unit from this's unit lists.
//----------------------------------------------------------------------

void Units::removeSuppressedUnit()
  {
  removeSuppressedUnit(units_top_);
  removeSuppressedUnit(units_bot_);
  }

//----------------------------------------------------------------------
// ***** This method is currently not being used.
// tokenize()
// Parse's strings into subwords suitable for unit analysis.
// Declared static in Units.h because it does not need to access a Units
// object.
// **** Needs to handle () just like stream output.
//----------------------------------------------------------------------

void Units::tokenize(const string& ustr, list<string>& tokens)
    {
  ISTRINGSTREAM ist(ustr.c_str());
  string token;
  char ch = 0;
  while (ist >> ch)
    {
    if (ch == '/')
      {
      tokens.push_back("/");
      }
    else if (ch == '1')
      {
      tokens.push_back("1");
      }
    else if (isalpha(ch))
      {
      ist.putback(ch);
      ist >> token;
      tokens.push_back(token);
      }
    else
      {
	Unit::UnitError("Bad unit string: " + ustr).throwMe();
      }
    }
  }

//----------------------------
// Protected internal methods
//----------------------------

//--------------------------------------------------
// reduceList(ulist,n)
//
// Reduce a unit list by the order n (see reduce())
//--------------------------------------------------

void Units::reduceList(list<Unit>& ulist, unsigned int n) 
  {
  //---------------------------------------
  // Reduce successively by dumping matches
  //---------------------------------------

  for (unsigned int i=2; i <= n; ++i)
    {
    for (list<Unit>::iterator p = ulist.begin();
         p != ulist.end(); ++p)
      {  // look for a match to each unit term
      Unit u = *p;  // unit term to match
      ++p;          // advance to start of search list
      list<Unit>::iterator pmatch = find(p,ulist.end(),u);
      --p;          // backup to original position
      if (pmatch == ulist.end())
        {  // no match - this is an error
          Unit::UnitError e("UnitVar reduction error; No match for " +
			    p->utp_->first, Unit::reduction_error);
	  e.throwMe();
        }
      else
        {  // found match, dump the match
        ulist.erase(pmatch);
        }
      }
    }
      
  }

//-----------------------------------
// tokensToUnits
//
// Convert token lists into unit list.
//-----------------------------------

// Need to rework some of this so that stream input is used for all processing
// so that I don't need to put the () processing into tokenize.
// tokensToUnits should accept numerator and denominator strings separately,
// and it should do the recursive call via toUnits.  The derived unit
// mapping should probably be strings instead of token lists.

void Units::tokensToUnits(list<string>& top_tokens, list<string>& bot_tokens,
  list<Unit>& utop, list<Unit>& ubot) 
  {
  // Process numerator tokens, expanding out derived units, and adding
  // the Unit terms to utop and ubot.

  for (list<string>::iterator p=top_tokens.begin();
       p != top_tokens.end(); ++p)
    {
    if (*p == "1") continue;  // skip over numerator '1'
    Unit::UTI utp = Unit::unit_table().find(*p);
    if (utp != Unit::unit_table().end())
      {  // found a matching unit name in the unit table
      utop.push_back(Unit(utp));
      }
    else
      { // No fundamental unit match, so look for a derived unit match.
      Unit::UDI udp = Unit::unit_derived().find(*p);
      if (udp != Unit::unit_derived().end())
        {  // found matching derived unit: recursively process it
        tokensToUnits(udp->second.top_tokens,udp->second.bot_tokens,utop,ubot);
        }
      else
        {  // No matching derived unit either, so this is an unknown unit.
         Unit::UnitError("Bad unit string: " + *p).throwMe();
        }
      }
    }
      
  // Process denominator tokens, expanding out derived units, and adding
  // the Unit terms to utop and ubot in inverted sense (ie., numerator of
  // denominator units go in the denominator, and denominator of denominator
  // units go in the numerator!)

  for (list<string>::iterator p=bot_tokens.begin();
       p != bot_tokens.end(); ++p)
    {
    Unit::UTI utp = Unit::unit_table().find(*p);
    if (utp != Unit::unit_table().end())
      {  // found a matching unit name in the unit table
      ubot.push_back(Unit(utp));
      }
    else
      { // No fundamental unit match, so look for a derived unit match.
      Unit::UDI udp = Unit::unit_derived().find(*p);
      if (udp != Unit::unit_derived().end())
        {  // found matching derived unit: recursively process it
        // Note reversal of utop,ubot order because these are denominator
        // units.
        tokensToUnits(udp->second.top_tokens,udp->second.bot_tokens,ubot,utop);
        }
      else
        {  // No matching derived unit either, so this is an unknown unit.
	  Unit::UnitError("Bad unit string: " + *p).throwMe();
        }
      }
    }

  }

//-------------------------------------------------------
// *****This method not used anymore.
// unitSplit
// Split token list into numerator and denominator parts.
//-------------------------------------------------------

// should not throw any exceptions
void Units::unitSplit(list<string>& tokens,
  list<string>& toptokens, list<string>& bottokens)
  {
  list<string>::iterator pdiv = find(tokens.begin(),tokens.end(),"/");
  if (pdiv == tokens.end())
    {  // no division symbol, so everything is in the numerator
    toptokens = tokens;
    return;
    }
  for (list<string>::iterator p = tokens.begin(); p != pdiv; ++p)
    {
    toptokens.push_back(*p);
    }
  ++pdiv;  // move to start of denominator list
  for (list<string>::iterator p = pdiv; p != tokens.end(); ++p)
    {
    if (*p == "(" || *p == ")") continue;  // skip any enclosing ()
    bottokens.push_back(*p);
    }
  }
  
//-------------------------------------------------------------------
// factorMatch
// Look for a match for u1 in the list unit_match_list.  When found,
// accumulate the conversion, and return an iterator to the matching
// element.  Return true if a match found, false otherwise.
// Should not throw any exceptions
//-------------------------------------------------------------------

bool Units::factorMatch(const Unit& u1,
  list<Unit>& unit_match_list, double& conv_factor,
  list<Unit>::iterator& u2p)
  {
  u2p=find_if(unit_match_list.begin(),unit_match_list.end(),CategoryEqual(u1));

  // If no match found, return without doing anything.
  if (u2p == unit_match_list.end()) return(false);

  // If the units within the category don't match, then a conversion
  // needs to be applied to the 2nd term (ie., *u2p) to match the
  // 1st term (u1).
  if (*u2p != u1)
    {
    conv_factor *= u2p->utp_->second.conv_factor /
                       u1.utp_->second.conv_factor;
    }
  return(true);
  }

//---------------------------------------------------------------------------
// streamUnitTermsIn(s,term_handling,uterms_list)
//
// Stream input of Unit terms into a list of strings.
// Each string in the list is a potential individual unit term.
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
// This function only handles simple lists of terms (eg., m m s s) but not
// numerator/denominator combinations which are handled by the << operator
// (with help from this protected function).
// The '/' and ')' characters are recognized as special, however,
// so that spaces are not required between the last numerator term and the '/',
// or between the last denominator term and the ')'.
// The '1' character is also recognized as special (for cases where there
// are no numerator units), and is a legitimate unit term.
// The last character read is returned because the std strstream appears to
// have a bug that prevents the last character on a stream from being putback.
// Thus, it needs to be checked separately.
// term_handling determines if individual unit terms will be checked against
// the list of known unit terms.
//---------------------------------------------------------------------------

// should not throw any exceptions
char Units::streamUnitTermsIn(istream& s, term_handlingE term_handling,
  list<string>& uterms_list)
  {
  string word;
  string ustr;  // accumulates the valid unit string
  string unit_term;

  //----------------------------------------------------------------------
  // Pull in whitespace delimited words and check each to see if it might
  // contribute to a valid unit string.
  //----------------------------------------------------------------------

  while (1)
    {
    streampos sp = s.tellg();  // remember current stream position
    if (!(s >> word)) break;   // no more input
    // need to clear the flags here or subsequent seekg commands will fail
    s.clear(std::ios::goodbit);

    string::size_type n1 = word.find('/');
    string::size_type n2 = word.find(')');
    if (n1 == string::npos && n2 == string::npos &&
        (term_handling == allow_any_terms || Unit::isUnit(word))
       )
      {  // simple unit term
      uterms_list.push_back(word);
      continue;
      }
    else if (n1 == string::npos && n2 == string::npos)
      {  // invalid unit term; this ends further unit input
      s.seekg(sp);
      break;  // done reading stream
      }
    else
      {  // Found a '/' or ')' so we have a possible unit term; xxx/.. or xxx)
      string::size_type n = min(n1,n2); // select the first occuring
      unit_term = word.substr(0,n);  // extract xxx before the first occuring
      if (term_handling == known_terms_only && !Unit::isUnit(unit_term))
        {  // invalid unit term; this ends further unit input
        if (unit_term.size() == 0)
          {  // the '/' or ')' is by itself, so reposition on it
          for (unsigned int i=0; i < word.length() - n; ++i)
            {  // put each denominator char back in stream buffer
            s.putback(word[word.length()-1-i]);
            }
          }
        else
          {  // reposition at start of invalid unit term
          s.seekg(sp);
          }
        }
      else
        {  // unit term before '/' or ')' is legitimate
        uterms_list.push_back(unit_term);  // final unit term
        // reposition on the '/' or ')'
        for (unsigned int i=0; i < word.length() - n; ++i)
          {  // put each denominator char back in stream buffer
          s.putback(word[word.length()-1-i]);
          }
// Statement below doesn't work because negative offsets are not allowed
// by the 3.4 compiler.
//        s.seekg(n - word.size(), std::ios::cur);
        // We need to clear the flags here just in case the stream ended right
        // after the last unit term which would leave the eof flag set, and
        // cause subsequent processing (of the denominator) to fail in an
        // obscure way!
        s.clear(std::ios::goodbit);
        }
      break;  // done reading stream
      }
    }

  if (uterms_list.empty())
    {  // no valid unit terms and not eof, so set failed
    s.clear(std::ios::failbit);
    }

  // Return last character read from stream
  return(word[word.size()-1]);

  }
      
//---------------------------------------------------------------------------
// streamUnitListsIn(s,term_handling,top_tokens,bot_tokens)
//
// Stream input of Unit lists.
// The numerator list is returned in top_tokens, and the denominator list
// is returned in bot_tokens.  The basic structure has to be n1/d1 or
// n1 n2 n3/(d1 d2 d3).
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
// term_handling determines if individual unit terms will be checked against
// the list of known unit terms.
//---------------------------------------------------------------------------

// should not throw any exceptions
void Units::streamUnitListsIn(istream& s, term_handlingE term_handling,
  list<string>& top_tokens, list<string>& bot_tokens) 
  {
  top_tokens.clear();
  bot_tokens.clear();

  // First read in the numerator terms
  char last_c = streamUnitTermsIn(s,term_handling,top_tokens);    

  if (s)
    {  // There is more input to process after the numerator terms
    // Check for a '/' which introduces denominator terms.
    // More than one denominator terms will be enclosed by ().
    char c = 0;
    streampos sp = s.tellg();  // remember position of '/' (if present)
    s.get(c);
    if (!s) return;
    if (c == '/')
      {  // we have a denominator
      c = ' ';
      while (isspace(c))
        {  // skip white space after '/'
        if (!(s >> c))
          {  // the '/' is the last character, so no units left
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        }

      if (c == '(')
        {  // we have a multi-term denominator inside ()
        last_c = streamUnitTermsIn(s,term_handling,bot_tokens);    
        if (!s)
          {  // no more input or invalid input
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        s >> c;
        if ((s && c != ')') || (!s && last_c != ')'))
          {
          // no more input or there is more input, but it isn't the needed ')'
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        }
      else
        {  // a simple denominator term
        s.seekg(sp+std::streamoff(1));                 // reposition after '/'
// Alternatives tried while porting to gcc 3.4.  Some may still work
//        s.seekg(sp+(int)1);                 // reposition after '/'
//        s.seekg(sp);                 // reposition at '/'
//        s.seekg(1,std::ios::cur);    // step over the '/'
//        s.get(c);                    // step over the '/'
//        if (!s) return;
        streamUnitTermsIn(s,term_handling,bot_tokens);    
        }
      }
    }
  
  }

//---------------------------------------------------------------------------
// streamUnitsIn(s)
//
// Stream input of Units (into *this).
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
//---------------------------------------------------------------------------

// should not throw any exceptions
void Units::streamUnitsIn(istream& s) 
  {
  list<string> top_tokens;
  list<string> bot_tokens;

  // Read in the numerator and denominator unit lists
  streamUnitListsIn(s,known_terms_only,top_tokens,bot_tokens);

  if (top_tokens.empty())
    {  // no valid unit string and not eof, so set failed
    s.clear(std::ios::failbit);
    }
  else
    {  // we have some legitimate unit string, so convert it
    // Here is where *this is actually set.
    tokensToUnits(top_tokens,bot_tokens,units_top_,units_bot_);
    }
  }

//------------------------------------------------------------------------
// removeSuppressedUnit(ulist)
//
// Remove all instances of the suppressed unit from the input unit list.
//------------------------------------------------------------------------

void Units::removeSuppressedUnit(list<Unit>& ulist)
  {
  remove(ulist.begin(),ulist.end(),Unit::suppressed_unit());
//  for (list<Unit>::iterator p = ulist.begin(); p != ulist.end(); ++p)
//   {
//    if (Unit::suppressed_unit() == *p) ulist.erase(p);
//    }
  }

//------------------------------------------------------------------------
// Initialize unit tables
// These static member functions provide construct on first use semantics.
//------------------------------------------------------------------------

// should not throw any exceptions
Unit::UNIT_TABLE& Unit::unit_table() 
  {
  static Unit::UNIT_TABLE* p = new Unit::UNIT_TABLE;
  return(*p);
  }

//should not throw any exceptions
Unit::UNIT_DERIVED& Unit::unit_derived() 
  {
  static Unit::UNIT_DERIVED* p = new Unit::UNIT_DERIVED;
  return(*p);
  }

Unit& Unit::suppressed_unit() 
  {
  // This assumes that the unit_table has been constructed already
  static Unit* p = new Unit(unit_table().find("rad"));
  return(*p);
  }

//----------------------------------------//
// I/O Supporting functions and operators //
//----------------------------------------//

//----------------------
// Stream output of Unit
//----------------------

ostream& operator<<(ostream& s, const Unit& u)
  {
  return s << u.utp_->first;
  }

//--------------------------------------------------------------
// Stream output of Units (string representation is written out)
//--------------------------------------------------------------

ostream& operator<<(ostream& s, const Units& u)
  {
  string ustr;
  u.toString(ustr);
  s << ustr;
  return s;
  }

//---------------------------------------------------------------------------
// Stream input of Units.
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
//---------------------------------------------------------------------------

istream& operator>>(istream& s, Units& u)
  {
  u.streamUnitsIn(s);
  return(s);
  }

//-------------------------
// Class UnitTable Methods
//-------------------------

// Construction

//---------------------------------------------------------------------------
// The default constructor is the only constructor.
// The default constructor is also the only method added by UnitTable.
// It serves to load data from static data arrays into the map from which
// UnitTable is derived.
// All other methods used to access the unit tables are normal map methods.
//---------------------------------------------------------------------------

UnitTable::UnitTable()
  {
  Unit base_unit;  // holds base unit for each category.
  bool base_unit_set = false;

  //-------------------------------------------------------
  // Loop over list of units (grouped by category)
  //-------------------------------------------------------

  for (unsigned int i=0; i < N_units_; ++i)
    {
    string unit_name = unit_name_[i];  // convert to c++ string

    // Build basic unit conversion structure using data from static tables.
    conversion c;
    c.conv_factor = conv_factor_[i];
    c.category_index = category_index_[i];
    c.category_name = category_name_[c.category_index];
    if (i > 0 && category_index_[i] != category_index_[i-1])
      {  // starting a new category (assuming categories grouped together)
      base_unit_set = false;
      }
    (*this)[unit_name] = c;  // add table entry

    // Set base unit
    if (c.conv_factor == 1.0)
      {
      base_unit_set = true;
      UTI utp = this->find(unit_name);
      if (utp == this->end())
        {
        Unit::UnitError e("Unit table error");
	e.throwMe();
        }
      base_unit = Unit(utp);
      }
    else if (!base_unit_set)
      {
      Unit::UnitError e("Base unit needs to be first in a category");
      e.throwMe();
      }

    // Put base unit in unit table
    UTI utp = this->find(unit_name);
    utp->second.base_unit = base_unit;
    }
  }

//-------------------------------------------------------
// Static data lists to build unit table mappings from.
//-------------------------------------------------------

const double UnitTable::conv_factor_[] =
  {
  1.0,1e-3,1e-5,1e-6,1e-9,1e-12,1e-13,
  1e-2,1e3,149598073,9.460753090819e12,3.0856780e13,
  1.0,60,3600,86400,31536000,0.001,1e-6,1e-9,
  1.0,1e-3,1e-6,1e-9,1e-12,1.66053873e-27,
  1.0,
  1.0,1e-3,1e-6,0.0174532925199,2.90888208665e-4,4.84813681108e-6,6.28318530718,
  1.0,1000,
  1.0,1e-3,1e-6,1e-9,1e3,1e6,1e9
  };
const char* UnitTable::unit_name_[] =
  {
  "km","m","cm","mm","um","nm","angstrom",
  "Dm","Mm","au","ly","parsec",
  "s","min","hr","day","yr","ms","us","ns",
  "kg","g","mg","ug","ng","amu",
  "K",
  "rad","mrad","urad","deg","arcmin","arcsec","revrad",
  "mol","kmol",
  "A","mA","uA","nA","KA","MA","GA"
  };
const int UnitTable::category_index_[] =
  {
  0,0,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,
  2,2,2,2,2,2,
  3,
  4,4,4,4,4,4,4,
  5,5,
  6,6,6,6,6,6,6
  };
const char* UnitTable::category_name_[] =
  {
  "length","time","mass",
  "temperature","angle","atomic_count","electric_current"
  };

//--------------------------------
// Class UnitDerivedTable Methods
//--------------------------------

// Construction

//---------------------------------------------------------------------------
// The default constructor is the only constructor.
// The default constructor is also the only method added by UnitDerivedTable.
// It serves to load data from static data arrays into the map from which
// UnitDerivedTable is derived.
// All other methods used to access the unit tables are normal map methods.
//---------------------------------------------------------------------------

UnitDerivedTable::UnitDerivedTable()
  {
  for (unsigned int i=0; i < N_derived_; ++i)
    {  // Loop over list of derived units
    string derived_unit_name = derived_unit_name_[i];  // convert to c++ string
    unit_tokens tokens;
    ISTRINGSTREAM ist(derived_unit_string_[i]);
    Units::streamUnitListsIn(ist, Units::allow_any_terms,
      tokens.top_tokens, tokens.bot_tokens);
    (*this)[derived_unit_name] = tokens;  // add to map
    }
  }

const char* UnitDerivedTable::derived_unit_name_[] =
  {
  "Hz","KHz","MHz","GHz",
  "hectare",
  "liter",
  "rps","rpm",
  "N",
  "dyne",
  "P","kPa","bar","mbar",
  "J","mJ","uJ","nJ",
  "erg","KJ","MJ","GJ",
  "W","mW","uW","nW","KW","MW","GW",
  "C",
  "GV","MV","KV","V","mV","uV","nV",
  "F","mF","uF","nF",
  "ohm","Kohm","Mohm",
  "S",
  "Wb",
  "T","mT","uT","nT","gauss","mgauss",
  "H"
  };
const char* UnitDerivedTable::derived_unit_string_[] =
  {
  "1/s","1/ms","1/us","1/ns",
  "Dm km",
  "1/(m m mm)",
  "revrad/s","revrad/min",
  "kg m/(s s s)",
  "g cm/(s s s)",
  "N/(m m)","N/(mm m)","N/(cm mm)","N/(m cm)",
  "kg m m/(s s)","kg mm m/(s s)","kg um m/(s s)","kg nm m/(s s)",
  "g cm cm/(s s)","kg m km/(s s)","kg km km/(s s)","kg km Mm/(s s)",
  "J/s","mJ/s","uJ/s","nJ/s","KJ/s","MJ/s","GJ/s",
  "A s",
  "GW/A","MW/A","KW/A","W/A","mW/A","uW/A","nW/A",
  "C/V","C/KV","C/MV","C/GV",
  "V/A","V/mA","V/uA",
  "1/ohm",
  "V s",
  "Wb/(m m)","Wb/(km m)","Wb/(km km)","Wb/(km Mm)","Wb/(km Dm)","mV s/(km Dm)",
  "Wb/A"
  };
