//----------------------------------------------------------------------------
// Units_impl2.h
//
// Class template member function definitions.
// These need to be in this .h file because g++ doesn't recognize the export
// keyword which would allow these definitions to be in a separate .C file.
// These need to be included after the regular UnitVar interface delcarations,
// but they are not considered part of the public interface, so they are
// separated into this file.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_impl2_h[] =
  "@(#) $Id: Units_impl2.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

//----------------------------
// UnitVar method definitions
//----------------------------

//--------------
// Construction
//--------------

template<class T_numtype>
UnitVar<T_numtype>::UnitVar()
  : a_(0) { }

//-------------------------------------------------------------------------
// UnitVar<T_numtype>(name)
//
// Construct from a name string.  If the name is a recognized constant,
// then the constant value (see below) is set along with the name.
// If the name is unknown, then only the name is set, while the UnitVar
// is initialized by default to 0.0 without units.
//-------------------------------------------------------------------------

template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const string& namestr)
  : name_(namestr)
  {
  if (namestr == "speed_light")
    {
    stringIn("2.997924580e8 m/s");
    }
  else if (namestr == "boltzmann_constant")
    {
    stringIn("1.3806503-e23 J/K");
    }
  else if (namestr == "gravitational_constant")
    {
    stringIn("6.673-e11 m m m/(s s kg)");
    }
  else if (namestr == "free_space_permittivity")
    {
    stringIn("8.854187818e-12 F/m");
    }
  }

template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const T_numtype& d)
  : a_(d) { }

template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const T_numtype& d, const std::string& ustr)
  throw(Unit::UnitError)
  : a_(d), u_(ustr) { }

template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const string& namestr, const T_numtype& d,
  const std::string& ustr)
  throw(Unit::UnitError)
  : name_(namestr), a_(d), u_(ustr) { }

// copy constructor
template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const UnitVar<T_numtype>& uv) throw()
  : name_(uv.name_), a_(uv.a_), u_(uv.u_) { }

// destructor
template<class T_numtype>
UnitVar<T_numtype>::~UnitVar()
  { }

//---------------------------
// Construct from binary file
//---------------------------

template<class T_numtype>
UnitVar<T_numtype>::UnitVar(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  // Order must match write below.
  file.read(a_);
  u_.read(file,usize);
  }

/**
template<class T_numtype>
template<class T_numtype2>
UnitVar<T_numtype>::UnitVar(const UnitVar<T_numtype2>& uvar)
  : a_(uvar.a_), u_(uvar.u_) { }
**/

// Testing

//----------------------------------------------------------
// selfTest
//
// Exercize the UnitVar type and verify proper operation.
// Actual work done in Units::selfTest.
//----------------------------------------------------------

template<class T_numtype>
bool UnitVar<T_numtype>::selfTest(const std::string& filename)
  {
  return(Units::selfTest(filename));
  }

//-------
// Setup
//-------

template<class T_numtype>
void UnitVar<T_numtype>::specifyUnitConversions(const std::string& filename)
  throw(Unit::UnitError)
  {
  Unit::specifyUnitConversions(filename);
  }

//-----------------
// Member Operators
//-----------------

// unary -
template<class T_numtype>
UnitVar<T_numtype>
  UnitVar<T_numtype>::operator-() const
  throw()
  {
  UnitVar<T_numtype> val = *this;
  val.negate();
  return(val);
  }

// +=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator+=(const T_numtype& arg)
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw("Unit mismatch performing += for UnitVar:" + name(),
      Unit::UnitError(Unit::mismatch));
    }
  a_ += arg;
  return(*this);
  }

// +=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator+=(const UnitVar<T_numtype>& arg)
  throw(Unit::UnitError)
  {
  // Match arg's units to *this's units (if possible).
  double c;
  try
    {
    c = u_.matchConversion(arg.u_);
    }
  catch(Unit::UnitError& e)
    {
    if (name() != "")
      {
      throw Unit::UnitError("[" + name() + "+=]" + e.msg, e.error_type);
      }
    else
      {
      throw;
      }
    }
  a_ += c * arg.a_;
  return(*this);
  }

// -=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator-=(const T_numtype& arg)
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw("Unit mismatch performing -= for UnitVar:" + name(),
      Unit::UnitError(Unit::mismatch));
    }
  a_ -= arg;
  return(*this);
  }

// -=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator-=(const UnitVar<T_numtype>& arg)
  throw(Unit::UnitError)
  {
  // Match arg's units to *this's units (if possible).
  double c;
  try
    {
    c = u_.matchConversion(arg.u_);
    }
  catch(Unit::UnitError& e)
    {
    if (name() != "")
      {
      throw Unit::UnitError("[" + name() + "-=]" + e.msg, e.error_type);
      }
    else
      {
      throw;
      }
    }
  a_ -= c * arg.a_;
  return(*this);
  }

// *=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator*=(const T_numtype& arg)
  throw(Unit::UnitError)
  {
  a_ *= arg;
  return(*this);
  }

// *=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator*=(const UnitVar<T_numtype>& arg)
  throw(Unit::UnitError)
  {
  // Merge arg's units with *this's units.
  double c = u_.factorMerge(arg.u_,Units::product);
  a_ *= c * arg.a_;
  return(*this);
  }

// /=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator/=(const T_numtype& arg)
  throw(Unit::UnitError)
  {
  a_ /= arg;
  return(*this);
  }

// /=
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator/=(const UnitVar<T_numtype>& arg)
  throw(Unit::UnitError)
  {
  // Merge arg's units with *this's units.
  double c = u_.factorMerge(arg.u_,Units::quotient);
  a_ *= c / arg.a_;
  return(*this);
  }

// Assignment
template<class T_numtype>
UnitVar<T_numtype>&
  UnitVar<T_numtype>::operator=(const UnitVar<T_numtype>& arg)
  throw()
  {
  // Self assignment is harmless.
  // Assignment leaves *this's name unchanged, but copies everything else.
  a_ = arg.a_;
  u_ = arg.u_;
  return(*this);
  }

//------------
// Predicates
//------------

// ==
template<class T_numtype>
bool UnitVar<T_numtype>::operator==(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  if (hasUnits() && b.hasUnits())
    {
    // Match b's units to this's units (if possible).
    double c;
    try
      {
      c = u_.matchConversion(b.u_);
      }
    catch(Unit::UnitError& e)
      {
      if (name() != "")
        {
        throw Unit::UnitError("[" + name() + "==]" + e.msg, e.error_type);
        }
      else
        {
        throw;
        }
      }
    return(a_ == c * b.a_);
    }
  else if (!hasUnits() && !b.hasUnits())
    {
    return(a_ == b.a_);
    }
  else
    {
    return(false);
    }
  }

// !=
template<class T_numtype>
bool UnitVar<T_numtype>::operator!=(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  return(!(UnitVar<T_numtype>::operator==(b)));
  }

// <
template<class T_numtype>
bool UnitVar<T_numtype>::operator<(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  if (hasUnits() && b.hasUnits())
    {
    // Match b's units to this's units (if possible).
    double c;
    try
      {
      c = u_.matchConversion(b.u_);
      }
    catch(Unit::UnitError& e)
      {
      if (name() != "")
        {
        throw Unit::UnitError("[" + name() + "<]" + e.msg, e.error_type);
        }
      else
        {
        throw;
        }
      }
    return(a_ < c * b.a_);
    }
  else if (!hasUnits() && !b.hasUnits())
    {
    return(a_ < b.a_);
    }
  else
    {
    return(false);
    }
  }

// >
template<class T_numtype>
bool UnitVar<T_numtype>::operator>(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  if (hasUnits() && b.hasUnits())
    {
    // Match b's units to this's units (if possible).
    double c;
    try
      {
      c = u_.matchConversion(b.u_);
      }
    catch(Unit::UnitError& e)
      {
      if (name() != "")
        {
        throw Unit::UnitError("[" + name() + ">]" + e.msg, e.error_type);
        }
      else
        {
        throw;
        }
      }
    return(a_ > c * b.a_);
    }
  else if (!hasUnits() && !b.hasUnits())
    {
    return(a_ > b.a_);
    }
  else
    {
    return(false);
    }
  }

// >=
template<class T_numtype>
bool UnitVar<T_numtype>::operator>=(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  return(!(UnitVar<T_numtype>::operator<(b)));
  }

// <=
template<class T_numtype>
bool UnitVar<T_numtype>::operator<=(const UnitVar<T_numtype>& b) const
  throw(Unit::UnitError)
  {
  return(!(UnitVar<T_numtype>::operator>(b)));
  }

// hasUnits
template<class T_numtype>
bool UnitVar<T_numtype>::hasUnits() const throw()
  {
  return(u_.hasUnits());
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

template<class T_numtype>
void UnitVar<T_numtype>::write(const FileMgr& file, int usize) const
  throw(Unit::UnitError, FileMgr::IoError)
  {
  file.write(a_);
  u_.write(file,usize);
  }

template<class T_numtype>
void UnitVar<T_numtype>::read(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  file.read(a_);
  u_.read(file,usize);
  }

//--------------------------------------------------------------------------
// readFloat(file,usize), writeFloat(file,usize)
//
// These methods read/write a float from/to the file regardless of the
// underlying type of UnitVar.  If the UnitVar is an int,
// then a conversion from/to float to/from int will occur.
//--------------------------------------------------------------------------

template<class T_numtype>
void UnitVar<T_numtype>::readFloat(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  float a;
  file.read(a);
  a_ = (T_numtype)a;    // convert float to underlying UnitVar type
  u_.read(file,usize);
  }

template<class T_numtype>
void UnitVar<T_numtype>::writeFloat(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  float a = (float)a_;  // convert underlying UnitVar type to float
  file.write(a);
  u_.write(file,usize);
  }

//--------------------------------------------------------------------------
// read(file,usize,ustr),write(file,usize,ustr)
//
// These read/write options will attempt to coerce the UnitVar to have the
// units in ustr (after being read from the file, or before being written
// to the file).
//--------------------------------------------------------------------------

template<class T_numtype>
void UnitVar<T_numtype>::read(const FileMgr& file, int usize,
  const std::string& ustr)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  file.read(a_);
  u_.read(file,usize);
  coerce(ustr);
  }

// See readFloat above
template<class T_numtype>
void UnitVar<T_numtype>::readFloat(const FileMgr& file, int usize,
  const std::string& ustr)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  float a;
  file.read(a);
  a_ = a;
  u_.read(file,usize);
  coerce(ustr);
  }

// See writeFloat above
template<class T_numtype>
void UnitVar<T_numtype>::writeFloat(const FileMgr& file, int usize,
  const std::string& ustr)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  // Convert units and convert underlying UnitVar type to float
  float a = (float)getInUnits(ustr);
  file.write(a);
  u_.write(file,usize);
  }

// streamIn
template<class T_numtype>
std::istream& UnitVar<T_numtype>::streamIn(std::istream& s)
  throw()
  {
  T_numtype a;
  if (!(s >> a)) return s;
  Units u;
  // Units don't have to be present, hence no return here.
  if (!(s >> u)) s.clear();
  a_ = a;
  u_ = u;
  return s;
  }

// streamOut
template<class T_numtype>
std::ostream& UnitVar<T_numtype>::streamOut(std::ostream& s)
  throw()
  {
  s << a_ << " " << u_;
  return s;
  }

//------------------------------------------------
// stringIn(uvar_str)
//
// Build from string (including value and units).
//------------------------------------------------

template<class T_numtype>
void UnitVar<T_numtype>::stringIn(const std::string& uvar_str)
  throw(Unit::UnitError)
  {
  ISTRINGSTREAM uvar_str_stream(uvar_str.c_str());
  if (!(uvar_str_stream >> a_))
    {
    throw Unit::UnitError("Bad unit number in string", Unit::formation);
    }
  uvar_str_stream >> u_;
  // EOF is ok, but other failures are not.
  // the bad() method is improperly implemented in cygwin-egcs so we have to
  // check good() and eof() instead.
  if (!(uvar_str_stream.good() || uvar_str_stream.eof()))
    {
    throw Unit::UnitError("Bad unit string in stringIn",
      Unit::formation);
    }
  }

//--------------------------
// Other arithmetic methods
//--------------------------

// pow

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::pow(const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  if (b.hasUnits())
    {
    throw Unit::UnitError("Exponent can't have units in " +
      name() + ".pow(" + b.name() + ")", Unit::mismatch);
    }
  return(pow(b.a_));  // call number pow()
  }

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::pow(const T_numtype& b)
  throw(Unit::UnitError)
  {
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
    throw Unit::UnitError(os.str(), Unit::mismatch);
    }
  else
    {  // Apply regular pow function (hopefully defined for T_numtype)
    a_ = std::pow(a_,b);
    }
  return(*this);
  }

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::pow(const int& n)
  throw()
  {  // number^integer
  a_ = std::pow(a_,n);
  if (n == 0)
    {  // A power of 0 means a unitless result.
    u_.clear();
    }
  else
    {  // Replicate the units n times.
    Units u(u_);
    for (int i=1; i < std::abs(n); ++i) u_.combine(u);
    }
  // If the power is negative then need to invert units
  if (n < 0) u_.invert();
  return(*this);
  }

// sqrt

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::sqrt()
  throw(Unit::UnitError)
  {
  try
    {
    u_.reduce(2);
    }
  catch(Unit::UnitError& e)
    {
    if (name() != "")
      {
      throw Unit::UnitError("[sqrt(" + name() + ")]" + e.msg, e.error_type);
      }
    else
      {
      throw;
      }
    }
  a_ = std::sqrt(a_);
  return(*this);
  }

// exp

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::exp()
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw Unit::UnitError("exp method needs unitless input in " +
      name() + ".exp()", Unit::mismatch);
    }
  a_ = std::exp(a_);
  return(*this);
  }

// log

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::log()
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw Unit::UnitError("log method needs unitless input in " +
      name() + ".log()", Unit::mismatch);
    }
  a_ = std::log(a_);
  return(*this);
  }

// abs

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::abs()
  throw()
  {
  a_ = fabs(a_);
  return(*this);
  }

// negate

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::negate()
  throw()
  {
  a_ = -a_;
  return(*this);
  }

// sin

template<class T_numtype>
T_numtype UnitVar<T_numtype>::sin() const
  throw(Unit::UnitError)
  {
  if (!hasUnits())
    {
    throw Unit::UnitError("sin(" + name() + ") needs units (rad)",
      Unit::mismatch);
    }
  return(std::sin(getInUnits("rad")));  // call number sin()
  }

// cos

template<class T_numtype>
T_numtype UnitVar<T_numtype>::cos() const
  throw(Unit::UnitError)
  {
  if (!hasUnits())
    {
    throw Unit::UnitError("cos(" + name() + ") needs units (rad)",
      Unit::mismatch);
    }
  return(std::cos(getInUnits("rad")));  // call number cos()
  }

// tan

template<class T_numtype>
T_numtype UnitVar<T_numtype>::tan() const
  throw(Unit::UnitError)
  {
  if (!hasUnits())
    {
    throw Unit::UnitError("tan(" + name() + ") needs units (rad)",
      Unit::mismatch);
    }
  return(std::tan(getInUnits("rad")));  // call number tan()
  }

// asin

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::asin()
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw Unit::UnitError("asin(" + name() + ") needs unitless input",
      Unit::mismatch);
    }
  a_ = std::asin(a_);
  coerce("rad");
  return(*this);
  }

// acos

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::acos()
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw Unit::UnitError("acos(" + name() + ") needs unitless input",
      Unit::mismatch);
    }
  a_ = std::acos(a_);
  coerce("rad");
  return(*this);
  }

// atan

template<class T_numtype>
UnitVar<T_numtype>& UnitVar<T_numtype>::atan()
  throw(Unit::UnitError)
  {
  if (hasUnits())
    {
    throw Unit::UnitError("atan(" + name() + ") needs unitless input",
      Unit::mismatch);
    }
  a_ = std::atan(a_);
  coerce("rad");
  return(*this);
  }

//---------------
// Other methods
//---------------

//----------------------------------------------------------------------------
// coerce(ustr)
//
// If the UnitVar has units, then they are forced to convert to the units
// specified by ustr. An exception is thrown if conversion is impossible
// because of a mismatch.
// If the UnitVar has no units, then it is given the specified units.
// The value of the UnitVar (T_numtype) is returned after any needed
// conversion has been applied.
//----------------------------------------------------------------------------

template<class T_numtype>
T_numtype UnitVar<T_numtype>::coerce(const std::string& ustr)
  throw(Unit::UnitError)
  {
  Units u(ustr);
  if (u_.hasUnits())
    {
    // Match units specified by ustr to *this's units (if possible).
    double c;
    try
      {
      c = u.matchConversion(u_);
      }
    catch(Unit::UnitError& e)
      {
      if (name() != "")
        {
        throw Unit::UnitError("[" + name() + ".coerce(" + ustr + ")]" + e.msg,
          e.error_type);
        }
      else
        {
        throw;
        }
      }
    a_ *= c;
    }
  u_ = u;  // match result's units with those specifed by ustr
  return(a_);
  }

//----------------------------------------------------------------------------
// getInUnits(ustr)
//
// The value of this UnitVar is returned after any needed
// conversion has been applied to express it in the units specified by ustr.
// An exception is thrown if conversion is impossible because of a mismatch.
//----------------------------------------------------------------------------

template<class T_numtype>
T_numtype UnitVar<T_numtype>::getInUnits(const std::string& ustr) const
  throw(Unit::UnitError)
  {
  Units u(ustr);
  // Match units specified by ustr to *this's units (if possible).
  double c;
  try
    {
    c = u.matchConversion(u_);
    }
  catch(Unit::UnitError& e)
    {
    if (name() != "")
      {
      throw Unit::UnitError("[" + name() + ".getInUnits(" + ustr + ")]" + e.msg,
        e.error_type);
      }
    else
      {
      throw;
      }
    }
  return(c*a_);
  }

//----------------------------------------------------------------------------
// getValue()
//
// This method returns the value regardless of any units present.
// It is intended for situations where the value is desired even though
// the unit category is not known and coerce can't be used.
// When possible, it is better to use coerce().
// This inspector method cannot fail.
//----------------------------------------------------------------------------

template<class T_numtype>
T_numtype UnitVar<T_numtype>::getValue() const throw()
  {
  return(a_);
  }

//----------------------------------------------------------------------------
// getUnits()
//
// This method returns a string representation of the units for this UnitVar.
// If no units are present, an empty string is returned.
// This inspector method cannot fail.
//----------------------------------------------------------------------------

template<class T_numtype>
string UnitVar<T_numtype>::getUnits() const throw()
  {
  string ustr;
  u_.toString(ustr);
  return(ustr);
  }

//----------------------------------------------------------------------------
// name()
//
// This method returns the name in a string.
// This inspector method cannot fail.
//----------------------------------------------------------------------------

template<class T_numtype>
string UnitVar<T_numtype>::name() const throw()
  {
  return(name_);
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
// Supporting Template functions
//------------------------------

// Binary operators for UnitVar.  Those that need access to the
// representation are friends of UnitVar.

// +
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator+(const UnitVar<T_numtype>& a, const T2& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result += b);
  }

// +
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator+(const T2& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result += b);
  }

// +
template<class T_numtype>
UnitVar<T_numtype>
  operator+(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result += b);
  }

// -
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator-(const UnitVar<T_numtype>& a, const T2& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result -= b);
  }

// -
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator-(const T2& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = b;
  return(result -= a);
  }

// -
template<class T_numtype>
UnitVar<T_numtype>
  operator-(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result -= b);
  }

// *
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator*(const T2& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result *= b);
  }

// *
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator*(const UnitVar<T_numtype>& a, const T2& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result *= b);
  }

// *
template<class T_numtype>
UnitVar<T_numtype>
  operator*(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result *= b);
  }

// /
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator/(const UnitVar<T_numtype>& a, const T2& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result /= b);
  }

// /
template<class T_numtype, class T2>
UnitVar<T_numtype>
  operator/(const T2& b, const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = b;
  return(result /= a);
  }

// /
template<class T_numtype>
UnitVar<T_numtype>
  operator/(const UnitVar<T_numtype>& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result /= b);
  }

// >>
template<class T_numtype>
std::istream& operator>>(std::istream& s, UnitVar<T_numtype>& uvar)
  {
  T_numtype a;
  if (!(s >> a)) return s;
  Units u;
  // Units don't have to be present, hence no return here.
  if (!(s >> u)) s.clear();
  uvar.a_ = a;
  uvar.u_ = u;
  return s;
  }

// <<
template<class T_numtype>
std::ostream& operator<<(std::ostream& s, const UnitVar<T_numtype>& uvar)
  {
  s << uvar.a_ << " " << uvar.u_;
  return s;
  }

//------------------------------------
// Other overloaded template functions
//------------------------------------

// pow
template<class T_numtype, class T_exptype>
UnitVar<T_numtype>
  pow(const UnitVar<T_numtype>& a, const T_exptype& b)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.pow(b));
  }

// pow
template<class T_numtype>
T_numtype
  pow(const T_numtype& a, const UnitVar<T_numtype>& b)
  throw(Unit::UnitError)
  {
  if (b.hasUnits())
    {
    throw Unit::UnitError("Exponent can't have units in pow(base,uvar)",
       Unit::mismatch);
    }
  return(pow(a,b.getValue()));
  }

// sqrt
template<class T_numtype>
UnitVar<T_numtype> sqrt(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.sqrt());
  }

// exp
template<class T_numtype>
T_numtype exp(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.exp().getValue());
  }

// log
template<class T_numtype>
T_numtype log(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.log().getValue());
  }

// abs
template<class T_numtype>
UnitVar<T_numtype>
  abs(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.abs());
  }

// sin
template<class T_numtype>
T_numtype sin(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.sin());
  }

// cos
template<class T_numtype>
T_numtype cos(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.cos());
  }

// tan
template<class T_numtype>
T_numtype tan(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.tan());
  }

// asin
template<class T_numtype>
UnitVar<T_numtype> asin(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.asin());
  }

// acos
template<class T_numtype>
UnitVar<T_numtype> acos(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.acos());
  }

// atan
template<class T_numtype>
UnitVar<T_numtype> atan(const UnitVar<T_numtype>& a)
  throw(Unit::UnitError)
  {
  UnitVar<T_numtype> result = a;
  return(result.atan());
  }

