//----------------------------------------------------------------------------
// Time.cpp
//
// This file contains method definitions for the Time handling classes
// These classes use the NAIF Spice toolkit to provide automatic handling
// of position and state vectors and the coordinate transformations between
// them.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_frame_c[] =
  "@(#) $Id: Time.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//---------------
// Other includes
//---------------

#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <string>
#include "Error.h"
#include "Time.h"
#include "Units.h"
#include "TemplateUtils.h"

using std::cout;
using std::endl;
using std::ios;
using std::setfill;
using std::setw;
using std::stringstream;

//-------------------------
// Methods for class Time
//-------------------------

//--------------
// Constructors
//--------------

//------------------------------------------------------
// Time()
//
// Default constructor initializes to zero ephemeris time.
//------------------------------------------------------

Time::Time()
  : UnitVar(0), valid_(false), encoded_sclk_set_(false), sclk_set_(false),
    encoded_sclk_(0), sclk_(0)
  {
  }

//------------------------------------------------------
// Time(utc_str)
//
// This constructor initializes from a utc time string.
//------------------------------------------------------

Time::Time(const string& utc_str)
  
  : UnitVar(0), valid_(false), encoded_sclk_set_(false), sclk_set_(false),
  encoded_sclk_(0), sclk_(0)
  {
  setUtc(utc_str);
  }

//-----------------------------------------------------------
// Time(sc_name,sclk)
//
// This constructor initializes from a sclk integer as
// supplied in RADAR telemetry.  The name of the spacecraft
// is given by sc_name.
//-----------------------------------------------------------

Time::Time(const string& sc_name, const unsigned int& sclk)
  
  : UnitVar(0), valid_(false), encoded_sclk_set_(false), sclk_set_(false), 
  encoded_sclk_(0), sclk_(0)
  {
  setSclk(sc_name,sclk);
  }

//--------------------
// Time(t)
//
// Copy constructor
//--------------------

Time::Time(const Time& t) // No exceptions
  : UnitVar(t), valid_(t.valid_), encoded_sclk_set_(t.encoded_sclk_set_),
    sclk_set_(t.sclk_set_), encoded_sclk_(t.encoded_sclk_), sclk_(t.sclk_)
  { }

//-------------------------------------------
// Time(u)
//
// Conversion constructor (UnitVar to Time)
//-------------------------------------------

Time::Time(const UnitVar& u) 
  : UnitVar(u), valid_(true), encoded_sclk_set_(false), sclk_set_(false), 
  encoded_sclk_(0), sclk_(0)
  {
  coerce("s");  // force error exception if not in time units
  }

// Destructor
Time::~Time() { }

//-----------
// Operators
//-----------

// Assignment
Time& Time::operator=(const Time& t) // No exceptions
  {
  // Self assignment is harmless
  UnitVar::operator=(t);  // handle base class assignment with its operator
  valid_ = t.valid_;
  encoded_sclk_set_ = t.encoded_sclk_set_;
  sclk_set_ = t.sclk_set_;
  encoded_sclk_ = t.encoded_sclk_;
  sclk_ = t.sclk_;
  return(*this);
  }

// unary -
Time Time::operator-() const 
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to negate invalid time."); e.throwMe();
    }
  Time val = *this;
  val.negate();
  val.encoded_sclk_set_ = false;
  val.sclk_set_ = false;
  return(val);
  }

// +=
Time& Time::operator+=(const Time& arg)
  
  {
  if(!valid_ || !arg.valid_)
    {
     ErrorMessage e("Attempt to add invalid time."); e.throwMe();
    }
  UnitVar::operator+=(arg);
  encoded_sclk_set_ = false;
  sclk_set_ = false;
  return(*this);
  }

// -=
Time& Time::operator-=(const Time& arg)
  
  {
  if(!valid_ || !arg.valid_)
    {
    ErrorMessage e("Attempt to subtract from invalid time."); e.throwMe();
    }
  UnitVar::operator-=(arg);
  encoded_sclk_set_ = false;
  sclk_set_ = false;
  return(*this);
  }

Time& Time::operator*=(const Uvar& arg)
  
  {
  if (!valid_)
    {
    ErrorMessage e("Attempt to multiply invalid time."); e.throwMe();
    }
  if (arg.hasUnits())
    {
    ErrorMessage e("Time can only multiply with unitless values"); e.throwMe();
    }
  UnitVar::operator*=(arg);
  encoded_sclk_set_ = false;
  sclk_set_ = false;
  return(*this);
  }

Time& Time::operator/=(const Uvar& arg)
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to divide invalid time."); e.throwMe();
    }
  if (arg.hasUnits())
    {
    ErrorMessage e("Time can only divide with unitless values"); e.throwMe();
    }
  UnitVar::operator/=(arg);
  encoded_sclk_set_ = false;
  sclk_set_ = false;
  return(*this);
  }

//---------------
// Predicates
//---------------

bool Time::valid(){ return(valid_);}

bool Time::operator==(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time."); e.throwMe();
    }
  return(UnitVar::operator==(arg));
  }

bool Time::operator!=(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time."); e.throwMe();
    }
  return(UnitVar::operator!=(arg));
  }

bool Time::operator>=(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  return(UnitVar::operator>=(arg));
  }

bool Time::operator<=(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  return(UnitVar::operator<=(arg));
  }

bool Time::operator<(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  return(UnitVar::operator<(arg));
  }

bool Time::operator>(const Time& arg) const
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  return(UnitVar::operator>(arg));
  }

//-----
// I/O
//-----

void Time::show() const // No exceptions
  {
  cout << "Time: " << getValue() << " " << getUnits() << endl;
  cout << "      " << "valid_ = " << valid_ << endl;
  cout << "      " << "encoded_sclk_set_ = " << encoded_sclk_set_ << endl;
  cout << "      " << "sclk_set_ = " << sclk_set_ << endl;
  cout << "      " << "encoded_sclk_ = " << encoded_sclk_ << endl;
  cout << "      " << "sclk_ = " << sclk_ << endl;
  }

//---------------
// Get/Set access
//---------------

//-------------------------------------------------------------------------
// setUtc(utc_str)
//
// This method sets from a UTC time string as defined by the Spice system.
//-------------------------------------------------------------------------

void Time::setUtc(const string& utc_str) // No exceptions
  {
  Uvar et=spice_utc_to_et(utc_str);
  UnitVar::operator=(et);
  valid_ = true;
  encoded_sclk_set_ = false;
  sclk_set_ = false;
  }

//-----------------------------------------------------------
// setSclk(sc_name,sclk)
//
// This method sets from a sclk integer as
// supplied in RADAR telemetry.  The name of the spacecraft
// is given by sc_name.
//-----------------------------------------------------------

void Time::setSclk(const string& sc_name, unsigned int sclk)
  {
  SpiceInt sc_id;

  // Get spice id for this body
  spice_target_id(sc_name, sc_id);

  Uvar et=spice_sclk_to_et(sc_id,sclk);

  UnitVar::operator=(et);
  sclk_ = sclk;
  valid_ = true;
  sclk_set_ = true;
  encoded_sclk_set_ = false;
  }

//-----------------------------------------------------------
// setEncodedSclk(sc_name,esclk)
//
// This method sets from an encoded sclk double.
// The name of the spacecraft is given by sc_name.
//-----------------------------------------------------------

void Time::setEncodedSclk(const string& sc_name, SpiceDouble esclk)
  {
  SpiceInt sc_id;

  // Get spice id for this body
  spice_target_id(sc_name, sc_id);

  *this=spice_encoded_sclk_to_et(sc_id,esclk);

//  UnitVar::operator=(et);
  encoded_sclk_ = esclk;
  valid_ = true;
  sclk_set_ = false;
  encoded_sclk_set_ = true;
  }

//-------------------------------------------------------------------------
// setEt(et)
//
// This method sets the indicated ephemeris time (secs).
//-------------------------------------------------------------------------

void Time::setEt(const Uvar& ephemtime)
  
  {
  UnitVar::operator=(ephemtime);
  valid_=true;
  sclk_set_ = false;
  encoded_sclk_set_ = false;
  coerce("s");
  }


//-----------------------------
// Get current UTC time
//-----------------------------

string Time::getCurrentUtcTime() 
{
  time_t unix_time = time(0);
  int seconds = gmtime(&unix_time) -> tm_sec;
  int minutes = gmtime(&unix_time) -> tm_min;
  int hours = gmtime(&unix_time) -> tm_hour;
  int month = gmtime(&unix_time) -> tm_mon; //0-11
  month += 1;
  int day_of_month = gmtime(&unix_time) -> tm_mday;
  int year = gmtime(&unix_time) -> tm_year;
  stringstream buf_stream(ios::app|ios::out);

  year += 1900;
  buf_stream << year << "-" <<
    std::setw(2) << std::setfill('0') << month << "-" <<
    std::setw(2) << std::setfill('0') << day_of_month << "T" <<
    std::setw(2) << std::setfill('0') << hours << ":" <<
    std::setw(2) << std::setfill('0') << minutes << ":" <<
    std::setw(2) << std::setfill('0') << seconds << ".000";
  return(buf_stream.str());
}  

//--------------------------------------------------
// getCurrentUtcDOYTime()
//
// Get current UTC time in ISOD (day-of-year) format
//--------------------------------------------------
string Time::getCurrentUtcDOYTime() 
{
  time_t unix_time = time(0);
  int seconds = gmtime(&unix_time) -> tm_sec;
  int minutes = gmtime(&unix_time) -> tm_min;
  int hours = gmtime(&unix_time) -> tm_hour;
  int day_of_year = gmtime(&unix_time) -> tm_yday + 1;
  int year = gmtime(&unix_time) -> tm_year;
  stringstream buf_stream(ios::app|ios::out);

  year += 1900;
  buf_stream << year << "-" <<
    std::setw(3) << std::setfill('0') << day_of_year << "T" <<
    std::setw(2) << std::setfill('0') << hours << ":" <<
    std::setw(2) << std::setfill('0') << minutes << ":" <<
    std::setw(2) << std::setfill('0') << seconds << ".000";
  return(buf_stream.str());
}  

//----------------------------------------------
//Get current local time
//----------------------------------------------
string Time::getCurrentLocalTime() 
  {
  time_t unix_time = time(0);
  int seconds = localtime(&unix_time) -> tm_sec;
  int minutes = localtime(&unix_time) -> tm_min;
  int hours = localtime(&unix_time) -> tm_hour;
  //int day_of_year=localtime(&unix_time) ->tm_yday;
  int month = localtime(&unix_time) ->tm_mon;//0-11
  month +=1;
  int day_of_month = localtime(&unix_time) -> tm_mday;
  int year = localtime(&unix_time)-> tm_year;
  year += 1900;
  string time_string = toStr(year)+"-"
    +toStr(month)+"-"+toStr(day_of_month)+"T"+toStr(hours)
    +":"+toStr(minutes)+":"+toStr(seconds)+".000";
  return(time_string);
  }
//------------------------------------------------------------------
// utc(format)
//
// Return time as a UTC time string using the indicated SPICE format.
// This method always uses a precision of 3 (fractional seconds).
//------------------------------------------------------------------

string Time::utc(const string& format) const 
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  string utc_str=spice_et_to_utc(et(),format);
  return(utc_str);
  }

//-----------------------------------------------------------
// esclk = encodedSclk(sc_name)
//
// Return time as SPICE compatible continuously encoded sclk.
// The underlying spice funtion used is sce2c_c.
// Encoded sclk is used with C-kernels.
//-----------------------------------------------------------

SpiceDouble Time::encodedSclk(const string& sc_name)
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  if (encoded_sclk_set_) return(encoded_sclk_);

  SpiceInt sc_id;

  // Get spice id for this body
  spice_target_id(sc_name,sc_id);

  encoded_sclk_=spice_et_to_encoded_sclk(sc_id,et());
  encoded_sclk_set_ = true;
  return(encoded_sclk_);
  }

//-----------------------------------------------------------
// sclk = Sclk(sc_name)
//
// Return time as Cassini Telemtry format (unsigned int) sclk.
// The underlying spice function used is scdecd_c
// Telemtry sclk is used for time segmenting L1A files
//-----------------------------------------------------------

unsigned int Time::sclk(const string& sc_name)
  
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  if (sclk_set_) return(sclk_);

  if (!encoded_sclk_set_) encodedSclk(sc_name);
  SpiceInt sc_id;

  // Get spice id for this body
  spice_target_id(sc_name,sc_id);


  sclk_=spice_decode_sclk(sc_id,encoded_sclk_);
  sclk_set_ = true;
  return(sclk_);
  }

//------------------------------------------------------------------
// et()
//
// Return ephemeris time in Uvar.
//------------------------------------------------------------------

const Uvar& Time::et() const 
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time"); e.throwMe();
    }
  return(*this);
  }

//------------------------------------------------------------------
// getEt(spice_et)
//
// Return ephemeris time in seconds using Spice compatible storage.
//------------------------------------------------------------------

void Time::getEt(SpiceDouble& spice_et) const
  {
  if(!valid_)
    {
    ErrorMessage e("Attempt to use invalid time");
    e.throwMe();
    }
  spice_et = get_in_base_units(*this);
  }

//--------------------------------------
// Supporting binary operators for Time
//--------------------------------------

/**
// +
Time operator+(Time t1, const Time& t2)
  
  {
  return(t1 += t2);
  }

// -
Time operator-(Time t1, const Time& t2)
  
  {
  return(t1 -= t2);
  }
**/

ostream& operator<<(ostream& s, const Time& t)
  {
  s << t.et();
  return(s);
  }

