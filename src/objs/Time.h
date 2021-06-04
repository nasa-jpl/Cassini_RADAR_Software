//==============================================================================
// Time.h
//
// This file contains the Time class declaration.
// The Time class provides a convenient interface to the various spice time
// formats and handles conversions between them.
// This header comment summarizes the interface.
// For details about a specific function, look at the declarations in this file.
//
// The NAIF SPICE toolkit forms the foundation of this time handling class.
// therefore the appropriate SPICE kernel files need to be loaded
// before any time conversions can be handled.  A static member function
// called spiceLoad is provided by the Frame class to load kernel files.
//
// Interface summary:
//
// class Time;
//   Defines a moment in time.
//
// Time Methods:
//
//   Construction:
//
//   Time();
//   Time(utc_str);   Initialize with UTC time string
//   Time(long sclk); Initialize with spacecraft clock count
//   Time(Uvar et);   Initialize with ephemeris time
//   Time(double et); Initialize with ephemeris time
//
//   Internal operators: +=, -=, *=, /=
//
//   Predicates:
//
//   Get/Set methods:
//
//   setUtc(utc_str);
//   setSclk(sclk);
//   setEt(et);
//   setCurrentUnixTime();
//   string utc(); Get UTC string representation
//   long sclk();  Get spacecraft clock count
//   double et();  Get ephemeris time
//   
//
// Other supporting functions:
//
// External Binary operators for Time: +, -, *, /
// ASCII stream I/O operators: <<, >>
//==============================================================================

#ifndef Time_H
#define Time_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include "Error.h"
#include "Units.h"

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//----------------------
// Forward declarations
//----------------------

class Time;

using std::string;
using std::ostream;
using std::istream;

//-------------------------
// Class Time declaration
//-------------------------

class Time : public UnitVar
  {
  public:

  enum errorE {unspecified, read_error, write_error};

  //--------------
  // Construction
  //--------------

  Time();
  Time(const string& utc_str);
  Time(const string& sc_name, const unsigned int& sclk);
  Time(const Time& t); // No exceptions // copy constructor
  Time(const UnitVar& u);  // conversion constructor
  ~Time();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  // Assignment
  Time& operator=(const Time& t); // No excpetions

  // unary -
  Time operator-() const;

  Time& operator+=(const Time& arg) ;
  Time& operator-=(const Time& arg) ;
  Time& operator*=(const Uvar& arg) ;
  Time& operator/=(const Uvar& arg) ;

  //------------
  // Predicates
  //------------

  bool valid();
  bool operator==(const Time& b) const ;
  bool operator!=(const Time& b) const ;
  bool operator>=(const Time& b) const ;
  bool operator<=(const Time& b) const ;
  bool operator>(const Time& b) const ;
  bool operator<(const Time& b) const ;

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const Time& v);
  friend istream& operator>>(istream& s, Time& v);
  void show() const; // NO exceptions

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  //----------------
  // Get/Set methods
  //----------------

  void setUtc(const string& utc_str); // No exceptions
  void setSclk(const string& sc_name, unsigned int sclk);
  void setEncodedSclk(const string& sc_name, SpiceDouble sclk);
  void setEt(const Uvar& ephemtime);
 
  string utc(const string& format) const;
  unsigned int sclk(const string& sc_name);
  SpiceDouble encodedSclk(const string& sc_name);
  const Uvar& et() const;
  // getEt returns with implicit units "s"
  void getEt(SpiceDouble& spice_et) const;
  static string getCurrentUtcTime();
  static string getCurrentUtcDOYTime();
  static string getCurrentLocalTime();
  private:

  // Internal representation

  bool valid_;
  bool encoded_sclk_set_;
  bool sclk_set_;
  double encoded_sclk_;
  unsigned int sclk_;
  
  };

//------------------------------
// Binary Operators for Time
//------------------------------

/**
// +
Time operator+(Time t1, const Time& t2);

// -
Time operator-(Time t1, const Time& t2);
**/

// >>
istream& operator>>(istream& s, Time& v);

// <<
ostream& operator<<(ostream& s, const Time& v);

#endif




