//==============================================================================
// Frame.h
//
// This file contains the interface classes and functions for automatic
// coordinate system handling and geometric calculation support.
// The public interface consists of the classes
// Frame, PositionVector, StateVector, Rotation, TargetGeom, and some
// associated functions.
// This header comment summarizes the interface.
// For details about a specific function, look at the declarations in this file
// or the function leader comment in the .cpp file.
//
// The NAIF SPICE toolkit forms the foundation of this automatic frame handling
// system, therefore the appropriate SPICE kernel files need to be loaded
// before any frame conversions can be handled.  A static member function
// called spiceLoad is provided to load kernel files. Another static member
// function config is provided to perform all necessary spice loads from config
// file.
//
// Class summary:
//
// class Frame;
//   Defines a coordinate system with a specific origin and orientation.
//
// class PositionVector;
//   Represents vectors in 3-D space.  A Frame and time will be associated with
//   each PositionVector so that frame transformations can be automatically
//   performed when needed by a vector operation.
//
// class DirectionVector;
//   Represents direction in 3-D space.  A Frame and time will be associated
//   with each DirectionVector so that frame transformations can be
//   automatically performed when needed by a vector operation.
//   DirectionVectors always have unit magnitude.
//
// class FloatVector;
//   Represents a direction with a magnitude (and possibly with units).
//   A Frame and time will be associated
//   with each FloatVector so that frame transformations can be
//   automatically performed when needed by a vector operation.
//   FloatVectors are appropriate for vectors like velocity and acceleration.
//
// class StateVector;
//   Represents an object state (position and velocity) in 3-D space.
//   A Frame and time will be associated with
//   each StateVector so that frame transformations can be automatically
//   performed when needed by a vector operation.
//
// class Rotation;
//   Contains data describing an ordered rotation about the axes of a
//   coordinate frame.
//
//==============================================================================

#ifndef Frame_H
#define Frame_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include "Config.h"
#include "SpiceZfc.h"

//----------------------
// Forward declarations
//----------------------

class Frame;
class GeomError;
class PositionVector;
class DirectionVector;
class FloatVector;
class StateVector;
class Rotation;

#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Array.h"


using std::string;
using std::istream;
using std::ostream;

//-------------------------
// Class Frame declaration
//-------------------------

class Frame
  {
  //----------------------------------------------------------------
  // Static Interface - Class services as opposed to object services
  //----------------------------------------------------------------

  public:


  enum GeometryModeE { READ_FILE, PRECOMPUTE, USE_SPICE };

  //---------
  // Testing
  //---------

  static bool selfTest();

  //------
  // Setup
  //------

  static void spiceLoad(const string& filename);
  static void spiceUnLoad(const string& filename);
  static void config(Config& cfg, bool need_ckernel=false);

  //----------------
  // Static services
  //----------------

  static void loadCkernelSet(Config& cfg, bool must_find=false);
  static void unloadCkernelSet();
  static void getCkernelList(list<string>& loaded_ckernel_filenames);
  static unsigned int numCkernelsLoaded();
  static void writeGeometryToFile();
  static void adjustSpecialGeometryArrays(Config& cfg);  
  static void setGeometryMode(GeometryModeE gm);
  static void testGeometryArrays();
  static void compareGeometryArraysToFile(const string& file);
  static void cleanUp();
  static void rotationMatrix(const Frame& frame_from, const Frame& frame_to,
    SpiceDouble et, SpiceDouble m[3][3]); //No exceptions

  //--------------
  // static services direct from intermediate geometry file
  //--------------
  static double directAccessTargetRadius(const Time t);
  static void directAccessRotMatB2T(const Time t,int bn,SpiceDouble m[3][3]);
  static void directAccessRotMatT2B(const Time t,int bn,SpiceDouble m[3][3]);
  static void directAccessSCPos(const Time t, double v[3]);
  static void directAccessSCVel(const Time t, double v[3]);
  static void directAccessSCAcc(const Time t, double v[3]);
  
  private:

  //-----------------------
  // Static implementation
  //-----------------------

  static list<string> loaded_ckernel_filenames_;


  // Id Handling for nonSpice frames
  static SpiceInt next_frame_id_;

  // Target dependent strings
  static string target_;
  static string target_frame_;

  // Commonly used spice frame and origin IDs set in config
  static SpiceInt target_frame_id_;
  static SpiceInt target_origin_id_;
  static SpiceInt cassini_origin_id_;
  static SpiceInt beam_frame_ids_[5];

  // Geometry Handling mode
  static GeometryModeE geom_mode_;

  // Strings Needed for READ_FILE and PRECOMPUTE modes
  static string target_geometry_filename_;
 
  // Arrays for storing entire set of Pre-computed geometry info
  static double* quat_; // indices are 1:time, 2:beam_number, 
                             // 3: w,x,y,z w=cos(angle/2) xyz=sin(angle/2)v
  static double*  d2qdt2_; // second derivative of rotmat wrt time
  static double*  sc_vel_;   // indices are 1:time 2:Xtar,Ytar,Ztar
                             // units are km/s
  static double*  sc_pos_;   // indices are 1:time 2:Xtar,Ytar,Ztar
                             // units are km
  static double* sc_acc_; // S/C acceleration
  static double* sc_jerk_; // derivative of acceleration w.r.t time


  // Current time to which geometry data has been interpolated
  // and flags for which data has thus far been interpolated
  static double current_time_; // time at which fields have most recently been
                              // interpolated
  static bool rotmat_available_[5]; // Has the current beams rotmat been
                                       // interpolated? ...
  static bool sc_pos_available_;
  static bool sc_vel_available_;
  
  // time limits for special geometry handling
  static double start_time_;
  static double time_step_;
  static int num_time_steps_;

  // Interploated geometry data arrays
  static double* rotmat_interp_; // indices are 1:beam_number, 
                             // 2: Xbeam, Ybeam, Zbeam, 3: Xtar, Ytar, Ztar
  static double*  sc_vel_interp_;   // index is Xtar,Ytar,Ztar
                             // units are km/s
  static double* sc_pos_interp_;   // index is Xtar,Ytar,Ztar
                             // units are km
  
  //-----------------------
  // Private static methods for Special Geometry Handling
  //------------------------

  static void allocateSpecialGeometryArrays();
  static double interpolate1D(double* array, double* d2arraydt2,
			      int start_idx, int skip, double t);
  static void interpolateQuaternion(double* qs, 
				    int start_idx, int skip, 
					     double t, double* q, int debug=0);
  static void printQuaternion(double q[4]);
  static void raiseQuaternionToPower(double q[4], double a);
  static void multiplyQuaternions(double q1[4], double q2[4], double q[4]);
  static void interpolateRotMat(int beam_number);
  static void interpolateSCVel();
  static void interpolateSCPos();
  static void updateTime(const Time& t);
  static void updateTime(double t_in_s);
  static void computeSpecialGeometryArraysFromSpice();
  static void computeDerivativeArrays();
  static void readGeometryFromFile();
  
  //-------------------------
  // Normal object interface
  //-------------------------

  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  enum indexE {X, Y, Z};

  //----------------------------------------------------------------------
  // Friend classes that need access to Frame's internal methods and data
  //----------------------------------------------------------------------

  friend class PositionVector;
  friend class StateVector;

  //--------------
  // Construction
  //--------------

  // computationally efficient constructors
  // preferable for low-level calls
  Frame();  // Default constructor
  Frame(const Frame& f); // const reference constructor
  Frame(int frame_id, int origin_id); // used for spice frames only
  Frame(const PositionVector& origin,
	DirectionVector x,
	DirectionVector y,
	DirectionVector z);
  Frame(DirectionVector x,
	DirectionVector y,
	DirectionVector z);
  Frame(const Frame& frame,
    PositionVector rel_origin);

  // legacy constructors (with string operations) 
  // ease of use for high-level calls
  Frame(const string& frame_name, const string& origin_name);
  Frame(const string& frame_name, const Frame& frame,
    PositionVector rel_origin);
 
  Frame(const string& frame_name, const PositionVector& origin,
	DirectionVector x,
	DirectionVector y,
	DirectionVector z);
  Frame(const string& frame_name,
	DirectionVector x,
	DirectionVector y,
	DirectionVector z);

  Frame(const string& frame_name, const PositionVector& rel_origin,
	const Rotation& rel_rot);
  Frame(const string& frame_name, const Frame& frame,
	const Rotation& rel_rot);

  //--------------
  // Predicates
  //--------------

  bool operator==(const Frame& b) const ; //No exceptions
  bool operator!=(const Frame& b) const ; //No exceptions
  bool isSpice() const ; //No exceptions

  //------------
  // Geometry speed-up non-static methods
  //-------------

  void setSpecialFrameNumber();
  //--------------
  // Other methods
  //--------------

  // do NOT overuse these routine as they are slow.
  string name() const;
  string spiceFrameName() const;

  // legacy ephemeris call (ease of use for high level calls)
  // throws error if corr!="NONE" 
  // We should perform any necessary relativistic corrections ourselves as
  // SPICE's version may not meet our specific need!
  void ephemeris(StateVector& sephem, const string& target, const Time& t,
    const string& corr) const;

  // stream-lined ephemeris call for low level calls
  // assumes correction="NONE"
  void ephemeris(StateVector& sephem, SpiceInt spice_target_id, 
		 double time_in_s) const;
  

  // These routines are not commonly used ... may be slow
  Array1D<double> quaternion(const DirectionVector& x,
    const DirectionVector& y, const DirectionVector& z) const;
  Array1D<double> quaternion(const Frame& f, const Time& t) const;
  void axialVectors(const Frame& f, const Time& t,
    DirectionVector& x, DirectionVector& y, DirectionVector& z) const;

  private:

  //-------------------------
  // Internal methods
  //-------------------------

  void setSpiceRelativeRotation(const Rotation& rot, const Frame& base_frame);
  void setSpiceRelativeRotation(DirectionVector x, 
				DirectionVector y,
				DirectionVector z, 
				const Frame& base_frame);
  void setSpiceRelativeOrigin(const PositionVector& rel_origin);
  void setDefiningSpiceName(const Frame& frame) ; //No exceptions
  void getNewFrameID();

  //-------------------------
  // Internal representation
  //-------------------------

  bool spice_frame_;                 // true if this frame is Spice frame
                                     // (NonSpice orign also causes this to 
                                     // be false)

  SpiceInt spice_frame_id_;          // int id of defining Spice frame
  SpiceInt frame_id_;                // int id of this frame  

  SpiceInt spice_origin_id_;         // int id of defining Spice origin
  SpiceDouble spice_xform_[3][3];    // xforms to the defining Spice frame
  SpiceDouble spice_invxform_[3][3]; // xforms from the defining Spice frame
  double rel_origin_[3];               // origin in defining Spice frame
                                       // in km
  
  // -- internal representation of for special frames which make use of 
  // -- precomputed geometry

  int special_frame_number_; // -1 not a special frame 
			     //  0 target body fixed
			     //  1-5 beam 1-5 frame


  int beam_number_;
  
  };

//----------------------------------
// Class PositionVector declaration
//----------------------------------

class PositionVector
  {
  public:

  enum indexE {X, Y, Z};

  //----------------------------------------------------------------------
  // Friend classes that need access to PositionVector's internals
  //----------------------------------------------------------------------

  friend class Frame;
  friend class DirectionVector;
  friend class FloatVector;
  friend class StateVector;

  //--------------
  // Construction
  //--------------

  PositionVector(); // No exceptions
  PositionVector(const string& name); // No exceptions
  PositionVector(const string& name, const Frame& frame, const Time& t,
    const Uvar& x, const Uvar& y, const Uvar& z); // No exceptions
  PositionVector(const Frame& frame, const Time& t,
    const Uvar& x, const Uvar& y, const Uvar& z); // No exceptions
  PositionVector(const Frame& frame, const double& t,
    const double& x, const double& y, const double& z); // No exceptions

  // take frame and time from r
  PositionVector(const PositionVector& r,
    const Uvar& x, const Uvar& y, const Uvar& z); // No exceptions

  // copy constructor
  PositionVector(const PositionVector& r); // No exceptions
  // type conversion constructors
  PositionVector(const DirectionVector& u); // No exceptions
  PositionVector(const FloatVector& u); // No exceptions
  ~PositionVector();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  // unary -
  PositionVector operator-() const ; //No exceptions

  PositionVector& operator+=(const PositionVector& arg);
  PositionVector& operator+=(const Uvar& arg);
  PositionVector& operator-=(const PositionVector& arg);
  PositionVector& operator-=(const Uvar& arg);
  PositionVector& operator*=(const Uvar& arg);
  PositionVector& operator*=(const double& arg);
  PositionVector& operator/=(const Uvar& a);
  PositionVector& operator/=(const double& arg);
  PositionVector& operator*=(const DirectionVector& arg);

  // commented out because result does not have km units and is NOT a
  // position
  // PositionVector& operator/=(const PositionVector& arg);
  // void invert();

  // Assignment
  PositionVector& operator=(const PositionVector& r) ; //No exceptions

  // Indexing
  Uvar& operator[](const indexE& i) ; //No exceptions
  const Uvar& operator[](const indexE& i) const ; //No exceptions
  double km(const indexE& i) const;
  void setkm(const indexE& i, double value);
 
  //------------
  // Predicates
  //------------

  bool operator==(PositionVector b) const;
  bool operator!=(PositionVector b) const;
  bool timeMatch(const PositionVector& b) const;
  bool timeMatch(const DirectionVector& b) const;
  bool timeMatch(const FloatVector& b) const;
  bool frameMatch(const PositionVector& b) const ; //No exceptions
  bool frameMatch(const DirectionVector& b) const ; //No exceptions
  bool frameMatch(const FloatVector& b) const ; //No exceptions
  bool isSpice() const ; //No exceptions

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const PositionVector& v);
  friend istream& operator>>(istream& s, PositionVector& v);
  friend DirectionVector operator/(const PositionVector& a, 
				   const PositionVector& b);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  //---------------
  // Other methods
  //---------------

  string name() const ; //No exceptions
  Time time() const ; //No exception
  void setTime(const Time& t) ; //No exceptions
  double timeInSeconds() const ; //No exceptions
  void setTimeInSeconds(double t) ; //No exceptions
  string frameName() const ; //No exceptions
  PositionVector& representInSpiceFrame();
  PositionVector& representIn(const Frame& frame);
  PositionVector& representIn(const PositionVector& r);
  PositionVector& representIn(const DirectionVector& u);
  PositionVector& representIn(const FloatVector& u);
  Uvar magnitude()  const;
  void scale(Uvar mag);
  Uvar angle(const PositionVector& v) const;
  double angleInRads(const PositionVector& v) const;
  PositionVector& rotateAbout(const DirectionVector& v, const Uvar& angle);

  private:

  // Internal methods
  static void frameRotate(const Frame& frame_from, const Frame& frame_to,
    SpiceDouble et, SpiceDouble v[3]) ; //No exceptions
  
  // Internal representation
  string* name_;
  Frame frame_;
  double t_; // in seconds  
  Uvar v_[3];

  };

//--------------------------------------
// Binary Operators for PositionVector.
//--------------------------------------

// + - * /
PositionVector operator+(PositionVector a, const PositionVector& b);
PositionVector operator+(PositionVector a, const Uvar& b);
PositionVector operator+(const Uvar& a, PositionVector b);
PositionVector operator-(PositionVector a, const PositionVector& b);
PositionVector operator-(PositionVector a, const Uvar& b);
PositionVector operator-(const Uvar& a, PositionVector b);
PositionVector operator*(PositionVector a, const DirectionVector& b);

// Commented out because result is not a position

// PositionVector operator/(PositionVector a, const PositionVector& b);
// PositionVector operator/(const Uvar& a, PositionVector b);


PositionVector operator/(PositionVector a, const Uvar& b);
PositionVector operator*(PositionVector a, const Uvar& b);
PositionVector operator*(const Uvar& a, PositionVector b);




// >>
istream& operator>>(istream& s, PositionVector& v);

// <<
ostream& operator<<(ostream& s, const PositionVector& v);

//---------------------------
// Other Supporting functions
//---------------------------

// Cross and dot products
FloatVector cross(const PositionVector& a, PositionVector b);
Uvar dot(const PositionVector& a, PositionVector b);

// Combinations of PositionVector and DirectionVector and FloatVector

Uvar dot(const PositionVector& a, DirectionVector b);
Uvar dot(const DirectionVector& a, PositionVector b);
Uvar dot(const PositionVector& a, FloatVector b);
Uvar dot(const FloatVector& a, PositionVector b);
Uvar dot(const FloatVector& a, DirectionVector b);
Uvar dot(const DirectionVector& a, FloatVector b);

FloatVector cross(const PositionVector& a, DirectionVector b);
FloatVector cross(const DirectionVector& a, PositionVector b);
FloatVector cross(const PositionVector& a, FloatVector b);
FloatVector cross(const FloatVector& a, PositionVector b);
FloatVector cross(const FloatVector& a, DirectionVector b);
FloatVector cross(const DirectionVector& a, FloatVector b);

//----------------------------------
// Class DirectionVector declaration
//----------------------------------

class DirectionVector
  {
  public:

  enum indexE {X, Y, Z};
  
  //----------------------------------------------------------------------
  // Friend classes that need access to DirectionVector's internals
  //----------------------------------------------------------------------
  friend class Frame;
  friend class PositionVector;
  friend class FloatVector;

  //--------------
  // Construction
  //--------------

  DirectionVector() ; //No exceptions
  DirectionVector(const string& name) ; //No exceptions
  DirectionVector(const Frame& frame, const double& t_in_s,
    const double& x, const double& y, const double& z) ; //No exceptions
  DirectionVector(const Frame& frame, const Time& t,
    const double& x, const double& y, const double& z) ; //No exceptions
  DirectionVector(const string& name, const Frame& frame, const Time& t,
    const double& x, const double& y, const double& z) ; //No exceptions
  DirectionVector(const string& name, const Frame& frame, const double& t,
    const double& x, const double& y, const double& z) ; //No exceptions
  DirectionVector(const string& name, const PositionVector& r);
  DirectionVector(const string& name, const DirectionVector& r); 
  //No exceptions

  DirectionVector(const string& name, const FloatVector& r);
  // copy constructors
  DirectionVector(const DirectionVector& r) ; //No exceptions  
  DirectionVector(const PositionVector& r);  // copy from PositionVector
  DirectionVector(const FloatVector& r);  // copy from FloatVector
  ~DirectionVector();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  // unary -
  DirectionVector operator-() const ; //No exceptions

  DirectionVector& operator+=(const DirectionVector& arg) ; //No exceptions
  DirectionVector& operator-=(const DirectionVector& arg) ; //No exceptions
  DirectionVector& operator*=(const DirectionVector& arg) ; //No exceptions
  DirectionVector& operator/=(const DirectionVector& arg) ; //No exceptions

  //DirectionVector& operator+=(const double& arg) ; //No exceptions
  //DirectionVector& operator-=(const double& arg) ; //No exceptions
  //DirectionVector& operator*=(const double& arg) ; //No exceptions
  //DirectionVector& operator/=(const double& arg) ; //No exceptions

  // Assignment
  DirectionVector& operator=(const DirectionVector& u) ; //No exceptions

  // Indexing
  double& operator[](const indexE& i) ; //No exceptions
  const double& operator[](const indexE& i) const ; //No exceptions

  //------------
  // Predicates
  //------------

  bool operator==(DirectionVector b) const;
  bool operator!=(DirectionVector b) const;
  bool timeMatch(const DirectionVector& b) const;
  bool timeMatch(const PositionVector& b) const;
  bool timeMatch(const FloatVector& b) const;
  bool timeMatch(const Time& t) const;
  bool timeMatch(const double& time_in_s) const;
  bool frameMatch(const DirectionVector& b) const ; //No exceptions
  bool frameMatch(const PositionVector& b) const ; //No exceptions
  bool frameMatch(const FloatVector& b) const ; //No exceptions
  bool isSpice() const ; //No exceptions

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const DirectionVector& v);
  friend istream& operator>>(istream& s, DirectionVector& v);

  //--------------------------
  // Get/Set methods
  //--------------------------

  void getPlanetodetic(Uvar& lat, Uvar& lon) const;
  void setPlanetodetic(const Uvar& lat, const Uvar& lon);
  void getAzimuthElevation(Uvar& azi, Uvar& elev) const;
  void setAzimuthElevation(const Uvar& azi, const Uvar& elev);
  void getRADEC(Uvar& ra, Uvar& dec) const;
  void setRADEC(const Uvar& ra, const Uvar& dec);
  void getSpherical(Uvar& theta, Uvar& phi) const;
  void setSpherical(const Uvar& theta, const Uvar& phi);

  void getPlanetodeticInRads(double& lat, double& lon) const;
  void setPlanetodeticInRads(const double& lat, const double& lon);
  void getAzimuthElevationInRads(double& azi, double& elev) const;
  void setAzimuthElevationInRads(const double& azi, const double& elev);
  void getRADECInRads(double& ra, double& dec) const;
  void setRADECInRads(const double& ra, const double& dec);
  void getSphericalInRads(double& theta, double& phi) const;
  void setSphericalInRads(const double& theta, const double& phi);

  //---------------
  // Other methods
  //---------------

  string name() const ; //No exceptions
  Time time() const ; //No exceptions
  void setTime(const Time& t) ; //No exceptions
  double timeInSeconds() const ; //No exceptions
  void setTimeInSeconds(double d) ; //No exceptions
  string frameName() const ; //No exceptions
  DirectionVector& representIn(const Frame& frame);
  DirectionVector& representIn(const DirectionVector& u);
  DirectionVector& representIn(const PositionVector& r);
  DirectionVector& representIn(const FloatVector& r);
  Uvar angle(DirectionVector v) const;
  Uvar angleInRads(DirectionVector v) const;
  DirectionVector& rotateAbout(const DirectionVector& v, const Uvar& angle);
  void invert();
  double* data();

  private:

  // Internal methods
  void unity_scale(bool allow_zero = false);

  // Internal representation
  string* name_;
  Frame frame_;
  double t_; // in seconds of ephemeris time
  SpiceDouble v_[3];

  };

//--------------------------------------
// Binary Operators for DirectionVector.
//--------------------------------------

// + - * /
DirectionVector operator+(DirectionVector a, const DirectionVector& b)
  ; //No exceptions

//DirectionVector operator+(DirectionVector a, const double& b)
//; //No exceptions
//DirectionVector operator+(const double a, const DirectionVector& b)
//; //No exceptions

DirectionVector operator-(DirectionVector a, const DirectionVector& b)
  ; //No exceptions

//DirectionVector operator-(DirectionVector a, const double& b)
//  ; //No exceptions

//DirectionVector operator-(const double a, const DirectionVector& b)
//; //No exceptions

DirectionVector operator*(DirectionVector a, const DirectionVector& b)
  ; //No exceptions

//DirectionVector operator*(DirectionVector a, const double& b)
//; //No exceptions

//DirectionVector operator*(const double a, const DirectionVector& b)
//; //No exceptions

DirectionVector operator/( DirectionVector a, const DirectionVector& b)
  ; //No exceptions

//DirectionVector operator/(DirectionVector a, const double& b)
//; //No exceptions

//DirectionVector operator/(const double a, const DirectionVector& b)
//; //No exceptions

// >>
istream& operator>>(istream& s, DirectionVector& v);

// <<
ostream& operator<<(ostream& s, const DirectionVector& v);

//---------------------------
// Other Supporting functions
//---------------------------

// Cross and dot products
DirectionVector cross(const DirectionVector& a, DirectionVector b);
double dot(const DirectionVector& a, DirectionVector b);

//----------------------------------
// Class FloatVector declaration
//----------------------------------

class FloatVector
  {
  public:

  enum indexE {X, Y, Z};

  //----------------------------------------------------------------------
  // Friend classes that need access to FloatVector's internals
  //----------------------------------------------------------------------

  friend class PositionVector;
  friend class DirectionVector;
  friend class StateVector;

  //--------------
  // Construction
  //--------------

  FloatVector() ; //No exceptions
  FloatVector(const string& name) ; //No exceptions
  FloatVector(const string& name, const Frame& frame, const Time& t,
    const Uvar& x, const Uvar& y, const Uvar& z) ; //No exceptions
  FloatVector(const Frame& frame, const double& t,
    const Uvar& x, const Uvar& y, const Uvar& z) ; //No exceptions
  FloatVector(const string& name, const PositionVector& r) ; //No exceptions
  FloatVector(const string& name, const FloatVector& r) ; //No exceptions
  FloatVector(const FloatVector& r) ; //No exceptions  // copy constructor
  FloatVector(const PositionVector& r) ; //No exceptions  // copy from PositionVector
  ~FloatVector();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  // unary -
  FloatVector operator-() const ; //No exceptions

  FloatVector& operator+=(const FloatVector& arg);
  FloatVector& operator-=(const FloatVector& arg);
  FloatVector& operator*=(const FloatVector& arg);
  FloatVector& operator/=(const FloatVector& arg);

  FloatVector& operator+=(const Uvar& arg);
  FloatVector& operator-=(const Uvar& arg);
  FloatVector& operator*=(const Uvar& arg);
  FloatVector& operator/=(const Uvar& arg);

  // Assignment
  FloatVector& operator=(const FloatVector& u) ; //No exceptions

  // Indexing
  Uvar& operator[](const indexE& i) ; //No exceptions
  const Uvar& operator[](const indexE& i) const ; //No exceptions
  
  //------------
  // Predicates
  //------------

  bool operator==(FloatVector b) const;
  bool operator!=(FloatVector b) const;
  bool timeMatch(const FloatVector& b) const;
  bool timeMatch(const PositionVector& b) const;
  bool timeMatch(const DirectionVector& b) const;
  bool timeMatch(const Time& t) const;
  bool timeMatch(const double& t) const;
  bool frameMatch(const FloatVector& b) const ; //No exceptions
  bool frameMatch(const PositionVector& b) const ; //No exceptions
  bool frameMatch(const DirectionVector& b) const ; //No exceptions
  bool isSpice() const ; //No exceptions

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const FloatVector& v);
  friend istream& operator>>(istream& s, FloatVector& v);

  //--------------------------
  // Get/Set methods
  //--------------------------

  //---------------
  // Other methods
  //---------------

  string name() const ; //No exceptions
  Time time() const ; //No exceptions
  void setTime(const Time& t) ; //No exceptions
  double timeInSeconds() const ; //No exceptions
  void setTimeInSeconds(double t_in_seconds) ; //No exceptions
  string frameName() const ; //No exceptions
  FloatVector& representIn(const Frame& frame);
  FloatVector& representIn(const FloatVector& u);
  FloatVector& representIn(const PositionVector& r);
  FloatVector& representIn(const DirectionVector& r);
  Uvar magnitude() const;
  void scale(Uvar mag);
  Uvar angle(const FloatVector& v) const;
  void invert();

  private:

  // Internal representation
  string* name_;
  Frame frame_;
  double t_;
  Uvar v_[3];

  };

//--------------------------------------
// Binary Operators for FloatVector.
//--------------------------------------

// + - * /
FloatVector operator+(FloatVector a, const FloatVector& b);
FloatVector operator+(FloatVector a, const Uvar& b);
FloatVector operator+(const Uvar a, const FloatVector& b);
FloatVector operator-(FloatVector a, const FloatVector& b);
FloatVector operator-(FloatVector a, const Uvar& b);
FloatVector operator-(const Uvar a, const FloatVector& b);
FloatVector operator*(const FloatVector& a, const FloatVector& b);
FloatVector operator*(FloatVector a, const Uvar& b);
FloatVector operator*(const Uvar a, const FloatVector& b);
FloatVector operator/(const FloatVector& a, const FloatVector& b);
FloatVector operator/(FloatVector a, const Uvar& b);
FloatVector operator/(const Uvar a, const FloatVector& b);

// >>
istream& operator>>(istream& s, FloatVector& v);

// <<
ostream& operator<<(ostream& s, const FloatVector& v);

//---------------------------
// Other Supporting functions
//---------------------------

// Cross and dot products
FloatVector cross(const FloatVector& a, FloatVector b);
Uvar dot(const FloatVector& a, FloatVector b);

//----------------------------------
// Class StateVector declaration
//----------------------------------

class StateVector
  {
  public:

  enum indexE {X, Y, Z};

  //----------------------------------------------------------------------
  // Friend classes that need access to StateVector's internals
  //----------------------------------------------------------------------

  friend class PositionVector;
  friend class Frame;

  //--------------
  // Construction
  //--------------

  StateVector() ; //No exceptions
  StateVector(const string& name) ; //No exceptions
  StateVector(const string& name, const PositionVector& r,
    const FloatVector& v);
  StateVector(const PositionVector& r, const FloatVector& v);
  StateVector(const StateVector& s) ; //No exceptions  // copy constructor
  ~StateVector();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  // unary -
  StateVector operator-() const ; //No exceptions

  StateVector& operator+=(const StateVector& arg);
  StateVector& operator-=(const StateVector& arg);
  StateVector& operator*=(const StateVector& arg);
  StateVector& operator/=(const StateVector& arg);

  StateVector& operator*=(const Uvar& arg);
  StateVector& operator/=(const Uvar& arg);

  // Assignment
  StateVector& operator=(const StateVector& s);

  //------------
  // Predicates
  //------------

  bool operator==(const StateVector& b) const;
  bool timeMatch(const StateVector& b) const;
  bool timeMatch(const Time& t) const;
  bool timeMatch(const double& t) const;

  //------------
  // Get and Set
  //------------

  PositionVector position() const ; //No exceptions
  FloatVector velocity() const ; //No exceptions
  Time time() const ; //No exceptions
  void setTime(const Time& t) ; //No exceptions
  double timeInSeconds() const ; //No exceptions
  void setTimeInSeconds(double t_in_s) ; //No exceptions
  void setState(const PositionVector& r, FloatVector v);

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const StateVector& v);
  friend istream& operator>>(istream& s, StateVector& v);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  //---------------
  // Other methods
  //---------------

  string name() const ; //No exceptions
  void invert();

  private:

  // Internal methods
  
  // Internal representation
  string* name_;
  Frame frame_;
  double t_;
  Uvar r_[3];
  Uvar v_[3];

  };

//--------------------------------------
// Binary Operators for StateVector.
//--------------------------------------

// + - * /
StateVector operator+(StateVector a1, const StateVector& a2);
StateVector operator-(StateVector a1, const StateVector& a2);
StateVector operator*(StateVector a1, const StateVector& a2);
StateVector operator*(StateVector a1, const Uvar& a2);
StateVector operator*(const Uvar& a1, StateVector a2);
StateVector operator/(StateVector a1, const StateVector& a2);
StateVector operator/(StateVector a1, const Uvar& a2);
StateVector operator/(const Uvar& a1, StateVector a2);

// >>
istream& operator>>(istream& s, StateVector& v);

// <<
ostream& operator<<(ostream& s, const StateVector& v);

//----------------------------------
// Class Rotation declaration
//----------------------------------

class Rotation
  {
  public:

  enum indexE {X, Y, Z};

  //--------------
  // Construction
  //--------------

  Rotation(const Uvar& r1, const Uvar& r2, const Uvar& r3,
    unsigned int o1, unsigned int o2, unsigned int o3); //No exceptions

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  Rotation& operator+=(const Rotation& arg);
  Rotation& operator-=(const Rotation& arg);
  Rotation& operator*=(const Rotation& arg);
  Rotation& operator/=(const Rotation& arg);

  bool operator==(const Rotation& b) const;

  //------------
  // Predicates
  //------------

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<<(ostream& s, const Rotation& v);
  friend istream& operator>>(istream& s, Rotation& v);

  //---------------
  // Other methods
  //---------------

  void rotatedAxes(DirectionVector& x, DirectionVector& y,
    DirectionVector& z) const;

  private:

  // Internal methods

  // Internal representation
  Uvar r1_,r2_,r3_;
  unsigned int o1_,o2_,o3_;

  };

//----------------------------------
// Class RotationalVelocity
//----------------------------------

class RotationalVelocity
  {
  public:

  enum indexE {X, Y, Z};

  //--------------
  // Construction
  //--------------

  RotationalVelocity() // No exceptions
    {
    }
  RotationalVelocity(const Frame& fsc, const Frame& ftarget,
		     const Time& t); 
  
  RotationalVelocity(const string&, const Frame& ref_frame,
		     const Time& t, Uvar x, Uvar y, Uvar z) // No exceptions
    {
      v_[0]=x;
      v_[1]=y;
      v_[2]=z;
      t_=t;
    }

  RotationalVelocity(const RotationalVelocity& rv) // No exceptions
    {
      v_[0]=rv.v_[X];
      v_[1]=rv.v_[Y];
      v_[2]=rv.v_[Z];;
      t_=rv.t_;
    }

  void setTime(const Time& t){ t_=t; };
  //---------
  // Testing
  //---------

  static bool selfTest();

  //---------------------
  //indexing
  //---------------------
  Uvar& operator[] (const indexE& i)//no exception
  {
    return(v_[i]);
  }
  
  private:
  Uvar v_[3];
  Time t_;
  };

//------------------------------------------------------
// Error handling class for Geometry related classes
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//------------------------------------------------------

class GeomError : public ErrorMessage
  {
  public:

  enum errorE {unspecified, scale_zero, time_mismatch, spice_nosolution,
               read_error, write_error, ambiguity_error, spice_bad_ckernel,
               spice_misc, spice_bad_target, spice_bad_target_radii};

  //--------------
  // Constructors
  //--------------
  GeomError(errorE err_type = unspecified,double et=0);
  GeomError(const string& emsg, errorE err_type = unspecified,double et=0);

  // Spice Error Handling Routine
  void initSpiceError();

  // Routine for self throw
  void throwMe();

  // Public type flag
  errorE error_type;
  double error_time; // ephemeris time in seconds

  protected:
  static bool spice_handling_enabled_;
  };

//-------------------------------
// Supporting external functions
//-------------------------------

void spice_surfpt(PositionVector& rsurf, bool& found, PositionVector observer,
  DirectionVector ulook, const PositionVector& target_radii);
void spice_ringplane_intercept(PositionVector& rint, bool& found,
  PositionVector observer, DirectionVector ulook,
  const PositionVector& target_radii);
void spice_nearpt(PositionVector& rnadir, Uvar& alt,
  PositionVector observer, const PositionVector& target_radii);
void spice_npedln(const Uvar& rx, const Uvar& ry, const Uvar& rz,
  const PositionVector& rpoint, DirectionVector udir,
  PositionVector& rnear, Uvar& dist);
void spice_target_value(const string& name, const string& value_name,
  Uvar& value);
void spice_target_id(const string& name, SpiceInt& id);
void spice_frame_id(const string& name, SpiceInt& id);
void spice_frame_name(SpiceInt id, string& name);
void spice_frame_name(SpiceInt id, char* name, int length);

void spice_target_radii(const string& name, const Frame& frame, const Time& t,
			SpiceInt id, 
			PositionVector& radii);
void spice_target_radii(const string& name, Uvec& radii);
Uvar spice_utc_to_et(const string& utc_str);
bool isUtc(const string& utc_str);
Uvar spice_sclk_to_et(SpiceInt sc_id,unsigned int sclk);
string spice_et_to_utc(const Uvar& et,const string& format);
double spice_et_to_encoded_sclk(SpiceInt sc_id, const Uvar& et);
Time spice_encoded_sclk_to_et(SpiceInt sc_id, double esclk);
unsigned int spice_decode_sclk(SpiceInt sc_id, double esclk);
void spice_write_ckernel(const string& filename, unsigned int Ncomchar, 
			 SpiceDouble begtime,
			 SpiceDouble endtime, SpiceInt inst_id, 
			 string ref_string, SpiceBoolean ang_vel_flag, 
			 string seg_id_string, int Nrecords, 
			 SpiceDouble spice_sclkdp[], 
			 SpiceDouble spice_quats[][4],
			 SpiceDouble spice_avvs[][3], 
			 SpiceDouble spice_start_x_time[]);


#endif







