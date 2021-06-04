// Class Methods for OblCylProj 

static const char rcs_id_projections_c[] =
  "@(#) $Id: Projections.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include<string>
#include"SARFunctions.h"
#include"Projections.h"
#include"Constants.h"
#include"DebugInfo.h"

using std::string;
using std::endl;

//--------------
// Constructors
//--------------

//------------------------------------------------------------------
// Trivial (standard) projection
//------------------------------------------------------------------
OblCylProj::OblCylProj()
  : is_trivial(true)
{
}

//------------------------------------------------------------------
// OblCylProj()
//
// Set up an oblique cylindrical projection in which the spacecraft
// position in s denotes the intersection of the equator and the
// prime meridian, and the spacecraft velocity direction defines the
// equator.  s is presumed to be the state of Cassini in the Target
// Body Fixed frame at the time of closest approach.  Note that the
// direction of the y-axis no longer depends on the look direction.
//------------------------------------------------------------------
OblCylProj::OblCylProj(StateVector& s)
{
  DebugInfo dbg("OblCylProj::OblCylProj");
  
  DirectionVector basis_x("", s.position());
  DirectionVector basis_z("", cross(basis_x, s.velocity()));
  DirectionVector basis_y("", cross(basis_z,basis_x));


  if(dbg.level){
    dbg.file << "Debug Output for OblCylProj constructor ..." << endl;
    dbg.file << "Closest approach S/C position direction vector=" << 
      DirectionVector("", s.position()) << endl;
    dbg.file << "Closest approach S/C velocity direction vector=" << 
      DirectionVector("", s.velocity()) << endl;
    dbg.file << "Basis X vector=" << basis_x << endl;
    dbg.file << "Basis Y vector=" << basis_y << endl;
    dbg.file << "Basis Z vector=" << basis_z << endl;
    dbg.file << "End Debug Output for OblCylProj constructor" << endl;
  }
  // Populate the rotation matrix
  rotation_matrix_[0][0] = basis_x[DirectionVector::X];
  rotation_matrix_[0][1] = basis_x[DirectionVector::Y];
  rotation_matrix_[0][2] = basis_x[DirectionVector::Z];
  rotation_matrix_[1][0] = basis_y[DirectionVector::X];
  rotation_matrix_[1][1] = basis_y[DirectionVector::Y];
  rotation_matrix_[1][2] = basis_y[DirectionVector::Z];
  rotation_matrix_[2][0] = basis_z[DirectionVector::X];
  rotation_matrix_[2][1] = basis_z[DirectionVector::Y];
  rotation_matrix_[2][2] = basis_z[DirectionVector::Z];

  // The Z basis vector gives the latitude (gamma_p) and longitude
  // (lambda_p) of the north pole in terms of the standard system.
  pole_latitude_ = atan2(rotation_matrix_[2][2],
    sqrt(pow(rotation_matrix_[2][0],2) + pow(rotation_matrix_[2][1],2)));
  pole_longitude_ = atan2(rotation_matrix_[2][1], rotation_matrix_[2][0]);

  // The X basis vector gives the reference latitude and longitude.
  reference_latitude_ = atan2(rotation_matrix_[0][2],
    sqrt(pow(rotation_matrix_[0][0],2) + pow(rotation_matrix_[0][1],2)));
  reference_longitude_ = atan2(rotation_matrix_[0][1], rotation_matrix_[0][0]);

  // Determine the pole rotation theta_p.  Use the equation for
  // tan(lambda_a + theta_p).  Set (gamma, lambda) to the reference
  // latitude and longitude; this point lies on the X basis vector
  // and thus lambda_a = 0.
  pole_rotation_ = atan2(cos(reference_latitude_)*sin(reference_longitude_ -
    pole_longitude_), sin(pole_latitude_)*cos(reference_latitude_)*
    cos(reference_longitude_ - pole_longitude_) -
    cos(pole_latitude_)*sin(reference_latitude_));

  is_trivial = false;
}

//------------------------------------------------------------------
// initBurstTransformation()
//
// Right now this does nothing
// But eventually it should be used to speed up Standard <> Oblique Cyl
// transformations within a single burst
//------------------------------------------------------------------
void 
OblCylProj::initBurstTransformation(){}

//------------------------------------------------------------------
// standardLatInRad()
//
// Given a latitude and longitude in the oblique system, return the
// corresponding latitude in the standard system.
//
// Use the equation for sin(gamma) from Section 2.6.2 of the BIDR SIS.
// gamma is the value being calculated, slatrad; gamma_a and lambda_a
// are the input parameters latrad and lonrad, respectively.
//------------------------------------------------------------------
double 
OblCylProj::standardLatInRad(double lonrad, double latrad) const{
  double slatrad;

  if(is_trivial){
    slatrad=latrad;
  }
  else{
    // asin() returns values between -pi/2 and pi/2
    slatrad = asin( sin(pole_latitude_)*sin(latrad) -
      cos(pole_latitude_)*cos(latrad)*cos(lonrad + pole_rotation_) );
  }
  return(slatrad);
}

//------------------------------------------------------------------
// standardLonInRad()
//
// Given a latitude and longitude in the oblique system, return the
// corresponding longitude in the standard system.
//
// Use the equation for tan(lambda - lambda_p) from Section 2.6.2 of the
// BIDR SIS.  lambda is the value being calculated, slonrad; gamma_a and
// lambda_a are the input parameters latrad and lonrad, respectively.
//
// The input parameters and return value are in radians.
//------------------------------------------------------------------
double 
OblCylProj::standardLonInRad(double lonrad, double latrad) const{
  double slonrad;

  if(is_trivial){
    slonrad=lonrad;
  }
  else{
    slonrad = pole_longitude_ +
      atan2(cos(latrad)*sin(lonrad + pole_rotation_),
      sin(pole_latitude_)*cos(latrad)*cos(lonrad + pole_rotation_) +
      cos(pole_latitude_)*sin(latrad));
  }
  return(boundLongitude(slonrad));
}

//------------------------------------------------------------------
// latInRad()
//
// Given a latitude and longitude in the standard system, return the
// corresponding latitude in the oblique system.
//
// The input parameters and return value are in radians.
//------------------------------------------------------------------
double 
OblCylProj::latInRad(double slonrad, double slatrad) const {
  double latrad;
  double s[3], obl[3];

  if(is_trivial){
    latrad=slatrad;
  }
  else{
    // Construct the position vector equivalent to (slatrad, slonrad)
    s[0] = cos(slatrad)*cos(slonrad);
    s[1] = cos(slatrad)*sin(slonrad);
    s[2] = sin(slatrad);
    rotate_vector(s, const_cast<double(*)[3]>(rotation_matrix_), obl);
    latrad = atan2(obl[2], sqrt(pow(obl[0],2) + pow(obl[1],2)));
  }
  return(latrad);
}

//------------------------------------------------------------------
// lonInRad()
//
// Given a latitude and longitude in the standard system, return the
// corresponding longitude in the oblique system.
//
// The input parameters and return value are in radians.
//------------------------------------------------------------------
double 
OblCylProj::lonInRad(double slonrad, double slatrad) const {
  double lonrad;
  double s[3], obl[3];

  if(is_trivial){
    lonrad=slonrad;
  }
  else{
    // Construct the position vector equivalent to (slatrad, slonrad)
    s[0] = cos(slatrad)*cos(slonrad);
    s[1] = cos(slatrad)*sin(slonrad);
    s[2] = sin(slatrad);
    rotate_vector(s, const_cast<double(*)[3]>(rotation_matrix_), obl);
    lonrad = atan2(obl[1], obl[0]);
  }
  return(boundLongitude(lonrad));
}

//------------------------------------------------------------------
// rotationMatrix()
//
// Return the 3x3 matrix for the transformation from standard to
// oblique coordinates.  Its rows are OBLIQUE_PROJ_X_AXIS_VECTOR,
// OBLIQUE_PROJ_Y_AXIS_VECTOR, and OBLIQUE_PROJ_Z_AXIS_VECTOR.
//------------------------------------------------------------------
void
OblCylProj::rotationMatrix(double m[3][3]){
  m[0][0] = rotation_matrix_[0][0];
  m[0][1] = rotation_matrix_[0][1];
  m[0][2] = rotation_matrix_[0][2];
  m[1][0] = rotation_matrix_[1][0];
  m[1][1] = rotation_matrix_[1][1];
  m[1][2] = rotation_matrix_[1][2];
  m[2][0] = rotation_matrix_[2][0];
  m[2][1] = rotation_matrix_[2][1];
  m[2][2] = rotation_matrix_[2][2];
}

//------------------------------------------------------------------
// poleLatitude()
//
// Return the value of OBLIQUE_PROJ_POLE_LATITUDE in degrees.
//------------------------------------------------------------------
double
OblCylProj::poleLatitude(){
  return(pole_latitude_ * radtodeg);
}

//------------------------------------------------------------------
// poleLongitude()
//
// Return the value of OBLIQUE_PROJ_POLE_LONGITUDE as a positive-west
// quantity in degrees.
//------------------------------------------------------------------
double
OblCylProj::poleLongitude(){
  if (pole_longitude_ == 0.0) {
    return(0.0);
  }
  return(360.0 - boundLongitude(pole_longitude_)*radtodeg);
}

//------------------------------------------------------------------
// poleRotation()
//
// Return the value of OBLIQUE_PROJ_POLE_ROTATION in degrees.
// Use a range of 0 to 360 degrees.
//------------------------------------------------------------------
double
OblCylProj::poleRotation(){
  return(boundLongitude(pole_rotation_)*radtodeg);
}

//------------------------------------------------------------------
// referenceLatitude()
//
// Return the value of REFERENCE_LATITUDE in degrees.
//------------------------------------------------------------------
double
OblCylProj::referenceLatitude(){
  return(radtodeg * reference_latitude_);
}

//------------------------------------------------------------------
// referenceLongitude()
//
// Return the value of REFERENCE_LONGITUDE as a positive-west
// quantity in degrees.
//------------------------------------------------------------------
double
OblCylProj::referenceLongitude(){
  if (reference_longitude_ == 0.0) {
    return(0.0);
  }
  return(360.0 - boundLongitude(reference_longitude_)*radtodeg);
}

//------------------------------------------------------------------
// boundLongitude()
//
// Recalculate angle d as the equivalent angle in the range
// [0, 2*pi) radians.
//------------------------------------------------------------------
double
OblCylProj::boundLongitude(double d) const {
  while (d < 0) {
    d += 2*pi;
  }
  while (d >= 2*pi) {
    d -= 2*pi;
  }
  return(d);
}
