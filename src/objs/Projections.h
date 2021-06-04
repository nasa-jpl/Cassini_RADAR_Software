#ifndef PROJECTIONS_H
#define PROJECTIONS_H
#include "Frame.h"

static const char rcs_id_projections_h[] =
  "@(#) $Id: Projections.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

class OblCylProj{

 public:
  OblCylProj(); // set up trivial (standard) projection

  // Set up an Oblique cylindrical projection
  // in which the subspacecraft position in s denotes the intersection of
  // the equation and the prime meridian, and the spacecraft velocity 
  // direction defines the equator. s is presumed to be the state of Cassini
  // in the Target Body Fixed frame at the time of closest approach.
  OblCylProj(StateVector& s);

  void initBurstTransformation(); // INCOMPLETE
  double standardLatInRad(double lonrad, double latrad) const;
  double standardLonInRad(double lonrad, double latrad) const;
  double latInRad(double slonrad, double slatrad) const;
  double lonInRad(double slonrad, double slatrad) const;
  void rotationMatrix(double m[3][3]);
  double poleLatitude();
  double poleLongitude();
  double poleRotation();
  double referenceLatitude();
  double referenceLongitude();

 private:
  bool is_trivial; // In this case the projection is the standard lat/lon
  double rotation_matrix_[3][3];
  double pole_latitude_;  // gamma_p in BIDR SIS
  double pole_longitude_; // lambda_p in BIDR SIS; stored internally in
                          // positive-east coordinates
  double pole_rotation_;  // theta_p in BIDR SIS
  double reference_latitude_;
  double reference_longitude_; // stored internally in positive-east coordinates
  double boundLongitude(double d) const;
};

#endif
