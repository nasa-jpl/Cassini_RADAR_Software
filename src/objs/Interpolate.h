//==========================================================//
// Copyright (C) 1997, California Institute of Technology.	//
// U.S. Government sponsorship acknowledged.				//
//==========================================================//

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <stdlib.h>

static const char rcs_id_interpolate_h[] =
	"@(#) $Id: Interpolate.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

//======================================================================
// DESCRIPTION
//		Interpolation functions
//		Polynomial fitting
//======================================================================

// for polint, n is the order+1
int polint(double xa[], double ya[], int n, double x, double* y);

// for polcoe, n is the order
int polcoe( double x[], double y[], int	n, double* cof);

// for polyval, N is the order
float polyval(float x, float a[], int N);

int  find_target(double x[2], double y[2], float z[2], double target_z,
         double* target_x, double* target_y);

int  get_quad_peak(double c[3], double* peak_location, double* peak_value);

int cubic_spline(double* x, double* y, int n, double yp1, double ypn, double* y2);
int interpolate_cubic_spline(double* xa, double* ya, double* y2a, int n,
			     double x, double* y);

int interpolate_cubic_spline(double x0, double xstep, 
			     double* ya, double* y2a, int n,
			     double x, double* y, int start_idx, 
			     int skip);

#endif



