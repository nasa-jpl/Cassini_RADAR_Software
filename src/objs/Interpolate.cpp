//==============================================================//
// Copyright (C) 1997-1998, California Institute of Technology. //
// U.S. Government sponsorship acknowledged.                    //
//==============================================================//

static const char rcs_id_interpolate_c[] =
    "@(#) $Id: Interpolate.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "Interpolate.h"

//--------//
// polint //
//--------//

int
polint(
	double		xa[],
	double		ya[],
	int			n,
	double		x,
	double*		y)
{
	//-----------------------------------//
	// allocate work arrays if necessary //
	//-----------------------------------//

	static double* c = NULL;
	static double* d = NULL;
	static int last_n = 0;

	if (n != last_n)
	{
		free(c);
		c = (double *)malloc(n * sizeof(double));
		if (c == NULL)
			return(0);

		free(d);
		d = (double *)malloc(n * sizeof(double));
		if (d == NULL)
			return(0);

		last_n = n;
		if (n == 0)
			return(1);
	}

	//-------------------------------------------//
	// find the index of the closest table entry //
	//-------------------------------------------//

	int ns = 0;
	double dif = fabs(x - xa[0]);
	for (int i = 0; i < n; i++)
	{
		double dift = fabs(x - xa[i]);
		if (dift < dif)
		{
			ns = i;
			dif = dift;
		}

		// initialize
		c[i] = ya[i];
		d[i] = ya[i];
	}

	double value = ya[ns--];
	for (int m = 1; m < n; m++)
	{
		for (int i = 0; i < n - m; i++)
		{
			double ho = xa[i] - x;
			double hp = xa[i+m] - x;
			double w = c[i+1] - d[i];
			double den = ho - hp;
			if (den == 0.0)
			{
				free(c);
				free(d);
				return(0);		// error
			}
			den = w / den;
			d[i] = hp * den;
			c[i] = ho * den;
		}

		double adj;
		if (2 * ns < (n - m - 1))
			adj = c[ns+1];
		else
			adj = d[ns--];
		value += adj;
	}

	*y = value;
	return(1);
}

//--------//
// polcoe //
//--------//

int
polcoe(
	double		x[],
	double		y[],
	int			n,
	double*		cof)
{
	//-----------------------------------//
	// allocate work arrays if necessary //
	//-----------------------------------//

	static double* s = NULL;
	static int last_n = 0;

	if (n != last_n)
	{
		free(s);
		s = (double *)malloc((n+1) * sizeof(double));
		if (s == NULL)
			return(0);

		last_n = n;
	}

	int i,j,k;
	double phi,ff,b;

	for (i=0; i <= n; i++)
	{
		s[i] = 0.0;
		cof[i] = 0.0;
	}

	s[n] = -x[0];
	for (i=1; i <= n; i++)
	{
		for (j=n-i; j <= n-1; j++)
		{
			s[j] -= x[i]*s[j+1];
		}
		s[n] -= x[i];
	}

	for (j=0; j <= n; j++)
	{
		phi = n+1;
		for (k=n; k >= 1; k--)
		{
			phi = k*s[k] + x[j]*phi;
		}
		ff = y[j]/phi;
		b = 1.0;
		for (k=n; k >= 0; k--)
		{
			cof[k] += b*ff;
			b = s[k] + x[j]*b;
		}
	}

	return(1);
}

//--------------------------------------------------------//
// Evaluate a polynomial function.
//
// y = a[N]*x^N + a[N-1]*x^(N-1) + ... + a[1]*x + a[0].
//
// y is the return value,
// a must be a vector of length N+1.
//--------------------------------------------------------//
 
float
polyval(
    float  x,
    float  a[],
    int    N)
{
    float value;
 
    value = a[N];
    for (int i= N-1; i >= 0; i--)
    {
        value = value*x + a[i];
    }
    return(value);
}

//-------------//
// find_target //
//-------------//
// Takes two (x, y, z) pairs and a z target value.  Returns the values
// of x and y that would correspond to the target z value.  Assumes
// linear relationship among all axes.

int
find_target(
    double   x[2],
    double   y[2],
    float    z[2],
    double   target_z,
    double*  target_x,
    double*  target_y)
{
    double delta_x = x[1] - x[0];
    double delta_y = y[1] - y[0];
    double delta_z = z[1] - z[0];

    double z_frac = (target_z - z[0]) / delta_z;
    *target_x = x[0] + z_frac * delta_x;
    *target_y = y[0] + z_frac * delta_y;
    return(1);
}

//---------------//
// get_quad_peak //
//---------------//
// finds the peak of a quadratic

int
get_quad_peak(
    double   c[3],
    double*  peak_location,
    double*  peak_value)
{
    if (c[2] == 0.0)
        return(0);
    *peak_location = -c[1] / (2.0 * c[2]);
    *peak_value = ((c[2] * (*peak_location)) + c[1]) * (*peak_location) + c[0];
    return(1);
}

//-------------------------------------------//
// Function cubic_spline                     //
// calculate 2nd derivatives                 //
// for use in cubic spline interpolation     //
//-------------------------------------------//

int cubic_spline(double* x, double* y, int n, double yp1, double ypn, 
		 double* y2){
  // CUBIC SPLINE INTERPOLATION FOR FUNCTION Y(X) GIVEN VALUES FOR THE
  //  DERIVATIVES YP1,YPN AT POINTS X1, XN.  RETURNS ARRAY Y2 OF VALUES FOR THE
  //  SECOND DERIVATIVE AT POINTS X_I.  "NATURAL CUBIC SPLINE" IS WHEN
  //  Y2(1)=Y2(N)=0, WHICH IS ACCOMPLISHED BY INPUTTING YP1,YPN > 1.d30.
  //  [CAN HAVE DIFFERENT BOUNDARY CONDITIONS ON DIFFERENT ENDS]

      double* u=new double[n];
      if(yp1 > 0.99e+30){
        y2[0]=0;
        u[0]=0;
      }
      else{
        y2[0]=-0.5;
        u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
      }
      for(int i=1;i< n-1; i++){
        double sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        double p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])
		  /(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;
      }
      double qn;
      if(ypn > 0.99e+30){
        qn=0.0;
        u[n-1]=0.0;
      }
      else{
        qn=0.5;
        u[n-1]=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
      }
      y2[n-1]=(u[n-1]-qn*u[n-2])/(qn*y2[n-2]+1.0);
      for(int k=n-2; k >= 0 ; k--){
        y2[k]=y2[k]*y2[k+1]+u[k];
      }

      delete[] u;
      return(1);
}
//---------------------------------------------//
// Function interpolate_cubic_spline           //
// Interpolate using a cubic spline            //
//---------------------------------------------//
int interpolate_cubic_spline(double* xa, double* ya, double* y2a, int n,
			     double x, double* y){
  //     DOES THE INTERPOLATION.
  //     XA AND YA TABULATE THE FUNCTION, Y2A IS THE OUTPUT OF SPLINE (AN
  //     ARRAY OF SECOND DERIVATIVE VALUES, X IS VALUE AT WHICH Y IS SOUGHT.
     
  int    klo=0;
  int    khi=n-1;
  while(khi-klo > 1){
    int k=(khi+klo)/2;
    if(xa[k]>x)
      khi=k;
    else
      klo=k;
  }
  double h=xa[khi]-xa[klo];
  if(h==0.0){
     fprintf(stderr,"interpolate_cubic_spline:  the xas must be distinct\n");
     exit(0);
  }
  double a=(xa[khi]-x)/h;
  double b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6;
  return(1);
}


//---------------------------------------------//
// Function interpolate_cubic_spline           //
// Interpolate using a cubic spline            //
//---------------------------------------------//
int interpolate_cubic_spline(double x0, double xstep, double* ya, 
			     double* y2a, int n, double x, double* y, 
			     int start_idx, int skip){
  //     DOES THE INTERPOLATION.
  //     XA AND YA TABULATE THE FUNCTION, Y2A IS THE OUTPUT OF SPLINE (AN
  //     ARRAY OF SECOND DERIVATIVE VALUES, X IS VALUE AT WHICH Y IS SOUGHT.
     
  double i=(x-x0)/xstep;
  int    klo=int(i);
  int    khi=klo+1;
  if(klo<0 || khi>n-1) return(0);
  double xlo=x0+klo*xstep;
  double xhi=xlo+xstep;


  double h=xstep;
  if(h==0.0){
     fprintf(stderr,"interpolate_cubic_spline:  the xas must be distinct\n");
     exit(0);
  }
  double a=(xhi-x)/h;
  double b=(x-xlo)/h;

  // apply skip;
  int jlo=start_idx+klo*skip;
  int jhi=jlo+skip;
  *y=a*ya[jlo]+b*ya[jhi]+((a*a*a-a)*y2a[jlo]+(b*b*b-b)*y2a[jhi])*(h*h)/6;
  return(1);
}

