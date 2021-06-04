//=============================================================================
// TemplateUtils.h
//
// This file contains utility template functions.
// Because these functions are templates, and gcc doesn't support the
// export keyword, the actual function definitions are in this file.
//
//=============================================================================

#ifndef TemplateUtils_H
#define TemplateUtils_H

#include <stdio.h>
#include <string>
#include <sstream>
#include "Error.h"
#include "Units.h"
#include "Frame.h"
#include "Array.h"
#include "Time.h"
#include "Config.h"

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;

//---------------------------------------------------------------------------
// str = toStr(arg)
//
// Convert the argument to a string using standard stream i/o.
//---------------------------------------------------------------------------

/** temporarily in Utils.h until I can change all the dependent modules...
template<class T>
string toStr(const T& arg)
  {    
  ostringstream ost;
  ost << arg;
  return(string(ost.str(),ost.pcount()));
  } 
**/

//-------------------------------------------------------------------------
// bracket_minimum(x1,x2,x3,funk)
//
// Bracket a minimum of funk between three points [x1,x2,x3] given an
// initial x1 and middle to search from in the downhill direction.
// funk is a function object.
//-------------------------------------------------------------------------

template <class F>
void bracket_minimum(Uvar& x1, Uvar& x2, Uvar& x3, F funk)
  throw(ErrorMessage)
  {
  Uvar f1 = funk(x1);
  Uvar f2 = funk(x2);
  if (f1 < f2)
    {
    swap(x1,x2);
    swap(f1,f2);
    }

  Uvar delta = x2 - x1;

  //--------------------------------------------------------------------------- 
  // Search for a 3rd point whose function value is more than f2 by successively
  // stepping away from x2 with jumps that grow by a factor of 2 with each
  // try.  Up to 100 tries are allowed.
  //--------------------------------------------------------------------------- 

  for (unsigned int i=0; i < 100; ++i)
    {
    x3 = x2 + delta;
    Uvar f3 = funk(x3);
    if (f3 > f2)
      {  // found a point higher so that x1,x2,x3 form a bracketing triple
      return;
      }
    delta *= 2;
    // Move x1,x2 downhill to follow x3 as it sweeps for a higher point
    // This makes a following search more efficient, and helps to prevent
    // an unusually large step to get f3 above the starting f2.
    x1 = x2;
    f1 = f2;
    x2 = x3;
    f2 = f3;
    }

  // If we made it to here, then no triple was found in the preceeding loop.
  throw
    ErrorMessage("bracket_minimum: Unable to bracket minimum within 100 steps");
  }

//---------------------------------------------------------------------
// golden_section_search(side1,middle,side2,x_tol,funk,final_x,final_y) 
//
// A general purpose, one-dimensional golden section search which
// finds the minimum of a function object (funk) bracketed in the interval
// [side1,middle,side2]
//---------------------------------------------------------------------

template <class F>
void golden_section_search(
    const Uvar& side1,
    const Uvar& middle,
    const Uvar& side2,
    const Uvar& x_tol,
    F   funk,
    Uvar&  final_x,
    Uvar&  final_y)
  {
  static const double golden_c = (3.0 - sqrt(5.0)) / 2.0;
//  static const double golden_r = 1.0 - golden_c;

  Uvar x0,x1,x2,x3;

  //--------------------------//
  // initialize search values //
  //--------------------------//

  x0 = side1;
  x3 = side2;
  if (fabs(side2 - side1) > fabs(middle - side1))
    {
    x1 = middle;
    x2 = middle + golden_c*(side2 - middle);
    }
  else
    {
    x2 = middle;
    x1 = middle - golden_c*(middle - side1);
    }

  Uvar f1 = funk(x1);
  Uvar f2 = funk(x2);

  while (fabs(x3 - x0) > x_tol)
    {
    if (f2 < f1)
      {
      x0 = x1;
      x1 = x2;
      x2 = x2 + golden_c * (x3 - x2);
      f1 = f2;
      f2 = funk(x2);
      }
    else
      {
      x3 = x2;
      x2 = x1;
      x1 = x1 - golden_c * (x1 - x0);
      f2 = f1;
      f1 = funk(x1);
      }
    }

    if (f1 < f2)
      {
      final_x = x1;
      final_y = f1;
      }
    else
      {
      final_x = x2;
      final_y = f2;
      }
  }

//-------------------------------------------------------------------------
// downhill_simplex(Dmat& p,Dvec& y,
//                  const unsigned int& ndim,
//                  const double& ftol,
//                  F funk, unsigned int& nfunk)
//
// miltidimensional minimization of the function funk(x) where x[1..ndim]
// is a vector in dim multidimensions.
//(see Numerical Recipes in C, Chapter 10 section 4)
//
// p: (ndim+1) * ndim matrix input.  Its ndim+1 rows are ndim-dimensional
//          vectors which are the vertices of the starting simplex
// y: (ndim+1) functional values at (ndim+1) vertices
// ndim: parameter dimension
// ftol: functional convergence tolerance to be achieved 
// On output, p and y will have been reset to ndim+1 new points all within
// ftol of a minimum function value, and nfunk gives the number of function
// evalulations taken.
//-------------------------------------------------------------------------

template <class F>
void downhill_simplex(Dmat& p,
		      Dvec& y, 
		      const unsigned int& ndim,
		      const double& ftol,
		      F funk, 
		      unsigned int& nfunk)
  {
    if (p.rows() != ndim+1) ErrorMessage("TemplateUtils.h::downhill_simples:p should be  a (ndim+1) X (ndim) matrix").throwMe();
    if (p.cols() != ndim) ErrorMessage("TemplateUtils.h::downhill_simples:p should be  a (ndim+1) X (ndim) matrix").throwMe();
    if (y.size() != ndim+1)  ErrorMessage("TemplateUtils.h::downhill_simples:y should have a size of ndim").throwMe();
    

    unsigned int mpoints = ndim+1;
    unsigned int ilo, ihi, inhi;
    double rtol,ysave,ytry;
    unsigned int Nmax = 5000;
    Dvec psum("psum",p.cols());

    //get psum
    for (unsigned int i = 0 ; i < p.cols();++i)
      {
	double sum = 0.0;
	for (unsigned int j = 0; j < p.rows();++j)  sum += p(j,i);
	psum(i) = sum;
      }

    nfunk = 0;

    for (unsigned int k = 0; k <Nmax;++k  )
      {
	//loop until k reaches Nmax
	ilo = 0;
	
	//determine which point is the highest, next highest, and lowest
	if (y(0) > y(1)) 
	  { 
	    ihi = 0; 
	    inhi = 1;
	  }
	else 
	  {
	    ihi = 1; 
	    inhi = 0;
	  }
	for (unsigned int i=0; i <mpoints;++i)
	  {
	    if (y(i) < y(ilo)) ilo = i;
	    if (y(i) > y(ihi))
	      {
		ihi = i;
		inhi = ihi;
	      }
	    else if (y(i) > y(inhi) && i != ihi) inhi = i;
	  }
	rtol = 2.0 * fabs(y(ihi)-y(ilo))/(fabs(y(ihi))+fabs(y(ilo)));

	//cout<<"ihi inhi ilo "<<ihi<<" "<<inhi<<" "<<ilo<<endl;
	//cout<<"y "<<y<<endl;
	//cout<<"nfunk "<<nfunk<<endl;
	
	if (rtol < ftol)
	  {
	    //break after save
	    double swap = y(0);
	    y(0) = y(ilo);
	    y(ilo) = swap;
	    for (unsigned int i = 0; i < ndim;++i)
	      {
		swap = p(0,i);
		p(0,i) = p(ilo,i);
		p(ilo,i) = swap;
	      }
	    break;
	  }
	if (nfunk > Nmax) ErrorMessage("Nmax exceeded").throwMe();
	nfunk +=2;

	//cout<<"reflection "<<endl;
	//reflect the simples from the high point
	ytry = amotry(p,y,psum,ndim,funk,ihi,-1.0);
	
	if (ytry <= y(ilo))
	  {
	    //cout<<"expansion "<<endl;
	    ytry = amotry(p,y,psum,ndim,funk,ihi,2.0);
	  }
	else if (ytry >=y(inhi))
	  {//worse reflection, do contraction
	    //cout<<"contraction "<<endl;
	    ysave = y(ihi);
	    ytry = amotry(p,y,psum,ndim,funk,ihi,0.5);
	    if (ytry >= ysave)
	      {//Can't seem to get rid of that high point, better contract around the lowest point
		for (unsigned int i =1; i<mpoints;++i)
		  {
		    if (i != ilo)
		      {
			for (unsigned j = 1; j < ndim;++j)
			  {
			    psum(j)= 0.5*(p(i,j)+p(ilo,j));
			    p(i,j) = psum(j);
			  }
			y(i) = funk(psum);
		      }
		  }
		nfunk += ndim;
		//sum over p
		for (unsigned int i = 0 ; i < p.cols();++i)
		  {
		    double sum = 0.0;
		    for (unsigned int j = 0; j < p.rows();++j)  sum += p(j,i);
		    psum(i) = sum;
		  }
	      }
	    else 
	      {
		--nfunk;
	      }
	  }
      }//begin new iteration
  }



template <class F>
double amotry(Dmat& p, Dvec& y, Dvec& psum, const unsigned int& ndim,
	    F funk, const unsigned int& ihi, const double& fac)
  {
    double fac1,fac2,ytry;
    Dvec ptry("ptry",p.cols());
    fac1 = (1.0 - fac)/double(ndim);
    fac2 = fac1 - fac;
    for (unsigned int j = 0; j < ndim;++j) 
      {
	ptry(j) = psum(j)*fac1-p(ihi,j)*fac2;
      }
    ytry = funk(ptry);
    if (ytry < y(ihi))
      {
	y(ihi) = ytry;
	for (unsigned int j = 0;j < ndim;++j)
	  {
	    psum(j)+= ptry(j) - p(ihi,j);
	    p(ihi,j) = ptry(j);
	  }
      }
    return(ytry);
  }
#endif 


