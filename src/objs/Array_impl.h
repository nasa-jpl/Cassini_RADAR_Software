//----------------------------------------------------------------------------
// Array_impl.h
//
// This file contains method definitions for the Array handling classes
// These classes use the NAIF Spice toolkit to provide automatic handling
// of position and state vectors and the coordinate transformations between
// them.
// A .h file is used here instead of .cpp because template classes need their
// method definitions in the same .h file (this file is included by Array.h).
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_array_impl_h[] =
  "@(#) $Id: Array_impl.h,v 11.5 2011/09/16 00:03:29 richw Exp $";

#include <string>
#include <complex>
#include <iostream>
#include <math.h>
#include "Array.h"
#include "Utils.h"
#include "Units.h"
//string toStr(ostringstream& arg);

//-------------------------------------------------------
// Methods for classes Array1D, Array2D
//-------------------------------------------------------

//--------------------------
// Methods for class Array1D
//--------------------------

//--------------
// Constructors
//--------------

template<class T>
Array1D<T>::Array1D(const string& name, unsigned int length = 1)
  : name_(name), v_(length)
  {
  if (length == 0)
    {
    ArrayError e("Can't make array: " + name + " with zero length");
    e.throwMe();
    }
  }

// copy constructor
template<class T>
Array1D<T>::Array1D(const Array1D<T>& a)
  : name_(a.name_), v_(a.v_)
  {
  for (unsigned int i=0; i < v_.size(); ++i)
    {
    v_[i] = a.v_[i];
    }
  }

// destructor
template<class T>
Array1D<T>::~Array1D()
  {
  }

//-----------
// Operators
//-----------

// unary -
template<class T>
Array1D<T> 
  Array1D<T>::operator-() const
  {
    Array1D<T> retval=*this;
    retval.negate();
    return(retval);
  }

// +=
template<class T>
Array1D<T>&
  Array1D<T>::operator+=(const Array1D<T>& arg)
  {
  if (arg.size() != size())
    {
    ArrayError e("Size mismatch adding arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] += arg.v_[i];
    }
  return(*this);
  }

// +=
template<class T>
Array1D<T>&
  Array1D<T>::operator+=(const T& arg)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] += arg;
    }
  return(*this);
  }

// -=
template<class T>
Array1D<T>&
  Array1D<T>::operator-=(const Array1D<T>& arg)
  {
  if (arg.size() != size())
    {
    ArrayError e("Size mismatch subtracting arrays " << name() <<
      " and " << arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] -= arg.v_[i];
    }
  return(*this);
  }

// -=
template<class T>
Array1D<T>&
  Array1D<T>::operator-=(const T& arg)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] -= arg;
    }
  return(*this);
  }

// *=
template<class T>
Array1D<T>&
  Array1D<T>::operator*=(const Array1D<T>& arg)
  {
  if (arg.size() != size())
    {
    ArrayError e("Size mismatch multiplying arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] *= arg.v_[i];
    }
  return(*this);
  }

// *=
template<class T>
Array1D<T>&
  Array1D<T>::operator*=(const T& arg)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] *= arg;
    }
  return(*this);
  }

// /=
template<class T>
Array1D<T>&
  Array1D<T>::operator/=(const Array1D<T>& arg)
  {
  if (arg.size() != size())
    {
    ArrayError e("Size mismatch dividing arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] /= arg.v_[i];
    }
  return(*this);
  }

// /=
template<class T>
Array1D<T>&
  Array1D<T>::operator/=(const T& arg)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] /= arg;
    }
  return(*this);
  }

// Subscripting
template<class T>
const T& Array1D<T>::operator()(unsigned int i) const
  {
  if (i >= size())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 1D array: " << name() << " with index: "
      << i;
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(v_[i]);
  }

template<class T>
T& Array1D<T>::operator()(unsigned int i)
  {
  if (i >= size())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 1D array: " << name() << " with index: "
      << i;
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(v_[i]);
  }

// Named method Subscripting - helps with gdb debugging
template<class T>
const T& Array1D<T>::elem(unsigned int i)
  const
  {
  if (i >= size())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 1D array: " << name() << " with index: "
      << i;
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(v_[i]);
  }

// Assignment
/**
template<class T>
template<class T2>
Array1D<T>& Array1D<T>::operator=(const Array1D<T2>& a) // No exceptions
  {
  resize(a.size());
  for (unsigned int i = 0; i < a.size(); ++i)
    {
    v_[i] = a(i);  // rely on assignment operator between T2 and T.
    }
  return(*this);
  }
**/

template<class T>
Array1D<T>& Array1D<T>::operator=(const Array1D& a) // No exceptions
  {
  if (v_ == a.v_) return *this;  // cover self assignment
  v_ = a.v_;  // transfer data, but keep name the same.
  return(*this);
  }

template<class T>
Array1D<T>& Array1D<T>::operator=(const T& a) // No exceptions
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] = a;
    }
  return(*this);
  }

//-------------------------
// I/O
//-------------------------

template<class T>
void Array1D<T>::show() const
  {
  std::cout << name() << "(" << size() << ")" << std::endl;
  int count = 1;
  for (unsigned int i = 0; i < size(); ++i)
    {
    if (i > 0) std::cout << ", ";
    std::cout << v_[i];
    count++;
    if (count >= 6)
      {
      count = 1;
      std::cout << std::endl;
      }
    }
  if (count != 1) std::cout << std::endl;
  }

template<class T>
void Array1D<T>::show(unsigned int i) const
  {
  std::cout << name() << "(" << i << "/" << size() << ") = " << elem(i) << std::endl;
  }

//-------------------------
// Other arithmetic methods
//-------------------------

//---------------------
// negate()
//
// Negate each element
//---------------------

template<class T>
Array1D<T>&
  Array1D<T>::negate()
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] = -v_[i];
    }
  return(*this);
  }

//-----------------------------------------
// pow(b)
//
// *this(i) = *this(i) ^ b
//-----------------------------------------

template<class T>
Array1D<T>&
  Array1D<T>::pow(const T& b)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] = pow(v_[i],b);
    }
  return(*this);
  }

//-----------------------------------------
// dft(b)
//-----------------------------------------
/*Array1D<double> dft(const Array1D<double>& b)
   {
      return(b);
   }
*/
//-----------------------------------------
// pow(b)
//
// *this(i) = *this(i) ^ b(i)
//-----------------------------------------

template<class T>
Array1D<T>&
  Array1D<T>::pow(const Array1D<T>& b)
  {
  if (b.size() != size())
    {
    ArrayError e("Size mismatch in " + name() + ".pow(" + b.name() + ")");
    e.throwMe();
    }
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] = pow(v_[i],b);
    }
  return(*this);
  }

//-----------------------------------------
// inverse_pow(base)
//
// Raise base to the power of each element
// ie., *this(i) = base ^ (*this(i))
//-----------------------------------------

template<class T>
Array1D<T>&
  Array1D<T>::inverse_pow(const T& base)
  {
  for (unsigned int i = 0; i < size(); ++i)
    {
    v_[i] = pow(base,v_[i]);
    }
  return(*this);
  }

//---------------
// Other methods
//---------------

//-------------------------------------------------------
// magnitude()
//
// Return the magnitude of *this vector (including units)
//-------------------------------------------------------

template<class T>
T Array1D<T>::magnitude() const 
  {
  T result = v_[0] * v_[0];
  for (unsigned int i=1; i < v_.size(); ++i)
    {
    result += v_[i] * v_[i];
    }
  return(sqrt(result));
  }

//-------------------------------------------------------
// min(unsigned int& idx)
//
// Return the minimum element of *this vector (including units)
// Store the index of the minimum element in idx
//-------------------------------------------------------

template<class T>
T Array1D<T>::min(unsigned int& idx) const 
  {
  T result = v_[0];
  idx=0;
  for (unsigned int i=1; i < v_.size(); ++i)
    {
    if(v_[i]<result)
      {
      result=v_[i];
      idx=i;
      }
    }
  return(result);
  }

//-------------------------------------------------------
// max(unsigned int& idx)
//
// Return the maximum element of *this vector (including units)
// Store the index of the maximum element in idx
//-------------------------------------------------------
template<class T>
T Array1D<T>::max(unsigned int& idx) const 
  {
  T result = v_[0];
  idx=0;
  for (unsigned int i=1; i < v_.size(); ++i)
    {
    if(v_[i]>result)
      {
      result=v_[i];
      idx=i;
      }
    }
  return(result);
  }
//-----------------
//  Original code
//-----------------
//template<class T>
//T Array1D<T>::max(unsigned int& idx) const 
//  {
//  return((-(*this)).min(idx));
//  }

//-------------------------------------------------------
// sum()
//
// Return sum of *this vector (including units)
//-------------------------------------------------------

template<class T>
T Array1D<T>::sum() const 
  {
  T result = v_[0];
  for (unsigned int i=1; i < v_.size(); ++i)
    {
    result+=v_[i];
    }
  return(result);
  }

//-------------------------------------------------------
// mean()
//
// Return mean of *this vector (including units)
//-------------------------------------------------------

template<class T>
T Array1D<T>::mean() const 
  {
  T result = v_[0];
  for (unsigned int i=1; i < v_.size(); ++i)
    {
    result+=v_[i];
    }
  return(result/v_.size());
  }

//-------------------------------------------------------
// rms()
//
// Return RMS of *this vector (including units)
//-------------------------------------------------------

template<class T>
T Array1D<T>::rms(unsigned int n) const 
  {
  T result = v_[0]*v_[0];
  for (unsigned int i=1; i < n; ++i)
    {
    result+=v_[i]*v_[i];
    }
  result=result/n;
  result=sqrt(result);
  return(result);
  }

//-------------------------------------------------
// scale(mag)
//
// Scale *this to have the indicated magnitude.
// Throw an exception if attempting to scale zero.
//-------------------------------------------------

template<class T>
void Array1D<T>::scale(T mag) 
  {
  T cur_mag = magnitude();
  if (cur_mag == 0) { 
    ArrayError e(ArrayError::scale_zero);
    e.throwMe();
  }
  for (unsigned int i=0; i < v_.size(); ++i)
    {
    v_[i] *= mag/cur_mag;
    }
  }

//-------------------------------------------------
// fabs()
// Make all negative elements of this vector positive
//-------------------------------------------------

template<class T>
void Array1D<T>::fabs() 
  {
  for (unsigned int i=0; i < v_.size(); ++i)
    {
    v_[i] = ::fabs(v_[i]);
    }
  }

//-------------------------------------------------
// name()
//
// Return the name of this array.
//-------------------------------------------------

template<class T>
string Array1D<T>::name() const // No exceptions
  {
  return(name_);
  }

//-------------------------------------------------
// setName()
//
// Set the name of this array.
//-------------------------------------------------

template<class T>
void Array1D<T>::setName(const string& newname) // No exceptions
  {
  name_ = newname;
  }

//-------------------------------------------------
// size()
//
// Return the number of elements.
//-------------------------------------------------

template<class T>
unsigned int Array1D<T>::size() const // No exceptions
  {
  return(v_.size());
  }

//-------------------------------------------------
// resize()
//
// Set the number of elements.
//-------------------------------------------------

template<class T>
void Array1D<T>::resize(unsigned int newsize)
  {
  if (newsize == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero length");
    e.throwMe();
    }
  v_.resize(newsize);
  }

//-----------------------------
// findInterpolatedIndices
// 
// computes floating point indices that would have
// the target value using linear interpolation
//
// For example= if A = { 0,10,20,30,0}
// A.findInterpolatedIndices(15) = {1.5,3.5}
//-------------------------------

template<class T>
Array1D<double> Array1D<T>::findInterpolatedIndices(const T& y)
 const 
{
  Array1D<T> tmp=*this-y;

  // count sign changes
  bool positive = (tmp(0)>0);
  unsigned int num_changes=0;
  for(unsigned int c=1; c<size(); c++){
    if(positive!=(tmp(c)>0))
    {
      num_changes++;
      positive=tmp(c)>0;
    }
  }
  if(num_changes==0)
    {
      ArrayError e("FindInterpolatedIndices:: Target not found in array "
		       + name());
      e.throwMe();
    }
  Array1D<double> idx("idx",num_changes);

  // linearly interpolate to find indices
  positive = (tmp(0)>0);
  unsigned int change_idx=0;
  
  for(unsigned int c=1; c<size(); c++){
    if(positive!=(tmp(c)>0))
    {
      
      idx(change_idx)=c-1;
      T dxdy=1/(tmp(c)-tmp(c-1));
      T dy=0-tmp(c-1);
      T dx=dxdy*dy;
      idx(change_idx)+=convert_to_double(dx);
      positive=(tmp(c)>0);
      change_idx++;
    }
  }
  return(idx);
}
//------------------------------

//------------------------------
// Supporting template functions
//------------------------------

//-------------------------------------------------------------------------
// Some notes about template ambiguities are in order here.
// The + and - operations are not defined for cases like (int + Uvec)
// or (Dvec - Uvar) because these are likely to be errors anyways.
// If the user knows that the Uvar has no units, and an operation with
// a non Uvar is desired, then convert_to_double should be used.
// The * and / operations do support some mixed case operations, but not
// all combinatorial possibilites are listed here.  If a new combination
// arises, then a specialized routine should be added for it.
// The long int routines currently allow an amibiguity.  If this arises,
// then more specializations will have to be explicitely written instead
// of relying on the template (see below).
//-------------------------------------------------------------------------

// +
template<class T, class T2>
Array1D<T>
  operator+(Array1D<T> a, const T2& b)
  {
  a += b;
  a.setName("+");
  return(a);
  }

// +
template<class T>
Array1D<T>
  operator+(const T& b, Array1D<T> a)
  {
  a += b;
  a.setName("+");
  return(a);
  }

// -
template<class T, class T2>
Array1D<T>
  operator-(Array1D<T> a, const T2& b)
  {
  a -= b;
  a.setName("-");
  return(a);
  }

// -
template<class T>
Array1D<T>
  operator-(const T& b, Array1D<T> a)
  {
  a.negate();
  a += b;
  a.setName("-");
  return(a);
  }

// *
template<class T, class T2>
Array1D<T>
  operator*(Array1D<T> a, const T2& b)
  {
  a *= b;
  a.setName("*");
  return(a);
  }

// *
template<class T>
Array1D<T>
  operator*(const T& b, Array1D<T> a)
  {
  a *= b;
  a.setName("*");
  return(a);
  }

// *
/** problem with this...
Uvec
  operator*(const Uvar& b, Array1D<double> a)
  {
  Uvec result("*",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i)*b;
    }
  return(result);
  }
**/

// *
// This will cause ambiguity if T == Array1D<double>
// Fortunately, it will cover other cases like int*Array1D<T> because
// int will convert to double.
template<class T>
Array1D<T>
  operator*(const double& b, Array1D<T> a)
  {
  a *= b;
  a.setName("*");
  return(a);
  }

// /
template<class T, class T2>
Array1D<T>
  operator/(Array1D<T> a, const T2& b)
  {
  a /= b;
  a.setName("/");
  return(a);
  }

// /
template<class T>
Array1D<T>
  operator/(const T& a, Array1D<T> b)
  {
  for (unsigned int i = 0; i < b.size(); ++i)
    {
    b(i) = a/b(i);
    }
  b.setName("/");
  return(b);
  }

// /
/** problem with this...
Uvec
  operator/(const Uvar& a, Array1D<double> b)
  {
  Uvec result("/",b.size());
  for (unsigned int i = 0; i < b.size(); ++i)
    {
    result(i) = a/b(i);
    }
  return(result);
  }
**/

// /
// This will cause ambiguity if T == Array1D<double>
template<class T>
Array1D<T>
  operator/(const double& a, Array1D<T> b)
  {
  for (unsigned int i = 0; i < b.size(); ++i)
    {
    b(i) = a/b(i);
    }
  b.setName("/");
  return(b);
  }

// ==
template<class T, class T2>
Array1D<bool>
  operator==(const Array1D<T>& a, const T2& b)
  
  {
  Array1D<bool> result("==",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i) == b;
    }
  return(result);
  }

// ==
template<class T, class T2>
Array1D<bool>
  operator==(const T2& b, const Array1D<T>& a)
  
  {
  Array1D<bool> result("==",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i) == b;
    }
  return(result);
  }

// ==
template<class T>
Array1D<bool>
  operator==(const Array1D<T>& a, const UnitVar& b)
  
  {
  Array1D<bool> result("==",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i) == b;
    }
  return(result);
  }

// ==
template<class T>
Array1D<bool>
  operator==(const UnitVar& b, const Array1D<T>& a)
  
  {
  Array1D<bool> result("==",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i) == b;
    }
  return(result);
  }

// ==
template<class T>
Array1D<bool>
  operator==(const Array1D<T>& a, const Array1D<T>& b)
  {
  if (a.size() != b.size())
    {
    ArrayError e("Size mismatch comparing(==) arrays " + a.name() +
      " and " + b.name());
    e.throwMe();
    }
  Array1D<bool> result("==",a.size());
  for (unsigned int i=0; i < a.size(); ++i)
    {
    result(i) = a(i) == b(i);
    }
  return(result);
  }

// <<
template<class T>
std::ostream& operator<<(std::ostream& s, const Array1D<T>& a)
  {
  for (unsigned int i=0; i < a.size(); ++i)
    {
    if (i != 0) s << " ";
    s << a(i);
    }
  return s;
  }

// >>
template<class T>
std::istream& operator>>(std::istream& s, Array1D<T>& a)
  {
  for (unsigned int i=0; i < a.size(); ++i)
    {
    s >> a(i);
    }
  return s;
  }

//------------------------------------------------------------
// c = pow(a,b)
//
// Raise the scalar value a to the powers in the Array1D b.
// The result is an Array1D with c(i) = a^(b(i)).
//------------------------------------------------------------

template<class T, class T2>
Array1D<T>
  pow(const T2& a, const Array1D<T>& b)
  
  {
  Array1D<T> result = b;
  result.inverse_pow(a);
  result.setName("pow");
  return(result);
  }

//------------------------------------------------------------
// c = pow(a,b)
//
// Raise the Array values to the scalar power b.
// The result is an Array1D with c(i) = a(i)^b.
// This should also handle the case where b is an array
// (type T2 is Array1D<type> in which case a different .pow()
// method would be invoked.
//------------------------------------------------------------

template<class T, class T2>
Array1D<T>
  pow(const Array1D<T>& a, const T2& b)
  
  {
  Array1D<T> result = a;
  result.pow(b);
  result.setName("pow");
  return(result);
  }

//----------------------------------------------------------
// selfTest
//
// Exercize the Array1D type and verify proper operation.
//----------------------------------------------------------

template<class T>
bool Array1D<T>::selfTest()
  {
  Uvec a("a",5);
  a = Uvar(3,"m");
  Uvar b(5,"m");
  a += b;
  a.show();
  b.show();
  if (a.size() != 5) return false;
  for (unsigned int i=0; i < a.size(); ++i)
    {
    if (a(i) != Uvar(8,"m")) return(false);
    }
  a(3) = Uvar(2,"m");
  Uvar c = a(3);
  Uvec d("d"),e("e");
  d = a - c;
  e = c - a;
  for (unsigned int i=0; i < d.size(); ++i)
    {
    if (d(i) != -e(i)) return(false);
    }

  Uvec array("array",5);
  array = Uvar(-8,"km");
  Uvec array_abs("array",5);
  array_abs = ::fabs(array);
  for (unsigned int i = 0; i <array.size();++i)
    {
    if (array(i).getInUnits("km") != -array_abs(i).getInUnits("km"))
	return (false);
    }

/**
  Dvec f("f",10);
  f.indexSet();
  ArrayIndices i = select(f >= 3 && f < 9));
  Dvec g("g",i.size());
  g = f(i);
  ArrayIndices i2(select(g > 4));
  Dvec h1("h1",i2.size());
  Dvec h2("h2",i2.size());
  h1 = f(i2);
  h2 = g(i2);
  std::cout << "h1=" << h1 << std::endl;
  std::cout << "h2=" << h2 << std::endl;
**/

  return(true);
  }

//--------------------------
// Methods for class Array2D
//--------------------------

//--------------
// Constructors
//--------------

template<class T>
Array2D<T>::Array2D(const string& name, unsigned int newrows = 1,
  unsigned int newcols = 1)
 
  : name_(name)
  {
  resize(newrows,newcols);
  }

// copy constructor
template<class T>
Array2D<T>::Array2D(const Array2D<T>& a)
  : name_(a.name_), m_(a.m_)
  {
  }

// destructor
template<class T>
Array2D<T>::~Array2D()
  {
  }

//-----------
// Operators
//-----------

// unary -
template<class T>
Array2D<T>
  Array2D<T>::operator-() const
  {
    Array2D<T> retval=*this;
    retval.negate();
    return(retval);
  }

// +=
template<class T>
Array2D<T>&
  Array2D<T>::operator+=(const Array2D<T>& arg)
  
  {
  if (arg.rows() != rows())
    {
    ArrayError e("Row size mismatch adding arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.cols() != cols())
    {
    ArrayError e("Column size mismatch adding arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] += arg.m_[i][j];
    }
  return(*this);
  }

// +=
template<class T>
Array2D<T>&
  Array2D<T>::operator+=(const T& arg)
  
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] += arg;
    }
  return(*this);
  }

// -=
template<class T>
Array2D<T>&
  Array2D<T>::operator-=(const Array2D<T>& arg)
  
  {
  if (arg.rows() != rows())
    {
    ArrayError e("Row size mismatch subtracting arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.cols() != cols())
    {
    ArrayError e("Column size mismatch subtracting arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] -= arg.m_[i][j];
    }
  return(*this);
  }

// -=
template<class T>
Array2D<T>&
  Array2D<T>::operator-=(const T& arg)
  
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] -= arg;
    }
  return(*this);
  }

// *=
template<class T>
Array2D<T>&
  Array2D<T>::operator*=(const Array2D<T>& arg)
  
  {
  if (arg.rows() != rows())
    {
    ArrayError e("Row size mismatch multiplying arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.cols() != cols())
    {
    ArrayError e("Column size mismatch multiplying arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] *= arg.m_[i][j];
    }
  return(*this);
  }

// *=
template<class T>
Array2D<T>&
  Array2D<T>::operator*=(const T& arg)
  
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] *= arg;
    }
  return(*this);
  }

// /=
template<class T>
Array2D<T>&
  Array2D<T>::operator/=(const Array2D<T>& arg)
  
  {
  if (arg.rows() != rows())
    {
    ArrayError e("Row size mismatch dividing arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.cols() != cols())
    {
    ArrayError e("Column size mismatch dividing arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] /= arg.m_[i][j];
    }
  return(*this);
  }

// /=
template<class T>
Array2D<T>&
  Array2D<T>::operator/=(const T& arg)
  
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] /= arg;
    }
  return(*this);
  }

// Subscripting
template<class T>
const T& Array2D<T>::operator()(unsigned int i, unsigned int j) const
 
  {
  if (i >= rows() || j >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with indices: "
      << "(" << i << "," << j << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j]);
  }

// Subscripting
template<class T>
T& Array2D<T>::operator()(unsigned int i, unsigned int j)
 
  {
  if (i >= rows() || j >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with indices: "
      << "(" << i << "," << j << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j]);
  }

// Named method Subscripting - helps with gdb debugging
template<class T>
const T& Array2D<T>::elem(unsigned int i, unsigned int j)
  const
  {
  if (i >= rows() || j >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with indices: "
      << "(" << i << "," << j << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j]);
  }

// Assignment
template<class T>
Array2D<T>& Array2D<T>::operator=(const Array2D& a) // No exceptions
  {
  if (m_ == a.m_) return *this;  // cover self assignment
  m_ = a.m_;  // transfer data, but keep name the same.
  return(*this);
  }

template<class T>
Array2D<T>& Array2D<T>::operator=(const T& a) // No exceptions
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] = a;
    }
  return(*this);
  }

//-------------------------
// Other access methods
//-------------------------

//-----------------------------------------------------
// v = getRow(irow)
//
// Get a copy of the indicated row from this 2D Array.
// Out of bounds row causes an exception.
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getRow(unsigned int irow) const
 
  {
  if (irow >= rows())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with row index: "
      << "(" << irow << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  Array1D<T> result("",cols());
  for (unsigned int i=0; i < cols(); ++i)
    {
    result(i) = m_[irow][i];
    }
  return(result);
  }


//-----------------------------------------------------
// v = getInterpRow(irow)
//
// Linearly interpolate a row from this 2D Array.
// Out of bounds row causes an exception.
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getInterpRow(double irow) const
 
  {
  if (irow > rows()-1 || irow < 0)
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with row index: "
      << "(" << irow << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  Array1D<T> result("",cols());
  for (unsigned int i=0; i < cols(); ++i)
    {
      unsigned int row1= (unsigned int) irow;
      unsigned int row2= row1+1;
      result(i) = m_[row1][i]*(row2-irow) + m_[row2][i]*(irow-row1);
    }
  return(result);
  }

//-----------------------------------------------------
// v = getCol(icol)
//
// Get a copy of the indicated col from this 2D Array.
// Out of bounds col causes an exception.
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getCol(unsigned int icol) const
 
  {
  if (icol >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with col index: "
      << "(" << icol << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  Array1D<T> result("",rows());
  for (unsigned int i=0; i < rows(); ++i)
    {
    result(i) = m_[i][icol];
    }
  return(result);
  }


//-----------------------------------------------------
// v = getInterpCol(icol)
//
// Linearly interpolate a column from this 2D Array.
// Out of bounds column causes an exception.
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getInterpCol(double icol) const
 
  {
  if (icol > cols()-1 || icol < 0)
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with col index: "
      << "(" << icol << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  Array1D<T> result("",rows());
  for (unsigned int i=0; i < rows(); ++i)
    {
      unsigned int col1= (unsigned int) icol;
      unsigned int col2= col1+1;
      result(i) = m_[i][col1]*(col2-icol) + m_[i][col2]*(icol-col1);
    }
  return(result);
  }

//-----------------------------------------------------
// v = getPartialCol(icol, irow1, irow2)
//
// Get a copy of the part of the indicated col from this 2D Array.
// The copied part is between irow1 and irow2 inclusive
// Out of bounds col or row causes an exception.
// irow1 > irow2 throws an exception
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getPartialCol(unsigned int icol, unsigned int irow1,
				     unsigned int irow2) const
 
  {
  if (icol >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with col index: "
      << "(" << icol << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  if (irow2 >= rows())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with row index: "
      << "(" << irow2 << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  if (irow1 > irow2)
    {
      string s="Bad row index order for Array2D::GetPartialCol.";
    ArrayError e(s,ArrayError::unspecified);
    e.throwMe();
    }
  Array1D<T> result("",irow2-irow1+1);
  for (unsigned int i=irow1; i <= irow2; ++i)
    {
    int j=i-irow1;  
    result(j) = m_[i][icol];
    }
  return(result);
  }

//-----------------------------------------------------
// v = getPartialRow(irow, icol1, icol2)
//
// Get a copy of the part of the indicated row from this 2D Array.
// The copied part is between icol1 and icol2 inclusive
// Out of bounds col or row causes an exception.
// icol2 < icol1 throws an exception.
//-----------------------------------------------------

template<class T>
Array1D<T> Array2D<T>::getPartialRow(unsigned int irow, unsigned int icol1,
				     unsigned int icol2) const
 
  {
  if (irow >= rows())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with row index: "
      << "(" << irow << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  if (icol2 >= cols())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 2D array: " << name() << " with col index: "
      << "(" << icol2 << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  if (icol1 > icol2)
    {
      string s="Bad row index order for Array2D::getPartialCol.";
      ArrayError e(s,ArrayError::unspecified);
      e.throwMe();
    }
  Array1D<T> result("",icol2-icol1+1);
  for (unsigned int i=icol1; i <= icol2; ++i)
    {
    int j=i-icol1;  
    result(j) = m_[irow][i];
    }
  return(result);
  }

//-------------------------
// Predicates
//-------------------------

//-----------------------------------------------------
// bool = match(a)
//
// Check for element by element match between *this and a.
// Dimensions must also match.
//-----------------------------------------------------

template<class T>
bool Array2D<T>::match(const Array2D<T>& a) const
  {
  if (a.rows() != rows() || a.cols() != cols())
    {
    return(false);
    }

  for (unsigned int i=0; i < a.rows(); ++i)
    {
    for (unsigned int j=0; j < a.cols(); ++j)
      {
      if (a.m_[i][j] != m_[i][j])
        {
        return(false);
        }
      }
    }
  return(true);
  }

//-------------------------
// Other arithmetic methods
//-------------------------

//---------------------
// negate()
//
// Negate each element
//---------------------

template<class T>
Array2D<T>&
  Array2D<T>::negate()
  {
  for (unsigned int i = 0; i < rows(); ++i)
  for (unsigned int j = 0; j < cols(); ++j)
    {
    m_[i][j] = -m_[i][j];
    }
  return(*this);
  }

//-----------------------
// Linear algebra methods
//-----------------------

//-------------------------------------------------------
// inv()
//
// Invert this matrix.
// An exception is thrown if this matrix is not square.
//-------------------------------------------------------

template<class T>
Array2D<T>&
  Array2D<T>::inv()
  {
  if (rows() != cols())
    {
    ArrayError e("Can't invert non-square array: " + name());
    e.throwMe();
    }

  Ivec indxc("indxc",rows());
  Ivec indxr("indxr",rows());
  Ivec ipiv("ipiv",rows());
  ipiv = 0;
  unsigned int irow,icol;

  for (unsigned int i=0; i < rows(); ++i)
    {
    T big = 0;
    for (unsigned int j=0; j < rows(); ++j)
      {
      if (ipiv(j) != 1)
        {
        for (unsigned int k=0; k < rows(); ++k)
          {
          if (ipiv(k) == 0)
            {
            if (::fabs(m_[j][k]) >= big)
              {
              big = ::fabs(m_[j][k]);
              irow = j;
              icol = k;
              }
            }
          else if (ipiv(k) > 1)
            {
            ArrayError e("Singular array: " + name());
	    e.throwMe();
            }
          }
        }
      }
    ipiv(icol) += 1;
    if (irow != icol)
      {
      for (unsigned int l=0; l < rows(); ++l) swap(m_[irow][l],m_[icol][l]);
      }
    indxr(i)=irow;
    indxc(i)=icol;
    if (m_[icol][icol] == 0.0)
      {
      ArrayError e("Singular array: " + name());
      e.throwMe();
      }
    T pivinv = 1.0/m_[icol][icol];
    m_[icol][icol] = 1.0;
    for (unsigned int l=0; l < rows(); ++l) m_[icol][l] *= pivinv;
    for (unsigned int ll=0; ll < rows(); ++ll)
      {
      if (ll != icol)
        {
        T dum = m_[ll][icol];
        m_[ll][icol] = 0.0;
        for (unsigned int l=0; l < rows(); ++l) m_[ll][l] -= m_[icol][l]*dum;
        }
      }
    }

  
  for (unsigned int l = rows()-1 ; l <rows(); --l)  
    {
    if (indxr(l) != indxc(l))
      {
      for (unsigned int k=0; k < rows(); ++k)
        {
        swap(m_[k][indxr(l)],m_[k][indxc(l)]);
        }
      }
    }

  return(*this);

  }

//---------------
// Other methods
//---------------

//-------------------------------------------------
// name()
//
// Return the name of this array.
//-------------------------------------------------

template<class T>
string Array2D<T>::name() const // No exceptions
  {
  return(name_);
  }

//-------------------------------------------------
// setName()
//
// Set the name of this array.
//-------------------------------------------------

template<class T>
void Array2D<T>::setName(const string& newname) // No exceptions
  {
  name_ = newname;
  }

//-------------------------------------------------
// size(rows,cols)
//
// Return the dimensions of this matrix.
//-------------------------------------------------

template<class T>
void Array2D<T>::size(unsigned int& Nrows, unsigned int& Ncols) const // No exceptions
  {
  Nrows = m_.size();
  Ncols = m_[0].size();
  }

//-------------------------------------------------
// r = rows()
//
// Return the number of rows of this matrix.
//-------------------------------------------------

template<class T>
unsigned int Array2D<T>::rows() const // No exceptions
  {
  return(m_.size());
  }

//-------------------------------------------------
// c = cols()
//
// Return the number of columns of this matrix.
//-------------------------------------------------

template<class T>
unsigned int Array2D<T>::cols() const // No exceptions
  {
  return(m_[0].size());
  }

//-------------------------------------------------
// resize()
//
// Set new dimensions for this matrix.
//-------------------------------------------------

template<class T>
void Array2D<T>::resize(unsigned int newrows, unsigned int newcols)
  {
  if (newrows == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero rows");
    e.throwMe();
    }
  if (newcols == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero cols");
    e.throwMe();
    }
  m_.resize(newrows);
  for (unsigned int i=0; i < newrows; ++i)
    {
    m_[i].resize(newcols);
    }
  }

//------------------------------
// Supporting template functions
//------------------------------

// +
template<class T, class T2>
Array2D<T>
  operator+(const Array2D<T>& a, const T2& b)
  
  {
  Array2D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T, class T2>
Array2D<T>
  operator+(const T2& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array2D<T>
  operator+(const Array2D<T>& a, const UnitVar& b)
  
  {
  Array2D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array2D<T>
  operator+(const UnitVar& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array2D<T>
  operator+(const Array2D<T>& a, const Array2D<T>& b)
  
  {
  Array2D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// -
template<class T, class T2>
Array2D<T>
  operator-(const Array2D<T>& a, const T2& b)
  
  {
  Array2D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// -
template<class T, class T2>
Array2D<T>
  operator-(const T2& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result.negate();
  result += b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array2D<T>
  operator-(const Array2D<T>& a, const UnitVar& b)
  
  {
  Array2D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array2D<T>
  operator-(const UnitVar& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result.negate();
  result += b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array2D<T>
  operator-(const Array2D<T>& a, const Array2D<T>& b)
  
  {
  Array2D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// *
template<class T, class T2>
Array2D<T>
  operator*(const Array2D<T>& a, const T2& b)
  
  {
  Array2D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T, class T2>
Array2D<T>
  operator*(const T2& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array2D<T>
  operator*(const Array2D<T>& a, const UnitVar& b)
  
  {
  Array2D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array2D<T>
  operator*(const UnitVar& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array2D<T>
  operator*(const Array2D<T>& a, const Array2D<T>& b)
  
  {
  Array2D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

/***
// * (matrix * vector)
template<class T>
Array1D<T>
  operator*(const Array2D<T>& a, const Array1D<T>& b)
  
  {
  Array1D<T> result("mv*",b.size());
  if (a.cols() != b.size())
    {
    ArrayError e("Matrix * vector size mismatch:
    e.throwMe();
    }
  for (unsigned int i=0; i < a.rows(); ++i)

  return(result);
  }
**/

// /
template<class T, class T2>
Array2D<T>
  operator/(const Array2D<T>& a, const T2& b)
  
  {
  Array2D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T, class T2>
Array2D<T>
  operator/(const T2& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array2D<T>
  operator/(const Array2D<T>& a, const UnitVar& b)
  
  {
  Array2D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array2D<T>
  operator/(const UnitVar& b, const Array2D<T>& a)
  
  {
  Array2D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array2D<T>
  operator/(const Array2D<T>& a, const Array2D<T>& b)
  
  {
  Array2D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// <<
template<class T>
std::ostream& operator<<(std::ostream& s, const Array2D<T>& a)
  {
  for (unsigned int i=0; i < a.rows(); ++i)
    {
    for (unsigned int j=0; j < a.cols(); ++j)
      {
      if (j != 0) s << " ";
      s << a(i,j);
      }
    s << std::endl;
    }
  return s;
  }

// >>
//
// Stream input.  This is freeform input using the already
// sized Array.

template<class T>
std::istream& operator>>(std::istream& s, Array2D<T>& a)
  {
  for (unsigned int i=0; i < a.rows(); ++i)
    {
    for (unsigned int j=0; j < a.cols(); ++j)
      {
      s >> a(i,j);
      }
    }
  return s;
  }

//----------------------------------------------------------
// selfTest
//
// Exercize the Array2D type and verify proper operation.
//----------------------------------------------------------

template<class T>
bool Array2D<T>::selfTest()
  {
  Umat a("a",5,4);
  a = Uvar(3,"m");
  Uvar b(5,"m");
  a += b;
  if (a.rows() != 5) return false;
  if (a.cols() != 4) return false;
  for (unsigned int i=0; i < a.rows(); ++i)
  for (unsigned int j=0; j < a.cols(); ++j)
    {
    if (a(i,j) != Uvar(8,"m")) return(false);
    }
  a(3,3) = Uvar(2,"m");
  Uvar c = a(3,3);
  Umat d("d"),e("e");
  d = a - c;
  e = c - a;
  for (unsigned int i=0; i < d.rows(); ++i)
  for (unsigned int j=0; j < d.cols(); ++j)
    {
    if (d(i,j) != -e(i,j)) return(false);
    }
  ofstream output_filestream("array2D_test.dat");
  if (!output_filestream)
    {
    ErrorMessage e("Can't open test file: array2D_test.dat");
    e.throwMe();
    }
  output_filestream << d;
  ifstream input_filestream("array2D_test.dat");
  if (!input_filestream)
    {
    ErrorMessage e("Can't open test file: array2D_test.dat");
    e.throwMe();
    }
  input_filestream >> e;
  if (!d.match(e)) return(false);
  
  // For some strange reason, the following declaration stops compiler
  // errors for the D3D declaration in Array3D<>selfTests.  This problem
  // does not occur for U3D!
  Dmat dd("dd");
  return(true);
  }

//--------------------------
// Methods for class Array3D
//--------------------------

//--------------
// Constructors
//--------------

template<class T>
Array3D<T>::Array3D(const string& name, unsigned int newN1 = 1,
  unsigned int newN2 = 1, unsigned int newN3 = 1)
 
  : name_(name)
  {
  resize(newN1,newN2,newN3);
  }

// copy constructor
template<class T>
Array3D<T>::Array3D(const Array3D<T>& a)
  : name_(a.name_), m_(a.m_)
  {
  }

// destructor
template<class T>
Array3D<T>::~Array3D()
  {
  }

//-----------
// Operators
//-----------

// unary -
template<class T>
Array3D<T>
  Array3D<T>::operator-() const
  {
    Array3D<T> retval=*this;
    retval.negate();
    return(retval);
  }

// +=
template<class T>
Array3D<T>&
  Array3D<T>::operator+=(const Array3D<T>& arg)
  
  {
  if (arg.sizeDim1() != sizeDim1())
    {
    ArrayError e("Dimension1 size mismatch adding arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.sizeDim2() != sizeDim2())
    {
    ArrayError e("Dimension2 size mismatch adding arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] += arg.m_[i][j][k];
    }
  return(*this);
  }

// +=
template<class T>
Array3D<T>&
  Array3D<T>::operator+=(const T& arg)
  
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] += arg;
    }
  return(*this);
  }

// -=
template<class T>
Array3D<T>&
  Array3D<T>::operator-=(const Array3D<T>& arg)
  
  {
  if (arg.sizeDim1() != sizeDim1())
    {
    ArrayError e("Dimension1 size mismatch subtracting arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.sizeDim2() != sizeDim2())
    {
    ArrayError e("Dimension2 size mismatch subtracting arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] -= arg.m_[i][j][k];
    }
  return(*this);
  }

// -=
template<class T>
Array3D<T>&
  Array3D<T>::operator-=(const T& arg)
  
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] -= arg;
    }
  return(*this);
  }

// *=
template<class T>
Array3D<T>&
  Array3D<T>::operator*=(const Array3D<T>& arg)
  
  {
  if (arg.sizeDim1() != sizeDim1())
    {
    ArrayError e("Dimension1 size mismatch multiplying arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.sizeDim2() != sizeDim2())
    {
    ArrayError e("Dimension2 size mismatch multiplying arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] *= arg.m_[i][j][k];
    }
  return(*this);
  }

// *=
template<class T>
Array3D<T>&
  Array3D<T>::operator*=(const T& arg)
  
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] *= arg;
    }
  return(*this);
  }

// /=
template<class T>
Array3D<T>&
  Array3D<T>::operator/=(const Array3D<T>& arg)
  
  {
  if (arg.sizeDim1() != sizeDim1())
    {
    ArrayError e("Dimension1 size mismatch dividing arrays " + name() +
      " and " + arg.name());
    e.throwMe();
    }
  if (arg.sizeDim2() != sizeDim2())
    {
    ArrayError e("Dimension2 size mismatch dividing arrays " + name() +
      " and "  + arg.name());
    e.throwMe();
    }
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] /= arg.m_[i][j][k];
    }
  return(*this);
  }

// /=
template<class T>
Array3D<T>&
  Array3D<T>::operator/=(const T& arg)
  
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] /= arg;
    }
  return(*this);
  }

// Subscripting
template<class T>
const T& Array3D<T>::operator()(unsigned int i, unsigned int j, unsigned int k)
  const
  {
  if (i >= sizeDim1() || j >= sizeDim2() || k >= sizeDim3())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 3D array: " << name() << " with indices: "
      << "(" << i << "," << j << "," << k << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j][k]);
  }

// Subscripting
template<class T>
T& Array3D<T>::operator()(unsigned int i, unsigned int j, unsigned int k)
 
  {
  if (i >= sizeDim1() || j >= sizeDim2() || k >= sizeDim3())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 3D array: " << name() << " with indices: "
      << "(" << i << "," << j << "," << k << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j][k]);
  }

// Named method Subscripting - helps with gdb debugging
template<class T>
const T& Array3D<T>::elem(unsigned int i, unsigned int j, unsigned int k)
  const
  {
  if (i >= sizeDim1() || j >= sizeDim2() || k >= sizeDim3())
    {
    OSTRINGSTREAM os;
    os << "Bounds error accessing 3D array: " << name() << " with indices: "
      << "(" << i << "," << j << "," << k << ")";
    ArrayError e(toStr(os),ArrayError::bounds_error);
    e.throwMe();
    }
  return(m_[i][j][k]);
  }

// Assignment
template<class T>
Array3D<T>& Array3D<T>::operator=(const Array3D& a) // No exceptions
  {
  if (m_ == a.m_) return *this;  // cover self assignment
  m_ = a.m_;  // transfer data, but keep name the same.
  return(*this);
  }

template<class T>
Array3D<T>& Array3D<T>::operator=(const T& a) // No exceptions
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] = a;
    }
  return(*this);
  }

//-------------------------
// Other access methods
//-------------------------

//-------------------------
// Other arithmetic methods
//-------------------------

//---------------------
// negate()
//
// Negate each element
//---------------------

template<class T>
Array3D<T>&
  Array3D<T>::negate()
  {
  for (unsigned int i = 0; i < sizeDim1(); ++i)
  for (unsigned int j = 0; j < sizeDim2(); ++j)
  for (unsigned int k = 0; k < sizeDim3(); ++k)
    {
    m_[i][j][k] = -m_[i][j][k];
    }
  return(*this);
  }

//-----------------------
// Linear algebra methods
//-----------------------

//---------------
// Other methods
//---------------

//-------------------------------------------------
// name()
//
// Return the name of this array.
//-------------------------------------------------

template<class T>
string Array3D<T>::name() const // No exceptions
  {
  return(name_);
  }

//-------------------------------------------------
// setName()
//
// Set the name of this array.
//-------------------------------------------------

template<class T>
void Array3D<T>::setName(const string& newname) // No exceptions
  {
  name_ = newname;
  }

//------------------------------------------------------
// size(N1,N2,N3)
//
// Return the sizes of the dimensions of this 3D array.
//------------------------------------------------------

template<class T>
void Array3D<T>::size(unsigned int& N1, unsigned int& N2, unsigned int& N3)
  const // No exceptions
  {
  N1 = m_.size();
  N2 = m_[0].size();
  N3 = m_[0][0].size();
  }

//-------------------------------------------------
// d1 = sizeDim1()
// d2 = sizeDim2()
// d3 = sizeDim3()
//
// Return the sizes of each dimension separately.
//-------------------------------------------------

template<class T>
unsigned int Array3D<T>::sizeDim1() const // No exceptions
  {
  return(m_.size());
  }

template<class T>
unsigned int Array3D<T>::sizeDim2() const // No exceptions
  {
  return(m_[0].size());
  }

template<class T>
unsigned int Array3D<T>::sizeDim3() const // No exceptions
  {
  return(m_[0][0].size());
  }

//-------------------------------------------------
// resize()
//
// Set new dimensions for this matrix.
//-------------------------------------------------

template<class T>
void Array3D<T>::resize(unsigned int newN1, unsigned int newN2,
  unsigned int newN3)
  {
  if (newN1 == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero for dim 1");
    e.throwMe();
    }
  if (newN2 == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero for dim 2");
    e.throwMe();
    }
  if (newN3 == 0)
    {
    ArrayError e("Can't make array: " + name() + " with zero for dim 3");
    e.throwMe();
    }
  m_.resize(newN1);
  for (unsigned int i=0; i < newN1; ++i)
    {
    m_[i].resize(newN2);
    for (unsigned int j=0; j < newN2; ++j)
      {
      m_[i][j].resize(newN3);
      }
    }
  }

//------------------------------
// Supporting template functions
//------------------------------

// +
template<class T, class T2>
Array3D<T>
  operator+(const Array3D<T>& a, const T2& b)
  
  {
  Array3D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T, class T2>
Array3D<T>
  operator+(const T2& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array3D<T>
  operator+(const Array3D<T>& a, const UnitVar& b)
  
  {
  Array3D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array3D<T>
  operator+(const UnitVar& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// +
template<class T>
Array3D<T>
  operator+(const Array3D<T>& a, const Array3D<T>& b)
  
  {
  Array3D<T> result = a;
  result += b;
  result.setName("+");
  return(result);
  }

// -
template<class T, class T2>
Array3D<T>
  operator-(const Array3D<T>& a, const T2& b)
  
  {
  Array3D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// -
template<class T, class T2>
Array3D<T>
  operator-(const T2& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result.negate();
  result += b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array3D<T>
  operator-(const Array3D<T>& a, const UnitVar& b)
  
  {
  Array3D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array3D<T>
  operator-(const UnitVar& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result.negate();
  result += b;
  result.setName("-");
  return(result);
  }

// -
template<class T>
Array3D<T>
  operator-(const Array3D<T>& a, const Array3D<T>& b)
  
  {
  Array3D<T> result = a;
  result -= b;
  result.setName("-");
  return(result);
  }

// *
template<class T, class T2>
Array3D<T>
  operator*(const Array3D<T>& a, const T2& b)
  
  {
  Array3D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T, class T2>
Array3D<T>
  operator*(const T2& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array3D<T>
  operator*(const Array3D<T>& a, const UnitVar& b)
  
  {
  Array3D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array3D<T>
  operator*(const UnitVar& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// *
template<class T>
Array3D<T>
  operator*(const Array3D<T>& a, const Array3D<T>& b)
  
  {
  Array3D<T> result = a;
  result *= b;
  result.setName("*");
  return(result);
  }

// /
template<class T, class T2>
Array3D<T>
  operator/(const Array3D<T>& a, const T2& b)
  
  {
  Array3D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T, class T2>
Array3D<T>
  operator/(const T2& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array3D<T>
  operator/(const Array3D<T>& a, const UnitVar& b)
  
  {
  Array3D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array3D<T>
  operator/(const UnitVar& b, const Array3D<T>& a)
  
  {
  Array3D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

// /
template<class T>
Array3D<T>
  operator/(const Array3D<T>& a, const Array3D<T>& b)
  
  {
  Array3D<T> result = a;
  result /= b;
  result.setName("/");
  return(result);
  }

/** some strange template problem with << and >>.
// <<
template<class T>
std::ostream& operator<<(std::ostream& s, const Array3D<T>& a)
  {
  for (unsigned int i=0; i < a.sizeDim1(); ++i)
    {
    for (unsigned int j=0; j < a.sizeDim2(); ++j)
      {
      for (unsigned int k=0; k < a.sizeDim3(); ++k)
        {
        if (k != 0) s << " ";
        s << a(i);
        }
      s << std::endl;
      }
    s << std::endl;
    }
  return s;
  }
**/

//----------------------------------------------------------
// selfTest
//
// Exercize the Array3D type and verify proper operation.
//----------------------------------------------------------

template<class T>
bool Array3D<T>::selfTest()
  {
  U3D a("a",5,4,3);
  a = Uvar(3,"m");
  Uvar b(5,"m");
  a += b;
  unsigned int N1,N2,N3;
  a.size(N1,N2,N3);
  if (N1 != 5) return false;
  if (N2 != 4) return false;
  if (N3 != 3) return false;
  for (unsigned int i=0; i < a.sizeDim1(); ++i)
  for (unsigned int j=0; j < a.sizeDim2(); ++j)
  for (unsigned int k=0; k < a.sizeDim3(); ++k)
    {
    if (a(i,j,k) != Uvar(8,"m")) return(false);
    }
  a(3,3,2) = Uvar(2,"m");
  Uvar c = a(3,3,2);
  U3D d("d"),e("e");
  d = a - c;
  e = c - a;
  for (unsigned int i=0; i < d.sizeDim1(); ++i)
  for (unsigned int j=0; j < d.sizeDim2(); ++j)
  for (unsigned int k=0; k < d.sizeDim3(); ++k)
    {
    if (d(i,j,k) != -e(i,j,k)) return(false);
    }
  // Below fails at compile time if Array2D::selfTest does not have
  // a Dmat declaration of some kind!
  D3D dd("dd");
  dd.resize(2,3,4);
  return(true);
  }

//---------------------------
// Other Supporting functions
//---------------------------

// min
template<class T>
T min(const Array1D<T>& v, unsigned int& idx)
  {
  T result = v.min(idx);
  return(result);
  }

// max
template<class T>
T max(const Array1D<T>& v, unsigned int& idx)
  {
  T result = v.max(idx);
  return(result);
  }

// sum
template<class T>
T sum(Array1D<T> v)
  {
  T result = v.sum();
  return(result);
  }

template<class T>
Array1D<T> fabs(Array1D<T> v)
  {
  Array1D<T> result = v;
  result.fabs();
  return(result);
  }
