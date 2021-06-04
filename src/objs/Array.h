//==============================================================================
// Array.h
//
// This file contains the interface classes and functions for a basic
// vector, matrix handling system.
// The public interface consists of the classes Array1D and Array2D.
// This header comment summarizes the interface.
// For details about a specific function, look at the declarations in this file.
//
// Interface summary:
//
// class Array1D;
//   Vector of data with bounds checking.
//
// class Array2D;
//   Matrix of data with bounds checking.
//
// Array1D Methods:
//
//   Construction:
//
//   Array1D(name,length);    Setup uninitialized vector and associated name.
//
//   Internal operators: +=, -=, *=, /=
//
//   Predicates:
//
//   I/O:
//   void write(file,usize); Write binary file with 3 options for units.
//
//   Other methods:
//   str = name();       Get name of this Array
//   scale(mag);         Scale array to have indicated magnitude
//
// Other supporting functions:
//
// External Binary operators: +, -, *, /, ==
// ASCII stream I/O operators: <<, >>
//==============================================================================

#ifndef Array_H
#define Array_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <sstream>

using std::string;
using std::istringstream;
using std::ostringstream;
using std::complex;

//----------------------
// Forward declarations
//----------------------

template <class T> class Array1D;
template <class T> class Array2D;
template <class T> class Array3D;
class ArrayError;

//--------------------------------------------------------------------
// The following typedef's are provided for commonly used Array types
//--------------------------------------------------------------------

#include "Units.h"

typedef Array1D<Uvar> Uvec;
typedef Array2D<Uvar> Umat;
typedef Array3D<Uvar> U3D;
typedef Array1D<double> Dvec;
typedef Array1D<char> Charvec;
typedef Array1D<unsigned char> CharUvec;
typedef Array2D<double> Dmat;
typedef Array3D<double> D3D;
typedef Array1D<float> Fvec;
typedef Array1D<float> fdata;
typedef Array2D<float> Fmat;
typedef Array3D<float> F3D;
typedef Array1D<complex<float> > CFvec;
typedef Array2D<complex<float> > CFmat;
typedef Array3D<complex<float> > CF3D;
typedef Array1D<complex<double> > CDvec;
typedef Array2D<complex<double> > CDmat;
typedef Array3D<complex<double> > CD3D;
typedef Array1D<int> Ivec;
typedef Array2D<int> Imat;
typedef Array3D<int> I3D;
typedef Array1D<unsigned short> SUvec;
typedef Array1D<unsigned long> LUvec;
typedef Array2D<unsigned long> LUmat;
typedef Array3D<unsigned long> LU3D;

#include "Error.h"
#include "Io.h"
#include "Utils.h"


//--------------------------------------------------------------------------
// IO declarations that are declared as friend functions, and therefore
// need to be declared here ahead of the class definitions
//--------------------------------------------------------------------------

// >>
template<class T>
std::istream& operator>>(std::istream& s, Array1D<T>& v);

// <<
template<class T>
std::ostream& operator<<(std::ostream& s, const Array1D<T>& v);

// >>
template<class T>
std::istream& operator>>(std::istream& s, Array2D<T>& v);

// <<
template<class T>
std::ostream& operator<<(std::ostream& s, const Array2D<T>& v);

//--------------------------
// Class Array1D declaration
//--------------------------

template<class T>
class Array1D
  {
  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;
  
  //--------------
  // Construction
  //--------------
  
  Array1D(const string& name, unsigned int length);
  Array1D(const Array1D& a);  // copy constructor
  ~Array1D();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  Array1D operator-() const;

  Array1D& operator+=(const Array1D& arg) ;
  Array1D& operator-=(const Array1D& arg) ;
  Array1D& operator*=(const Array1D& arg) ;
  Array1D& operator/=(const Array1D& arg) ;

  Array1D& operator+=(const T& arg) ;
  Array1D& operator-=(const T& arg) ;
  Array1D& operator*=(const T& arg) ;
  Array1D& operator/=(const T& arg) ;

//  Array1D operator/(const T& arg) ;

//  bool operator==(const Array1D& b) const;

  // Subscripting
  const T& operator()(unsigned int i) const;
  T& operator()(unsigned int i);
  // elem is provided to aid access in gdb - it is the same as const operator()
  const T& elem(unsigned int i) const;

  // Assignment
  template <class T2>
    Array1D& operator=(const Array1D<T2> a) // No excepts
    {
    resize(a.size());
    for (unsigned int i = 0; i < a.size(); ++i)
      {
      v_[i] = a(i);  // rely on assignment operator between T2 and T.
      }
    return(*this);
    }
  Array1D& operator=(const Array1D& a); // NO excepts
  Array1D& operator=(const T& a); // NO excepts

  //------------
  // Predicates
  //------------

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  void write(const FileMgr& file, int usize = -1) const ;
  void read(const FileMgr& file, int usize = -1) ;
  // For the read option below, the elements are coerced to have the units
  // contained in ustr after being read from the file.
  void read(const FileMgr& file, int usize, const string& ustr);
  void show() const;
  void show(unsigned int i) const;

  friend ostream& operator<< <>(ostream& s, const Array1D& v);
  friend istream& operator>> <>(istream& s, Array1D& v);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  Array1D& negate();
  Array1D& pow(const T& b);
  Array1D& pow(const Array1D<T>& b);
  Array1D& inverse_pow(const T& base);

  //---------------
  // Other methods
  //---------------

  T magnitude() const ; 
  T min(unsigned int& idx) const;
  T max(unsigned int& idx) const;
  T sum() const;
  T mean() const;
  T rms(unsigned int n=size()) const; 
  void fabs() ;
  void scale(T mag) ;
  string name() const; // NO excepts
  void setName(const string& newname); // NO excepts
  unsigned int size() const; // NO excepts
  void resize(unsigned int newsize);
  Array1D<double> findInterpolatedIndices(const T& y) 
    const;
  private:

  string name_;
  vector<T> v_;

  };

//--------------------------
// Class Array2D declaration
//--------------------------

template<class T>
class Array2D
  {
  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  //--------------
  // Construction
  //--------------

  Array2D(const string& name, unsigned int rows, unsigned int cols)
    ;
  Array2D(const Array2D& a);  // copy constructor
  ~Array2D();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  Array2D operator-() const;

  Array2D& operator+=(const Array2D& arg) ;
  Array2D& operator-=(const Array2D& arg) ;
  Array2D& operator*=(const Array2D& arg) ;
  Array2D& operator/=(const Array2D& arg) ;

  Array2D& operator+=(const T& arg) ;
  Array2D& operator-=(const T& arg) ;
  Array2D& operator*=(const T& arg) ;
  Array2D& operator/=(const T& arg) ;

  bool operator==(const Array2D& b) const;

  // Subscripting
  const T& operator()(unsigned int i, unsigned int j) const ;
  T& operator()(unsigned int i, unsigned int j) ;
  // elem is provided to aid access in gdb - it is the same as const operator()
  const T& elem(unsigned int i, unsigned int j) const
    ;

  // Assignment
  Array2D& operator=(const Array2D& a); // NO excepts
  Array2D& operator=(const T& a); // NO excepts

  //----------------------
  // Other access methods
  //----------------------

  Array1D<T> getRow(unsigned int irow) const ;
  Array1D<T> getInterpRow(double irow) const ;
  Array1D<T> getCol(unsigned int icol) const ;
  Array1D<T> getInterpCol(double irow) const ;
  Array1D<T> getPartialCol(unsigned int icol, unsigned int irow1,
			   unsigned int irow2) const ;
  Array1D<T> getPartialRow(unsigned int irow, unsigned int icol1,
			   unsigned int icol2) const ;


  //------------
  // Predicates
  //------------

  bool match(const Array2D& a) const;

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

  friend ostream& operator<< <>(ostream& s, const Array2D& v);
  friend istream& operator>> <>(istream& s, Array2D& v);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  Array2D& negate();

  //-------------------------
  // Linear algebra methods
  //-------------------------

  Array2D& inv();

  //---------------
  // Other methods
  //---------------

  string name() const; // NO excepts
  void setName(const string& newname); // NO excepts
  void size(unsigned int& Nrows, unsigned int& Ncols) const; // NO excepts
  unsigned int rows() const; // NO excepts
  unsigned int cols() const; // NO excepts
  void resize(unsigned int newrows, unsigned int newcols);

  private:

  string name_;
  vector<vector<T> > m_;

  };

//--------------------------
// Class Array3D declaration
//--------------------------

template<class T>
class Array3D
  {
  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  //--------------
  // Construction
  //--------------

  Array3D(const string& name, unsigned int N1, unsigned int N2, unsigned int N3)
    ;
  Array3D(const Array3D& a);  // copy constructor
  ~Array3D();

  //---------
  // Testing
  //---------

  static bool selfTest();

  //----------------
  // Operators
  //----------------

  Array3D operator-() const;

  Array3D& operator+=(const Array3D& arg) ;
  Array3D& operator-=(const Array3D& arg) ;
  Array3D& operator*=(const Array3D& arg) ;
  Array3D& operator/=(const Array3D& arg) ;

  Array3D& operator+=(const T& arg);
  Array3D& operator-=(const T& arg);
  Array3D& operator*=(const T& arg);
  Array3D& operator/=(const T& arg);

  bool operator==(const Array3D& b) const;

  // Subscripting
  const T& operator()(unsigned int i, unsigned int j, unsigned int k) const;
  T& operator()(unsigned int i, unsigned int j, unsigned int k);
  // elem is provided to aid access in gdb - it is the same as const operator()
  const T& elem(unsigned int i, unsigned int j, unsigned int k) const;

  // Assignment
  Array3D& operator=(const Array3D& a); // No exceptions
  Array3D& operator=(const T& a); // No exceptions

  //----------------------
  // Other access methods
  //----------------------

  //------------
  // Predicates
  //------------

  //-------------------------------------------------------------
  // I/O
  //-------------------------------------------------------------

//  friend ostream& operator<< <>(ostream& s, const Array3D& v);
//  friend istream& operator>> <>(istream& s, Array3D& v);

  //--------------------------
  // Other arithmetic methods
  //--------------------------

  Array3D& negate();

  //-------------------------
  // Linear algebra methods
  //-------------------------

  //---------------
  // Other methods
  //---------------

  string name() const; // No exceptions
  void setName(const string& newname); // No exceptions
  void size(unsigned int& N1, unsigned int& N2, unsigned int& N3) const; // No exceptions
  unsigned int sizeDim1() const; // No exceptions
  unsigned int sizeDim2() const; // No exceptions
  unsigned int sizeDim3() const; // No exceptions
  void resize(unsigned int newN1, unsigned int newN2, unsigned int newN3);

  private:

  string name_;
  vector<vector<vector<T> > > m_;

  };

//--------------------------------
// Class ArrayIndices declaration
//--------------------------------

class ArrayIndices
  {
  public:

  //--------------
  // Construction
  //--------------

  ArrayIndices();

  private:

  //-------------------------
  // Internal representation
  //-------------------------

  unsigned int level_;
  vector<unsigned int>* i_;

  };

//-----------------------------------------------------------------
// Error handling class for Array1D, 2D, 3D (including definitions)
// The error type can be obtained from the public variable
// error_type, or by reading the error message in msg.
//-----------------------------------------------------------------

//---------------------
// Error handling class
//---------------------

class ArrayError : public ErrorMessage
  {
  public:

  enum errorE {unspecified, scale_zero, read_error, write_error, bounds_error};

  // Constructors

  ArrayError(errorE err_type = unspecified) // no exceptions
    : error_type(err_type)
    {
    if (error_type == unspecified)
      msg = "Unspecified Array Error";
    else if (error_type == scale_zero)
      msg = "Array tried to scale a zero vector";
    else if (error_type == read_error)
      msg = "Array read error";
    else if (error_type == write_error)
      msg = "Array write error";
    else if (error_type == bounds_error)
      msg = "Array bounds error";
    }

  ArrayError(const string& emsg, errorE err_type = unspecified) // no excepts
    : error_type(err_type)
    {
    msg = emsg;
    }
  void throwMe() { throw *this;}
  // Public type flag
  errorE error_type;
  };

//------------------------------
// Binary Operators for Arrays.
//------------------------------

/**
// + - * /
template<class T>
Array1D<T> operator+(const Array1D<T>& a, const Array1D<T>& b);
template<class T>
Array1D<T> operator-(const Array1D<T>& a, const Array1D<T>& b);
template<class T>
Array1D<T> operator*(const Array1D<T>& a, const Array1D<T>& b);
template<class T>
Array1D<T> operator/(const Array1D<T>& a, const Array1D<T>& b);
**/

//---------------------------
// Other Supporting functions
//--------------------------

template<class T>
T min(const Array1D<T>& v, unsigned int& idx);

template<class T>
T max(const Array1D<T>& v, unsigned int& idx);

template<class T>
T sum(Array1D<T> v);

template<class T>
Array1D<T>  fabs(Array1D<T> v);

ArrayIndices select(const Array1D<bool>& b);

//-----------------------------------------------------------------------------
// The remaining declarations and required definitions are placed in a separate
// file.  They are not considered part of the normal interface.
//-----------------------------------------------------------------------------

#include "Array_impl.h"

#endif


