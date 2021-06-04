//=============================================================================
// MFile.h
//
// This file contains the MFile class declarations.  The MFile class 
// handles outputs to a MATLAB readable script file
//
// Interface summary:
//
//-----------------------------------------------------------------------
#ifndef MFILE_H
#define MFILE_H

#include<stdio.h>
#include<complex>
#include<string>
using std::string;
using std::complex;

class MFile
{
 public:
  
  // constructor and destructor
  MFile();
  ~MFile();
  
  // open and close
  void close();
  void open(const char* s);
  void flush();

  // output a comment
  void comment(const string& s);
  
  // routines for outputing set statements (arrays and scalars of various
  // types
  void set(const string & left, float v, int d1=0, int d2=0);
  void set(const string & left, complex<float>* v, int len, int d1=0, 
	   int d2=0);
  void set(const string & left, float* v, int len, int d1=0, int d2=0);
  void leftHandSide(const string & left, int d1=0, int d2=0);

private:

  void dieOnUnopenedFile(const string& funcname);

FILE* fp;

};
#endif
