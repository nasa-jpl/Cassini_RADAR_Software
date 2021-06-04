//=============================================================================
// Io.h
//
// This file contains the FileMgr class declarations and some i/o supporting
// functions.  The FileMgr class provides a convenient handler for heterogenous
// binary data files that may need to move from big-endian systems to
// little-endian systems or vica-versa.
//
// Interface summary:
//
// FileMgr methods:
//
//   Construction:
//
//   FileMgr(filename,mode); Construct opening filename in indicated mode
//     Because the constructor opens a file, a copy-constructor, assignment
//     operator, and destructor are also defined.  The copy-constructor and
//     assignment operator make sure that one file is not opened twice, and
//     the destructor closes the file.
//
//   State control:
//
//   machineType();  -- to be removed
//   setByteOrderReversed();
//   setByteOrderNormal();
//
//   File control and status:
//
//   rewind();
//   name();
//   eof();
//   getPosition();
//   setPosition(pos);
//
//   I/O:
//
//   read(arg);           Read argument from binary file based on argument type
//   write(arg);          Write argument to binary file based on argument type
//   readKnown(arg);      Read a known value and set byte-flipping status
//
// Error handling: FileMgr::IoError
//
// Supporting functions:
//
//   byte_reverse(d);             Reverse byte order of d
//   byte_int(d,order,type,b);    Convert bytes in b to int returned in d
//   byte_int(dvec,order,type,b); Convert bytes in b to int returned in dvec
//   i = bitget(ii,b1,b2);        Extract bit positions b1..b2 from ii
//   binary_read(file,d);         function to access method read
//
//=============================================================================

#ifndef Io_H
#define Io_H

#include <stdio.h>
#include <string>
#include <vector>
#include "Error.h"
class UnitVar;


using std::string;
using std::vector;
using std::swap;

int chrbitget(const unsigned char  c, const int b1, const int b2); // temporary use to get bit information

class FileMgr;
template<class T> void byte_reverse(T& d);

//-------------------------------
// File Management               
//-------------------------------

class FileMgr
  {
  //----------------
  // private methods
  //----------------
  private:

  FileMgr(const FileMgr& fp);
  FileMgr& operator= (const FileMgr& fp);

  //----------------
  // public methods
  //----------------
  public:

  enum orderE {little_endian, big_endian, unknown};
  enum typeE {unsigned_int, signed_int};

  //---------------------
  // Error handling class
  //---------------------

//  class IoError;

  class IoError : public ErrorMessage
    {
    public:
      IoError(const string& filename, const string& smsg) // No exceptions
      {
      msg = "File Access Error for: " + filename + "\n  " + smsg;
      }
      void throwMe(){ throw *this; }
    };

  // Constructors/Destructors
  FileMgr(); // default constructor
  FileMgr(const string& name, const string& attrib);
  ~FileMgr();

  // routines for closing and reopening a file
  void reopen(const string& name, const string& attrib);
  void close();


  // Other Methods
  bool eof() const;
  int getPosition() const;
  void setPosition(int position);
  orderE machineType() const;
  void setByteOrderReversed();
  void setByteOrderNormal();
  void rewind();
  string name() const;

  //--------------------------
  // String read/write methods
  //--------------------------

  void read(string& s) const;
  void write(const string& s) const;
  void readXYZ(double& x, double& y, double& z);
  void writeXYZ(const double& x, const double& y, const double& z);

  // differs from Uvar read/write in that non base units I/O is enabled
  // despite performance hit
  //float method
  void readFloatInUnits(UnitVar& uv, double conversion_factor, const string& ustr);
  void readDoubleInUnits(UnitVar& uv, double conversion_factor, const string& ustr);

  void writeFloatInUnits(const UnitVar& uv, const string& ustr);
  void writeDoubleInUnits(const UnitVar& uv, const string& ustr);


  //--------------------------
  // read a line from an ASCII file
  //--------------------------

  void readLine(string& s) const;

  //--------------------------
  // read the next string from an ASCII file
  // ignoring whitespace
  //--------------------------

  string parseNextString() const;

  //----------------------------------------
  // Generic single datum read/write methods
  //----------------------------------------

  // read
  template<class T>
  void read(T& d) const
    {
    if(!p_) {
      FileMgr::IoError(filename_,"Error reading from unopened file").throwMe();
    }
    int n = fread(&d,sizeof(d),1,p_);
    if (n != 1){
      FileMgr::IoError e(filename_,"Read Error");
      e.throwMe();
    }
    if (byte_reversal_) byte_reverse(d);
    return;
    }

  // write
  template<class T>
    void write(const T& d) const
    {
    if(!p_) {
      IoError(filename_,"Error writing to unopened file").throwMe();
    }
    int n = fwrite(&d,sizeof(d),1,p_);
    if (n != 1) { 
      IoError e(filename_,"Write Error");
      e.throwMe();
    } 
   return;
    }

  // read a vector
  template<class T>
  void read(vector<T>& d) const
    {
    if(!p_) {
      IoError(filename_,"Error reading from unopened file").throwMe();
    }
    for (typename vector<T>::iterator p = d.begin(); p != d.end(); ++p)
      {  // read each element
      read(*p);
      }
    }

  //write a vector
  template <class T>
  void write(vector<T>& d) const
  {
    if(!p_) {
      IoError(filename_,"Error writing to unopened file").throwMe();
    }
    for (typename vector<T>::iterator p=d.begin(); p !=d.end(); ++p)
      {//write each element
	write(*p);
      }
  }

  //----------------------------------------------------------------
  // Read a known datum and set byte order reversal and machine type
  // appropriately.
  //----------------------------------------------------------------

  // readKnown
  template<class T>
  void readKnown(const T& d)
    {
    if(!p_) {
      IoError(filename_,"Error reading from unopened file").throwMe();
    }
    int pos=getPosition();
    vector<unsigned char> b(sizeof(d));
    read(b);
    T d1,d2;
    byte_int(d1,big_endian,unsigned_int,b);
    byte_int(d2,little_endian,unsigned_int,b);
    if (d1 == d) machine_type_ = big_endian;
    else if (d2 == d) machine_type_ = little_endian;
    else IoError(filename_,"Could not match known datum").throwMe();

    setPosition(pos);
    T d0;
    read(d0);
    if (d0 != d)
      {
      byte_reverse(d0);
      if (d0 == d) setByteOrderReversed();
      else IoError(filename_,"Could not match known datum (2)").throwMe();
      }
    }

  protected:

  // Internal representation

  FILE *p_;
//  std::fstream *p_;  need binary i/o on streams to work before I use this.
  bool byte_reversal_;
  orderE machine_type_;
  string filename_;

  };

//------------------------------
// Supporting template functions
//------------------------------

//-------------------------------------------------------------
// byte_reverse(d)
//
// Reverse the byte order of the datum d (passed by reference).
// This is intended for simple types like int, double, etc.
//-------------------------------------------------------------

template<class T>
void byte_reverse(T& d)
  {
  char* p = (char *)&(d);
  for (unsigned int i=1; i <= sizeof(T)/2; i++)
    {
    swap(*p, *(p + sizeof(T) - (i*2-1)));
    ++p;
    } 
  }

//----------------------------------------------------------------------------
// byte_int(d,order,type,b)
//
// Convert data from a byte array into other integer forms with possible byte
// order swapping.
// This is intended for simple types like int, double, etc.
//
// d = converted integer data
// order = big_endian or little_endian
// type = 'signed' or 'unsigned'
// b = vector of byte data
//----------------------------------------------------------------------------

template<class Tint_type>
void byte_int(Tint_type& d, FileMgr::orderE order, FileMgr::typeE type,
  const vector<unsigned char>& b, unsigned int start_pos = 0,
  unsigned int num_elems = sizeof(Tint_type))
  {
  if (sizeof(Tint_type) < num_elems)
    ErrorMessage("Error: return integer too small in byte_int").throwMe();
  if (start_pos + num_elems > b.size())
    ErrorMessage("Error: Invalid subrange in byte_int").throwMe();

  d = 0;

  if (order == FileMgr::big_endian)
    {
    int n = 0;
    vector<unsigned char>::const_iterator p = b.begin() + start_pos;
    for (unsigned int i=0; i < num_elems; ++i)
      {
      d += (*p) << n;  // add in byte contribution based on position
      n += 8;
      ++p;
      }
    }
  else if (order == FileMgr::little_endian)
    {
    int n = 0;
    vector<unsigned char>::const_iterator p = b.begin() + start_pos +
      num_elems - 1;
    for (unsigned int i=0; i < num_elems; ++i)
      {
      d += (*p) << n;  // add in byte contribution based on position
      n += 8;
      --p;
      }
    }
  else
    {
      ErrorMessage e("Error: Invalid machine type in byte_int");
      e.throwMe();
    }

  if (sizeof(Tint_type) > num_elems)
    {  // destination is larger than the source byte array
    Tint_type maxint = 1;
    maxint <<= 8*num_elems;
    if (type == FileMgr::signed_int && d >= maxint/2)
      {
      d = -(maxint - d);  // fix sign appropriately
      }
    }
  }

//---------------------------------------------------------
// byte_int(dvec,order,type,b)
//
// Like byte_int above, but generates a vector of results.
//---------------------------------------------------------

template<class Tint_type>
void byte_int(vector<Tint_type>& d,
  FileMgr::orderE order, FileMgr::typeE type,
  const vector<unsigned char>& b, unsigned int start_pos = 0,
  unsigned int num_elems = sizeof(Tint_type))
  {
  if (num_elems % sizeof(Tint_type) != 0)
    ErrorMessage("Error: non-integral number of elements in byte_int").throwMe();

  unsigned int Nints = num_elems/sizeof(Tint_type);
  d.clear();
  for (unsigned int i=0; i < Nints; ++i)
    {
    Tint_type d1;
    byte_int(d1,order,type,b,start_pos + i*sizeof(Tint_type));
    d.push_back(d1);
    }
  }

//--------------------------------------------------------------------------
// i = bitget(ii,b1,b2)
//
// Extract the bits from bit position b1..b2 inclusive (zero-offset) from
// the integer ii.  The bits are returned as the least
// significant bits of an integer of the same type as ii (all others zeroed).
//--------------------------------------------------------------------------

template<class Tint_type>
Tint_type bitget(Tint_type ii, unsigned int b1, unsigned int b2)
  {
  if (sizeof(Tint_type) > sizeof(int))
    ErrorMessage("Integer too long for bitget").throwMe();
  if (b2 < b1) ErrorMessage("Bad position range in bitget").throwMe();
//  cout << ii << " " << b1 << " " << b2 << endl;

  if (b2 < 8*sizeof(Tint_type) - 1)
    {  // There are bits above b2.
    // Turn off bits below b2 by shifting them out.
    Tint_type d = ii >> b2+1;
//    cout << d << " ";
    d <<= b2+1;
//    cout << d << endl;

    // Now subtract off the bits above b2.
    ii -= d;
//    cout << ii << endl;
    }

  // Right shift by b1 leaving bits b1..b2 in the least sig. bit position.
  ii >>= b1;
//  cout << ii << endl;

  return(ii);
  }




//----------------------------------------//
// I/O Supporting functions and operators //
//----------------------------------------//

// binary_read(file,d)
//
// Read a simple data type using a FileMgr object.
// This function is provided to support templatized objects that need to
// read an object that may be simple or more complex.  The simple cases
// are handled by this templatized version of binary_read.  FileMgr already
// knows how to handle the simple data types (float, int, etc).
// More complex types will either need an explicit specialization of
// binary_read, or FileMgr will need to understand them.
// Alternatively, just provide the more complex types with their own binary
// read method and don't use binary_read at all.

//template<class T>
//void binary_read(FileMgr& file, T& d)
//  {
//  file.ReadSingle(d);
//  }

#endif
