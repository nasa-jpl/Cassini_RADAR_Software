#include <stdio.h>
#include <iostream>
#include <string>
#include "Io.h"
#include "Units.h"

//------------------------//
// Class FileMgr Methods //
//------------------------//

// Constructors/Destructors
FileMgr::FileMgr() 
  : byte_reversal_(false), machine_type_(unknown),
    filename_("")
  {
  p_ = NULL;
  }

FileMgr::FileMgr(const string& name, const string& attrib) 
  : byte_reversal_(false), machine_type_(unknown),
    filename_(name)
  {
  p_ = fopen(name.c_str(),attrib.c_str());
  if (!p_) IoError(name,"Open Error").throwMe();
  }

// Copy constructor
FileMgr::FileMgr(const FileMgr& fp)
  {
    IoError e(fp.filename_,"Do NOT use copy constructor with FileMgr");
    e.throwMe();
  }

FileMgr::~FileMgr() { if (p_) fclose(p_); }

// Assignment operator
FileMgr& FileMgr::operator= (const FileMgr& fp)
  {
    IoError e(fp.filename_, "Do NOT use assignment operator with FileMgr");
    e.throwMe();
    return(*this);
  }

void FileMgr::close(){
  if (p_) fclose(p_);
  p_=NULL;
}

void FileMgr::reopen(const string& name, const string& attrib) 
  {
    close();
    p_ = fopen(name.c_str(),attrib.c_str());
    if (!p_) IoError(name,"Open Error").throwMe();  
  }
//--------------
// Other Methods
//--------------

FileMgr::orderE FileMgr::machineType() const { return(machine_type_); }
void FileMgr::setByteOrderReversed() { byte_reversal_ = true; }
void FileMgr::setByteOrderNormal() { byte_reversal_ = false; }
void FileMgr::rewind() { 
  if(p_) setPosition(0);
  else IoError(filename_,"rewind:File not open").throwMe();
}
string FileMgr::name() const { return(filename_); }

//---------------------------------------------
// eof()
// Check for end of file on the next read.
// Throws an exception if an i/o error occurs.
//---------------------------------------------

bool FileMgr::eof() const 
  {
  if(!p_) {
    IoError(filename_,"Error checking for end of unopened file").throwMe();
  }
  int c = fgetc(p_);
  bool e = feof(p_);
  if (e) return(e);
  if (ungetc(c,p_) == EOF)
    {
    IoError(filename_,"Error checking for end of file").throwMe();
    }
  return(e);
  }

//---------------------------------------------------------
// getPosition()
// Returns the integer offset of the current file position.
// Throws an exception if an i/o error occurs.
//---------------------------------------------------------

int FileMgr::getPosition() const 
  {
  if(!p_) {
    IoError(filename_,"Cannot get position in unopened file").throwMe();
  }
  int position = ftell(p_);
  if (position < 0)
    {
    IoError e(filename_,"Error getting current file position");
    e.throwMe();
    }
  return(position);
  }

//--------------------------------------------------------------
// setPosition()
// Sets the current file position to the supplied integer offset.
// The offset is with respect to the start of the file.
// Throws an exception if an i/o error occurs.
//--------------------------------------------------------------

void FileMgr::setPosition(int position) 
  {
  if(!p_) {
    IoError(filename_,"Cannot set position in unopened file").throwMe();
  }
  if (fseek(p_,position,0) != 0)
    {
    IoError e(filename_,"Error setting current file position");
    e.throwMe();
    }
  }

//------------------------------------------------
// read
// Read in characters to fill the supplied string.
//------------------------------------------------

void FileMgr::read(string& s) const 
  {
  if(!p_) {
    IoError(filename_,"Error reading unopened file").throwMe();
  }
  for (unsigned int i=0; i < s.length(); i++)
    {
    char c;
    int n = fread(&c,sizeof(char),1,p_);
    if (n != 1) IoError(filename_,"Read Error").throwMe();
    s[i] = c;
    }
  }

// read XYZ
void
FileMgr::readXYZ(double & x, double& y, double& z)
{
  read(x);
  read(y);
  read(z);
}

// write XYZ
void
FileMgr::writeXYZ(const double & x, const double& y, const double& z)
{
  write(x);
  write(y);
  write(z);
}

//-------------------------------------------------------
// write
// Write out characters from the supplied string.
//-------------------------------------------------------

void FileMgr::write(const string& s) const 
  {
  if(!p_) {
    IoError(filename_,"Error writing to unopened file").throwMe();
  }
  for (unsigned int i=0; i < s.length(); i++)
    {
    char c = s[i];
    int n = fwrite(&c,sizeof(char),1,p_);
    if (n != 1) IoError(filename_,"Write Error").throwMe();
    }
  }
//---------------------------------
// readInUnits, writeFloatInUnits  
// I/O routines that can handle non-base unit values
// and take the resultant performance hit on the chin
//----------------------------------

// ustr= Units to write in
void
FileMgr::writeFloatInUnits(const Uvar& uv, const string& ustr){
  float value= uv.getInUnits(ustr);
  write(value);
}

// ustr= Units to write in
void
FileMgr::writeDoubleInUnits(const Uvar& uv, const string& ustr){
  double value= uv.getInUnits(ustr);
  write(value);
}


//-------------------------------------------------------
// ustr = base units, non-auto unit modes will throw exception otherwise
// conversion_factor  is the conversion from input units to
// base units
//----------------------------------------------------------
void
FileMgr::readFloatInUnits(Uvar& uv, double conversion_factor,
const string& ustr){
  float value= 0;
  read(value);
  uv=Uvar(value*conversion_factor,ustr);
}

void
FileMgr::readDoubleInUnits(Uvar& uv, double conversion_factor,
const string& ustr){
  double  value= 0;
  read(value);
  uv=Uvar(value*conversion_factor,ustr);
}

void FileMgr::readLine(string& s) 
const 
{
  if(!p_) {
    IoError(filename_,"Error reading from unopened file").throwMe();
  }
  bool end_line_found =false;
  for (unsigned int i=0; i < s.length(); i++)
    {
    char c;
    int n = fread(&c,sizeof(char),1,p_);
    if (n != 1) IoError(filename_,"Read Line Error").throwMe();
    if(c=='\n')
      {
	end_line_found=true;
	break;
      }
    s[i] = c;
    if(eof()) break;
    }
  if(!end_line_found && !eof())
    {
      IoError(filename_,"String too short to read line").throwMe();
    }
}

//------------------------------------------------
// parse a string from an ascii file ignoring white space
//------------------------------------------------

string FileMgr::parseNextString() const 
  {
    if(!p_) {
      IoError(filename_,"Error reading from unopened file").throwMe();
    }
    string s="";
    bool string_complete=false;
    bool string_begun=false;
    while (!string_complete)
    {
      char c;
      int n = fread(&c,sizeof(char),1,p_);
      if (n != 1) IoError(filename_,"Read Error").throwMe();
      switch(c){
      case ' ':
      case '\n':
      case '\t':
	if (string_begun) string_complete=true;
	break;
      default :
	string_begun=true;
        s.append(1,c);
      }
    }
    return(s);
  }

//-------------------------------
// Supporting ordinary functions
//-------------------------------

//___________________________________________________________________
// i = chrbitget(const unsigned char c, unsigned int b1, unsigned int b2
//
// Extract the bits  (8-bit) from bit position b1..b2 
// inclusive (zero-offset) from the unsigned character c. 
//  

int chrbitget(const unsigned char c,const  int b1, const  int b2)
{
  unsigned char c1, c2;
  if (sizeof(c) > sizeof(unsigned char)) {
    ErrorMessage e("char too big for cbitbet");
    e.throwMe();
  }
  if (b2 < b1) {ErrorMessage("Bad position range in bitget ").throwMe();}
  if (b2 >7 || b1 < 0) {
    ErrorMessage e("Bad initial or final bit position");
    e.throwMe();
  }

  if (b2 <= 7)
  {
    c1 = c; 
    c1 >>=b2 +1;
    c1 <<= b2 +1;
  }
  c2 = c - c1;
  c2 >>= b1;
  int i = c2;
  return (i); 
  
}








