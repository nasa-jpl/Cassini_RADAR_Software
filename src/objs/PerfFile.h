//-----------------------------------------------
//PerfFile.h
// This file contains the PerfFile class declaration.
// This class provides an interface to handle performance calculation data
//  read/write.
// This class takes an  input of ambiguityFNN(four nearest neighbor) class
// object and writes thermal snr and amb ratio for both single and multilook
// using the Titan body fixed frame. 
// 
//
// class PerfFile;
//
//
//   Construction:
//  PerfFile(const string& filename,const string& filetype) 
//
//-----------------------
//supporting functions
//----------------------
// void load(const AmbiguityFNN& amb);
//
//--------------------------
//I/O
//--------------------------
//void readRecord(AmbiguityFNN& amb)
//void writeRecord() throw(ErrorMessage);
//--------------------------------------------------------



#ifndef PerfFile_H
#define PerfFile_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class PerfFile;
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "Ambiguity.h"
#include "Io.h"

using std::string;


//-------------------------
// Class Ieb declaration
//-------------------------
class PerfFile
  {
  public: 
    
    
    //--------------
    // Construction
    //-------------- 
    PerfFile(const string& filename,const string& mode);
    ~PerfFile();//destructor
    
    //-----------------------
    //supporting functions
    //----------------------
    void close();
    
    //--------------------------
    //I/O
    //--------------------------
    void readRecord(AmbiguityFNN& amb) ;
    void writeRecord(const AmbiguityFNN& amb);
  
    bool eof();
   
   

  private:
    //---------------------------
    //internal representation
    //---------------------------
    bool file_opened_;
    FileMgr file_;
    string filetype_;
  
  };



#endif






