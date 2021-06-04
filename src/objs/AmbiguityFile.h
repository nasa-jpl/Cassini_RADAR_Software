//-----------------------------------------------
//AmbiguityFile.h
// This file contains the AmbiguityFile class declaration.
// The AmbiguityFile class provides an interface to handle 
// file containing ambiguity geometry calculation.
//
//
// 
//
//
// class AmbiguityFile;
//
//
//   Construction:
//  AmbiguityFile(const string& filename,const string& filetype) 
//
//-----------------------
//supporting functions
//----------------------
// void load(const Ambiguity& amb);
//
//--------------------------
//I/O
//--------------------------
//void readRecord(Ambiguity& amb) 
//void writeRecord() 
//--------------------------------------------------------



#ifndef AmbiguityFile_H
#define AmbiguityFile_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class AmbiguityFile;
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Ambiguity.h"
#include "Io.h"

using std::string;


//-------------------------
// Class Ieb declaration
//-------------------------
class AmbiguityFile
  {
  public: 
    
   
    //--------------
    // Construction
    //-------------- 
    AmbiguityFile(const string& filename,const string& mode);
    ~AmbiguityFile();
    
    //-------------
    //{get,set} method
    //------------


    //--------------------
    //Clear method
    //-------------------
    

    //--------------------
    //Supporting functions
    //--------------------
    void load(const Ambiguity& amb);

    //-----------------------------
    //I/O
    //----------------------------
    void readRecord(Ambiguity& amb);
    void writeRecord();
    void close();
    bool eof();


  private:
    //---------------------------
    //internal representation
    //---------------------------
    bool ambiguity_data_loaded_;
    bool file_opened_;

    FileMgr file_;
    string filetype_;
   
    Ambiguity amb_;
  };


#endif






