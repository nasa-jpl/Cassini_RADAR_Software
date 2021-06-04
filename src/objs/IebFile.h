//-----------------------------------------------
//IebFile.h
// This file contains the IebFile class declaration.
// The IebFile class provides an interface to handle IEB record file such as
// writing and reading
//
// This class will be shared with CASSINI_RADAR_RMSS group (GH). 
//
//
// class IebFile;
//
//
//   Construction:
//  IebFile(const string& filename,const string& filetype) throw(ErrorMessage)
//
//-----------------------
//supporting functions
//----------------------
// void load(const Ieb& ieb);

//--------------------------
//I/O
//--------------------------
//void readRecord(const Time& t, Ieb& ieb) throw(ErrorMessage);
//void writeRecord() throw(ErrorMessage);
//--------------------------------------------------------



#ifndef IebFile_H
#define IebFile_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class IebFile;
#include <vector>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Ieb.h"
#include "Io.h"

using std::string;


//-------------------------
// Class Ieb declaration
//-------------------------
class IebFile
  {
  public: 
    
    //---------------------
    //typedef's and enums
    //---------------------
    typedef std::vector<unsigned short> sdata;//short integer data
    
    //--------------
    // Construction
    //-------------- 
    IebFile(const string& filename,const string& mode);
    ~IebFile();//destructor
    
    

    //-----------------------
    //supporting functions
    //----------------------
    void load(const Ieb& ieb);
    
    //--------------------------
    //I/O
    //--------------------------
    void readAllIebRecords();
    void writeAllIebRecords() throw(ErrorMessage);
    Ieb getIeb(const Time& t) ;

    //--------------------------------
    //adding and deleting ieb record at time t
    //-------------------------------
    void addIeb(const Ieb& ieb);
    void deleteIeb(const Ieb& ieb);

    //------------------------
    //return the size of ieb records, start time, and end time
    //------------------------
    Time getStartTime();
    Time getEndTime();
    unsigned int getNumberofRecords();

    //------------------
    //size of slow and field fields:constant fixed
    //-----------------
    static const unsigned int Nslowfield = 10;
    static const unsigned int Nfastfield = 7;
    
  private:
    //---------------------------
    //internal representation
    //---------------------------
    bool slowfastfield_loaded_;
    bool file_opened_;
    bool ieb_all_written_;
    bool ieb_all_read_;
    //---------------------
    //slow and fast fields
    //---------------------
    sdata  slowfield_;   
    sdata  fastfield_;   

    
    FileMgr file_;
    string filetype_;   
    vector<Ieb> ieb_list_;  
  };


#endif






