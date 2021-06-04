//--------------------------------------------------------------------
// Ckernel.h
//
// This files contains the interface class and functions to generate
// quaternions based on attitude and time information
//
//void LoadIvd() 
//  Load Ivd data such as number of data, time, attitude information
//  the other for minus z-axis.  UTC time will be stored as Ephemeris time
//void GenerateQuats() 
//  From ephemeris time and attitude of x and mz, generate quaternions
//void WriteData() 
//  Write data into a ckernel file
//----------------------------------------------------------------------

#ifndef Ckernel_H
#define Ckernel_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "Array.h"
#include "Ivd.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Frame.h"
#include "Time.h"
#include "Constants.h"
#include "SpiceUsr.h"
#include "Ivd.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

class Ckernel
  {
  // Static interface
  public:

  static void validTimeRange(const string& filename, const string& sc_name,
    Time& start_time, Time& end_time);

  // Object interface
 
  public:

  // Constructors  
  Ckernel();

  // I/O
  void WriteData(string& ckernel_file) ;

  // Data handling
  void LoadIvd(const Ivd& ivd) ;
  void GenerateQuats() ;

  //------------------------
  //Ivd public data fields
  //------------------------

  Dmat quats,avvs;
  Dmat attitude_x,attitude_mz;
  Uvec et_x, et_mz;
  Dvec sclkdp_x,sclkdp_mz;  

  Time start_x,end_x,start_mz,end_mz;
  string head_x,base_x,head_mz,base_mz; 
  unsigned int N_dat;
 
  //------------------------
  // Internal representation
  //------------------------

  private:

  const static unsigned int maxrec = 20000;  
  bool ivd_data_loaded_;
  bool quats_data_loaded_;
 
  SpiceInt INST_;
  string REF_;
  string CENTERNAME_;
  string SCNAME_;
  string SEGID_;
  SpiceBoolean AVFLAG_;  
  };
  
#endif





