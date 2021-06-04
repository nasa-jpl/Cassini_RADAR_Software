//--------------------------------------------------------------------
//--------------------------------------------------------------------
//Ivd.h
//
//This files contains the interface class and functions for reading two
// Ivd files to generate quaternions that will be used to generate 
// a ckernel file with Data Type 3 (see ck.req in spice documents)
//
//In Data Type 3, we do not calculate the angular velocity vector and
// set all the values to 0.0.
//
//Two Ivd files contain attitude information of two axis (x and -z)
// as a function of UTC time. 
//
//
//
//
//void ReadDataFiles() throw(ErrorMessage)
//  Read in time and attitude information from two IVD files, 
//  one for x-axis and the other for minus z-axis.  
//  UTC time will be stored as Ephemeris time

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef Ivd_H
#define Ivd_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "Array.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Frame.h"
#include "Time.h"
#include "Constants.h"
#include "SpiceUsr.h"
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;



class Ivd
  {
  public:

  // Constructors  
  Ivd();

  // Predicates
  // no Predicates

  // I/O
  // no I/O

  // Data handling
  void ReadDataFile(ifstream& ivd_file) throw(ErrorMessage);
  void GenerateQuats() throw(ErrorMessage);
 
  // Other methods
  // no Other methods

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

  bool ivd_data_loaded;
 
  private:

  const static unsigned int maxrec_;  

  
 
  string SCNAME_;
  //------------------------
  // Internal representation
  //------------------------
  };

  
  
#endif





