//=============================================================================
// RadiometerData.h
//
// This file contains the RadiometerData class declarations.
// The RadiometerData class provides a convenient package to hold vectors
// of radiometer data taken from a L1A file.
//
// Interface summary:
//
// RadiometerData methods:
//
//   Construction:
//
//   RadiometerData(L1B_filename); Construct from data in L1A file
//   RadiometerData(L1B);          Construct from data in L1A file
//
// RadiometerData public data:
//
// Radiometer data are stored in vectors of UnitVar's which are publicly
// accessible.
//
//=============================================================================

#ifndef RADIOMETERDATA_H
#define RADIOMETERDATA_H

#include <string>
#include <vector>
#include "Array.h"
#include "Units.h"
#include "Error.h"
#include "L1B.h"

using std::vector;
using std::string;

class RadiometerData
  {
  public:

//  typedef UnitVar<double> Uvar;
//  typedef Array1D<Uvar> Uvec;

  //--------------
  // Constructors
  //--------------

  RadiometerData(L1B& l1b,Time t1=Time(), Time t2=Time()) 
    throw(ErrorMessage);

  
  //--------------------
  // Public data vectors
  //--------------------

  // Raw data
  Uvec cnt_rl,cnt_nd,cnt_radio,rad;
  Uvec csr,r_mode,ctbe;
  Array1D<unsigned int> sclk;
  Uvec brst; 
  Uvec hip,cip,rip,bpd;

  Uvec pul,tro,rwd;
  Uvec pri;

  Uvec be1tmp,be2tmp,be3tmp,be4tmp,be5tmp;
  Uvec wgb1t1,wgb3t1,wgb3t2,wgb3t3,wgb5t1;
  Uvec mratmp,mruttm,nsdtmp,rlotmp;

  // Processed data
  unsigned int record_count;
  Time epoch;
  Array1D<Time> time;
  Ivec beam_number;
  Uvec ncnt_rl, ncnt_nd, ncnt_radio;

  // Values which should be const static but cannot be
  Uvar misc_delay;
  Uvar delta_tau;
  int offset;
  
  
  private:

  
  //-------------------
  // Internal variables
  //-------------------

  //------------------
  // Private Methods
  //------------------

  void processTime();
  void extractBeamNumbers();
  void normalize();

  };

#endif





