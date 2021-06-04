//-------------------------------------------------------------------
// InstructDensity.h
//
// This file contains the instruction density methods.
// The InstructDensity class provides following processes:
//   (1) Set number of commands:    InstructDensity::setDensity()
//   (2) Calculate valid time:      InstructDensity::ComputeValTime()                  

#ifndef InstructDensity_H
#define InstructDensity_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class Baq;
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"

using std::string;


//-------------------------
// Class InstructDensity declaration
//-------------------------
class InstructDensity
  {
  public:
    InstructDensity();

    void setNCommand(const unsigned int& Nc);
    Uvar computeValTime(Uvar& t);
    void setNoNavErrors();
    void setNavErrors();
    Uvar validTime(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);
    Uvar validTime_inc(const string& target_name,
                   const unsigned int& m,
                   const Time& t,
                   const Uvar& pri,
                   const Uvar& pul,
                   const Uvar& tro,
                   const Uvar& cbw,
                   const Uvar& rcv,
                   const Uvar& freq,
                   const Uvar& Inc_angle_bound);
    Uvar validTime_num(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);
    Uvar validTime1(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);
    Uvar validTime2(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);
    Uvar getRWDmargin(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);
    Uvar getCSFmargin(const string& target_name,
		   const unsigned int& m, 
		   const Time& t, 
		   const Uvar& pri,
		   const Uvar& pul, 
		   const Uvar& tro,
                   const Uvar& cbw, 
		   const Uvar& rcv, 
		   const Uvar& freq);

  private:
    Uvar T_min, T_max;
    unsigned int N_min, N_max, Nc_;
    double c0, c1, c2, c3, c4;
    Uvar a0, a1, a2, a3;
    Uvar dX_, dY_, dZ_;
    Uvar dh_, nav_dt_, ds_, dazi_, dele_;
    Uvar pRWD_pS_,pRWD_pazi_,pCSF_pazi_,pRWD_pele_,pCSF_pele_;
    Uvar dRWD_dt_, dCSF_dt_, dRWD2_dt2_, dCSF2_dt2_;
    bool nav_error_set_;

    void setLongitudinalShift(const string& target_name, const Time& t);
    void setHorizontalShift(const string& target_name, const Time& t);
    void setAltitudeShift(const string& target_name, const Time& t);
    Uvar PartialRWDvsShift(const string& target_name, const Time& t);
    void getPartialRWDandCSFvsBeamPointing(const string& target_name, const Time& t, const Uvar& freq);
    void getDerivativeOfRWDandCSF(const string& target_name, const Time& t, const Uvar& freq);
    Uvar getValidTime(const string& target_name, const Time& t, const Uvar& freq, const Uvar& RWDmargin, const Uvar& CSFmargin);
  };
#endif

















































