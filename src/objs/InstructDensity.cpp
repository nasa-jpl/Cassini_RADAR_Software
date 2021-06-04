//-------------------------------------------------------------------------------------------
// InstructDensity.cpp
//
// This file contains the instruction density methods.
// The InstructDensity class provides following processes:
//   (1) Set zero navigation errors:    InstructDensity::setNoNavErrors()
//   (2) Set navigation erros:          InstructDensity::setNavErrors()
//   (3) Calculate valid time:          InstructDensity::validTime()   
//   (4) Calculate partial derivatives for RWD and CSF with respect to beam pointing angles
//                                      InstructDensity::getPartialRWDandCSFvsBeamPointing()
//   (5) Calculate RWD changes due to spacecraft position errors
//                                      InstructDensity::PartialRWDvsShift() 
//   (6) Calculate the chanage of the (x,y,z) coordinates due the spacecraft position error
//       along flight path              InstructDensity::setLongitudinalShift()      
//   (7) Calculate the chanage of the (x,y,z) coordinates due the spacecraft altitude error
//                                      InstructDensity::setAltitudeShift()      
//   (8) Calculate the chanage of the (x,y,z) coordinates due the spacecraft sideslip error
//                                      InstructDensity::setHorizontalShift()   
//   (9) Calculate time-derivatives of RWD and CSF
//                                      InstructDensity::getDerivativeOfRWDandCSF()  
// 
// Reference: Memo yz_2003_06_27_valid_time
//------------------------------------------------------------------------------------------

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "InstructDensity.h"
#include "TargetGeom.h"
#include "Frame.h"
#include "Io.h"

using std::string;
// add by Gary on 07/08/2004:
using std::cout;
using std::cerr;
using std::endl;
using std::terminate;

//-------------------------
// InstructDensity Methods
//-------------------------
//-------------------------
// Constructors
//-------------------------
InstructDensity::InstructDensity()
{
  nav_error_set_ = false;

  // the following initials are not used in the updated code.
  N_min = 130;
  N_max = 480;

  T_min = Uvar(-20.0*60.0,"s");
  T_max = Uvar(20.0*60.0,"s");

  c0 =  1.592172323e-10;
  c1 = -2.8402088562e-7;
  c2 =  2.0006061533e-4;
  c3 = -0.067552944;
  c4 =  9.29002570;

  a3 = Uvar(4.0,"s");      
  a2 = 1.0 / Uvar(1.2040e5,"s"); 
  a1 = Uvar(1.1e2,"s");  
  a0 = Uvar(30.0,"s");      

}
//---------------------------------------------------------------------------
// This method defined the ZERO navigation errors for the valid time 
// calculation.
// 
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
void InstructDensity::setNoNavErrors()
{

  dh_ = Uvar(0.,"km");
  nav_dt_ = Uvar(0.,"s");
  ds_ = Uvar(0.,"km");
  dazi_ = Uvar(0.,"rad");
  dele_ = Uvar(0.,"rad");
  nav_error_set_ = true;
}

//---------------------------------------------------------------------------
// This method defined the navigation errors for the valid time calculation.
// 
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
void InstructDensity::setNavErrors()
{

  dh_ = Uvar(5.,"km");
  nav_dt_ = Uvar(5.,"s");
  ds_ = Uvar(5.,"km");
  dazi_ = Uvar(0.2*atan(1.)/45,"rad");
  dele_ = Uvar(0.2*atan(1.)/45,"rad");
  nav_error_set_ = true;
}

//---------------------------------------------------------------------------
// This method is to calculate the IEB valid time by using hybrid approach -
// First order analytical approach for SAR and ALT.
// Numerical approach for SCAT.
//
// Input parameters:
// (1) Radar mode                        m
// (2) Epoch time                        t
// (3) Number of pulse                   pul
// (4) Transmission/receiving offset     tro 
// (5) Chirp bandwidth                   cbw
// (6) Receiving frequency window size   rcv
// (7) Radar carrier frequency           freq
//
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
Uvar InstructDensity::validTime(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
//cout << "ghdebug: ID: in SAR mode" << endl;
//printf ("ghdebug:ID:in SAR mode");
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf, val_time;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
	{
	  for (unsigned int i_azi = 0; i_azi < 2; i_azi++)	    
	    {
	      
	      Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
	      Uvar ele = Uvar(0,"rad");
	      Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
	      DirectionVector look("bore",fbeam,t,0,0,1);
	      look.setAzimuthElevation(azi,ele);
	      TargetGeom tg(t);
	      tg.setTarget(target_name,ftitan);
	      
	      tg.setState(sc_state);
	      tg.setLookDirection(look); 
	      
	      minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);
	      
	    }
	}
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;
      
      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;
      
      // Calculate time-derivatives of RWD and CSF
      getDerivativeOfRWDandCSF(target_name, t, freq);
      
      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps) 
			- fabs(pRWD_pazi_) - fabs(pRWD_pele_));
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));
      
      //cout << "GHdebug: rwdmargin=" << rwd_margin ;
      //cout << "  fabs(pRWD_pt)=" << fabs(pRWD_pt) ;
      //cout << "  fabs(pRWD_ph)=" << fabs(pRWD_ph) ;
      //cout << "  fabs(pRWD_ps)="<< fabs(pRWD_ps)<< "  fabs(pRWD_pazi_)="<< fabs(pRWD_pazi_) << "  fabs(pRWD_pele_)" << fabs(pRWD_pele_) << endl;

      if(( m == 0 ) || ( m == 8 )) // SCAT mode
	{
	  // ignore negative margines
	  if(rwd_margin < 0.) rwd_margin = 0.;
	  if(csf_margin < 0.) csf_margin = 0.;
      
	  val_time = getValidTime(target_name, t, freq, rwd_margin, csf_margin);
	}
      else
	{
	  // Calculat valid times for RWD and CSF
	  
	  if(dRWD_dt_ != Uvar(0.,"s")/Uvar(1.,"s"))
	    val_time_rwd = rwd_margin/fabs(dRWD_dt_);
	  else
	    val_time_rwd = Uvar(0,"s");
	  
	  //cout<<csf_margin<<"  "<<dCSF_dt_<<endl;
	  if(dCSF_dt_ != Uvar(0.,"1/s")/Uvar(1.,"s"))
	    val_time_csf = csf_margin/fabs(dCSF_dt_);
	  else
	    val_time_csf = Uvar(0,"s");
	  //cout<<val_time_rwd<<"  "<<val_time_csf<<endl;

	  if(val_time_rwd < val_time_csf)
	    val_time = val_time_rwd;
	  else
	    val_time = val_time_csf;
	}
    }
  else
    {
      val_time= Uvar(0,"s");      
    }
 
    return(val_time);
}


// ******************** Added by Gary on 03/07/2007 to add Incidence Angle for Scat *****

//---------------------------------------------------------------------------
// This method is to calculate the IEB valid time by using hybrid approach -
// First order analytical approach for SAR and ALT.
// Numerical approach for SCAT.
//
// Input parameters:
// (1) Radar mode                        m
// (2) Epoch time                        t
// (3) Number of pulse                   pul
// (4) Transmission/receiving offset     tro
// (5) Chirp bandwidth                   cbw
// (6) Receiving frequency window size   rcv
// (7) Radar carrier frequency           freq
// (8) Max Incidence Angle for Scat mod  inc_angle
//
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
Uvar InstructDensity::validTime_inc(const string& target_name,
                                const unsigned int& m,
                                const Time& t,
                                const Uvar& pri,
                                const Uvar& pul,
                                const Uvar& tro,
                                const Uvar& cbw,
                                const Uvar& rcv,
                                const Uvar& freq,
				const Uvar& Inc_angle_bound)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;

      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3);

      Uvar inc_angle = tg.incidenceAngle();
      //Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");// GH 03/15/07.  Now passing this parameter in directly
      
      if(inc_angle>Inc_angle_bound)
        less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
//cout << "ghdebug: ID: in SAR mode" << endl;
//printf ("ghdebug:ID:in SAR mode");
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5;
    }
  else
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf, val_time;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);

      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);

      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);

      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);

      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
        {
          for (unsigned int i_azi = 0; i_azi < 2; i_azi++)
            {

              Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
              Uvar ele = Uvar(0,"rad");
              Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
              DirectionVector look("bore",fbeam,t,0,0,1);
              look.setAzimuthElevation(azi,ele);
              TargetGeom tg(t);
              tg.setTarget(target_name,ftitan);

              tg.setState(sc_state);
              tg.setLookDirection(look);

              minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);

            }
        }
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;

      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;

      // Calculate time-derivatives of RWD and CSF
      getDerivativeOfRWDandCSF(target_name, t, freq);

      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps)
                        - fabs(pRWD_pazi_) - fabs(pRWD_pele_));
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));

      //cout << "GHdebug: rwdmargin=" << rwd_margin ;
      //cout << "  fabs(pRWD_pt)=" << fabs(pRWD_pt) ;
      //cout << "  fabs(pRWD_ph)=" << fabs(pRWD_ph) ;
      //cout << "  fabs(pRWD_ps)="<< fabs(pRWD_ps)<< "  fabs(pRWD_pazi_)="<< fabs(pRWD_pazi_) << "  fabs(pRWD_pele_)" << fabs(pRWD_pele_) << endl;

      if(( m == 0 ) || ( m == 8 )) // SCAT mode
        {
          // ignore negative margines
          if(rwd_margin < 0.) rwd_margin = 0.;
          if(csf_margin < 0.) csf_margin = 0.;

          val_time = getValidTime(target_name, t, freq, rwd_margin, csf_margin);
        }
      else
        {
          // Calculat valid times for RWD and CSF

          if(dRWD_dt_ != Uvar(0.,"s")/Uvar(1.,"s"))
            val_time_rwd = rwd_margin/fabs(dRWD_dt_);
          else
            val_time_rwd = Uvar(0,"s");

          //cout<<csf_margin<<"  "<<dCSF_dt_<<endl;
          if(dCSF_dt_ != Uvar(0.,"1/s")/Uvar(1.,"s"))
            val_time_csf = csf_margin/fabs(dCSF_dt_);
          else
            val_time_csf = Uvar(0,"s");
          //cout<<val_time_rwd<<"  "<<val_time_csf<<endl;

          if(val_time_rwd < val_time_csf)
            val_time = val_time_rwd;
          else
            val_time = val_time_csf;
        }
    }
  else
    {
      val_time= Uvar(0,"s");
    }

    return(val_time);
}

// **************************************************************************************

//---------------------------------------------------------------------------
// This method is to calculate the IEB valid time by using numerical approach.
//
// Input parameters:
// (1) Radar mode                        m
// (2) Epoch time                        t
// (3) Number of pulse                   pul
// (4) Transmission/receiving offset     tro 
// (5) Chirp bandwidth                   cbw
// (6) Receiving frequency window size   rcv
// (7) Radar carrier frequency           freq
//
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
Uvar InstructDensity::validTime_num(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf, val_time;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
	{
	  for (unsigned int i_azi = 0; i_azi < 2; i_azi++)	    
	    {
	      
	      Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
	      Uvar ele = Uvar(0,"rad");
	      Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
	      DirectionVector look("bore",fbeam,t,0,0,1);
	      look.setAzimuthElevation(azi,ele);
	      TargetGeom tg(t);
	      tg.setTarget(target_name,ftitan);
	      
	      tg.setState(sc_state);
	      tg.setLookDirection(look); 
	      
	      minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);
	      
	    }
	}
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;
      
      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;
      
      // Calculate time-derivatives of RWD and CSF
      getDerivativeOfRWDandCSF(target_name, t, freq);
      
      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps) 
			- fabs(pRWD_pazi_) - fabs(pRWD_pele_));
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));

      // ignore negative margines
      if(rwd_margin < 0.) rwd_margin = 0.;
      if(csf_margin < 0.) csf_margin = 0.;
      
      val_time = getValidTime(target_name, t, freq, rwd_margin, csf_margin);
    }
  else
    {
      val_time = Uvar(0,"s");      
    }
  
  
    return(val_time);
}

//---------------------------------------------------------------------------
// This method is to calculate the IEB valid time by using the first order 
// analytical approach.
//
// Input parameters:
// (1) Radar mode                        m
// (2) Epoch time                        t
// (3) Number of pulse                   pul
// (4) Transmission/receiving offset     tro 
// (5) Chirp bandwidth                   cbw
// (6) Receiving frequency window size   rcv
// (7) Radar carrier frequency           freq
//
// Reference: Memo yz_2003_06_27_valid_time
//---------------------------------------------------------------------------
Uvar InstructDensity::validTime1(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
	{
	  for (unsigned int i_azi = 0; i_azi < 2; i_azi++)	    
	    {
	      
	      Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
	      Uvar ele = Uvar(0,"rad");
	      Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
	      DirectionVector look("bore",fbeam,t,0,0,1);
	      look.setAzimuthElevation(azi,ele);
	      TargetGeom tg(t);
	      tg.setTarget(target_name,ftitan);
	      
	      tg.setState(sc_state);
	      tg.setLookDirection(look); 
	      
	      minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);
	      
	    }
	}
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;
      
      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;
      
      // Calculate time-derivatives of RWD and CSF
      getDerivativeOfRWDandCSF(target_name, t, freq);
      
      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps) 
			- fabs(pRWD_pazi_) - fabs(pRWD_pele_));
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));
      
      // Calculat valid times for RWD and CSF
      
      if(dRWD_dt_ != Uvar(0.,"s")/Uvar(1.,"s"))
	val_time_rwd = rwd_margin/fabs(dRWD_dt_);
      else
	val_time_rwd = Uvar(0,"s");
	
      //cout<<csf_margin<<"  "<<dCSF_dt_<<endl;
      if(dCSF_dt_ != Uvar(0.,"1/s")/Uvar(1.,"s"))
	val_time_csf = csf_margin/fabs(dCSF_dt_);
      else
	val_time_csf = Uvar(0,"s");
      //cout<<val_time_rwd<<"  "<<val_time_csf<<endl;
      
      // The following commented code is to solve the valid time from tbe combination of
      // the first and second derivatives of RWD and CSF. Result shows that it is more 
      // sensitive to noise, thus we remain in the calculation by using the more stable 
      // first order approximation only.

      /*if(dRWD2_dt2_ == Uvar(0.,"s")/Uvar(1.,"s")/Uvar(1.,"s"))
	{
	  if(dRWD_dt_ == Uvar(0.,"s")/Uvar(1.,"s"))
	    //ErrorMessage("Both first and second order derivatives of RWD are zero.").throwMe();
	    val_time_rwd = 0;
	  else
	    val_time_rwd = fabs(rwd_margin/dRWD_dt_);	    
	}
      else
	{
	  Uvar val_time_rwd1 = 
               fabs((-dRWD_dt_ + sqrt(fabs(pow(dRWD_dt_,2.)+2*dRWD2_dt2_*rwd_margin)))/dRWD2_dt2_);
	  Uvar val_time_rwd2 = 
               fabs((-dRWD_dt_ - sqrt(fabs(pow(dRWD_dt_,2.)+2*dRWD2_dt2_*rwd_margin)))/dRWD2_dt2_);
	  val_time_rwd = val_time_rwd1;
	  if(val_time_rwd > val_time_rwd2)
	    val_time_rwd = val_time_rwd2;
	}
      
      if(dCSF2_dt2_ == Uvar(0.,"1/s")/Uvar(1.,"s")/Uvar(1.,"s"))
	{
	  if(dCSF_dt_ == Uvar(0.,"1/s")/Uvar(1.,"s"))
	    //ErrorMessage("Both first and second order derivatives of CSF are zero.").throwMe();
	    val_time_csf = 0;
	  else
	    val_time_csf = fabs(csf_margin/dCSF_dt_);	    
	}
      else
	{
	  Uvar val_time_csf1 = 
               fabs((-dCSF_dt_ + sqrt(fabs(pow(dCSF_dt_,2.)+2*dCSF2_dt2_*csf_margin)))/dCSF2_dt2_);
	  Uvar val_time_csf2 = 
               fabs((-dCSF_dt_ - sqrt(fabs(pow(dCSF_dt_,2.)+2*dCSF2_dt2_*csf_margin)))/dCSF2_dt2_);
	  val_time_csf = val_time_csf1;
	  if(val_time_csf > val_time_csf2)
	    val_time_csf = val_time_csf2;
	} */
    }
  else
    {
      val_time_rwd = Uvar(0,"s");
      val_time_csf = Uvar(0,"s");      
    }
  
  //cout<<val_time_rwd<<" "<<val_time_csf<<endl;

  if(val_time_rwd < val_time_csf)
    return(val_time_rwd);
  else
    return(val_time_csf);
}

//---------------------------------------------------------------------------
// This method is to calculate the IEB valid time by using the second order
// analytical approach
//---------------------------------------------------------------------------
Uvar InstructDensity::validTime2(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
	{
	  for (unsigned int i_azi = 0; i_azi < 2; i_azi++)	    
	    {
	      
	      Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
	      Uvar ele = Uvar(0,"rad");
	      Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
	      DirectionVector look("bore",fbeam,t,0,0,1);
	      look.setAzimuthElevation(azi,ele);
	      TargetGeom tg(t);
	      tg.setTarget(target_name,ftitan);
	      
	      tg.setState(sc_state);
	      tg.setLookDirection(look); 
	      
	      minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);
	      
	    }
	}
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;
      
      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;
      
      // Calculate time-derivatives of RWD and CSF
      getDerivativeOfRWDandCSF(target_name, t, freq);
      
      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps) 
			- fabs(pRWD_pazi_) - fabs(pRWD_pele_));
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));
      
      // Calculat valid times for RWD and CSF
      
      if(dRWD_dt_ != Uvar(0.,"s")/Uvar(1.,"s"))
	val_time_rwd = rwd_margin/fabs(dRWD_dt_);
      else
	val_time_rwd = Uvar(0,"s");
	
      //cout<<csf_margin<<"  "<<dCSF_dt_<<endl;
      if(dCSF_dt_ != Uvar(0.,"1/s")/Uvar(1.,"s"))
	val_time_csf = csf_margin/fabs(dCSF_dt_);
      else
	val_time_csf = Uvar(0,"s");
      //cout<<val_time_rwd<<"  "<<val_time_csf<<endl;
      
      // The following commented code is to solve the valid time from tbe combination of
      // the first and second derivatives of RWD and CSF. Result shows that it is more 
      // sensitive to noise, thus we remain in the calculation by using the more stable 
      // first order approximation only.

      if(dRWD2_dt2_ == Uvar(0.,"s")/Uvar(1.,"s")/Uvar(1.,"s"))
	{
	  if(dRWD_dt_ == Uvar(0.,"s")/Uvar(1.,"s"))
	    //ErrorMessage("Both first and second order derivatives of RWD are zero.").throwMe();
	    val_time_rwd = 0;
	  else
	    val_time_rwd = fabs(rwd_margin/dRWD_dt_);	    
	}
      else
	{
	  Uvar val_time_rwd1 = 
               fabs((-dRWD_dt_ + sqrt(fabs(pow(dRWD_dt_,2.)+2*dRWD2_dt2_*rwd_margin)))/dRWD2_dt2_);
	  Uvar val_time_rwd2 = 
               fabs((-dRWD_dt_ - sqrt(fabs(pow(dRWD_dt_,2.)+2*dRWD2_dt2_*rwd_margin)))/dRWD2_dt2_);
	  val_time_rwd = val_time_rwd1;
	  if(val_time_rwd > val_time_rwd2)
	    val_time_rwd = val_time_rwd2;
	}
      
      if(dCSF2_dt2_ == Uvar(0.,"1/s")/Uvar(1.,"s")/Uvar(1.,"s"))
	{
	  if(dCSF_dt_ == Uvar(0.,"1/s")/Uvar(1.,"s"))
	    //ErrorMessage("Both first and second order derivatives of CSF are zero.").throwMe();
	    val_time_csf = 0;
	  else
	    val_time_csf = fabs(csf_margin/dCSF_dt_);	    
	}
      else
	{
	  Uvar val_time_csf1 = 
               fabs((-dCSF_dt_ + sqrt(fabs(pow(dCSF_dt_,2.)+2*dCSF2_dt2_*csf_margin)))/dCSF2_dt2_);
	  Uvar val_time_csf2 = 
               fabs((-dCSF_dt_ - sqrt(fabs(pow(dCSF_dt_,2.)+2*dCSF2_dt2_*csf_margin)))/dCSF2_dt2_);
	  val_time_csf = val_time_csf1;
	  if(val_time_csf < val_time_csf2)
	    val_time_csf = val_time_csf2;
	} 
    }
  else
    {
      val_time_rwd = Uvar(0,"s");
      val_time_csf = Uvar(0,"s");      
    }
  
  //cout<<val_time_rwd<<" "<<val_time_csf<<endl;

  /*if(val_time_rwd < val_time_csf)
    return(val_time_rwd);
  else
  return(val_time_csf);*/

  return(dRWD_dt_);
}

//---------------------------------------------------------------------------
// Calculate RWD margin
//---------------------------------------------------------------------------
Uvar InstructDensity::getRWDmargin(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      
      // Modify margins by navigation and beam pointing errors
      rwd_margin = (rwd_margin - fabs(pRWD_pt) - fabs(pRWD_ph) - fabs(pRWD_ps) 
			- fabs(pRWD_pazi_) - fabs(pRWD_pele_));

      if(rwd_margin < 0) rwd_margin = 0.;

    }
  else
    {
      rwd_margin = Uvar(0,"s");    
    }
  
  return(rwd_margin);
}

//---------------------------------------------------------------------------
// Calculate CSF margin
//---------------------------------------------------------------------------
Uvar InstructDensity::getCSFmargin(const string& target_name,
				const unsigned int& m, 
                                const Time& t, 
                                const Uvar& pri, 
				const Uvar& pul, 
				const Uvar& tro,
				const Uvar& cbw, 
				const Uvar& rcv,
				const Uvar& freq)
{

  if(!nav_error_set_)  ErrorMessage("Navigation errors haven't been set yet.").throwMe();

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  Uvar rwd_margin, csf_margin;
  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;
  bool less_60deg;
  unsigned int Nb_s, Nb_e, Nb;

  less_60deg = true;
  if(( m == 0 ) || ( m == 8 )) // SCAT mode
    {  
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1;
 
      DirectionVector look_b3("bore",fbeam,t,0,0,1);
      TargetGeom tg(t);
      tg.setTarget(target_name,ftitan);
      tg.setState(sc_state);
      tg.setLookDirection(look_b3); 
      
      Uvar inc_angle = tg.incidenceAngle();
      Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
      if(inc_angle>Inc_angle_bound)
	less_60deg = false;
    }
  else if(( m == 1)||( m == 9)) // ALT mode
    {
      rwd_margin = 3*pri;
      Nb_s = 3;
      Nb_e = 3;
      Nb = 1; 
    }
  else if(( m == 2)||( m == 3)||( m == 10)||( m == 11)) // SAR mode
    {
      rwd_margin = 0.02*(pul*pri+tro);
      Nb_s = 1;
      Nb_e = 5;
      Nb = 5; 
    }
  else 
    ErrorMessage("no such radar mode!").throwMe();

  Uvar val_time_rwd, val_time_csf;
  if(less_60deg)
    {
      // RWD change due to the spacecraft altitude error
      setAltitudeShift(target_name, t);
      Uvar pRWD_ph = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the navigation time (position error along flight path) error
      setLongitudinalShift(target_name, t);
      Uvar pRWD_pt = PartialRWDvsShift(target_name, t);
      
      // RWD change due to the spacecraft horizontal position error
      setHorizontalShift(target_name, t);
      Uvar pRWD_ps = PartialRWDvsShift(target_name, t);
      
      // RWD and CSF change dut to azimuth pointing error
      getPartialRWDandCSFvsBeamPointing(target_name, t, freq);
      
      // Get CSF margin
      Uvec minmax_doppler("minmax_doppler",Nb*2);
      for (unsigned int n = Nb_s; n <= Nb_e; n++)
	{
	  for (unsigned int i_azi = 0; i_azi < 2; i_azi++)	    
	    {
	      
	      Uvar azi = Uvar((i_azi-0.5)*0.4*atan(1.)/45,"rad");
	      Uvar ele = Uvar(0,"rad");
	      Frame fbeam("CASSINI_RADAR_"+toStr(n),"Cassini");
	      DirectionVector look("bore",fbeam,t,0,0,1);
	      look.setAzimuthElevation(azi,ele);
	      TargetGeom tg(t);
	      tg.setTarget(target_name,ftitan);
	      
	      tg.setState(sc_state);
	      tg.setLookDirection(look); 
	      
	      minmax_doppler((n-Nb_s)*2+i_azi) = tg.doppler(lambda);
	      
	    }
	}
      Uvar dop_min, dop_max;
      min_max(dop_min,dop_max,minmax_doppler);
      Uvar dfc_all = (dop_min + dop_max)/2;
      
      csf_margin = ((rcv - cbw)/2. + dfc_all - dop_max)/2.;
      
      // Modify margins by navigation and beam pointing errors
      csf_margin = (csf_margin - fabs(pCSF_pazi_) - fabs(pCSF_pele_));

      if(csf_margin < 0) csf_margin = 0.;
      
    }
  else
    {
      csf_margin = Uvar(0,"1/s");    
    }
  
  return(csf_margin);
}

//--------------------------------------------------------------------------------
// This method is to calculate the changes of RWD and CSF due to the beam pointing
// errors (azimuth and elevation angles) 
//--------------------------------------------------------------------------------
void InstructDensity::getPartialRWDandCSFvsBeamPointing(const string& target_name, const Time& t, const Uvar& freq)
{

  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state_1("sc_state");
  ftitan.ephemeris(sc_state_1,"Cassini",t,"NONE");

  DirectionVector look_1("bore",fbeam,t,0,0,1);
  TargetGeom tg_1(t);
  tg_1.setTarget(target_name,ftitan);
  tg_1.setState(sc_state_1);
  tg_1.setLookDirection(look_1); 
  
  Uvar DOP_1 = tg_1.doppler(lambda);
  Uvar RWD_1 = 2*tg_1.range()/speed_light;
  
  // Change azimuth angle
  DirectionVector look_2("bore",fbeam,t,0,0,1);
  look_2.setAzimuthElevation(dazi_,Uvar(0.0,"rad"));
  TargetGeom tg_2(t);
  tg_2.setTarget(target_name,ftitan);
  tg_2.setState(sc_state_1);
  tg_2.setLookDirection(look_2);

  Uvar DOP_2 = tg_2.doppler(lambda);
  Uvar RWD_2 = 2*tg_2.range()/speed_light;

  pCSF_pazi_ = -(DOP_2 - DOP_1);
  pRWD_pazi_ = RWD_2 - RWD_1;

  // Change elevation angle
  DirectionVector look_3("bore",fbeam,t,0,0,1);
  look_3.setAzimuthElevation(Uvar(0.0,"rad"),dele_);
  TargetGeom tg_3(t);
  tg_3.setTarget(target_name,ftitan);
  tg_3.setState(sc_state_1);
  tg_3.setLookDirection(look_3);

  Uvar DOP_3 = tg_3.doppler(lambda);
  Uvar RWD_3 = 2*tg_3.range()/speed_light;

  pCSF_pele_ = -(DOP_3 - DOP_1);
  pRWD_pele_ = RWD_3 - RWD_1;

  //cout<<pRWD_pele_<<" "<<pCSF_pele_<<endl;
  //cout<<pRWD_pazi_<<" "<<pCSF_pazi_<<endl;

}

//---------------------------------------------------------------------------------
// This method is to calculate the RWD changes due to the spacecraft position error
// dX, dY, dZ, in the Titan fix coordinates.
//---------------------------------------------------------------------------------
Uvar InstructDensity::PartialRWDvsShift(const string& target_name, const Time& t)
{

  Uvar speed_light = Uvar("speed_light");

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state_1("sc_state");
  ftitan.ephemeris(sc_state_1,"Cassini",t,"NONE");

  DirectionVector look_1("bore",fbeam,t,0,0,1);
  TargetGeom tg_1(t);
  tg_1.setTarget(target_name,ftitan);
  tg_1.setState(sc_state_1);
  tg_1.setLookDirection(look_1); 
  
  Uvar RWD_1 = 2*tg_1.range()/speed_light;
  
  Uvar Xsc = sc_state_1.position()[PositionVector::X] + dX_;
  Uvar Ysc = sc_state_1.position()[PositionVector::Y] + dY_;
  Uvar Zsc = sc_state_1.position()[PositionVector::Z] + dZ_;
  
  Uvar Vxsc = sc_state_1.velocity()[FloatVector::X];
  Uvar Vysc = sc_state_1.velocity()[FloatVector::Y];
  Uvar Vzsc = sc_state_1.velocity()[FloatVector::Z];

  PositionVector p("p",ftitan, t, Xsc, Ysc, Zsc);
  FloatVector v("v",ftitan, t, Vxsc, Vysc, Vzsc);
  StateVector sc_state_2("sci",p,v);
 
  DirectionVector look_2("bore",fbeam,t,0,0,1);
  TargetGeom tg_2(t);
  tg_2.setTarget(target_name,ftitan);
  tg_2.setState(sc_state_2);
  tg_2.setLookDirection(look_2);

  
  Uvar RWD_2 = 2*tg_2.range()/speed_light;

  //cout<<RWD_2<<" "<<RWD_1<<endl;

  return(RWD_2 - RWD_1);

}

//---------------------------------------------------------------------------------
// This method is used to calculate the changes of the (x,y,z) coordinates of the 
// spacecraft due to the position error along the flight path.
//
// The parameter passed in this method is 
// (1) Epoch time     t
//---------------------------------------------------------------------------------
void InstructDensity::setLongitudinalShift(const string& target_name, const Time& t)
{
  Frame ftitan("IAU_"+target_name,target_name);
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");

  Uvar Vxsc = sc_state.velocity()[FloatVector::X];
  Uvar Vysc = sc_state.velocity()[FloatVector::Y];
  Uvar Vzsc = sc_state.velocity()[FloatVector::Z];

  Uvar Vsc = sqrt(pow(Vxsc,2.0)+pow(Vysc,2.0)+pow(Vzsc,2.0));
  dX_ = Vxsc*nav_dt_;
  dY_ = Vysc*nav_dt_;
  dZ_ = Vzsc*nav_dt_;
  //cout<<sqrt(pow(dX_,2.0)+pow(dY_,2.0)+pow(dZ_,2.0))<<endl;
}

//---------------------------------------------------------------------------------
// This method is used to calculate the changes of the (x,y,z) coordinates of the 
// spacecraft due to the sideslip error.
//
// The parameter passed in this method is 
// (1) Epoch time     t
//---------------------------------------------------------------------------------
void InstructDensity::setHorizontalShift(const string& target_name, const Time& t)
{
  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");
  DirectionVector look("bore",fbeam,t,0,0,1);
  TargetGeom tg(t);
  tg.setTarget(target_name,ftitan);
  tg.setState(sc_state);
  tg.setLookDirection(look); 

  DirectionVector Vsc = sc_state.velocity();
  double vx = Vsc[DirectionVector::X];
  double vy = Vsc[DirectionVector::Y];
  double vz = Vsc[DirectionVector::Z];

  DirectionVector P_n = tg.nadir();
  double px = P_n[DirectionVector::X];
  double py = P_n[DirectionVector::Y];
  double pz = P_n[DirectionVector::Z];

  double sx = vy*pz - vz*py;
  double sy = vz*px - vx*pz;
  double sz = vx*py - vy*px;

  double s = sqrt(pow(sx,2.0)+pow(sy,2.0)+pow(sz,2.0));

  dX_ = sx/s*ds_;
  dY_ = sy/s*ds_;
  dZ_ = sz/s*ds_;

}

//---------------------------------------------------------------------------------
// This method is used to calculate the changes of the (x,y,z) coordinates of the 
// spacecraft due to the altitude error.
//
// The parameter passed in this method is 
// (1) Epoch time     t
//---------------------------------------------------------------------------------
void InstructDensity::setAltitudeShift(const string& target_name, const Time& t)
{
  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");
  DirectionVector look("bore",fbeam,t,0,0,1);
  TargetGeom tg(t);
  tg.setTarget(target_name,ftitan);
  tg.setState(sc_state);
  tg.setLookDirection(look); 

  DirectionVector Vsc = sc_state.velocity();
  double vx = Vsc[DirectionVector::X];
  double vy = Vsc[DirectionVector::Y];
  double vz = Vsc[DirectionVector::Z];

  DirectionVector P_n = tg.nadir();
  double px = P_n[DirectionVector::X];
  double py = P_n[DirectionVector::Y];
  double pz = P_n[DirectionVector::Z];

  double sx = vy*pz - vz*py;
  double sy = vz*px - vx*pz;
  double sz = vx*py - vy*px;

  double hx = vy*sz - vz*sy;
  double hy = vz*sx - vx*sz;
  double hz = vx*sy - vy*sx;  

  double h = sqrt(pow(hx,2.0)+pow(hy,2.0)+pow(hz,2.0));

  dX_ = hx/h*dh_;
  dY_ = hy/h*dh_;
  dZ_ = hz/h*dh_;

}
//---------------------------------------------------------------------------------
// This method is to calculate the first and second time-derivatives of RWD and CSF
// The parameters passed in this method are 
//
// (1) Epoch time          t
// (2) Carrier frequency   freq
//---------------------------------------------------------------------------------
void InstructDensity::getDerivativeOfRWDandCSF(const string& target_name, const Time& t, const Uvar& freq)
{

  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  // RWD and CSF at time t
  DirectionVector look0("bore",fbeam,t,0,0,1);
  TargetGeom tg0(t);
  tg0.setTarget(target_name,ftitan);
  tg0.setState(sc_state);
  tg0.setLookDirection(look0); 
  Uvar RWD0 = 2*tg0.range()/speed_light;
  Uvar dop_freq0 = tg0.doppler(lambda);

  Uvar dt = Uvar(1,"s");

  // RWD and CSF at time t - dt
  ftitan.ephemeris(sc_state,"Cassini",t-dt,"NONE");
  DirectionVector look1("bore",fbeam,t-dt,0,0,1);
  TargetGeom tg1(t-dt);
  tg1.setTarget(target_name,ftitan);
  tg1.setState(sc_state);
  tg1.setLookDirection(look1); 
  Uvar RWD1 = 2*tg1.range()/speed_light;
  Uvar dop_freq1 = tg1.doppler(lambda);

  // RWD and CSF at time t + dt
  ftitan.ephemeris(sc_state,"Cassini",t+dt,"NONE");
  DirectionVector look2("bore",fbeam,t+dt,0,0,1);
  TargetGeom tg2(t+dt);
  tg2.setTarget(target_name,ftitan);
  tg2.setState(sc_state);
  tg2.setLookDirection(look2); 
  Uvar RWD2 = 2*tg2.range()/speed_light;
  Uvar dop_freq2 = tg2.doppler(lambda);

  double key = 1.;
  if((RWD1==Uvar(0.0,"s"))||(RWD2==Uvar(0.0,"s")))
    key = 0.;

  // The first order time-derivatives
  dRWD_dt_ = key*(RWD2 - RWD1)/(2*dt);
  dCSF_dt_ = -key*(dop_freq2 - dop_freq1)/(2*dt);

  // The second order time-derivatives
  dRWD2_dt2_ = key*(RWD2 + RWD1 -2.0*RWD0)/dt/dt;
  dCSF2_dt2_ = -key*(dop_freq2 + dop_freq1 - 2.0*dop_freq0)/dt/dt;

}

//---------------------------------------------------------------------------------
// calculate valid time base on numerical searching
//---------------------------------------------------------------------------------
Uvar InstructDensity::getValidTime(const string& target_name, const Time& t, const Uvar& freq,
     const Uvar& RWDmargin, const Uvar& CSFmargin)
{

  Uvar speed_light = Uvar("speed_light");
  Uvar lambda=speed_light/freq;

  Frame ftitan("IAU_"+target_name,target_name);
  Frame fbeam("CASSINI_RADAR_"+toStr(3),"Cassini");
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");


  // RWD and CSF at time t
  DirectionVector look0("bore",fbeam,t,0,0,1);
  TargetGeom tg0(t);
  tg0.setTarget(target_name,ftitan);
  tg0.setState(sc_state);
  tg0.setLookDirection(look0); 
  Uvar RWD0 = 2*tg0.range()/speed_light;
  Uvar dop_freq0 = tg0.doppler(lambda);

  Uvar dt = Uvar(1,"s");

  // forward search
  Uvar dRWD = RWD0 - RWD0;
  Uvar dCSF = dop_freq0 - dop_freq0;
  Uvar vt_rwd_f = Uvar(0.,"s");
  Uvar vt_rwd_b = Uvar(0.,"s");
  Uvar vt_csf_f = Uvar(0.,"s");
  Uvar vt_csf_b = Uvar(0.,"s");
  Uvar valid_time, valid_time_rwd, valid_time_csf;
  Uvar Inc_angle_bound = Uvar(60.*atan(1.)/45.,"rad");
  
  if(tg0.foundSurfaceIntercept())
    {
      // forward search 
      while(dRWD < RWDmargin)
	{
	  vt_rwd_f += dt;  
	  ftitan.ephemeris(sc_state,"Cassini",t + vt_rwd_f,"NONE");
	  DirectionVector look1("bore",fbeam,t + vt_rwd_f,0,0,1);
	  TargetGeom tg1(t + vt_rwd_f);
	  tg1.setTarget(target_name,ftitan);
	  tg1.setState(sc_state);
	  tg1.setLookDirection(look1);
 
	  Uvar RWD1 = 2*tg1.range()/speed_light; 

	  Uvar inc_angle = tg1.incidenceAngle();
	  if(inc_angle>Inc_angle_bound)
	    RWD1 = 0;

	  dRWD = fabs(RWD1 - RWD0);
	}
      //cout<<"vt_rwd_f: "<<vt_rwd_f<<endl;

      // backward search 
      Uvar dRWD = RWD0 - RWD0;
      while(dRWD < RWDmargin)
	{
	  vt_rwd_b -= dt;  
	  ftitan.ephemeris(sc_state,"Cassini",t + vt_rwd_b,"NONE");
	  DirectionVector look1("bore",fbeam,t + vt_rwd_b,0,0,1);
	  TargetGeom tg1(t + vt_rwd_b);
	  tg1.setTarget(target_name,ftitan);
	  tg1.setState(sc_state);
	  tg1.setLookDirection(look1); 

	  Uvar RWD1 = 2*tg1.range()/speed_light; 

	  Uvar inc_angle = tg1.incidenceAngle();
	  if(inc_angle>Inc_angle_bound)
	    RWD1 = 0;

	  dRWD = fabs(RWD1 - RWD0);
	}
      //cout<<"vt_rwd_b: "<<vt_rwd_b<<endl;

      valid_time_rwd = fabs(vt_rwd_f);
      if(valid_time_rwd > fabs(vt_rwd_b))
	valid_time_rwd = fabs(vt_rwd_b);

      //------------------------------------------      
      // valid time for CSF
      //------------------------------------------
      // forward search 
      while((dCSF < CSFmargin)&&(vt_csf_f < valid_time_rwd))
	{
	  vt_csf_f += dt;  
	  ftitan.ephemeris(sc_state,"Cassini",t + vt_csf_f,"NONE");
	  DirectionVector look1("bore",fbeam,t + vt_csf_f,0,0,1);
	  TargetGeom tg1(t + vt_csf_f);
	  tg1.setTarget(target_name,ftitan);
	  tg1.setState(sc_state);
	  tg1.setLookDirection(look1); 

	  Uvar dop_freq1 = tg1.doppler(lambda);

	  Uvar inc_angle = tg1.incidenceAngle();
	  if(inc_angle>Inc_angle_bound)
	    dop_freq1 = 0;

	  dCSF = fabs(dop_freq1 - dop_freq0);
	}
      //cout<<"vt_csf_f: "<<vt_csf_f<<endl;

      // backward search 
      Uvar dCSF = dop_freq0 - dop_freq0;
      while((dCSF < CSFmargin)&&(fabs(vt_csf_b)<valid_time_rwd))
	{
	  vt_csf_b -= dt;  
	  ftitan.ephemeris(sc_state,"Cassini",t + vt_csf_b,"NONE");
	  DirectionVector look1("bore",fbeam,t + vt_csf_b,0,0,1);
	  TargetGeom tg1(t + vt_csf_b);
	  tg1.setTarget(target_name,ftitan);
	  tg1.setState(sc_state);
	  tg1.setLookDirection(look1); 

	  Uvar dop_freq1 = tg1.doppler(lambda);

	  Uvar inc_angle = tg1.incidenceAngle();
	  if(inc_angle>Inc_angle_bound)
	    dop_freq1 = 0;

	  dCSF = fabs(dop_freq1 - dop_freq0);
	}
      //cout<<"vt_csf_b: "<<vt_csf_b<<endl;

      valid_time_csf = fabs(vt_csf_f);
      if(valid_time_csf > fabs(vt_csf_b))
	valid_time_csf = fabs(vt_csf_b);
   
      // compare with valid time for RWD and CSF, choose the minimum
      valid_time = valid_time_rwd;
      if(valid_time > valid_time_csf)
	valid_time = valid_time_csf;

    }
  else
    {
      valid_time = 0;
      //cout<<"RWD0: "<<RWD0<<endl;
      //cout<<"DOP0: "<<dop_freq0<<endl;
    }
  
  return(valid_time);
}

//----------------------------------------------------------------------------
// This method is to set the number of command for initialization
// This method was developed for the SAR mode of T7 flyby only, and will not
// be used in the updated class.
//----------------------------------------------------------------------------
void InstructDensity::setNCommand(const unsigned int& Nc)
{

  if(( Nc >= N_min ) && ( Nc <= N_max ))
    Nc_ = Nc;
  else 
    ErrorMessage("Number is commands is out of range!").throwMe();
   
}

//------------------------------------------------------------------------------
// This method is to calculate the valid time for the SAR mode of T7 flyby only.
// It will not be used in the updated class.
//------------------------------------------------------------------------------
Uvar InstructDensity::computeValTime(Uvar& t)
{
  Uvar SarModeDuration = Uvar(960,"s");
  if(( t < T_min ) || ( t > T_max ))
    ErrorMessage("The time is out of range!").throwMe();

  // Calculate scaling coefficient based on Nc
  double c = c0*pow(2*Nc_,4.) + c1*pow(2*Nc_,3.) + c2*pow(2*Nc_,2.) + c3*(2*Nc_) + c4;

  // Calculate valid time without offset
  Uvar y1 = a0*exp(-pow(t/a1,2.)) + a2*pow(t,2.);

  // Add offset to y1 scaled by c
  Uvar valid_time;
  if(Nc_ <= 240)
    valid_time = c * y1 + a3;
  else
    valid_time = SarModeDuration/Nc_;
  
  return(valid_time);
}
