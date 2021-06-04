//--------------------------------------------------------------------
//--------------------------------------------------------------------
// IebProfile.h
//  The tool contains the interface calss and functions for reading
//  IEB Division  parameters from config file and providing division
//  parameters at a requested time.  
//  
//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef IebProfile_H
#define IebProfile_H

#include <string>
#include <map>
#include "Array.h"
#include "Units.h"
#include "Error.h"
#include "Constants.h"



class IebProfile
  {
  public:

    // Constructors  
    IebProfile(Config& cfg);

 
    //setup
    void config(Config& cfg);

    //--------------
    // Predicates
    //--------------
    bool validDivision(const string& suffix) const;
    
    //------------------------
    // Special service methods
    //------------------------
    void divisionCheck(const string& suffix) const;
  
    //---------------------------------------------------
    //Obtain sar or scat prf value at the request time t
    //---------------------------------------------------
    Uvar getSarPrf(const Uvar& epoch_relative_time);
    Uvar getScatPrf(const Uvar& epoch_relative_time);

    //--------------------------------------------------------------
    // Apply duty cycle limit set by cfg keywords dutycycle_limit_*
    //--------------------------------------------------------------

    double limitDutyCycle(const double& dutycycle, const Uvar& prf);

    //----------------------------
    //get Time information 
    //------------------------------  
    Time getTriggerTime();
    Uvar getFirstDivisionStartTime();
    Uvar getLastDivisionEndTime();

    //---------------------
    //getND and RL
    //----------------------
    Uvar getND();
    Uvar getMaxND();
    Uvar getRL();
    Uvar getMaxRL();

    //-----------------------------------------
    //IEB Division supportin functions
    //-------------------------------------------
    bool foundDivision(const Uvar& epoch_relative_time) const;
    string getDivisionName(const Uvar& epoch_relative_time) const;
    string getDivisionMode(const string& div_name) const;
    Uvar getDivisionTimestep(const string& div_name) const;
    unsigned int getDivisionBem(const string& div_name) const;
    unsigned int getDivisionBaq(const string& div_name) const;
    unsigned int getDivisionCsr(const string& div_name) const;
    double  getDivisionNoiseBitSetting(const string& div_name) const;
    //Uvar  getDivisionBpd(const string& div_name) const;
    double getDivisionDutycycle(const string& div_name) const;
    Uvar getDivisionPrf(const string& div_name) const ;
    unsigned int getDivisionNpulse(const string& div_name) const;
    unsigned int getDivisionN_bursts_in_flight(const string& div_name) const;
    unsigned int getDivisionPercentBW(const string& div_name) const;
    string getDivisionAutorad(const string& div_name) const;
    Uvar getDivisionRip(const string& div_name) const;
    double getDivisionDatarate(const string& div_name) const;
    int getDivisionTro(const string& div_name) const;    
    
  private:

    
    //cfg filename
    string cfg_filename_;

    //ieb trigger time
    Time trigger_time_;

    //--------------------------------
    //parameter used for constructing divisions
    //--------------------------------
    mutable Config::numbers div_start_time_map_;
    mutable Config::numbers div_end_time_map_;
    mutable Config::numbers div_time_step_map_;
    mutable Config::strings div_mode_map_;
    mutable Config::numbers div_bem_map_;
    mutable Config::numbers div_baq_map_;
    mutable Config::numbers div_csr_map_;
    mutable Config::numbers div_dutycycle_map_;
    mutable Config::numbers div_noise_bit_setting_map_;
    mutable Config::numbers div_bpd_map_;
    mutable Config::numbers div_prf_map_;
    mutable Config::numbers div_percentBW_map_;
    mutable Config::numbers div_pulse_map_;
    mutable Config::numbers div_bursts_in_flight_map_;
    mutable Config::strings div_autorad_map_;
    mutable Config::numbers div_rip_map_;
    mutable Config::numbers div_rate_map_;
    mutable Config::numbers div_tro_map_; 
    //parameters needed for polynomial fit of sar prf
    Dvec sar_prf_time_polyfit_;
    // parameters needed for dutycycle limits
    Uvec dutycycle_limit_prf_;
    Dvec dutycycle_limit_d_;
    // parameters storing the Scat PRF profile
    string target_name_;
    Uvec scat_prf_;
    Uvec scat_thetai_;
    
    //ND and RL
    Uvar ND_;
    Uvar RL_;
    Uvar Max_ND_;
    Uvar Max_RL_;
  };

  
  
#endif





