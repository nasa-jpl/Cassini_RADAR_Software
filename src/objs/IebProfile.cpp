//-------------------------------------------
//IebProfile.cpp
//  This file contains method definitions for IebProfile class.
//  IebProfile provides support for constructing IebList by conating
// all division parameters and providing time-dependent  outputs
//-------------------------------------------

#include <string>
#include <map>
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "IebProfile.h"
#include "Plot.h"
#include "Utils.h"
#include "TargetGeom.h"
#include "Frame.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;




//----------------------------
//  Methods for Class IebProfile   
//----------------------------

//--------------
// Constructos: there is only one method
//--------------
IebProfile::IebProfile(Config& cfg)
  :cfg_filename_(cfg.filename()),
   sar_prf_time_polyfit_("sar prf polyfit coeff"),
   dutycycle_limit_prf_("duty_cycle_limit_prf"),
   dutycycle_limit_d_("duty_cycle_limit_d"),
   scat_prf_("scat_prf_"),
   scat_thetai_("scat_thetai_")
//old:: :cfg_filename_(cfg.filename()),trigger_time_(cfg.str("IEB_Trigger_time")),sar_prf_time_polyfit_("sar prf polyfit coeff")
  {
    config(cfg);
    trigger_time_.setSclk (cfg.str("spacecraft"), cfg.getInt("IEB_Trigger") ) ;
    trigger_time_ = trigger_time_ + cfg.num("IDAPT_delta");
    // Add difference between the command in the PEF going to CDS and RADAR execution
    trigger_time_ = trigger_time_ + cfg.num("CDS_Cmd_Delay");
  }

//--------------------------
//Setup: config method-> Read all division parameters
//----------------------------
void IebProfile::config(Config& cfg)
  {
    //ND and RL
    ND_ = cfg["ND_Rip"];
    RL_ =cfg["RL_Rip"];
    Max_ND_ = cfg["Max_ND_Rip"];
    Max_RL_ = cfg["Max_RL_Rip"];

    if( (Max_ND_ - ND_) < Uvar(0.999*mstos,"s"))
      {
	cout<<"max nd and nd "<< Max_ND_ <<" "<<ND_<<endl;
	ErrorMessage("IebProfile.cpp::config method: Max ND should be at least 1 ms larger than normal ND").throwMe(); 
      }
    if( (Max_RL_ - RL_)<Uvar(4.999*mstos,"s")) 
      {
	cout<<"max rl and rl "<< Max_RL_<<" "<<RL_<<endl;
	ErrorMessage("IebProfile.cpp::config method: max RL should be at least 5 ms larger than normal RL value ").throwMe();
      }
    //---------------------------------
    //Collect time information to construct IEB division
    //-----------------------------------
    cfg.suffixSet("div_start_time_",div_start_time_map_);
    cfg.suffixSet("div_end_time_",div_end_time_map_);
    cfg.suffixSet("div_time_step_",div_time_step_map_);
    cfg.suffixSet("div_mode_",div_mode_map_);
    cfg.suffixSet("div_bem_",div_bem_map_);
    cfg.suffixSet("div_baq_",div_baq_map_);
    cfg.suffixSet("div_csr_",div_csr_map_);
    cfg.suffixSet("div_dutycycle_",div_dutycycle_map_);
    cfg.suffixSet("div_noise_bit_setting_",div_noise_bit_setting_map_);
    cfg.suffixSet("div_prf_",div_prf_map_);
    cfg.suffixSet("div_number_of_pulses_",div_pulse_map_);
    cfg.suffixSet("div_n_bursts_in_flight_",div_bursts_in_flight_map_);
    cfg.suffixSet("div_percent_of_BW_",div_percentBW_map_);
    cfg.suffixSet("div_auto_rad_",div_autorad_map_);
    cfg.suffixSet("div_max_data_rate_",div_rate_map_);
    cfg.suffixSet("div_rip_",div_rip_map_);
    cfg.suffixSet("div_tro_",div_tro_map_);
    //debugging purpose
    //make sure start time is always ealier than end time
    Config::numbers::const_iterator pstart=div_start_time_map_.begin();
    Config::numbers::const_iterator pend=div_end_time_map_.begin();
    Config::numbers::const_iterator pstep=div_time_step_map_.begin();
    
    //check start time is always earler than end time inside the same div.
    //if start time is later than end time inside the same div., throw error
    for( pstart=div_start_time_map_.begin();
	 pstart != div_start_time_map_.end()
	   && pend !=div_end_time_map_.end();
	 ++pstart,++pend)
      {
	if(pstart->second > pend->second) ErrorMessage("Flyby.cpp:Config(): Inside a divsion, start time is later thane end time").throwMe();
      }
    //check division time step > 0 s
    for(pstep = div_time_step_map_.begin(); pstep !=div_time_step_map_.end();
	++pstep)
      {
	if(pstep ->second < Uvar(0,"s")) ErrorMessage("Flyby.cpp::division time step should be larger than 0 s").throwMe();
      }
    
    pstart=div_start_time_map_.begin();
    pend=div_end_time_map_.begin();
    pstep=div_time_step_map_.begin();
    //check map size: all should have the same size! otherwise, throw error
    if(div_start_time_map_.size() != div_end_time_map_.size()
       || div_start_time_map_.size() != div_mode_map_.size()
       ||div_start_time_map_.size() != div_time_step_map_.size()
       || div_start_time_map_.size() !=div_bem_map_.size()
       || div_start_time_map_.size() != div_baq_map_.size()
       || div_start_time_map_.size() !=div_csr_map_.size()
       || div_start_time_map_.size() !=div_dutycycle_map_.size()
       || div_start_time_map_.size() != div_prf_map_.size()
       || div_start_time_map_.size() != div_pulse_map_.size()
       || div_start_time_map_.size() != div_bursts_in_flight_map_.size()
       || div_start_time_map_.size() != div_percentBW_map_.size()
       || div_start_time_map_.size() !=div_autorad_map_.size()
       || div_start_time_map_.size()!=div_rip_map_.size()
       || div_start_time_map_.size()!=div_rate_map_.size()
       || div_start_time_map_.size() != div_noise_bit_setting_map_.size()
       || div_start_time_map_.size() != div_tro_map_.size())
      {
	cout<<"start  "<<div_start_time_map_.size()<<endl;
	cout<<" end   "	    << div_end_time_map_.size()<<endl;
	cout<<" time step "   <<div_time_step_map_.size()<<endl;
	cout<<" mode " 	    <<div_mode_map_.size()<<endl;
	cout<<" bem "	    <<div_bem_map_.size()<<endl;
	cout<<" baq " 	    <<div_baq_map_.size()<<endl;
	cout<<" csr "	    <<div_csr_map_.size()<<endl;
	cout<<" dutycycle "	    <<div_dutycycle_map_.size()<<endl;
	cout<<" prf " 	    <<div_prf_map_.size()<<endl;
	cout<<" pulse "   <<div_pulse_map_.size()<<endl;
	cout<<" bursts in flight " <<div_bursts_in_flight_map_.size()<<endl;
	cout<<" percentBW "  <<div_percentBW_map_.size()<<endl;
	cout<<" autorad "   <<div_autorad_map_.size()<<endl;
	cout<<" rip "    <<div_rip_map_.size()<<endl;
	cout<<" rate"   <<div_rate_map_.size()<<endl;
	cout<<" noise bit setting "  <<div_noise_bit_setting_map_.size()<<endl;
        cout<<" div tro map "<< div_tro_map_.size()<<endl;
	ErrorMessage("IebProfile.cpp: Division size mismatch").throwMe();
      }
    
    
    //debugging
    //print all the division parameters on the screen
    /*
      Config::strings::const_iterator pmode=div_mode_map_.begin();
      Config::numbers::const_iterator pbem=div_bem_map_.begin();
      Config::numbers::const_iterator  pbaq=div_baq_map_.begin();
      Config::numbers::const_iterator  pcsr=div_csr_map_.begin();
      Config::numbers::const_iterator  pdutycycle=div_dutycycle_map_.begin();
      Config::numbers::const_iterator pprf=div_prf_map_.begin();
      Config::numbers::const_iterator  ppulse=div_pulse_map_.begin();
      Config::strings::const_iterator  pautorad=div_autorad_map_.begin();
      Config::numbers::const_iterator  prip=div_rip_map_.begin();
      Config::numbers::const_iterator  prate=div_rate_map_.begin(); 
      for( pstart=div_start_time_map_.begin();
      pstart != div_start_time_map_.end()
      && pend != div_end_time_map_.end()
      && pmode !=div_mode_map_.end()
      && pbem !=div_bem_map_.end()
      && pbaq !=div_baq_map_.end()
      && pcsr !=div_csr_map_.end()
      && pdutycycle !=div_dutycycle_map_.end()
      && pprf !=div_prf_map_.end()
      && ppulse != div_pulse_map_.end()
      && pautorad !=div_autorad_map_.end()
      && prip!=div_rip_map_.end()
      && prate != div_rate_map_.end();
      ++pstart,++pend,++pmode,++pbem,++pbaq,++pcsr,
      ++pdutycycle,++pprf,++ppulse,++pautorad,++prip,++prate)
      
      {
      cout<<pstart->first<<" "<<pstart->second<<endl;
      cout<<pend->first<<" "<<pend->second<<endl;
      cout<<pmode->first<<" "<<pmode->second<<endl;
      cout<<pbem->first<<" "<<pbem->second<<endl;
      cout<<pbaq->first<<" "<<pbaq->second<<endl;
      cout<<pcsr->first<<" "<<pcsr->second<<endl;
      cout<<pdutycycle->first<<" "<<pdutycycle->second<<endl;
      cout<<pprf->first<<" "<<pprf->second<<endl;
      cout<<ppulse->first<<" "<<ppulse->second<<endl;
      cout<<pautorad->first<<" "<<pautorad->second<<endl;
      cout<<prip->first<<" "<<prip->second<<endl;
      cout<<prate->first<<" "<<prate->second<<endl;
      cout<<endl;
      }  
    */

   
    Uvar constant_prf_offset=Uvar(0,"Hz");
    if(cfg.keywordExists("constant_prf_offset")){
      cout<<"config has constant prf offset "<<endl;
      constant_prf_offset=cfg["constant_prf_offset"];
      cout<<"prf prfile will have a constant offset of "<<constant_prf_offset.getInUnits("Hz")<<" Hz with respect to ideal profile "<<endl;
    }
    
    //--------------------------------
    //Use polynomial to fit SAR prf(KHz) vs time(min) 
    //----------------------------------
    double units_in_min= 60.0;
    double units_in_thousand=1000.0;
    
    Uvar read_keyword;
    Uvec prf("prf values"),time("time prf");
    Dvec prf_value_in_KHz("prf in khz"),time_in_min("time in min");
    string prf_selection_option = cfg.str("prf_profile_option");
    if (prf_selection_option =="constant_prf")
      {
	Uvar const_prf = cfg["constant_prf_value"];
	sar_prf_time_polyfit_.resize(1);
	sar_prf_time_polyfit_(0) = const_prf.getInUnits("Hz")/1000.0;
	cout<<"constant prf will be used "<<endl;
      }
    else if (prf_selection_option =="polynomial_fit")
      {
	read_keyword=cfg["number_of_time_prf_variations"];
	unsigned int N_time_prf = (unsigned int) (read_keyword.getInUnits(""));
	time.resize(2*N_time_prf);
	prf.resize(2*N_time_prf); 
	time_in_min.resize(2*N_time_prf);
	prf_value_in_KHz.resize(2*N_time_prf);
	
	if (N_time_prf < 2)
	  {
	    ErrorMessage("For poly fit, number_of_altitude_prf_variations should be larger than 1").throwMe();
	  }
	for (unsigned int i = 0; i < N_time_prf;++i)
	  {
	    
	    if(cfg.keywordExists("time_prf01")){//new format
	      time(N_time_prf-i-1) = cfg["time_prf"+toStr(i+1,2)];
	      time_in_min(N_time_prf-i-1) =time(N_time_prf-i-1).getInUnits("s")/units_in_min;
	      
	      prf(N_time_prf-i-1) = cfg["prf"+toStr(i+1,2)];
	    }
	    else{//old format
	      time(N_time_prf-i-1) = cfg["time_prf"+toStr(i+1)];
	      time_in_min(N_time_prf-i-1) =time(N_time_prf-i-1).getInUnits("s")/units_in_min;
	      
	      prf(N_time_prf-i-1) = cfg["prf"+toStr(i+1)];
	    }

	    prf(N_time_prf-i-1) = prf(N_time_prf-i-1)+constant_prf_offset;
	    if( prf(N_time_prf-i-1)<Uvar(0,"Hz")) ErrorMessage("constant prf offset is too high so that some of adjusted PRF values are belows 0").throwMe();
	    prf_value_in_KHz(N_time_prf-i-1) = prf(N_time_prf-i-1).getInUnits("Hz")/units_in_thousand;
	  }
	for (unsigned int i = N_time_prf; i < 2*N_time_prf;++i)
	  {
	    time(i) = -time(2*N_time_prf-i-1);
	    time_in_min(i)=time(i).getInUnits("s")/units_in_min;
	    prf(i)= prf(2*N_time_prf-i-1);
	    prf_value_in_KHz(i)= prf(i).getInUnits("Hz")/units_in_thousand;
	  }
	//debug
	//cout<<"time "<<time<<endl;
	//cout<<"prf "<<prf_value<<endl;
	//cout<<"time in min "<<time_in_min<<endl;
	//cout<<"prf in KHz "<<prf_value_in_KHz<<endl;
	
	read_keyword = cfg["time_prf_polynomial_power"];
	unsigned int prf_poly_power = (unsigned int) read_keyword.getInUnits("");
	prf_poly_power = prf_poly_power+1;  
	sar_prf_time_polyfit_.resize(prf_poly_power);
	sar_prf_time_polyfit_ = 0.0;  
	if (prf_poly_power > N_time_prf)
	  {
	    ErrorMessage("Reduce number of polynomial power to fit the data").throwMe();
	  }
	//remove unit
	polyfit(sar_prf_time_polyfit_,prf_value_in_KHz,time_in_min); 
	cout<<"Finish finding polyfit coefficents for prf profile"<<endl;
      }
    else
      {
	ErrorMessage("Unknown prf selection option was chosen").throwMe();
      }
    //debugging: if want to show the fit, turn on plot.show("x")
    // plot prf vs time  when polyfit option is chosen
    if (prf_selection_option =="polynomial_fit")
      {
	Plot plot;
	plot.addXY(time,"min",prf,"Hz",line("none"),sym("circle","red",2));
	Uvec prf_fit("prf fit",prf.size());
	for (unsigned int i = 0; i < time.size();++i)
	  {
	    double time_in_min= time(i).getInUnits("s")/units_in_min;
	    double prf_in_KHz= 0;
	    for (unsigned int j = 0; j < sar_prf_time_polyfit_.size();++j)
	      {
		prf_in_KHz += sar_prf_time_polyfit_(j)*pow(time_in_min,j);
	      }
	    prf_fit(i) = Uvar(prf_in_KHz * 1000.0,"Hz");
	  } 
	plot.addXY(time,"min",prf_fit,"Hz",line("solid","black",2),sym("none"));
	plot.setTitle("prf vs alt fit ");
	//GH: need to remove this plot popping up
	//plot.show("x");
      }

    //-----------------------------------------------------------------
    // Dutycycle limit setup
    // These arrays store a set of dutycycle limits that can be applied
    // to a user supplied dutycycle given the corresponding PRI.
    // See the method limitDutyCycle() in this class.
    //-----------------------------------------------------------------

    int N_dutycycle_limit =
      convert_to_int(cfg["dutycycle_limit_number_points"]);
    dutycycle_limit_prf_.resize(N_dutycycle_limit);
    dutycycle_limit_d_.resize(N_dutycycle_limit);
    for (int i=0; i < N_dutycycle_limit; ++i)
      {
      dutycycle_limit_prf_(i) = cfg["dutycycle_limit_prf_" + toStr(i+1)];
      dutycycle_limit_d_(i) =
        convert_to_double(cfg["dutycycle_limit_d_" + toStr(i+1)]);
      }

    //-----------------------------------------------------------------
    // Scat PRF array setup.
    // These arrays store a set of incidence angle,prf pairs that specify
    // a PRF profile for scatterometer divisions that set div PRF to zero.
    // See the method getScatPrf() in this class.
    //-----------------------------------------------------------------

    target_name_ = cfg.str("target");
    int N_scat_prf =
      convert_to_int(cfg["scat_prf_number_points"]);
    scat_prf_.resize(N_scat_prf);
    scat_thetai_.resize(N_scat_prf);
    for (int i=0; i < N_scat_prf; ++i)
      {
      scat_prf_(i) = cfg["scat_prf_" + toStr(i+1)];
      scat_thetai_(i) = cfg["scat_thetai_" + toStr(i+1)];
      }

    //---------------------------
    //dutycycle polyfit
    //this part will be turned off on the assumption that
    // we are going to use a constant dutycycle for sar
    //-------------------------------
    /*
    string dutycycle_selection_option = cfg.str("dutycycle_profile_option");
    Uvec altitude_dutycycle("altitude for dutycycle");
    Dvec dutycycle_value("duty cycle values"); 
    Uvec dutycycle_plot("for plot purpose");
    if (dutycycle_selection_option =="constant_dutycycle")
    {	
    sar_dutycycle_polyfit_.resize(1);
    sar_dutycycle_polyfit_(0)= cfg["constant_dutycycle_value"].getInUnits("");
    //cout<<"constant dutycycle will be used "<<endl;
    }
    else if (dutycycle_selection_option =="polynomial_fit")
    {
    read_keyword=cfg["number_of_altitude_dutycycle_variations"];
    unsigned int N_alt_dutycycle = (unsigned int)read_keyword.getInUnits("");
    altitude_dutycycle.resize(N_alt_dutycycle);
    dutycycle_value.resize(N_alt_dutycycle); 
    Dvec altitude_dutycycle_in_thousandkm("altitude in km",N_alt_dutycycle);
    dutycycle_plot.resize(N_alt_dutycycle);
    if (N_alt_dutycycle < 2) 
    {
    ErrorMessage("For poly fit, number_of_altitude_dutycycle_variations should be larger than 1").throwMe();
    }
    for (unsigned int i = 0; i < N_alt_dutycycle;++i)
    {
    altitude_dutycycle(i) = cfg["alt_dutycycle"+toStr(i+1)];
    altitude_dutycycle_in_thousandkm(i) = altitude_dutycycle(i).getInUnits("km")/units_in_thousand;
    dutycycle_value(i) = cfg["dutycycle"+toStr(i+1)].getInUnits("");
    dutycycle_plot(i) = cfg["dutycycle"+toStr(i+1)];
    }
    for (unsigned int i = 0; i <N_alt_dutycycle-1; i++)
    {
    if (altitude_dutycycle(i) > altitude_dutycycle(i+1))
    {
    ErrorMessage("Reorder altitude from low altitude to high altitude").throwMe();
    }    
    }
    read_keyword = cfg["alt_dutycycle_polynomial_power"];
    unsigned int dutycycle_poly_power = int(read_keyword.getInUnits(""));
    dutycycle_poly_power = dutycycle_poly_power+1;  
    sar_dutycycle_polyfit_.resize(dutycycle_poly_power);         
    sar_dutycycle_polyfit_ = 0.0;  
    if (dutycycle_poly_power > N_alt_dutycycle)
    {
    throw ErrorMessage("Reduce number of polynomial power to fit the data");
    }
    //remove unit
    polyfit(sar_dutycycle_polyfit_,dutycycle_value,altitude_dutycycle_in_thousandkm);   
    cout<<"Finished finding polyfit coefficient for duty cycle profile"<<endl;
    }
    else
    {
    ErrorMessage("Unknown dutycycle selection option was chosen").throwMe();
    }
    
    //debugging purpose
    //To see the dutycycle profile, turn on plot.show("x")
    // plot alt vs prf when polyfit option is chosen
    if (dutycycle_selection_option =="polynomial_fit")
    {
    Plot plot;
    plot.addXY(altitude_dutycycle,"km",dutycycle_plot," ",line("none"),sym("circle","red",2));
    Uvec dutycycle_fit("dutycycle fit",dutycycle_plot.size());
    for (unsigned int i = 0; i < altitude_dutycycle.size();++i)
    {
    double alt_in_1000km = altitude_dutycycle(i).getInUnits("km")/1000.0;
    double dutycycle= 0;
    for (unsigned int j = 0; j < sar_dutycycle_polyfit_.size();++j)
    {
    dutycycle += sar_dutycycle_polyfit_(j)*pow(alt_in_1000km,j);
    }
    dutycycle_fit(i) = dutycycle;
    }
    plot.addXY(altitude_dutycycle,"km",dutycycle_fit," ",line("solid","black",2),sym("none"));
    plot.setTitle("dutycycle vs alt fit ");
    //plot.show("x");
    }
    */
    
  }//end of config setup method

//------------------------------
//Predicates
//---------------------------------
bool IebProfile::validDivision(const string& suffix) const
 {
   Config::numbers::const_iterator p = div_start_time_map_.find(suffix);
   return( p !=div_start_time_map_.end());
 }

//--------------------------------
//Check there is any valid division
//----------------------------------
void IebProfile::divisionCheck(const string& suffix) const
  {
    if(!validDivision(suffix))
      {
	ErrorMessage e("Division suffix (" + suffix +
		       ") not present in " +cfg_filename_);
	e.throwMe();
      }
  }

//-----------------------------------------
// get  time information
//------------------------------------------
Time IebProfile::getTriggerTime()
  {
    return(trigger_time_);
  }

Uvar IebProfile::getFirstDivisionStartTime()
  {
    Config::numbers::const_iterator pstart=div_start_time_map_.begin();
    Uvar start_time= pstart->second;
   
    for( pstart=div_start_time_map_.begin();
	 pstart != div_start_time_map_.end();
	 ++pstart)
      {
	if(start_time > pstart->second) start_time = pstart->second;
      }
    return(start_time);
  }

Uvar IebProfile::getLastDivisionEndTime()
  {
    Config::numbers::const_iterator pend=div_end_time_map_.begin();
    Uvar end_time = pend->second;
    for( pend=div_end_time_map_.begin();
	 pend != div_end_time_map_.end();
	 ++pend)
      {
	if(end_time < pend->second) end_time = pend->second;
      }
    return(end_time);
  }
//-----------------------------------
//getND
//------------------------------------
Uvar IebProfile::getND()
  {
    return(ND_);
  }
Uvar IebProfile::getMaxND()
  {
    return(Max_ND_);
  }
//-----------------------------
//getRL
//--------------------------------
Uvar IebProfile::getRL()
  {
    return(RL_);
  }
Uvar IebProfile::getMaxRL()
  {
    return(Max_RL_);
  }

//----------------------------------------------
//Obtain sar prf value at the request time t
//------------------------------------------------
Uvar IebProfile::getSarPrf(const Uvar& epoch_relative_time)
  {
    if(epoch_relative_time < Uvar(-60*60,"s") 
       || epoch_relative_time > Uvar(60*60,"s"))
     cout << "Warning: IebProfile.cpp::getSarPrf: outside sar time, from -30 to 30 min" << endl;
	// ErrorMessage("IebProfile.cpp::getSarPrf: outside sar time, from -30 to 30 min").throwMe();

    if(epoch_relative_time < Uvar(-30*60,"s") 
       || epoch_relative_time > Uvar(30*60,"s")){
      cout<<"#------------------------------------#"<<endl;
      cout<<"#  You are trying to run SAR outside #"<<endl;
      cout<<"# +-30 min w.r.t. the closest approach time #"<<endl;
      cout<<"#  ----------------------------------#"<<endl;
      cout<<endl;
    }
      

    double time_in_min= epoch_relative_time.getInUnits("s")/60.0;
    double prf_in_KHz= 0;
    for (unsigned int j = 0; j < sar_prf_time_polyfit_.size();++j)
      {
	prf_in_KHz += sar_prf_time_polyfit_(j)*pow(time_in_min,j);
      }
    return( Uvar(prf_in_KHz * 1000.0,"Hz"));
  }

//---------------------------------------------------------------------
// Apply limits listed in the cfg file to the input dutycycle value
// using the input PRF.
//-----------------------------------------------------------------------

double  IebProfile::limitDutyCycle(const double& dutycycle, const Uvar& prf)
  {
  double r_duty = dutycycle;
  int N = dutycycle_limit_prf_.size();
  for (int i=0; i < N-1; ++i)
    {  // locate requested prf within prf limit array and apply check
    if (prf >= dutycycle_limit_prf_(i) && prf < dutycycle_limit_prf_(i+1))
      {
      if (dutycycle > dutycycle_limit_d_(i)) r_duty = dutycycle_limit_d_(i);
      return(r_duty);
      }
    }
  // requested prf is at or above the final limit, apply last check
  if (dutycycle > dutycycle_limit_d_(N-1)) r_duty = dutycycle_limit_d_(N-1);
  return(r_duty);
  }

//---------------------------------------------------------------------
// Apply limits listed in the cfg file to the input dutycycle value
// using the input PRF.
//-----------------------------------------------------------------------

Uvar IebProfile::getScatPrf(const Uvar& t)
  {
  Uvar r_prf;
  // Do geometry calculations to get incidence angle
  StateVector sc_state;
  Frame target_frame("IAU_" + target_name_, target_name_);
  target_frame.ephemeris(sc_state,"Cassini",t,"NONE");
  Frame fbeam("CASSINI_RADAR_3","Cassini");
  DirectionVector boresight("boresight",fbeam,t,0,0,1);
  TargetGeom tg(t);
  tg.setState(sc_state);
  tg.setLookDirection(boresight);
  tg.setTarget(target_name_,target_frame);
  Uvar beam_inc_angle = tg.incidenceAngle();
  int N = scat_prf_.size();
  if (tg.foundSurfaceIntercept())
    {
    for (int i=0; i < N-1; ++i)
      {  // locate requested thetai within thetai limit array 
      if (beam_inc_angle >= scat_thetai_(i) &&
          beam_inc_angle < scat_thetai_(i+1))
        {
        r_prf = scat_prf_(i);
        return(r_prf);
        }
      }
    }
  // requested thetai is at or above the final thetai, so return final prf
  // Or, no intercept was found, so return the final prf
  r_prf = scat_prf_(N-1);
  return(r_prf);
  }

//------------------------------------
//IEB division supporting functions
//---------------------------------------
bool IebProfile::foundDivision(const Uvar& epoch_relative_time) const
  {
    bool found = false;
    Config::numbers::const_iterator pend=div_end_time_map_.begin();  
    for(Config::numbers::const_iterator pstart=div_start_time_map_.begin();
	pstart     != div_start_time_map_.end() 
	  && pend  != div_end_time_map_.end();
	  ++pstart,++pend)
      {
	if(pstart->second<= epoch_relative_time 
	   && epoch_relative_time < pend->second)
	  {
	    found = true;
	    break;
	  }
      }
    return(found);
  }

//----------------------------------
//IEB division supporting functions
//-----------------------------------
string IebProfile::getDivisionName(const Uvar& epoch_relative_time) const
  { 
    Config::numbers::const_iterator pend=div_end_time_map_.begin();
    Config::strings::const_iterator pmode=div_mode_map_.begin();
    if(!foundDivision(epoch_relative_time)) 
      ErrorMessage("IebProfile.cpp::No division").throwMe();
    string return_mode="";//default
    for(Config::numbers::const_iterator pstart=div_start_time_map_.begin();
	pstart     != div_start_time_map_.end() 
	  && pend  != div_end_time_map_.end()
	  && pmode !=div_mode_map_.end()
	  ;++pstart,++pend,++pmode)
      {
	if(pstart->second<= epoch_relative_time 
	   && epoch_relative_time < pend->second)
	  {
	    //cout<<pstart->second<<" "<<pend->second<<endl;
	    return_mode=pstart->first;
	    break;
	  }
      }
    return(return_mode);
  }

string IebProfile::getDivisionMode(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_mode_map_[div_name]);
  }
Uvar IebProfile::getDivisionTimestep(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_time_step_map_[div_name]);
  }
unsigned int IebProfile::getDivisionBem(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int) round_double(div_bem_map_[div_name].getInUnits("")));
  }
unsigned int IebProfile::getDivisionBaq(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int) round_double(div_baq_map_[div_name].getInUnits("")));
  }
unsigned int IebProfile::getDivisionCsr(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int)round_double(div_csr_map_[div_name].getInUnits("")));
  }
double  IebProfile::getDivisionNoiseBitSetting(const string& div_name) const
  {
    divisionCheck(div_name);
    return( div_noise_bit_setting_map_[div_name].getInUnits(""));
  }
//Uvar  IebProfile::getDivisionBpd(const string& div_name) const
//  {
//    divisionCheck(div_name);
//    return(div_bpd_map_[div_name]);
//  }
double IebProfile::getDivisionDutycycle(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_dutycycle_map_[div_name].getInUnits(""));
  }
Uvar IebProfile::getDivisionPrf(const string& div_name) const 
  {
    divisionCheck(div_name);
    return(div_prf_map_[div_name]);
  }
unsigned int IebProfile::getDivisionNpulse(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int)round_double(div_pulse_map_[div_name].getInUnits("")));
  }
unsigned int IebProfile::getDivisionN_bursts_in_flight(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int)round_double(div_bursts_in_flight_map_[div_name].getInUnits("")));
  }
unsigned int IebProfile::getDivisionPercentBW(const string& div_name) const
  {
    divisionCheck(div_name);
    return((unsigned int)round_double(div_percentBW_map_[div_name].getInUnits("")));
  }
string IebProfile::getDivisionAutorad(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_autorad_map_[div_name]);
  }
Uvar IebProfile::getDivisionRip(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_rip_map_[div_name]);
  }
double IebProfile::getDivisionDatarate(const string& div_name) const
  {
    divisionCheck(div_name);
    return(div_rate_map_[div_name].getInUnits(""));
  }
int IebProfile::getDivisionTro(const string& div_name) const
  {
   divisionCheck(div_name);
    return(round_double( div_tro_map_[div_name].getInUnits("")));
  }
