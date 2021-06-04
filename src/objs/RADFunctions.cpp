#include <vector>
#include"RADFunctions.h"
#include "Units.h"
#include "Radiometer.h"
#include "TemplateUtils.h"

using std::cout;
using std::endl;


Uvar  CalibrateTsysandYFactor(const unsigned int& i_beam_rec, 
			      double& Y_factor,
			      const Uvar& avg_time, 
			      const Uvar& Tcoldsky,
			      const vector<Uvar>& beam_rel_time,
			      const vector<Uvar>& beam_ncnt_rl, 
			      const vector<Uvar>& beam_ncnt_radio,
			      const vector<Uvar>& beam_rlotmp, 
			      const vector<bool>& is_beam_coldsky)

  { 
    unsigned int Nbeam_record = beam_rel_time.size();
    if(!is_beam_coldsky[i_beam_rec]) ErrorMessage("Try to calibration Tsys while pointing to target").throwMe();
    vector<Uvar> coldsky_ncnt_rl;
    vector<Uvar> coldsky_ncnt_radio;
    vector<Uvar> coldsky_rlotmp;
    coldsky_ncnt_rl.clear();
    coldsky_ncnt_radio.clear();
    coldsky_rlotmp.clear();
    
    //backward data collection
    for (unsigned int i_cd_b=i_beam_rec;i_cd_b!=0;i_cd_b--)
      {//whithin a minute while pointing cold sky
	if( (beam_rel_time[i_cd_b]<
	     (beam_rel_time[i_beam_rec]-avg_time))
	    || !is_beam_coldsky[i_cd_b]) break;
	coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_b]);
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_b]);
	coldsky_rlotmp.push_back(beam_rlotmp[i_cd_b]);
      }
  
    //forward data collection
    for (unsigned int i_cd_f=i_beam_rec;i_cd_f<Nbeam_record;i_cd_f++)
      {//within a minimute while pointing cold sky
	if( (beam_rel_time[i_cd_f]>
	     (beam_rel_time[i_beam_rec]+avg_time))
	    || !is_beam_coldsky[i_cd_f]) break;
	coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_f]);
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_f]);
	coldsky_rlotmp.push_back(beam_rlotmp[i_cd_f]);
      }
    unsigned int Ncoldsky_record = coldsky_ncnt_rl.size();
    if(Ncoldsky_record!= coldsky_ncnt_radio.size()||
       Ncoldsky_record!=coldsky_rlotmp.size())
      ErrorMessage("container size mismatch: rl radio rlotmp").throwMe();
    Uvar avg_coldsky_ncnt_rl, avg_coldsky_ncnt_radio,avg_coldsky_rlotmp;
    for(unsigned int i_avg=0;i_avg<Ncoldsky_record;i_avg++)
      {
	avg_coldsky_ncnt_rl+=coldsky_ncnt_rl[i_avg];
	avg_coldsky_ncnt_radio+=coldsky_ncnt_radio[i_avg];
	avg_coldsky_rlotmp+=coldsky_rlotmp[i_avg];
      }
    avg_coldsky_ncnt_rl/=double(Ncoldsky_record);
    avg_coldsky_ncnt_radio/=double(Ncoldsky_record);
    avg_coldsky_rlotmp/=double(Ncoldsky_record);
    
    Y_factor=(avg_coldsky_ncnt_radio/avg_coldsky_ncnt_rl).getInUnits("");
    return((Tcoldsky - Y_factor*avg_coldsky_rlotmp)/(Y_factor-1.0));
  }



Uvar  CalibrateTsysUsingConstGain(const unsigned int& i_beam_rec,
				  const Uvar& gain_sys,
				  const Uvar& avg_time, 
				  const Uvar& Tcoldsky,
				  const vector<Uvar>& beam_rel_time,
				  const vector<Uvar>& beam_ncnt_radio,
				  const vector<bool>& is_beam_coldsky) 
  {
    unsigned int Nbeam_record = beam_rel_time.size();
    if(!is_beam_coldsky[i_beam_rec]) ErrorMessage("Try to calibration Tsys while pointing to target").throwMe();
    
    vector<Uvar> coldsky_ncnt_radio;
    coldsky_ncnt_radio.clear();
      
    //backward data collection
    for (unsigned int i_cd_b=i_beam_rec;i_cd_b!=0;i_cd_b--)
      {//whithin a minute while pointing cold sky
	if( (beam_rel_time[i_cd_b]<
	     (beam_rel_time[i_beam_rec]-avg_time))
	    || !is_beam_coldsky[i_cd_b]) break;
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_b]);
      }
    
    //forward data collection
    for (unsigned int i_cd_f=i_beam_rec;i_cd_f<Nbeam_record;i_cd_f++)
      {//within a minimute while pointing cold sky
	if( (beam_rel_time[i_cd_f]>
	     (beam_rel_time[i_beam_rec]+avg_time))
	    || !is_beam_coldsky[i_cd_f]) break;
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_f]);
      }
    
    unsigned int Ncoldsky_record = coldsky_ncnt_radio.size();
    
    Uvar avg_coldsky_ncnt_radio;
    for(unsigned int i_avg=0;i_avg<Ncoldsky_record;i_avg++)
      {
	avg_coldsky_ncnt_radio+=coldsky_ncnt_radio[i_avg];
      }
    avg_coldsky_ncnt_radio/=double(Ncoldsky_record);
    Uvar Tsys = avg_coldsky_ncnt_radio/gain_sys -Tcoldsky;
    return(Tsys );    
  }


Uvar  CalibrateGsysUsingConstTsys(const unsigned int& i_beam_rec,
				  const Uvar& Tsys,
				  const Uvar& avg_time, 
				  const Uvar& Tcoldsky,
				  const vector<Uvar>& beam_rel_time,
				  const vector<Uvar>& beam_ncnt_radio,
				  const vector<bool>& is_beam_coldsky) 
  {
    unsigned int Nbeam_record = beam_rel_time.size();
    if(!is_beam_coldsky[i_beam_rec]) ErrorMessage("Try to calibration system gain  while pointing to target").throwMe();
    
    vector<Uvar> coldsky_ncnt_radio;
    coldsky_ncnt_radio.clear();
       
    //backward data collection
    for (unsigned int i_cd_b=i_beam_rec;i_cd_b!=0;i_cd_b--)
      {//whithin a minute while pointing cold sky
	if( (beam_rel_time[i_cd_b]<
	     (beam_rel_time[i_beam_rec]-avg_time))
	    || !is_beam_coldsky[i_cd_b]) break;
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_b]);
      }
	     
    //forward data collection
    for (unsigned int i_cd_f=i_beam_rec;i_cd_f<Nbeam_record;i_cd_f++)
      {//within a minimute while pointing cold sky
	if( (beam_rel_time[i_cd_f]>
	     (beam_rel_time[i_beam_rec]+avg_time))
	    || !is_beam_coldsky[i_cd_f]) break;
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_f]);
      }
    
    unsigned int Ncoldsky_record = coldsky_ncnt_radio.size();
    
    Uvar avg_coldsky_ncnt_radio;
    for(unsigned int i_avg=0;i_avg<Ncoldsky_record;i_avg++)
      {
	avg_coldsky_ncnt_radio+=coldsky_ncnt_radio[i_avg];
      }
    avg_coldsky_ncnt_radio/=double(Ncoldsky_record);
    
    Uvar Gsys = avg_coldsky_ncnt_radio/(Tsys +Tcoldsky);
    return(Gsys );    
  }


//-----------------------------------
//compute average radiometer count and RL count
//-------------------------------------
void  computeAvg_radio_rl_rlotmp(const unsigned int& i_beam_rec,
			  const Uvar& avg_time, 
			  const vector<Uvar>& beam_rel_time,
			  const vector<Uvar>& beam_ncnt_radio,
			  const vector<Uvar>& beam_ncnt_rl,
			  const vector<Uvar>& beam_rlotmp,
			  const vector<bool>& is_beam_coldsky,
			  Uvar& avg_cnt_radio,
			  Uvar& avg_cnt_rl,
			  Uvar& avg_rlotmp) 
  {
    unsigned int Nbeam_record = beam_rel_time.size();
    if(!is_beam_coldsky[i_beam_rec]) ErrorMessage("Try to calibration system gain  while pointing to target").throwMe();
    vector<Uvar> coldsky_ncnt_rl;
    vector<Uvar> coldsky_ncnt_radio;
    vector<Uvar> coldsky_rlotmp;
   
    coldsky_ncnt_rl.clear();
    coldsky_ncnt_radio.clear();
    coldsky_rlotmp.clear();
    
    
    //backward data collection
    for (unsigned int i_cd_b=i_beam_rec;i_cd_b!=0;i_cd_b--)
      {//whithin a minute while pointing cold sky
	if( (beam_rel_time[i_cd_b]<
	     (beam_rel_time[i_beam_rec]-avg_time))
	    || !is_beam_coldsky[i_cd_b]) break;
	coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_b]);
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_b]);
	coldsky_rlotmp.push_back(beam_rlotmp[i_cd_b]);
      }
  
    //forward data collection
    for (unsigned int i_cd_f=i_beam_rec;i_cd_f<Nbeam_record;i_cd_f++)
      {//within a minimute while pointing cold sky
	if( (beam_rel_time[i_cd_f]>
	     (beam_rel_time[i_beam_rec]+avg_time))
	    || !is_beam_coldsky[i_cd_f]) break;
	coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_f]);
	coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_f]);
	coldsky_rlotmp.push_back(beam_rlotmp[i_cd_f]);
      }
    unsigned int Ncoldsky_record = coldsky_ncnt_rl.size();
    if(Ncoldsky_record!= coldsky_ncnt_radio.size()||
       Ncoldsky_record!=coldsky_rlotmp.size())
      
      ErrorMessage("container size mismatch: rl radio rlotmp wgbtmp").throwMe();

    avg_cnt_radio = Uvar(0,"1/s");
    avg_cnt_rl=Uvar(0,"1/s");
    avg_rlotmp=Uvar(0,"K");
    for(unsigned int i_avg=0;i_avg<Ncoldsky_record;i_avg++)
      {
	avg_cnt_rl+=coldsky_ncnt_rl[i_avg];
	avg_cnt_radio+=coldsky_ncnt_radio[i_avg];
	avg_rlotmp+=coldsky_rlotmp[i_avg];
      }
  avg_cnt_radio/=double(Ncoldsky_record);
  avg_cnt_rl/=double(Ncoldsky_record);
  avg_rlotmp/=double(Ncoldsky_record);
  }

//----------------------------------
//compute system gain by using cold load (RL)
//----------------------------------
void compute_GT_using_RL(const unsigned int& i_beam_rec,
			 const double& frontend_loss_figure ,
			 const Uvar& avg_time,
			 const Uvar& Tcoldsky,
			 const vector<Uvar>& beam_rel_time,
			 const vector<Uvar>& beam_ncnt_radio,
			 const vector<Uvar>& beam_ncnt_rl,
			 const vector<Uvar>& beam_rlotmp,
			 const vector<Uvar>& beam_wgbtmp,
			 const vector<bool>& is_beam_coldsky,
			 Uvar& system_gain,
			 Uvar& receiver_temp){
  
  
  unsigned int Nbeam_record = beam_rel_time.size();
  if(!is_beam_coldsky[i_beam_rec]) ErrorMessage("Try to calibration system gain  while pointing to target").throwMe();
  vector<Uvar> coldsky_ncnt_rl;
  vector<Uvar> coldsky_ncnt_radio;
  vector<Uvar> coldsky_rlotmp;
  vector<Uvar> coldsky_wgbtmp;
  coldsky_ncnt_rl.clear();
  coldsky_ncnt_radio.clear();
  coldsky_rlotmp.clear();
  coldsky_wgbtmp.clear();

  //backward data collection
  for (unsigned int i_cd_b=i_beam_rec;i_cd_b!=0;i_cd_b--)
    {//whithin a minute while pointing cold sky
      if( (beam_rel_time[i_cd_b]<
	   (beam_rel_time[i_beam_rec]-avg_time))
	  || !is_beam_coldsky[i_cd_b]) break;
      coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_b]);
      coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_b]);
      coldsky_rlotmp.push_back(beam_rlotmp[i_cd_b]);
      coldsky_wgbtmp.push_back(beam_wgbtmp[i_cd_b]);
    }
  
  //forward data collection
  for (unsigned int i_cd_f=i_beam_rec;i_cd_f<Nbeam_record;i_cd_f++)
    {//within a minimute while pointing cold sky
      if( (beam_rel_time[i_cd_f]>
	   (beam_rel_time[i_beam_rec]+avg_time))
	  || !is_beam_coldsky[i_cd_f]) break;
      coldsky_ncnt_rl.push_back(beam_ncnt_rl[i_cd_f]);
      coldsky_ncnt_radio.push_back(beam_ncnt_radio[i_cd_f]);
      coldsky_rlotmp.push_back(beam_rlotmp[i_cd_f]);
      coldsky_wgbtmp.push_back(beam_wgbtmp[i_cd_f]);
    }
  unsigned int Ncoldsky_record = coldsky_ncnt_rl.size();
  if(Ncoldsky_record!= coldsky_ncnt_radio.size()||
     Ncoldsky_record!=coldsky_rlotmp.size() ||
     Ncoldsky_record!=coldsky_wgbtmp.size())
    ErrorMessage("container size mismatch: rl radio rlotmp wgbtmp").throwMe();

  Uvar avg_coldsky_ncnt_rl, avg_coldsky_ncnt_radio;
  Uvar avg_coldsky_rlotmp, avg_coldsky_wgbtmp;

  for(unsigned int i_avg=0;i_avg<Ncoldsky_record;i_avg++)
    {
      avg_coldsky_ncnt_rl+=coldsky_ncnt_rl[i_avg];
      avg_coldsky_ncnt_radio+=coldsky_ncnt_radio[i_avg];
      avg_coldsky_rlotmp+=coldsky_rlotmp[i_avg];
      avg_coldsky_wgbtmp+=coldsky_wgbtmp[i_avg];
    }
  avg_coldsky_ncnt_rl/=double(Ncoldsky_record);
  avg_coldsky_ncnt_radio/=double(Ncoldsky_record);
  avg_coldsky_rlotmp/=double(Ncoldsky_record);
  avg_coldsky_wgbtmp/=double(Ncoldsky_record);

  //------------------------
  //rl_counts = system_gain *(rl_temp + receiver_only_temp)
  //sky_counts = system_gain/noise_figure *
  //                ( sky_temp + (noise_figure-1)*waveguid_temp 
  //                    + noise_figure*receiver_only_temp)
  //--------------------------

 
  Uvar Tp =avg_coldsky_wgbtmp;
  
  Radiometer_coldsky rad(avg_coldsky_ncnt_radio, Tcoldsky,
		 avg_coldsky_ncnt_rl, avg_coldsky_rlotmp,
		 Tp,frontend_loss_figure);
  
  rad.computeGainandReceiverTemp(system_gain, receiver_temp);

  /*
  cout<<"avg coldsky count "<< avg_coldsky_ncnt_radio<<endl;
  cout<<"cold sky temp "<< Tcoldsky<<endl;
  cout<<"rl count "<< avg_coldsky_ncnt_rl<<endl;
  cout<<"rl temp "<< avg_coldsky_rlotmp<<endl;
  cout<<"frontend loss "<< frontend_loss_figure<<endl;
  cout<<"gain and receiver temp "<< system_gain<<" "<<receiver_temp<<endl;
  Uvar ratio = (avg_coldsky_rlotmp+receiver_temp);
  ratio /=(receiver_temp+ (1.0-1/frontend_loss_figure)*Tp+Tcoldsky/frontend_loss_figure);
  cout<<"ratio "<< ratio<<endl;
  cout<<"count ratio "<< avg_coldsky_ncnt_rl / avg_coldsky_ncnt_radio<<endl;
  cout<<"system temperature "<< frontend_loss_figure*receiver_temp +(frontend_loss_figure-1.0)*Tp<<endl;
  */
  
}


//-------------------------
//Compute time-averaged RL counts
//---------------------------
void compute_time_averaged_rl(const Uvar& avg_time,
			      const vector<Uvar>& beam_ncnt_rl,
			      const vector<Uvar>& beam_rel_time,
			      vector<Uvar>&  avg_beam_ncnt_rl){
  
  unsigned int Nbeam_record = beam_rel_time.size();
  avg_beam_ncnt_rl.resize(Nbeam_record);
  
  for(unsigned int i=0;i<Nbeam_record;++i){

    vector<Uvar> rl_record;
    //backward data collection
    for (unsigned int i_cd_b=i;i_cd_b!=0;i_cd_b--){
      //whithin a minute while pointing cold sky
      if( beam_rel_time[i_cd_b]<  (beam_rel_time[i]-avg_time)) break;
      rl_record.push_back(beam_ncnt_rl[i_cd_b]);
    }
    //forward search
    //forward data collection
    for (unsigned int i_cd_f=i;i_cd_f<Nbeam_record;i_cd_f++){
      //within a minimute while pointing cold sky
      if( beam_rel_time[i_cd_f]> (beam_rel_time[i]+avg_time)) break;
      rl_record.push_back(beam_ncnt_rl[i_cd_f]);
    }
    
    //take average
    Uvar sum=Uvar(0,"1/s");
    for(unsigned int k=0;k<rl_record.size();++k)
      sum= sum + rl_record[k];

    avg_beam_ncnt_rl[i]=sum/double(rl_record.size());
  }
}
