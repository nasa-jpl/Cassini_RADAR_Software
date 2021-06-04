//-----------------------------------------------------------------------------
// DopplerCentroid.cpp
// This is class will not be used 
//-----------------------------------------------------------------------------

#include <string>

//----------------------------
//Forward declaration
//----------------------------
#include <fstream>
#include <iostream>
#include <iomanip>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "DopplerCentroid.h"
#include "Plot.h"
#include "Constants.h"
#include "TargetGeom.h"
#include "Utils.h"
#include "Ivd.h"
#include "Ckernel.h"
using std::string;
using std::cout;
using std::endl;
using std::setw;

//-------------------------
// DopplerCentroid Methods
//-------------------------
//-------------------------
// Constructors
//-------------------------
DopplerCentroid::DopplerCentroid() 
  :look_direction_set_(false),
   target_radius_set_(false),
   azimuth_steering_angle_set_(false),
   wavelength_set_(false),
   fitting_parameter_set_(false),
   l_chisq_(true),
   i_paramest_("i_paramest"),
   i_usedata_("used data"),
   i_asaa_("i_asaa"),
   yaw_rad_("each beams yaw"),
   pitch_rad_("each beams rad"),
   roll_rad_("roll rad"),
   position_("position"),
   velocity_("velocity"),
   beam_record_indicator_("load beam record"),
   beam_first_range_("beam_first range"),
   beam_delta_range_("beam delta range"),
   azimuth_angles_("azimuth angles"),
   azimuth_bias_angles_("azimuth bias angles"),
   r_vecin_("r_vecin"),
   r_vobs_("r_vobs"),
   r_cov_("r_cov"),
   r_a_("r_a"),
   r_at2_("r_at2"),
   r_u_("r_u"),
   r_v_("r_v"),
   r_w_("r_w"),
   r_chisq_("r_chisq"),
   r_cvm_("r_cmv")
  {
    Nbeam_ = 5;
    i_rd_=1;//only one observation
    i_fp_=18;//18 parameters
    NPP_=2+Nbeam_;//yaw + pitch + 5 beams
    rscale_ = 100000.0;
    MP_ = 10000;//max number of data point
    i_yaw_=0;//first parameter
    i_pitch_=1;//second parameter

    yaw_rad_.resize(Nbeam_);
    pitch_rad_.resize(Nbeam_);
    roll_rad_.resize(Nbeam_);
    position_.resize(Nbeam_,3);
    velocity_.resize(Nbeam_,3);
    beam_record_indicator_.resize(Nbeam_);
    beam_record_indicator_=0;//no beam is loaded
    beam_first_range_.resize(Nbeam_);
    beam_delta_range_.resize(Nbeam_);
    azimuth_angles_.resize(Nbeam_);
    azimuth_bias_angles_.resize(Nbeam_);

    //default
    azimuth_bias_angles_ = Uvar(0,"rad");
    //set azimuth bias angles based on oneway beam gain
    /*
    azimuth_bias_angles_(0)= Uvar(-0.0077*degtorad,"rad");//beam 1
    azimuth_bias_angles_(1)= Uvar( 0.0401*degtorad,"rad");//beam 2
    azimuth_bias_angles_(2)= Uvar(-0.0215*degtorad,"rad");//beam 3
    azimuth_bias_angles_(3)= Uvar( 0.0354*degtorad,"rad");//beam 4
    azimuth_bias_angles_(4)= Uvar(-0.0058*degtorad,"rad");//beam 5
    */
    
    //set azimuth bias angles based on averaged power distribution
    azimuth_bias_angles_(0)= Uvar(-0.014*degtorad,"rad");//beam 1
    azimuth_bias_angles_(1)= Uvar( 0.031*degtorad,"rad");//beam 2
    azimuth_bias_angles_(2)= Uvar(-0.0257*degtorad,"rad");//beam 3
    azimuth_bias_angles_(3)= Uvar( 0.026*degtorad,"rad");//beam 4
    azimuth_bias_angles_(4)= Uvar(-0.013*degtorad,"rad");//beam 5
    

    fsc_=Frame("CASSINI_SC_COORD","Cassini");//cassini frame
    href_ = 0.0;
    x_direction_to_Px_=0;//no direction set
  
    r_vecin_.resize(i_fp_,MP_);
    r_vobs_.resize(i_rd_,MP_);
    r_cov_.resize(i_rd_,i_rd_,MP_);
    r_a_.resize(NPP_);
    r_at2_.resize(NPP_);
    r_u_.resize(NPP_,NPP_);
    r_v_.resize(NPP_,NPP_);
    r_w_.resize(NPP_);
    r_chisq_.resize(i_rd_,MP_);
    r_cvm_.resize(NPP_,NPP_);
    i_usedata_.resize(i_rd_); i_usedata_ = 1;
    i_paramest_.resize(NPP_);
    i_asaa_.resize(NPP_);
    for(unsigned int i=0;i<Nbeam_;++i) i_asaa_(i)= 2 + i;
    r_cov_=1.0;
  }

//look direction
void DopplerCentroid::setLookDirection(const string& look_direction)
 {
   if(look_direction_set_) ErrorMessage("DopplerCentroid.cpp:: look direction has been already set").throwMe();
   if(look_direction=="Right" || look_direction=="right") 
     look_direction_=-1.0;
   else if(look_direction=="Left" || look_direction=="left")
     look_direction_= 1.0;
   else
     ErrorMessage("DopplerCentroid.cpp:: invalid look direction").throwMe();
   look_direction_set_ = true;
 } 

//set target
void DopplerCentroid::setTarget(const string& target)
  {
    if(target_radius_set_) ErrorMessage("DopplerCentroid.cpp:: target radius has been already set").throwMe();
    TargetGeom tg;
    tg.setTarget(target);
    radius_ = tg.radius().getInUnits("m");
    ftarget_=Frame("IAU_"+target,target);
    fj2000_=Frame("J2000",target);
    target_radius_set_ = true;
  } 

//set wavelength
void DopplerCentroid::setWavelength(const Uvar& lambda)
  {
    if(wavelength_set_) ErrorMessage("DopplerCentroid.cpp:: wavelength  has been already set").throwMe();
    lambda_ = lambda.getInUnits("m");
    wavelength_set_ = true;
  }

//compute beam's azimuth steering angle
void DopplerCentroid::computeAzimuthSteeringAnglesAtEpoch(const Time& epoch_time)
  {
    Time t=epoch_time;
    for(unsigned int i_beam=0;i_beam<Nbeam_;++i_beam){
   
      Frame fbeam("CASSINI_RADAR_"+toStr(i_beam+1),"Cassini");
      SpiceDouble et;
      t.getEt(et);
      
      //state vector
      StateVector sc_state;
      //TCN coordinate
      ftarget_.ephemeris(sc_state,"Cassini",t,"NONE");
      DirectionVector sc_z(" ",sc_state.position());
      DirectionVector sc_y(" ",cross(sc_z,sc_state.velocity()));
      DirectionVector sc_x(" ",cross(sc_y,sc_z));
      
      
      //SPACE CRAFT FRAME
      DirectionVector x("",fsc_,t,1,0,0);
      DirectionVector y("",fsc_,t,0,1,0);
      DirectionVector z("",fsc_,t,0,0,1);
      
      //When spacecraft rotated w.r.t. z-axis by 180 degree
      double z_rotation = dot(sc_x, x);
      if(z_rotation < 0)
	{
	  //cout<<"z rotation "<<endl;
	  x = -x;
	  y= -y;
	  x_direction_to_Px_ = -1;
	} 
      else
	{
	  x_direction_to_Px_ = 1;
	  //cout<<"no z rotation "<<endl;
	}
      Frame FSC("new sc frame",sc_state.position(),x,y,z);
      //transformation matrix
      
      //from beam to FSC transformation
      SpiceDouble FSC_to_beam[3][3];	    
      Frame::rotationMatrix(FSC,fbeam,et,FSC_to_beam);//from FSC to beam    
      SpiceDouble* beam_angles;
      beam_angles = new SpiceDouble[3];
      m2eul_c(FSC_to_beam,1,2,3,beam_angles, beam_angles+1,beam_angles+2);
      //cout<<FSC_to_beam[0][0]<< " "<<FSC_to_beam[0][1]<< " "
      //<<FSC_to_beam[0][2]<< " "<<endl;
      //cout<<FSC_to_beam[1][0]<< " "<<FSC_to_beam[1][1]<< " "
      //<<FSC_to_beam[1][2]<< " "<<endl;
      //cout<<FSC_to_beam[2][0]<< " "<<FSC_to_beam[2][1]<< " "
      //<<FSC_to_beam[2][2]<< " "<<endl;
      //cout<<endl;
      //cout<<"FSC to beam "<<  beam_angles[0]*radtodeg<<" "
      //<< beam_angles[1]*radtodeg<<" "
      //  <<beam_angles[2]*radtodeg<<endl;
      azimuth_angles_(i_beam)= beam_angles[1];
     
      //add bias to the angle
      azimuth_angles_(i_beam) += azimuth_bias_angles_(i_beam).getInUnits("rad");

      //no z rotation, switch the sign
      if(z_rotation < 0.0) azimuth_angles_(i_beam) = -azimuth_angles_(i_beam);


     
     
      //double anlge_in_degree = azimuth_angles_(i_beam)*radtodeg;
      //azimuth_angles_(i_beam) = double (round_double(anlge_in_degree*100.0))/100.0*degtorad;
      cout<<"azimuth bias angle of beam  "<<i_beam+1<<" "
	  << azimuth_angles_(i_beam)*radtodeg<<endl;
      delete[] beam_angles;     
    }    
    azimuth_steering_angle_set_ = true;
  }
//set fitting parameter
void DopplerCentroid::setFittingParameters(const bool& fit_yaw,
			     const bool& fit_pitch)
  {
    
    i_paramest_=0;
    if(fit_yaw) i_paramest_(i_yaw_) = 1;
    if(fit_pitch) i_paramest_(i_pitch_)=1;
    if(i_paramest_.sum()==0) 
      ErrorMessage("DopplerCentroid: no fitting paramter is set").throwMe();
		    
    fitting_parameter_set_ = true;
  }


void 
DopplerCentroid
::loadCassini5BeamRecords(const vector<StateVector>& sc_state_record,
		    const vector<unsigned int>& beam_record_indicator,
		    const vector<unsigned int>& range_doppler_record_length,
		    const Umat& range_values,
		    const Umat& doppler_values)

{
  if( (beam_record_indicator.size() != sc_state_record.size()) ||
      (beam_record_indicator.size() != range_doppler_record_length.size()) ||
      (beam_record_indicator.size() != Nbeam_))
    ErrorMessage("DopplerCentroid.cpp: container size mismatch").throwMe();
  if(range_values.rows()!=Nbeam_ || doppler_values.rows()!=Nbeam_) 
    ErrorMessage("DopplerCentroid.cpp:: range and doppler size does not match to number of beams used").throwMe();

 


  for(unsigned int i=0; i<sc_state_record.size();++i){
    unsigned int beam_id= beam_record_indicator[i];
    
    //beam load indicator: 1 loaded 0 not loaded
    beam_record_indicator_(beam_id)++;
    
    //keep beam3 time
    if( (beam_id+1)==3) t_beam3_ = sc_state_record[i].time();
   
    computeTRMatrix(sc_state_record[i], beam_id);
        
    //assume regular interval and increasing order
    Uvar deltar= range_values(beam_id,1)-range_values(beam_id,0);
    if(range_values(beam_id,1) < range_values(beam_id,0)) 
      ErrorMessage("DopplerCentroid:: range is not arranged in increasing order").throwMe();
    
    beam_first_range_(beam_id) = range_values(beam_id,0).getInUnits("m");
    beam_delta_range_(beam_id)= deltar.getInUnits("m");

    for(unsigned int j=0; j<range_doppler_record_length[i]-1;++j){
      if( fabs( range_values(beam_id,j+1)-range_values(beam_id,j)-deltar)>deltar/100.0){
	cout<<j+1 <<"th record is not good "<<endl;
	cout<<range_values(beam_id,j+1)<<" "<<range_values(beam_id,j)<<" "<<deltar<<endl;
	cout<<"diff "<< range_values(beam_id,j+1) - range_values(beam_id,j) <<endl;
	ErrorMessage("DopplerCentroid::irregular range bins ").throwMe();
      }
    }
    
    for(unsigned int j=0; j< range_doppler_record_length[i];++j){   
      record_beam_id_.push_back(beam_id);
      record_time_.push_back(sc_state_record[i].time());
      range_.push_back(range_values(beam_id,j).getInUnits("m"));
      doppler_.push_back(doppler_values(beam_id,j).getInUnits("Hz")); 
    }
  }
}

/*
void 
DopplerCentroid::loadRangeDopplerData(const StateVector& sc_state,
				      const unsigned int& beam_id,
				      const Uvec& range,
				      const Uvec& doppler)
{
  if(range.size()!=doppler.size()) 
    ErrorMessage("range and doppler have different sizes").throwMe();
  if(beam_record_indicator_(beam_id)==1) 
    ErrorMessage("DopplerCentroid:: beam data has been already loaded").throwMe();

  beam_record_indicator_(beam_id)++;//0 means beam not used
  
  if( (beam_id+1)==3) t_beam3_ = sc_state.time();//keep beam3 time

  computeTRMatrix(sc_state, beam_id);

  //assume regular interval and increasing order
  Uvar deltar= range(1)-range(0);
  if(range(1) < range(0)) 
    ErrorMessage("DopplerCentroid:: range is not arranged in increasing order").throwMe();
  

  for(unsigned int i=0; i<range.size()-1;++i)
    if( fabs( range(i+1)-range(i)-deltar)>deltar/100.0)
      {
	cout<<range(i+1)<<" "<<range(i)<<" "<<deltar<<endl;
	cout<<"diff "<< range(i+1) - range(i) <<endl;
	ErrorMessage("DopplerCentroid::irregular range bins ").throwMe();
      }
  

  beam_first_range_(beam_id) = range(0).getInUnits("m");
  beam_delta_range_(beam_id)= deltar.getInUnits("m");

  for(unsigned int i=0;i<range.size();++i){
    record_beam_id_.push_back(beam_id);
    record_time_.push_back(sc_state.time());
    range_.push_back(range(i).getInUnits("m"));
    doppler_.push_back(doppler(i).getInUnits("Hz"));}
  
}

void 
DopplerCentroid::loadRangeDopplerDataInkm_Hz(const StateVector& sc_state,
					     const unsigned int& beam_id,
					     const Dvec& range_in_km,
					     const Dvec& doppler_in_Hz)
{
  if(range_in_km.size()!=doppler_in_Hz.size()) ErrorMessage("range and doppler have different sizes").throwMe(); 
  if(beam_record_indicator_(beam_id)==1) 
    ErrorMessage("DopplerCentroid:: beam data has been already loaded").throwMe();
  
  beam_record_indicator_(beam_id)++;//0 means beam not used
  
  if( (beam_id+1)==3) t_beam3_ = sc_state.time();//keep beam3 time
  
  //compute Transformation Matrix
  computeTRMatrix(sc_state, beam_id);
  
  double deltar= range_in_km(1)-range_in_km(0);
  if(range_in_km(1) < range_in_km(0))
    ErrorMessage("DopplerCentroid:: range is not arranged in increasing order").throwMe();
  

  for(unsigned int i=0; i<range_in_km.size()-1;++i)
    if( (range_in_km(i+1)-range_in_km(i))!=deltar) ErrorMessage("DopplerCentroid::irregular range bins ").throwMe();

 
  beam_first_range_(beam_id) = range_in_km(0)*1000.0;
  beam_delta_range_(beam_id)= deltar*1000;

  for(unsigned int i=0;i<range_in_km.size();++i){
    record_beam_id_.push_back(beam_id);
    record_time_.push_back(sc_state.time());
    range_.push_back(range_in_km(i)*1000.0);
    doppler_.push_back(doppler_in_Hz(i));}
    
}

*/


void DopplerCentroid::computeYawPitch(Uvar& yaw,
				      Uvar& delta_yaw,
				      Uvar& pitch,
				      Uvar& delta_pitch,
				      const bool& show_plot)
{
  //see whether fitting parameters are set
  if(!fitting_parameter_set_) 
    ErrorMessage("DopplerCentroid::No fitting paramter is set").throwMe();
  //check container size
  if( (record_beam_id_.size() != record_time_.size()) &&
      (record_time_.size() != range_.size()) &&
      (range_.size() != doppler_.size())){
    cout<<"container size beam_id, time, range, doppler "<<endl;
    cout<<record_beam_id_.size()<<" "<< record_time_.size()<<" "<< range_.size()
	<<" "<<doppler_.size()<<endl;
    ErrorMessage("DopplerCentroid.cpp:: container size mismatch").throwMe();
  }
  //check transformation matrix is set
  if(!azimuth_steering_angle_set_) ErrorMessage("DopplerCentroid::beam's azimuth and steering angle is not set").throwMe();
  
  //check 0 data record
  if(record_beam_id_.size()==0) 
    ErrorMessage("No data has been loaded").throwMe();
  
  //check how many beams are being used
  for(unsigned int i=0;i<Nbeam_;++i){
    if(beam_record_indicator_(i)==0) 
      ErrorMessage("beam "+toStr(i+1)+" is not loaded").throwMe();
    if(beam_record_indicator_(i) > 1)
      ErrorMessage("beam "+toStr(i+1) +" is loaded more than once").throwMe();
  }
  
  

  //generate range and doppler data based on yaw, pitch and range
 
  for(unsigned int i=0;i<range_.size();++i)
    {
      r_sch_[0]= position_(record_beam_id_[i],0);
      r_sch_[1]= position_(record_beam_id_[i],1);
      r_sch_[2]= position_(record_beam_id_[i],2);

      r_schvel_[0]=velocity_(record_beam_id_[i],0);
      r_schvel_[1]=velocity_(record_beam_id_[i],1);
      r_schvel_[2]=velocity_(record_beam_id_[i],2);


      dop_derivatives(radius_,
			 r_sch_,
			 r_schvel_,
			 yaw_rad_(record_beam_id_[i]),
			 pitch_rad_(record_beam_id_[i]),
			 azimuth_angles_(record_beam_id_[i]),
			 href_,
			 range_[i],
			 lambda_,
			 look_direction_,
			 rscale_,
			  r_f_,
			 r_dfdy_,
			 r_dfdp_,
			 r_dfdvs_,
			 r_dfdvc_,
			 r_dfdvh_,
			 r_dfdhref_,
			 r_dfdwvl_,
			 r_dfdasa_);
     
      doppler_ck_.push_back(r_f_);
      r_vobs_(0,i)= doppler_[i]- r_f_;//store the difference
      r_cov_(0,0,i)= 1.0;//variance
      r_vecin_(0,i) = radius_;
      r_vecin_(1,i) = r_sch_[0];
      r_vecin_(2,i) = r_sch_[1];
      r_vecin_(3,i) = r_sch_[2];//altitude 
      r_vecin_(4,i) = r_schvel_[0];
      r_vecin_(5,i) = r_schvel_[1];
      r_vecin_(6,i) = r_schvel_[2];
      r_vecin_(7,i) = yaw_rad_(record_beam_id_[i]);
      r_vecin_(8,i) = pitch_rad_(record_beam_id_[i]);
      r_vecin_(9,i) = href_;//set to 0
      r_vecin_(10,i) = beam_first_range_(record_beam_id_[i]) ;//first range bin
      r_vecin_(11,i) = beam_delta_range_(record_beam_id_[i]);//delta range
      r_vecin_(12,i) = lambda_;
      r_vecin_(13,i) = rscale_;
      r_vecin_(14,i) = record_beam_id_[i];//beam number
      r_vecin_(15,i) = look_direction_;//look direction(-1 for right, 1 for left) 
      r_vecin_(16,i) =azimuth_angles_(record_beam_id_[i]);
      r_vecin_(17,i) = (unsigned int) 
	round_double( (range_[i]- beam_first_range_(record_beam_id_[i]))
		      /beam_delta_range_(record_beam_id_[i]));    
       
    }
  
  //debug
  //cout<<"each beam's start range "<< beam_first_range_<<endl;
  //cout<<"each beams's delta range "<< beam_delta_range_<<endl;
  //construct covariant matrix


  //svdvecfit
 svdvecfit(range_.size(),i_rd_,i_fp_,r_vecin_,r_vobs_,
	   r_cov_,NPP_,r_a_,r_at2_,r_u_,r_v_,r_w_,
	   r_chisq_,l_chisq_,
	   i_paramest_,i_usedata_); 
 yaw_bias_rad_ =  r_at2_(i_yaw_)/rscale_;
 pitch_bias_rad_ = r_at2_(i_pitch_)/rscale_;
 //svdvar
 svdvar(r_v_, r_w_, r_cvm_);

 //debug
 //cout<<" Initial yaw"<<yaw_rad_*radtodeg;
 //cout<<" Yaw bias "<<yaw_bias_rad_*radtodeg;
 //cout<<" Initial pitch "<< pitch_rad_*radtodeg;
 //cout<<" Pitch bias "<<pitch_bias_rad_*radtodeg;
 //cout<<" intial roll angle "<< roll_rad_*radtodeg;
 //cout<<endl;
  
 //return variables
 yaw= Uvar(yaw_rad_(2),"rad");//beam3's yaw
 delta_yaw= Uvar(yaw_bias_rad_,"rad");
  
 pitch=Uvar(pitch_rad_(2),"rad");//beam3's pitch
 delta_pitch=Uvar(pitch_bias_rad_,"rad");

 //---------------------------------------------
 //store result for ivd file generation
 //need to add - sign because of direction of transformation from FTCN to FSC 
 //----------------------------------------------------
 
 //---------------------------
 //Time
 //---------------------------
 if(beam_record_indicator_(2)==0)
   {//beam3 is not loaded, error at this time
     ErrorMessage("DopplerCentroid:: beam3 record is not avail ").throwMe();
   }
 beam3_time_record_.push_back(t_beam3_);


 //-------------------------
 //Yaw Pitch Roll: + transformation from FSC to FTCN
 // To make a proper transitin from FTCN to FSC, reverse sign
 //-------------------------

 //---------------------------------
 //yaw::beam3 ckernel + range_doppler fit
 //----------------------------------
 yaw_record_ftcn_to_fsc_.push_back( -Uvar(yaw_rad_(2)+yaw_bias_rad_,"rad"));
 //----------------------------------
 //pitch::beam 3 ckernel + range_doppler fit
 //----------------------------------
 pitch_record_ftcn_to_fsc_.push_back( -Uvar(pitch_rad_(2)+pitch_bias_rad_,"rad"));
 //---------------------------
 //Roll
 //--------------------------
 roll_record_ftcn_to_fsc_.push_back( -Uvar(roll_rad_(2),"rad"));//beam 3 roll from ckernel
  
 //-------------------------------------------
 //if the user want to see the fitting result
 //-----------------------------------------
 if(show_plot)
   {
     for(unsigned int i_beam =0;i_beam<Nbeam_;i_beam++){
       cout<<"beam number "<< i_beam+1<<endl;
       cout<<"Position "<<  position_.getRow(i_beam)<<endl;
       cout<<"velocity "<<velocity_.getRow(i_beam)<<endl;
       cout<<"yaw "<< yaw_rad_(i_beam)*radtodeg<<endl;
       cout<<"pitch "<< pitch_rad_(i_beam)*radtodeg<<endl;
       cout<<"azimuth "<<azimuth_angles_(i_beam)*radtodeg<<endl;
       cout<<"-------------------------------------------"<<endl;
     }
 
     cout<<"----------------------------------------------------- "<<endl;
     cout<<"fitting results "<<endl;
     cout<<"yaw bias in deg "<<yaw_bias_rad_*radtodeg<<endl;
     cout<<"pitch bias in deg "<<pitch_bias_rad_*radtodeg<<endl;
     cout<<"----------------------------------------------------- "<<endl;
     for(unsigned int i=0;i<range_.size();++i){
       r_sch_[0]= position_(record_beam_id_[i],0);
       r_sch_[1]= position_(record_beam_id_[i],1);
       r_sch_[2]= position_(record_beam_id_[i],2);
       
       r_schvel_[0]=velocity_(record_beam_id_[i],0);
       r_schvel_[1]=velocity_(record_beam_id_[i],1);
       r_schvel_[2]=velocity_(record_beam_id_[i],2);
       dop_derivatives(radius_,
		       r_sch_,
		       r_schvel_,
		       yaw_rad_(record_beam_id_[i])+yaw_bias_rad_,
		       pitch_rad_(record_beam_id_[i])+pitch_bias_rad_,
		       azimuth_angles_(record_beam_id_[i]),
		       href_,
		       range_[i],
		       lambda_,
		       look_direction_,
		       rscale_,
		       r_f_,
		       r_dfdy_,
		       r_dfdp_,
		       r_dfdvs_,
		       r_dfdvc_,
		       r_dfdvh_,
		       r_dfdhref_,
		       r_dfdwvl_,
		       r_dfdasa_);
       if(fabs(r_f_)<10e6){
	 new_doppler_.push_back(r_f_);//store expected doppler
       }
       else{
	 cout<<"unusual fitting "<< r_f_<<endl;
	 new_doppler_.push_back(0);
       }
       
       diff_doppler_.push_back(r_f_ - doppler_[i]);
       //cout<<"fitting "<< range_[i]<<" "<< new_doppler_.back()<<endl;
     }   
     displayFitting();
   }
   
 //reset after computation
 reset();
 
 
}



void 
DopplerCentroid::generateCkernelWithBias(const Time& start_time,
					 const Time& end_time,
					 const Uvar& time_step,
					 const Uvar& yaw_bias_from_fsc_to_ftcn,
				       const Uvar& pitch_bias_from_fsc_to_ftcn,
					 string& ck_filename)
{

    //construct spacecraft attitude using
    //beam3_time_record_
    //yaw_record_ftcn_to_fsc_
    //pitch_record_ftcn_to_fsc_
    //roll_record_ftcn_to_fsc_
  /*
    if(beam3_time_record_.size()==0) 
      ErrorMessage("DopplerCentroid::no record to generate ivd file").throwMe();
    if( (beam3_time_record_.size() != yaw_record_ftcn_to_fsc_.size()) &&
	(beam3_time_record_.size() != pitch_record_ftcn_to_fsc_.size()) &&
	(beam3_time_record_.size() != roll_record_ftcn_to_fsc_.size()))
      ErrorMessage("DopplerCentroid:: time, yaw, pitch, roll container size mismatch").throwMe();
  */

  if(end_time < start_time) ErrorMessage("end time is earlier than start time ").throwMe();
  if(time_step<= Uvar(0,"s")) ErrorMessage("time step less than 0").throwMe();
  string ivd_filename="doppler_centroid.ivd";
  ofstream  ivd(ivd_filename.c_str());
  if (!ivd) throw ErrorMessage("Can't create file" + ivd_filename);
  
  
  unsigned int Ntime = (unsigned int) round_double( ( (end_time - start_time)/time_step).getInUnits(""));
  Dmat x_direction("x_position"),y_direction("y_direction"),z_direction("z_direction");
  x_direction.resize(Ntime,3);
  y_direction.resize(Ntime,3);
  z_direction.resize(Ntime,3);
  DirectionVector x,y,z;
  Uvar zero_deg = Uvar(0,"rad");
  
  Uvar yaw0, pitch0, roll0;

  Time t;
  for(unsigned int i=0; i<Ntime;++i){
    
    t = start_time;
    t += time_step * i;

    StateVector sc_state;
    ftarget_.ephemeris(sc_state,"Cassini",t,"NONE");
    computeYaw_Pitch_Roll(sc_state,ftarget_,fsc_,yaw0,pitch0,roll0);//rotation from FSC to FTCN

   
    //compute yaw and pitch angles
    DirectionVector sc_z(" ",sc_state.position());
    DirectionVector sc_y(" ",cross(sc_z, sc_state.velocity()));
    DirectionVector sc_x(" ",cross(sc_y,sc_z));
    Frame ftcn("tnc",sc_state.position(),sc_x,sc_y,sc_z);//ftcn coordinate	
   
    //from FTCN to FSC: 1-2-3 rotation: YAW first
    Rotation azimuth_rot(zero_deg,zero_deg,-(yaw0+yaw_bias_from_fsc_to_ftcn),
			 0,1,2);
    Frame fazimuth("azimuth_rot",ftcn,azimuth_rot);

    //Rotation around y_axis with pitch bias:: PITCH second
    Rotation pitch_rotation(zero_deg,-(pitch0+pitch_bias_from_fsc_to_ftcn),zero_deg,0,1,2);
    Frame fpitch("pitch rotation ", fazimuth, pitch_rotation);
    
    //Rotation around x_axis:: ROLL
    Rotation look_rot(-roll0,zero_deg,zero_deg,0,1,2);
    Frame fsc("space_craft",fpitch,look_rot);
    
    x = DirectionVector("x",fsc,t,1,0,0);
    y = DirectionVector("y",fsc,t,0,1,0);
    z = DirectionVector("z",fsc,t,0,0,1);
    if(x_direction_to_Px_==-1)
      {//spacecraft x direction to -vx
	x = -x;
	y=-y;
      }
    x.representIn(fj2000_);
    y.representIn(fj2000_);
    z.representIn(fj2000_);
    x_direction(i,0) = x[DirectionVector::X];
    x_direction(i,1) = x[DirectionVector::Y];
    x_direction(i,2) = x[DirectionVector::Z];
    
    z_direction(i,0) = z[DirectionVector::X];
    z_direction(i,1) = z[DirectionVector::Y];
    z_direction(i,2) = z[DirectionVector::Z];
    //cout<<"output using setAzimuthIncidence"<<endl;
    //cout<<"X-direction "<<x<<endl;
    //cout<<"Y-direction "<<y<<endl;
    //cout<<"Z-direction "<<z<<endl;
    
  }//loop over time
  


    //------------------------------
    //Current product date
    //-----------------------------
    string product_create_date=Time::getCurrentLocalTime();
   

    //-----------------------
    //According to the IVD SIS
    // only two specifications in the header are important
    //1.  SUBTYPE = IVD;
    //2.  PARAMETER_API_VERSION_ID = 1.0
    // we need to know more about what the version id means
    //-----------------------
    ivd<<"CCSD3ZF0000100000001NJPL3KS0L015MSASPARM"<<endl;
    ivd<<"MISSION_NAME = CASSINI"<<endl;
    ivd<<"MISSION_ID = *"<<endl;
    ivd<<"SPACECRAFT_NAME = *"<<endl;
    ivd<<"SPACECRAFT_ID = *"<<endl;
    ivd<<"SPACECRAFT_SUBSYSTEM = CASSINI"<<endl;
    ivd<<"DATA_SET_ID = ivd "<<endl;
    ivd<<"FILE_NAME = *"<<endl;
    ivd<<"OBJECT_NAME = *"<<endl;
    ivd<<"PRODUCER_ID = *"<<endl;
    ivd<<"PRODUCER_FULL_NAME = *"<<endl;
    ivd<<"PRODUCER_PROGRAM = *"<<endl;
    ivd<<"PRODUCT_CREATION_TIME = "<<product_create_date<<endl;
    ivd<<"FILE_RECORDS = *"<<endl;
    ivd<<"SPACECRAFT_CLOCK_START_COUNT = *"<<endl;
    ivd<<"SPACECRAFT_CLOCK_STOP_COUNT = *"<<endl;
    ivd<<"EARTH_RECEIVE_START_TIME = *"<<endl;
    ivd<<"EARTH_RECEIVE_STOP_TIME = *"<<endl;
    ivd<<"SPACECRAFT_EVENT_START_TIME = *"<<endl;
    ivd<<"SPACECRAFT_EVENT_STOP_TIME = *"<<endl;
    ivd<<"RECORD_CREATION_START_TIME = *"<<endl;
    ivd<<"RECORD_CREATION_STOP_TIME = *"<<endl;
    ivd<<"DESCRIPTION = *"<<endl;
    ivd<<"INPUT_FILES = *"<<endl;
    ivd<<"HISTORY = *"<<endl;
    ivd<<"FILE_STATUS = *"<<endl;
    ivd<<"SEQUENCE_ID = *"<<endl;
    ivd<<"SUBTYPE = IVD"<<endl;
    ivd<<"PARAMETER_API_VERSION_ID = 1.0"<<endl;
    ivd<<"APPLICABLE_START_TIME = "<<start_time.utc("ISOD")<<endl;
    ivd<<"APPLICABLE_STOP_TIME  = "<<end_time.utc("ISOD")<<endl;
    ivd<<"CCSD$$MARKERMSASPARMNJPL3IF0037600000001"<<endl;
    
    //----------------------
    //writing plus x position
    //----------------------
    
    
    ivd<<"START:      "<<"'"<< start_time.utc("ISOD")<<"'"<<endl;
    ivd<<"END:        "<<"'"<< end_time.utc("ISOD")<<"'"<<endl;
    ivd<<"HEAD:       "<<"'TITAN_IVD_RADAR_X'"<<endl;
    ivd<<"BASE:       "<<"'CASSINI'"<<endl;
    ivd<<" "<<endl;
    for (unsigned int i = 0 ; i < Ntime; i++)
      {
	t = start_time;
	t += time_step * i;
       
	ivd<<"Time:    "<<"'"<<t.utc("ISOD")<<"'"<<endl;
	ivd.precision(10);
	ivd<<"Position:  "<<"  "<<x_direction(i,0)<<setw(5)<<"  "<<x_direction(i,1)<<setw(5)<<"  "<<x_direction(i,2)<<endl;
	ivd<<" "<<endl;     
      }
    
    //--------------------
    // writing -z position
    //--------------------
    ivd<<"START:      "<<"'"<< start_time.utc("ISOD")<<"'"<<endl;
    ivd<<"END:        "<<"'"<< end_time.utc("ISOD")<<"'"<<endl;
    ivd<<"HEAD:       "<<"'TITAN_IVD_RADAR_MZ'"<<endl;
    ivd<<"BASE:       "<<"'CASSINI'"<<endl;
    ivd<<" "<<endl;
    for (unsigned int i = 0 ; i < Ntime; i++)
      {
	t = start_time;
	t += time_step * i;
	ivd<<"Time:    "<<"'"<<t.utc("ISOD")<<"'"<<endl;
	ivd.precision(10);
	ivd<<"Position:  "<<"  "<<-z_direction(i,0)<<setw(5)<<"  "<<-z_direction(i,1)<<setw(5)<<"  "<<-z_direction(i,2)<<endl;
	ivd<<" "<<endl; 
      }
    ivd.close();


   

    ifstream ivd_file(ivd_filename.c_str());
    if(!ivd_file) throw ErrorMessage("Can't opne file"+ivd_filename);
    Ivd ivd_dop;
    Ckernel ckernel;
 
    ivd_dop.ReadDataFile(ivd_file);
    ckernel.LoadIvd(ivd_dop);
    ckernel.GenerateQuats();
    ckernel.WriteData(ck_filename);
    ivd_file.close();
  }

//------------------------------
//generate 
//---------------------------
void 
DopplerCentroid::generateCkernelTimeVaryingBias(const Time& start_time,
						const Time& end_time,
						const Time& epoch_time,
						const Uvar& time_step,
						const Dvec& yaw_coeff,
						const Dvec& pitch_coeff,
						const vector<Uvar>& ca_time,
						const Dvec& ca_diff_coeff,
						string& ck_filename){
  if(end_time < start_time) ErrorMessage("end time is earlier than start time ").throwMe();
  if(time_step<= Uvar(0,"s")) ErrorMessage("time step less than 0").throwMe();
  string ivd_filename="doppler_centroid.ivd";
  ofstream  ivd(ivd_filename.c_str());
  if (!ivd) throw ErrorMessage("Can't create file" + ivd_filename);
  
  
  unsigned int Ntime = (unsigned int) round_double( ( (end_time - start_time)/time_step).getInUnits(""));
  Dmat x_direction("x_position"),y_direction("y_direction"),z_direction("z_direction");
  x_direction.resize(Ntime,3);
  y_direction.resize(Ntime,3);
  z_direction.resize(Ntime,3);
  DirectionVector x,y,z;
  Uvar zero_deg = Uvar(0,"rad");
  
  Uvar yaw0, pitch0, roll0;

  Time t;
 
  double epoch_relative_time;
  double yp, yy;
  Uvar pitch_bias_from_fsc_to_ftcn;
  Uvar yaw_bias_from_fsc_to_ftcn;

  for(unsigned int i=0; i<Ntime;++i){
    
    t = start_time;
    t += time_step * i;

    StateVector sc_state;
    ftarget_.ephemeris(sc_state,"Cassini",t,"NONE");
    computeYaw_Pitch_Roll(sc_state,ftarget_,fsc_,yaw0,pitch0,roll0);//rotation from FSC to FTCN

   
    //compute yaw and pitch angles
    DirectionVector sc_z(" ",sc_state.position());
    DirectionVector sc_y(" ",cross(sc_z, sc_state.velocity()));
    DirectionVector sc_x(" ",cross(sc_y,sc_z));
    Frame ftcn("tnc",sc_state.position(),sc_x,sc_y,sc_z);//ftcn coordinate	
 
    epoch_relative_time =  (t - epoch_time).getInUnits("min");

    //---------------
    //fit yaw and pitch angles: global fitting
    //---------------
    yy = yaw_coeff(0);
    yp= pitch_coeff(0);
    for(unsigned int k=1;k<yaw_coeff.size();++k){
      yy += pow(epoch_relative_time, k)*yaw_coeff(k);
    }

    for(unsigned int k=1;k<pitch_coeff.size();++k){
      yp += pow(epoch_relative_time, k)*pitch_coeff(k);
    }
    pitch_bias_from_fsc_to_ftcn=Uvar(yp*degtorad,"rad");
    yaw_bias_from_fsc_to_ftcn=Uvar(yy*degtorad,"rad");


    //-----------------
    //special fitting around CA
    //------------------
    if(ca_time.size()!=0){
      if( (t-epoch_time) >=  ca_time[0] &&
	  (t-epoch_time) <= ca_time[ca_time.size()-1]) {
	double yp=ca_diff_coeff(0);
	for(unsigned int k=1;k<ca_diff_coeff.size();++k)
	  yp += pow(epoch_relative_time,k)*ca_diff_coeff(k);
	pitch_bias_from_fsc_to_ftcn += Uvar(yp*degtorad,"rad");
      }
    }
     
    //---------------------
    //Construct transformation
    //-----------------------


    //from FTCN to FSC: 1-2-3 rotation: YAW first
    Rotation azimuth_rot(zero_deg,zero_deg,-(yaw0+yaw_bias_from_fsc_to_ftcn),
			 0,1,2);
    Frame fazimuth("azimuth_rot",ftcn,azimuth_rot);

    //Rotation around y_axis with pitch bias:: PITCH second
    Rotation pitch_rotation(zero_deg,-(pitch0+pitch_bias_from_fsc_to_ftcn),zero_deg,0,1,2);
    Frame fpitch("pitch rotation ", fazimuth, pitch_rotation);
    
    //Rotation around x_axis:: ROLL
    Rotation look_rot(-roll0,zero_deg,zero_deg,0,1,2);
    Frame fsc("space_craft",fpitch,look_rot);
    
    x = DirectionVector("x",fsc,t,1,0,0);
    y = DirectionVector("y",fsc,t,0,1,0);
    z = DirectionVector("z",fsc,t,0,0,1);
    if(x_direction_to_Px_==-1)
      {//spacecraft x direction to -vx
	x = -x;
	y=-y;
      }
    x.representIn(fj2000_);
    y.representIn(fj2000_);
    z.representIn(fj2000_);
    x_direction(i,0) = x[DirectionVector::X];
    x_direction(i,1) = x[DirectionVector::Y];
    x_direction(i,2) = x[DirectionVector::Z];
    
    z_direction(i,0) = z[DirectionVector::X];
    z_direction(i,1) = z[DirectionVector::Y];
    z_direction(i,2) = z[DirectionVector::Z];
    //cout<<"output using setAzimuthIncidence"<<endl;
    //cout<<"X-direction "<<x<<endl;
    //cout<<"Y-direction "<<y<<endl;
    //cout<<"Z-direction "<<z<<endl;
    
  }//loop over time
  


    //------------------------------
    //Current product date
    //-----------------------------
    string product_create_date=Time::getCurrentLocalTime();
   

    //-----------------------
    //According to the IVD SIS
    // only two specifications in the header are important
    //1.  SUBTYPE = IVD;
    //2.  PARAMETER_API_VERSION_ID = 1.0
    // we need to know more about what the version id means
    //-----------------------
    ivd<<"CCSD3ZF0000100000001NJPL3KS0L015MSASPARM"<<endl;
    ivd<<"MISSION_NAME = CASSINI"<<endl;
    ivd<<"MISSION_ID = *"<<endl;
    ivd<<"SPACECRAFT_NAME = *"<<endl;
    ivd<<"SPACECRAFT_ID = *"<<endl;
    ivd<<"SPACECRAFT_SUBSYSTEM = CASSINI"<<endl;
    ivd<<"DATA_SET_ID = ivd "<<endl;
    ivd<<"FILE_NAME = *"<<endl;
    ivd<<"OBJECT_NAME = *"<<endl;
    ivd<<"PRODUCER_ID = *"<<endl;
    ivd<<"PRODUCER_FULL_NAME = *"<<endl;
    ivd<<"PRODUCER_PROGRAM = *"<<endl;
    ivd<<"PRODUCT_CREATION_TIME = "<<product_create_date<<endl;
    ivd<<"FILE_RECORDS = *"<<endl;
    ivd<<"SPACECRAFT_CLOCK_START_COUNT = *"<<endl;
    ivd<<"SPACECRAFT_CLOCK_STOP_COUNT = *"<<endl;
    ivd<<"EARTH_RECEIVE_START_TIME = *"<<endl;
    ivd<<"EARTH_RECEIVE_STOP_TIME = *"<<endl;
    ivd<<"SPACECRAFT_EVENT_START_TIME = *"<<endl;
    ivd<<"SPACECRAFT_EVENT_STOP_TIME = *"<<endl;
    ivd<<"RECORD_CREATION_START_TIME = *"<<endl;
    ivd<<"RECORD_CREATION_STOP_TIME = *"<<endl;
    ivd<<"DESCRIPTION = *"<<endl;
    ivd<<"INPUT_FILES = *"<<endl;
    ivd<<"HISTORY = *"<<endl;
    ivd<<"FILE_STATUS = *"<<endl;
    ivd<<"SEQUENCE_ID = *"<<endl;
    ivd<<"SUBTYPE = IVD"<<endl;
    ivd<<"PARAMETER_API_VERSION_ID = 1.0"<<endl;
    ivd<<"APPLICABLE_START_TIME = "<<start_time.utc("ISOD")<<endl;
    ivd<<"APPLICABLE_STOP_TIME  = "<<end_time.utc("ISOD")<<endl;
    ivd<<"CCSD$$MARKERMSASPARMNJPL3IF0037600000001"<<endl;
    
    //----------------------
    //writing plus x position
    //----------------------
    
    
    ivd<<"START:      "<<"'"<< start_time.utc("ISOD")<<"'"<<endl;
    ivd<<"END:        "<<"'"<< end_time.utc("ISOD")<<"'"<<endl;
    ivd<<"HEAD:       "<<"'TITAN_IVD_RADAR_X'"<<endl;
    ivd<<"BASE:       "<<"'CASSINI'"<<endl;
    ivd<<" "<<endl;
    for (unsigned int i = 0 ; i < Ntime; i++)
      {
	t = start_time;
	t += time_step * i;
       
	ivd<<"Time:    "<<"'"<<t.utc("ISOD")<<"'"<<endl;
	ivd.precision(10);
	ivd<<"Position:  "<<"  "<<x_direction(i,0)<<setw(5)<<"  "<<x_direction(i,1)<<setw(5)<<"  "<<x_direction(i,2)<<endl;
	ivd<<" "<<endl;     
      }
    
    //--------------------
    // writing -z position
    //--------------------
    ivd<<"START:      "<<"'"<< start_time.utc("ISOD")<<"'"<<endl;
    ivd<<"END:        "<<"'"<< end_time.utc("ISOD")<<"'"<<endl;
    ivd<<"HEAD:       "<<"'TITAN_IVD_RADAR_MZ'"<<endl;
    ivd<<"BASE:       "<<"'CASSINI'"<<endl;
    ivd<<" "<<endl;
    for (unsigned int i = 0 ; i < Ntime; i++)
      {
	t = start_time;
	t += time_step * i;
	ivd<<"Time:    "<<"'"<<t.utc("ISOD")<<"'"<<endl;
	ivd.precision(10);
	ivd<<"Position:  "<<"  "<<-z_direction(i,0)<<setw(5)<<"  "<<-z_direction(i,1)<<setw(5)<<"  "<<-z_direction(i,2)<<endl;
	ivd<<" "<<endl; 
      }
    ivd.close();


   

    ifstream ivd_file(ivd_filename.c_str());
    if(!ivd_file) throw ErrorMessage("Can't opne file"+ivd_filename);
    Ivd ivd_dop;
    Ckernel ckernel;
 
    ivd_dop.ReadDataFile(ivd_file);
    ckernel.LoadIvd(ivd_dop);
    ckernel.GenerateQuats();
    ckernel.WriteData(ck_filename);
    ivd_file.close();
 
}

Uvar 
DopplerCentroid::computeDopplerAtZeroAzimuth(
			  const StateVector& sc_state_target_frame,
			  const unsigned int& beam_id,
			  const Uvar& range)
{
  if(!look_direction_set_ ||
     !target_radius_set_ ||
     !azimuth_steering_angle_set_||
     !wavelength_set_) ErrorMessage("can not compute").throwMe();

  computeTRMatrix(sc_state_target_frame,beam_id);
  //get postion
  r_sch_[0]= position_(beam_id,0);
  r_sch_[1]= position_(beam_id,1);
  r_sch_[2]= position_(beam_id,2);
       
  //get velocity in sch
  r_schvel_[0]=velocity_(beam_id,0);
  r_schvel_[1]=velocity_(beam_id,1);
  r_schvel_[2]=velocity_(beam_id,2);
  //range in meter
  double r_range = range.getInUnits("km")*1000.0;
 
  dop_derivatives(radius_,
		  r_sch_,
		  r_schvel_,
		  yaw_rad_(beam_id),
		  pitch_rad_(beam_id),
		  azimuth_angles_(beam_id),
		  href_,
		  r_range,
		  lambda_,
		  look_direction_,
		  rscale_,
		  r_f_,
		  r_dfdy_,
		  r_dfdp_,
		  r_dfdvs_,
		  r_dfdvc_,
		  r_dfdvh_,
		  r_dfdhref_,
		  r_dfdwvl_,
		  r_dfdasa_);


  //if(beam_id==1){
  //cout<<"dop class: radius sch schvel yaw pitch azimuht "<< radius_ <<endl;
  //for(unsigned int i=0;i<3;++i) cout<<r_sch_[i]<<" ";
  //cout<<endl;
  //for(unsigned int i=0;i<3;++i) cout<<r_schvel_[i]<<" ";
  //cout<<endl;
  //cout<<yaw_rad_(beam_id)*radtodeg<<" "<<pitch_rad_(beam_id)*radtodeg<<" "<<azimuth_angles_(beam_id)*radtodeg<<endl;
  //cout<<"doppler "<< r_f_<<endl;
  //}
  
  return(Uvar(r_f_,"Hz"));


}

void
DopplerCentroid::computeDopplerAtZeroAzimuth(
			  const StateVector& sc_state_target_frame,
			  const unsigned int& beam_id,
			  const Uvec& ranges,
			  Uvec& dopplers)
{
  if(!look_direction_set_ ||
     !target_radius_set_ ||
     !azimuth_steering_angle_set_||
     !wavelength_set_) ErrorMessage("can not compute").throwMe();
  if(dopplers.size() != ranges.size()) dopplers.resize(ranges.size());

  computeTRMatrix(sc_state_target_frame,beam_id);
  //get postion
  r_sch_[0]= position_(beam_id,0);
  r_sch_[1]= position_(beam_id,1);
  r_sch_[2]= position_(beam_id,2);
       
  //get velocity in sch
  r_schvel_[0]=velocity_(beam_id,0);
  r_schvel_[1]=velocity_(beam_id,1);
  r_schvel_[2]=velocity_(beam_id,2);
  
  //range in meter
  for(unsigned int i=0;i<ranges.size();++i){
    double r_range = ranges(i).getInUnits("km")*1000.0;
    dop_derivatives(radius_,
		    r_sch_,
		    r_schvel_,
		    yaw_rad_(beam_id),
		    pitch_rad_(beam_id),
		    azimuth_angles_(beam_id),
		    href_,
		    r_range,
		    lambda_,
		    look_direction_,
		    rscale_,
		    r_f_,
		    r_dfdy_,
		    r_dfdp_,
		    r_dfdvs_,
		    r_dfdvc_,
		    r_dfdvh_,
		    r_dfdhref_,
		    r_dfdwvl_,
		    r_dfdasa_);
    
    dopplers(i)=Uvar(r_f_,"Hz");
  }
}


 void DopplerCentroid::getDopDerivativeYPAZ(Uvar& dfdy, Uvar& dfdp, Uvar& dfdasa)
{
  dfdy = Uvar(r_dfdy_,"Hz/rad")*rscale_;
  dfdp = Uvar(r_dfdp_,"Hz/rad")*rscale_;
  dfdasa = Uvar(r_dfdasa_,"Hz/rad")*rscale_;
}
//----------------------------------
//Private function declarations
//
//----------------------------------


//----------------------------------
//Update with every new data set
//------------------------------------
 void DopplerCentroid::computeTRMatrix(const StateVector& sc_state,
				const unsigned int& beam_id)
  {
    
    t_ypaz_ =sc_state.time();
    //ftcn(or SCH) coordinate definition
    //compute yaw and pitch angles
    DirectionVector sc_z(" ",sc_state.position().representIn(ftarget_));
    DirectionVector sc_y(" ",cross(sc_z,
				   sc_state.velocity().representIn(ftarget_)));
    DirectionVector sc_x(" ",cross(sc_y,sc_z));
    Frame ftcn("tnc",sc_state.position(),sc_x,sc_y,sc_z);//ftcn coordinate

    //SPACE CRAFT FRAME
    DirectionVector x("",fsc_,t_ypaz_,1,0,0);
    DirectionVector y("",fsc_,t_ypaz_,0,1,0);
    DirectionVector z("",fsc_,t_ypaz_,0,0,1);
    //When spacecraft rotated w.r.t. z-axis by 180 degree
    double z_rotation = dot(sc_x, x);
    if(z_rotation < 0)
      {
	x = -x;
	y= -y;
      }
    Frame FSC("new sc frame",sc_state.position(),x,y,z);
    SpiceDouble fsc_to_ftcn[3][3];
    SpiceDouble et;
    t_ypaz_.getEt(et);
    Frame::rotationMatrix(FSC,ftcn,et,fsc_to_ftcn);//from fsc to ftcn:3-2-1
    SpiceDouble* Spice_angle;
    Spice_angle = new SpiceDouble[3];    
    m2eul_c(fsc_to_ftcn,3,2,1,Spice_angle, Spice_angle+1,Spice_angle+2);
    //keep yaw and pitch
    yaw_rad_(beam_id)= *Spice_angle;
    pitch_rad_(beam_id)=*(Spice_angle+1);
    roll_rad_(beam_id)= *(Spice_angle+2);
    delete[] Spice_angle;

   
   
    position_(beam_id,0)=0;
    position_(beam_id,1)=0;
    position_(beam_id,2)= sc_state.position().magnitude().getInUnits("m")
      -radius_;//altitude
    
    //velocity in tcn
    FloatVector v=sc_state.velocity().representIn(ftcn);
    
    velocity_(beam_id,0)=dot(v,sc_x).getInUnits("m/s");
    velocity_(beam_id,1)=dot(v,sc_y).getInUnits("m/s");
    velocity_(beam_id,2)=dot(v,sc_z).getInUnits("m/s");
    
    //debug
    //in TCN, v_y component should be close to 0
    if(velocity_(beam_id,1) > 0.1){   
      cout<<"direction x "<< sc_x <<endl;
      cout<<"direction y "<< sc_y<<endl;
      cout<<"direction z "<<sc_z<<endl;
      cout<<"beam number "<< beam_id+1<<endl;
      cout<<"velocity in ftarget "<<v.representIn(ftarget_)<<endl;
      cout<<"dot (v, sc_x ) "<< dot(v,sc_x)<<" "<< dot(v,sc_y)<<" "<<dot(v,sc_z)<<endl;
      cout<<"velocity  in tcn"<< velocity_.getRow(beam_id)<<endl;
      cout<<"position "<<position_.getRow(beam_id)<<endl;
      cout<<"yaw pitch roll "<< yaw_rad_(beam_id)*radtodeg<< " "<< pitch_rad_(beam_id)*radtodeg<<" "<< roll_rad_(beam_id)*radtodeg<<endl;
      cout<<"time "<<t_ypaz_.utc("ISOD")<<endl;
      ErrorMessage("DopplerCentroid.cpp:: v in TCN has non-zero y component ").throwMe();
    }
  }

void
DopplerCentroid::displayFitting()
{

  vector<double> x;
  for(unsigned int i=0;i<range_.size();++i) x.push_back(range_[i]/1000.0);

  Plot a;
  a.addXY(x,doppler_,line("none"),sym("circle","red",0.25));
  a.addXY(x,doppler_ck_,line("none"),sym("circle","black",1));
  a.addXY(x,new_doppler_,line("none"),sym("circle","green",2));
  a.setTitle("r-data, black-ckernel, green-fit ");
  a.show("x");
  Plot b;
  b.addXY(x,diff_doppler_,line("none"),sym("circle","black",1));
  b.addLegend("residual doppler");
  b.show("x");
}


void
DopplerCentroid::reset()
{
  //at the end of computation, clear the container
  yaw_rad_=0.0;
  pitch_rad_=0.0;
  roll_rad_=0.0;
  position_=0.0;
  velocity_=0.0;
  beam_record_indicator_ = 0;
  beam_first_range_=0.0;
  beam_delta_range_=0.0;

  record_beam_id_.clear();
  record_time_.clear();
  range_.clear();
  doppler_.clear();
  doppler_ck_.clear();
  new_doppler_.clear();
  diff_doppler_.clear();
 
  r_vecin_=0.0;
  r_vobs_=0.0;
  r_a_=0.0;
  r_at2_=0.0;
  r_u_=0.0;
  r_v_=0.0;
  r_w_=0.0;
  r_chisq_=0.0;
  r_cvm_=0.0;
  wavelength_set_ =false;
}



void track_rel_int_doppler(const vector<unsigned int>& sab_record,
			   const vector<Uvar>& prf_record,
			   const vector<double>& fract_record,
			   const vector<Uvar>& geom_doppler,
			   const vector<int>& geom_int_doppler,
			   vector<int>& rel_int_doppler,
			   const bool& show_plot)
{
  
  try{
    //size check
    if(prf_record.size() != sab_record.size() ||
       prf_record.size() != fract_record.size() ||
       prf_record.size() != geom_doppler.size() ||
       prf_record.size() != geom_int_doppler.size() )
      ErrorMessage("size mismatch").throwMe();

 
    //count relative int
    rel_int_doppler.push_back(0);//first value


    //cycled variables
    int int_prf_summed=0;
    //double look_backward, look_forward;
   
    Uvar delta_int;
    //temporary containers
    vector<Uvar> diff_int_geom, diff_int_meas;
    diff_int_geom.clear();
    diff_int_meas.clear();
    diff_int_geom.push_back(0);
    diff_int_meas.push_back(0);


    for(unsigned int ii=1;ii<sab_record.size();++ii){
      diff_int_geom.push_back(geom_int_doppler[ii]-geom_int_doppler[ii-1]);
      
      delta_int = -(fract_record[ii] * prf_record[ii] - fract_record[ii-1]*prf_record[ii-1]);
      delta_int -= geom_int_doppler[ii]*(prf_record[ii]-prf_record[ii-1]);
      delta_int += geom_doppler[ii]- geom_doppler[ii-1];

      delta_int /= prf_record[ii-1];

      int_prf_summed = int_prf_summed + int(round_double(delta_int.getInUnits("")));
      rel_int_doppler.push_back(int_prf_summed);
    }
 

    if(sab_record.size() != rel_int_doppler.size())
      ErrorMessage("size mismatch").throwMe();


    if(show_plot){
      
      
      for(unsigned int ii=1;ii<sab_record.size();++ii)
	diff_int_meas.push_back(rel_int_doppler[ii]-rel_int_doppler[ii-1]);
      
      
      //plot input record
      vector<Uvar> x_axis; 
      vector<Uvar> y_axis;
      for(unsigned int i=0;i<sab_record.size();++i)
	x_axis.push_back(sab_record[i]);
      
      Plot a;
      a.addXY(x_axis,"",prf_record,"KHz",line("none"),sym("circle","red",1));
      a.setTitle("sab record vs prf ");
      a.show("x");
      
      for(unsigned int i=0;i<fract_record.size();++i)
	y_axis.push_back(fract_record[i]);
      Plot b;
      b.addXY(x_axis,"",y_axis,"",line("none"),sym("circle","red",1));
      b.setTitle("sab record vs measured fract ");
      b.show("x");
      
      Plot c;
      c.addXY(x_axis,"",geom_doppler,"KHz",line("none"),sym("circle","red",1));
      c.setTitle("sab record vs predicted boresight doppler ");
      c.show("x");
      
      y_axis.clear();
      for(unsigned int i=0;i<geom_int_doppler.size();++i)
	y_axis.push_back(geom_int_doppler[i]);
      Plot d;
      d.addXY(x_axis,"",y_axis,"",line("none"),sym("circle","red",1));
      d.setTitle("sab record vs predicted int doppler ");
      d.show("x");
      
      Plot e;
      e.addXY(x_axis,"",diff_int_geom,"",line("none"),sym("circle","red",0.5));
      e.addXY(x_axis,"",diff_int_meas,"",line("none"),sym("square","blue",1));
      e.setTitle("diff int geom-red, meas-blue");
      e.show("x");
    }
  }
  catch(ErrorMessage& e){
    cerr<<"Error: "<<e.msg<<endl;
  }
}





