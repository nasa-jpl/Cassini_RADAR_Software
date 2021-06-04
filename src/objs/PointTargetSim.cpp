//----------------------------------------------------------------------------
// PointTargetSim.cpp
//
// This file contains method definitions for the PointTargetSim class.
//----------------------------------------------------------------------------//----------------------
// Configuration Control
//----------------------

static const char rcs_id_cassini_sim_c[] =
  "@(#) $Id: PointTargetSim.cpp,v 11.7 2012/08/15 22:15:27 bstiles Exp $";

#include "PointTargetSim.h"
#include "Constants.h"
#include "config_keywords.h"
#include "Units.h"
#include "Array.h"
#include "SARFunctions.h"
#include "Distributions.h"
#include "DebugInfo.h"
#include<iostream>
#include<time.h>

//#define DEBUG

using std::cout;
using std::endl;
using std::cerr;




//---------------
// Constructors
//---------------

PointTargetSim::PointTargetSim(Config& cfg, ParamSourceE p)
  : CassiniSim(cfg), target_beyond_horizon_(false), ideal_tracking_on_(false),
    compute_echo_amplitude_(false), 
    simulate_backscatter_(false),
    simulate_area_(false),
    area_is_relative_(false),
    ideal_gain_setting_(false),
    use_config_bpd_(false),
    l1b_template(NULL),l1b(NULL), ieb_file(NULL),
    beam_order_idx(0), target_lat("target_lat",1), 
    target_lon("target_lon",1), target_height("target_height",1),
    target_ampl("target_ampl",1),
    target_type_id(0), target_azim("target_azim",1),
    target_elev("target_elev",1), target_beam("target_beam",1),
    target_alt("target_alt",1), 
    target_inbound("target_inbound",1), 
    target_special("target_special",1), radar_param_source(p),
    flyby(cfg), burst_no(0), data_take_id(0),
    check_this_burst(false), perform_single_burst_check(true)
  {
  config(cfg);
  }
  
PointTargetSim::~PointTargetSim(){
  if(l1b_template){
    delete l1b_template;
    l1b_template=NULL;
  }
  if(l1b){
    delete l1b;
    l1b=NULL;
  }
  if(ieb_file){
    delete ieb_file;
    ieb_file=NULL;
  }
  diag_file.close();
}
void PointTargetSim::config(Config& cfg){
  double sc=get_time(); // needed for manual profiling

  DebugInfo dbg("PointTargetSim::config");

  string diag_filename = cfg.str("PTS_DIAG_FILENAME");
  diag_file.open(diag_filename.c_str());

  //----------------------------------------------
  // Open and read IebFile if necessary
  //----------------------------------------------
  if(radar_param_source==IEB){
    string ieb_filename = cfg.str("ieb_file");
    ieb_delta_trigger_time = cfg["IEB_delta_trigger"];
    if(ieb_file) delete ieb_file;
    ieb_file=new IebList(cfg);
    ieb_file->readAllIebRecords(ieb_filename,"r");
  } 

  //----------------------------------------------
  // Open and read Input L1BFile if necessary
  //----------------------------------------------
  else if(radar_param_source==L1BFILE){
    string l1a_a_filename = cfg.str(L1A_PASSIVE_MODE_FILENAME);
    if(l1b_template) delete l1b_template;
    l1b_template = new L1B(l1a_a_filename,"rb","passive");
    l1b_template->readHeader();
 
  }


  //-----------------------------------------------
  // CONFIGURE START AND END TIMES
  //-----------------------------------------------

  int time_config_method = cfg.getInt("SIM_TIME_CONFIG_METHOD");
  string start_time_string, end_time_string;
  if(time_config_method==0){
    ErrorMessage e("PointTargetSim:: Time config method 0 is obsolete");
    e.throwMe();
  }
  switch(time_config_method){
  case 1: // Use times referenced to the epoch from Flyby
    start_time=flyby.startTime()+flyby.epochTime();
    end_time=flyby.endTime()+flyby.epochTime();
    break;
  case 2: // use ieb or l1bfile times
    if(radar_param_source == IEB){
      start_time=ieb_file->getFirstValidSlowFastFieldTime()+Uvar(5,"s");
      end_time=ieb_file->getEndTime()-Uvar(5,"s");
      t_=start_time; // needed to get parameters from first IEB
    }
    else if(radar_param_source == L1BFILE){
      // gives a little margin to ascertain that all records fall within
      // the range
      start_time=l1b_template->getStartTime()-Uvar(2,"s");
      end_time=l1b_template->getEndTime()+Uvar(2,"s");
    }
    else{
      ErrorMessage e("Error::PointTargetSim Cannot Configure Times");
      e.throwMe();
    }

    break;
  default:
    ErrorMessage e("Error::PointTargetSim Cannot Configure Times");
    e.throwMe();
  }

  //-----------------------------------------------------------
  // For IEB case make sure start and end times are within IEB
  //-----------------------------------------------------------
  if(radar_param_source == IEB){
    Time ieb_start_time=ieb_file->getFirstValidSlowFastFieldTime();
    Time ieb_end_time=ieb_file->getEndTime();
    
    if(ieb_end_time<end_time){
	cerr << "Warning end_time beyond end of IEB, truncating ..." << endl;
	end_time=ieb_end_time-Uvar(5,"s");
    }
    if(ieb_start_time>start_time){
      cerr 
	<< "Warning start_time before first slow/fast field, truncating ..." 
	<< endl;
      start_time=ieb_start_time+Uvar(5,"s");
    }
  }

  //----------------------------------------------------------------
  // CONFIGURE SIMULATION OPTIONS
  //----------------------------------------------------------------

  // configure whether or not (1/0) to compute echo amplitude
  // if not amplitude is automatically set to 1.0.

  target_info_in_diag_file=(bool)cfg.getInt("PTSIM_PRODUCE_FULL_DIAG_FILE");
  compute_echo_amplitude_=(bool)cfg.getInt("COMPUTE_ECHO_AMPLITUDE");

  use_random_phase_shift_=(bool)cfg.getInt("USE_RANDOM_PHASE_SHIFT");

  // This is only used for point target grids !!!
  // No need for random shifting of individual point targets
  // and such shifting would tend to screw up location checks
  use_random_location_shift_=(bool)cfg.getInt("USE_RANDOM_LOCATION_SHIFT");


  simulate_backscatter_=(bool)cfg.getInt("SIMULATE_BACKSCATTER");

  if( simulate_backscatter_ ){
    muhleman_k1_= cfg["Muhleman_backscatt_model_k1"].getInUnits("");
    muhleman_k2_= cfg["Muhleman_backscatt_model_k2"].getInUnits("");
  }

  // configure whether or not (1/0) to simulate thermal noise
  simulate_noise_=(bool)cfg.getInt("SIMULATE_THERMAL_NOISE");

  // configure whether or not (1/0) to simulate ADC quantization
  simulate_quantization_=(bool)cfg.getInt("SIMULATE_QUANTIZATION");

  // configure whether or not (1/0) to simulate baq
  simulate_baq_=(bool)cfg.getInt("SIMULATE_BAQ");

  // configure whether or not (1/0) to compute an ideal
  // attenuator  gain setting for each burst
  // if not gain is comnfigured
  ideal_gain_setting_=(bool)cfg.getInt("IDEAL_GAIN_SETTING");

  // configure whether or not to output
  // debugging info for a particular burst.
  // set burst number to check if desired.
  // Output goes to diag file.
  perform_single_burst_check=(bool)cfg.getInt("PERFORM_SINGLE_BURST_CHECK"); 
  if(perform_single_burst_check){
    check_burst_idx=cfg.getInt("CHECK_BURST_NO");
    if(check_burst_idx < 1){
      ErrorMessage 
	e("PointTargetSim:config bad Check Burst Number should be >= 1");
      e.throwMe();
    }
  }
  // configure number of samples to skip between delay/doppler/gain 
  // calculations
  sample_skip=cfg.getInt("SAMPLE_SKIP");

  // configure whether or not (1/0) to perform ideal doppler/range tracking
  ideal_tracking_on_=(bool)cfg.getInt("IDEAL_TRACKING_ON");

  // configure maximum angle to target for which an echo will be computed
  // This refers to the angle betwwen the boresight and the target
  // GAIN outside this angle is assumed to be zero.
  max_target_angle=cfg["MAX_TARGET_ANGLE"];
  
  //-------------------------------------------------------
  // Configure output file
  //-------------------------------------------------------

  string l1b_a_filename = cfg.str(L1B_ACTIVE_MODE_FILENAME);
  if(l1b) delete l1b;
  l1b = new L1B(l1b_a_filename,"wb","active");
  l1b->config(cfg);
  l1b->writeHeader();

  //--------------------------------------------
  // Read burst period if needed
  //--------------------------------------------
  use_config_bpd_=(bool)cfg.getInt("USE_CONFIG_BPD");
  if(use_config_bpd_ || radar_param_source==CONFIG){
    l1b->bpd=cfg["BURST_PERIOD"];
  }

  //-----------------------------------------------
  // Read in offset for anti-aliasing filter
  //-----------------------------------------------

  filter_centered=(bool)cfg.getInt("ASSUME_FILTER_CENTERED");

  if(!filter_centered){
    filter_start_freq_offset=cfg["FILTER_START_FREQ_OFFSET"];
  }
  //--------------------------------------------------------
  // Configure Point Target Location Info
  //--------------------------------------------------------

  // Determine name of target body
  target = cfg.str("target");

  // Determine number of point targets
  num_targets=cfg.getInt("NUM_TARGETS");

  // If num target is zero no point target configuration is necessary
  if(num_targets!=0){
    only_write_targeted_bursts=(bool)cfg.getInt("ONLY_WRITE_TARGETED_BURSTS");
    // Determine point target type
    // 0= LAT/LON
    // 1= S/C Altitude, beam number, azimuth and elevation angles
    //    beam number determines frame but all beams see all point targets
    // 2= Beam number and azimuth and elevation angles (every burst)
    //    each beam only sees point targets for that beam
    target_type_id = cfg.getInt("TARGET_TYPE");

    
    // determine amplitudes of point targets
    // bias to model if simulate_backscatter_ is set
    if(target_type_id<=2){
      target_ampl.resize(num_targets); 
 
      for (unsigned int i=0;i<num_targets;i++){
	  target_ampl(i)=cfg["POINT_TARGET_AMPLITUDE"+toStr(i+1)];
      }
      
    }
    
    if(target_type_id==6){
      target_ampl.resize(1);
      target_lat.resize(1);
      target_height.resize(1);
      target_lon.resize(1);
      target_special.resize(1);
      if(num_targets != 1){
	cerr << "Setting num_taregst to 1 for nadir target type id (6)" << endl;
      }
      num_targets=1;
    }
    double gain; // USED if TARGET_TYPE_ID =3
    
    // VARIABLES USED if TARGET_TYPE_ID =4,5
    int grid_size1, grid_size2; 
    Uvar grid_res1, grid_res2, tampl, tht;
    Uvar elev_min,azi_min,lat_min,lon_min;

    bool htused=false;
    bool image_file_used=cfg.keywordExists("POINT_TARGET_GRID_IMAGE_FILE");
    bool height_file_used=cfg.keywordExists("POINT_TARGET_GRID_HEIGHT_IMAGE_FILE");
    FILE* imagefp=NULL;
    string imfilename;

    FILE* htimfp=NULL;
    string htimfilename;
      
    switch(target_type_id){
      
    case 0: // Read target latitude and longitude if necessary
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_height.resize(num_targets);
      target_special.resize(num_targets);
      for (unsigned int i=0;i<num_targets;i++){
        htused=cfg.keywordExists("POINT_TARGET_HEIGHT"+toStr(i+1));
        if(htused){
	  target_height(i)=cfg["POINT_TARGET_HEIGHT"+toStr(i+1)];
	}
        else{
	  target_height(i)=Uvar(0,"km");
	}
	target_lat(i)=cfg["POINT_TARGET_LAT"+toStr(i+1)];
	target_lon(i)=cfg.convertWestLongitude("POINT_TARGET_LON"+toStr(i+1));
	target_special(i)=1; // all targets are debugged
         
      }
      break;
      
    case 1: // Read S/C altitude and position in beam and convert to lat/lon
      target_height.resize(num_targets);
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_alt.resize(num_targets);
      target_beam.resize(num_targets);    
      target_azim.resize(num_targets);    
      target_elev.resize(num_targets);    
      target_inbound.resize(num_targets);
      target_special.resize(num_targets);    
      for (unsigned int i=0;i<num_targets;i++){
        target_height(i)=Uvar(0,"km");
	target_alt(i)=cfg["POINT_TARGET_SC_ALT"+toStr(i+1)];
	target_beam(i)=cfg.getInt("POINT_TARGET_BEAM"+toStr(i+1));
	target_azim(i)=cfg["POINT_TARGET_AZIMUTH"+toStr(i+1)];
	target_elev(i)=cfg["POINT_TARGET_ELEVATION"+toStr(i+1)];
	target_inbound(i)=cfg.getInt("POINT_TARGET_INBOUND"+toStr(i+1));
	altBAEToLatLon(i);
	target_special(i)=1; // all targets are debugged
      }
      break;
      
    case 2: // Read Beam number, Azimuth, and Elevation angles
      target_height.resize(num_targets);
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_beam.resize(num_targets);
      target_azim.resize(num_targets);    
      target_elev.resize(num_targets);    
      target_special.resize(num_targets); 
   
      for (unsigned int i=0;i<num_targets;i++){
	target_height(i)=Uvar(0,"km");
	target_beam(i)=cfg.getInt("POINT_TARGET_BEAM"+toStr(i+1));
	target_azim(i)=cfg["POINT_TARGET_AZIMUTH"+toStr(i+1)];
	target_elev(i)=cfg["POINT_TARGET_ELEVATION"+toStr(i+1)];
	target_special(i)=1; // All targets are debugged
      }
      break;
    case 3: // Set Point Targets to 3dB two-way extrema of each beam
      num_targets=25;
      target_height.resize(num_targets);
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_beam.resize(num_targets);
      target_azim.resize(num_targets);    
      target_elev.resize(num_targets); 
      target_ampl.resize(num_targets);
      target_special.resize(num_targets);

  
      // all amplitudes are the same
      for (unsigned int i=0;i<num_targets;i++){
	target_height(i)=Uvar(0,"km");
	target_ampl(i)=cfg["POINT_TARGET_AMPLITUDE1"];
      }
      
      // read in two-way gain contour for the point targets
      gain=cfg.getDouble("TARGET_TWO_WAY_GAIN_DB");

      // assign azimuth and elevations and beam numbers
      for (unsigned int i=0;i<5;i++){
	int tn=i*5;
	Uvec azim=beam_[i].elevAndGainToAzim(0,-gain/2);
	Uvec elev=beam_[i].azimAndGainToElev(0,-gain/2);

	// boresight
	target_beam(tn)=i+1;
	target_azim(tn)=Uvar(0,"deg");
	target_elev(tn)=Uvar(0,"deg");
	target_special(tn)=1; // Only boresight targets are debugged
	tn++;

	//top
	target_beam(tn)=i+1;
	target_azim(tn)=Uvar(0,"deg");     
	target_elev(tn)=elev(1);
	target_special(tn)=0;
	tn++;

	//bottom
	target_beam(tn)=i+1;
	target_azim(tn)=Uvar(0,"deg");     
	target_elev(tn)=elev(0);
	target_special(tn)=0;
	tn++;

	//left
	target_beam(tn)=i+1;
	target_azim(tn)=azim(0);    
	target_elev(tn)=Uvar(0,"deg");   
	target_special(tn)=0;
	tn++;
	
	//right
	target_beam(tn)=i+1;
	target_azim(tn)=azim(1);    
	target_elev(tn)=Uvar(0,"deg"); 
	target_special(tn)=0; 
 
      }
      // reassign target type id it is 2 from here on out.
      target_type_id=2;
      break;

      // grid of point target in azimuth and elevation for all beams
    case 4:
      
      // set area calculation Booleans
      simulate_area_=true;
      area_is_relative_=true;

      // get grid size 
      grid_size1= cfg.getInt("POINT_TARGET_GRID_NUM_AZIM");
      grid_size2= cfg.getInt("POINT_TARGET_GRID_NUM_ELEV");
      
      // get grid resolution
      grid_res1 = cfg["POINT_TARGET_GRID_AZIM_RES"];
      grid_res2 = cfg["POINT_TARGET_GRID_ELEV_RES"];
      grid_azim_res=grid_res1;
      grid_elev_res=grid_res2;

      // get minimum azimuth and elevation
      azi_min = cfg["POINT_TARGET_GRID_MIN_AZIMUTH"];
      elev_min = cfg["POINT_TARGET_GRID_MIN_ELEVATION"];

      // computer center of grid for later use.
      grid_azim_center=azi_min+grid_azim_res*grid_size1;
      grid_elev_center=elev_min+grid_elev_res*grid_size2;

      // compute number of point targets in grid
      num_targets=grid_size1*grid_size2*5;
      target_height.resize(num_targets);
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_special.resize(num_targets);
      target_beam.resize(num_targets);
      target_azim.resize(num_targets);    
      target_elev.resize(num_targets); 
      target_ampl.resize(num_targets);

      tampl=cfg["POINT_TARGET_AMPLITUDE1"];
      for (int b=0;b<5;b++){

	for(int a=0;a<grid_size1;a++){
	  for(int e=0;e<grid_size2;e++){
	    int tn=grid_size1*grid_size2*b+grid_size2*a+e;
	    target_height(tn)=Uvar(0,"km");	    
	    // set beam and position in beam for each point target
	    target_beam(tn)=b+1;
            Uvar atweak=0, etweak=0;
            if(use_random_location_shift_){
	      atweak=grid_res1*drand48()-grid_res1/2;
	      etweak=grid_res2*drand48()-grid_res2/2;
	    }
	    target_azim(tn)=azi_min+a*grid_res1+atweak;
	    target_elev(tn)=elev_min+e*grid_res2+etweak;
	    
	    // set amplitude bias

	    // add correction factor for portion of beam footprint "covered"
            // by each point target (assumes each azim/elev pixel has equal 
	    // area.
            target_ampl(tn)=tampl;

            // backscatter if desired) and footprint area are computed later;
	    
            // special debugging flags are set to 0 for every target
	    // except the one in the center of the grid            
            if(a==grid_size1/2 && e==grid_size2/2){
	      target_special(tn)=1;
	    }
            else target_special(tn)=0;
	  }
	}
      }
      
      // reassign target type id it is 2 from here on out.
      target_type_id=2;
      break;
      
      // grid of point targets in latitude and longitude
    case 5:


      // set up area simulation booleans
      simulate_area_=true;
      area_is_relative_=false;

      // get grid size 
      grid_size1= cfg.getInt("POINT_TARGET_GRID_NUM_LAT");
      grid_size2= cfg.getInt("POINT_TARGET_GRID_NUM_LON");
      
      // get grid resolution
      grid_res1 = cfg["POINT_TARGET_GRID_LAT_RES"];
      grid_res2 = cfg["POINT_TARGET_GRID_LON_RES"];
      
      
      // get minimum azimuth and elevation
      lat_min = cfg["POINT_TARGET_GRID_MIN_LATITUDE"];
      lon_min = cfg.convertWestLongitude("POINT_TARGET_GRID_WESTERNMOST_LONGITUDE");
      
      // compute number of point targets in grid
      num_targets=grid_size1*grid_size2;
      target_height.resize(num_targets);
      target_lat.resize(num_targets);
      target_lon.resize(num_targets);
      target_ampl.resize(num_targets);
      target_special.resize(num_targets);

      tampl=cfg["POINT_TARGET_AMPLITUDE1"];
      htused=cfg.keywordExists("POINT_TARGET_HEIGHT1");
      if(htused)
        tht=cfg["POINT_TARGET_HEIGHT1"];
      else
	tht=Uvar(0,"km");
      if(image_file_used){
	imfilename=cfg.str("POINT_TARGET_GRID_IMAGE_FILE");
	imagefp=fopen(imfilename,"r");
      }


      if(height_file_used){
	htimfilename=cfg.str("POINT_TARGET_GRID_HEIGHT_IMAGE_FILE");
	htimfp=fopen(htimfilename,"r");
      }

      for(int a=0;a<grid_size1;a++){
        // debug hack
        printf("% d of %d point target rows configured\n",a,grid_size1);
	for(int e=0;e<grid_size2;e++){
	  int tn=grid_size2*a+e;
	  Uvar lattweak=0, lontweak=0;
	  if(use_random_location_shift_){
	    lattweak=grid_res1*drand48()-grid_res1/2;
	    lontweak=grid_res2*drand48()-grid_res2/2;
	  }
	  // compute position for each point target
	  target_lat(tn)=lat_min+a*grid_res1;
	  target_lon(tn)=lon_min+e*grid_res2;
	  
	  // set amplitudes identically
          // unless image file is used then assign value from file
	  target_ampl(tn)=tampl;
          target_height(tn)=tht;
	  if(image_file_used){
	    float val=0.0;
	    if(!fread(&val,sizeof(float),1,imagefp)){
	      ErrorMessage e("Cannot read point target grid image file "+imfilename);
	      e.throwMe();
	    }
	    target_ampl(tn)*=val;
	  }

	  if(height_file_used){
	    float val=0.0;
            // floating point values in fiule are in meters
	    if(!fread(&val,sizeof(float),1,htimfp)){
	      ErrorMessage e("Cannot read point target grid heightimage file "+htimfilename);
	      e.throwMe();
	    }
	    target_height(tn)=Uvar(val,"m");
	  }

	  // add area correction factor 
	  // compute corners of grid cell surrounding point target
	  // since we are only interested in the area
	  // incidentals like aspherical body shape or longitude sign 
	  // convention are negligible
          // the effect of surface height is also ignored
        
	  Uvar pixel_area=getLatLonPixelAreaOnSphere(target_lat(tn),
						       target_lon(tn),
						       grid_res1,
						       grid_res2);
	  target_ampl(tn)*=pixel_area;

          target_lat(tn)+=lattweak;
          target_lon(tn)+=lontweak;

	  // special debugging flags are set to 0 for every target
	  // except the one in the center of the grid            
	  if(a==grid_size1/2 && e==grid_size2/2){
	    target_special(tn)=1;
	  }
	  else target_special(tn)=0;
	}
      }
      if(image_file_used) fclose(imagefp);
      if(height_file_used) fclose(htimfp);
     
      // reassign target type id it is 0 from here on out.
      target_type_id=0;
      break;
    case 6: // single point target at nadir for every burst
      target_special(0)=1;
      break;
    default:
      ErrorMessage e("Error:PointTargetSim:BadTargetType");
      e.throwMe();
    }
  }
  // set target_type_id to 0 for num_targets==0 case
  else{
    target_type_id=0;
  }

  // convert target lat lons from Planetographic to Planetodetic
  // for target_type_id =0 (or 5 which has already been reset to 0 by  now)
  
  double radx=default_target_radii.km(PositionVector::X);
  double rady=default_target_radii.km(PositionVector::Y);
  double radz=default_target_radii.km(PositionVector::Z);
  double rad=sqrt(radx*radx+rady*rady+radz*radz)/sqrt(3);
  if(target_type_id==0){
    for(unsigned int tt=0;tt<num_targets;tt++){
      double lat_centric=get_in_base_units(target_lat(tt));
      double lon_centric=get_in_base_units(target_lon(tt));
      double zc=sin(lat_centric)*rad;
      double xc=cos(lon_centric)*cos(lat_centric)*rad;
      double yc=sin(lon_centric)*cos(lat_centric)*rad;
      double denom=xc*xc/(radx*radx)+yc*yc/(rady*rady)+zc*zc/(radz*radz);
      double k=sqrt(1/denom);
      double x=k*xc;
      double y=k*yc;
      double z=k*zc;
      double lat=asin(z/radz);
      double lon=atan2(y/rady,x/radx);
      target_lat(tt)=Uvar(lat,"radians");
      target_lon(tt)=Uvar(lon,"radians");
    }
  }
   
  //--------------------------------------------------------
  // SET UP RADAR CONSTANTS
  //--------------------------------------------------------

  // configure SAR calibration coefficient
  // need to change so it can access calibration coefficients for other modes
  // and change the code that uses it so it can distinguish between them.
  sar_cal_coeff=cfg["squared_deviation_of_system_noise_input_at_SARL"];

  // configure SAR system noise temperature
  Tsys=cfg["Tsys"];

  // configure shift of down converted signal 
  // i.e., the frequency of the down converted return at input to the ADC
  // if there was no chirp or doppler shift
  updown_shift=-cfg["stalo_frequency"];


  // configure transmit power
  transmit_power=cfg["Pt"];
  

  //-----------------------------------------------------------------
  // Complete configuration of parameters from appropriate source
  //-----------------------------------------------------------------
  
  // Source=CONFIG
  if(radar_param_source == CONFIG)  getRadarParams(cfg);

  // Source=L1BFILE

  //  if Source=IEB or L1BFILE No need to do anything Radar Parameters 
  //  are set for each burst

  // handle special debug settings 
  DebugInfo dbgsp_radar_params("Special_Radar_Parameters_Debug");
  if(dbgsp_radar_params.level){
    if(dbgsp_radar_params.level & 0x01) special_debug_pulse_duty=cfg["INTRA_BURST_DUTY_CYCLE"];
    if(dbgsp_radar_params.level & 0x02) special_debug_attr_gain_bias_db=cfg.getInt("DEBUG_ATTN_GAIN_BIAS_IN_DB");
    if(dbgsp_radar_params.level & 0x04) special_tro_in_pri=cfg.getInt("RX_WINDOW_EXTRA_PRI");
  }

  // needed for manual profiling
  double se=get_time();
  pts_config_time+=se-sc;

  if(dbg.level){
    Time et=flyby.epochTime();
    dbg.file<<"Debugging PointTargetSim::config start_time " << start_time -et
	    <<" end_time " << end_time -et << endl;
  }
}

//----------------------------------------------
// Routine uses a combination of IEB info in a L1A file and/or
// (Eventually RMSS output should be used)
// and config file info to perform a point target simulation
// writes out to a L1B file.
//---------------------------------------------------

void PointTargetSim::run(Config& cfg){
  // handle special debug settings 
  DebugInfo dbgsp_radar_params("Special_Radar_Parameters_Debug");
  DebugInfo dbgsp("Special_Calibration_Debug");

  // Write comments to diag file
  diag_file.comment("########################################");
  diag_file.comment("Burst number indexed arrays");
  diag_file.comment("t:= minutes from closest approach");
  diag_file.comment("beam_no: beam number (1-5)");
  diag_file.comment("");
  diag_file.comment("-------------------------------------------------");
  diag_file.comment("The following fields are relevant for TARGET_TYPE=3");
  diag_file.comment("and an appropriate TARGET_TWO_GAIN_DB value ....");
  diag_file.comment("start_tm: Time margin (s) at beginning of Rx window");
  diag_file.comment("end_tm: Time margin (s) at end of Rx window");
  diag_file.comment("start_fm: Freq margin (Hz) at start of analog filter");
  diag_file.comment("end_fm: Freq margin (Hz) at end of analog filter");
  diag_file.comment("start_pm: Margin in PRI beginning of Rx window");
  diag_file.comment("start_pm: Margin in PRI  at end of Rx window");
  diag_file.comment("(negative start_pm and end_pm indicate missed pulses)");
  diag_file.comment("-------------------------------------------------");
  diag_file.comment("pulse_duty: Pulse (Short term) duty cycle percentage");
  diag_file.comment("burst_duty: Burst (Long term) duty cycle percentage");
  diag_file.comment("tro: Number of extra PRI's in receive window.");
  diag_file.comment("pul: Number of pulses transmitted");
  diag_file.comment("rwd: Receive window delay in s.");
  diag_file.comment("scalt: S/C altitude at time t in km");
  diag_file.comment("min_dop: Min doppler (Hz) for all valid point targets");
  diag_file.comment("max_dop: Max doppler (Hz) for all valid point targets");
  diag_file.comment("min_rtt: Min round trip (s) for all valid point targets");
  diag_file.comment("max_rtt: Max round trip (s)for all valid point targets");
  diag_file.comment("rms0: RMS of signal prior to adding thermal noise");
  diag_file.comment("min0: Min signal value prior to adding thermal noise");
  diag_file.comment("max0: Max signal value prior to adding thermal noise");
  diag_file.comment("rms1: RMS of signal prior to quantization/truncation");
  diag_file.comment("min1: Min signal value prior to quant/trunc");
  diag_file.comment("max1: Max signal value prior to quant trunc");
  diag_file.comment("rms2: RMS of signal prior to baq encode/decode");
  diag_file.comment("min2: Min signal value prior to baq encode/decode");
  diag_file.comment("max2: Max signal value prior to baq encode decode");
  diag_file.comment("rms3: RMS of signal after baq encode/decode");
  diag_file.comment("min3: Min signal value after baq encode/decode");
  diag_file.comment("max3: Max signal value after baq encode decode");
  diag_file.comment("##############################################");
  diag_file.comment("");
  diag_file.comment("##############################################");
  diag_file.comment("Target Number indexed arrays                  ");
  diag_file.comment("Includes data from most recent burst in which");
  diag_file.comment("the target is valid.");
  diag_file.comment("targetang: Angle between boresight and target");
  diag_file.comment("fdop: doppler frequency in Hz");
  diag_file.comment("range: range to target in km");
  diag_file.comment("height: surface height of target above ref sphere in m");
  diag_file.comment("scale: Amplitude of echo from point target");
  diag_file.comment("lat: geodetic target latitude in degrees");
  diag_file.comment("lon: geodetic target longitude in degrees");
  diag_file.comment("##############################################");
  diag_file.comment("");
     

  // manual profiling
  double run_clock_start=get_time();

  //--------------------------
  //Need beam frame for measurement geometry calculation for each burst
  //------------------------
  Umat azim_1way3dB_ellipse_fit("azimuth ellipse fit",5,4);
  Umat elev_1way3dB_ellipse_fit("azimuth ellipse fit",5,4);
  Umat azim_2way3dB_ellipse_fit("azimuth ellipse fit",5,4);
  Umat elev_2way3dB_ellipse_fit("azimuth ellipse fit",5,4);
  azim_1way3dB_ellipse_fit=Uvar(0,"rad");
  elev_1way3dB_ellipse_fit=Uvar(0,"rad");
  azim_2way3dB_ellipse_fit=Uvar(0,"rad");
  elev_2way3dB_ellipse_fit=Uvar(0,"rad");

  //make all 5 beams and store
  for (unsigned int i = 0; i < 5;i++)
    {
      Uvec azi1_fit("azi fit",4), elev1_fit(" ",4);    
      beam_[i].computeBestFitEllipse(azi1_fit,elev1_fit,-3.0);//one way 3dB  
      for(unsigned int j=0;j<4;++j)
	{
	  azim_1way3dB_ellipse_fit(i,j)= azi1_fit(j);
	  elev_1way3dB_ellipse_fit(i,j)=elev1_fit(j);
	}    
      Uvec azi2_fit("azi fit",4), elev2_fit(" ",4);   
      beam_[i].computeBestFitEllipse(azi2_fit,elev2_fit,-1.5);//two way 3dB   
      for(unsigned int j=0;j<4;++j)
	{
	  azim_2way3dB_ellipse_fit(i,j)= azi2_fit(j);
	  elev_2way3dB_ellipse_fit(i,j)=elev2_fit(j);
	}
    }
 

  // check for L1B NULL (closed L1B file) most likely to occur
  // if run is called twice.
  if(l1b==NULL){
    ErrorMessage e("Cannot run PointTargetSim twice");
  }

 
  // set Time to be the start_time read from the config file

  setTime(start_time);
    
  // If L1B_FILE skip to start time in input file
  // and set up ending and beginning record numbers
  
  unsigned int start_record_num=0,end_record_num=0, input_record_no=0;
  if(radar_param_source==L1BFILE){
    start_record_num=
      l1b_template->sclkToRecord(start_time.sclk("Cassini"));
    end_record_num=
      l1b_template->sclkToRecord(end_time.sclk("Cassini"));
    if(start_record_num>=l1b_template->recordCount() || 
       end_record_num == 0){
      ErrorMessage e("PointTargetSim::Times out of range for Input L1B file");
      e.throwMe();
    }
    l1b_template->skipRecord(start_record_num);
    input_record_no=start_record_num;
  }


  // loop through times Burst Period is extracted from L1A file and
  // used to increment time.
  // Loop continues until end_time (read from config file) is reached


  while(t_<=end_time){
    // read next input L1B record break if end record_num reached
    if(radar_param_source==L1BFILE){
      if(input_record_no>=end_record_num) break;
      l1b->copyPassive(*l1b_template);
      setTime(l1b->t+l1b->transmit_time_offset); 
      input_record_no++;

      // if the L1A-L1B conversion fails due to bad
      // geometry  throw an error
      if(!l1b->goodGeometry()){
	char msg[60];
	sprintf(msg,"Bad L1A-L1B convert: quality_flag=%d",(int)l1b->quality_flag);
	ErrorMessage e(msg);
	e.throwMe();
      }
      // configure starting frequency of anti-alias filter
      // this is read from the config file until the precise formulation
      // is determined we assume this is -adc/2
      if(filter_centered){
	filter_start_freq=-l1b->adc/2+(l1b->adc/2-l1b->rc_bw)/2;
      }
      else{
	filter_start_freq=-l1b->adc/2.0+filter_start_freq_offset;
      }

      // close but no cigar
      // chirp_rate=l1b->csq*l1b->slow_cfs/l1b->csd/(1+l1b->csq);
      setChirpRate();

      if(!l1b->isCal()) updateRadarParams();

      // fake out num bursts in flight greater than 1
      // Currently this only works for L1BFILE and not IEB param source mode
      
      if(l1b->num_bursts_in_flight>1){
	l1b->rwd+=(l1b->num_bursts_in_flight-1)*l1b->bpd;
	l1b->num_bursts_in_flight=1;
      }
    }

    else if (radar_param_source == CONFIG){
      updateRadarParams();
    }

    // IEB
    else{
      getRadarParamsFromIeb();
      if(!l1b->isCal()) updateRadarParams();
    }
    
    // case for special radar parameter debug settings
    if(dbgsp_radar_params.level){
       int bit[32];
       int flag=dbgsp_radar_params.level;
       for(int i=0;i<3;i++){
	 bit[i]=flag & 0x00001;
	 flag = flag >> 1;
       }
       if(bit[0]){
         // change duty cycle but keep chirp bandwidth the same!!!
         // also keep chirp step duration the same!
         Uvar chirpbw=(l1b->csq)*l1b->slow_cfs;
	 l1b->chirp_length = l1b->pri * special_debug_pulse_duty;
	 l1b->csq=(unsigned int)
	   floor(get_in_base_units(l1b->chirp_length/l1b->csd)+0.5);
	 l1b->csq-=1;
	 l1b->slow_cfs=chirpbw/l1b->csq;
         setChirpRate();
       }
       if(bit[1]){
	 l1b->at1_db-=special_debug_attr_gain_bias_db;
	 l1b->at3_db-=special_debug_attr_gain_bias_db;
	 l1b->at4_db-=special_debug_attr_gain_bias_db;
	 l1b->at1 = cag_table_.gainLinearScale((int)l1b->at1_db);
	 l1b->at3 = cag_table_.gainLinearScale((int)l1b->at3_db);
	 l1b->at4 = cag_table_.gainLinearScale((int)l1b->at4_db);
	 // atx_each values are not set because the reverse mapping
         // nominal attenuation value -> mask is not implemented and
	 // is ambiguous, RMSS uses some table which is unavailable here.
       }
       if(bit[2]){
         Uvar old_tro=l1b->tro;
	 l1b->tro=special_tro_in_pri*l1b->pri;
	 l1b->ctrx=l1b->pul+special_tro_in_pri; 
	 Uvar receive_window_length=l1b->ctrx*l1b->pri;
	 Uvar num_samples=receive_window_length*l1b->adc;
	 double num_samps_int=floor(num_samples.getInUnits("")+0.5);

	 // Assumes pri and chirp_length do not need to be updated to
         // insure integer number of samples (no check if wrong)
	 //         if(num_samples.getInUnits("")!=num_samps_int){
	 //  ErrorMessage e("PointTargetSim:: special radar params case causes
	 //  Non integer number of samples.");
	 //  e.throwMe();
	 //  }

         // recompute rwd and Nradar_data to account for change in tro
	 l1b->rwd+=old_tro/2-l1b->tro/2;
         l1b->Nradar_data=(unsigned int)num_samps_int;
         l1b->checkNradar_data(); // make sure Nradar_data is in valid 	 
       }
    }
    // skip any calibration or radiometer only bursts
    if(l1b->isCal() || l1b->isPassive()){
      if(radar_param_source!=L1BFILE){
	setTime(t_+Time(l1b->bpd));
      }
      continue;
    }

    // store burst start time (t_ is going to be modified during processing)
    burst_start_time=t_;

 
    // Set up single burst checking capability
    if(perform_single_burst_check && (int)burst_no==check_burst_idx){
      check_this_burst=true;
    }
    else{
      check_this_burst=false;
    }

    // initialize echo
    // old data which is beyond the current value of Nradar_data is not
    // reset to zero, but that should not matter.
    for(unsigned int c=0;c<l1b->Nradar_data;c++) l1b->radar_data(c)=0;

    // initialize margin values to ridiculously large numbers for each burst
    start_time_margin_in_s=HUGE_VAL;
    end_time_margin_in_s=HUGE_VAL;
    start_freq_margin_in_Hz=HUGE_VAL;
    end_freq_margin_in_Hz=HUGE_VAL;

    // intialize round_trip_time min and max values to extreme numbers
    min_rtt_in_s=HUGE_VAL;
    max_rtt_in_s=-HUGE_VAL;
    min_dop_in_Hz=HUGE_VAL;
    max_dop_in_Hz=-HUGE_VAL;


    if(ideal_tracking_on_){
      idealDopplerRangeTrack();
    }

    computeDoubleParams();

    // For each point target

    
    // for(int c=0;c<10;c++) // HACK :: UNCOMMENT FOR PROFILING

    float average_echo_ampl=0.0;
    int num_nonzero_pt_echoes=0;
    bool target_found = false;
    for(unsigned int i=0; i<num_targets; i++){

      if(target_type_id==2) // point target specified by (beam azim elev)
	{
          // For TARGET TYPE 2 each beam only sees point targets with
          // that beam number
	  if(l1b->beam_number != (unsigned int)target_beam(i)) continue;

          // reference time for point target 
          // is mid-way between center of range window
          // and center of receive window,
          Time target_time = burst_start_time+
	    (l1b->rwd + l1b->ctrx*l1b->pri/2 + l1b->pul*l1b->pri/2)/2;

          // get target lat/lon
          timeBAEToLatLon(target_time,i);
	}
      else if(target_type_id==3){
       ErrorMessage e("Target Type ID =3 not set to 2 after initialization.");
       e.throwMe();
      }
      else if(target_type_id==6){
	// reference time for point target 
          // is mid-way between center of range window
          // and center of receive window,
          Time target_time = burst_start_time+
	    (l1b->rwd + l1b->ctrx*l1b->pri/2 + l1b->pul*l1b->pri/2)/2;
          setTime(target_time);
	  tg_.nadirLatLon(target_lat(i),target_lon(i));
          
      }
      bool target_in_range = checkTargetInRange(i);

      // nadir target assumed always in range!

      if(target_type_id == 6 ) target_in_range=true;
      if(target_in_range){
	target_found=true;
	float val=addEcho(i);
	if(val>0.0){
	  average_echo_ampl+=val;
          num_nonzero_pt_echoes++;
	}
      }
    }
    average_echo_ampl/=num_nonzero_pt_echoes;
    // special calibration debug outputs
    if(dbgsp.level){
      dbgsp.file << endl << "PointTargetSim::run  BurstNo:" << burst_no
		 << " BeamNo:" << l1b->beam_number
		 << " AverageEchoAmpl:"  << average_echo_ampl
		 << " NumNonzeroPointTargetEchoes:" << num_nonzero_pt_echoes
		 << endl << endl;
    }
    // Add Noise 
    unsigned int dummy=0;
    double rms0=l1b->radar_data.rms(l1b->Nradar_data);
    double max0=l1b->radar_data.max(dummy);
    double min0=l1b->radar_data.min(dummy);

    if (simulate_noise_) addNoise();

    double rms1=l1b->radar_data.rms(l1b->Nradar_data);
    double max1=l1b->radar_data.max(dummy);
    double min1=l1b->radar_data.min(dummy);


    // Multiply by attenuator gain AND IF
    // simulate_quantization_ is set 
    // convert to integer in range -128 to 127 
    convertToByte();

    double rms2=l1b->radar_data.rms(l1b->Nradar_data);
    double max2=l1b->radar_data.max(dummy);
    double min2=l1b->radar_data.min(dummy);

    // BAQ encode/decode 

    if (simulate_baq_) simulateBaq();


    double rms3=l1b->radar_data.rms(l1b->Nradar_data);
    double max3=l1b->radar_data.max(dummy);
    double min3=l1b->radar_data.min(dummy);
    

    // compute all measurement geometry before writing record
    l1b->computeMeasurementGeometry(azim_1way3dB_ellipse_fit,
				 elev_1way3dB_ellipse_fit,
				 azim_2way3dB_ellipse_fit,
				 elev_2way3dB_ellipse_fit);
    // write L1B record
    // If only_write_target_bursted is set
    // and no target was found this is skipped
    if(target_found || !only_write_targeted_bursts)
      l1b->writeRecord();

    // write outputs to diag_file


   // compute time of burst for output to diagnostic file
    float rel_bst_in_min
      =(burst_start_time-flyby.epochTime()).getInUnits("min");
    diag_file.set("t",rel_bst_in_min,burst_no);
    diag_file.set("beam_no",(float)l1b->beam_number,burst_no);
    diag_file.set("start_tm",start_time_margin_in_s,burst_no);
    diag_file.set("end_tm",end_time_margin_in_s,burst_no);
    diag_file.set("start_fm",start_freq_margin_in_Hz,burst_no);
    diag_file.set("end_fm",end_freq_margin_in_Hz,burst_no);
    float start_pri_margin=start_time_margin_in_s/pri_in_s;
    float end_pri_margin=end_time_margin_in_s/pri_in_s;
    diag_file.set("start_pm",start_pri_margin,burst_no);
    diag_file.set("end_pm",end_pri_margin,burst_no);
    float pulse_duty=100*chirp_length_in_s/pri_in_s;
    float bpd_in_s;
    if(radar_param_source!=IEB)
      bpd_in_s=l1b->bpd.getInUnits("s");
    else 
      bpd_in_s=ieb.getBpd().getInUnits("s");
    float burst_duty=100*l1b->pul*chirp_length_in_s/bpd_in_s;  
    diag_file.set("pulse_duty",pulse_duty,burst_no);
    diag_file.set("burst_duty",burst_duty,burst_no);     
    diag_file.set("tro",(l1b->tro/l1b->pri).getInUnits(""),burst_no);
    diag_file.set("pul",l1b->pul,burst_no);
    diag_file.set("rwd",l1b->rwd.getInUnits("s"),burst_no);
    diag_file.set("min_dop",min_dop_in_Hz,burst_no);
    diag_file.set("max_dop",max_dop_in_Hz,burst_no);
    diag_file.set("min_rtt",min_rtt_in_s,burst_no);
    diag_file.set("max_rtt",max_rtt_in_s,burst_no);
    diag_file.set("rms0",rms0, burst_no);
    diag_file.set("min0",min0, burst_no);
    diag_file.set("max0",max0, burst_no);
    diag_file.set("rms1",rms1, burst_no);
    diag_file.set("min1",min1, burst_no);
    diag_file.set("max1",max1, burst_no);
    diag_file.set("rms2",rms2, burst_no);
    diag_file.set("min2",min2, burst_no);
    diag_file.set("max2",max2, burst_no);
    diag_file.set("rms3",rms3, burst_no);
    diag_file.set("min3",min3, burst_no);
    diag_file.set("max3",max3, burst_no);
    setTime(burst_start_time);
    Uvar scalt=altitude();
    diag_file.set("scalt",scalt.getInUnits("km"),burst_no);
    // set next time
    if(radar_param_source!=L1BFILE){
       setTime(burst_start_time+Time(l1b->bpd));
    }
    
  }
  
  // rewrite L1B header
  l1b->rewriteHeader();

  // manual profiling
  double run_clock_end=get_time();
  pts_run_time+=run_clock_end-run_clock_start;
}


bool PointTargetSim::checkTargetInRange(const int target_idx){

  DebugInfo dbg("PointTargetSim::checkTargetInRange");
  // halfway between middle of receive and transmit windows
  Uvar midoffset = 0.5*(l1b->rwd+0.5*l1b->ctrx*l1b->pri+0.5*l1b->pul*l1b->pri);
  Time midpoint=burst_start_time+midoffset;
  setBeamTime(l1b->beam_number,midpoint);
  DirectionVector boresight_look=tg_.lookDirection();
  
  setBeamTimeLatLonHeight(l1b->beam_number,midpoint,
			  target_lat(target_idx),target_lon(target_idx),
			  target_height(target_idx));

  // check for beyond horizon case
  if(!tg_.foundSurfaceIntercept()) return(false);
 
  DirectionVector target_look=tg_.lookDirection(); 
  Uvar ang=target_look.angle(boresight_look);

  // Right now the diag_file output 
  // overwrites the previous burst with the same target_idx
  // This works well is each beam is seen once and target type is 2 or 3
  // but it needs to be made more general 
  if(target_info_in_diag_file){
    diag_file.set("targetang",ang.getInUnits("deg"),target_idx+1);
  }
  // Purposefully uses same reference time as SAR processor rather than
  // computing the more accurate values from getDelayDopplerAndScale
  setBeamTimeLatLonHeight(l1b->beam_number,midpoint,
			  target_lat(target_idx),target_lon(target_idx),
			  target_height(target_idx));
  
  StateVector st=tg_.state();
  PositionVector p = st.position();
  FloatVector v = st.velocity();
  //cout.precision(20);
  //cout << l1b->sab_counter << ": Position " << p << " Velocity"  << v  
  //     << " Time: " << st.time() << endl;
  Uvar range=tg_.range();
  fdop=tg_.doppler(lambda_chirp); // special relativistic Doppler
  fdop_in_Hz=get_in_base_units(fdop);

  PositionVector target_pos=tg_.surfaceIntercept();
  Uvar height = target_pos.magnitude()-tg_.radius();  

  // Output for DEBUG level 1 outputs everything that goes into doppler computation
  if(dbg.level==1){
      dbg.file.precision(20);
      dbg.file << "  lambda_chirp=" << lambda_chirp.getInUnits("km")
      << ", look_x = " << target_look[DirectionVector::X]
      << ", look_y = " << target_look[DirectionVector::Y]
      << ", look_z = " << target_look[DirectionVector::Z]
      << ", vel_x = " << v[FloatVector::X].getInUnits("km/s") << " km/s,"
      << "vel_y = " << v[FloatVector::Y].getInUnits("km/s") << " km/s,"
      << "vel_z = " << v[FloatVector::Z].getInUnits("km/s") << " km/s,"
      << ", fdop=" << fdop_in_Hz << endl;
  } // end debug level = 1


  if(target_info_in_diag_file){
    diag_file.set("fdop",fdop_in_Hz,target_idx+1);
    diag_file.set("range",range.getInUnits("km"),target_idx+1);
    diag_file.set("scale",scale,target_idx+1);
    diag_file.set("lat",target_lat(target_idx).getInUnits("deg"),target_idx+1);
    diag_file.set("lon",target_lon(target_idx).getInUnits("deg"),target_idx+1);
    diag_file.set("height",height.getInUnits("m"),target_idx+1);
  }
  if(ang<max_target_angle) return(true);
  else return(false);
}

// sets Receive Window delay and chirp start frequency so
// beam boresight is as close to center of time and frequency gate as 
// possible 
void PointTargetSim::idealDopplerRangeTrack()
{
  // manual profiling 
  double cs=get_time();


  // set time to burst_start_time compute delay
  setBeamTime(l1b->beam_number,l1b->t);
  Uvar delay=tg_.range()*2/speed_light;


  //let the center of the transmit window return as closes as possible
  // to center of recieve gate;
  Uvar transmit_center=l1b->pul*l1b->pri/2;
  Uvar receive_center=transmit_center+delay;
  Uvar receive_window=l1b->ctrx*l1b->pri;
  Uvar ideal_rwd=receive_center - receive_window/2;
  l1b->rwd=floor((ideal_rwd/l1b->pri).getInUnits("") + 0.5);
  l1b->rwd*=l1b->pri;

  //let the center frequency be as close as possible to the center of 
  // receive freq window.
  // this ignores the small effect of varying csf on doppler
  Uvar center_freq=filter_start_freq + l1b->rc_bw/2;
  Uvar start_freq=center_freq-chirp_rate*l1b->chirp_length/2;
  Uvar ideal_csf=start_freq-tg_.doppler(lambda_)-updown_shift;
  Uvar csf_quanta=Uvar(457.764,"Hz");
  l1b->fast_csf=floor((ideal_csf/csf_quanta).getInUnits("") + 0.5);
  l1b->fast_csf*=csf_quanta;


  // Print out ideal tracking results
#ifdef DEBUG
  cout << "Tracking:lat:" << tg_.lat().getInUnits("deg") << " lon:" 
       << tg_.lon().getInUnits("deg") << " fast_csf:" << l1b->fast_csf 
       << " rwd:" << l1b->rwd << endl;
  cout << "Tracking:center_freq=" << center_freq <<" start_freq=" << 
    start_freq << endl;
  cout << "Tracking ideal_csf=" << ideal_csf << " updown_shift=" <<
    updown_shift << " doppler=" << tg_.doppler(lambda_) << endl;
#endif

  //sanity checks
  if(l1b->rwd+l1b->ctrx*l1b->pri > l1b->bpd){
    ErrorMessage e("PointTargetSim Receive window extends past end of burst");
    e.throwMe();
  }
  if(l1b->pul*l1b->pri > l1b->rwd){
    ErrorMessage e("PointTargetSim Transmit window extends into RX window");
    e.throwMe();
  }

  // manual profiling
  double ce=get_time();
  pts_ideal_track_time+=ce-cs;

}

// This routine is called for each point target
// it computes the returned signal sample by sample
// then it adds the corresponding contribution to the radar_data_buffer
float PointTargetSim::addEcho(const int target_idx)
{
    DebugInfo dbg("PointTargetSim::addEcho");
    if(dbg.level)
      dbg.file << "PointTargetSim::addEcho burst number " << burst_no << " Target "  << target_idx << endl;
    // initialize return value
    float average_ampl=0.0;
    int num_valid_points=0;

    // Special debugging routine coordinated among L1I, SARProcParams and
    // PointTargetSim objects
    // Each bit in dbgsp.level turns on/off (0/1)  
    // a particular calibration parameter.
    // OFF means set parameter to unity
    // 
    // Parameter Flags bits from LSB to MSB
    // (from L1I::calibrate) bits 0-5
    // Xcal, g2, area, range^4, cag, norm_value,  ...  
    //
    // (from SARProcParams::Xcal . i.e., factors in Xcal) 
    // bits 6-9
    // Xconstant_==Pt*lambda*lambda*X_FACTOR_CONSTANT_CORRECTION/64/pi/pi/pi,
    // Xbeam_, Xtimevar_, Xmode_==sar_cal_coeff
    //
    // (from PointTargetSim::addEcho)
    // bits 10-12
    // scale, area_correction, sigma0_value==ValueFromBackscatterModel
    // 
    // (from PointTargetSim::convertToByte)
    // bit 13
    // gain == attenuator_gain
    // 
    // (from PointTargetSim::getDelayDopplerAndScale)
    // bits 14-18
    // target_ampl==AMPLFROMFILE*footprint_fraction, 
    // Xconstant==transmit_power*lambda_*lambda_/(64*pi*pi*pi),
    // g2,range^4,2*sar_cal_coeff (OFF sets sar_cal_coeff=0.5)
    DebugInfo dbgsp("Special_Calibration_Debug");


    // manual profiling
    double cs=get_time();

    // compute area correction associated with point target
    double area_correction=1.0; 


    if(simulate_area_ && area_is_relative_){
 
      // This is only computed for target_type=4 (elev,azim) point target grid
      // compute area correction for a point target at the center of the grid
      // with the prescribed resolution in azim and elevation
      Uvar az1=grid_azim_center - grid_azim_res/2.0;
      Uvar el1=grid_elev_center - grid_elev_res/2.0;
      Uvar az2=grid_azim_center + grid_azim_res/2.0;
      Uvar el2=grid_elev_center + grid_elev_res/2.0;  
 
      // set up "corner points" for rectangle in beam frame
      PositionVector p[4];
      int b=l1b->beam_number;
      Time t=burst_start_time;
      Uvar target_radius=(default_target_radii.magnitude())/sqrt(3.0);
  
      setBeamTimeAzimuthElevation(b,t,az1,el1);
      p[0]=tg_.surfaceIntercept();
      setBeamTimeAzimuthElevation(b,t,az1,el2);
      p[1]=tg_.surfaceIntercept();
      setBeamTimeAzimuthElevation(b,t,az2,el1);
      p[2]=tg_.surfaceIntercept();
      setBeamTimeAzimuthElevation(b,t,az2,el2);
      p[3]=tg_.surfaceIntercept();

      Uvar dummy1,dummy2,dummy3;
      Uvar area=getSphericalTriangleArea(p[0],p[1],p[2],target_radius,
				   dummy1,dummy2,dummy3);
      area+=getSphericalTriangleArea(p[1],p[2],p[3],target_radius,
				   dummy1,dummy2,dummy3);

      area_correction=get_in_base_units(area);
    }

    // compute backscatter model value if desired
    double sigma0_value=1.0;
    if(simulate_backscatter_){
      setBeamTimeLatLonHeight(l1b->beam_number,burst_start_time,
			      target_lat(target_idx),target_lon(target_idx),
			      target_height(target_idx));
      Uvar inc=tg_.incidenceAngle();
      sigma0_value=muhleman_backscatter(muhleman_k1_,muhleman_k2_,inc);
    }

    // compute phase shift of point target return
    double phi0=0.0;
    if(use_random_phase_shift_){
      phi0=Uniform(pi,0).GetNumber();
    }
    // compute start of pass band
    double pass_band_start=filter_start_freq_in_Hz; // assumes that anti-aliasing filter
                            // is ideal and passes frequencies from 
                            // -ADC/2 to -ADC/2 + the receiver bandwidth
                            // of the ADC.

    // compute end of pass band
    double pass_band_end=pass_band_start+rc_bw_in_Hz;


    // compute sampling interval 
    double sample_interval=1/adc_in_Hz;

    // pulse_width
    double  pulse_width=chirp_length_in_s; 
    
  
    int pulsenum=0;

    
    // all times relative to transmit window start (except tp)
    double te_start=rwd_in_s;
    double te=0;
    double tp;
    double f0=	fast_csf_in_Hz + 
      updown_shift_in_Hz;
    double next_delay=0,next_fdop=0,next_scale=0;
    double ddelay=0,dfdop=0,dscale=0;
    double cr= chirp_rate_in_Hz_per_s;
    double f_start,f_inst;
    delay_in_s=0;
    fdop_in_Hz=0;
    scale=0;
    float value=0;

    int first_active_sample=1;
    for(unsigned int c=0;c<l1b->Nradar_data;c++){

      // compute time te at middle of sample
      te=te_start+c*sample_interval; 

      // compute round_trip_time

      if(c%sample_skip==0 || c==l1b->Nradar_data-1){
      
	// compute derivatives if necessary
        if(sample_skip!=1 && c!=l1b->Nradar_data-1){
          double of,od,os;
          // Make sure initial values are obtained
          if(c==0){
	    getDelayDopplerAndScale(te,target_idx);
	    next_fdop=fdop_in_Hz;
	    next_delay=delay_in_s;
	    next_scale=scale;
	  }
          // Save previously calculated values for current time instant
          of=next_fdop;
          od=next_delay;
          os=next_scale;

          // Compute values for next sample skip value
	  double big_interval=sample_interval*sample_skip;
	  getDelayDopplerAndScale(te+big_interval,
				  target_idx);
	  next_fdop=fdop_in_Hz;
          next_delay=delay_in_s;
	  next_scale=scale;
          
          // Swap current value back in.
          fdop_in_Hz=of;
          // if it is going beyond the horizon pre-empt it!
          if(next_delay==0) delay_in_s=0;
          else delay_in_s=od;
          scale=os;
	  
          // Compute derivatives
          dfdop=(next_fdop-fdop_in_Hz)/big_interval;
          ddelay=(next_delay-delay_in_s)/big_interval;
	  dscale=(next_scale-scale)/big_interval;
	} // end case in which derivatives need to be calculated

	// last sample OR sample_skip=1 special case
	else  getDelayDopplerAndScale(te,target_idx);

	if(delay_in_s > max_rtt_in_s) max_rtt_in_s=delay_in_s;
	if(delay_in_s < min_rtt_in_s) min_rtt_in_s=delay_in_s;
	if(fdop_in_Hz > max_dop_in_Hz) max_dop_in_Hz=fdop_in_Hz;
	if(fdop_in_Hz < min_dop_in_Hz) min_dop_in_Hz=fdop_in_Hz;

      } // end full calculation performed this interval case
      else{
        // approximate using derivatives
	fdop_in_Hz+=dfdop*sample_interval;
        // if beyond horizon at current time assume it stays so until
        // next delay calculation
	if(delay_in_s!=0) delay_in_s+=ddelay*sample_interval;
        scale+=dscale*sample_interval;
      }

      // check for beyond horizon case (delay=0), 
      // or zero gain case scale=0.0
      // in either case this point target
      // does not contribute to this sample
      if (scale==0.0) continue;
      if (delay_in_s==0.0) continue;
      


      // compute transmit_time
      double tt=te-delay_in_s;


      // compute time from start and number of most recent pulse;
      pulsenum=int(floor(tt/pri_in_s));
      tp=tt-pulsenum*pri_in_s;

      // compute time margins and update global values if necessary

      if(c==0){
	float margin=-(tp+pulsenum*pri_in_s);
	if(margin < start_time_margin_in_s)  start_time_margin_in_s=margin;
      }

      if(c==l1b->Nradar_data-1){
	float margin=(tp+pulsenum*pri_in_s)-pri_in_s*l1b->pul;
	if(margin < end_time_margin_in_s)  end_time_margin_in_s=margin;
      }

      // if tp is inside pulse then compute echo;
      if(tp>0 && tp<pulse_width && pulsenum < (int)l1b->pul && pulsenum>=0){
	// using chirp rate and chirp starting frequency compute value of
	// waveform and instantaneous frequency at transmit time tt
	//ignoring chirp quantization for now




	// compute instantaneous freq
        f_start=f0+fdop_in_Hz;
	f_inst=f_start+cr*tp;


	// compute frequency margins and update global values if necessary

	float margin=f_inst-pass_band_start;
	if(margin < start_freq_margin_in_Hz)  start_freq_margin_in_Hz=margin;

	margin=pass_band_end-f_inst;
	if(margin < end_freq_margin_in_Hz)  end_freq_margin_in_Hz=margin;


	// No signal if instanteous freq is outside of passband
        // THIS IS AN APPROXIMATION
        if(f_inst >= pass_band_start && f_inst <= pass_band_end){


	  float ampl=scale*sqrt(area_correction*sigma0_value);

          // Special_Calibration_Debug Case
          if(dbgsp.level){
	    float s2=scale;
            float a2=area_correction;
            float sig2=sigma0_value;
	    if(0x0400 & dbgsp.level) s2=1.0;
	    if(0x0800 & dbgsp.level) a2=1.0;
            if(0x1000 & dbgsp.level) sig2=1.0;
	    ampl=s2*sqrt(a2*sig2);
    
            // writes output to debug file for first active sample only
	    // this is only written when Special_Calibration_Debug has a
            // a nonzero debug level. To get this output w/o disabling
            // any calibration parameters use a nonzero flag value with the
            // 19 least significant bits set to zero i.e., 1048576
	    if(first_active_sample && target_special(target_idx)){
	      first_active_sample=0;
	      dbgsp.file << endl << "PointTargetSim::addEcho "
			 << "Burst number " << burst_no 
			 << " Beam number " << l1b->beam_number << endl
			 << " CalParams:  Echo amplitude ="
		         << ampl << " Area correction=" 
			 << area_correction 
			 << " km^2  Sigmao value=" << sigma0_value 
			 << "  scale=" << scale << endl
                         << " Location:"
			 << "  range=" 
			 << get_in_base_units((delay*speed_light)/2)
			 << " km   fdop=" << get_in_base_units(fdop)/1000 
			 << " kHz " 
			 << endl; 
	    }
	  }


          average_ampl+=ampl;
          num_valid_points++;

          
	  value=get_point_target_echo(tp,te,updown_shift_in_Hz,
				      carrier_freq_in_Hz,
				      delay_in_s,
				      fast_csf_in_Hz,cr,
				      ampl,
				      phi0);
	  l1b->radar_data(c)+=value;

	}


        // Debug level 1 output calibartion info for each point target and each burst
	if(dbg.level==1 && c==l1b->Nradar_data/2 ){
	  Uvar Xconstant=
	    2*sar_cal_coeff*transmit_power*lambda_chirp*lambda_chirp
	    /64/pi/pi/pi;


// Debugging instruments
	  dbg.file << "PointTargetSim::addEcho "
		   << "  Area correction=" << area_correction 
	           << "  Sigmao value=" << sigma0_value 
	           << endl 
                   << "  scale=" << scale
		   << "  Lamda_chirp=" << lambda_chirp 
		   << "  range=" << (delay*speed_light)/2
		   << "  fdop=" << fdop 
	           << endl 
	           << " target_azim=" << target_azim(target_idx)
	           << " target_elev=" << target_elev(target_idx)
	           << " Target Amplitude=" << target_ampl(target_idx)
	           << " Target Index=" << target_idx 
		   << endl 
		   << " Area x sigma0 = " << 
	                target_ampl(target_idx)*area_correction*sigma0_value
		   << " 2*sar_cal*Pt*lamda^2/64pi^3 =" << Xconstant
	           << " G^2/R^4=" << scale/target_ampl(target_idx)/Xconstant
                   << " Beam number =" << l1b-> beam_number 
		   << endl << endl;
	}


        // Output for DEBUG level 2 keeps track of delay and doppler for each sample and each point target
	if(dbg.level==2){
	  if(c%sample_skip==0){
	    dbg.file.precision(20);
	    dbg.file << "sample " << c << "  te=" << te << ", delay=" << delay_in_s 
	          << ",tp=" << tp << ", updown=" <<  updown_shift_in_Hz
                  << ", carrier_freq=" <<  carrier_freq_in_Hz 
	          << ", fast_csf=" << fast_csf_in_Hz
	          << ", cr=" << cr << ", scale=" << scale 
		  << ", fdop=" << fdop_in_Hz << ", f_inst=" << f_inst
	          << ", value=" << value << endl;
          }
	} // end debug level = 2

      
      }



                                       
    }

    average_ampl/=num_valid_points;
    // manual profiling
    double ce=get_time();
    pts_addecho_time+=ce-cs;
    return(average_ampl);
}

void PointTargetSim::addNoise()
{
  noise_var=boltzmann_constant*Tsys*l1b->rc_bw*sar_cal_coeff;
  Gaussian g(noise_var.getInUnits(""),0.0);
  for(unsigned int c=0;c<l1b->Nradar_data;c++){
    float noise_val=g.GetNumber();
    l1b->radar_data(c)+=noise_val;
  }
}

void PointTargetSim::convertToByte(){

  // Special Calibration Debugging see comments in addEcho()
  DebugInfo dbgsp("Special_Calibration_Debug");
  float gain=0;
  if(ideal_gain_setting_){

    // In this case an idealized gain value is used
    // so that the maximum/minimum values barely fit within 
    // a signed char.
    // Attenuator settings in the l1b are set accordingly.

    float m=fabs(l1b->radar_data(0));
    for(unsigned int c=1;c<l1b->Nradar_data;c++){
      float x=fabs(l1b->radar_data(c));
      if(x>m) m=x;
    }
    gain=127/m;
    double lossdB=-10*log10(gain);
    int nominal_loss_dB= int(ceil(lossdB));
    
    if(nominal_loss_dB < 0 ){
      string s = 
	"PointTargetSim::convertToByte Ideal attenuator loss < 0 dB\n ";
	s=s+"Increase point Target amplitude by " + toStr(-nominal_loss_dB);
        s=s+" dB\n or disable simulate_quantization_ or ideal_gain_setteing_"; 
      ErrorMessage e(s);
      e.throwMe();
    }

    if(nominal_loss_dB > 72 ){
      string s = 
	"PointTargetSim::convertToByte Ideal attenuator loss < 0 dB\n ";
	s=s+"Decrease point Target amplitude by " + toStr(72-nominal_loss_dB);
        s=s+" dB\n or disable simulate_quantization_ or ideal_gain_setting_."; 
      ErrorMessage e(s);
      e.throwMe();
    }

    
    float gain2=1;

    // make sure the calibrated gain is less than the maximal gain
    // may require a couple iterations due to difference between nominal
    // gain and calibrated gain from CAGTable.
    while(gain2>gain && nominal_loss_dB <= 37){
      l1b->at1_db=(unsigned)nominal_loss_dB*2; // gain is power to power!!!
      l1b->at3_db=(unsigned)nominal_loss_dB*2;
      l1b->at4_db=(unsigned)nominal_loss_dB*2;
      l1b->at1 = cag_table_.gainLinearScale((int)l1b->at1_db);
      l1b->at3 = cag_table_.gainLinearScale((int)l1b->at3_db);
      l1b->at4 = cag_table_.gainLinearScale((int)l1b->at4_db);
      // atx_each is not assigned needs table used by RMSS

      gain2=l1b->computeCalibratedAttenuatorGain();
      nominal_loss_dB++;
    }
    gain=gain2; // Attenuator gain is power to power like everything
                // else
  }
  else {
    // gain is determined from attenuator settings in L1B and calibrated
    // attenuator gain (CAG) table
    gain=l1b->computeCalibratedAttenuatorGain();
  }

  if(dbgsp.level){
    if (0x2000 & dbgsp.level) gain =1.0;

    // To get this output w/o disabling
    // any calibration parameters use a nonzero flag value with the
    // 19 least significant bits set to zero i.e., 1048576
    dbgsp.file << endl << "PointTargetSim::convertToByte gain =" << gain 
	       << endl;
  }
  
  for(unsigned int c=0;c<l1b->Nradar_data;c++){
    // gain is power to power liek everything else
    l1b->radar_data(c)*=sqrt(gain); 
    if(simulate_quantization_){
      l1b->radar_data(c)=floor(l1b->radar_data(c))+0.5;
      if(l1b->radar_data(c) > 127.5) l1b->radar_data(c)=127.5;
      if(l1b->radar_data(c) < -127.5) l1b->radar_data(c)=-127.5;
    }
  }
}

void PointTargetSim::simulateBaq(){
  Uvar tro_pri = l1b->tro/l1b->pri;
  int tro_int=int(floor(tro_pri.getInUnits("")+0.5));
  baq.setParams(l1b->pri,l1b->pul,l1b->baq_mode,l1b->adc,tro_int);
  Ivec Thresh("Thresh");
  Charvec tmp("tmp",l1b->Nradar_data);
  Fvec tmp_float("tmp_float",l1b->Nradar_data);
  for(unsigned int c=0;c<l1b->Nradar_data;c++) 
    tmp(c)=(char)l1b->radar_data(c);
  baq.compuThreshold(tmp, Thresh);

  Charvec words("words_name");
  baq.Encode_Nbit(tmp, Thresh, words);

  baq.Decode_Nbit(words, Thresh, tmp_float);
  for(unsigned int c=0;c<l1b->Nradar_data;c++) 
    l1b->radar_data(c)=tmp_float(c);  
}

void PointTargetSim::getDelayDopplerAndScale(double t_receive, 
					const int target_idx){

  static int last_target_idx=-1;
  if(!configured_)
    {
      ErrorMessage e("CassiniSim not configured");
      e.throwMe();
    }
  
  // Special Calibration Debugging see comments in addEcho()
  DebugInfo dbgsp("Special_Calibration_Debug");

  setBeamTimeLatLonHeight(l1b->beam_number,burst_start_time+Uvar(t_receive,"s"),
			  target_lat(target_idx),target_lon(target_idx), 
			  target_height(target_idx));

  // Beyond horizon special case
  if(!tg_.foundSurfaceIntercept()){
    delay=Uvar(0,"s");
    delay_in_s=0;
    return;
  }

  receive_range = tg_.range();
  delay = 2*receive_range/speed_light;
  delay_in_s=get_in_base_units(delay);



  // compute antenna_gain on receive
  double antenna_gain_receive;
  if(compute_echo_amplitude_){
    look_direction=tg_.lookDirection();
    look_direction=
      look_direction.representIn(beam_[cur_beam_].getFrame(dummy_time));
    look_direction.getAzimuthElevation(azimuth,elevation);
    antenna_gain_receive=beam_[cur_beam_].bilinear(azimuth,elevation);
  }
  // Outside of antenna pattern special case
  if(antenna_gain_receive==0){
    scale=0;
    return;
  }

  // Doppler is approximate
  // computed as average of receive time and transmit time doppler 
  // using carrier frequency only (ignores chirp)
  // This is sufficient as it is only used to estimate anti-aliasing filter
  // effect which is also an approximation.
  fdop=tg_.doppler(lambda_chirp);
  

  // iterate once to get reasonably accurate time
  delay_previous_iteration = delay;
  delay_difference=delay;
  while(delay_difference>=Uvar(10e-14,"s")){
    double tt= t_receive-delay_in_s;
    setBeamTimeLatLonHeight(l1b->beam_number,burst_start_time+Uvar(tt,"s"),
			    target_lat(target_idx),target_lon(target_idx),
			    target_height(target_idx));

    // Beyond horizon special case
    if(!tg_.foundSurfaceIntercept()){
      delay=Uvar(0,"s");
      delay_in_s=0;
      return;
    }
    transmit_range = tg_.range();
    delay=(receive_range+transmit_range)/speed_light;
    delay_difference=delay-delay_previous_iteration;

    delay_previous_iteration=delay;
    delay_in_s=get_in_base_units(delay);
  }

  fdop+=tg_.doppler(lambda_chirp);
  fdop/=2;
  fdop_in_Hz=get_in_base_units(fdop);
  
  double antenna_gain_transmit;
  // compute antenna_gain on transmit
  if(compute_echo_amplitude_){
    look_direction=tg_.lookDirection();
    look_direction=
      look_direction.representIn(beam_[cur_beam_].getFrame(dummy_time));
    look_direction.getAzimuthElevation(azimuth,elevation);
    antenna_gain_transmit=beam_[cur_beam_].bilinear(azimuth,elevation);

    
    dXdA=transmit_power*lambda_*lambda_*
      antenna_gain_receive*antenna_gain_transmit/
      (64*pi*pi*pi*transmit_range*transmit_range*
       receive_range*receive_range);

    power_receive=dXdA*target_ampl(target_idx);

    // constant 2 accounts for fact that sar_cal_coeff is the ratio
    // of RSS data number to power recieve
    // but scale is the amplitude of a sinusoid 
    scale_uvar=sqrt(2*sar_cal_coeff*power_receive);

    // Special_Calibration_Debug case (see addEcho for more details)
    if(dbgsp.level){
      float ta=get_in_base_units(target_ampl(target_idx));
      float Xconstant=get_in_base_units(transmit_power*lambda_*lambda_/
					(64*pi*pi*pi));
      float g2=get_in_base_units(antenna_gain_receive*antenna_gain_transmit);
      float r4=get_in_base_units(transmit_range*transmit_range*
				 receive_range*receive_range);
      float cal_coeff=2*get_in_base_units(sar_cal_coeff);

      if(0x4000 & dbgsp.level) ta=1.0;
      if(0x8000 & dbgsp.level) Xconstant=1.0;
      if(0x10000 & dbgsp.level){
	float maxgain =beam_[cur_beam_].getMaxGain();
	maxgain*=maxgain;
        if(g2<maxgain/2.0) g2=0.0;
	else g2=1.0;
      }
      if(0x20000 & dbgsp.level) r4=1.0;
      if(0x40000 & dbgsp.level) cal_coeff=1.0;

      float pr=ta*Xconstant*g2/r4;
      scale_uvar=sqrt(cal_coeff*pr);

 
      // writes output to debug file for first scale computation in
      // special point target
      // this is only written when Special_Calibration_Debug has a
      // a nonzero debug level. To get this output w/o disabling
      // any calibration parameters use a nonzero flag value with the
      // 19 least significant bits set to zero i.e., 1048576
      if(last_target_idx!=target_idx && target_special(target_idx)){
	last_target_idx=target_idx;
	dbgsp.file << endl << "PointTargetSim::getDelayDopplerAndScale "
		   << "target_ampl =" << ta
		   << ", Xconstant =" << Xconstant << endl
		   << " g2 =" << g2
		   << ", range^4 =" << r4
		   << " km^4,  2*sar_cal_coeff=" << cal_coeff << endl;
      }

    } // end of Special_Calibration_Debug case
    scale=get_in_base_units(scale_uvar);



  }
 
  else{
    scale=1.0;
  }

}

void
PointTargetSim::computeDoubleParams()
{
 chirp_rate_in_Hz_per_s=chirp_rate.getInUnits("Hz/s");
 rc_bw_in_Hz=l1b->rc_bw.getInUnits("Hz");
 filter_start_freq_in_Hz=filter_start_freq.getInUnits("Hz");
 adc_in_Hz=l1b->adc.getInUnits("Hz");
 chirp_length_in_s=l1b->chirp_length.getInUnits("s");
 carrier_freq_in_Hz=(speed_light/lambda_).getInUnits("Hz");
 rwd_in_s=l1b->rwd.getInUnits("s");
 fast_csf_in_Hz=l1b->fast_csf.getInUnits("Hz");
 updown_shift_in_Hz=updown_shift.getInUnits("Hz");
 pri_in_s=l1b->pri.getInUnits("s");
}

void PointTargetSim::altBAEToLatLon(const int target_idx){

  // Get target altitude and tolerance
  Uvar af=target_alt(target_idx);
  Uvar alt_tolerance=Uvar(0.001,"km");

  // Determine time and altitude ranges for search
  Time t1,t2;
  t1=flyby.lowestAltitudeTime();

  if(target_inbound(target_idx)){
    if(t1 > end_time) t1=end_time;
    t2=start_time;
  }
  else {
    if(t1 < start_time) t1=start_time;
    t2=end_time;
  }

  setTime(t2);
  Uvar a2=altitude();
  setTime(t1);
  Uvar a1=altitude();

  // Check for out of range condition
  if ( a1 > af && a2 > af || a1 < af && a2 < af ){
    OSTRINGSTREAM os;
    os << "Error:PointTargetSim target altitude (" << af
       << ") is outside altitude range ("<< a1 <<"," << a2 <<")";
    ErrorMessage e(toStr(os));
    e.throwMe();
  }


  // Intial altitude, Time guess and altitude error
  Uvar ai=a1;
  Time ti=t1;
  Uvar alt_error=fabs(ai-af);

  //-----------------------------------------------------
  // Binary search to find time corresponding to altitude
  //-----------------------------------------------------
  while(alt_error>alt_tolerance){
    // compute altitude && time guess
    ti=t1+((af-a1)/(a2-a1))*(t2-t1);
    setTime(ti);
    ai=altitude();

    // compute error
    alt_error=fabs(ai-af);

    // compute new bounds
    if((af>a1 && af<ai)||(af>ai && af<a1)){
      a2=ai;
      t2=ti;
    } 
    else{
      a1=ai;
      t1=ti;
    }
  }

  // Set target altitude to value from search so that simulation is exact
  target_alt(target_idx)=ai;

  // Convert beam, azim, and elevation to lat lon coordinates
  timeBAEToLatLon(ti,target_idx);
  
}

void PointTargetSim::timeBAEToLatLon(const Time& t, const int target_idx)
{
  setTime(t);
  Frame beamframe=beam_[target_beam(target_idx)-1].getNominalFrame();
  DirectionVector look("look",beamframe,t,0,0,0);
  look.setAzimuthElevation(target_azim(target_idx),target_elev(target_idx));
  tg_.setLookDirection(look);
  target_lat(target_idx)=tg_.lat();
  target_lon(target_idx)=tg_.lon();
}

void PointTargetSim::getRadarParams(Config& cfg){
  data_take_id=cfg.getInt("data_take_number");
  l1b->sync=SYNC_VAL;
  l1b->rc_bw=cfg["RECEIVER_BANDWIDTH"];
  l1b->csr=0;
  l1b->r_mode=cfg.getInt("PTSIM_RADAR_MODE"); // set radar mode 0 = SCAT, 1= ALT, 2= SARL, 3= SARH


  l1b->csd=Uvar(1e-09,"s"); // refects continuous chirp assumption
  l1b->adc=cfg["ADC_SAMPLE_RATE"];

  if(filter_centered){
    filter_start_freq=-l1b->adc/2+(l1b->adc/2-l1b->rc_bw)/2;
  }
  else{
    filter_start_freq=-l1b->adc/2.0+filter_start_freq_offset;
  }

  l1b->baq_mode=cfg.getInt("BAQ_MODE");


  l1b->pri=1/(cfg["PRF"]);
  unsigned int tro_in_pri=cfg.getInt("RX_WINDOW_EXTRA_PRI");
  l1b->tro=tro_in_pri*l1b->pri;
  l1b->chirp_length=l1b->pri*cfg["INTRA_BURST_DUTY_CYCLE"];

  double max_duty_cycle = cfg.getDouble("TOTAL_DUTY_CYCLE");
  Uvar ratio=l1b->bpd/l1b->chirp_length;
  l1b->pul=(long unsigned int)floor(max_duty_cycle*ratio.getInUnits(""));
  l1b->ctrx=l1b->pul+tro_in_pri;

  // configure whether or not (1/0) to compute an ideal
  // attenuator  gain setting for each burst
  // if not gain is comnfigured
  if(!ideal_gain_setting_) {
    l1b->at1_db=cfg.getInt("ATTENUATOR_LOSS_1");
    l1b->at3_db=cfg.getInt("ATTENUATOR_LOSS_3");
    l1b->at4_db=cfg.getInt("ATTENUATOR_LOSS_4");
    l1b->at1 = cag_table_.gainLinearScale((int)l1b->at1_db);
    l1b->at3 = cag_table_.gainLinearScale((int)l1b->at3_db);
    l1b->at4 = cag_table_.gainLinearScale((int)l1b->at4_db);
    // atx_each not assigned needs table used by RMSS to set masks
  }

  // compute number of samples taken
  Uvar receive_window_length=l1b->ctrx*l1b->pri;
  
  Uvar num_samples=receive_window_length*l1b->adc;
  double num_samps_int=floor(num_samples.getInUnits("")+0.5);
  if(num_samples.getInUnits("")!=num_samps_int){
    num_samples=num_samps_int;
    l1b->tro/=l1b->pri;
    l1b->chirp_length/=l1b->pri;
    l1b->pri=num_samples/l1b->adc/l1b->ctrx;
    l1b->tro*=l1b->pri;
    l1b->chirp_length*=l1b->pri;
    if(DebugInfo::allWarnings){
      cerr << "Warning PointTargetSim::getRadarParams Non integer number of samples." ;
      cerr << endl << "Changing PRF to: "<< (1/(l1b->pri)).getInUnits("KHz") << " KHz" << endl;
    }
  }

  l1b->Nradar_data=(unsigned int)(num_samples.getInUnits(""));

  if(l1b->Nradar_data > 32768){
    ErrorMessage e("CONFIG mode produced too many samples per burst, try reducing TOTAL_DUTY_CYCLE");
    e.throwMe();
  }

  l1b->csq=(unsigned int)floor((l1b->chirp_length/l1b->csd).getInUnits(""));
  l1b->csq-=1;

  if(!ideal_tracking_on_){
    l1b->fast_csf=cfg["CHIRP_START_FREQ"];
    l1b->rwd=cfg["RX_WINDOW_DELAY_PRI"]*l1b->pri;
  }

  // will not use this chirp rate exactly
  chirp_rate=cfg["CHIRP_RATE"];

  l1b->slow_cfs=chirp_rate*(l1b->csq+1)*l1b->csd/l1b->csq;
  if(DebugInfo::allWarnings)
    cerr << "Warning not using exact config file chirp rate =" << chirp_rate;
  setChirpRate();
  if(DebugInfo::allWarnings)
    cerr <<" but instead chirp_rate = " << chirp_rate << endl;

  //-----------------------
  // Beam Number Order
  //-----------------------

  //  Load  first beam order value 
  beam_order[0] = (unsigned int) cfg.getInt("BEAM_ORDER_1");

  unsigned int i = 2;
  while (i<=5)
    {
      OSTRINGSTREAM os;
      os << "BEAM_ORDER_" << i;
      if (cfg.keywordExists(toStr(os)))
        {
	  beam_order[i-1]=(unsigned)cfg.getInt(toStr(os));
        }
      else
        {
	  break;  // only load a sequential set of names
        }
      // Next beam_order is set to zero
      // if a value is available for it in the config file
      // it will be overwritten
      beam_order[i]=0;
      ++i;
    }
}

void PointTargetSim::getRadarParamsFromIeb(){
  DebugInfo dbg("PointTargetSim::getRadarParamsFromIeb");
  if(dbg.level){
    Time et =flyby.epochTime();
    dbg.file << "PointTargetSim::getRadarParamsFromIeb: start_time " 
	     << start_time -et << " end_time " << end_time -et
             << " t_ " << t_ -et << " delta_trigger " << ieb_delta_trigger_time
	     << " lastIEBtime " << ieb_file->getEndTime() - et << endl;
  }
  ieb = ieb_file->getIeb(t_-ieb_delta_trigger_time);
  l1b->csr=ieb.getCsr();
  l1b->r_mode=ieb.getMod();
  data_take_id=ieb.getDtn();
  l1b->sync=SYNC_VAL;
  if(!use_config_bpd_){
    l1b->bpd=ieb.getBpd();
  }
  l1b->rc_bw=ieb.getRcv();
  l1b->adc=ieb.getAdc();

  // configure starting frequency of anti-alias filter
  // this is read from the config file until the precise formulation
  // is determined we assume this is -adc/2

  if(filter_centered){
    filter_start_freq=-l1b->adc/2+(l1b->adc/2-l1b->rc_bw)/2;
  }
  else{
    filter_start_freq=-l1b->adc/2.0+filter_start_freq_offset;
  }



  l1b->csd= ieb.getCsd(); // refects continuous chirp assumption

  l1b->baq_mode= ieb.getBaq();


  l1b->pri = ieb.getPri();
  int tro_in_pri=ieb.getTroInPriUnits();
  l1b->tro=tro_in_pri*l1b->pri;


  l1b->chirp_length=ieb.getTaup();


  l1b->pul=ieb.getPul();
  l1b->ctrx=l1b->pul+tro_in_pri;

  if(!ideal_gain_setting_) 
    { 
      l1b->at1_db=ieb.getAt1();
      l1b->at3_db=ieb.getAt3();
      l1b->at4_db=ieb.getAt4();
      l1b->at1 = cag_table_.gainLinearScale((int)l1b->at1_db);
      l1b->at3 = cag_table_.gainLinearScale((int)l1b->at3_db);
      l1b->at4 = cag_table_.gainLinearScale((int)l1b->at4_db);
      ieb.getAttenuationBitMask(l1b->at1_each,l1b->at3_each,l1b->at4_each);
    }

  // compute number of samples taken
  l1b->ctrx=tro_in_pri+l1b->pul;
  Uvar receive_window_length=l1b->ctrx*l1b->pri;
  Uvar num_samples=receive_window_length*l1b->adc;
  int num_samps_int= int(num_samples.getInUnits("")+0.5);

  l1b->Nradar_data=num_samps_int;


  l1b->csq= ieb.getCsq();
  l1b->slow_cfs = ieb.getCfs();

  if(!ideal_tracking_on_){
    l1b->fast_csf=ieb.getChirpStartFrequency();
    l1b->rwd=ieb.getRwd();
  }

  setChirpRate();


  //-----------------------
  // Beam Number Order
  //-----------------------

  //  Load  beam ordering 
  int bem=ieb.getBem();

  int num_beams=0;
  for(int c=0;c<5;c++){
    if(bem & 1){
      beam_order[num_beams]=c+1;
      num_beams++;
    }
    bem = bem >> 1;
  }
  
  beam_order[num_beams]=0; // Add trailing zero to beam order

  if(num_beams<1 && !l1b->isCal()){
    ErrorMessage e("Error:PointTargetSim: No beams enabled for normal ops");
    e.throwMe();
  }
}

void PointTargetSim::setChirpRate(){
  chirp_rate=l1b->csq*l1b->slow_cfs*l1b->adc/
    floor(get_in_base_units(l1b->chirp_length*l1b->adc+0.5));
}
void PointTargetSim::updateRadarParams(){
  double cs=get_time();

  if(radar_param_source != L1BFILE){
    burst_no++;
    l1b->createRecord(t_);
    l1b->record_id=1000000*data_take_id+burst_no;
    // if the L1B record creation fails due to bad
    // geometry  throw an error
    if(!l1b->goodGeometry()){
      char msg[60];
      sprintf(msg,"Bad L1A-L1B convert: quality_flag=%d",(int)l1b->quality_flag);
      ErrorMessage e(msg);
      e.throwMe();
    }
    if(beam_order[beam_order_idx]==0) beam_order_idx=0;
    l1b->beam_number=beam_order[beam_order_idx];
    beam_order_idx++;
    l1b->sab_counter=burst_no;
  }
  else{
    burst_no=l1b->record_id%1000000;
  }
  // update lambda_chirped wavelength for center of chirp
  Uvar chirp_freq=speed_light/lambda_; // carrier
  chirp_freq+=l1b->fast_csf; // add chirp start freq
  chirp_freq+=l1b->slow_cfs*(l1b->csq-1)/2; // add half of chirp bandwidth
  lambda_chirp=speed_light/chirp_freq; // convert to wavelength
  
  double ce=get_time();
  pts_update_radar_params_time+=ce-cs;

}


double
PointTargetSim::getPulseTimeAndNumber(double t_return_in_s, 
				      const Dvec& delay_in_s,
				      unsigned int& pulsenum, 
				      int first_sample)
{
  double cs=get_time();
  double t_in_pulse;
  if(first_sample) pulsenum=0;
  double pulse_end_time;
  double offset=pri_in_s+pulsenum*pri_in_s;
  while(pulsenum<l1b->pul){
    // pulse end time = rtt of pulse (delay(pulsenum) 
    // + pulse transmit time (pulsenum+1)*pri
    pulse_end_time=delay_in_s(pulsenum)+offset;

    // if return time is before end of pri # pulsenum
    // compute time in pulse and break
    if(t_return_in_s < pulse_end_time){
      // time in pulse is return_time - pulse start time
      // calling routine checks if greater than pulse_width 
      t_in_pulse = t_return_in_s-pulse_end_time+pri_in_s;
      break;
    }

    // otherwise increment pulsenum
    else{
      pulsenum++;
      offset+=pri_in_s;
    }
  }

  // if outside the time range  of the returned echo
  if(pulsenum>l1b->pul || t_in_pulse < 0 ) t_in_pulse=2*pri_in_s; 
      // return value guaranteed to be > pulse_width

  // return time from start of pulse
  double ce=get_time();
  pts_getPTAN_time+=ce-cs;
  return(t_in_pulse);
}










