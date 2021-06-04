static const char rcs_id_sarprocparams_c[] =
  "@(#) $Id: SARProcParams.cpp,v 11.9 2012/11/28 19:01:24 bstiles Exp $";

#include <strings.h>
#include"SARProcParams.h"
#include "DebugInfo.h"
#include "Utils.h"
#include "Constants.h"
#include "BurstData.h"
#include "L1I.h"
#include "Distributions.h"
using std::endl;
using std::cerr;
using std::cout;

SARProcParams::SARProcParams(Config& cfg, int max_sab_counter)
  : data_noise_out_fp(NULL),
    quality_override(NULL), doppler_offset(NULL), range_offset(NULL),
    Xtimevar_coeff_("Xtimevar_coeffs",1),
    fdopc_coeff_("fdopc_coeffs",1,1), rangec_coeff_("rangec_coeffs",1,1),
    flyby_(cfg), gain_cutoff_coeff_("gain_cutoff_coeffs",1,1),
    num_bursts_processed(0),
    num_bursts_skipped(0),
    min_sab_counter_processed(0),
    max_sab_counter_processed(0),
    num_bursts_outside_region(0),
    num_bursts_forced_bad(0),
    num_bursts_forced_good(0),
    num_bursts_wrong_beam(0),
    num_bursts_bad_ckernel(0),
    num_bursts_bad_ephemeris(0),
    num_bursts_off_limb(0),
    num_bursts_bad_downlink(0),
    num_bursts_wrong_rmode(0),
    num_calibration_bursts(0),
    num_bursts_no_pulses(0)

{
  // Read constant portions of Xfactor from config file
  Xcorr_constant=cfg["X_FACTOR_CONSTANT_CORRECTION"];
  transmit_power=cfg["Pt"];
  carrier_frequency = cfg["carrier_frequency"];
  lambda=Uvar("speed_light")/carrier_frequency;
 
  azimuth_res_coeff=cfg["AZIMUTH_RESOLUTION_COEFFICIENT"];
  range_res_coeff=cfg["RANGE_RESOLUTION_COEFFICIENT"];

  Xconstant_=transmit_power*lambda*lambda*Xcorr_constant/(64*pi*pi*pi);

  // Read in mode dependent portions of Xfactor from config file
  // Xmode is indexed by the r_mode value in the SAB
  // This probably should be handled by an enum, but integer assignments
  // to enums are platform dependent, so the conversion from r_mode to
  // enum must be performed explicitly anyway.

  // Compressed Scatteromet mode will not be handled correctly, but
  // it would not be worthwhile to run that data through the sar processor 
  // anyway

  // The mode dependent calibration term is the same with or w/o auto-gain
  // and it is 0 for radiometer only, calibration, and spare modes
  // Scatterometer mode

  Xmode_[0]=cfg["squared_deviation_of_system_noise_input_at_ALTL"];

  // Altimeter mode
  Xmode_[1]=cfg["squared_deviation_of_system_noise_input_at_ALTH"];

  // SAR low res  mode
  Xmode_[2]=cfg["squared_deviation_of_system_noise_input_at_SARL"];

  // SAR high res mode
  Xmode_[3]=cfg["squared_deviation_of_system_noise_input_at_SARH"];

  Xmode_[4]=Xmode_[5]=Xmode_[6]=Xmode_[7]=0; // Passive mode and calibration modes

  // Scatterometer mode with auto gain
  Xmode_[8]=Xmode_[0];

  // Altimeter mode with auto gain
  Xmode_[9]=Xmode_[1];

  // SAR low res mode with auto gain
  Xmode_[10]=Xmode_[2];

  // SAR high res mode with auto gain
  Xmode_[11]=Xmode_[3];

  // Spare modes
  Xmode_[12]=Xmode_[13]=Xmode_[14]=Xmode_[15]=0;


  // Scatterometer mode
  Tsys_[0]=cfg["SYSTEM_TEMPERATURE_ALTL"];

  // Altimeter mode
  Tsys_[1]=cfg["SYSTEM_TEMPERATURE_ALTH"];

  // SAR low res  mode
  Tsys_[2]=cfg["SYSTEM_TEMPERATURE_SARL"];

  // SAR high res mode
  Tsys_[3]=cfg["SYSTEM_TEMPERATURE_SARH"];

  Tsys_[4]=Tsys_[5]=Tsys_[6]=Tsys_[7]=0; // Passive mode and calibration modes

  // Scatterometer mode with auto gain
  Tsys_[8]=Tsys_[0];

  // Altimeter mode with auto gain
  Tsys_[9]=Tsys_[1];

  // SAR low res mode with auto gain
  Tsys_[10]=Tsys_[2];

  // SAR high res mode with auto gain
  Tsys_[11]=Tsys_[3];

  // Spare modes
  Tsys_[12]=Tsys_[13]=Tsys_[14]=Tsys_[15]=0;

  // Time/temperature varying calibration not performed by
  // SAR processor, reads value at closest approach
  // and polynomial coefficients
  // so for right now read a constant from the config file

  // Read Polynomial coefficients for time-varying X
  
  Xtimevar_poly_order_=cfg.getInt("XFACTOR_TIME_POLYNOMIAL_ORDER");
  Xtimevar_coeff_.resize(Xtimevar_poly_order_+1);

  // for now set constant portion of Xtimevar to zero
  if(Xtimevar_poly_order_==0){
    Xtimevar_coeff_(0)=1;
  }
  else{
    for(int c=0;c<=Xtimevar_poly_order_;c++){
      Xtimevar_coeff_(c)=cfg["XFACTOR_TIME_POLYNOMIAL_COEFF_"+toStr(c)];
    }
  }
  // initialize Xtimevar to constant portion
  Xtimevar_=Xtimevar_coeff_(0);

  // Read beam dependent portion of Xfactor from config file
  for(int b=0;b<5;b++)
    Xbeam_[b]=cfg["X_FACTOR_BEAM_"+toStr(b+1)];

  // Read in noise subtraction offset and initialize parameters
   noise_subtraction_on = (bool) cfg.getInt("ENABLE_NOISE_SUBTRACTION");
   if(!noise_subtraction_on){
     fprintf(stderr,"Warning: Processing without Noise Subtraction!!\n");
     fprintf(stdout,"Warning: Processing without Noise Subtraction!!\n");
     
   }
   quant_noise_scale=1.0;
   quant_noise_offset=0;
   thermal_noise_offset=0;

  // Read in SAR processor options from config file
  single_beam_mode = (bool) cfg.getInt("SINGLE_BEAM_MODE");
  if(single_beam_mode)
    beam_number=cfg.getInt("SARPROC_BEAM_NUMBER");
  process_SAR_mode_only = (bool) cfg.getInt("PROC_SAR_MODE_ONLY");
  if(cfg.keywordExists("PROC_SCAT_MODE_ONLY")){
    process_scat_mode_only = (bool) cfg.getInt("PROC_SCAT_MODE_ONLY");
  }

  else{
    process_scat_mode_only = false;
  }
  if(cfg.keywordExists("PROC_ALT_MODE_ONLY")){
    process_alt_mode_only = (bool) cfg.getInt("PROC_ALT_MODE_ONLY");
  }
  else{
    process_alt_mode_only = false;
  }
  if(cfg.keywordExists("PROC_SARL_MODE_ONLY")){
    process_SARL_mode_only = (bool) cfg.getInt("PROC_SARL_MODE_ONLY");
  }
  else{
    process_SARL_mode_only = false;
  }
  if(cfg.keywordExists("PROC_SARH_MODE_ONLY")){
    process_SARH_mode_only = (bool) cfg.getInt("PROC_SARH_MODE_ONLY");
  }
  else{
    process_SARH_mode_only = false;
  }


  bool SL=process_SARL_mode_only;
  bool SH=process_SARH_mode_only;
  bool S=process_SAR_mode_only;
  bool R=process_scat_mode_only;
  bool A=process_alt_mode_only;
  bool bad_mode=(SL&&(SH||R||A)) || (SH&&(R||A)) || (R&&(A||S)) || (A && S);
  if(bad_mode){
    ErrorMessage e("Fatal Error: Bad SAR processor radar mode selection");
    e.throwMe();
  }
  strip_processor_on = (bool) cfg.getInt("STRIP_PROC_ON");   
  azimuth_range_orthogonalization_on= 
    (bool) cfg.getInt("AZIMUTH_RANGE_ORTHOGONALIZATION_ON");
  azimuth_deramp_processing_on= 
    (bool) cfg.getInt("AZIMUTH_DERAMP_PROCESSING_ON");

  calibration_on=
    (bool) cfg.getInt("SAR_CALIBRATION_ON");

  if(cfg.keywordExists("FORCE_BIDR_REGION")){
    force_bidr_region=
      (bool) cfg.getInt("FORCE_BIDR_REGION");
    forced_first_line=cfg.getInt("FORCED_FIRST_BIDR_LINE");
    forced_last_line=cfg.getInt("FORCED_LAST_BIDR_LINE");
    cerr << "Warning Implementing Special BIDR Region forcing ... " << endl;
    cerr << "warning(cont) Only lines " << forced_first_line 
	 << " through " << forced_last_line << " will be outputted." << endl;
  }
  else force_bidr_region=false;

  // Set maximal relative gain correction to apply 
  // I.e. for all points below maximal_gain_corr dB down from the peak of
  // the oneway gain pattern the 2*maximal_gain_corr dB or correction 
  // will be  applied in the calibration and no more. This is used to avoid
  // garbage leaking into the edges of the useable area through the sinc
  // interpolator, If this number is too high we can get this garbage.
  // If it is too low, features near the edge will be darkened and slightly
  // squeezed together (due to losing energy on one side of the sinc).
  maximal_gain_corr=cfg.getDouble("MAXIMAL_ONEWAY_GAIN_CORRECTION_IN_DB");
  
  if(maximal_gain_corr<0){
    ErrorMessage e("Fatal Error: MAXIMAL_ONEWAY_GAIN_CORRECTION_IN_DB is negative.");
    e.throwMe();
  }
  if(maximal_gain_corr<5){
    cerr << "Warning: MAXIMAL_ONEWAY_GAIN_CORRECTION_IN_DB < 5.0 dB" <<endl;
  }
  // Read Polynomial coefficients for time-varying Gain Cut-off Threshold
  // for computing usable extent
  // gain cutoff polynomial computes one-way gain cut-off in dB
  // but getGain2Thresh() method returns two-way threshold linear scale
  // the model is a polynomial in time from closest approach
  gain_cutoff_enabled=(bool)cfg.getInt("GAIN_CUTOFF_ENABLED");
  if(gain_cutoff_enabled){
    gain_cutoff_poly_order_=cfg.getInt("GAIN_CUTOFF_POLYNOMIAL_ORDER");
    gain_cutoff_coeff_.resize(5,gain_cutoff_poly_order_+1);
    for(int b=0;b<5;b++){
      for(int c=0;c<=gain_cutoff_poly_order_;c++){
	gain_cutoff_coeff_(b,c)=cfg.getDouble("GAIN_CUTOFF_BEAM_" +
					      toStr(b+1)+ 
				  "_POLYNOMIAL_COEFF_"+toStr(c));
      }
    }
  }
  azimuth_width_cutoff_enabled=(bool)cfg.getInt("AZIMUTH_WIDTH_CUTOFF_ENABLED");
  if(azimuth_width_cutoff_enabled){
    azimuth_width_percent=(float)cfg.getDouble("AZIMUTH_WIDTH_PERCENT");
  }

  full_useable_calc_enabled=(bool)cfg.getInt("FULL_USABILITY_CALC_ENABLED");

  remove_doppler_phase_ramp=(bool)cfg.getInt("REMOVE_DOPPLER_PHASE_RAMP");

  enable_data_dependent_cal=(bool)cfg.getInt("ENABLE_DATA_DEPENDENT_CAL");

  // Configure SARAmb object if desired
  if(full_useable_calc_enabled){
    if(gain_cutoff_enabled){
      ErrorMessage e("GAIN_CUTOFF_ENABLED and FULL_USABILITY_CALC_ENABLE cannot both be set.");
      e.throwMe();
    }
    use_multilook_ambiguity=(bool)cfg.getInt("use_multilook_ambiguity");
    amb_sar.config(cfg);
    amb_sar.setTarget(default_target_name
		      ,Frame(default_target_frame_spice_id,
			     default_target_spice_id));
    Frame ftrack=flyby_.trackFrame();
    amb_sar.setTrackFrame(ftrack);
  }
  string sigtypestr;
  if(cfg.keywordExists("SAR_PROC_SIGNAL_TYPE")){
    sigtypestr=cfg.str("SAR_PROC_SIGNAL_TYPE");
    if(sigtypestr=="THEOR_NOISE"){
      signal_type=THEOR_NOISE;
    }
    else if(sigtypestr=="DATA_NOISE"){
      signal_type=DATA_NOISE;
      string filename=cfg.str("DATA_NOISE_OUTFILE");
      if(filename!="none" && filename!="NONE" && filename!="None"){
	data_noise_out_fp=fopen(filename,"w");
      }
    }
    else signal_type=NOMINAL;

 
  }
  else{
    signal_type=NOMINAL;
  }
  simBaqForNoise=false;
  if(signal_type!=NOMINAL){
    cerr << "Warning:: Fudging data to create noise only case of type "
	 << sigtypestr << endl;
    cout << "Warning:: Fudging data to create noise only case of type "
	 << sigtypestr << endl;
    simBaqForNoise=(bool)cfg.getInt("SIM_BAQ_FOR_NOISE");
  }
  string region_str=cfg.str("SAR_PROC_SEGMENT_TYPE");
  if(region_str=="FULL") region_type_=FULL;
  else if(region_str=="SAB_COUNTER"){
    region_type_=BURSTNO;
    burstno_start=cfg.getInt("SAR_PROC_START_SAB");
    burstno_end=cfg.getInt("SAR_PROC_END_SAB");
  }
  else if(region_str=="ABSOLUTE_TIME"){
    region_type_=ABSTIME;
    abstime_start=cfg.getTime("SAR_PROC_START_TIME");
    abstime_end=cfg.getTime("SAR_PROC_END_TIME");
  }
  else if(region_str=="TIME_FROM_CLOSEST_APPROACH"){
    region_type_=CLATIME;
    reltime_start=cfg.getTime("SAR_PROC_START_TIME");
    reltime_end=cfg.getTime("SAR_PROC_END_TIME");
  }
  else if(region_str=="TIME_FROM_EPOCH"){
    region_type_=EPOCHTIME;
    reltime_start=cfg.getTime("SAR_PROC_START_TIME");
    reltime_end=cfg.getTime("SAR_PROC_END_TIME");
  }
  else if(region_str=="TIME_FROM_TRIGGER"){
    region_type_=TRIGGERTIME;
    reltime_start=cfg.getTime("SAR_PROC_START_TIME");
    reltime_end=cfg.getTime("SAR_PROC_END_TIME");
  }
  else if(region_str=="LONLAT" || region_str=="LATLON"){
    region_type_=LATLON;
    lat_south=cfg["SAR_PROC_MIN_LAT"];
    lat_north=cfg["SAR_PROC_MAX_LAT"];
    elon_west=cfg.convertWestLongitude("SAR_PROC_WESTERNMOST_LON");
    elon_east=cfg.convertWestLongitude("SAR_PROC_EASTERNMOST_LON");
  }
  else{
    ErrorMessage e("SARProcParams::config bad SAR_PROC_SEGMENT_TYPE "
		   +region_str);
    e.throwMe();
  }

  // allocate and initialize quality override array;
  num_sabs=max_sab_counter; // allows sab_counter to start from 0
  if(num_sabs<=0){
    ErrorMessage e("SARProcParams Attempt to allocate quality_override array with 0 or negative size="+toStr(num_sabs)+"\nMaybe sab_counter in last record of L1B is bad.");
    e.throwMe();
  }
  quality_override=(char*)make_array(sizeof(char),1,num_sabs);

  if(quality_override==NULL){
    ErrorMessage e("SARProcParams insufficent memory to allocate quality_override array (size="+toStr(num_sabs)+")\nMaybe sab_counter in last record of L1B is bad.");
    e.throwMe();
  }

  for(int c=0;c<num_sabs;c++) quality_override[c]=-1;

  // Read DATA_QUALITY_OVERRIDE_FILE if it exists
  if(cfg.keywordExists("DATA_QUALITY_OVERRIDE_FILE")){
    string file=cfg.str("DATA_QUALITY_OVERRIDE_FILE");
    if(strcasecmp("none",file.c_str())!=0) 
      readQualityOverrideFile(file);
  }


  // use old linear approximation of doppler and doppler rate variation
  // with range
  linear_dopran_approx=false;
  if(cfg.keywordExists("SARPROC_LINEAR_DOPPLER_VS_RANGE")){
    linear_dopran_approx=(bool)cfg.getInt("SARPROC_LINEAR_DOPPLER_VS_RANGE");
    if(linear_dopran_approx){
      cerr << "Warning: Using Linear approximation for doppler and doppler rate variation over range." << endl;
    }
  }

  // Set up configured doppler centroid if desired
  // for now it reads a polynomial for doppler centroid
  // and sets the derivative with respect to range to zero
  // This mode should be used for debugging only so it prints a warning.
  // to stderr and stdout

  set_doppler_centroid=false;
  set_burst_range_dop=false;
  if(cfg.keywordExists("BURST_RANGE_DOPPLER_OFFSETS_FILE")){
    string brfile=cfg.str("BURST_RANGE_DOPPLER_OFFSETS_FILE");
    set_doppler_centroid=true;
    set_burst_range_dop=true;
    doppler_offset=(float*)malloc(sizeof(float)*num_sabs);
    range_offset=(float*)malloc(sizeof(float)*num_sabs);
    for(int c=0;c<num_sabs;c++){
      range_offset[c]=-1000;
    }
    readBurstRangeAndDopplerFile(brfile);
  }

  if(cfg.keywordExists("SET_DOPPLER_CENTROID")){
    set_doppler_centroid=(bool)cfg.getInt("SET_DOPPLER_CENTROID");
    if(set_doppler_centroid){

      // Read Polynomial coefficients for time-varying fdopc and rangec
  
      fdopc_poly_order_=cfg.getInt("DOPPLER_CENTROID_POLYNOMIAL_ORDER");
      fdopc_coeff_.resize(5,fdopc_poly_order_+1);
      fdopc_mean=cfg["DOPPLER_CENTROID_MEAN"];
      fdopc_std=cfg["DOPPLER_CENTROID_STD"];
      time_mean=cfg["ESTRD_TIME_MEAN"];
      time_std=cfg["ESTRD_TIME_STD"];
      for(int b=0;b<5;b++){
	
	for(int c=0;c<=fdopc_poly_order_;c++){
	  fdopc_coeff_(b,c)=cfg["DOPPLER_CENTROID_BEAM_" +
					      toStr(b+1)+ 
				"_POLYNOMIAL_COEFF_"+toStr(c)];
	}
      }


      rangec_poly_order_=cfg.getInt("RANGE_CENTROID_POLYNOMIAL_ORDER");
      rangec_coeff_.resize(5,rangec_poly_order_+1);
      rangec_mean=cfg["RANGE_CENTROID_MEAN"];
      rangec_std=cfg["RANGE_CENTROID_STD"];

     for(int b=0;b<5;b++){
      for(int c=0;c<=rangec_poly_order_;c++){
	rangec_coeff_(b,c)=cfg["RANGE_CENTROID_BEAM_" +
			       toStr(b+1)+ 
			       "_POLYNOMIAL_COEFF_"+toStr(c)];
      }
     }
    }
  }

     
  
}



void SARProcParams::topoInit(float oneway_gain_cutoff, float prf_percent){
  noise_subtraction_on=1;
  single_beam_mode=0;
  strip_processor_on=1;
  simBaqForNoise=0;
  azimuth_width_cutoff_enabled=true;
  azimuth_width_percent=prf_percent;
  gain_cutoff_enabled=true;
  gain_cutoff_poly_order_=0;
  gain_cutoff_coeff_.resize(5,1);
  for(int b=-0;b<5;b++)gain_cutoff_coeff_(b,0)=oneway_gain_cutoff;
  full_useable_calc_enabled=false;
  if(set_burst_range_dop){
    free(doppler_offset);
    free(range_offset);
  }
  
  // Range centroid processing needs to be off but doppler centroid
  // processing may still be allowed BWS Oct 2 2012
  if(set_doppler_centroid==true){
    rangec_poly_order_=0;
    rangec_mean=Uvar(0,"km");
    rangec_std==Uvar(1,"km");
    for(int b=0;b<5;b++){
      rangec_coeff_(b,0)=0;
    }
  }
  
  set_burst_range_dop=false;
}


// Routine for estimating noise and filling up the buffer with noise
// The operation of the routine is controlled by signal_type
// If NOMINAL nothing is done
// If THEOR_NOISE, a Theoretical Gaussian Noise is produced from Tsys
// If DATA_NOISE, uses a precomputed noise_var_estimate (in dn)

 void SARProcParams::replaceSignalWithNoise(int radar_mode, L1I& b, float noise_var_est) const
{
  if(signal_type==NOMINAL) return;

  // compute noise variance estimate analytically
  // otherwise used input value
  else if (signal_type==THEOR_NOISE){
    Uvar Pn=boltzmann_constant*Tsys_[radar_mode]*b.rc_bw;
    double cag=b.computeCalibratedAttenuatorGain();
    noise_var_est=get_in_base_units(Pn*Xmode_[radar_mode])*cag;
    
  }

  Gaussian g(noise_var_est,0);
  for(unsigned int c=0;c<b.Nradar_data;c++){
    b.radar_data[c]=g.GetNumber();
  }  

} 

// Computes estimate of qunatization noise floor, thermal noise floor,
// and scaling constant required to preserve unity gain after removal
// of quantization noise floor
void SARProcParams::estimateNoiseSubtractionParameters(L1I& b, RasList* raslist){
  DebugInfo dbg ("SARProcParams::estimateNoiseSubtractionParameters");
  // compute sample variance of smallest possible signal in Dn
  // to get quantization noise floor.
  // probably should be done block by block
  // for now one quant_noise and thermal_noise estimate for the whole echo
  // is used.
  if(dbg.level) dbg.file << "SARProcParams:estimateNoiseSubtractionParameters" 
			 << endl ;
  if(!raslist){
    cerr << "SARProcParams::estimateNoiseSubtractionParameters Bad raslist pointer" << endl;
    exit(1);
  }

  // Estimate thermal noise
  float var=b.rms_radar_data;
  var*=var;
  float thermal_fudge_factor=0.95;
  double sarproc_gain=b.getSARProcGain();
  Uvar Pn=boltzmann_constant*Tsys_[b.r_mode]*b.rc_bw;
  double cag=b.computeCalibratedAttenuatorGain();
  if(dbg.level) dbg.file << "   CAG computed" << endl;
  float thermal_var=get_in_base_units(Pn*Xmode_[b.r_mode])*cag;
  thermal_noise_offset=thermal_var*sarproc_gain*thermal_fudge_factor;
   
 
  quant_noise_scale=1.0;
  quant_noise_offset=0.0;

  // If BAQ is used estimate quantization noise
  // we used the commented out code
  // for ENCELADUS SAR where nadir is in the main lobe
  // and thus quantization noise is poorly modeled
  // the exclusion for hisar calibrate was left commented out for RHEA SAR
  // if(b.baq_mode!=5 && !b.use_HISAR_calibrate){
  bool target_not_enceladus=(strcasecmp(default_target_name.c_str(),"enceladus"))!=0;
  
  if(b.baq_mode!=5 & target_not_enceladus){
    float quant_var=b.minimumBAQThresholdVariance(raslist);
    float baq_outp_rat=b.estimateBaqOutputRatio(raslist); 
    float true_var=var/baq_outp_rat;
    float true_SNR=(true_var-thermal_var)/thermal_var;
    if(true_SNR<0) true_SNR=0;
    float adj_quant_var=(0.5*quant_var+true_SNR*quant_var)/(1+true_SNR);	  
    quant_noise_offset=sarproc_gain*adj_quant_var;	  
    quant_noise_scale=baq_outp_rat*(1-adj_quant_var/true_var);
 }
  else{
    cout << "Warning BAQ contrast correction is turned OFF" << endl;
  }
  
  // quant noise gain is sarproc_gain (at least for high SNR)
  // thermal noise gain is 2*sarproc_gain
 

  if(dbg.level) dbg.file << "   ThOff " << thermal_noise_offset << " QOff "
			 << quant_noise_offset << " QScale " 
			 << quant_noise_scale << endl;
}
 
Uvar SARProcParams::NoiseVarianceToSystemTemperature(float noise_var,L1I& b) const
{
    double cag=b.computeCalibratedAttenuatorGain();
    Uvar retval = noise_var/(Xmode_[b.r_mode]*cag*b.rc_bw*boltzmann_constant);
    return(retval);
}
//  Main routine to determine whether or not to process a burst
//  called by BIDR::BIDR for the purpose of setting grid boundary
//  and compiled burst skip statistics
//  This routine resets the overrides so that the main SAR processor
//  can call a simplified version.
//   
//
//  First quality overrides are checked if 0 skipped if 1 process
//  if -1 proceed with other checks
//  Others checks are : 
//  1) check BAD SAB bit in engineering quality flag and/or missing geometry
//  2) check beam_number if single beam mode
//  3) check radar mode and calibration mode (processSARonly flag used)
//     to determine whether or not to process other radar modes
//  4) region restrictions are employed
//   
bool SARProcParams::skipBurst(const OblCylProj& proj, const BurstData& bd){

  DebugInfo dbg("SARProcParams::skipBurst");
  
  
  int si=bd.sab_counter-1;
  if(si >= num_sabs || si< 0){
    ErrorMessage e("SARProcParams::skipBurst: sab_counter " + toStr(si+1) +
		   " outside range [1," + toStr(num_sabs) + "]");
    e.throwMe();
  }
  

  int qo=quality_override[si];

  // Burst marked BAD in QUALITY_OVERRIDE_FILE
  if(qo==0){
    num_bursts_skipped++;
    num_bursts_forced_bad++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Quality_Override=BAD" << endl;

    return(1); 
  }

  // Burst marked GOOD or 4DB (and thus GOOD) in QUALITY_OVERRIDE_FILE
  if(qo>0){
    num_bursts_processed++;
    num_bursts_forced_good++;
    if(min_sab_counter_processed==0) min_sab_counter_processed=si+1;
    max_sab_counter_processed=si+1;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " processed, Quality_Override=GOOD" << endl;
    return(0); 
  }

  // Skip beams not processed 
  if(single_beam_mode && (int)bd.beam_number != beam_number){ 
    num_bursts_wrong_beam++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Beam deselected" << endl;
    return(1);
  }

  // Skip bursts outside of desired processing region
  if(!inRegion(proj,bd)){
    num_bursts_outside_region++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Outside Region" << endl;
    return(1);
  }

  // Skip bursts with bad geometry or corrupted SABs and WARN!!
  if(bd.isCal()){
    quality_override[si]=0;
    num_calibration_bursts++;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Calibration burst" << endl;
    return(1);
  }

  // Skip Non-SAR bursts if desired
  if(process_SAR_mode_only && !bd.isSAR()){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Non-SAR burst" << endl;
    return(1);
  }


  // Skip all but SARL mode bursts if desired
  if(process_SARL_mode_only && (bd.adc < Uvar(900000,"Hz") || bd.adc > Uvar(1200000,"Hz"))){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Non-SARL burst" << endl;
    return(1);
  }

  // Skip all but SARH mode bursts if desired
  if(process_SARH_mode_only && (bd.adc < Uvar(1200000,"Hz") || bd.adc > Uvar(2100000,"Hz"))){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Non-SARH burst" << endl;
    return(1);
  }

  // Skip all but ALTH mode bursts if desired
  if(process_alt_mode_only && (bd.adc < Uvar(4000000,"Hz") || bd.adc > Uvar(15000000,"Hz"))){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Non-ALT burst" << endl;
    return(1);
  }

  // Skip all ALTH mode bursts unless we are specifically using that mode (It would bomb anyway due to array sizing!)
  if(!process_alt_mode_only && bd.adc > Uvar(4000000,"Hz")){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, ALTH burst" << endl;
    return(1);
  }

  // Skip all but SCAT mode bursts if desired
  if(process_scat_mode_only && (bd.adc < Uvar(100000,"Hz") || bd.adc > Uvar(300000,"Hz"))){
    num_bursts_wrong_rmode++;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, Non-SCAT burst" << endl;
    return(1);
  }

  // The rest of the burst skip criteria are RARE and therefore generate
  // warnings

  // skip data with no pulses received
  if(bd.num_pulses_received==0){
    num_bursts_no_pulses++;
    cerr << "Warning SAB counter "<<si+1<< " skipped, no pulses received." << endl;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, no pulses received" << endl;
    return(1);
  }

  // skip data with bad ckernel
  if(bd.quality_flag & 0x01){
    num_bursts_bad_ckernel++;
    cerr << "Warning SAB counter "<<si+1<< " skipped, ckernel gap" << endl;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, ckernel gap" << endl;
    return(1);
  }

  // skip data with bad ephemeris
  if(bd.quality_flag & 0x02){
    num_bursts_bad_ephemeris++;
    cerr << "Warning SAB counter "<<si+1<< " skipped,  ephemeris gap" << endl;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, ephemeris gap" << endl;
    return(1);
  }

  // skip data corrupted during downlink
  if(bd.quality_flag & 0x20){
    num_bursts_bad_downlink++;
    cerr << "Warning SAB counter "<<si+1<< " skipped,  bad downlink" << endl;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, bad downlink" << endl;
    return(1);
  }

  // skip data off planet or on the limb
  if(bd.science_qual_flag & 0x0180){
    num_bursts_off_limb++;
    cerr << "Warning SAB counter "
	 <<si+1<< " skipped,  off planet or on limb." << endl;
    quality_override[si]=0;
    num_bursts_skipped++;
    if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " skipped, off planet or on limb" << endl;
    return(1);
  }

  
  // Process all other bursts
  num_bursts_processed++;
  if(min_sab_counter_processed==0) min_sab_counter_processed=si+1;
  max_sab_counter_processed=si+1;
  if(dbg.level) dbg.file << "SPP:skipBurst SAB #"<<si+1<< " processed normally" << endl;
  return(0);
}
			      
void SARProcParams::reportSkippedBursts(ofstream& afs){
  afs << endl;
  afs << "Burst skip statistics Report" << endl;
  afs << "-----------------------------------------------------------" << endl;
  afs << "Total Number of Bursts processed: " << num_bursts_processed << endl;
  afs << "In sab_counter range ("<< min_sab_counter_processed << ","
      << max_sab_counter_processed << ")" << endl;
  afs << "Number of Bursts processed due to quality override:" 
      << num_bursts_forced_good << endl;
  afs << endl;
  afs << "Total Number of burst skipped: " << num_bursts_skipped << endl;
  afs << "Number of bursts skipped due to quality override: " << num_bursts_forced_bad << endl;
  afs << "Number of bursts skipped outside region: " << num_bursts_outside_region << endl;
  afs << "Number of bursts skipped due to deselected beam: " << num_bursts_wrong_beam << endl;
  afs << "Number of bursts skipped due to ckernel_gap: " << num_bursts_bad_ckernel << endl;
  afs << "Number of bursts skipped due to ephemeris gap: " << num_bursts_bad_ephemeris << endl;
  afs << "Number of bursts skipped due to footprint off planet: " << num_bursts_off_limb << endl;
  afs << "Number of bursts skipped due to bad downlink: " << num_bursts_bad_downlink << endl;
  afs << "Number of bursts skipped due to no pulses received: " << num_bursts_no_pulses << endl;
  afs << "Number of calibration bursts skipped: " << num_calibration_bursts << endl;
  afs << "Number of nonSAR bursts skipped: " << num_bursts_wrong_rmode << endl;
}

bool SARProcParams::skipBurst(int sab_counter){
  int si=sab_counter-1;
  int qo=quality_override[si];

  // Burst marked BAD previously
  if(qo==0){
    return(1); 
  }

  // Burst marked GOOD previously
  return(0); 

}

SARProcParams::~SARProcParams()
{

  if(data_noise_out_fp){
    fclose(data_noise_out_fp);
    data_noise_out_fp=NULL;
  }
  // free quality_override array
  if(quality_override!=NULL){
   free_array((void*)(quality_override),1,num_sabs);
  }
  if(doppler_offset!=NULL) free(doppler_offset);
  if(range_offset!=NULL) free(range_offset);

}

bool SARProcParams::isUseable(double fdop_in_Hz, double range_in_km)
{
  return(amb_sar.isGoodBin(range_in_km, fdop_in_Hz,use_multilook_ambiguity));
}
double SARProcParams::getAmbigRatio(double fdop_in_Hz, double range_in_km)
{
  return(amb_sar.getAmbigRatio(range_in_km, fdop_in_Hz,use_multilook_ambiguity));
}
void  SARProcParams::initUseableCalculation(const BurstData& data,
					    const StateVector& sc_state,
					    const Uvar& bore_range,
					    const Uvar& bore_doppler,
					    double lambda_in_km,
					    double look_repeat_time_in_s)
					   
{
  DebugInfo dbg("SARProcParams::initUseableCalculation");
   amb_sar.clear();
   amb_sar.setLambda(lambda_in_km);
   Time t=sc_state.time();
   Uvar range_width=data.pri*speed_light/2.0;
   Uvar prf=1/data.pri;
   Uvar chirp_bandwidth=data.csq*data.slow_cfs;
   amb_sar.setTime(t);
   amb_sar.setState(sc_state);
   amb_sar.setProcessWindow(bore_range, bore_doppler,
			    -range_width/2.0,
			    range_width/2.0,
			    -prf/2.0,
			    prf/2.0,
			    prf);
   amb_sar.setBeamNumber(data.beam_number);
   amb_sar.calAmbGeometry();

   // Compute ground range and ground azimuth resolutions
   Frame beam_frame(beam_frame_spice_id[data.beam_number-1],
		   cassini_spice_id);  
   FloatVector velocity = sc_state.velocity();
   velocity.representIn(beam_frame);
   DirectionVector velocity_dir = velocity;
   double et;
   t.getEt(et);
   DirectionVector bore(beam_frame,et,0.0,0.0,1.0);
   Uvar speed = velocity.magnitude();
   Uvar Vst_bore = velocity.magnitude();
   Vst_bore *= sqrt(1.0 - pow(dot(velocity_dir,bore),2));
   Uvar receive_window = double(data.num_pulses_received) * data.pri;
   Uvar x_res = lambda * bore_range/(2.0 * Vst_bore * receive_window);
   Uvar rg_res = speed_light/(2.0 * chirp_bandwidth)/sin(data.act_incidence_angle);


   // Compute NoisePower Pn and X0
   Uvar Nr = data.chirp_length * chirp_bandwidth; 
   Uvar Ni = double(data.num_pulses_received)*Nr;
   Uvar X0 = transmit_power * Ni/pow(4.0*pi,3);
   Uvar Pn=boltzmann_constant * Tsys_[data.r_mode] * data.rc_bw;
   if(!use_multilook_ambiguity) amb_sar.calUsableArea(X0,Pn,x_res,rg_res);
   else{
     amb_sar.calMultilookUsableArea(X0,Pn,x_res,rg_res,look_repeat_time_in_s,velocity,data.num_pulses_received);
   }
   if(dbg.level){
     Uvar seconds_since_start=t-flyby_.startTime();
     Uvar time_from_closest_approach=t-flyby_.epochTime();
     dbg.file << "Time from start of file = " 
	      << seconds_since_start.getInUnits("s") << " s" << endl;
     dbg.file << "Time closest approach = " 
	      << time_from_closest_approach.getInUnits("s") << " s" << endl;
     dbg.file << "Beam Number = " << data.beam_number << endl;

     // print out tables of various amb_sar intermediate values
     dbg.file << "amb_ratio table" << endl;

     int Ngrid=amb_sar.getRangeGridsize();
    
     dbg.file.precision(5);
     dbg.file << "          ";
     for(int d=0;d<(int)Ngrid;d++){
       dbg.file.width(10);
       dbg.file << amb_sar.doppler[d];
     }
     dbg.file << endl;
     for(int r=0;r<(int)Ngrid;r++){
       dbg.file.width(10);
       dbg.file << amb_sar.range[r];
       for(int d=0;d<(int)Ngrid;d++){
	 dbg.file.width(10);
	 dbg.file << amb_sar.amb_ratio_SL[r][d];
       }
       dbg.file << endl;
     }
     dbg.file << endl << endl;




     dbg.file << "good_bin table" << endl;
     dbg.file.precision(5);

     dbg.file << "          ";
     for(int d=0;d<(int)Ngrid;d++){
       dbg.file.width(10);
       dbg.file << amb_sar.doppler[d];
     }
     dbg.file << endl;
     for(int r=0;r<(int)Ngrid;r++){
       dbg.file.width(10);
       dbg.file << amb_sar.range[r];
       for(int d=0;d<(int)Ngrid;d++){
	 dbg.file.width(10);
	 dbg.file << amb_sar.good_bin[r][d];
       }
       dbg.file << endl;
     }
     dbg.file << endl << endl;
     
     
     dbg.file << "snr table" << endl;
     dbg.file.precision(5);
     dbg.file << "          ";
     for(int d=0;d<(int)Ngrid;d++){
       dbg.file.width(10);
       dbg.file << amb_sar.doppler[d];
     }
     dbg.file << endl;
     for(int r=0;r<(int)Ngrid;r++){
       dbg.file.width(10);
       dbg.file << amb_sar.range[r];
       for(int d=0;d<(int)Ngrid;d++){
	 dbg.file.width(10);
	 dbg.file << amb_sar.thermal_snr_SL[r][d];
       }
       dbg.file << endl;
     }
     dbg.file << endl << endl;

     
     
     dbg.file << "One way gain table" << endl;
     dbg.file.precision(5);
     dbg.file << "          ";
     for(int d=0;d<(int)Ngrid;d++){
       dbg.file.width(10);
       dbg.file << amb_sar.doppler[d];
     }
     dbg.file << endl;
     for(int r=0;r<(int)Ngrid;r++){
       dbg.file.width(10);
       dbg.file << amb_sar.range[r];
       for(int d=0;d<(int)Ngrid;d++){
	 dbg.file.width(10);
	 dbg.file << amb_sar.oneway_beam_gaindB[r][d];
       }
       dbg.file << endl;
     }
     dbg.file << endl << endl;
   }
}

float SARProcParams::getGain2Thresh(const Time& time, int beam_num,int sab_counter){
  Uvar time_from_closest_approach=time-flyby_.epochTime();
  double t=time_from_closest_approach.getInUnits("s");
  float gaindb=0;
  for(int c=gain_cutoff_poly_order_;c>0;c--){
    gaindb+=gain_cutoff_coeff_(beam_num-1,c);
    gaindb*=t;
  }
  gaindb+=gain_cutoff_coeff_(beam_num-1,0); // One way gain in db
  // Force 4 db gain cut-off if specified in quality override file
  if(quality_override[sab_counter-1]==4) gaindb=4;
  gaindb*=2; // Two way loss in dB (3 means 3 dB down from the two-way peak)

  return(pow(10.0,-0.1*gaindb)); // convert to linear scale and return  
}
const Uvar&  
SARProcParams::powerToDataNumber(int radar_mode)
{
  return(Xmode_[radar_mode]);
}

const Uvar&  
SARProcParams::getTsys(int radar_mode)
{
  return(Tsys_[radar_mode]);
}
const Uvar&  
SARProcParams::Xcal(const Time& time,
			  int beam_num, int radar_mode)
{

  Uvar time_from_closest_approach;
  DebugInfo dbg("SARProcParams::updateXcal");
    
  // special debugging routine coordinated among L1I, SARProcParams and
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

  if(Xtimevar_poly_order_!=0){

    Uvar ti=1.0;
    Xtimevar_=0;
 
    // compute time from closest approach
    time_from_closest_approach=time-flyby_.epochTime();

    // compute timevarying correction
    for(int c=0;c<=Xtimevar_poly_order_;c++){
      Xtimevar_+=ti*Xtimevar_coeff_(c);
      ti*=time_from_closest_approach;
    }
  }

  // Case for Special_Calibration_Debug
  if(dbgsp.level){
    if ( 0x0040 & dbgsp.level) Xconstant_=1.0;
    if ( 0x0080 & dbgsp.level) Xbeam_[beam_num-1] = 1.0;
    if ( 0x0100 & dbgsp.level) Xtimevar_= 1.0;
    if ( 0x0200 & dbgsp.level) Xmode_[radar_mode]=1.0;
  }

  // compute complete Xcal
  Xcal_=Xconstant_*Xbeam_[beam_num-1]*Xtimevar_*Xmode_[radar_mode];

  // Debugging code
  if(dbg.level){
    dbg.file << "SARProcParams::updateXcal debugging message .. " << endl;
    dbg.file << "    Time (UTC) "<< time.utc("ISOD")
             << "    from closest approach: "
	     << time_from_closest_approach << ", Beam: " << beam_number
	     << ", Radar_Mode: " << radar_mode << endl;
    dbg.file << "Xconstant_: "<< Xconstant_ <<",Xbeam_: " 
	     << Xbeam_[beam_num-1]
	     << ",Xtimevar_: " << Xtimevar_ << ",Xmode_:"<< Xmode_[radar_mode]
	     << ",Xcal_: " << Xcal_;
    
  }
  return(Xcal_);
}


bool SARProcParams::inRegion(const OblCylProj& proj, const BurstData& bd){

  /*** get values from BurstData ***/
  int bn=bd.record_id%1000000;
  Time t=bd.t;
  Uvar ct=bd.time_from_closest_approach;
  Uvar et=bd.time_from_epoch;
  Uvar tt=t - trigger_time_;
  Uvar lat=bd.act_centroid_lat;
  Uvar lon=bd.act_centroid_lon;


  string msg;
  
  switch(region_type_){
  case FULL:
    return(true);
    break;
  case BURSTNO:
    return(bn<=burstno_end && bn>=burstno_start);
    break;
  case ABSTIME:
    return(t<=abstime_end && t>=abstime_start);
    break;
  case CLATIME:
   return(ct<=reltime_end && ct>=reltime_start);
    break;
  case EPOCHTIME:
   return(et<=reltime_end && et>=reltime_start);
    break;
  case TRIGGERTIME:
   return(tt<=reltime_end && tt>=reltime_start);
    break;
  case LATLON:
    if(lat>=lat_south && lat<=lat_north){
      // choose lon boundaries where elon_west < elon_east <= elon_west+ 2pi
      
      while(elon_west>elon_east)  elon_west-=Uvar(2*pi,"rad");
      while(elon_east>elon_west+Uvar(2*pi,"rad")) elon_east-=Uvar(2*pi,"rad");

      // to avoid floating point error problem make very small 
      // or very large ranges
      // equal to full range and thus burst is in region
      if(elon_east-elon_west<Uvar(0.00001,"rad")) return(true);
      if(elon_east-elon_west>Uvar(2*pi-0.00001,"rad")) return(true);

      // check if lon is in range
      while(lon>elon_east)lon-=Uvar(2*pi,"rad");
      while(lon<elon_west)lon+=Uvar(2*pi,"rad");
      if(lon>elon_east) return(false);
      else return(true);
    }
    else return(false);
    break;
  case LINENO:
    msg="LINENO region type not yet implemented";
    break;
  default:
    msg="INVALID region type (SHOULD NEVER HAPPEN)";
    break;
  }
  ErrorMessage e("SARProcParams::inRegion "+msg);
  e.throwMe();
  return(false);
}

void 
SARProcParams::setDopplerProfile(const Time& time, int beam_num,
				 Uvar& slope, Uvar& doppler_boresight, 
				 Uvar& range_boresight,
				 Uvar& rate_slope,
				 Uvar& doppler_rate_boresight,
				 int sab_counter) const
{
  DebugInfo dbg("SARProcParams::setDopplerProfile");
  if(dbg.level) dbg.file << "Start debugging SARProcParams::setDopplerProfile"
			 << endl;

  if(set_burst_range_dop){
    slope=0;
    rate_slope=0;
    doppler_rate_boresight=0;
    if(sab_counter==-1){
      ErrorMessage e("SAB COUNTER not set in call to setDopplerProfile");
      e.throwMe();
    }
    doppler_boresight=doppler_offset[sab_counter];
    range_boresight=range_offset[sab_counter];
    if(range_boresight<-999){
      ErrorMessage e("Range Offset not computable (SAB missing in BurstRangeDop file?)");
      e.throwMe();
    }
  }
  else{
    Uvar t=time-flyby_.epochTime();
    if(dbg.level) dbg.file << "Time from Closest approach = " << t << endl;

    t-=time_mean;
    t/=time_std;
    if(dbg.level) dbg.file << "Norm time " << t <<  endl;

    // set nominal doppler and doppler rate ddop/dt
    doppler_boresight=0;
    for(int c=fdopc_poly_order_;c>0;c--){
      doppler_boresight+=fdopc_coeff_(beam_num-1,c);
      doppler_boresight*=t;
      if(c>1){
	doppler_rate_boresight+=fdopc_coeff_(beam_num-1,c)*c;
	doppler_rate_boresight*=t;
      }
      else{
	doppler_rate_boresight+=fdopc_coeff_(beam_num-1,1);
      }
    }

    doppler_boresight+=fdopc_coeff_(beam_num-1,0); 

    if(dbg.level) dbg.file << "Norm doppler correction " << doppler_boresight <<  endl;
    
    doppler_boresight*=fdopc_std;
    doppler_boresight+=fdopc_mean;
    if(dbg.level) dbg.file << "Doppler correction " << doppler_boresight <<  endl;
    // set nominal range  
    range_boresight=0;
    
    
    
    for(int c=rangec_poly_order_;c>0;c--){
      
      range_boresight+=rangec_coeff_(beam_num-1,c);
      range_boresight*=t;
    }
    
    range_boresight+=rangec_coeff_(beam_num-1,0);
    if(dbg.level) dbg.file << "Norm range correction " << range_boresight <<  endl;
    range_boresight*=rangec_std;
    range_boresight+=rangec_mean;
    if(dbg.level) dbg.file << "range correction " << range_boresight <<  endl;
    // cout << "range correction " << range_boresight << endl;
    // set range derivatives to zero .....
    slope=0;
    rate_slope=0;
  }
}
//------------------------------------------------------------------
// getClosestApproachState()
//
// Return state vector in Target Body Fixed frame at time of closest
// approach
//------------------------------------------------------------------
StateVector SARProcParams::getClosestApproachState() {

  StateVector s;

  // Construct Target Body Fixed frame
  Frame tbff(default_target_frame_spice_id,default_target_spice_id);

  // Get time and then state vector at time of closest approach
  Time lowest_altitude_time = flyby_.lowestAltitudeTime();
  tbff.ephemeris(s,"Cassini",lowest_altitude_time,"NONE");
  return(s);
}

//------------------------------------------------------------------
// LookDirection()
//
// Return look direction ("Right" or "Left") at time of closest
// approach
//------------------------------------------------------------------
string SARProcParams::LookDirection() {
  return(flyby_.LookDirection());
}

// reads a quality override file
// file contains string constants and integers separated by white space
// "ALL_BAD" string constant must be the first string in the file it
// indicates that all sabs are treated as bad unless otherwise indicated.
//
// "ALL_GOOD" string constant must be the first string in the file it
// indicates that all sabs are treated as good unless otherwise indicated
//
// If neither of these are present in the file the default behavior is to
// apply the standard mode, region, and data quality checks in the SAR 
// processor
//
// "GOOD_RANGE" indicates that the next two white space separated integers in
// the file define an inclusive range of sab_counters for GOOD SABs
//
// "BAD_RANGE" indicates that the next two white space separated integers in
// the file define an inclusive range (min, max ) of sab_counters for BAD SABs
// 
// "4DB_RANGE" is a special case. It uses an inclusive sab counter range as do
// GOOD_RANGE and BAD_RANGE. It sets to quality_override flag to the special 
// value of 4. SABs in this range are processed with a GAIN CUTOFF of 4 dB 
// (one way) to define
// the processing window instead of the nominal manner. All the SABs in
// this range will forcibly included as if they were in a GOOD_RANGE. To force
// exclusion of a subset of the range append a BAD_RANGE later in the file.
// 
// -1 for maximum means extend range to maximum sab_counter value
// 
// Overlapping ranges are disallowed as are setting a range 
// to the default value
// i.e. (ALL_BAD + BAD_RANGE)
//
void SARProcParams::readQualityOverrideFile(string file)
{
  FileMgr f(file,"r");
  char default_setting=-1;
  int line=0;
  while(!f.eof()){
    string keyword="";
    try{
      keyword=f.parseNextString();
    }
    // parseNextString returns a error if there is only whitespace on a line
    // I'll go ahead and ignore that.
    catch(ErrorMessage e){
      if(line==0){
	cerr << "Warning: empty QualityOverrideFile " << file << endl;
      }
      break;
    }
    if(f.eof() && keyword==""){
      if(line==0){
	cerr << "Warning: empty QualityOverrideFile " << file << endl;
      }
      break;
    }

    // ALL_BAD keyword
    if(keyword=="ALL_BAD"){
      // Bad Format Case
      if (line!=0){
	ErrorMessage e("Bad format in QualityOverrideFile " + file +
		       " on line " + toStr(line+1) + 
		       "\nALL_BAD is only allowed on first line.");
	e.throwMe();	       
      }
      else{
        default_setting=0;
        for(int c=0;c<num_sabs;c++) quality_override[c]=default_setting;
      }
    }

    else if(keyword=="ALL_GOOD"){
      // Bad Format Case
      if (line!=0){
	ErrorMessage e("Bad format in QualityOverrideFile " + file +
		       " on line " + toStr(line+1) + 
		       "\nALL_GOOD is only allowed on first line.");
	e.throwMe();	       
      }

      else{
        default_setting=1;
        for(int c=0;c<num_sabs;c++) quality_override[c]=default_setting;
      }
    }

    // GOOD_RANGE and BAD_RANGE keywords
    else if(keyword=="GOOD_RANGE" || keyword=="BAD_RANGE" || keyword=="4DB_RANGE"){
      char val=1;
      if(keyword=="BAD_RANGE") val=0;

      else if(keyword=="4DB_RANGE") val=4; // Way to force 4 db gain cutoff
      // Only works if gain cut-off is enabled in config file
      // Should probably only be used when gain cut-off is nominally a constant
      // 5 dB value

      string mins=f.parseNextString();
      string maxs=f.parseNextString();
      int min_sab=atoi(mins.c_str());
      int max_sab=atoi(maxs.c_str());

      // BAD FORMAT case:
      if(val==default_setting){
	ErrorMessage e("Bad format in QualityOverrideFile " + file +
		       " on line " + toStr(line+1) + 
		       "\nGOOD_RANGE disallowed when ALL_GOOD set and"
		       +"\nBAD_RANGE disallowed when ALL_BAD set.");
	e.throwMe();
      }

      // BAD FORMAT case:
      if(min_sab>max_sab || max_sab > num_sabs || min_sab < 0){
	if(max_sab!=-1){
	  ErrorMessage e("Bad format in QualityOverrideFile " + file +
			 " on line " + toStr(line+1) + 
			 "\nInvalid SAB range ("+toStr(min_sab)+","
			 +toStr(max_sab)+")");
	  e.throwMe();
	}

        // max_sab == -1 case
	else{
	  max_sab=num_sabs;
	}
      }

      // Set range to appropriate value
      for(int c=min_sab-1;c<max_sab;c++){
        // BAD FORMAT CASE
	if(quality_override[c]!=default_setting){
	  ErrorMessage e("Bad format in QualityOverrideFile " + file +
			 " on line " + toStr(line+1) + 
			 "\n SAB range ("+toStr(min_sab)+","
			 +toStr(max_sab)+") Overlaps previous range.");
	  e.throwMe();
	}
	quality_override[c]=val;
      }
    }
    else{
      ErrorMessage e("Bad format in QualityOverrideFile " + file +
			 " on line " + toStr(line+1) +
		     "\n Bad keyword "+keyword+" found.");
      e.throwMe();
    }
    line++;
  }
  f.close();
}

// reads Burst by Burst Doppler and Range File
void
SARProcParams::readBurstRangeAndDopplerFile(string file){
  FileMgr f(file,"r");
  int line=0;
  while(!f.eof()){
    string str="";
    try{
      str=f.parseNextString(); // skip time
    }
    // parseNextString returns a error if there is only whitespace on a line
    // I'll go ahead and ignore that.
    catch(ErrorMessage e){
      if(line==0){
	cerr << "Warning: empty QualityOverrideFile " << file << endl;
      }
      break;
    }
    if(f.eof() && str==""){
      if(line==0){
	cerr << "Warning: empty QualityOverrideFile " << file << endl;
      }
      break;
    }
    str=f.parseNextString(); // SAB
 
    int sab=atoi(str.c_str());
    str=f.parseNextString(); // skip beam
    str=f.parseNextString(); // skip range index
    str=f.parseNextString(); // skip doppler index
    str=f.parseNextString(); // range
    
    range_offset[sab]=atof(str.c_str());
    str=f.parseNextString();  // doppler
    doppler_offset[sab]=atof(str.c_str());
    line++;

    
  }
  f.close();
  // set unset (ie, Beams 1,2,4,5) sabs if necessary
  int last_good=-1;
  for(int c=0;c<num_sabs;c++){
    if(range_offset[c]>-999){

      // fill in early sabs
      if(last_good==-1){
	for(int d=c-1;d>=0;d--){
	  range_offset[d]=range_offset[c];
          doppler_offset[d]=doppler_offset[c];
	}
      }
      // fill in last few sabs (works for 1 2 3 4 5 case)
      else if(!single_beam_mode){
	for(int d=c-2;d<=c-1;d++){
	  range_offset[d]=range_offset[c];
          doppler_offset[d]=doppler_offset[c];
	}
      }
      // interpolate between two values
      else{
	float width=c-last_good;
	for(int d=last_good+1;d<c;d++){
	  float c0=(c-d)/width;
	  float c1=(d-last_good)/width;
	  range_offset[d]=c0*range_offset[last_good]+c1*range_offset[c];
          doppler_offset[d]=c0*doppler_offset[last_good]+c1*doppler_offset[c];
	}
      }
      last_good=c;      
    }
    else{
      range_offset[c]=range_offset[last_good];
      doppler_offset[c]=doppler_offset[last_good];
    }
  }
}
