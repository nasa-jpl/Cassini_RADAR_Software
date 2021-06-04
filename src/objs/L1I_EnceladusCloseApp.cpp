static const char rcs_id_l1i_c[] =
  "@(#) $Id: L1I_EnceladusCloseApp.cpp,v 11.1 2017/10/18 17:57:23 bstiles Exp $";

#include "L1I.h"
#include "Error.h"
#include "SARFunctions.h"
#include "RangeDopplerFile.h"
#include "config_keywords.h"
#include "Constants.h"
#include "DebugInfo.h"
#include "BIDR.h"
#include <string>

// This sould be set to 0 for nominal operations
// set it to 1 for ql1 if there appears to be a very BAD attitude error

using std::cout;
using std::endl;
using std::cerr;


//-----------------
// Methods for L1I
//-----------------

//--------------
// Constructors
//--------------


L1I::L1I(const std::string& filename, const std::string& mode)
  : 
  BurstData(filename,mode,"active"), 
  data_complete(false),
  sigma0fp(NULL), is0fp(NULL), Xfp(NULL), Xambfp(NULL), dcfp(NULL), rfp(NULL), mffp(NULL),
  g2fp(NULL),areafp(NULL),rangefp(NULL), fdopfp(NULL), estrdfp(NULL),
  azimfp(NULL),elevfp(NULL),pixdfp(NULL), rspecfp(NULL),
  use_upper_band(false),
  use_linear_chirp(false),
  use_trivial_pulse_segmentation_(false),
  perform_single_burst_check(false),
  check_this_burst(false)
  
{
  
  header_size_=L1I_HEADER_LENGTH;
  record_size_=L1I_RECORD_LENGTH;

  // append parameters 
  parameters_.appendParameter("rms_radar_data",sizeof(float),
			      Parameter::FLOAT,(void*)(&rms_radar_data));
  parameters_.appendParameter("num_range_bins",num_range_bins);
  parameters_.appendParameter("num_pulses_received",num_pulses_received);
  parameters_.appendParameter("window_start_delay",window_start_delay);
  parameters_.appendParameter("window_start_idx",window_start_idx);


  // Allocate arrays
  allocateArrays();
}

void

L1I::zeroArrays(){

  // zeros all arrays
  // SHOULD BE CALLED BEFORE READING DATA FROM LBDR
  complex<float> zc(0,0);
  for(int i=0;i<L1I_RANGE_WIDTH;i++){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      dechirped[i][j]=zc;
      raw_sar_data[i][j]=zc;
      sigma0[i][j]=zc;
      Xfactor[i][j]=0.0;
      area_buf[i][j]=0.0;
      g2_buf[i][j]=0.0;
      dop_buf[i][j]=0.0;
      azim_buf[i][j]=0.0;
      elev_buf[i][j]=0.0;
    }

    fdop_centroid[i]=0.0;
    fdop_rate[i]=0.0;
    range[i]=0.0;
    matched_filter[i]=zc;
    fft_matched_filter[i]=zc;
    interp_range_buf[i]=zc;    
  }

  for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
    
    interp_doppler_buf[j]=zc;
    azimuth_ref_func[j]=zc;
  }

  for(int c=0;c<GENERIC_ARRAY_WIDTH;c++){
    radar_data[c]=0.0;
    tmp_complex_array[c]=0.0;
    tmp_complex_array2[c]=0.0;
  }
  for(int c=0;c<FFT_INTERP_SIZE;c++){
    fftinterp_fdop_centroid[c]=0;
  }
  for(int c=0;c<FFT_INTERP_SIZE;c++){
    for(int d=0;d<FFT_INTERP_SIZE;d++){
      interp_sigma0[c][d]=zc;
    }
  }
}			     

void
L1I::allocateArrays()
{
  DebugInfo dbg("L1I::allocateArrays");

  long unsigned int bytes=0;
  // allocate 2-D arrays
  dechirped=(complex<float>**)make_array(sizeof(complex<float>),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(complex<float>),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  

  raw_sar_data=(complex<float>**)make_array(sizeof(complex<float>),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(complex<float>),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);  

  sigma0=(complex<float>**)make_array(sizeof(complex<float>),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(complex<float>),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  interp_sigma0=(complex<float>**)make_array(sizeof(complex<float>),2,
					 FFT_INTERP_SIZE,FFT_INTERP_SIZE);
  bytes+=size_array(sizeof(complex<float>),2,FFT_INTERP_SIZE,FFT_INTERP_SIZE);

  Xfactor=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  Xambig=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  ests0=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  ests0ambig=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  area_buf=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  g2_buf=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  dop_buf=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  azim_buf=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);

  elev_buf=(float**)make_array(sizeof(float),2,
					 L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  bytes+=size_array(sizeof(float),2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);


  // allocate 1-D arrays
  radar_data = new float[GENERIC_ARRAY_WIDTH];
  bytes+=GENERIC_ARRAY_WIDTH*sizeof(float);

  fdop_centroid = new double[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(double);

  fftinterp_fdop_centroid = new double[FFT_INTERP_SIZE];
  bytes+=FFT_INTERP_SIZE*sizeof(double);

  fdop_rate = new double[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(double);

  range = new double[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(double);

  azimuth_ref_func = new complex<float>[L1I_DOPPLER_WIDTH];
  bytes+=L1I_DOPPLER_WIDTH*sizeof(complex<float>);

  matched_filter = new complex<float>[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(complex<float>);

  fft_matched_filter = new complex<float>[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(complex<float>);

  tmp_complex_array = new complex<float>[GENERIC_ARRAY_WIDTH];
  bytes+=GENERIC_ARRAY_WIDTH*sizeof(complex<float>);

  tmp_complex_array2 = new complex<float>[GENERIC_ARRAY_WIDTH];
  bytes+=GENERIC_ARRAY_WIDTH*sizeof(complex<float>);

  interp_range_buf = new complex<float>[L1I_RANGE_WIDTH];
  bytes+=L1I_RANGE_WIDTH*sizeof(complex<float>);

  interp_doppler_buf = new complex<float>[L1I_DOPPLER_WIDTH];
  bytes+=L1I_DOPPLER_WIDTH*sizeof(complex<float>);

  if(dbg.level){
    dbg.file << "L1I allocates " << bytes 
	     << " bytes of array + " << sizeof(L1I) << " bytes of object = ";
    bytes+=sizeof(L1I);
    dbg.file << bytes << " bytes Total" << endl;
  }

}  

L1I::~L1I(){

  // Delete output arrays
  free_array((void*)dechirped,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  free_array((void*)raw_sar_data,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  free_array((void*)sigma0,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH);
  free_array((void*)interp_sigma0,2,FFT_INTERP_SIZE,FFT_INTERP_SIZE);
  free_array((void*)Xfactor,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)Xambig,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)ests0,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)ests0ambig,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)area_buf,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)g2_buf,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)dop_buf,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)azim_buf,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 
  free_array((void*)elev_buf,2,L1I_RANGE_WIDTH,L1I_DOPPLER_WIDTH); 

  // Delete intermediate 1-D arrays
  delete[] radar_data;
  delete[] fdop_centroid;
  delete[] fftinterp_fdop_centroid;
  delete[] fdop_rate;
  delete[] range;
  delete[] azimuth_ref_func;
  delete[] matched_filter;
  delete[] fft_matched_filter;
  delete[] tmp_complex_array;
  delete[] tmp_complex_array2;
  delete[] interp_range_buf;
  delete[] interp_doppler_buf;


  // close intermediate data files if necessary
  if(sigma0fp) fclose(sigma0fp);
  if(Xfp) fclose(Xfp);
  if(Xambfp) fclose(Xambfp);
  if(dcfp) fclose(dcfp);
  if(rspecfp) fclose(rspecfp);
  if(rfp) fclose(rfp);
  if(mffp) fclose(mffp);
  if(g2fp) fclose(g2fp);
  if(areafp) fclose(areafp);
  if(rangefp) fclose(rangefp);
  if(fdopfp) fclose(fdopfp);
  if(estrdfp) fclose(estrdfp);
  if(azimfp) fclose(azimfp);
  if(elevfp) fclose(elevfp);
  if(pixdfp) fclose(pixdfp);
  if(is0fp) fclose(is0fp);
} 

void L1I::config(Config& cfg)
{
  // set up cassini sim -- trying to get rid of this
  cassini_sim.config(cfg);

  // set up access to target radii as doubles
  default_target_radii_in_km[0]=default_target_radii[PositionVector::X].getInUnits("km");
  default_target_radii_in_km[1]=default_target_radii[PositionVector::Y].getInUnits("km");
  default_target_radii_in_km[2]=default_target_radii[PositionVector::Z].getInUnits("km");

  if(cfg.keywordExists("USE_QL1_DOPRAN_CORRECTION")){
    quick_look_dopran=(bool)cfg.getInt("USE_QL1_DOPRAN_CORRECTION");
    if(quick_look_dopran){
      fprintf(stderr,"WARNING You are using QL1 Doppler and Range Corrections\n");
      fprintf(stderr,"NEVER us this for final processing.\n");
      fprintf(stderr,"All locations will be messed up by 1-2 km\n");
    }
  }
  else{
    quick_look_dopran=false;
  }



  if(cfg.keywordExists("ORTHORECTIFY")){
    orthorectify=(bool)cfg.getInt("ORTHORECTIFY");
    if(orthorectify){
      string tmapfile=cfg.str("ORTHO_TOPOMAP_FILE");
      ortho_tmap.config(cfg,tmapfile);
      ortho_height_tolerance=cfg["ORTHO_HEIGHT_TOLERANCE"];
      fprintf(stderr,"L1I::config Orthorectification On with height tolerance %g km\n",ortho_height_tolerance.getInUnits('m'));
    }
    
  }
  else{
    orthorectify=false;
  }


 if(cfg.keywordExists("USE_HISAR_CALIBRATION")){
    HISAR_ambig_mode = 0;
    use_HISAR_calibrate=(bool)cfg.getInt("USE_HISAR_CALIBRATION");
    if(use_HISAR_calibrate){
      signal_to_amb_thresh=cfg["signal_to_amb_ratiodB"].getInUnits("");
      min_oneway_gain=cfg["min_oneway_gaindB_wrt_peak"].getInUnits("");
      exclude_mirror=(bool)cfg.getInt("HISAR_EXCLUDE_MIRROR_AMBIGUITY");
      nadir_exclude_ang=cfg["HISAR_NADIR_TRACK_EXCLUSION_ANGLE"];
      HISAR_ambig_mode = cfg.getInt("HISAR_CALIBRATION_AMBIG_MODE");
 
      fprintf(stderr,"L1I::config HISAR Calibration On with SignalToAMb threshold %g dB with ambigmode = %d\n",signal_to_amb_thresh,HISAR_ambig_mode);
      fprintf(stdout,"L1I::config HISAR Calibration On with SignalToAMb threshold %g dB with ambigmode = %d\n",signal_to_amb_thresh,HISAR_ambig_mode);
      signal_to_amb_thresh=pow(10,0.1*signal_to_amb_thresh);
      min_oneway_gain=pow(10,0.1*min_oneway_gain);
    }
    
  }
  else{
    use_HISAR_calibrate=false;
  }

  if(cfg.keywordExists("NADIR_SAR_MODE_ENABLED")){
    nadir_mode=(bool)cfg.getInt("NADIR_SAR_MODE_ENABLED");
  }
  else{
    nadir_mode=false;
  }
  // The REPLACE_X Mode is used to output arbitrary parameters to the X
  // backplane. As a side effect the new X parameter is the tie breaker for
  // choosing which beam to process.
  if(cfg.keywordExists("L1I_REPLACE_X")){
    replace_X_mode=(bool)cfg.getInt("L1I_REPLACE_X");
    if(replace_X_mode){
      string rxp_str=cfg.str("L1I_REPLACE_X_PARAM");
      if(rxp_str=="RANGE_RES") replace_X_param=RANGE_RES;
      else if(rxp_str=="AZIMUTH_RES") replace_X_param=AZIMUTH_RES;
      else if(rxp_str=="LOOKVEL_ANGLE") replace_X_param=LOOKVEL_ANGLE;
      else if(rxp_str=="RES_ELEMENT_ANGLE") replace_X_param=RES_ELEMENT_ANGLE;
      else if(rxp_str=="AMBIGUITY_RATIO") replace_X_param=AMBIG_RATIO;
      else if(rxp_str=="AMBIG_RATIO") replace_X_param=AMBIG_RATIO;
      else if(rxp_str=="ORTHO_LOCAL_RADIUS") replace_X_param=ORTHO_LOCAL_RAD;
      else if(rxp_str=="ORTHO_LOCAL_RAD") replace_X_param=ORTHO_LOCAL_RAD;
      else if(rxp_str=="ORTHO_HEIGHT") replace_X_param=ORTHO_HEIGHT;
      else if(rxp_str=="NONE") replace_X_param=NONE;
      else{
	cerr << "BAD L1I_REPLACE_X_PARAM" << endl;
        exit(1);
      }
      cerr << "Warning L1I_REPLACE_X Activated; X will be fudged." << endl;
      cout << "Warning L1I_REPLACE_X Activated; X will be fudged." << endl;
    }
  }
  else replace_X_mode=false;

  // Set up sinc or FFT interpolation filters
  use_sinc_interpolation_ = cfg.getInt("USE_SINC_INTERPOLATION");
  use_fft_interpolation_ = cfg.getInt("USE_FFT_INTERPOLATION");
  if (use_sinc_interpolation_ && use_fft_interpolation_){
    cerr << "Warning both SINC and FFT interploation enabled ..." << endl;
    cerr << "   disabling SINC." << endl;
    cout << "Warning both SINC and FFT interploation enabled ..." << endl;
    cout << "   disabling SINC." << endl;
    use_sinc_interpolation_=0;
  }
  if (use_sinc_interpolation_)
  {
    if (! sinc_table_range.configSet())
    {
      sinc_table_range.config(cfg, "RANGE");
    }
    if (! sinc_table_doppler.configSet())
    {
      sinc_table_doppler.config(cfg, "DOPPLER");
    }
  }

  // Setup beams
  bool read_beam_pattern;
  string beam_pattern_source = cfg.str("beam_pattern_source");
  if (beam_pattern_source == "file")
    {
    read_beam_pattern = true;
    }
  else
    {
    read_beam_pattern = false;
    }

  beam_.resize(5);
  for (unsigned int i=1; i <= 5; ++i)
    {
    beam_[i-1] = Beam(i,cfg,read_beam_pattern);
    Uvar elev=0;
    Uvar azim=beam_[i-1].maxGainLineToAzim(elev);
    beam_azioff_rad[i-1]=get_in_base_units(azim);
    }

  

  // read in instrument specific parameters --- need to move to SARProcParams
  updown_shift=-cfg["stalo_frequency"];
  filter_centered=(bool)cfg.getInt("ASSUME_FILTER_CENTERED");
  if(!filter_centered){
    filter_start_freq_offset=cfg["FILTER_START_FREQ_OFFSET"];
  }
  transmit_power=cfg["Pt"];


  // read in sar processing params -- need to move to SARProcParams
  range_ref_window_param=cfg.getDouble("RANGE_REF_WINDOW_PARAM");
  azimuth_ref_window_param=cfg.getDouble("AZIMUTH_REF_WINDOW_PARAM");
  DC_null_bins=cfg.getInt("DC_NULL_BINS");
  if(DC_null_bins < 0 || (DC_null_bins > 0 && DC_null_bins%2==0)) {
    ErrorMessage e("L1I::config Bad DC_null_Bins should be (0,1,3,5,...)");
    e.throwMe();
  }
  use_upper_band=(bool)cfg.getInt("USE_UPPER_BAND");
  use_linear_chirp=(bool)cfg.getInt("SARPROC_USE_LINEAR_CHIRP");
  perform_single_burst_check=(bool)cfg.getInt("PERFORM_SINGLE_BURST_CHECK");

  output_scatterometer_info=(bool)cfg.getInt("SARPROC_OUTPUT_SCATTEROMETER_INFO");

  pulse_spread_in_pri=cfg.getDouble("PULSE_SPREADING_IN_PRI");
  if(perform_single_burst_check){
    check_burst_idx=cfg.getInt("CHECK_BURST_NO");
    if(check_burst_idx < 1){
      ErrorMessage e("L1I:config bad Check Burst Number should be >= 1");
      e.throwMe();
    }
    check_burst_filename=cfg.str("CHECK_BURST_FILENAME");
    cbf.open(check_burst_filename.c_str());
  }



  // use stupid pulse segmentation
  if(cfg.keywordExists("USE_TRIVIAL_PULSE_SEGMENTATION")){
    use_trivial_pulse_segmentation_=(bool)cfg.getInt("USE_TRIVIAL_PULSE_SEGMENTATION");
    // Warn user
    if(use_trivial_pulse_segmentation_){
      cout << "Warning Using Trivial Pulse Sementation (For Debug Only ..)"
	   << endl;
      cout << "Warning Using Trivial Pulse Sementation (For Debug Only ..)"
	   << endl;
    }
  }
  else{
    use_trivial_pulse_segmentation_=false;
  }
  // set up intermediate data files if they exist
  if(cfg.keywordExists("L1I_EST_RANGE_DOPPLER_OFFSET_FILE")){
    string filename=cfg.str("L1I_EST_RANGE_DOPPLER_OFFSET_FILE");
    estrdfp=fopen(filename,"w");
    estrd_nskip=cfg.getInt("ESTRD_SKIP");
    estrd_dopwid=cfg.getInt("ESTRD_DOPPLER_WIDTH");
    estrd_ranwid=cfg.getInt("ESTRD_RANGE_WIDTH");
    estrd_iskip=0;
  }

  if(cfg.keywordExists("L1I_SIGMA0_FILE"))
  {
    string filename=cfg.str("L1I_SIGMA0_FILE");
    sigma0fp=fopen(filename,"w");    
  }

  if(cfg.keywordExists("L1I_INTERP_SIGMA0_FILE"))
  {
    string filename=cfg.str("L1I_INTERP_SIGMA0_FILE");
    is0fp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_XFACTOR_FILE"))
  {
    string filename=cfg.str("L1I_XFACTOR_FILE");
    Xfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_XAMBIG_FILE"))
  {
    string filename=cfg.str("L1I_XAMBIG_FILE");
    Xambfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_DECHIRPED_FILE"))
  {
    string filename=cfg.str("L1I_DECHIRPED_FILE");
    dcfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_RANGE_SPECTRUM_FILE"))
  {
    string filename=cfg.str("L1I_RANGE_SPECTRUM_FILE");
    rspecfp=fopen(filename,"w");    
  }


  if( cfg.keywordExists("L1I_RAW_SAR_FILE"))
  {
    string filename=cfg.str("L1I_RAW_SAR_FILE");
    rfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_MATCHED_FILTER_FILE"))
  {
    string filename=cfg.str("L1I_MATCHED_FILTER_FILE");
    mffp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_GAIN_SQUARED_FILE"))
  {
    string filename=cfg.str("L1I_GAIN_SQUARED_FILE");
    g2fp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_AREA_FILE"))
  {
    string filename=cfg.str("L1I_AREA_FILE");
    areafp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_RANGE_FILE"))
  {
    string filename=cfg.str("L1I_RANGE_FILE");
    rangefp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_DOPPLER_FILE"))
  {
    string filename=cfg.str("L1I_DOPPLER_FILE");
    fdopfp=fopen(filename,"w");    
  }
  
  if( cfg.keywordExists("L1I_AZIMUTH_FILE"))
  {
    string filename=cfg.str("L1I_AZIMUTH_FILE");
    azimfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_ELEVATION_FILE"))
  {
    string filename=cfg.str("L1I_ELEVATION_FILE");
    elevfp=fopen(filename,"w");    
  }

  if( cfg.keywordExists("L1I_PIXEL_DUMP_FILE") ){
    string filename=cfg.str("L1I_PIXEL_DUMP_FILE");
    pixdfp=fopen(filename,"w"); 
    pixd_noise_factor=cfg.getDouble("L1I_PIXEL_DUMP_NOISE_FACTOR");   
  }
}

//--------------
// I/O
//--------------


//------------------------------------------------------------------
// writeHeader()
//
// Write a L1I header to the current L1I file in binary.
// This method defines the L1I Header format.
// Takes a num_records parameter
// various overloads follow
//------------------------------------------------------------------

void L1I::writeHeader(int num_records) 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1IError e("Can't write to input file " + filename_, L1I::write_error);
      e.throwMe();
    }

  // Should not call this routine twice
  if (header_handled_)
    {
      L1IError e("Attempt to write header twice", L1I::write_error);
      e.throwMe();
    }
  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",num_records,record_size_);
  header=headbuf;
  if ( header.length() != header_size_){
    char msg[100];
    sprintf(msg,"Incorrect header size (%d bytes) should be %d",
	    header.length(),header_size_);
    L1IError e(msg,L1I::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  }

void L1I::writeHeader(){
  writeHeader(record_count_);
}

void L1I::writeHeader(BurstData& l1b){
  writeHeader(l1b.recordCount());
}

void L1I::rewriteHeader(){
int position = file_.getPosition();
  file_.rewind();
  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",record_count_,record_size_);
  header=headbuf;
  if ( header.length() != header_size_){
    char msg[100];
    sprintf(msg,"Incorrect header size (%d bytes) should be %d",
	    header.length(),header_size_);
    L1IError e("Header length is incorrect",L1I::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  file_.setPosition(position);
}

//------------------------------------------------------------------
// readHeader()
//
// Read a L1I header from the current L1I file in binary.
// This method defines the L1I Header format.
//------------------------------------------------------------------

void L1I::readHeader() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1IError e("Can't read from output file " + filename_,
	       L1I::read_error);
      e.throwMe();
    }
  if (header_handled_)
    {
      L1IError e("Attempt to read header twice",
	       L1I::read_error);
      e.throwMe();
    }
  checkEndedness();
  string header;
  header.resize(L1I_HEADER_LENGTH);
  file_.read(header);
  string token=get_next_token(header," /n/t"); // skip text
  token=get_next_token(header," /n/t"); // skip text
  
  token=get_next_token(header," /n/t"); // get num records
  record_count_=atoi(token.c_str());
  records_counted_=true;

  token=get_next_token(header," /n/t"); // skip text
  token=get_next_token(header," /n/t"); // skip text

  token=get_next_token(header," /n/t"); // check record size
  unsigned int recsize=(unsigned int) atoi(token.c_str());
  if(recsize!=record_size_){
    L1IError e("Record size incorrect in " + filename_, L1I::read_error);
    e.throwMe();
  }
  header_handled_ = true;
  }

//------------------------------------------------------------------
// writeRecord()
//
// Write the current L1I record to the end of the current L1I file in binary.
// This method defines the L1I record format.
// It is the responsibility of the user to ensure that valid data is present
// before calling this routine.  No tracking of data status is performed.
//------------------------------------------------------------------

void L1I::writeRecord() 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1IError e("Can't write to input file " + filename_, L1I::write_error);
      e.throwMe();
    }
  if (!header_handled_)
    {
      L1IError e("Can't write L1I record until header is written",
		 L1I::write_error);
      e.throwMe();
    }
  if (!data_complete)
    {
      L1IError e("Can't write L1I record until data is complete",
		 L1I::write_error);
      e.throwMe();
    }

  DebugInfo dbg("L1I::writeSARData");

  if(file_.getPosition() > MAX_FILE_LENGTH-(int)record_size_) 
    createNextFile();


  int start_position=file_.getPosition();
  writeSARData();
    
  int end_position=file_.getPosition();
  int len_written=end_position-start_position;
  if((unsigned int)len_written!=record_size_ && !dbg.level){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_written,record_size_);
    L1IError(msg,L1I::write_error).throwMe();
  }
  record_count_++;
  }


//------------------------------------------------------------------
// readRecord()
//
// Read a L1I record from the current L1I file.
//------------------------------------------------------------------

void L1I::readRecord() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1IError e("Can't read from output file " + filename_,
	       L1I::read_error);
//
// create next file
      e.throwMe();
    }  if (!header_handled_)
    {
      L1IError e("Can't read L1I record until header is read",
	       L1I::read_error);
      e.throwMe();
    } 

  // This check has the side effect of moving to the next  file if necessary
  if(eof())
    {
      L1IError e("Unexpected EOF in L1I file "+filename_,
	       L1I::read_error);
      e.throwMe();
    } 


 
  int start_position=file_.getPosition();
  
  readPassiveSABData();  
  readGeometry(); 
  readSARData();
  
  
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,record_size_);
    L1IError e(msg,L1I::read_error);
    e.throwMe();
  }
  data_complete=true;
  data_read_=true;  
  }




void L1I::resetRecord(){
  quality_flag=0;
  data_read_=false;
  data_complete=false;
}




// Parameter access methods

double L1I::getSpecialParam(const char* name){
  double retval;
  
  retval=BurstData::getSpecialParam(name);
  return(retval);
}
void L1I::enableSpecialParam(const char* name)
{
  BurstData::enableSpecialParam(name);
}
void L1I::disableSpecialParam(const char* name)
{
  BurstData::disableSpecialParam(name);
}

 

void L1I::readSARData()
{
  ErrorMessage e("L1I::ReadSARData is incomplete");
  e.throwMe();

  DebugInfo dbg("L1I::readSARData");
  if(dbg.level){
    file_.read(rms_radar_data);
    file_.read(num_range_bins);
    file_.read(num_pulses_received);

    unsigned int pos;
    for(unsigned int i=0;i<num_range_bins;i++){
      for(unsigned int j=0;j<num_pulses_received;j++){
	file_.read(raw_sar_data[i][j]);
      }
      pos=file_.getPosition();
      file_.setPosition(pos+(L1I_DOPPLER_WIDTH-num_pulses_received)*sizeof(float));
    }
    pos=file_.getPosition();
    unsigned int offset=L1I_DOPPLER_WIDTH*
      (L1I_RANGE_WIDTH-num_range_bins)*sizeof(float);
    file_.setPosition(pos+offset);
    
    file_.read(window_start_idx);
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      file_.read(matched_filter[i]);
    }
    
    for(unsigned int i=0;i<num_range_bins;i++){
      for(unsigned int j=0;j<num_pulses_received;j++){
	file_.read(dechirped[i][j]);
      }
      pos=file_.getPosition();
      file_.setPosition(pos+(L1I_DOPPLER_WIDTH-num_pulses_received)*2*sizeof(float));
    }
    pos=file_.getPosition();
    offset=L1I_DOPPLER_WIDTH*
     (L1I_RANGE_WIDTH-num_range_bins)*2*sizeof(float);
    file_.setPosition(pos+offset);
  }
}
  

void L1I::writeSARData()
{
  DebugInfo dbg("L1I::writeSARData");
  file_.write(num_range_bins);
  file_.write(Nfft_az);
  file_.write(num_pulses_received);  
  if(dbg.level){
    file_.write(rms_radar_data);
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      for(unsigned int j=0;j<L1I_DOPPLER_WIDTH;j++){
	file_.write(raw_sar_data[i][j]);
      }
    }
    file_.write(window_start_idx);
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      file_.write(matched_filter[i]);
    }
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      for(unsigned int j=0;j<L1I_DOPPLER_WIDTH;j++){
	file_.write(dechirped[i][j]);
      }
    }

  }
  for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
    for(unsigned int j=0;j<L1I_DOPPLER_WIDTH;j++){
      file_.write(sigma0[i][j]);
    }
  }
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      for(unsigned int j=0;j<L1I_DOPPLER_WIDTH;j++){
	file_.write(range[i]);
      }
    }
    
    for(unsigned int i=0;i<L1I_RANGE_WIDTH;i++){
      for(unsigned int j=0;j<L1I_DOPPLER_WIDTH;j++){
	double dop=absDoppler(i,j);
	file_.write(dop);
      }
    }

}

double L1I::absDoppler(int range_i, int doppler_j){
  double retval=fdop_centroid[range_i]+relDoppler(doppler_j);
  return(retval);
}

double L1I::relDoppler(int j){
  double bw=prf_in_Hz;
  double step=bw/Nfft_az;
  double dop=j*step-bw/2;
  return(dop);
}


void L1I::computeWindows()
{
  DebugInfo dbg("L1I::computeWindows");
  // for now arbitrarily set nominal spread(two-sided)  to X PRI 
  // (where X is set by config file variable PULSE_SPREADING_IN_PRI.
  // Actually the maximal useful overlap for segmentation should
  // be 1 pulse_width (two-sided) 
  nominal_pulse_spread=pulse_spread_in_pri*pri;

  // Commented out to allow windows to overlap
  // if(chirp_length+nominal_pulse_spread > pri){
  //  L1IError e("Windows Overlap.",internal_error);
  //  e.throwMe();
  //}

  // A special configurable debug mode in which pulse sgementation merely
  // splits the whole receive window into overlapping segments
  if(use_trivial_pulse_segmentation_){
    window_start_delay=0;
    num_pulses_received=ctrx;
    first_pulse_idx=0;
    window_start_idx=0;
  }
  else{

    // cout << " The real nominal delay is .."<< nominal_delay << endl;



    // OLD way range centroid estimate effects geo-location
    if(quick_look_dopran)
     window_start_delay=nominal_delay-nominal_pulse_spread/2-rwd + rc_time_shift;

    //NEW way only gain is effected by range centroid change
    else
      window_start_delay=nominal_delay-nominal_pulse_spread/2-rwd;
   

    
    // Omit missed and partial pulses at beginning
    num_pulses_received=pul;
    first_pulse_idx=0;
    while(window_start_delay<0 && num_pulses_received > 0){
      num_pulses_received--;
      window_start_delay+=pri;
      first_pulse_idx++;
    }
    
    if(num_pulses_received==0) return;

    
    window_start_idx=(unsigned int)((window_start_delay*adc).getValue()+0.5);  
  }

  // force even window_start_idx
  if(window_start_idx%2==1){
    window_start_idx++;
    window_start_delay+=1/adc;
  }
  // Parameters for sizing range compression arrays
  num_samples_per_window=(unsigned int)
    (((chirp_length+nominal_pulse_spread)*adc).getValue()+0.5);
  if(num_samples_per_window%2!=0) num_samples_per_window++;
  num_range_bins=num_samples_per_window/2;
  
  window_step=(unsigned int)((pri*adc).getValue()+0.5);
  if(window_step%2!=0 && DebugInfo::allWarnings)
    cerr << "Warning: Odd number of samples per PRI" << endl;

  // compute extent of central pri in window
  start_range_bin=(num_samples_per_window-window_step)/4;
  end_range_bin=start_range_bin+window_step/2-1;

  // Omit missed or partial pulses at the end
  unsigned int end_sample=window_start_idx+
    (num_pulses_received-1)*window_step
    +num_samples_per_window-1;


  if(window_start_idx + num_samples_per_window -1 > Nradar_data)
    {
      num_pulses_received=0;
    }
  else{
    while(end_sample>Nradar_data){
      num_pulses_received--;
      end_sample=window_start_idx+
	(num_pulses_received-1)*window_step
	+num_samples_per_window-1;
    }

  }


  Nfft_r = (int)get_next_power_of_2((unsigned int)num_samples_per_window);
  if(check_this_burst){
    cbf.set("pulse_spread",nominal_pulse_spread.getInUnits("s"));
    cbf.set("start_delay",window_start_delay.getInUnits("s"));
    cbf.set("start_idx",window_start_idx);
    cbf.set("window_len",num_samples_per_window);
    cbf.set("Nfft_r",Nfft_r);
    cbf.set("window_step",window_step);
    cbf.set("num_pulses",num_pulses_received);
    cbf.set("num_data",Nradar_data);
    cbf.set("end_idx",end_sample);
    cbf.set("pri",pri.getInUnits("s"));
  }

  if(dbg.level==1){
    dbg.file << "Start Debugging L1I::computeWindows" << endl;
    dbg.file << "pulse_spread " << nominal_pulse_spread.getInUnits("s") << endl;
    dbg.file << "start_delay " << window_start_delay.getInUnits("s") << endl;
    dbg.file << "start_idx " << window_start_idx << endl;
    dbg.file << "window_len " << num_samples_per_window << endl;
    dbg.file << "Nfft_r " << Nfft_r << endl;
    dbg.file << "window_step" << window_step << endl;
    dbg.file << "num_pulses" << num_pulses_received << endl;
    dbg.file << "num_data" << Nradar_data << endl;
    dbg.file << "end_idx" << end_sample << endl;
    dbg.file << "pri" << pri.getInUnits("s") << endl;
    dbg.file << "End Debugging L1I::computeWindows" << endl;
  }
 


  // Compute effective duty cycle of pulse-segmented 
  // uncompressed image

  num_signal_samples=0;
  int noise_samples=0;
  int one_side_spread_samples=
    (int)((adc*nominal_pulse_spread).getValue()/2+0.5);
  int one_pulse_samples=(int)((chirp_length*adc).getValue()+0.5);
  int one_gap_samples=window_step-one_pulse_samples;
  
  num_signal_samples=one_pulse_samples*num_pulses_received;
  noise_samples=one_gap_samples*(num_pulses_received)
    +2*one_side_spread_samples;
  effective_duty_cycle=(float)num_signal_samples/(float)(num_signal_samples+
						     noise_samples);

 
} 



/****
// weird HACKed routine for forcing the data to shift within useable area
void L1I::computeFakeWindows()
{

  int fake_npr=pul;
  int fake_fpi=0;
  DebugInfo dbg("L1I::computeWindows");
  // for now arbitrarily set nominal spread(two-sided)  to X PRI 
  // (where X is set by config file variable PULSE_SPREADING_IN_PRI.
  // Actually the maximal useful overlap for segmentation should
  // be 1 pulse_width (two-sided) 
  nominal_pulse_spread=pulse_spread_in_pri*pri;

  // Commented out to allow windows to overlap
  // if(chirp_length+nominal_pulse_spread > pri){
  //  L1IError e("Windows Overlap.",internal_error);
  //  e.throwMe();
  //}

  // A special configurable debug mode in which pulse sgementation merely
  // splits the whole receive window into overlapping segments
  if(use_trivial_pulse_segmentation_){
    window_start_delay=0;
    num_pulses_received=ctrx;
    first_pulse_idx=0;
    window_start_idx=0;
  }
  else{

    // cout << " The real nominal delay is .."<< nominal_delay << endl;

    // HACK ALERT
    // window_start_delay=nominal_delay-nominal_pulse_spread/2-rwd;

    window_start_delay=nominal_delay-nominal_pulse_spread/2-rwd;  
    
    // Omit missed and partial pulses at beginning
 
    while(window_start_delay<0 && fake_npr > 0){
      fake_npr--;
      window_start_delay+=pri;
      fake_fpi++;
    }
    
    if(fake_npr==0) return;

    
    window_start_idx=(unsigned int)((window_start_delay*adc).getValue()+0.5);  
  }

  // force even window_start_idx
  if(window_start_idx%2==1){
    window_start_idx++;
    window_start_delay+=1/adc;
  }
 
  
  // Omit missed or partial pulses at the end
  unsigned int end_sample=window_start_idx+
    (fake_npr-1)*window_step
    +num_samples_per_window-1;


  if(window_start_idx + num_samples_per_window -1 > Nradar_data)
    {
      fake_npr=0;
    }
  else{
    while(end_sample>Nradar_data){
      fake_npr--;
      end_sample=window_start_idx+
	(fake_npr-1)*window_step
	+num_samples_per_window-1;
    }

  }

} 
 
***/

void L1I::computeMatchedFilterHensley()
{

  //    initialize matched filter and zero pad to end
  for(int i=0;i<Nfft_r;i++){
    matched_filter[i]=complex<float>(0.0,0.0);
  }     
      
  //    compute the chirp bandwidth and center frequency
  double chirp_bandwidth = (csq)*slow_cfs_in_Hz;  
  double r_f0 = updown_shift_in_Hz+fast_csf_in_Hz+chirp_bandwidth/2;
  r_f0=r_f0+nominal_doppler_centroid.getInUnits("Hz");

  // chirp phase at middle of pulse
  double r_phase0 = 0;

  //------------------------------------------------------------------- 
  //     generate the time domain signal 
  //------------------------------------------------------------------- 
  int i_npts=0;   
  if(!use_linear_chirp){
    dig_chirp_hensley(chirp_length_in_s,chirp_bandwidth,csd_in_s,r_f0,
       r_phase0,adc_in_Hz,i_npts,matched_filter);
  }
  else{
    dig_chirp_hensley(chirp_length_in_s,chirp_bandwidth,
    chirp_length_in_s,r_f0,
    r_phase0,adc_in_Hz,i_npts,matched_filter);
  }
  


  //-----------------------------------------------------------------     
  //     compute scale of reference function 
  //-----------------------------------------------------------------     
      
  double fftgain=0;
  tmp_complex_array[0]=complex<float>(1,0);
  for(int i=1;i<Nfft_r;i++) tmp_complex_array[i]=complex<float>(0,0);

  fft(tmp_complex_array,tmp_complex_array2,Nfft_r);
  ifft(tmp_complex_array2,tmp_complex_array,Nfft_r/2);
  fftgain=real(tmp_complex_array[0]);
  // cout << "fftg1 " << fftgain << endl;

  //-------------------------------------------------------------
  //    Apply window to range reference function
  //-------------------------------------------------------------


  double a=range_ref_window_param;
  double b=1-a;   


  for (int i = 0;  i < i_npts; i++){
    double r_win = a - b * cos((2.* pi*i)/(i_npts-1));
    matched_filter[i] *= r_win;
  }


         
  //------------------------------------------------------
  // calculate fft of range reference function
  //-------------------------------------------------------       
  fft(matched_filter,fft_matched_filter,Nfft_r); 

  double sum2=0,sum4=0;
  for(int i=0;i<Nfft_r;i++){
    float dv=abs(fft_matched_filter[i]*fft_matched_filter[i]);
    sum4+=dv*dv;
    sum2+=dv;
  }  

  double scale=sqrt(sum2)/(fftgain*sqrt(sum4));
 
  // scale reference function
  for( int i=0; i< i_npts; i++ ){
    matched_filter[i]*=scale;
  }     
  for( int i=0; i< Nfft_r; i++ ){
    fft_matched_filter[i]*=scale;
  }     
  
  //----------------------------------------------------    
  //    zero out DC from spectrum
  //----------------------------------------------------
  
  if(DC_null_bins != 0){       
    for(int i = -DC_null_bins/2; i<= DC_null_bins/2; i++){
      double wgt = 0.5 - 0.5 * cos((float)i/(DC_null_bins/2)*pi);
      int ind = (i+Nfft_r)%Nfft_r;
      fft_matched_filter[ind] *= wgt;
    }
  
    //-------------------------------------------------------
    //    Zero out mid bandwidth DC term in offset video data
    //-------------------------------------------------------

    for(int i = -DC_null_bins/2; i<= DC_null_bins/2; i++){
      double wgt = 0.5 - 0.5 * cos((float)i/(DC_null_bins/2)*pi);
      int ind = i+Nfft_r/2;
      fft_matched_filter[ind] *= wgt;
    }
  } // end if DC_null_bins nonzero

}

void
L1I::computeDoubleParams()
{
 rc_bw_in_Hz=rc_bw.getInUnits("Hz");
 adc_in_Hz=adc.getInUnits("Hz");
 chirp_length_in_s=chirp_length.getInUnits("s");
 rwd_in_s=rwd.getInUnits("s");
 fast_csf_in_Hz=fast_csf.getInUnits("Hz");
 updown_shift_in_Hz=updown_shift.getInUnits("Hz");
 pri_in_s=pri.getInUnits("s");
 prf_in_Hz=1/pri_in_s;
 slow_cfs_in_Hz=slow_cfs.getInUnits("Hz");
 csd_in_s=csd.getInUnits("s");
 lambda_in_km=cassini_sim.lambda().km();
 lambda_chirp_in_km=lambda_chirp.km();
 ground_impact_time.getEt(ground_impact_time_in_s);
 nominal_doppler_centroid_in_Hz=nominal_doppler_centroid.getInUnits("Hz");

  Uvar nr=cassini_sim.altitude();
  nadir_range_in_km=nr.getInUnits("km");
  Uvar range_window=pri*speed_light/2;
  range_window_in_km=range_window.getInUnits("km");
  Uvar ndop=cassini_sim.nadirDoppler();
  nadir_fdop_in_Hz=ndop.getInUnits("Hz");
  nadir_fdop_in_Hz*=lambda_chirp_in_km/lambda_in_km;
}

void L1I::computeDopplerAndRange(const SARProcParams& spp)
{
  DebugInfo dbg("L1I::computeDopplerAndRange");
  if(dbg.level){
    dbg.file << "Start Debug L1I::computeDopplerAndRange ..." << endl;
    dbg.file << "Beam AZI OFF in degrees="
	     <<beam_azioff_rad[beam_number-1]*180/pi << endl;
    dbg.file << "Range [look] fdop_centroid fdop_rate [look_from_dop]" << endl;
    dbg.file.precision(20);
  } 
  // Constructs intermediate Uvars once
  static Uvar td,r,nom_range,fdc,fdcrate;
  int mult=2; // 2 real samples per range bin

  // get doppler and range steps as doubles
  doppler_step_=prf_in_Hz/Nfft_az;

  // speed of light * half of the time between complex samples
  range_step_=get_in_base_units(speed_light)/adc_in_Hz;  

  if(spp.azimuth_range_orthogonalization_on) cerr << "Warning AZIMUTH_RANGE_ORTOGONALIZATION is ON!!!" << endl;
  

  // HACK range_boresight holds the delta range value
  nom_range=nominal_delay*speed_light/2;
   
  // For each equi-range line compute doppler centroid and rate
  for(unsigned int k=0;k<num_range_bins;k++){
    // compute range

    // OLD WAY range centroid estimate effects geo-location
    if(quick_look_dopran)
       td=mult*k/adc+window_start_idx/adc+rwd-pri*first_pulse_idx-rc_time_shift;    // New way only gain is effected by range centroid estimate
    else
      td=mult*k/adc+window_start_idx/adc+rwd-pri*first_pulse_idx;
    r=td*speed_light/2;
    range[k]=r.km();
    
    // compute nominal doppler and range indices first time through
    if(k==0){
      nominal_range_idx=int(get_in_base_units((nominal_delay-td)
					      /(mult/adc))   +0.5);

      if(nominal_range_idx > (signed int) num_range_bins){
	if(!use_trivial_pulse_segmentation_){
	  ErrorMessage e("Error Nominal range index is outside of range array");
	  e.throwMe();
	}
	else{
	  nominal_range_idx=num_range_bins/2;
	}
      }
    }


    // compute doppler centroid vs. range

    double look[3];
    int look_good; // return value from getFasTLookVector
                   // if it is 0 then range is off planet use 
                   // nominal_doppler centroid and 0 rate
    
    if(!spp.linear_dopran_approx){
      // WARNING! WARNING! WARNING!
      // Using a constant azim offset
      // fails when elev direction and range direction are
      // much different, so far this only occurs with high SAR near the pole
      // In this case it is best to set azimuth_range_orthogonalization_on=0
      look_good=getFastLookVector(range[k],beam_azioff_rad[beam_number-1],look);
    }
    if(spp.azimuth_range_orthogonalization_on){
      // special debugging case in which doppler centroids are specified
      // as polynomials in Config file 

      if(spp.linear_dopran_approx){
	fdc=nominal_doppler_centroid+ddopdrange*(r-nom_range);
	fdop_centroid[k]=get_in_base_units(fdc);
      }
      else if(!look_good){
	fdop_centroid[k]=get_in_base_units(nominal_doppler_centroid);
      }
      else{
	fdop_centroid[k]=getFastDopplerCentroid(look);
      }
    }
    else{
      fdop_centroid[k]=get_in_base_units(nominal_doppler_centroid);
    }
    // compute doppler rate
    if(spp.azimuth_deramp_processing_on){
      if(spp.set_doppler_centroid || spp.linear_dopran_approx){
	fdcrate=nominal_doppler_rate+ddopratedrange*(r-nom_range);
	fdop_rate[k]=get_in_base_units(fdcrate);
      }
      else if(!look_good){
	fdop_rate[k]=0;
      }
      else{
	fdop_rate[k]=getFastDopplerRate(look,range[k],fdop_centroid[k]);
      }
    }
    else{
      fdop_rate[k]=0;
    }
    if(dbg.level==1){
      if(k==0){
	dbg.file <<"range[0]=" << range[0] 
		 <<" range_step=" << range_step_
		 <<" doppler_step=" << doppler_step_ << endl;
      }

      dbg.file << "Range step " << k <<": fdop[0]=" 
<< fdop_centroid[k]-prf_in_Hz/2
	       << endl;
    }
    if(dbg.level==2){
      double look2[3],otherlook[3];
      fast_doppler_range_to_TBF_look(fdop_centroid[k],range[k],target_radius_,scpos_,scvel_,lambda_special,look2,otherlook);

      double cosangle1=dot(boresight_in_tbf_,look2);
      double cosangle2=dot(boresight_in_tbf_,otherlook);
      for(int c=0;c<3;c++){
	if(cosangle1<cosangle2) look2[c]=otherlook[c];
      }
      dbg.file << range[k] <<"  [" << look[0] << " " << look[1] << " "
	       << look[2] << "] " << fdop_centroid[k] << " "
	       << fdop_rate[k] << " [" << look2[0] << " " << look2[1]
	       << " " << look2[2] << "]" << endl;
    }
    // compute range component of Xfactor;
    for(int j=0;j<Nfft_az;j++){
      Xfactor[k][j]=1.0/(range[k]*range[k]*range[k]*range[k]);
    }

  }

  // fix range error due to doppler variation with range
  // (right now this creates an irregular grid which we will fix later)
  /**
  double chirp_samps=floor(chirp_length_in_s*adc_in_Hz+0.5);
  double chirp_bw=csq*slow_cfs_in_Hz;
  double cr=chirp_bw/(chirp_samps/adc_in_Hz);
    // range[k] err = -range_step*(fdop_centroid[k]-nominal_dop_centroid)
    //                / change in chirp over one sample period (complex) 
  for(unsigned int k=0;k<num_range_bins;k++){
    cout << chirp_samps << " " << chirp_bw << " " << cr << " " << range[k]
	 << " " << range_step_ << " ";
    double ddop= fdop_centroid[k]-get_in_base_units(nominal_doppler_centroid);
    double dchirp = 2*cr/adc_in_Hz;
     range[k]+= range_step_*ddop/dchirp;
     cout << range[k] << " " << ddop << " " << dchirp << endl;
  }
  **/
  if(dbg.level){
   dbg.file << "End Debug L1I::computeDopplerAndRange ..." << endl;
  }

}

double
 L1I::getFastDopplerCentroid(double look[3]){


  // dop=2/lamb * dot(v,l)
  // classical doppler
  // double dop=scvel_[0]*look[0]+scvel_[1]*look[1]+scvel_[2]*look[2];
  // dop*=2/lambda_chirp_in_km;


  // Doppler due to special relativity
  double vr = (scvel_[0]*look[0]+scvel_[1]*look[1]+scvel_[2]*look[2]);
  double dop=2*vr*chirp_freq_in_Hz/(speed_light_in_kps-vr);
  
  //#define DEBUG
#ifdef DEBUG
  cout.precision(20);
  cout << "  lambda_chirp=" << lambda_chirp_in_km
                     << ", look_x = " << look[0]
                     << ", look_y = " << look[1]
                     << ", look_z = " << look[2]
		     << ", vel_x = " << scvel_[0]<< " km/s,"
		     << "vel_y = " << scvel_[1] << " km/s,"
		     << "vel_z = " << scvel_[2]<< " km/s,"
		  << ", fdop=" << dop << endl;
#endif
  return(dop);
}

double
L1I:: getFastDopplerRate(double look[3],double range_in_km, double fdop){
  //dfdt=2/lamb  * [ dot(a,l)-|v|^2/range] + f^2*lambda/2range
  double adotl=scacc_[0]*look[0]+scacc_[1]*look[1]+scacc_[2]*look[2];
  double dfdt=(2/lambda_chirp_in_km)*
    (adotl-scspeed_*scspeed_/range_in_km) +
    (lambda_chirp_in_km/2/range_in_km)*fdop*fdop;
  return(dfdt);
}


// returns look vector in beam frame given range and azimuth
int
L1I::getFastLookVector(double range_in_km,float azim_rad,double look[3]){

  
  // --- Nadir pointing frame in titan body coordinates ---
  double px = rotmat_[2][0];
  double py = rotmat_[2][1];
  double pz = rotmat_[2][2];
  
  // determine unit vector toward nadir in target body fixed coordinates
  // This code makes a Spherical body assumption
  double P0[3];
  double r0=sqrt(scpos_[0]*scpos_[0]+scpos_[1]*scpos_[1]+scpos_[2]*scpos_[2]);
  P0[0]=-scpos_[0]/r0;
  P0[1]=-scpos_[1]/r0;
  P0[2]=-scpos_[2]/r0;
  
  //--------------------------
  // define beam direction
  //--------------------------
  // compute y-component of P0 in beam frame
  double pyb=P0[0]*rotmat_[1][0]+P0[1]*rotmat_[1][1]+P0[2]*rotmat_[1][2];
  int beam_direction = 1;
  if(pyb>0) beam_direction = -1;

  //cout<<"beam direction "<<beam_direction<<endl;
 
  double cx = P0[0];
  double cy = P0[1];
  double cz = P0[2];
  
  double b0 = sqrt(pow(cy*pz-py*cz,2.)+pow(cz*px-pz*cx,2.)+pow(cx*py-px*cy,2.));
  double bx = (cy*pz-py*cz)/b0;
  double by = (cz*px-pz*cx)/b0;
  double bz = (cx*py-px*cy)/b0;
  
  double ax = by*cz-cy*bz;
  double ay = bz*cx-cz*bx;
  double az = bx*cy-cx*by;
  
  // Matrix product rotmat_ * Trans(T_T2r)
  double m11 = rotmat_[0][0]*ax +  rotmat_[0][1]*ay +  rotmat_[0][2]*az;
  double m12 = rotmat_[0][0]*bx +  rotmat_[0][1]*by +  rotmat_[0][2]*bz;
  double m13 = rotmat_[0][0]*cx +  rotmat_[0][1]*cy +  rotmat_[0][2]*cz;
  double m21 = rotmat_[1][0]*ax +  rotmat_[1][1]*ay +  rotmat_[1][2]*az;
  double m22 = rotmat_[1][0]*bx +  rotmat_[1][1]*by +  rotmat_[1][2]*bz;
  double m23 = rotmat_[1][0]*cx +  rotmat_[1][1]*cy +  rotmat_[1][2]*cz;
  double m31 = rotmat_[2][0]*ax +  rotmat_[2][1]*ay +  rotmat_[2][2]*az;
  double m32 = rotmat_[2][0]*bx +  rotmat_[2][1]*by +  rotmat_[2][2]*bz;
  double m33 = rotmat_[2][0]*cx +  rotmat_[2][1]*cy +  rotmat_[2][2]*cz;  
    
 
  if ((range_in_km<=r0-target_radius_)||
      (range_in_km>=sqrt(r0*r0-target_radius_*target_radius_)))
    {
      if(DebugInfo::allWarnings){
	cerr<< "Warning L1I::getFastLookVector: Range is either too large or too small" << endl;
      }
      return(0);
    }
    
  
  double xi = (r0*r0 + range_in_km*range_in_km - target_radius_*target_radius_) / (2.*range_in_km*r0);
  //cout << "   xi = "<<xi<<endl;
  
  //----- Calculate elevation ------
  double azi0=beam_azioff_rad[beam_number-1];
  double a = (m11-m31*tan(azi0))*sqrt(1-xi*xi);
  double b = (m12-m32*tan(azi0))*sqrt(1-xi*xi);
  double c = (m13-m33*tan(azi0))*xi;
  //cout << " a= " << a << " b=" <<b << " c=" << c << endl;
  
  double sol = a*a + b*b -c*c;
  if (sol >= 0)
    {
      double sin_phi = (double(beam_direction)*a*sqrt(sol) - b*c)/(a*a + b*b);
      double cos_phi = -(b*sin_phi+c)/a;
      //cout<<"sin_phi: "<<sin_phi<<" cos_phi: "<<cos_phi<<endl;
      double sin_ele0 = (m21*sqrt(1-xi*xi)*cos_phi + m22*sqrt(1-xi*xi)*sin_phi + m23*xi);
      double ele0=asin(sin_ele0);
      // compute look vector in beam frame
      double look_bf[3];
      look_bf[0]=sin(azi0)*cos(ele0);
      look_bf[1]=sin(ele0);
      look_bf[2]=cos(azi0)*cos(ele0);
      rotate_vector(look_bf,rotmat_rev_,look);
      
    }
  else{
    return(0);
  }
  return(1);
}

void L1I::computeBurstParams(const SARProcParams& spp, ofstream& afs, RasList* raslist)
{ 
  DebugInfo dbg("L1I::computeBurstParams");
  if(perform_single_burst_check){
    if (burst_no==check_burst_idx){
      check_this_burst=true;
    }
    else check_this_burst=false;
  }


  // add transmit offset into time t so that all geometry computations
  // are made at the appropriate time
  t+=transmit_time_offset;

  // need to modify this so fiter_start_freq_offset does not depend on
  // anti-aliasing filter bandwidth (or at least Config file parameter does
  // NOT
  if(filter_centered){
    filter_start_freq=-adc/2+(adc/2-rc_bw)/2;
  }
  else{
    filter_start_freq=-adc/2+filter_start_freq_offset;
  }



  // Oct 24 2010 - Bryan Stiles
  // If a rerouted chirp burst is being processed,
  // fix it to be processable. This should only occur if
  // we have marked that burst as "GOOD" in the quality override file.
  if(csr==4) {
    fixRerouteChirpBurst();
  }

  //--- Update cassini_sim -----//
  // use time halfway between transmit window center and receive window
  // center to get nominal_doppler_centroid and nominal_delay
  Uvar receive_window=ctrx*pri;
  Uvar transmit_window=pul*pri;
  Uvar transmit_center=transmit_window/2;
  Uvar receive_center=rwd+receive_window/2;
  Uvar time_offset=(transmit_center+receive_center)/2;

  // Choice of ground impact time has a very small impact on range
  // Even if it is off several pri ( 1 ms) if will be at most 3 m off.
  ground_impact_time=t+time_offset;
  


  
  cassini_sim.setBeamTime(beam_number,ground_impact_time);



  // Estimate the wavelength lambda_chirp at midpoint of each pulse for
  // use in accurately estimating doppler  

  Uvar chirp_freq=speed_light/cassini_sim.lambda(); // carrier
  chirp_freq+=fast_csf; // add chirp start freq
  chirp_freq+=slow_cfs*(csq-1)/2; // add half of chirp bandwidth
  lambda_chirp=speed_light/chirp_freq; // convert to wavelength
  
  // Compute intermediate arrays of geometry data for use in later calculations
  initGeometryCalculations();





  // Special Debugging Option
  // set doppler centroid values and range values from config file
 
  rc_time_shift=0;
  ac_dop_shift=0;
  rc_range_shift_in_km=0;
  ac_dop_shift_in_Hz=0;
  if(spp.set_doppler_centroid){
    
    Uvar dnomrange;
    Uvar ddopc;
    spp.setDopplerProfile(ground_impact_time, beam_number, ddopdrange, 
			  ddopc, dnomrange,
			  ddopratedrange,nominal_doppler_rate,sab_counter);

    // shifts in range and azimuth compression to avoid pointing/ephem error
    rc_time_shift=2*dnomrange/speed_light;
    // cout << "rc_time_shift " << rc_time_shift << endl;
    // cout << "pri_in_s " << pri_in_s << endl;
    
    rc_time_shift=Uvar(fmod(get_in_base_units(rc_time_shift),get_in_base_units(pri)),"s");
    // cout << "fmod(rc_time_shift,pri_in_s)" << rc_time_shift << endl;
    rc_range_shift_in_km=get_in_base_units(rc_time_shift*speed_light/2);
    ac_dop_shift=ddopc;
    ac_dop_shift_in_Hz=get_in_base_units(ddopc);
  }

  //  Old linear doppler vs. range doppler profile estimate
  //  only used in special debug case
  else if(spp.linear_dopran_approx){
    cassini_sim.getDopplerProfile(lambda_chirp,
			
	  ddopdrange, nominal_doppler_centroid,
				  nominal_range, ddopratedrange, 
				  nominal_doppler_rate);
  }

  
  // New version is in computeDopplerAndRange

  // set nominal delay
  nominal_delay=nominal_range*2/speed_light;


  if(check_this_burst){
    cbf.set("nom_del",nominal_delay.getInUnits("s"),burst_no);
    cbf.set("nom_dop",nominal_doppler_centroid.getInUnits("Hz"),burst_no);
  }

 

  // convert Uvar and Time parameters to doubles as needed
  // All such parameters need to be have already been computed at this point
  computeDoubleParams();



  //compute pulse segmentation parameters
  //Nominal range spreading is computed using a particular limiting antenna
  //gain (Should be output of SAR preproc Gl(time to closest approach).
  

  // cout << "nominal_delay " << nominal_delay << endl; 
  computeWindows();

  // Special noise only processing case
  // THIS GETS RID OF THE SIGNAL !!!!!
  if(spp.signal_type!=SARProcParams::NOMINAL){
    double est_noise_var=0.0;
    double min_var=0;
    double minmean=0;
    // compute est_noise_variance from data if desired
    if(spp.signal_type==SARProcParams::DATA_NOISE){
      // compute noises variance from PRI at edge of window
      // chooses side of window which is most removed from signal
      
      int i1;
      

      int end_idx=window_start_idx+(num_pulses_received-1)*window_step
	+num_samples_per_window-1;

      int dnns=(int)(pri_in_s*adc_in_Hz+0.5);

      // closer to start case
      if(first_pulse_idx>0 || (Nradar_data-end_idx)>window_start_idx) 
	i1=Nradar_data-dnns;
      // closer to end case
      else i1=0;
      

      float dnmean=0;      
      for(int c=i1;c<i1+dnns;c++){
	//fprintf(spp.data_noise_out_fp,"%d %g\n", c,radar_data[c]);
	est_noise_var+=radar_data[c]*radar_data[c];
	dnmean+=radar_data[c];
      }
      //est_noise_var-=dnmean*dnmean/dnns;
      //fprintf(spp.data_noise_out_fp,"&");

      est_noise_var/=dnns;
      dnmean/=dnns;

      if(spp.simBaqForNoise){
	for(int c=i1;c<i1+dnns;c++){
	  radar_data[c]=0;
	}
	BAQEncodeAndDecode(raslist);

	for(int c=i1;c<i1+dnns;c++){

	  //fprintf(spp.data_noise_out_fp,"%d %g\n", c,radar_data[c]);
	  min_var+=radar_data[c]*radar_data[c];
	  minmean+=radar_data[c];
	}

	//fprintf(spp.data_noise_out_fp,"&");
	min_var/=dnns;
	est_noise_var-=min_var;
      }
      Uvar est_tsys=spp.NoiseVarianceToSystemTemperature(est_noise_var,*this);
      if(spp.data_noise_out_fp){
	fprintf(spp.data_noise_out_fp,"%d %d %g %g %g %g %g\n",(int)sab_counter,
		(int) beam_number,
		get_in_base_units(rc_bw),est_noise_var,
		get_in_base_units(est_tsys), min_var, dnmean);
        fflush(spp.data_noise_out_fp);
      }
    }
    spp.replaceSignalWithNoise(r_mode,*this,est_noise_var);
    if(spp.simBaqForNoise && spp.signal_type==SARProcParams::THEOR_NOISE){
      
      BAQEncodeAndDecode(raslist);  
    }
  }
  if(dbg.level>0){
    dbg.file << "L1I::computeBurstParams is complete." << endl;
    dbg.file << "Nominal Doppler Centroid=" << nominal_doppler_centroid;
    dbg.file << ",Nominal Delay=" << nominal_delay << endl;    
    dbg.file << "Nominal Doppler Rate=" << nominal_doppler_rate;
    dbg.file << ", Doppler Range Slope=" << ddopdrange;
    dbg.file << ", Doppler Rate Range Slope=" << ddopratedrange << endl;
    dbg.file << "Num_range_bins=" << num_range_bins;
    dbg.file << ", Num_pulses_received=" << num_pulses_received << endl;
  }
  
  // return if no pulses found
  if(num_pulses_received==0) return;

  // checks for sizing error
  if(num_samples_per_window>L1I_RANGE_WIDTH ||
     num_pulses_received > L1I_DOPPLER_WIDTH){
    L1IError e("L1I bad sizes DopWidth=" + toStr(num_pulses_received) +
	       + "Range Width="  
	       + toStr(num_samples_per_window),
	       internal_error);
    e.throwMe();
  }  

  //get size for azimuth FFTs 
  //size of range fft computed in computeWindows
  Nfft_az=(int)get_next_power_of_2((unsigned int) num_pulses_received);
  if(Nfft_az > L1I_DOPPLER_WIDTH){
    L1IError e("L1I bad size DopWidth=" + toStr(Nfft_az),
	       internal_error);
    e.throwMe();
  }

  // compute Doppler, dDopplerdTime  and Range for each range bin
  computeDopplerAndRange(spp); 



  // compute effective chirp rate for use in correcting range error 
  // due to doppler
  double chirp_samps=floor(chirp_length_in_s*adc_in_Hz+0.5);
  double chirp_bw=csq*slow_cfs_in_Hz;
  chirp_rate_in_Hz_per_s=chirp_bw/(chirp_samps/adc_in_Hz);
  
  // compute correction coefficient to account for error in range compression
  // due to doppler
  // range(true)=range(SAR)+RDCC*(fdop-nominal_doppler_centroid)
  range_doppler_correction_coeff=range_step_
    /(2*chirp_rate_in_Hz_per_s/adc_in_Hz);

  // initialize RangeDoppler to AzimuthElevation transformation in beam frame
  initAzimElevTransformation();  



  afs << " Computed Burst Parameters ... " << endl
      << " Nominal Boresight ( Range = "
      << get_in_base_units(nominal_delay*speed_light)/2 << " km, "
    "Doppler = " <<  get_in_base_units(nominal_doppler_centroid)/1000
      << "kHz) " << endl
      << " " << num_samples_per_window*num_pulses_received 
      << " samples processed in receive window." << endl
      << " effective duty cycle for pulse segmented window is:"
      <<  effective_duty_cycle*100 <<"%" << endl;
}

float L1I::estimateBaqInputRatio()
{
  double  prebaq_absmean=0.0;
  double  prebaq_absmean_ends=0.0;
  unsigned int  offset=(unsigned int)(get_in_base_units(pri*adc*8)+0.5);
  for(unsigned int c=0;c<Nradar_data;c++){
    if(c<offset || c>Nradar_data-offset-1){
      prebaq_absmean_ends+=fabs(radar_data[c]);
    }
    else{
      prebaq_absmean+=fabs(radar_data[c]);
    }
  }
  prebaq_absmean/=(Nradar_data-offset*2);
  prebaq_absmean_ends/=offset*2;
  float baq_inp_rat=(float)prebaq_absmean_ends/prebaq_absmean;
  return(baq_inp_rat);
}
float L1I::numBlankPulsesForBAQ(){
        Uvar tro_pri = tro/pri;
	int tro_int=int(floor(tro_pri.getInUnits("")+0.5));
	float num_first_blank=get_in_base_units((nominal_delay-rwd)/pri);
        if(num_first_blank<0) num_first_blank=0;
	if(num_first_blank>8) num_first_blank=8;
        Uvar end_rec_window=rwd+(tro_int+pul)*pri;
        Uvar end_last_pulse=nominal_delay+pul*pri;
	float num_last_blank=get_in_base_units((end_rec_window-end_last_pulse)/pri);
        if(num_last_blank<0) num_last_blank=0;
	if(num_last_blank>8) num_last_blank=8;

	float num_blank_pris=num_first_blank+num_last_blank;
        return(num_blank_pris);
}
 
float L1I::estimateBaqOutputRatio(RasList* raslist){
  int use_hfb=0;
  float y;
  if(use_hfb){
    Uvar tro_pri = tro/pri;
    int tro_int=int(floor(tro_pri.getInUnits("")+0.5));
    baq.setParams(pri,pul,baq_mode,adc,tro_int);
    Ivec Thresh("Thresh");
    // Get Thresholds from rasfile
    raslist->findRecord(sab_counter);
    raslist->decodeBaqThresholds();
    Thresh.resize(raslist->baq_decoded.size());
    for(unsigned int c=0;c<raslist->baq_decoded.size();c++)
      Thresh(c)=raslist->baq_decoded(c);
    y=estimateBaqOutputRatio(Thresh);
  }
  else y=estimateBaqOutputRatio();
  return(y);

}


float L1I::estimateBaqOutputRatio(Ivec& Thresh){
 

  float frac=fracHighBit(Thresh);
  float numblank=numBlankPulsesForBAQ();
  float y= FHBToBAQBias(frac,numblank);
  return(y);
}


float L1I::estimateBaqOutputRatio(const L1B& l1b){
  // returns 1 for 8 to 8
  if(l1b.baq_mode==5) return(1);

  pri=l1b.pri;
  adc=l1b.adc;
  Nradar_data=l1b.Nradar_data;
  baq_mode=l1b.baq_mode;
  for(unsigned int c=0;c<Nradar_data;c++) radar_data[c]=l1b.radar_data(c);
  float y=estimateBaqOutputRatio();
  return(y);
}

float L1I::estimateBaqOutputRatio(){
 


  float x=estimateBaqInputRatio();
  if(x<0.77) x=0.77;
  if(x>1.026) x=1.026;
  float y=-2.3936+5.7853*x-2.43*x*x;

  return(y);
}


float L1I::FHBToBAQBias(float fhb, float numblank){
  static float a0[8]={1,1,1,1,1,1,1,1};
  static float a1[8]={0,0,0,0,0,0,0,0};
  static float a2[8]={0,0,0,0,0,0,0,0};
  static float a3[8]={0,0,0,0,0,0,0,0};

  float maxfhb = 0.5;
  float minfhb = 0.3;
  if(fhb>maxfhb) fhb=maxfhb;
  if(fhb<minfhb) fhb=minfhb;
  int i=(int)(numblank-0.5);
  if(i>7) i=7;
  if(i<0) i=0;
  float y=a0[i]+a1[i]*fhb+a2[i]*fhb*fhb+a3[i]*fhb*fhb*fhb;
  return(y);
}

float L1I::fracHighBit(Ivec& Thresh){
  // Right Now this only works for 8 to 2 BAQ mode
  if(Thresh.size()!=24) return(0);
  int nhighbits=0;
  int NSamplesPerBlock=(int)(get_in_base_units(pri*adc)/24 +0.5);
  int NSamplesPerPri=(int)(get_in_base_units(pri*adc)+0.5);
  for(unsigned int c=0;c<Nradar_data;c++){
    int pulseno=c/NSamplesPerPri;
    int blockno=(c-pulseno*NSamplesPerPri)/NSamplesPerBlock;
    if(blockno==24)blockno=23;
    if(fabs(radar_data[c])> Thresh(blockno)/2) nhighbits++;
  }
  return((float)nhighbits/(float)Nradar_data);
}
void L1I::rangeCompressHensley(const SARProcParams& spp, ofstream& afs)
{

  
  // cout << "Starting index " << window_start_idx << endl;
  // cout << "Start delay " << window_start_delay << endl;

  DebugInfo dbg("L1I::rangeCompressHensley");
  if(check_this_burst) 
    cbf.comment("Using Hensley's range compression ...");

  if(check_this_burst){
      cbf.comment("Number of Samples in each window");
      cbf.set("num_samples_per_window",num_samples_per_window);
      cbf.comment("Number of Pulses Received");
      cbf.set("num_pulses_received",num_pulses_received);
  }

  // compute reference function (matched filter)
  computeMatchedFilterHensley();

  // compute mean square of input values for use in normalization
  mean_square_inputs=0.0;
  float sum_square_overlapped_inputs=0;
  int num_samples_used=
    (num_pulses_received-1)*window_step+num_samples_per_window;
  sum_square_dechirp=0.0;
  sum_square_fulldechirp=0.0;


  for(int i=0;i<num_samples_used;i++){
    float value=radar_data[window_start_idx+i];
    mean_square_inputs+=value*value;
  }
  // for each pulse perform range compression
  for(unsigned int i=0;i<num_pulses_received;i++){

    int pulse_offset=(int)(window_start_idx+i*window_step);
    // copy and zero-pad signal 
    for(int j=0;j<(int)num_samples_per_window;j++){
      int sample_offset=pulse_offset+j;
      float real_value=radar_data[sample_offset];
      tmp_complex_array[j]=complex<float>(real_value,0.0);
      sum_square_overlapped_inputs+=real_value*real_value;
    }
    for(int j=(int)num_samples_per_window;j<Nfft_r;j++){
      tmp_complex_array[j]=complex<float>(0,0);
    }
    float input_energy=0.0;
    for(int j=0;j<Nfft_r;j++) 
      input_energy+=abs(tmp_complex_array[j])*abs(tmp_complex_array[j]);
    if(check_this_burst){
      cbf.comment("Offset into echo array");
      cbf.set("pulse_offset",pulse_offset,i+1);
      cbf.comment("");
      cbf.comment("Zero-pad echo data in window");
      cbf.set("hrcdata1",tmp_complex_array,Nfft_r,i+1,-1);
    }

    
    // perform FFT on signal
    fft(tmp_complex_array,tmp_complex_array2,Nfft_r);

    if(check_this_burst){
      cbf.comment("");
      cbf.comment("Perform FFT on echo data");
      cbf.set("hrcdata2",tmp_complex_array2,Nfft_r,i+1,-1);
    }

    if(!use_upper_band){
      //----------------------------------------------------------
      // Use Lower sideband
      //----------------------------------------------------------

      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Using lower sideband ...");
      }

      // multiply frequency domain signal and conjugate of reference function
      for(int j=0; j<Nfft_r/2; j++){ 
	tmp_complex_array2[j]*=conj(fft_matched_filter[j]);
      }

      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Multiply in frequency domain by conj(ref_func)");
	cbf.set("hrcdata3",tmp_complex_array2,Nfft_r,i+1,-1);
        cbf.set("mf_fft",fft_matched_filter,Nfft_r/2,i+1,-1);
      }      

      // baseband (rotate spectra)
      for(int j=0;j<Nfft_r/4;j++) 
	tmp_complex_array2[j+Nfft_r/2]=tmp_complex_array2[j];
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Rotate Spectra");
	cbf.set("hrcdata4",tmp_complex_array2,Nfft_r,i+1,-1);
      }      

      // inverse FFT 
      ifft(&tmp_complex_array2[Nfft_r/4],tmp_complex_array,Nfft_r/2);
     
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Inverse FFT");
	cbf.set("hrcdata5",tmp_complex_array,Nfft_r/2,i+1,-1);
      }      
      
      // copy to storage array
      for(int j=0;j<Nfft_r/2;j++){
	dechirped[j][i]=tmp_complex_array[j]; 
	float val=abs(tmp_complex_array[j]);
        val*=val;
	sum_square_fulldechirp+=val;
        if(j<(int)num_range_bins) sum_square_dechirp+=val;
      }
    }
    else{
      //------------------------------------------------------------------
      // Use Upper sideband
      //------------------------------------------------------------------
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Using lower sideband ...");
      }

      // multiply frequency domain signal and conjugate of reference function
      for(int j=Nfft_r/2; j<Nfft_r; j++){ 
	tmp_complex_array2[j]*=conj(fft_matched_filter[j]);
      }
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Multiply in frequency domain by conj(ref_func)");
	cbf.set("hrcdata3",tmp_complex_array2,Nfft_r,i+1,-1);
        cbf.set("mf_fft",fft_matched_filter,Nfft_r/2,i+1,-1);
      }      

      // baseband (rotate spectra)
      for(int j=Nfft_r/4;j<Nfft_r/2;j++) 
	tmp_complex_array2[i]=tmp_complex_array2[i+Nfft_r/2];


      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Rotate Spectra");
	cbf.set("hrcdata4",tmp_complex_array2,Nfft_r,i+1,-1);
      }      
      
      // inverse FFT 
      ifft(&tmp_complex_array2[Nfft_r/4],tmp_complex_array,Nfft_r/2);
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Inverse FFT");
	cbf.set("hrcdata5",tmp_complex_array,Nfft_r/2,i+1,-1);
      }      
      
      // copy to storage array
      
      for(int j=0;j<Nfft_r/2;j++){
	dechirped[j][i]=tmp_complex_array[j]; 
	float val=abs(tmp_complex_array[j]);
        val*=val;
	sum_square_fulldechirp+=val;
        if(j>=start_range_bin && j <= end_range_bin) sum_square_dechirp+=val;
      }
    }
    float output_energy=0.0;
    for(int j=0;j<Nfft_r/2;j++) 
      output_energy+=abs(dechirped[i][j])*abs(dechirped[i][j]);
    if(dbg.level)
      dbg.file << "RangeCompress Pulse " << i << " Input Energy= " 
	       << input_energy << " Output Energy= " << output_energy
	       << endl;
  }
  rangeCompressGain=sum_square_dechirp/mean_square_inputs;
  float rangeTruncationLoss=sum_square_fulldechirp/sum_square_dechirp;
  overlap_ratio=sum_square_overlapped_inputs/mean_square_inputs;
  mean_square_inputs/=(num_samples_used);
  afs << " Range Compressed. Mean Square Inputs is: " << mean_square_inputs
      << ", Segmentation Overlap Energy Ratio is:" << overlap_ratio << endl;
  afs << " rangeCompressGain is: " << rangeCompressGain << endl;
  afs << " range truncation loss is:" <<  rangeTruncationLoss << endl;
}

//---------------------------------------------------------------------
// L1I::azimuthCompress
// Multiply range compressed arrays by complex exponential to perform 
// deramp and doppler centering for each iso-range line then perform
// forward FFT
//---------------------------------------------------------------------
void
L1I::azimuthCompress(const SARProcParams& spp, ofstream& afs)
{
  DebugInfo dbg("L1I::azimuthCompress");
  float dopshift=get_in_base_units(ac_dop_shift);
//------------------------------
// azimuth compress for each range line
//------------------------------
  for(unsigned int k=0;k<num_range_bins;k++){
    


    // compute reference time for first window in seconds
    float t0=pri_in_s/2;
    
    
    float freq_offset;
    
    // OLD WAY doppler centroid estimate effects geo-location
    if(quick_look_dopran)
      freq_offset=2*pi*(fdop_centroid[k]+ dopshift -prf_in_Hz/2);
    
    // NEW way doppler centroid estimate only effects gain
    else
      freq_offset=2*pi*(fdop_centroid[k]-prf_in_Hz/2);

    float freq_rate=pi*fdop_rate[k];

    // compute reference function
    for(unsigned int p=0;p<num_pulses_received;p++){
      float t=t0+pri_in_s*p;
      float r_phase= -(freq_rate*t+freq_offset)*t;
      // Azimuth Ref function is complex conjugate of the doppler "chirp"
      // Multiplying by this is the same and correlation in the frequency
      // domain
      azimuth_ref_func[p]=complex<float>(cos(r_phase),std::sin(r_phase));
    }



    //-----------------------------------------------------------------     
    //     compute scale of reference function 
    //-----------------------------------------------------------------     
      
    double fftgain=0;
    tmp_complex_array[0]=complex<float>(1,0);
    for(int i=1;i<Nfft_az;i++) tmp_complex_array[i]=complex<float>(0,0);

    fft(tmp_complex_array,tmp_complex_array2,Nfft_az);

    double sfft=0;
    for(int i=0;i<Nfft_az;i++) sfft+=abs(tmp_complex_array2[i])*abs(tmp_complex_array2[i]);
    fftgain=sfft;

    // compute smoothing window in freq domain
    float a=azimuth_ref_window_param;
    float b=1-a;  
    


  
    // multiply reference by windowing function
    
    

    double sum2=0,sum4=0;
    for(unsigned int p=0;p<num_pulses_received;p++){
      tmp_complex_array[p] = a - b * cos((2.* pi*p)/num_pulses_received);
      azimuth_ref_func[p] *= tmp_complex_array[p]; 
      float dv=abs(azimuth_ref_func[p]*azimuth_ref_func[p]);
      sum4+=dv*dv;
      sum2+=dv;
    }
    
    double scale=sqrt(sum2)/sqrt(sum4*fftgain);

    // scale ref func
    for(unsigned int p=0;p<num_pulses_received;p++){
      azimuth_ref_func[p]*=scale;
    }
    float input_energy=0.0;
    for(unsigned int p=0;p<num_pulses_received;p++) 
      input_energy+=abs(dechirped[k][p])*abs(dechirped[k][p]);
    // get array of dechirped signal for range bin
    // take complex conjugate for multiplying with azimuth ref function
    for(unsigned int p=0;p<num_pulses_received;p++){
      tmp_complex_array[p]=conj(dechirped[k][p]);
    } 



    // zero pad signal
    for(unsigned int p=num_pulses_received;p<(unsigned int)Nfft_az;p++){
      tmp_complex_array[p]=0.0;
    }  

  
    // multiply reference function by signal
    for(unsigned int p=0;p<num_pulses_received;p++){
      tmp_complex_array[p]*=azimuth_ref_func[p];
    }
  

    if(check_this_burst){
      cbf.comment("");
      cbf.comment("Multiply dechirped data by azimuth reference function");
      cbf.set("acdata1",tmp_complex_array,Nfft_az,-1,k+1);
      cbf.set("azref",azimuth_ref_func,Nfft_az,-1,k+1);
    }
    // forward fft
    fft(tmp_complex_array,tmp_complex_array2, Nfft_az); 

   if(check_this_burst){
      cbf.comment("");
      cbf.comment("Take FFT");
      float acscale=(float)scale;
      float acfftgain=(float)fftgain;
      cbf.set("acscale",&acscale,1,1,k+1);
      cbf.set("acfftgain",&acfftgain,1,1,k+1);
      cbf.set("acdata2",tmp_complex_array2,Nfft_az,-1,k+1);
    }
    // store values in array
    for(int p=0;p<Nfft_az;p++){
      raw_sar_data[k][p]=tmp_complex_array2[p];
    }

    float output_energy=0.0;
    for(int p=0;p<Nfft_az;p++) 
      output_energy+=abs(raw_sar_data[k][p])*abs(raw_sar_data[k][p]);

    if(dbg.level)
      dbg.file << "AzimtuhCompress Range Bin  " << k << " Input Energy= " 
	       << input_energy << " Output Energy= " << output_energy
	       << endl;
  }
  // compute sum square of compressed pixels for use in normalization
  image_sum=0.0;
  for(int k=0;k<(int)num_range_bins;k++){
    for(int p=0;p<(int)Nfft_az;p++){
      if(k>=start_range_bin && k<=end_range_bin){
	float detected_value=abs(raw_sar_data[k][p]);
	detected_value*=detected_value;
	image_sum+=detected_value;
      }
    }
  } 
  azimuthCompressGain=image_sum/sum_square_dechirp;
  afs << " Azimuth Compressed. Sum of detected image pixels is " 
      << image_sum << endl;
  afs << "  azimuthCompressionGain = " << azimuthCompressGain << endl;
}

void
L1I::computeUsableExtent(const SARProcParams& spp)
{

  // assign fdop_min and fdop_max
  fdop_max=0.8*prf_in_Hz;
  fdop_min=-fdop_max;
 
  // assign range_min and range_max
  float range_width=(-nominal_pulse_spread*speed_light/2).km();
  range_min=-range_width/2;
  range_max=-range_width/2;  
}

void
L1I::removeDopplerPhaseRamp(const SARProcParams& spp, ofstream& afs){
  DebugInfo  dbg("L1I::removeDopplerPhaseRamp");
  float w;  // theoretical doppler phase ramp;
  w=-pi*(num_pulses_received-1)/Nfft_az;

  float phr; // empirical range phase ramp
  float pha; // empirical doppler phase ramp; 
  if(spp.remove_doppler_phase_ramp){
    // compute theoretical phase ramp (in rad/pix from range error and 
    // transmit frequency

    // compute xgoto-lchange in range per pixel in doppler direction
    // (this value is non-zero because of doppler induced error in range 
    // compression


 

    // compute empirical range phase ramp
    complex<float> cr(0.0,0.0);
    for(unsigned int i = 1; i<num_range_bins;i++){
      for(int j = 1;j<Nfft_az;j++){
        if(Xfactor[i][j]!=0.0)
	  cr = cr + sigma0[i][j]*conj(sigma0[i-1][j]);
      }
    }


    phr = atan2(imag(cr),real(cr));

    
    // compute empirical azimuth phase ramp
    complex<float> ca(0.0,0.0);
    for(unsigned int i = 1; i<num_range_bins;i++){
      for(int j = 1;j<Nfft_az;j++){
	if(Xfactor[i][j]!=0.0)
	  ca = ca + sigma0[i][j]*conj(sigma0[i][j-1]);
      }
    }

    
    pha = atan2(imag(ca),real(ca));

    

    // precompute doppler phase ramp corrections
    for(int j=0;j<Nfft_az;j++){
      tmp_complex_array[j]=complex<float>(cos(j*pha),-sin(j*pha));
    }

    // remove range and doppler phase ramps

    for(unsigned int i=0;i<num_range_bins;i++){
      tmp_complex_array2[i]=complex<float>(cos(i*phr),-sin(i*phr));
      for(int j=0;j<Nfft_az;j++){
	sigma0[i][j]*=tmp_complex_array[j]*tmp_complex_array2[i];
      }
    }
  }
  else{
    for(unsigned int i=0;i<num_range_bins;i++){
      for(int j=0;j<Nfft_az;j++){
	sigma0[i][j]*=complex<float>(cos(j*w),-sin(j*w));
      }
    }      
  }
  if(dbg.level){
    dbg.file << "L1I::removeDopplerPhaseRamp SabCounter:" << sab_counter
	     << " Theoretical doppler phase ramp = " << w << " rad/pix " 
	     << endl;
    if(spp.remove_doppler_phase_ramp){
      dbg.file <<  " Empirical doppler phase ramp = " << pha << endl 
	       <<  " Empirical range phase ramp = " << phr << endl;
    }
  }
}

float L1I::getSARProcGain()
{
  double retval=0.0;
  for(int c=0;c<Nfft_r/2;c++){
    float v=abs(fft_matched_filter[c]);
    retval+=v*v;
  }
  retval/=Nfft_r/2;
  retval*=2*(float)num_pulses_received/(float)Nfft_az;
  return((float)retval);
}

int L1I::getGoodSigma0s(float* val, int maxlen){
  int retval=0; 
   for(int i=0; i<(int)num_range_bins; i++){
      for(int j=0;j<Nfft_az;j++){
	if(Xfactor[i][j]!=0){
          if(retval>=maxlen){
	     fprintf(stderr,"L1I::getGoodSigma0s::Memory allocation error: more than %d good pixels\n",maxlen);
	    exit(2);
	  }
	  val[retval]=abs(raw_sar_data[i][j]);
	  val[retval]*=val[retval];
	  retval++;
	}
      }
    }
   return(retval);
}


float
L1I::correlateRawWithX(){
  float corr=0;
  float energy_raw=0;
  float energy_X=0;
  float mean_raw=0;
  float mean_X=0;
  int nr=end_range_bin-start_range_bin+1;
  for(int i=start_range_bin;i<=end_range_bin;i++){
    for(int j = 0; j<Nfft_az;j++){
      double val=abs(raw_sar_data[i][j]);
      val*=val;
      mean_raw+=val;
      mean_X+=Xfactor[i][j];
    }
  }
  mean_raw=mean_raw/(Nfft_az*nr);
  mean_X=mean_X/(Nfft_az*nr);
  for(int i=start_range_bin;i<=end_range_bin;i++){
    for(int j = 0; j<Nfft_az;j++){
      double val=abs(raw_sar_data[i][j]);
      val*=val;
      val=val-mean_raw;
      double val2=Xfactor[i][j]-mean_X;
      corr+=val*val2;
      energy_raw+=val*val;
      energy_X+=val2*val2;
    }    
  }
  corr/=sqrt(energy_X*energy_raw);
  return(corr);
}

int
// tolerance is a HARD_CODED HACK 0.1 is typical 0.01 is for closest approach
// Enceladus SAR
L1I::HISARcalibrate(SARProcParams& spp, double radii[3], double  freq_shift, int ambigmode, double tol=0.01){
  
  // ------------
  // Routine loops over antenna pattern space
  // determines  1) which range doppler pixel each azim elev pixel contributes to
  // and 2) whether it is in the processing window and thus contributes to X OR
  // outside the processing window and thus contributing to Xamb
  // This routine is called instead of calibrate and does not have many of the bells and whistles available
  // to calibrate, i.e., no ability to adjust calibration for attitude error (mo Doppler and Range correction)
  // only works for triaxial body NOT for arbitrary topography.
  //
  // Routine has two modes of operation 
  // 1) ambigmode = 0 (Xfactor is set to zero when Xambig/Xfactor exceeds the configured threshold) 
  // 2) ambigmode = 1 ( Xfactor is set to Xfactor+Xambig to account for all ambiguous contributions)
  // 3) ambigmode = 2 same as 1 but a bcakscatter model is incorporated into X
  // this third  mode is most useful when we are doing a search to determine radii and/or lambda_offset
  // 4) ambigmode=3 Xfactor is set to zero when Xambig*ests0ambig/Xfactor/ests0 > thresh

  //---------------------------
  // Determine which target is in use
  // in order to set backscatter model
  // if desired
  //----------------------------


  corrTypeE backscat_mod_type;
  if(strcasecmp("Titan",default_target_name.c_str())==0){
    backscat_mod_type=HAGHAG;
  }
  else if(strcasecmp("Enceladus",default_target_name.c_str())==0){
    backscat_mod_type=ENCELADUS_WYE1;
  }
  else if(strcasecmp("Rhea",default_target_name.c_str())==0){
    backscat_mod_type=RHEA_WYE1;
  }
  else{
    backscat_mod_type=INVALID_CORR;
  }

  if(ambigmode>=2 && backscat_mod_type==INVALID_CORR){
    cerr << "Error: L1I::HISARcalibrate AmbigMode=2 but target " << default_target_name << " has no assigned backscatter model " << endl;
    exit(1);
  }

  //------------------
  // Toss burst if echo buffer is all zeros
  //-----------------
  if(mean_square_inputs==0) return(0);

  //------------------------------------------
  // compute constant part of X and Xambig
  //-----------------------------------------
  float norm_value = (float) 1.0/(float)(num_signal_samples);
  float cag=computeCalibratedAttenuatorGain();
  Uvar Xcal=spp.Xcal(t,beam_number,r_mode);
  float Xconst=get_in_base_units(Xcal*cag/norm_value);


  //---------------------------------------
  // Initialize Xfactor amd Xambig
  //---------------------------------------
  for(unsigned int i=0;i<num_range_bins;i++){
    for(int j=0;j<Nfft_az;j++){
      Xfactor[i][j]=0;
      Xambig[i][j]=0;
      ests0[i][j]=0;
      ests0ambig[i][j]=0;
    }
  }

  //------------------------------------------------
  // Determine bounds for azimuth and elevation loop
  // Presume twice the nominal beam widths is OK
  //------------------------------------------------
  double amax=1.5*get_in_base_units(beam_[beam_number-1].getAzimuthWidthOneWay());
  double amin=-amax;
  double emax=1.5*get_in_base_units(beam_[beam_number-1].getElevationWidthOneWay());
  double emin=-emax;
  
  
  //-----------------------------------------------
  // Determine azimuth and elevation resolution 
  // and number of integration steps
  // for now ~100 m by default (near boresight since bins are equal angle this becomes much coarser near limb)
  //-----------------------------------------------


  double u_bore[3]={0,0,1};
  double u_tbf[3];
  double de=0.00001;
  double da=0.00001;
  double z=sqrt(1-de*de);
  double u_bore_de[3]={0,de,z};
  double u_tbf_de[3];
  double u_bore_da[3]={da,0,z};
  double u_tbf_da[3];
  rotate_vector(u_bore,rotmat_rev_,u_tbf); // compute boresight vector in target body fixed frame
  rotate_vector(u_bore_de,rotmat_rev_,u_tbf_de); // compute vector de away from boresight in elevation in target body fixed frame
  rotate_vector(u_bore_da,rotmat_rev_,u_tbf_da); // compute vector da away from boresight in azimuth in target body fixed frame

  double pos_bore[3];
  double pos_bore_de[3];
  double pos_bore_da[3];

  // compute surface locations on triaxial body
  double dummyr;
  int found=get_surface_intercept_triaxial(scpos_,u_tbf,radii,pos_bore,dummyr);
  found = found && get_surface_intercept_triaxial(scpos_,u_tbf_de,radii,pos_bore_de,dummyr);
  found = found && get_surface_intercept_triaxial(scpos_,u_tbf_da,radii,pos_bore_da,dummyr);

  if(!found){
    // for now if boresight is off target throw out burst
    return(0);
  }

  float g2thresh=beam_[beam_number-1].getMaxGain();
  g2thresh*=g2thresh*min_oneway_gain*min_oneway_gain;
  
  float nadir_exclude_ang_float = get_in_base_units(nadir_exclude_ang);
  double ares=tol*da/vector_dist(pos_bore,pos_bore_da);
  double eres=tol*de/vector_dist(pos_bore,pos_bore_de);
  int Nazim=(int)((amax-amin)/ares +0.5);
  int Nelev=(int)((emax-emin)/eres +0.5);

  ares=(amax-amin)/Nazim;
  eres=(emax-emin)/Nelev;

  //------------------------------
  // Loop in azimuth over elevation integration bins
  //-------------------------------
  double look[3], look_tbf[3], pos[3];
 
  for(int a=0;a<Nazim;a++){
    for(int e=0;e<Nelev;e++){
      double azim=amin+a*ares;
      double elev=emin+e*eres;
      bool gaintoosmall=false;
      bool goodside=true;
      //-------------------
      // Determine range and Doppler for center of bin
      //-------------------
      double r, fdop;
      look[1]=sin(elev);
      look[0]=sin(azim)*cos(elev);
      look[2]=sqrt(1-look[0]*look[0]-look[1]*look[1]);


      rotate_vector(look,rotmat_rev_,look_tbf);
      if(!get_surface_intercept_triaxial(scpos_,look_tbf,radii,pos,r)){
	continue; //skip bins with centers off the limb;
      }
      if(exclude_mirror) goodside=goodSideOfNadir(pos,nadir_exclude_ang_float);
      double frequency=speed_light_in_kps/lambda_chirp_in_km+freq_shift;
      double vr=scvel_[0]*look_tbf[0]+scvel_[1]*look_tbf[1]+scvel_[2]*look_tbf[2];
      fdop=2*vr*frequency/(speed_light_in_kps-vr);

      //--------------------
      // Determine antenna gain for center of bin
      //---------------------
      double g2=beam_[beam_number-1].bilinear(azim,elev);
      g2*=g2;

      if(g2<g2thresh) gaintoosmall=true;
      //--------------------
      // Determine surface area of bin
      //--------------------
      float area=getBinAreaOnTriaxialBody(radii,azim,ares,elev,eres);
      
      //---------------------
      // Compute dX for bin
      //---------------------
      float r4=r*r*r*r;
      float dX=Xconst*g2*area/r4;
      float ds0=0;
      //---------------------
      // Include backscatter model in dX if ambigmode==2
      //---------------------
      if(ambigmode==2){
        double unorm[3]={pos[0]/radii[0],pos[1]/radii[1],pos[2]/radii[2]};
        float inc=acos(-dot(look_tbf,unorm));
	float s0=theorBackscatter(inc,backscat_mod_type);
        dX*=s0;
      }
      if(ambigmode==3){
	double unorm[3]={pos[0]/radii[0],pos[1]/radii[1],pos[2]/radii[2]};
        float inc=acos(-dot(look_tbf,unorm));
	float s0=theorBackscatter(inc,backscat_mod_type);
        ds0+=s0;
      }
      //----------------------
      // Determine SAR pixel indices that include energy of bin
      // and whether bin in in processing window
      // Update Xambig and Xfactor
      // accordingly 
      //---------------------

       // account for doppler induced error in range
       // converts true range value to range compression range estimate
       r-=range_doppler_correction_coeff*(fdop-nominal_doppler_centroid_in_Hz);
       bool inprocwindow=true;
       float range_min=range[0];
       float relrange=(r-range_min)/range_step_;
       int i = int(floor(relrange +0.5));


       // check for range_index out of bounds rotate indices by mutiples of processor window widths to
       // compute ambiguous energy if necessary
       while(i<start_range_bin){
	 inprocwindow=false;
         i+=window_step/2;
       }
       while(i> end_range_bin){
	 inprocwindow=false;
         i-=window_step/2;
       }

       double reldop=fdop-nominal_doppler_centroid_in_Hz;
       double mindop=-prf_in_Hz/2;
       float jf= (reldop-mindop)/doppler_step_;
       int j=(int)(floor(jf +0.5));    
       jf=jf-j;

       

       while(j<0){
	 inprocwindow=false;
         j%=Nfft_az;
         j+=Nfft_az;
       }
       while(j> Nfft_az-1){
	 inprocwindow=false;
         j%=Nfft_az;
       }
       // add dX to multiple Doppler bins to account for actual Doppler resolution (as opposed to bin resolution)
       // we multiply by the inverse of the resolution factor to make this an unbiased adjustment
       float inv_resolution_factor = (float)num_pulses_received/(float)Nfft_az;
       float resolution_factor = 1/inv_resolution_factor;
       dX=dX*inv_resolution_factor;
       ds0=ds0*inv_resolution_factor;
       int minj= (int)floor(j+jf-resolution_factor/2+0.5);
       int maxj= (int)floor(j+jf+resolution_factor/2+0.5);
       for(int j2=minj;j2<maxj;j2++){
         if(j2<0) j2=j2+Nfft_az;
         if(j2>Nfft_az-1) j2=j2-Nfft_az;
	 if(inprocwindow) Xfactor[i][j2]+=dX;
	 else Xambig[i][j2]+=dX;
       
       
	 if(ambigmode==3){
	   if(inprocwindow) ests0[i][j2]+=ds0;
	   else ests0ambig[i][j2]+=ds0;
	 }

	 // handle gaintoosmall case
	 // Odd way to do it but this should force the ambiguity to signal ratio to fail the threshold check
	 // thereby assuring that poor gain and mirror ambig pixels are excluded (may also exclude limb pixels, but that's good too.)
	 if(inprocwindow && (gaintoosmall || !goodside) && (ambigmode==0 || ambigmode==3)) {
	   Xambig[i][j2]+=dX;
	 }
       } // end of j2 loop which accounts for blurring of X due to true doppler resolution


    } // end elevation loop
  } // end azimuth loop


  
  //--------------------------
  // Set outputs based on ambigmode
  //--------------------------
  for(unsigned int i=0;i<num_range_bins;i++){
    for(int j=0;j<Nfft_az;j++){
      if(ambigmode==1 || ambigmode==2) Xfactor[i][j]+=Xambig[i][j];

      if(Xfactor[i][j]==0) sigma0[i][j]=complex<float>(0,0);
      else sigma0[i][j]=raw_sar_data[i][j]/(float)sqrt(Xfactor[i][j]);

      if(ambigmode==3){
	double ratio=Xfactor[i][j]*ests0[i][j]/(Xambig[i][j]*ests0ambig[i][j]);
	if(ratio<signal_to_amb_thresh || Xfactor[i][j]==0){
	  Xfactor[i][j]=0;
	}
      }
      if(ambigmode==0){
	double ratio=Xfactor[i][j]/Xambig[i][j];
	if(ratio<signal_to_amb_thresh || Xfactor[i][j]==0 ){
	  Xfactor[i][j]=0;
	}
      }
    }
  }


  //-------------------------------
  // Write out intermediate results if desired
  //-------------------------------

  if(sigma0fp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(sigma0[i][j]),sizeof(complex<float>),1,sigma0fp);
      }
    }
  }




  if(Xfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(Xfactor[i][j]),sizeof(float),1,Xfp);
      }
    }
  }

  if(Xambfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(Xambig[i][j]),sizeof(float),1,Xambfp);
      }
    }
  }

  if(dcfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(dechirped[i][j]),sizeof(complex<float>),1,dcfp);
      }
    }
  }

  if(rspecfp!=NULL){

    for(int i=0;i<Nfft_r*2;i++){ 
      tmp_complex_array[i]=0;
    }     
    for(unsigned int j=0;j<num_pulses_received;j++){
      for(int i=0;i<Nfft_r;i++){
	tmp_complex_array[i]=dechirped[i][j];	
      }
      fft(tmp_complex_array,tmp_complex_array2,Nfft_r);
      for(int i=0;i<Nfft_r/2;i++){
	tmp_complex_array[i+Nfft_r]+=abs(tmp_complex_array2[Nfft_r/2+i])*abs(tmp_complex_array2[Nfft_r/2+i]);
	tmp_complex_array[i+3*Nfft_r/2]+=abs(tmp_complex_array2[i])*abs(tmp_complex_array2[i]);
      }
    }
    float v;
    for(int i=0;i<Nfft_r;i++){
      v=abs(tmp_complex_array[i+Nfft_r])/num_pulses_received;
      fprintf(rspecfp,"%g %g\n",i*(adc_in_Hz/Nfft_r/2)-adc_in_Hz/4,
	      v);
    }
    fprintf(rspecfp,"&\n");
    complex<float> zc(0,0);
    for(int i=0;i<Nfft_r*2;i++){
      tmp_complex_array[i]=zc;
      tmp_complex_array2[i]=zc;
    }
  }
  
  float corrected_range;
  float range_correction;
  if(rangefp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
        if(i>=(int)num_range_bins || j>=(int)Nfft_az) corrected_range=0;
	else{
	  range_correction=dop_buf[i][j]-nominal_doppler_centroid_in_Hz;
	  range_correction*=range_doppler_correction_coeff;
	  corrected_range=(float)(range[i]+range_correction);
	}
	fwrite(&(corrected_range),sizeof(float),1,rangefp);
      }
    }
  }

  if(rfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(raw_sar_data[i][j]),sizeof(complex<float>),1,rfp);
      }
    }
  }
  float dop;
  if(fdopfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
        if(i<start_range_bin || i>end_range_bin || j>Nfft_az)
	  dop=absDoppler(i,j);
	else dop = 0;
	fwrite(&dop,sizeof(float),1,fdopfp);
      }
    }
  }

  return(1);
}

void
L1I::calibrate(SARProcParams& spp, ofstream& afs, float height_in_km)
{
  // static variable ony used when orthorectify is on
  static float dh=ortho_height_tolerance.getInUnits("km");

  if(use_HISAR_calibrate){
    HISARcalibrate(spp,default_target_radii_in_km,0,HISAR_ambig_mode);
    // add values to standard info line
    spp.cout_line="";
    spp.cout_line+=toStr(sab_counter);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(get_in_base_units(time_from_closest_approach));
    spp.cout_line+=" ";
    double lat=get_in_base_units(act_centroid_lat)*radtodeg;
    double lon=get_in_base_units(act_centroid_lon)*radtodeg;
    lon=360.0-lon; // positive west convention output to match other outputs
    while(lon>360)lon-=360;
    while(lon<0)lon+=360;
    spp.cout_line+=toStr(lat);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(lon);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(beam_number);
    spp.cout_line+=" ";
    return;
  }
  if(quick_look_dopran)
    fprintf(stderr,"WARNING You are using QL1 Doppler and Range Corrections\n");
  fft_interpd_=false;
  if(spp.calibration_on){

    // compute coefficents to get range_res and azimuth_res if X is replaced
     float chirp_bw=get_in_base_units(slow_cfs*(csq-1));
     float range_res=get_in_base_units(speed_light)/(2*chirp_bw);
     float doppler_res=1/(num_pulses_received*pri_in_s);

     float rcoeff=get_in_base_units(spp.range_res_coeff*range_res/range_step_);
     float acoeff=get_in_base_units(spp.azimuth_res_coeff*doppler_res/doppler_step_);

    // Special functionality for noise only processing 
    // to "beat down the noise on noise estimates" hah

    if(spp.signal_type!=SARProcParams::NOMINAL){
      cerr << "Warning: performing noise only processing" << endl;
      // average raw_sar_data
      double average=0.0;
      for(unsigned int i=0;i<num_range_bins;i++){
	for(int j=0;j<Nfft_az;j++){
	  double val=abs(raw_sar_data[i][j]);
	  
	  average+=val*val;
	}
      }
      average=sqrt(average/(num_range_bins*Nfft_az));
      for(unsigned int i=0;i<num_range_bins;i++){
	for(int j=0;j<Nfft_az;j++){
	  raw_sar_data[i][j]=complex<float>((float)average,0.0);
	}
      }
    }

    // compute normalization value needed so that detected pixel values are
    // sum to the received power in square data numbers
    //    float norm_value=mean_square_inputs*overlap_ratio/image_sum
    //  /effective_duty_cycle;
    float norm_value;

    float datdepnv=1/(azimuthCompressGain*rangeCompressGain*num_signal_samples);//    float constnv=(float)Nfft_r/(float)(2.0*num_signal_samples*Nfft_az);
    float constnv = (float) 1.0/(float)(num_signal_samples);

    
    // Nfft_az term should be computed using fftgain method in azimuthCompress
    // to allow fft library usage if desired
    if(spp.enable_data_dependent_cal){
      norm_value=datdepnv;
    }

    else{
      norm_value=constnv;
    }
    /**** Commented out calibration debugging stuff
    cout << "datadepnorm_val " << datdepnv << " constnorm_val " << constnv
	 << endl;
    cout << "AzimuthCompressGain " << azimuthCompressGain << endl;
    cout << "RangeCompressGain " << rangeCompressGain << endl;
    cout << "Nfft_range " << Nfft_r << endl;
    cout << "Nfft_azim " << Nfft_az << endl;
    cout << "num_signal_samples " << num_signal_samples << endl;
    *****/
    float cag=1.0;
    DebugInfo computeGain2_dbg("L1I::computeGain2");

    // calibrate routine debugging
    // level =1 standard calibration debug
    // level =2 usability computation debugging 
    DebugInfo dbg("L1I::calibrate");
 
    if(dbg.level)
      dbg.file << "azimuthCompressGain=" << azimuthCompressGain << " rangeCompressGain=" 
	       << rangeCompressGain << " num_signal_samples=" << num_signal_samples << endl;
    dbg.file << "Nfft_az=" << Nfft_az << endl;
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

    // compute attenuator gain 
    // NOW like everything else this is power to power
    cag=computeCalibratedAttenuatorGain();
    Uvar Xcal=spp.Xcal(t,beam_number,r_mode);
    Uvar pow2dat=spp.powerToDataNumber(r_mode);
    // Get gain threshold to apply to compute useable extent
    float g2max=0; // peak two-way gain
    float g2thresh_rel=0; // useable area gain threshold with respect to max gain (linear scale)
    float g2thresh=0; // absolute useable area gain threshold  (linear_scale)

    g2max=beam_[beam_number-1].getMaxGain(); // peak one-way gain
    g2max*=g2max; // peak two-way gain

    // compute absolute minimum gain (maximal relative gain correction)
    //  to apply 
    // I.e. for all points below maximal_gain_corr dB down from the peak of
    // the oneway gain pattern the 2*maximal_gain_corr dB or correction 
    // will be  applied in the calibration and no more. This is used to avoid
    // garbage leaking into the edges of the useable area through the sinc
    // interpolator, If this number is too high we can get this garbage.
    // If it is too low, features near the edge will be darkened and slightly
    // squeezed together (due to losing energy on one side of the sinc).
    double abs_min_gain= g2max/pow(10.0,0.2*spp.maximal_gain_corr);
    
    if(spp.gain_cutoff_enabled){
      // specified linear scale two-way gain thresholds
      // the region of gain higher than this threshold defines the
      // useable area
      g2thresh_rel=spp.getGain2Thresh(ground_impact_time, beam_number,
				      sab_counter);
      g2thresh=g2thresh_rel*g2max;
    }
    else if ( spp.full_useable_calc_enabled){
      int num_beams=5; // this is a hack but it should be OK unless beam 3 
      // only SAR gets used in which case it should be conservative
      double look_repeat_time_in_s=bpd.getInUnits("s")*num_beams; 
      spp.initUseableCalculation(*this,scstate_,
				 nominal_delay*speed_light/2.0-Uvar(rc_range_shift_in_km,"km"),
				 nominal_doppler_centroid-Uvar(ac_dop_shift_in_Hz,"Hz"),lambda_chirp_in_km,
				 look_repeat_time_in_s);
    }
  


    float alt;
    if(dbg.level==1){
      alt=sqrt(scpos_[0]*scpos_[0]+scpos_[1]*scpos_[1]
		     +scpos_[2]*scpos_[2]);
      alt=alt-target_radius_;
      dbg.file << "L1I::calibrate: Burst Number " << burst_no 
	       << " Altitude " << alt << " km" 
	       << "Beam " << beam_number <<  endl;
    }

    // Values for testing Useability Calculation
    float max_bad_g2rel=-20, min_good_g2rel=0, g2rel;
    int num_bad=0,num_good=0;
    float mbdop,mbran,mgdop,mgran;





    for(unsigned int i=0;i<num_range_bins;i++){
      for(int j=0;j<Nfft_az;j++){
	// at this point Xfactor should be 1/R^4;
	
	
	// May take dop and range as input instead of ij

	double doppler_in_hz=absDoppler(i,j); 
	double range_in_km=range[i];     
	float azdeg, eldeg;
        float g2;
        // Old Way Only Gain is effected.doppler and range centroid estimates
        // effect geo-location


        // if orthorectification is on find correct height in km.
        // the drawback to doing it this way is that the same height is still used for the whole pixel
        // BWS Oct 18 2010
        if(orthorectify){
	  float lastheight = height_in_km;
          double lat,lon;
          bool topomap_inbounds=false;
	  while(1){
	    bool onbody=dopplerRangeHeightToLatLon(doppler_in_hz,range_in_km,height_in_km,lat,lon);
         
            // try again with zero height if location was not on body 
            if(!onbody){
              height_in_km=0;
	      onbody=dopplerRangeHeightToLatLon(doppler_in_hz,range_in_km,height_in_km,lat,lon);	      
	    }

            // if still not on body drop out with height = 0
            // gain and area calculations will compute 0 values later on.
            if(!onbody) break;


	    height_in_km=ortho_tmap.getInterpolatedHeightInKm(lat,lon,topomap_inbounds);
	    if(!topomap_inbounds){
	      fprintf(stderr,"L1I::calibrate Bad TopoMap for Orthorectification (need map that covers whole imaging region)\n");
	      exit(1);
	    }

            // break when height tolerance achieved
	    if(fabs(height_in_km-lastheight)<dh) break;

            lastheight=height_in_km; // update previous height estimate holder

	  } // end height search while loop
	} // end orthorectification on case

        if(quick_look_dopran)
	  g2=computeGain2(doppler_in_hz,range_in_km, computeGain2_dbg,
			      azdeg,eldeg,height_in_km); 
        // New Way Only Gain is effected.
        else
	  g2=computeGain2(doppler_in_hz-ac_dop_shift_in_Hz,
			  range_in_km-rc_range_shift_in_km, computeGain2_dbg,
			      azdeg,eldeg,height_in_km); 

        if(g2<abs_min_gain) g2=abs_min_gain;

        float area;
	area=getPixelAreaOnSphere(doppler_in_hz,range_in_km,height_in_km);
	// for debugging only
	//#define AREA_DEBUG
#ifdef AREA_DEBUG
	// compute u
        double shouldbe_radius;
	double normal[3];
	float inc_angle;
	double u[3], u1[3], u2[3];
        float area_est;
	// get two possible looks in TBF coordinates
	// most parameters to this function were computed in 
	// initGeometryCalculations.
	// If no such look exists retrun 0 gain
	if(! fast_doppler_range_to_TBF_look(doppler_in_hz,range_in_km,
					    target_radius_+height_in_km, scpos_, scvel_, 
					    lambda_special, u1,u2)){
	  area_est=0;
	}
 
	else{
	  // pick the look closest in angle to the boresight
	  // compares dot products (highest wins)
	  // boresight in TBF was precomputed in initGeometryCalculations.
	  // If the SAR image intersects nadar this procedure will fail because a
	  // particular doppler/range bin will contain energy from BOTH mirror looks.
	  // Of course the one closest to the boresight will have the larger gain.
	  double cosangle1=dot(boresight_in_tbf_,u1);
	  double cosangle2=dot(boresight_in_tbf_,u2);
	  for(int c=0;c<3;c++){
	    if(cosangle1>cosangle2) u[c]=u1[c];
	    else u[c]=u2[c];
	  }
	  // compute normal
	  double sum=0;
	  for(int c=0;c<3;c++){
	    normal[c]=scpos_[c]+range_in_km*u[c];
	    sum+=normal[c]*normal[c];
	  }
          
	  sum=sqrt(sum);
	  for(int c=0;c<3;c++){
	    normal[c]/=sum;
	  }
          shouldbe_radius=sum;
          // just for fun make the look vector point up
          u[0]=-u[0];
	  u[1]=-u[1];
	  u[2]=-u[2];
	  inc_angle=acos(dot(u,normal));
	  
	  area_est= range_step_*range_in_km*lambda_chirp_in_km*doppler_step_/
	  (sin(inc_angle) *2 *scspeed_*
	   sqrt(1-pow((lambda_chirp_in_km*doppler_in_hz/2/scspeed_),2)));
	}
        FILE* adfp=fopen("areadebug.txt","a");
	fprintf(adfp,"%g %g %g %g %g %g %g %g %g %g %g\n",area,area_est,inc_angle*180/pi,area_est/area,u[0],u[1],u[2],normal[0],normal[1],normal[2],shouldbe_radius);
	fclose(adfp);

#endif        
        // output info about range and azimuth resolution
	if(dbg.level==1){
          // only output info for central pixel
	  if(i==num_range_bins/2 && j==Nfft_az/2){
	    float chirp_bw=get_in_base_units(slow_cfs*(csq-1));
	    float range_res=get_in_base_units(speed_light)/(2*chirp_bw);
	    float doppler_res=1/(num_pulses_received*pri_in_s);

	    dbg.file << "L1I::calibrate " 
		     << "pixel size (range_km,doppler_hz)="
		     <<"(" << range_step_ << "," << doppler_step_ << ")"
		     << endl;
	 
	    dbg.file << "L1I::calibrate " 
		     << "Effective resolution (range_km,doppler_hz)="
		     <<"(" << range_res << "," << doppler_res << ")"
		     << endl;
            dbg.file << "L1I::calibrate " 
		     << "pixel size on ground(km,km)="
		     <<"(" << range_ground_width_ << "," 
		     << doppler_ground_width_ << ")"
		     << endl;
            float doppler_ground_res=doppler_ground_width_ *
	      doppler_res/doppler_step_;
	    float range_ground_res=range_ground_width_*range_res/range_step_;
            dbg.file << "L1I::calibrate " 
		     << "Effective resolution on ground(km,km)="
		     <<"(" << range_ground_res 
		     << "," 
		     << doppler_ground_res
		     << ")"
		     << endl;
	    dbg.file << "L1I::calibrate " << "area="<< area << " (km km)" 
		     <<endl;

            // compute acute angle of pixel parallelogram
            // 90 if pixel is rectangular or spherical target yields an
            // area greater than the rectangle
            float ratio=area/doppler_ground_width_/range_ground_width_;	   
	    float para_ang;

	    if(ratio>1.0) 
	      {
		para_ang=90;
	      }
	    else{
	      para_ang=asin(ratio)*radtodeg;
	    }
	    dbg.file << "L1I::calibrate " 
                     << "Acute angle in pixel parallelogram is "
		     << para_ang << " degrees" << endl;

            dbg.file<< "L1I::calibrate(Summary_for_XMGR) "
		    << burst_no << " " << beam_number << " "
		    << alt << " " << para_ang << " "
	            << range_ground_res << " " << doppler_ground_res << " "
		    << area << " " << g2 << " " << datdepnv << " " << constnv 
		    << " " << rangeCompressGain << " " << azimuthCompressGain 
	            << " " << pul << " " << num_pulses_received
		    << endl << endl;
	  }
	}

        // Case for Special_Calibration_Debug
	if(dbgsp.level){
	  if ( 0x0001 & dbgsp.level) Xcal=1.0;
	  if ( 0x0002 & dbgsp.level){
	    if(g2<g2max/2.0) g2=0.0;
	    else g2=1.0;
            g2thresh=g2thresh_rel;
	  }
	  if ( 0x0004 & dbgsp.level) area = 1.0;
	  if ( 0x0008 & dbgsp.level) Xfactor[i][j] = 1.0;
	  if ( 0x0010 & dbgsp.level) cag = 1.0;
	  if ( 0x0020 & dbgsp.level) norm_value = 1.0;
	}


	// All echo samples are zero case (will never happen with real data -- we hope)
	if(mean_square_inputs==0){
	  Xfactor[i][j]=0;
          sigma0[i][j]=complex<float>(0.0,0.0);
	}
        else{
	  // cag is not squared because it is a power gain like 
          // everything else;


	  Xfactor[i][j]*=get_in_base_units(Xcal*g2*area)*cag
	    /norm_value;

          // HACK to turn off gain
	  //Xfactor[i][j]*=get_in_base_units(Xcal*area)*cag/norm_value;

	  sigma0[i][j]=raw_sar_data[i][j]/(float)sqrt(Xfactor[i][j]);

          //REPLACE X if desired
          if(replace_X_mode){
            double sinangle,cosangle;
	    switch(replace_X_param){

	      // June 23 2005 corrected NONE case to avoid applying XFACTOR
              // Bug fix was not tested....
	    case NONE:
	      sigma0[i][j]=raw_sar_data[i][j];
	      Xfactor[i][j]=1.0;
	      
	      break;
	    case RANGE_RES:
	      Xfactor[i][j]=range_ground_width_*rcoeff;
	      break;
	    case AZIMUTH_RES:
	      Xfactor[i][j]=doppler_ground_width_*acoeff;
	      break;
	    case LOOKVEL_ANGLE:
	      cosangle=fabs(doppler_in_hz*lambda_chirp_in_km/(2*scspeed_));

	      if(cosangle>1) cosangle=1;
	      Xfactor[i][j]=acos(cosangle)*180/pi;
	      break;
	    case RES_ELEMENT_ANGLE:
	      sinangle=area/(doppler_ground_width_*range_ground_width_);
              if(sinangle>1) sinangle=1;
	      Xfactor[i][j]=asin(sinangle)*180/pi;
	      break;
            case AMBIG_RATIO: // not in dB
	      Xfactor[i][j]=spp.getAmbigRatio(doppler_in_hz,range_in_km);
              break;
            case ORTHO_LOCAL_RAD:
              Xfactor[i][j]=target_radius_+height_in_km;
           case ORTHO_HEIGHT:
	     Xfactor[i][j]=height_in_km*1000; // output in m
	    }
	  }
	}

        // assign special output arrays

	//when this line is not commented out the area output file
	// takes the calibrated, detected s0 value w/o gain correction
//#define SPECIAL_SH_OUTPUT

#ifdef SPECIAL_SH_OUTPUT
	area_buf[i][j]=abs(sigma0[i][j]);
	area_buf[i][j]*=area_buf[i][j];
	area_buf[i][j]*=g2;
#else
        area_buf[i][j]=area;
#endif
        g2_buf[i][j]=g2;
        dop_buf[i][j]=doppler_in_hz;
	azim_buf[i][j]=azdeg;
	elev_buf[i][j]=eldeg;

        // Truncate doppler pixels outside percentage of PRF if
        // azimuth_width_cutoff_enabled;
        if(spp.azimuth_width_cutoff_enabled){
	  float thresh1=(100.0-spp.azimuth_width_percent)*0.005;
	  float thresh2=1-thresh1;
	  if(float(j)/float(Nfft_az-1)<thresh1)
	    Xfactor[i][j]=0;
    	  if(float(j)/float(Nfft_az-1)>thresh2)
	    Xfactor[i][j]=0;
	}

        // Truncate image beyond usable extent if gainCutoffEnabled
        if(spp.gain_cutoff_enabled){
	  if(g2<g2thresh){
	    Xfactor[i][j]=0.0;
	  }
	}  
        else if(spp.full_useable_calc_enabled){
          g2rel=10*log10(g2/g2max);
	  if (!(spp.isUseable(doppler_in_hz-ac_dop_shift_in_Hz,range_in_km-rc_range_shift_in_km))){
	    Xfactor[i][j]=0.0;
	    num_bad++;
	    if(g2rel>max_bad_g2rel){
	      max_bad_g2rel=g2rel;
	      mbdop=doppler_in_hz;
	      mbran=range_in_km;
	    }
	  }
	  else{
	    num_good++;
	    if(g2rel<min_good_g2rel){
	      min_good_g2rel=g2rel;
	      mgdop=doppler_in_hz;
	      mgran=range_in_km;
	    }
	  }   
	}
        // output scaling for boresight pixel to annotation file
	if((int)j==Nfft_az/2 && (int)i==nominal_range_idx){
	  float r4=get_in_base_units(Xcal*g2*area)*cag/norm_value
	    /Xfactor[i][j];

	  afs << " Calibrated Burst Data ..." << endl;
	  afs << " Total Xfactor = " << Xfactor[i][j] 
	      << "= Xcal * G^2 * Area * cag / norm_value / Range^4" << endl;
	  afs << " Xcal = " << get_in_base_units(Xcal) 
              << ", Two-Way Antenna Gain= (G^2) " 
	      << get_in_base_units(g2) << endl;
          afs << " Pixel area = " << get_in_base_units(area) << " sq. km, " 
	      << "Squared Calibrated Attenuator Voltage Gain (cag) = " 
	      << cag << endl;
	  afs << " Normalization Value = " << norm_value 
	      << ", Range^4 = " << r4 << " km^4" << endl;
	}	
	if(dbg.level>=2){
	  dbg.file << "L1I::calibrate (i,j)=(" 
		   << i <<"," << j <<"), g2=" << g2
		   <<  ", area=" << area << ", Xcal=" << Xcal
		   << ",Xfactor=" << Xfactor[i][j] << ", sigma0="
		   << sigma0[i][j] << "beam="
		   << beam_number << endl;
	}
      } // end of azimuth bin loop
    } // end of range bin loop

    
    // output 11 x 11 central detected sigma0 pixels to annotation file
    int min_range=nominal_range_idx-5;
    int max_range=nominal_range_idx+5;
    int min_dop=Nfft_az/2-5;
    int max_dop=Nfft_az/2+5;

    if(spp.full_useable_calc_enabled){
      
      afs << "Full Useability Calculation is Enabled." << endl;
      afs << "Found " << num_good << " useable pixels and " 
	  << num_bad << "unuseable ones." << endl;
      if(num_bad){
	afs << "Maximum One way Gain Contour with unuseable data is " 
	    << max_bad_g2rel/2 << " dB." << endl;
	spp.amb_sar.pointReport(afs,mbran,mbdop,spp.use_multilook_ambiguity);
      }
      if(num_good){
	afs << "Minimum One way Gain Contour with useable data is  "
	    << min_good_g2rel/2 << " dB" << endl;
	spp.amb_sar.pointReport(afs,mgran,mgdop,spp.use_multilook_ambiguity);
      }
    }
    if(min_range<0 || min_dop <0 
       || max_dop>Nfft_az-1 || max_range>(int)num_range_bins-1){
      afs << "!Error: Boresight 11x11 images out of range." << endl;
    }
    else{
      afs << " Producing 11x11 pixel boresight sigma0 dB image." << endl;
      afs << " Left to right = (" << relDoppler(min_dop) 
	  << "," << relDoppler(max_dop) << ") Hz" << endl;
      afs << " Top to Bottom = (" << range[min_range]
	  << "," << range[max_range] << ") km" << endl;
      
      afs.precision(3);
      for(int i=min_range;i<=max_range;i++){
	afs << "    ";
	for(int j=min_dop; j<= max_dop ; j++){
	  float val=abs(sigma0[i][j]);
          afs.width(8);
	  afs << 10*log10(val*val);
	}
	afs << endl;
      }
      afs.precision(5);
    }


    // Output average s0 = integrated power / integrated X  
    // Integration performed for all X greater than Xboresight - 3 dB
    float Xtot=0, Ptot=0, Xmax=0;
    int numpix=0;
    for(int i=0; i<(int)num_range_bins; i++){
      for(int j=0;j<Nfft_az;j++){
	if(Xfactor[i][j]>Xmax) Xmax=Xfactor[i][j];
      }
    }
    float Xthresh=Xmax/2.0;

    for(int i=0; i<(int)num_range_bins; i++){
      for(int j=0;j<Nfft_az;j++){
	if(Xfactor[i][j]<Xthresh) continue;
	Xtot+=Xfactor[i][j];
        float val=abs(raw_sar_data[i][j]);
	Ptot+=val*val;
	numpix++;
      }
    }

    if(output_scatterometer_info){
      total_echo_energy=(Ptot*num_pulses_received*chirp_length*norm_value)/(pow2dat*cag);
      x_factor=(Xtot*spp.quant_noise_scale*num_pulses_received*chirp_length*norm_value)/(pow2dat*cag);
      noise_echo_energy=(spp.quant_noise_scale*spp.thermal_noise_offset+spp.quant_noise_offset)*numpix;
      noise_echo_energy*=(num_pulses_received*chirp_length*norm_value)/(cag*pow2dat*spp.quant_noise_scale);
      sigma0_uncorrected=(Ptot-spp.quant_noise_offset*numpix
			  -spp.thermal_noise_offset*numpix
			  *spp.quant_noise_scale)/Xtot/spp.quant_noise_scale;
  
      Uvar SNR= (total_echo_energy-noise_echo_energy)/noise_echo_energy;
      float neffpix= (float)numpix*(float)num_pulses_received/(float)Nfft_az;
      sigma0_uncorrected_std=(1+ 2/SNR+ 1/SNR/SNR)/neffpix;
      Uvar s0ne=noise_echo_energy/x_factor;
      sigma0_uncorrected_std=(sigma0_uncorrected_std*sigma0_uncorrected*sigma0_uncorrected);
      sigma0_uncorrected_std=MAX(0.04*s0ne*s0ne,sigma0_uncorrected_std);      
      sigma0_uncorrected_std=sqrt(sigma0_uncorrected_std);
      sigma0_corrected=Uvar(incAngleCorrect(sigma0_uncorrected.getInUnits(""),act_incidence_angle.getInUnits("rad"),0,0,HAGHAG));

 // clear Scat data invalid bit (bit 3)
 	science_qual_flag=science_qual_flag & 0xfffffff7;
    }

    float sig0_ave=10*log10(Ptot/Xtot);

    afs << " Sigma0 Average with 3 dB two-way footprint is: "
	<< sig0_ave << endl;
    afs << " Total Power in footprint is: "
	<< Ptot << " (data numbers squared)" << endl;
    afs << " Total X in footprint is: "
	<< Xtot << endl;
    afs << " Number of pixels in footprint is: " << numpix << endl;

    afs << "Apparent image overlap factor is:" << image_sum/Ptot << endl;
    
    // add values to standard info line
    spp.cout_line+=toStr(sab_counter);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(get_in_base_units(time_from_closest_approach));
    spp.cout_line+=" ";
    double lat=get_in_base_units(act_centroid_lat)*radtodeg;
    double lon=get_in_base_units(act_centroid_lon)*radtodeg;
    lon=360.0-lon; // positive west convention output to match other outputs
    while(lon>360)lon-=360;
    while(lon<0)lon+=360;
    spp.cout_line+=toStr(lat);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(lon);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(beam_number);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(sig0_ave);
    spp.cout_line+=" ";
    spp.cout_line+=toStr(10*log10(Xtot));
    spp.cout_line+=" ";
    spp.cout_line+=toStr(num_pulses_received);
    spp.cout_line+=" ";
  }
  // No calibration case; set Xfactor to one; pass raw data to sigma0
  else{
    for(unsigned int i=0;i<num_range_bins;i++){
      for(int j=0;j<Nfft_az;j++){
	Xfactor[i][j]=1;
	sigma0[i][j]=raw_sar_data[i][j];
      }
    }
  }


  data_complete=true;


  // estimate range and doppler offsets and output if desired
  if(estrdfp!=NULL){
    estrd_iskip++;
    if(estrd_iskip%estrd_nskip==0){
    int rangemax_i, dopmax_j;
    float maxE=0, E=0;
    float drangec,ddopc;


    // loop over offsets
    for(int io=-estrd_ranwid/2;io<estrd_ranwid/2;io++){
      for(int jo=-estrd_dopwid/2;jo<estrd_dopwid/2;jo++){
	E=0;
	for(int i=start_range_bin;i<end_range_bin;i++){
	  int i2=i+io;
	  if(i2<0 || i2>L1I_RANGE_WIDTH-1) continue;
	  for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
	    int j2=(j+jo);
	    if(j2<0)j2+=L1I_DOPPLER_WIDTH;
            if(j2>L1I_DOPPLER_WIDTH-1) j2-=L1I_DOPPLER_WIDTH;
	    float dv=abs(raw_sar_data[i2][j2]);
	    dv*=dv*Xfactor[i][j];
	    E+=dv;
	  }
	}
	if(E>maxE){
	  maxE=E;
	  rangemax_i=io;
	  dopmax_j=jo;
	}
      }
    }
    drangec=rangemax_i*range_step_;


    // March 16 2006 modifying this to account for the doppler centroid shift
    // due to the range shift 
    // ORIGINAL VERSION ddopc=dopmax_j*doppler_step_;
    ddopc=dopmax_j*doppler_step_+dop_buf[nominal_range_idx+rangemax_i][0]-dop_buf[nominal_range_idx][0];
    if(maxE!=0){
      fprintf(estrdfp,"%g %d %d %d %d %g %g\n", get_in_base_units(time_from_closest_approach),
	      (int) sab_counter, 
	      (int) beam_number, rangemax_i, dopmax_j, drangec,ddopc);
      fflush(estrdfp);

    }
    }
  }
  // write out intermediate results if desired

  // remove doppler phase ramp empirically if spp boolean is set
  // analytically otherwise
  removeDopplerPhaseRamp(spp,afs);
  
  if(sigma0fp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(sigma0[i][j]),sizeof(complex<float>),1,sigma0fp);
      }
    }
  }




  if(Xfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(Xfactor[i][j]),sizeof(float),1,Xfp);
      }
    }
  }

  if(dcfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(dechirped[i][j]),sizeof(complex<float>),1,dcfp);
      }
    }
  }

  if(rspecfp!=NULL){

    for(int i=0;i<Nfft_r*2;i++){ 
      tmp_complex_array[i]=0;
    }     
    for(unsigned int j=0;j<num_pulses_received;j++){
      for(int i=0;i<Nfft_r;i++){
	tmp_complex_array[i]=dechirped[i][j];	
      }
      fft(tmp_complex_array,tmp_complex_array2,Nfft_r);
      for(int i=0;i<Nfft_r/2;i++){
	tmp_complex_array[i+Nfft_r]+=abs(tmp_complex_array2[Nfft_r/2+i])*abs(tmp_complex_array2[Nfft_r/2+i]);
	tmp_complex_array[i+3*Nfft_r/2]+=abs(tmp_complex_array2[i])*abs(tmp_complex_array2[i]);
      }
    }
    float v;
    for(int i=0;i<Nfft_r;i++){
      v=abs(tmp_complex_array[i+Nfft_r])/num_pulses_received;
      fprintf(rspecfp,"%g %g\n",i*(adc_in_Hz/Nfft_r/2)-adc_in_Hz/4,
	      v);
    }
    fprintf(rspecfp,"&\n");
    complex<float> zc(0,0);
    for(int i=0;i<Nfft_r*2;i++){
      tmp_complex_array[i]=zc;
      tmp_complex_array2[i]=zc;
    }
  }
  
  float corrected_range;
  float range_correction;
  if(rangefp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
        if(i>=(int)num_range_bins || j>=(int)Nfft_az) corrected_range=0;
	else{
	  range_correction=dop_buf[i][j]-nominal_doppler_centroid_in_Hz;
	  range_correction*=range_doppler_correction_coeff;
	  corrected_range=(float)(range[i]+range_correction);
	}
	fwrite(&(corrected_range),sizeof(float),1,rangefp);
      }
    }
  }

  if(rfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(raw_sar_data[i][j]),sizeof(complex<float>),1,rfp);
      }
    }
  }

  if(g2fp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(g2_buf[i][j]),sizeof(float),1,g2fp);
      }
    }
  }

 if(areafp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(area_buf[i][j]),sizeof(float),1,areafp);
      }
    }
  }

 if(fdopfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(dop_buf[i][j]),sizeof(float),1,fdopfp);
      }
    }
  }

  if(azimfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(azim_buf[i][j]),sizeof(float),1,azimfp);
      }
    }
  }

  if(elevfp!=NULL){
    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
	fwrite(&(elev_buf[i][j]),sizeof(float),1,elevfp);
      }
    }
  }
  if(mffp!=NULL){
    for(int i=0;i<L1I_RANGE_WIDTH;i++){
      fwrite(&(fft_matched_filter[i]),sizeof(complex<float>),1,mffp);
    }    
  }

    
}

void
L1I::fftinterp(){
  int nfftdop, nfftrange,numinrange, numindop;
  
  // determine rectangular processing window to interpolate from

  int mini=(int)num_range_bins,maxi=0;
  int minj=Nfft_az, maxj=0;
 
  for(unsigned int i = 1; i<num_range_bins;i++){
    for(int j = 1;j<Nfft_az;j++){
      if(Xfactor[i][j]!=0.0){
	if ((int)i>maxi) maxi=(int)i;
        if ((int)i<mini) mini=(int)i;
	if (j>maxj) maxj=j;
        if (j<minj) minj=j;	
      }	
    }
  }

  fftinterp_range_min=range[mini];
  fftinterp_reldop_min=doppler_step_*minj-prf_in_Hz/2;
  numindop=maxj-minj+1;
  numinrange=maxi-mini+1;
  if(numindop<=0 || numinrange<=0){
    cerr << "Fatal error Sab number " << sab_counter 
         << " computed negative (or zero) fftinterp window size" << endl;
    exit(1);
  }

  complex<float> zc = complex<float>(0.0,0.0);
  for(int i=mini;i<maxi;i++){
    for(int j=minj;j<maxj;j++){
      if(Xfactor[i][j]!=0) interp_sigma0[i-mini][j-minj]=sigma0[i][j];
      else interp_sigma0[i-mini][j-minj]=zc;
    }
  }
  nfftdop=get_next_power_of_2(numindop);   
  nfftrange=get_next_power_of_2(numinrange);  
  int Sdop=FFT_INTERP_SIZE/nfftdop;
  int Srange=FFT_INTERP_SIZE/nfftrange;
  fftinterp_range_step=range_step_/Srange;
  fftinterp_range_step=range_step_/Srange;
  fftinterp_doppler_step=doppler_step_/Sdop;
  
  // perform 2d fft on input window
  for(int i=0;i<nfftdop;i++){
    for(int j=0;j<nfftrange;j++){
      tmp_complex_array[j]=interp_sigma0[i][j];
    }
    fft(tmp_complex_array,tmp_complex_array2,nfftrange);
    for(int j=0;j<nfftrange;j++){
      interp_sigma0[i][j]=conj(tmp_complex_array2[j]);
    }
  }
  for(int j=0;j<nfftrange;j++){
    for(int i=0;i<nfftdop;i++){
      tmp_complex_array[i]=interp_sigma0[i][j];
    }
    fft(tmp_complex_array,tmp_complex_array2,nfftdop);
    for(int i=0;i<nfftdop;i++){
      interp_sigma0[i][j]=conj(tmp_complex_array2[i]);
    }
  }
  // zeropad fft
  float zeropad_loss=Sdop*Srange;
  for(int i=0;i<nfftrange/2;i++){
    for(int j=0;j<nfftdop/2;j++){ 
      interp_sigma0[i][FFT_INTERP_SIZE-j-1]=interp_sigma0[i][nfftdop-j-1];
      interp_sigma0[i][nfftdop-j-1]=zc;
      interp_sigma0[FFT_INTERP_SIZE-i-1][FFT_INTERP_SIZE-j-1]=interp_sigma0[nfftrange-i-1][nfftdop-j-1];
      interp_sigma0[nfftrange-i-1][nfftdop-j-1]=zc;
      interp_sigma0[FFT_INTERP_SIZE-i-1][j]=interp_sigma0[nfftrange-i-1][j];
      interp_sigma0[nfftrange-i-1][j]=zc;
    }
  } 
  // ifft
  for(int i=0;i<FFT_INTERP_SIZE;i++){
    for(int j=0;j<FFT_INTERP_SIZE;j++){
      tmp_complex_array[j]=interp_sigma0[i][j];
    }
    ifft(tmp_complex_array,tmp_complex_array2,FFT_INTERP_SIZE);
    for(int j=0;j<FFT_INTERP_SIZE;j++){
      interp_sigma0[i][j]=conj(tmp_complex_array2[j]);
    }
  }
  for(int j=0;j<FFT_INTERP_SIZE;j++){
    for(int i=0;i<FFT_INTERP_SIZE;i++){
      tmp_complex_array[i]=interp_sigma0[i][j];
    }
    ifft(tmp_complex_array,tmp_complex_array2,FFT_INTERP_SIZE);
    for(int i=0;i<FFT_INTERP_SIZE;i++){
      interp_sigma0[i][j]=conj(tmp_complex_array2[i])*zeropad_loss;
    }
  }
  

  // eliminate extrapolated samples
  fftinterp_num_doppler_bins=(nfftdop-1)*Sdop +1;
  fftinterp_num_range_bins=(nfftrange-1)*Srange +1;

  // bilinearly interpolate fdop_centroid
  for(int i=0;i<fftinterp_num_range_bins;i++){
    int i0=i/Srange;
    int i1=i0+1;
    float ifloat= (float)i/(float)Srange;
    float c0=i1-ifloat;
    float c1=ifloat-i0;
    fftinterp_fdop_centroid[i]=c0*fdop_centroid[i0+mini]+c1*fdop_centroid[i1+mini];
  }
  // output interpolated s0s to file if desired
  if(is0fp!=NULL){
    for(int j=0;j<FFT_INTERP_SIZE;j++){
      for(int i=0;i<FFT_INTERP_SIZE;i++){
	fwrite(&(interp_sigma0[i][j]),sizeof(complex<float>),1,is0fp);
      }
    }
  }
  fft_interpd_=true;
}

void
L1I::pixelDump(OblCylProj& proj){
  static int looknumber=1;
  static int first_pixel=1;
  static double xoff,yoff;

    for(int j=0;j<L1I_DOPPLER_WIDTH;j++){
      for(int i=0;i<L1I_RANGE_WIDTH;i++){
        // skip unusable pixels
        
        // perform correction to get exact range to pixel accounts for range
        // error due to doppler
	double range_correction=dop_buf[i][j]-nominal_doppler_centroid_in_Hz;
	range_correction*=range_doppler_correction_coeff;
	double corrected_range=(float)(range[i]+range_correction);
	if(Xfactor[i][j]==0) continue;
        if(first_pixel){
	  first_pixel=0;
          // compute reference x and y

          getOblCylLatLon(proj,corrected_range,dop_buf[i][j],&xoff,&yoff);
	  cout.precision(20);
	  cout << "Xoff " << xoff*radtodeg << " Yoff " << yoff*radtodeg << endl;
	}
        // compute oblique latitude
        // compute oblique longitude of center pixel
        double x,y;
        getOblCylLatLon(proj,corrected_range,dop_buf[i][j],&x,&y);
	x=x-xoff;
        y=y-yoff;
        x=x*radtodeg;
        y=y*radtodeg; 
        // unwrap y if necessary
        if(y>180) y=-(360-y);
        else if(y<-180) y=360+y;
 
        // compute radar mode number
        int rmn;
	if(adc>Uvar(900000,"Hz") && adc<Uvar(1100000,"Hz")) rmn=1;
	else if(adc>Uvar(1900000,"Hz") && adc<Uvar(2100000,"Hz")) rmn=0;
	else if (adc>Uvar(4900000,"Hz") && adc<Uvar(5100000,"Hz")) rmn=-1;
	else rmn=-2;

        // detect backsactter value
        double s0=abs(sigma0[i][j]);
        s0*=s0;
	// estimate snr
        float snr=s0*Xfactor[i][j]/pixd_noise_factor;
        
        // print to file
        fprintf(pixdfp,"%g %g %g %g %d %d %g %g\n",x,y,elev_buf[i][j],azim_buf[i][j],looknumber,rmn,s0,snr);
      }
    } 
    // update look number
    looknumber++; 
}

void 
L1I::setBadSARData(){
  // Incomplete
  ErrorMessage e("SAR data is bad .....");
  e.throwMe();
}

void
L1I::fixRerouteChirpBurst(){
  if(csr!=4){
    cerr << "Warning: attempted to call L1I::fixRerouteChirpBurst on data with csr!=4" << endl;
    cerr << "Ignoring ..." << endl;
  }
  beam_number=3; // set a an actual beam number

  // fake a recieve window delay that adds the nominal boresight range to
  // the -2 ms delay actually used for rerouted chirp. This should be the rerouted
  // chirp "point target" in the center of the range window.

  Uvar origrwd=rwd;
  Uvar lastrwd=-1000;

  Uvar apparent_range;
  while(fabs(rwd-lastrwd)>Uvar(1e-6,"s")){ 
    lastrwd=rwd;
    Uvar receive_window=ctrx*pri;
    Uvar transmit_window=pul*pri;
    Uvar transmit_center=transmit_window/2;
    Uvar receive_center=rwd+receive_window/2;
    Uvar time_offset=(transmit_center+receive_center)/2;
    ground_impact_time=t+time_offset;
    cassini_sim.setBeamTime(beam_number,ground_impact_time);
    Uvar rtt=(cassini_sim.boresightRange()*2)/speed_light;
    rwd=origrwd+rtt;
    apparent_range=(rtt*speed_light)/2.0;
  }

    cout << "L1I:fixRerouteChirpBurst: the apparent range for the reroute chirp in SAB " << sab_counter << " should be " << apparent_range.getInUnits("km") << " km" << endl;

    cerr << "L1I:fixRerouteChirpBurst: the apparent range for the reroute chirp in SAB " << sab_counter << " should be " << apparent_range.getInUnits("km") << " km" << endl;

    // fake a chirp start frequency that put 0 Hz in middle of Frequency window
    Uvar fdop;
    Uvar fcsforig=fast_csf;
    Uvar lastfcsf=100e6;
    while(fabs(fast_csf-lastfcsf)> Uvar(0.1,"Hz")){
      lastfcsf=fast_csf;
      PositionVector dummy;
      Uvar chirp_freq=speed_light/cassini_sim.lambda(); // carrier
      chirp_freq+=fast_csf; // add chirp start freq
      chirp_freq+=slow_cfs*(csq-1)/2; // add half of chirp bandwidth
      Uvar lchirp=speed_light/chirp_freq; // convert to wavelength
  
      fdop=cassini_sim.getDoppler(lchirp);
      fast_csf=fcsforig-fdop;
    }   
    

    cout << "L1I:fixRerouteChirpBurst: the apparent Doppler for the reroute chirp in SAB " << sab_counter << " should be " << fdop.getInUnits("Hz") << " Hz" << endl;

    cerr << "L1I:fixRerouteChirpBurst: the apparent Doppler for the reroute chirp in SAB " << sab_counter << " should be " << fdop.getInUnits("Hz") << " Hz" << endl;

    // recompute footprint geometry using new rwd
    // SAR processor uses this to determine processing window in OblCyl

    Umat azim_1way3dB_ellipse_fit("azimuth ellipse fit",5,4);
    Umat elev_1way3dB_ellipse_fit("azimuth ellipse fit",5,4);
    Umat azim_2way3dB_ellipse_fit("azimuth ellipse fit",5,4);
    Umat elev_2way3dB_ellipse_fit("azimuth ellipse fit",5,4);
    azim_1way3dB_ellipse_fit=Uvar(0,"rad");
    elev_1way3dB_ellipse_fit=Uvar(0,"rad");
    azim_2way3dB_ellipse_fit=Uvar(0,"rad");
    elev_2way3dB_ellipse_fit=Uvar(0,"rad");
  
    for (unsigned int i = 0; i < 5;i++){
      Uvec azi1_fit("azi fit",4), elev1_fit(" ",4);    
      beam_[i].computeBestFitEllipse(azi1_fit,elev1_fit,-3.0);//one way 3dB  
      for(unsigned int j=0;j<4;++j){
	azim_1way3dB_ellipse_fit(i,j)= azi1_fit(j);
	elev_1way3dB_ellipse_fit(i,j)=elev1_fit(j);
      }    
      Uvec azi2_fit("azi fit",4), elev2_fit(" ",4);   
      beam_[i].computeBestFitEllipse(azi2_fit,elev2_fit,-1.5);//two way 3dB   
      for(unsigned int j=0;j<4;++j){
	azim_2way3dB_ellipse_fit(i,j)= azi2_fit(j);
	elev_2way3dB_ellipse_fit(i,j)=elev2_fit(j);
      }
    }
    
    // to make computeMeasurementGeometry work we also need to fudge csr 
    csr=0;
    computeMeasurementGeometry(azim_1way3dB_ellipse_fit,
		     elev_1way3dB_ellipse_fit,
		     azim_2way3dB_ellipse_fit,
		     elev_2way3dB_ellipse_fit);
  
}

// This routine is INCOMPLETE
// For now it only returns one value per burst but it should return one value
// per pixel computed from slon and slat
float 
L1I::getIncidenceAngle(float& slat, float& slon){

  //-------------------------------------------------------
  // compute incidence angle on surface of reference sphere
  //--------------------------------------------------------
  double pos[3], look[3], normal[3];
  // convert lat/lon to position in target body fixed coordinates on the
  // reference sphere (and unit surface normal vector)
  float theta=pi/2-slat;
  normal[0]= sin(theta)*cos(slon); 
  normal[1]= sin(theta)*sin(slon);
  normal[2]= cos(theta);
  for(int c=0;c<3;c++) pos[c]=normal[c]*target_radius_;

  // compute look unit vector; scpos_ was already computed during burst
  // processing
  double sum=0;
  for(int c=0;c<3;c++) {
    look[c]=scpos_[c]-pos[c];
    sum+=look[c]*look[c];
  }
  sum=sqrt(sum);
  for(int c=0;c<3;c++) {
    look[c]/=sum;
  }

  // compute incidence angle
  float inc=acos(dot(look,normal));
  if(inc<0 || inc> pi/2 +0.001){
    fprintf(stderr,"getIncidenceAngle::Weird incidence angle %g , look =[%g %g %g] normal=[%g %g %g]",
	    inc,(float)look[0],(float)look[1],(float)look[2],(float)normal[0],(float)normal[1],(float)normal[2]);
  }
  return(inc);
}


bool
L1I::initInterpolate(double fdop_in_Hz, double range_in_km){

  // get minimum range (range[0] is actually the midpoint of the first pixel)
  // so nearest indices are determined by truncation rather than rounding
  // no +0.5 in floor( ... )
  float range_min=range[0];
  float relrange=(range_in_km-range_min)/range_step_;
  current_range_idx_=int(floor(relrange));
    
  
  // check for range_index out of bounds
  if(current_range_idx_<0 || current_range_idx_> (int)num_range_bins-1) 
    return(false);
  
  double reldop=fdop_in_Hz-fdop_centroid[current_range_idx_];
  double mindop=-prf_in_Hz/2;
  current_doppler_idx_=(int)(floor((reldop-mindop)/doppler_step_));

  // check for doppler index out of bounds
  if(current_doppler_idx_<0 || current_doppler_idx_> Nfft_az-1) 
    return(false); 

  // check for zero Xfactor out of useable area case
  if(Xfactor[current_range_idx_][current_doppler_idx_]==0.0) return(false);

  if (use_sinc_interpolation_)
    {
      // check for range index out of bounds for sinc interpolation
      if (! sinc_table_range.inBounds(range_in_km, range_min, range_step_,
				      0, num_range_bins-1)) return(false);
      
      // linearly interpolate fdop_centroid to range_in_km and recalculate reldop
      if (current_range_idx_ < (int)num_range_bins-1)
	{
	  reldop=fdop_in_Hz-(fdop_centroid[current_range_idx_] +
			     (fdop_centroid[current_range_idx_+1]-fdop_centroid[current_range_idx_])*
			     (relrange-floor(relrange)));
	}
      // check for doppler index out of bounds for sinc interpolation
      if (! sinc_table_doppler.inBounds(reldop, mindop, doppler_step_,
					  0, Nfft_az-1)) return(false);
    }

  if(use_fft_interpolation_){
    if(!fft_interpd_) fftinterp();
    float fft_relrange=(range_in_km-fftinterp_range_min)/fftinterp_range_step;
    fft_range_idx_=int(floor(fft_relrange));
    
    // check for range_index out of bounds
    if(fft_range_idx_<0 || fft_range_idx_> (int)fftinterp_num_range_bins-1) 
      return(false);
   
    double fft_reldop=fdop_in_Hz-fftinterp_fdop_centroid[fft_range_idx_];
    fft_doppler_idx_=(int)(floor((fft_reldop-fftinterp_reldop_min)
				 /fftinterp_doppler_step));

    // check for doppler_index out of bounds
    if(fft_doppler_idx_<0 || fft_doppler_idx_> (int)fftinterp_num_doppler_bins-1) 
      return(false);
  }
  return(true);
}

void  
L1I::initAzimElevTransformation(){
  // If speedup is needed this routine will initialize
  // the taylor series expansions for dop/range to azim/elev transformations

  // for now do nothing

  // order of taylor series needs to be 
  // determined and coefficient arrays allocated
  // in config
  
  // computeAzimElevDerivatives();
  // compute Taylor coefficients

}

void 
L1I::computeSARAncillaryData(const SARProcParams& spp, OblCylProj& proj){

  // clear SAR Anciallary data invalid bit (bit 9)
  science_qual_flag=science_qual_flag & 0xfffffdff;

  // compute effective resolutions
  float chirp_bw=get_in_base_units(slow_cfs*(csq-1));
  float range_res=get_in_base_units(speed_light)/(2*chirp_bw);
  float doppler_res=1/(num_pulses_received*pri_in_s);


  // call pixel area routine to get range_ground_width_and doppler_ground_width_
  float doppler_in_hz=get_in_base_units(nominal_doppler_centroid);
  float range_in_km=get_in_base_units(nominal_delay*speed_light/2);
  float dummy_area;
  dummy_area=getPixelAreaOnSphere(doppler_in_hz,range_in_km);

  // theoretical effective resolution
  float doppler_ground_res=doppler_ground_width_ *
	      doppler_res/doppler_step_;
  float range_ground_res=range_ground_width_*range_res/range_step_;

  // multiply resolutions by empirical quantity from config file
  sar_azimuth_res=Uvar(doppler_ground_res,"km")*spp.azimuth_res_coeff;
  sar_range_res=Uvar(range_ground_res,"km")*spp.range_res_coeff;

  // compute BIDR coordinates
  double slat=get_in_base_units(act_centroid_lat);
  double slon=get_in_base_units(act_centroid_lon);
  double lat=proj.latInRad(slon,slat);
  double lon=proj.lonInRad(slon,slat);
  sar_centroid_bidr_lat=Uvar(lat,"rad");
  sar_centroid_bidr_lon=Uvar(lon,"rad");
}

// for now this does nothing if necessary it will be used to speed up
// transformation between projected lat lon and range/doppler coordinates
void  
L1I::initLatLonTransformation(const SARProcParams& spp, OblCylProj& proj){
  // for now we use this to dump pixel info if desired
  if(pixdfp!=NULL){
    pixelDump(proj);
  }
}

float L1I::sigma0Interpolate(){
  // get each weight using  SincTable::weight()
  // and compute weighted sum
  // and detect s0 complex==> real
  complex<float>s0;
  if(!use_fft_interpolation_)s0=sigma0[current_range_idx_][current_doppler_idx_];
  else s0=interp_sigma0[fft_range_idx_][fft_doppler_idx_];
  float detected_s0=abs(s0);
  detected_s0*=detected_s0;
  return(detected_s0);
}

float L1I::sigma0Interpolate(double fdop_in_Hz, double range_in_km){

  int k;
  complex<float> s0;

  if (! use_sinc_interpolation_)
  {

    if(!use_fft_interpolation_)s0=sigma0[current_range_idx_][current_doppler_idx_];

    // Use nearest neighbor instead
    else s0=interp_sigma0[fft_range_idx_][fft_doppler_idx_];

  }
  else
  {
    double range_min = range[0];
    double mindop = -prf_in_Hz/2;
    double reldop;

    // Note that only the elements of interp_range_buf between indices
    // min_range_idx_ and max_range_idx_ are used in the interpolations.
    min_range_idx_ = sinc_table_range.minIndex(range_in_km, range_min,
      range_step_, 0);
    max_range_idx_ = sinc_table_range.maxIndex(range_in_km, range_min,
      range_step_, num_range_bins-1);
    // linearly interpolate fdop_centroid to range_in_km and calculate reldop
    if (current_range_idx_ < (int)num_range_bins-1)
    {
      float relrange=(range_in_km-range_min)/range_step_;
      reldop = fdop_in_Hz-(fdop_centroid[current_range_idx_] +
	(fdop_centroid[current_range_idx_+1]-fdop_centroid[current_range_idx_])*
	(relrange-floor(relrange)));
    }
    else
    {
      reldop = fdop_in_Hz - fdop_centroid[current_range_idx_];
    }
    for (k = min_range_idx_; k <= max_range_idx_; k++)
    {
      interp_range_buf[k] = sinc_table_doppler.interpolate(reldop,
	mindop, doppler_step_, 0, Nfft_az-1, &sigma0[k][0]);
    }
    // Interpolate interp_range_buf to range_in_km.
    s0 = sinc_table_range.interpolate(range_in_km, range_min,
      range_step_, 0, num_range_bins-1, interp_range_buf);
  }

  float detected_s0=abs(s0);
  detected_s0*=detected_s0;
  return(detected_s0);
}
bool L1I::nadirAmb(double fdop_in_Hz, double range_in_km){
  // account for doppler induced error in range
  // converts true range value to range compression range estimate
  double nadir_r=nadir_range_in_km-range_doppler_correction_coeff*(nadir_fdop_in_Hz-nominal_doppler_centroid_in_Hz);
  float dd=fmod(fabs(nadir_fdop_in_Hz-fdop_in_Hz),prf_in_Hz);
  float dr=fmod(fabs(nadir_r-range_in_km),range_window_in_km);
  bool retval=false;
  if(dd<10*doppler_step_ && dr < 10*range_step_) retval=true;
  return(retval);
}
float L1I::XInterpolate(){
  // get each weight using  SincTable::weight()
  // and compute weighted sum
  float X;
  if (! use_sinc_interpolation_)
  {
    X=Xfactor[current_range_idx_][current_doppler_idx_];
  }
  else
  {
    // INCOMPLETE -- still uses nearest neighbor
    X=Xfactor[current_range_idx_][current_doppler_idx_];
  }
  return(X);
}


// routine to compute bounds for each burst in Latitude and Longitude
// for now extrapolate 3 dB two-way axes by 1.5 to get latlon box
void L1I::reportBoundary(const OblCylProj& proj){
  double lat[5];
  double lon[5];
  double slat[5];
  double slon[5];

  // centroid
  slat[0]=get_in_base_units(act_centroid_lat);
  slon[0]=get_in_base_units(act_centroid_lon);
  
  // ellipse points
  slat[1]=get_in_base_units(act_ellipse_pt1_lat);
  slon[1]=get_in_base_units(act_ellipse_pt1_lon);

  slat[2]=get_in_base_units(act_ellipse_pt2_lat);
  slon[2]=get_in_base_units(act_ellipse_pt2_lon);

  slat[3]=get_in_base_units(act_ellipse_pt3_lat);
  slon[3]=get_in_base_units(act_ellipse_pt3_lon);

  slat[4]=get_in_base_units(act_ellipse_pt4_lat);
  slon[4]=get_in_base_units(act_ellipse_pt4_lon);
  
  // compute ellipse points projected into OblCyl
  cerr << "L1I::boundaryReport (act_ellipse points and boundary in OblCyl"
       << endl;
  for(int i=0;i<5;i++){
    lat[i]=proj.latInRad(slon[i],slat[i]);
    lon[i]=proj.lonInRad(slon[i],slat[i]);

    if(i==0) cerr << "centroid_lat="<< lat[i]*radtodeg << ",lon=" <<
	       lon[i]*radtodeg << endl;
    else cerr << "ellipse pt#" << i <<" lat="<<lat[i]*radtodeg<< ",lon=" <<
	   lon[i]*radtodeg << endl;
  }  
  cerr << "lonmin=" << lonlat_bounds[0]*radtodeg << ", lonmax=" << lonlat_bounds[1]*radtodeg << endl;
  cerr << "latmin=" << lonlat_bounds[2]*radtodeg << ", latmax=" << lonlat_bounds[3]*radtodeg << endl;
}

// routine to compute bounds for each burst in Latitude and Longitude
// for now extrapolate 3 dB two-way axes by 1.5 to get latlon box
void L1I::computeLatLonBounds(const OblCylProj& proj){
  double lat[5];
  double lon[5];
  double dlat[5];
  double dlon[5];
  double slat[5];
  double slon[5];

  // if footprint is off the limb the routine defaults to the entire
  // sphere
  if(science_qual_flag & 0x0180){
    lonlat_bounds[0]=-pi; // LONMIN
    lonlat_bounds[1]=+pi; // LONMAX
    lonlat_bounds[2]=-pi/2; // LATMIN
    lonlat_bounds[3]=+pi/2; // LATMAX
    return;
  }

  // centroid
  slat[0]=get_in_base_units(act_centroid_lat);
  slon[0]=get_in_base_units(act_centroid_lon);
  
  // ellipse points
  slat[1]=get_in_base_units(act_ellipse_pt1_lat);
  slon[1]=get_in_base_units(act_ellipse_pt1_lon);

  slat[2]=get_in_base_units(act_ellipse_pt2_lat);
  slon[2]=get_in_base_units(act_ellipse_pt2_lon);

  slat[3]=get_in_base_units(act_ellipse_pt3_lat);
  slon[3]=get_in_base_units(act_ellipse_pt3_lon);

  slat[4]=get_in_base_units(act_ellipse_pt4_lat);
  slon[4]=get_in_base_units(act_ellipse_pt4_lon);
  
  // compute ellipse points projected into OblCyl
  for(int i=0;i<5;i++){
    lat[i]=proj.latInRad(slon[i],slat[i]);
    lon[i]=proj.lonInRad(slon[i],slat[i]);
  }  

  // compute corners from ellipse points (fails close to pole but that's OK
  // once oblique cylindrical is set up )

  // compute delta lat,lons from center to ellipse points
  for(int i=1;i<5;i++){
    dlat[i]=lat[0]-lat[i];
    dlon[i]=lon[0]-lon[i];
    if(dlon[i]>pi) dlon[i]=dlon[i]-2*pi;
    if(dlon[i]<-pi) dlon[i]=2*pi+dlon[i];
  }

  // use delta to compute parallelogram corners (store in lat lon)
  int k=0;
  for(int i=1;i<3;i++){
    for(int j=3;j<5;j++){
      k++;
      lat[k]=lat[0]+dlat[i]+dlat[j];
      lon[k]=lon[0]+dlon[i]+dlon[j];
    }
  }
  // get max lat and lon distances from boresight for corners
  float maxdlat=0.0;
  float maxdlon=0.0;
  for(int i=0;i<4;i++){
    float dlat=fabs(lat[0]-lat[i+1]);
    float dlon=fabs(lon[0]-lon[i+1]);
    if(dlon>pi) dlon=2*pi-dlon;
    if(dlat>maxdlat) maxdlat=dlat;
    if(dlon>maxdlon) maxdlon=dlon;
  }

  // set up bounds with 2X margin over 3 dB two-way footprint

  float BURST_PAD_FACTOR=4.0;
  if(adc<Uvar(900000,"Hz") || adc > Uvar(4000000,"Hz")) BURST_PAD_FACTOR = 4.0;
  // gives scat mode high sar more real estate to work with

  lonlat_bounds[0]=lon[0]-BURST_PAD_FACTOR*maxdlon; // LONMIN
  lonlat_bounds[1]=lon[0]+BURST_PAD_FACTOR*maxdlon; // LONMAX
  lonlat_bounds[2]=lat[0]-BURST_PAD_FACTOR*maxdlat; // LATMIN
  lonlat_bounds[3]=lat[0]+BURST_PAD_FACTOR*maxdlat; // LATMAX

}

bool
L1I::goodSideOfNadir(double slon, double slat){

  // compute x,y,z coordinates of pixel position (slon,slat)
  double pos[3], normal[3];
  double theta=pi/2-slat;
  normal[0]= sin(theta)*cos(slon); 
  normal[1]= sin(theta)*sin(slon);
  normal[2]= cos(theta);
  for(int c=0;c<3;c++) pos[c]=normal[c]*target_radius_;
 
  return(goodSideOfNadir(pos,0));
}

bool
L1I::goodSideOfNadir(double pos[3], double eang){
  double dbore[3], dpos[3], mnorm[3];
  // compute x,y,z coordinates of nadir (assumes sphere).
  double dboremag=0,dposmag=0,mnormmag=0;


  // take the cross product of the velocity vector and the s/c position vector to get the normal
  // to the mirror plane
  mnorm[0]=scpos_[1]*scvel_[2]-scpos_[2]*scvel_[1];
  mnorm[1]=scpos_[2]*scvel_[0]-scpos_[0]*scvel_[2];
  mnorm[2]=scpos_[0]*scvel_[1]-scpos_[1]*scvel_[0];
  
  // compute units vectors to boresight and pos from spacecraft and normalize mnorm unit vector
  for(int c=0;c<3;c++){
    dpos[c]=pos[c]-scpos_[c];
    dbore[c]=borepos_[c]-scpos_[c];
    dposmag+=dpos[c]*dpos[c];
    dboremag+=dbore[c]*dbore[c];
    mnormmag+=mnorm[c]*mnorm[c];
  }
  mnormmag=sqrt(mnormmag);
  dposmag=sqrt(dposmag);
  dboremag=sqrt(dboremag);

  for(int c=0;c<3;c++){
    mnorm[c]/=mnormmag;
    dpos[c]/=dposmag;
    dbore[c]/=dboremag;
  }

  // compute dot products of dbore and dpos with mnorm
  double dotbore=0, dotpos=0;
  for(int c=0;c<3;c++){
    dotbore+=dbore[c]*mnorm[c];
    dotpos+=dpos[c]*mnorm[c];
  }

  // check to see if signs of dot products are the same
  double ind=dotbore*dotpos; // negative is a mirror ambiguity, positive is not

  
  if(ind>0 & eang==0){ return(true);} // return true if not a mirror ambiguities and there no exclusion angle about nadir
  else if(ind>0){ // if not not a mirror ambiguity
    float ang=fabs(asin(dotpos));
    if(ang<eang) return(false); // return false if within eang degrees of mirror plane
    else return(true); // return true if outside eang degrees of mirror plane
  }
  // return false for mirror ambiguities
  else return(false);  
}

// convert from Planetocentric LatLon to Doppler and Range.
// Three modes 
// Nominal (Titan SAR) mode assumes planet is a sphere.
// Orthorectification mode (Uses a topomap to get the surface height
// which is looked up by Planetocentric lat lon
// Triaxial Mode (Uses the triaxial shape for body) when HISARcalibrate is on

void
L1I::latLonToDopplerAndRange(double slon, 
			     double slat, 
			     double& fdop_in_Hz,
			     double& range_in_km)
{
 
  /**** 
 static TargetGeom tg;
 tg.reset(ground_impact_time);
 tg.setTarget();
 tg.setState(scstate_); // sc_state_ computed in computeBurstParams
 tg.setLatLon(Uvar(slat,"rad"),Uvar(slon,"rad"));
 fdop_in_Hz=get_in_base_units(tg.doppler(lambda_chirp));
 range_in_km=get_in_base_units(tg.range());
  *****/
 

  //-------------------------------------------------------
  // compute doppler and range from coordinates of reference sphere
  //--------------------------------------------------------
  
  double pos[3], look[3], normal[3];

  // convert lat/lon to position in target body fixed coordinates on the
  // reference sphere (and unit surface normal vector)
  double theta=pi/2-slat;
  normal[0]= sin(theta)*cos(slon); 
  normal[1]= sin(theta)*sin(slon);
  normal[2]= cos(theta);
  if(!orthorectify && !use_HISAR_calibrate){
    for(int c=0;c<3;c++) pos[c]=normal[c]*target_radius_;
  }
  // triaxial case
  else if(use_HISAR_calibrate){
    double center[3]={0,0,0};
    double dummyr;
    if(! get_surface_intercept_triaxial(center,normal,default_target_radii_in_km,pos,dummyr)){
	cerr << "Weird error in L1I::latLonToDopplerRange that shouldn't happen" 
	<< endl;
	exit(2);
      }
  }
  // orthorectification case
  else{
    bool topo_inbounds=false;
    float orthoheight=ortho_tmap.getInterpolatedHeightInKm(slat,slon,topo_inbounds);
    if(!topo_inbounds){
      fprintf(stderr,"L1I::latLonToDopplerAndRange Bad TopoMap for Orthorectification (need to cover whole imaging area)\n");
      exit(1);
    }
    for(int c=0;c<3;c++) pos[c]=normal[c]*(target_radius_+orthoheight);
    //for(int c=0;c<3;c++) pos[c]=normal[c]*(target_radius_+orthoheight);
  } // end orthrectification
  // compute look unit vector; scpos_ was already computed during burst
  // processing
  double sum=0;
  for(int c=0;c<3;c++) {
    look[c]=pos[c]-scpos_[c];
    sum+=look[c]*look[c];
  }
  range_in_km=sqrt(sum);

  for(int c=0;c<3;c++) {
    look[c]/=range_in_km;
  }

  // compute dot product of look vector and s/c velocity vector for
  // use in doppler computation;
  double vrange_bore=dot(look,scvel_);

  // Classical doppler
  // fdop_in_Hz=2.0*vrange_bore/lambda_chirp_in_km;

  //special_relativity
  fdop_in_Hz=2.0*vrange_bore*chirp_freq_in_Hz/
    (speed_light_in_kps-vrange_bore);
 
  // account for doppler induced error in range
  // converts true range value to range compression range estimate
  range_in_km-=range_doppler_correction_coeff*
    (fdop_in_Hz-nominal_doppler_centroid_in_Hz);
    

}
bool L1I::dopplerRangeHeightToLatLon(double fdop_in_Hz,
				     double range_in_km, double height_in_km,double& slat, double& slon){
  double u[3], u1[3], u2[3], pos[3];

  // get two possible looks in TBF coordinates
  // most parameters to this function were computed in 
  // initGeometryCalculations.
  // If no such look exists return 0 gain
  if(! fast_doppler_range_to_TBF_look(fdop_in_Hz,range_in_km,
	   target_radius_+height_in_km, scpos_, scvel_, lambda_special, u1,u2)){
    return(false);
  }
 

  // pick the look closest in angle to the boresight
  // compares dot products (highest wins)
  // boresight in TBF was precomputed in initGeometryCalculations.
  // If the SAR image intersects nadar this procedure will fail because a
  // particular doppler/range bin will contain energy from BOTH mirror looks.
  // Of course the one closest to the boresight will have the larger gain.
  double cosangle1=dot(boresight_in_tbf_,u1);
  double cosangle2=dot(boresight_in_tbf_,u2);
 
  double r=target_radius_+height_in_km;
  for(int c=0;c<3;c++){
    if(cosangle1>cosangle2) u[c]=u1[c];
    else u[c]=u2[c];

    pos[c]=scpos_[c]+range_in_km*u[c];
    pos[c]=pos[c]/r;
  }
  slat=asin(pos[2]);
  slon=atan2(pos[1],pos[0]);
  
  return(true);
}
// BAQ simulation used in noise only processing
void L1I::BAQEncodeAndDecode(RasList* raslist){
  Uvar tro_pri = tro/pri;
  int tro_int=int(floor(tro_pri.getInUnits("")+0.5));
  baq.setParams(pri,pul,baq_mode,adc,tro_int);
  Ivec Thresh("Thresh");
  Charvec tmp("tmp",Nradar_data);
  Fvec tmp_float("tmp_float",Nradar_data);
  for(unsigned int c=0;c<Nradar_data;c++) 
    tmp(c)=(char)radar_data[c];
  
  // Get Thresholds from rasfile
  raslist->findRecord(sab_counter);
  raslist->decodeBaqThresholds();

  Charvec words("words_name");
  Thresh.resize(raslist->baq_decoded.size());
  for(unsigned int c=0;c<raslist->baq_decoded.size();c++)
    Thresh(c)=raslist->baq_decoded(c);
  
  baq.Encode_Nbit(tmp,Thresh, words);
  
  baq.Decode_Nbit(words, Thresh , tmp_float);
  for(unsigned int c=0;c<Nradar_data;c++) 
    radar_data[c]=tmp_float(c);  
}

float L1I::minimumBAQThresholdVariance(RasList* raslist){
  Uvar tro_pri = tro/pri;
  int tro_int=int(floor(tro_pri.getInUnits("")+0.5));
  baq.setParams(pri,pul,baq_mode,adc,tro_int);
  Ivec Thresh("Thresh");
  // Get Thresholds from rasfile
  raslist->findRecord(sab_counter);
  raslist->decodeBaqThresholds();
  Thresh.resize(raslist->baq_decoded.size());
  for(unsigned int c=0;c<raslist->baq_decoded.size();c++)
    Thresh(c)=raslist->baq_decoded(c);
  return(baq.minVariance(Thresh));
}
//-------------------------------------------------------------------------
// routine for computing two-way antenna gain from doppler and range
//-------------------------------------------------------------------------
float
L1I::computeGain2(double doppler_in_hz, double range_in_km, 
		  const DebugInfo& dbg, float& azim_deg, float& elev_deg,
		  double height_in_km)
{
  //---------------------------------------------------------------
  // fast gain computation
  //---------------------------------------------------------------
  double u[3], u1[3], u2[3];

  // get two possible looks in TBF coordinates
  // most parameters to this function were computed in 
  // initGeometryCalculations.
  // If no such look exists retrun 0 gain

  // utilize approximate but fast special relativity correction
  
  if(! fast_doppler_range_to_TBF_look(doppler_in_hz,range_in_km,
	   target_radius_+height_in_km, scpos_, scvel_, lambda_special, u1,u2)){
    return(0.0);
  }
 

  // pick the look closest in angle to the boresight
  // compares dot products (highest wins)
  // boresight in TBF was precomputed in initGeometryCalculations.
  // If the SAR image intersects nadar this procedure will fail because a
  // particular doppler/range bin will contain energy from BOTH mirror looks.
  // Of course the one closest to the boresight will have the larger gain.
  double cosangle1=dot(boresight_in_tbf_,u1);
  double cosangle2=dot(boresight_in_tbf_,u2);
  for(int c=0;c<3;c++){
    if(cosangle1>cosangle2) u[c]=u1[c];
    else u[c]=u2[c];
  }

  // Transform look into beamframe
  // rotmat_ was precomputed in initGeometryCalculations()
  double u_bf[3]={0,0,0};
  rotate_vector(u,rotmat_,u_bf);

  //Get azimuth and elevation in radians
  double elev=asin(u_bf[1]);
  double azim=asin(u_bf[0]/cos(elev));

  
  //Compute two-way gain
  // beam was assigned in computeBurstParams()
  azim_deg = (float) azim*radtodeg;
  elev_deg = (float) elev*radtodeg;

  double gain2=beam_[beam_number-1].bilinear(azim,elev);
  gain2*=gain2; // approximate two-way gain as square of one-way gain at 
                // ground impact time

  //-------------------------------------------------------------------
  // For non-zero debug levels re-compute gain the slow way and compare
  //-------------------------------------------------------------------

  if(dbg.level){

    // construct TargetGeom using doppler and range and get look vector
    // A unique look vector is determined by 
    // TargetGeom::setRangeDopplerInTargetFrame
    // uses ground_impact_time and lambda_chirp computed in 
    // computeBurstParameters()
    TargetGeom tg;
    tg.reset(ground_impact_time);
    tg.setTarget();
    tg.setState(scstate_); // sc_state_ computed in computeBurstParams
    tg.setRangeDopplerInTargetFrame(Uvar(range_in_km,"km"),
				    Uvar(doppler_in_hz,"Hz"),lambda_chirp);
    DirectionVector look=tg.lookDirection();


    // Represent Look In Beam Frame

    look.representIn(Frame(beam_frame_spice_id[beam_number-1],
			   cassini_spice_id));


    // Get azimuth and elevation
    Uvar azim_uv,elev_uv;
    look.getAzimuthElevation(azim_uv,elev_uv);

    // Compute two-way gain
    double gain_check2=beam_[beam_number-1].bilinear(azim_uv,elev_uv);
    gain_check2*=gain_check2;

    double percent_error=100*(gain_check2-gain2)/gain_check2;

    // compute debugging values
    double u_chk[3], ubf_chk[3],cosangle_chk;
    ubf_chk[0]=look[DirectionVector::X];
    ubf_chk[1]=look[DirectionVector::Y];
    ubf_chk[2]=look[DirectionVector::Z];
    look.representIn(Frame(default_target_frame_spice_id,
			   default_target_spice_id));
    u_chk[0]=look[DirectionVector::X];
    u_chk[1]=look[DirectionVector::Y];
    u_chk[2]=look[DirectionVector::Z];
    cosangle_chk=dot(boresight_in_tbf_,u_chk);
    

    int debug_level=dbg.level;
    if (debug_level>3) debug_level=3; // allows use of switch statement 
                                      // modify if additional levels are
                                      // desired
    
    switch(debug_level){
    case 3:   
      dbg.file << "----------------------------------------------------"
	       << endl;
      dbg.file << "L1I::ComputeGain2: azim(deg):" << azim*180/pi 
	       << " elev(deg):" << elev*180/pi 
	       << " cosangle1:" << cosangle1
	       << " cosangle2:" << cosangle2 << endl;
      dbg.file << "L1I::ComputeGain2: azim_chk(deg):" 
	       << azim_uv.getInUnits("deg")
	       << " elev_chk(deg):" 
	       << elev_uv.getInUnits("deg")
	       << " cosangle_chk:" << cosangle_chk << endl;
      dbg.file << "L1I::ComputeGain2: ulook_target_body_fixed(x,y,z):"
	       << "(" << u[0] <<"," << u[1] << "," << u[2] << ")" << endl;
      dbg.file << "L1I::ComputeGain2: Check ulook_target_body_fixed(x,y,z):"
	       << "(" << u_chk[0] <<"," << u_chk[1] << "," << u_chk[2] << ")" 
	       << endl;
      dbg.file << "L1I::ComputeGain2: ulook_beam_frame(x,y,z):"
	       << "(" << u_bf[0] <<"," << u_bf[1] << "," << u_bf[2] << ")" 
	       << endl;
      dbg.file << "L1I::ComputeGain2: Check ulook_beam_frame(x,y,z):"
	       << "(" << ubf_chk[0] <<"," << ubf_chk[1] << "," << ubf_chk[2] 
	       << ")" 
	       << endl;
      dbg.file << "L1I::ComputeGain2: boresight_target_body_fixed(x,y,z):"
	       << "(" << boresight_in_tbf_[0] <<"," 
	       << boresight_in_tbf_[1] << "," << boresight_in_tbf_[2] << ")" 
	       << endl;
    case 2:
      dbg.file << "L1I::ComputeGain2: Range:" << range_in_km << " Doppler:"
	       << doppler_in_hz << " GainPercentError:" << percent_error
	       << " Gain " << gain2 << " GainCheck" << gain_check2
	       << endl;
 
    case 1:
      if(percent_error >1){
      dbg.file << "Warning ! percent gainerror exceeded in L1I::ComputeGain2:"
	       << endl << "      Range:" << range_in_km << " Doppler:"
	       << doppler_in_hz << " GainPercentError:" << percent_error
	       << endl;
      }
      break;
    default:
      break;
    }
  }

  return(gain2);
}


//--------------------------------
// Routine returns area of a bin in azimuth and elevation in square km
// inputs are center in azim and elev (in Beam Frame) and resolution in each dimension (all in radians)
//--------------------------------
float L1I::getBinAreaOnTriaxialBody(double radii[3], double azim0, double azim_res, double elev0, double elev_res){
  // This routine uses the two triangle surface area approximation method.
  
  // First we compute the location of the four corners of the pixel.
  // This part will not operate correctly if the doppler/range image 
  // intersects the
  // ground track.
  float area=0;
  
  double u[4][3], ubeam[3];
  double pos[4][3];
  int k=0;
  for(int i=0;i<2;i++){
    double a=azim0+(i-0.5)*azim_res;
    for(int j=0;j<2;j++){
      double e=elev0+(j-0.5)*elev_res;

      // compute look in target body fixed for each corner
      ubeam[1]=sin(e);
      ubeam[0]=sin(a)*cos(e);
      ubeam[2]=sqrt(1-ubeam[1]*ubeam[1]-ubeam[0]*ubeam[0]);

      rotate_vector(ubeam,rotmat_rev_,u[k]);
      k++;
    }
  }

 
  // get surface points from looks
  double r;
  for(int c=0;c<4;c++){
    double found=get_surface_intercept_triaxial(scpos_,u[c],radii,pos[c],r);
    if(!found) return(0); // returns zero area for bin if a corner is off the body
  }
  
 
  // Then we compute the area of the two triangles(assumes surface area is approximately planar)

  area+=getFlatTriangleArea(pos[0],pos[1],pos[2]);
  area+=getFlatTriangleArea(pos[1],pos[2],pos[3]);
  return(area);
}

// returns pixel area in square km
float L1I::getPixelAreaOnSphere(double doppler, double range, double height_in_km){
  // This routine uses the two triangle surface area approximation method.
  
  // First we compute the location of the four corners of the pixel.
  // This part will not operate correctly if the doppler/range image 
  // intersects the
  // ground track.
  float area=0;
  
  double u[4][3], u1[3], u2[3];
  int corner_idx=0;
  for(int i=0;i<2;i++){
    double d=doppler+(i-0.5)*doppler_step_;
    for(int j=0;j<2;j++){
      double r=range+(j-0.5)*range_step_;

      // compute mirror looks to each corner point in RangeDoppler space
      // If looks are not found set area to zero.
      if(!fast_doppler_range_to_TBF_look(d,r,target_radius_+height_in_km,
	scpos_, scvel_, lambda_special, u1,u2)){
	if(DebugInfo::allWarnings){
	  cerr << "Warning L1I::getPixelAreaOnSphere: Pixel not on surface!!!";
	}
	return(0.0);
      }

      // pick the look closest in angle to the boresight
      // compares dot products (highest wins)
      // boresight in TBF was precomputed in initGeometryCalculations.
      // If the SAR image intersects nadar this procedure will fail because a
      // particular doppler/range bin will contain energy from BOTH mirror looks.
      // Of course the one closest to the boresight will have the larger gain.
      double cosangle1=dot(boresight_in_tbf_,u1);
      double cosangle2=dot(boresight_in_tbf_,u2);
      for(int c=0;c<3;c++){
        // pick look direction among mirrors
	if(cosangle1>cosangle2) u[corner_idx][c]=u1[c];
	else u[corner_idx][c]=u2[c];

	// get surface point from lookdir
	u[corner_idx][c]*=r;
	u[corner_idx][c]+=scpos_[c];
      }
      corner_idx++;
    }
  }

  // Then we compute the area of the two triangles inscribed on the sphere and
  // sum.
  double s01,s12,s20,s23,s31;
  area+=getSphericalTriangleArea(u[0],u[1],u[2],target_radius_+height_in_km,
				 s01,s12,s20);
  area+=getSphericalTriangleArea(u[1],u[2],u[3],target_radius_+height_in_km,
				 s12,s23,s31);
  doppler_ground_width_=(s20+s31)/2;
  range_ground_width_=(s01+s23)/2;
  return(area);
}


// returns pointers to OblCyl Lat and Lon

void L1I::getOblCylLatLon(OblCylProj& proj, double range, double doppler, double* lat, double* lon){
  // First compute standard lat and lon
  // compute mirror looks to each corner point in RangeDoppler space
  // If looks are not found set area to zero.
  double u1[3],u2[3];
  if(!fast_doppler_range_to_TBF_look(doppler,range,target_radius_,
	scpos_, scvel_, lambda_special, u1,u2)){
    cerr << "Error L1I::getOblCylLatLon: Point not on surface!!!";
    exit(0);
  }
  
  // pick the look closest in angle to the boresight
  // compares dot products (highest wins)
  // boresight in TBF was precomputed in initGeometryCalculations.
  // If the SAR image intersects nadar this procedure will fail because a
  // particular doppler/range bin will contain energy from BOTH mirror looks.
  // Of course the one closest to the boresight will have the larger gain.
  double cosangle1=dot(boresight_in_tbf_,u1);
  double cosangle2=dot(boresight_in_tbf_,u2);
  for(int c=0;c<3;c++){
    // pick look direction among mirrors
    if(cosangle1<cosangle2) u1[c]=u2[c];
    
    // get surface point (unit vector) from lookdir
    u1[c]*=range;
    u1[c]+=scpos_[c];
    u1[c]/=target_radius_;
  }

  //Compute standard lat and lon from xyz
  double theta,slon,slat;
  theta=acos(u1[2]);
  slon=atan2(u1[1],u1[0]);
  slat=pi/2-theta;

  // Transform to Oblique cylindrical
  
  *lat= proj.latInRad(slon,slat);
  *lon= proj.lonInRad(slon,slat);

  /***** DEBUG_TOOLS
  static int call_no=0;

  printf("%g %g %g %g %g %g %g %g %g\n",u1[0],u1[1],u2[2],range,doppler,slon,slat,*lon,*lat);
  fflush(stdout);
  if(call_no) exit(0);
  call_no++;
  ****/
}

//---------------------------------------------------
// routine for computing intermediate values needed by
// computeGain2 and computePixelAreaOnSphere
// getFastDopplerCentroid and getFastDopplerRate
// getFastLookVector
// This routine is only called once per burst in order
// to enhance computational efficiency.
//----------------------------------------------------

 void L1I::initGeometryCalculations()
{
  //---------------------------------------------------
  // perform geometry initializations
  //---------------------------------------------------
  DebugInfo dbg("L1I::initGeometryCalculations");
  if(dbg.level){
     dbg.file << "Start Debugging L1I::initGeometryCalculations Burst #"
	      << sab_counter << endl;
 }

  // #define BORESIGHT_FROM_INTGEOMFILE
  #ifndef BORESIGHT_FROM_INTGEOMFILE
  // compute boresight ground location in TBF using preprocessor data

    cassini_sim.boresightPositionInKm(borepos_);


 
  // Radius of target body
  // for sphere tangent to boresight centered at body center

  target_radius_=sqrt(borepos_[0]*borepos_[0]+borepos_[1]*borepos_[1]+
		      borepos_[2]*borepos_[2]);

  // Compute state vector -- trying to get rid of CassiniSim dependency
  Frame f(default_target_frame_spice_id,default_target_spice_id);
  ground_impact_time.getEt(ground_impact_time_in_s);

  f.ephemeris(scstate_,cassini_spice_id,ground_impact_time_in_s);
  
  // Estimate S/C acceleration using +- .01 s time interval
  // SLOW_ALERT: Possible minor time bottleneck
  f.ephemeris(scstate_tmp1_,cassini_spice_id,ground_impact_time_in_s+0.01);
  f.ephemeris(scstate_tmp2_,cassini_spice_id,ground_impact_time_in_s-0.01);
  scacc_vec_=(scstate_tmp1_.velocity()-scstate_tmp2_.velocity())
    /Uvar(0.02,"s");
  scacc_[0]=get_in_base_units(scacc_vec_[FloatVector::X]);
  scacc_[1]=get_in_base_units(scacc_vec_[FloatVector::Y]);
  scacc_[2]=get_in_base_units(scacc_vec_[FloatVector::Z]);

  // Spacecraft position 
  PositionVector pos=scstate_.position();
  scpos_[0]=pos.km(PositionVector::X);
  scpos_[1]=pos.km(PositionVector::Y);
  scpos_[2]=pos.km(PositionVector::Z);

  // Spacecraft velocity
  FloatVector vel=scstate_.velocity();
  scvel_[0]=get_in_base_units(vel[FloatVector::X]);
  scvel_[1]=get_in_base_units(vel[FloatVector::Y]);
  scvel_[2]=get_in_base_units(vel[FloatVector::Z]);  
  scspeed_=sqrt(scvel_[0]*scvel_[0]+scvel_[1]*scvel_[1]+scvel_[2]*scvel_[2]);



  // Boresight of beam in Target Body Fixed (TBF)
  Frame bf(beam_frame_spice_id[beam_number-1],cassini_spice_id);
  Frame tbf(default_target_frame_spice_id,default_target_spice_id);
  DirectionVector b(bf,ground_impact_time_in_s,0,0,1);
  b.representIn(tbf);
  boresight_in_tbf_[0]=b[DirectionVector::X];
  boresight_in_tbf_[1]=b[DirectionVector::Y];
  boresight_in_tbf_[2]=b[DirectionVector::Z];

  // compute lambda_chirp and range to boresight;
  lambda_chirp_in_km=get_in_base_units(lambda_chirp);
  chirp_freq_in_Hz=get_in_base_units(speed_light/lambda_chirp);

  double vr=(scvel_[0]*boresight_in_tbf_[0]+scvel_[1]*boresight_in_tbf_[1]+scvel_[2]*boresight_in_tbf_[2]);
  // special lambda value used to approximate the effects of special relativity in fast computations 
  // used for calls to fast_doppler_range_to_TBF_look
  lambda_special=lambda_chirp_in_km*speed_light_in_kps/(speed_light_in_kps - vr);

  double range_in_km=sqrt((scpos_[0]-borepos_[0])*(scpos_[0]-borepos_[0])+
    (scpos_[1]-borepos_[1])*(scpos_[1]-borepos_[1])+
    (scpos_[2]-borepos_[2])*(scpos_[2]-borepos_[2]));

  nominal_range=Uvar(range_in_km,"km");
  if(nadir_mode) nominal_range=cassini_sim.altitude();
  nominal_delay=nominal_range*2/speed_light;

  // compute rotation matrix from TBF to beam frame
  tbf.rotationMatrix(tbf,bf,ground_impact_time_in_s,rotmat_);

  // compute rotation matrix from beam frame to TBF
  bf.rotationMatrix(bf,tbf,ground_impact_time_in_s,rotmat_rev_);


  // assign nominal_doppler_centroid to doppler at nominal_range and
  // the appropriate azimuth offset
  // assumes boresight is on the surface
  double look[3];
  // If you cannot find the appropriate azimuth offset (unusual case
  // in which elevation is nearly perpendicular to range
  if(!getFastLookVector(range_in_km,beam_azioff_rad[beam_number-1],look)){  
    for(int c=0;c<3;c++) {
      look[c]=boresight_in_tbf_[c];
    }
  }
  nominal_doppler_centroid=Uvar(getFastDopplerCentroid(look),"Hz");
#endif
  // NEW NEW NEW test test test
  // let's do this using the intermediate geometry data file only
  // ASSUMES TARGET IS SPHERICAL !!!!!!!!
  // We need rotmat rotmat_rev_ look nominal_doppler_centroid
  // borepos_ target_radius_ scstate_ scacc_ scpos_ scvel_
  // boresight_in_tbf_ lambda_chirp_in_km  nominal_range nominal_delay
#ifdef BORESIGHT_FROM_INTGEOMFILE
  // set ground_impact_time
  ground_impact_time.getEt(ground_impact_time_in_s);
  
  // directly access geometry values from Frame object globals
  target_radius_= Frame::directAccessTargetRadius(ground_impact_time);
  Frame::directAccessRotMatB2T(ground_impact_time,beam_number,rotmat_rev_);
  Frame::directAccessRotMatT2B(ground_impact_time,beam_number,rotmat_);
  Frame::directAccessSCPos(ground_impact_time,scpos_);
  Frame::directAccessSCVel(ground_impact_time,scvel_);
  scspeed_=sqrt(scvel_[0]*scvel_[0]+scvel_[1]*scvel_[1]+scvel_[2]*scvel_[2]);
  double posmag=sqrt(scpos_[0]*scpos_[0]+scpos_[1]*scpos_[1]
		     +scpos_[2]*scpos_[2]);
  Frame::directAccessSCAcc(ground_impact_time,scacc_);

  // compute lambda_chirp 
  lambda_chirp_in_km=get_in_base_units(lambda_chirp);

  // compute boresight lookvector
  double look_inbeamframe[3]={0,0,1};
  mxv_c(rotmat_rev_,look_inbeamframe,boresight_in_tbf_);
  // here we assume a sphere  to compute the range
  double B=2*dot(boresight_in_tbf_,scpos_);
  double C=posmag*posmag - target_radius_*target_radius_;
  double check=B*B-4*C;
  double range_in_km;
  if(check<0) range_in_km=0;
  else{
    double r1 =(-B+sqrt(check))/2;
    double r2 =(-B-sqrt(check))/2;
    range_in_km=MIN(r1,r2);
  }
  if(range_in_km<=0){
    cerr << "Warning BAD Boresight Range " << range_in_km << " " 
	 << "SAB:" << sab_counter << endl;
  }
  // compute boresight ground position
  borepos_[0]=scpos_[0]+range_in_km*boresight_in_tbf_[0];
  borepos_[1]=scpos_[1]+range_in_km*boresight_in_tbf_[1];
  borepos_[2]=scpos_[2]+range_in_km*boresight_in_tbf_[2];

  nominal_range=Uvar(range_in_km,"km");
  if(nadir_mode)nominal_range=cassini_sim.altitude();
  nominal_delay=nominal_range*2/speed_light;
  
  // assign nominal_doppler_centroid to doppler at nominal_range and
  // the appropriate azimuth offset
  // assumes boresight is on the surface
  double look[3];
    
  Frame tbf(default_target_frame_spice_id,default_target_spice_id);
  PositionVector pos(tbf,ground_impact_time_in_s,scpos_[0],scpos_[1],
		     scpos_[2]);
  FloatVector vel(tbf,ground_impact_time_in_s,scvel_[0],scvel_[1],scvel_[2]);
  scstate_=StateVector(pos,vel);
  getFastLookVector(range_in_km,beam_azioff_rad[beam_number-1],look);
  nominal_doppler_centroid=Uvar(getFastDopplerCentroid(look),"Hz");
#endif
  // END of NEW NEW test test test
 if( dbg.level){
   dbg.file.precision(20);
   
   dbg.file << "borepos_=[" << borepos_[0] << " " << borepos_[1] << " "
	    << borepos_[2] <<"] km" << endl;
   dbg.file << "scpos_=[" << scpos_[0] << " " << scpos_[1] << " "
	    << scpos_[2] <<"] km" << endl;
   dbg.file << "scvel_=[" << scvel_[0] << " " << scvel_[1] << " "
	    << scvel_[2] <<"] km/s" << endl;
   dbg.file << "scacc_=[" << scacc_[0] << " " << scacc_[1] << " "
	    << scacc_[2] <<"] km/s" << endl;
   dbg.file << "boresight_in_tbf_=[" << boresight_in_tbf_[0] << " " 
	    << boresight_in_tbf_[1] << " "
	    << boresight_in_tbf_[2] <<"]" << endl;
   dbg.file << "nominal_doppler_centroid=" << nominal_doppler_centroid
	    << " Hz nominal_delay=" << nominal_delay << " s" << endl;
   dbg.file << "nominal range=" << nominal_range << " km" << endl;
   dbg.file << "ground_impact_time=" << ground_impact_time_in_s << " s"
	    << endl; 
   dbg.file << "rotmat_=" << endl;
   for(int i=0;i<3;i++){
     for(int j=0;j<3;j++){
       dbg.file << rotmat_[i][j] << " ";
     }
     dbg.file << endl;
   }
   dbg.file << "rotmat_rev_=" << endl;
   for(int i=0;i<3;i++){
     for(int j=0;j<3;j++){
       dbg.file << rotmat_rev_[i][j] << " ";
     }
     dbg.file << endl;
   }
   dbg.file << "End Debugging L1I::initGeometryCalculations Burst #"
	      << sab_counter << endl;
  }
  //---------------------------------------------------
  // perform Area calculation initializations if necessary
  //---------------------------------------------------
  }



//------------------------
// routines for getting the true unadulterated range and Doppler
// associated with a pixel (Range and Doppler centroid offsets removed)
// used to produce range and doppler backplanes for stereo processing
//------------------------
double L1I::rawRange(double fdop_in_Hz, double range_in_km){
  double rawdop=rawDoppler(fdop_in_Hz);
  double range_correction=rawdop-nominal_doppler_centroid_in_Hz;
  range_correction*=range_doppler_correction_coeff;
  // Old way  range centroid effect geo-location
  // double rawrange=range_in_km+range_correction+rc_range_shift_in_km;
  // New way  range centroid effects gain only
  double rawrange=range_in_km+range_correction;
  return(rawrange);
}
double L1I::rawDoppler(double fdop_in_Hz){
  double rawdop;
  // OLD way doppler centroid estimate effects geolocation
  if(quick_look_dopran)
     rawdop=fdop_in_Hz+ac_dop_shift_in_Hz;
  // New way doppler centroid estiamte only effects Gain
  else
    rawdop=fdop_in_Hz;
  return(rawdop);
}

// routine for getting estimated area in square km of a resolution pixel
// used by BIDR in TOPO mode
// for now it does not vary with location within burst
// and assumes rectangular pixels
double L1I::resolutionArea(SARProcParams& spp){
  double chirp_bw=csq*slow_cfs_in_Hz;
  float range_res=get_in_base_units(speed_light)/(2*chirp_bw);
  float doppler_res=1/(num_pulses_received*pri_in_s);
  float rcoeff=get_in_base_units(spp.range_res_coeff*range_res/range_step_);
  float acoeff=get_in_base_units(spp.azimuth_res_coeff*doppler_res/doppler_step_);
  return(rcoeff*acoeff*range_ground_width_*doppler_ground_width_);
}
// Output to Range And Doppler File (Doppler Centroid Tracker format)

void L1I::outputRangeDoppler(RangeDopplerFile& rdf){

rdf.loadBurstData(t,act_geom_time_offset, sab_counter,beam_number-1, 1.0/pri);
Uvec r("",num_range_bins);
Dvec fd("",num_range_bins);
Uvar nom_range=nominal_delay*speed_light/2.0;

// Adjustment factor included to account for track_doppler erroneously
// ignoring chirp start frequency
//double adjust_dop=get_in_base_units(lambda_chirp/cassini_sim.lambda());
 double adjust_dop=1.0; // parameter that can be used to account for known
                        // doppler errors --- see above 

// doppler adjustment applied here
double nom_dop=get_in_base_units(nominal_doppler_centroid*adjust_dop);


int intdop=(int)floor(nom_dop*pri_in_s+0.5);
double fdave=0.0;

for(unsigned int c=0;c<num_range_bins;c++){
	r(c)=Uvar(range[c],"km");
        // doppler adjustment applied here
	fd(c)=fdop_centroid[c]*adjust_dop-intdop/pri_in_s;
        fd(c)*=pri_in_s;
        fdave+=fd(c);
}
fdave/=num_range_bins;

// output fdnom instead of fdave to compare with Doppler Tracker geometry
double fdnom= nom_dop*pri_in_s-intdop;

rdf.loadComputedRangeDoppler(num_range_bins,r,fd,nom_range,fdnom,intdop);
rdf.saveLoadedData();
}
