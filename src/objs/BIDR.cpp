static const char rcs_id_bidr_c[] =
  "@(#) $Id: BIDR.cpp,v 11.11 2013/10/15 23:02:19 bstiles Exp $";

#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include "BIDR.h"
#include "config_keywords.h"
#include "Constants.h"
#include "SimpleArray.h"
#include "Utils.h"

using std::cin;
using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::ios;
using std::setfill;
using std::setprecision;
using std::setw;
using std::sort;
using std::stringstream;

//-------------------------
// Static initializations
//-------------------------

string BIDR::pds_version_id_ = "PDS3";
string BIDR::record_type_ = "FIXED_LENGTH";

string BIDR::data_set_id_ = "CO-SSA-RADAR-5-BIDR-V1.0";
string BIDR::data_set_name_ = "CASSINI ORBITER SSA RADAR 5 BIDR V1.0";
string BIDR::producer_institution_name_ = "JET PROPULSION LABORATORY";
string BIDR::producer_id_ = "JPL";
string BIDR::instrument_host_name_ = "CASSINI ORBITER";
string BIDR::instrument_host_id_ = "CO";
string BIDR::instrument_name_ = "CASSINI RADAR";
string BIDR::instrument_id_ = "RADAR";
string BIDR::mission_name_ = "CASSINI-HUYGENS";
string BIDR::sample_type_[BIDRTYPECOUNT] =
{
  "PC_REAL", // INC
  "UNSIGNED INTEGER", // BEAMMASK
  "PC_REAL", // S0_UNCORR
  "PC_REAL", // S0_CORR
  "PC_REAL", // LAT
  "PC_REAL", // LON
  "INTEGER", // START_BURST_NUM     -- nondeliverable product
  "INTEGER", // END_BURST_NUM       -- nondeliverable product
  "UNSIGNED INTEGER", // NUM_LOOKS
  "UNSIGNED INTEGER", // S0_CORR_DB -- nondeliverable product (from JPL)
  "PC_REAL",  // DOPPLER
  "PC_REAL",  // RANGE
  "INTEGER",  // NOMINAL_BURST_NUM
  "PC_REAL",  // S0NSQT_UNCORR
  "PC_REAL",  // SO_STD
  "PC_REAL",  // SO_NOISE_EQUIV
  ""  // UNUSED
};
int BIDR::sample_bits_[BIDRTYPECOUNT] =
{
  BITS_PER_BYTE * sizeof(float), // INC
  BITS_PER_BYTE * sizeof(unsigned char),  // BEAMMASK
  BITS_PER_BYTE * sizeof(float), // S0_UNCORR
  BITS_PER_BYTE * sizeof(float), // S0_CORR
  BITS_PER_BYTE * sizeof(float), // LAT
  BITS_PER_BYTE * sizeof(float), // LON
  BITS_PER_BYTE * sizeof(int), // START_BURST_NUM
  BITS_PER_BYTE * sizeof(int), // END_BURST_NUM
  BITS_PER_BYTE * sizeof(unsigned char), // NUM_LOOKS
  BITS_PER_BYTE * sizeof(unsigned char),  // S0_CORR_DB
  BITS_PER_BYTE * sizeof(float), // DOPPLER
  BITS_PER_BYTE * sizeof(float), // RANGE
  BITS_PER_BYTE * sizeof(int), // NOMINAL_BURST_NUM
  BITS_PER_BYTE * sizeof(float), // S0NSQT_UNCORR
  BITS_PER_BYTE * sizeof(float), // S0_STD
  BITS_PER_BYTE * sizeof(float), // S0_NOISE_EQUIV
  0  // UNUSED
};

double BIDR::scaling_factor_ = 1.0;
double BIDR::offset_ = 0.0;
string BIDR::missing_constant_[BIDRTYPECOUNT] =
{
  ISIS_NULL, // INC
  "0",       // BEAMMASK
  ISIS_NULL, // S0_UNCORR
  ISIS_NULL, // S0_CORR
  ISIS_NULL, // LAT
  ISIS_NULL, // LON
  "0", // START_BURST_NUM
  "0", // END_BURST_NUM
  "0", // NUM_LOOKS
  "0", // S0_CORR_DB
  ISIS_NULL, // DOPPLER
  ISIS_NULL, // RANGE
  "0", // NOMINAL_BURST_NUM
  ISIS_NULL, // S0NSQT_UNCORR
  ISIS_NULL, // S0_STD
  ISIS_NULL, // S0_NOISE_EQUIV
  ""  // UNUSED
};
string BIDR::note_[BIDRTYPECOUNT] =
{
"The data values in this file are the incidence angles between the \
look vector and the surface normal at each pixel.",
"The data values in this file identify the beam(s) used to produce \
the backscatter values for each pixel.",
"The data values in this file are Synthetic Aperture Radar (SAR) \
normalized backscatter cross-section values.  The values are \
physical scale (not in dB) and have not been corrected for \
incidence-angle effects. Biases due to noise have not been \
removed.",
"The data values in this file are Synthetic Aperture Radar (SAR) \
normalized backscatter cross-section values.  The values are \
physical scale (not in dB) and have been corrected for \
incidence-angle effects and biases due to thermal and quantization \
noise have been removed.  The raw backscatter values have been \
multiplied by the function f(I), where I is the incidence angle \
and f(I) = 0.2907/(f1(I)+f2(I)+f3(I)), for \
f1(I)=2.8126*(cos(I)^4+893.9677*sin(I)^2)^(-1.5), \
f2(I)=0.5824*(cos(I)^4+34.1366*sin(I)^2)^(-1.5), \
and f3(I)=0.3767*cos(I)^1.9782.",
"The data values in this file are the ordinary latitudes \
in the body fixed coordinate system for Titan.",
"The data values in this file are the ordinary longitudes (positive \
west) in the body fixed coordinate coordinate system for Titan.",
"The data values in this file are the indices of the starting (first) \
burst used to produce the backscatter values for each pixel.",
"The data values in this file are the indices of the ending (last) \
burst used to produce the backscatter values for each pixel.",
"The data values in this file are the number of looks used to produce \
the backscatter values for each pixel.",
"",
"The data values in this file is the Doppler freq in Hz corresponding \
to the pixel location in the nominal burst.",
"The data values in this file is the Range in km corresponding \
to the pixel location in the nominal burst.",
"The data values in this file are the indices of the nominal (center) \
burst used to produce the backscatter values for each pixel.",
"The data values in this file are Synthetic Aperture Radar (SAR) \
normalized backscatter cross-section values.  The values are \
physical scale (not in dB) and have not been corrected for \
incidence-angle effects. Biases due to thermal and quantization \
noise have been removed.",
"The data values in this file are standard deviations of \
Synthetic Aperture Radar (SAR) normalized backscatter cross-section \
noise subtracted values w/o incidence angle correction.  \
The values are physical scale (not in dB).",
"The data values in this file are noise equivalent sigma-0 values \
associated with Synthetic Aperture Radar (SAR) normalized backscatter \
cross-section noise subtracted values w/o incidence angle correction.  \
The values are physical scale (not in dB).",
""
};

string BIDR::map_projection_type_ = "OBLIQUE CYLINDRICAL";
string BIDR::first_standard_parallel_ = "N/A";
string BIDR::second_standard_parallel_ = "N/A";
string BIDR::positive_longitude_direction_ = "WEST";
double BIDR::proj_center_latitude_ = 0.0;
double BIDR::proj_center_longitude_ = 0.0;
int BIDR::line_first_pixel_ = 1;
int BIDR::sample_first_pixel_ = 1;
double BIDR::map_projection_rotation_ = 90.0;
string BIDR::coordinate_system_type_ = "BODY-FIXED ROTATING";
string BIDR::coordinate_system_name_ = "PLANETOGRAPHIC";

//--------------
// Constructors
//--------------


BIDR::BIDR(const string& prefix, Config& cfg, L1B& l1b,
  SARProcParams& spp, const string& mode, BIDRModeE proc_mode=NORMAL)
  :  procMode(proc_mode),
     inc_file(NULL),
     beammask_file(NULL),
     s0_uncorr_file(NULL),
     s0_corr_file(NULL),
     lat_file(NULL),
     lon_file(NULL),
     start_burst_num_file(NULL),
     end_burst_num_file(NULL),
     num_looks_file(NULL),
     X_file(NULL),
     std_file(NULL),
     s0ne_file(NULL),
     beam_deltas_file(NULL),
     dop_file(NULL),
     range_file(NULL),
     nom_burst_num_file(NULL),
     topoplus_file(NULL),
     num_lines_of_output(0),
     num_looks_required(0),
     max_looks_allowed(100000000),
     num_looks_to_skip(0),
     enable_max_looks(0),
     maximum_latitude_(-90),
     minimum_latitude_(90),
     maximum_longitude_(-180),
     minimum_longitude_(180),
     first_longitude_(-1000),
     auto_overwrite_bidr_(false),
     check_pds_string_lengths_(true),
     use_config_data_take_number_(false)
     
  {
  overlap_set=false;
  initProcMode();
  if((procMode==TOPO || procMode==TOPO_PLUS) && !spp.noise_subtraction_on){
    ErrorMessage e("BIDR Fatal Error: Topo mode requires noise_subtraction_on");
    e.throwMe();
  }
  // set up bad value float by specific bits.
  int val=BAD_VALUE_HEX;
  float* ptr=(float*) &val;
  BAD_VALUE=*ptr;

  DebugInfo dbg("BIDR::BIDR");

  // For now if mode is not write throw an error
  if(mode!="w"){
    ErrorMessage e("For now the BIDR constructor only works for write mode");
    e.throwMe();
  }

  if (cfg.keywordExists("SAR_PROC_AUTO_OVERWRITE_BIDR"))
  {
    auto_overwrite_bidr_ = (bool) cfg.getInt("SAR_PROC_AUTO_OVERWRITE_BIDR");
  }
  openFiles(prefix,mode,spp.noise_subtraction_on);

  if (cfg.keywordExists("CHECK_PDS_STRING_LENGTHS"))
  {
    check_pds_string_lengths_ = (bool) cfg.getInt("CHECK_PDS_STRING_LENGTHS");
  }
  if (!check_pds_string_lengths_)
  {
    // Display a warning message.  Normally, this option should be enabled.
    cerr << "Warning:  PDS string length limit checks have been disabled." <<
      endl;
    cerr << "To enable the checks, set CHECK_PDS_STRING_LENGTHS to 1 in the " <<
      "config file." << endl;
  }

  // Determine whether of not oblique cylindrical will be used
  use_oblique_cyl=(bool)cfg.getInt("USE_OBLIQUE_CYLINDRICAL");

  // beam feathering
  feather_beams=(bool)cfg.getInt("ENABLE_BEAM_FEATHERING");
  if(feather_beams && spp.single_beam_mode){
    cerr<<"Warning Beam Feathering enabled in Single Beam Mode Ignoring ..."
	<< endl;
    feather_beams=false;
  }
  
  dominant_beam3=false;
  if(cfg.keywordExists("SAR_PROC_DOMINANT_BEAM_3")){
    dominant_beam3=(bool)cfg.getInt("SAR_PROC_DOMINANT_BEAM_3");
  }

  mark_nadir_amb=false;
  if(cfg.keywordExists("BIDR_MARK_NADIR")){
    mark_nadir_amb=(bool)cfg.getInt("BIDR_MARK_NADIR");
  }


  output_beam_deltas=false;
  beam_deltas_file=NULL;
  if(cfg.keywordExists("SAR_PROC_BEAM_DELTAS_FILE")){
    output_beam_deltas=true;
    string tmpstr=cfg.str("SAR_PROC_BEAM_DELTAS_FILE");
    beam_deltas_file=fopen(tmpstr,"w");
  }

  // Determine Incidence Angle Correction Model
  string incmod_string=cfg.str("INCIDENCE_ANGLE_CORRECTION_MODEL");
  if(incmod_string=="VENUS" || incmod_string=="Venus"){
    inc_correction_model=VENUS;
    venus_k1=0.0188;
    venus_k2=0.111;
    venus_s045=muhleman_backscatter(venus_k1,venus_k2,Uvar(45,"deg"));
  }

  // Correct name for this backscatter model is DBLINE041104 as it was 
  // produced on Nov 4 2004. Original incorrect date is allowed for
  // backward compatibility.
  else if(incmod_string=="dBLine041004" || incmod_string=="DBLINE041004"
	  || incmod_string=="dBLine041104" || incmod_string=="DBLINE041104"){
    inc_correction_model=DBLINE041104;
    dbline_offset=4.1667;
    dbline_slope=-0.2333*180/pi;
    dbline_s045=pow(10,-0.6333);
  }
  else if(incmod_string=="RKIRK01TA"){
    inc_correction_model=RKIRK01TA;
    invsin_coeff=0.0814;
    invsin_s045=0.0814/sin(45*pi/180);
  }
  else if(incmod_string=="HAGHAG"){
    inc_correction_model=HAGHAG;
  }
  else if(incmod_string=="ENCELADUS_WYE1"){
    inc_correction_model=ENCELADUS_WYE1;
  }
  else if(incmod_string=="RHEA_WYE1"){
    inc_correction_model=RHEA_WYE1;
  }
  else if(incmod_string=="DIONE_WYE1"){
    inc_correction_model=DIONE_WYE1;
  }
  else if(incmod_string=="INVALID" || incmod_string=="Invalid"){
    inc_correction_model=INVALID_CORR;
  }
  else{
    ErrorMessage e("BIDR::config BAD INCIDENCE_ANGLE_CORRECTION_MODEL "+ incmod_string);
    e.throwMe();
  }
  if(inc_correction_model!=HAGHAG){
    fprintf(stderr,"Warning Using Obsolete INC Correction MODEL\n");
    fprintf(stdout,"Warning Using Obsolete INC Correction MODEL\n");
  }
  // if so set up projection accordingly
  if(use_oblique_cyl){
    StateVector s=spp.getClosestApproachState();
    proj=OblCylProj(s);
  }
  look_direction_=spp.LookDirection();

  // if not oblique cylindrical use default (trivial) Projection (do nothing here)
  
  // determine number of looks required to image
  num_looks_required=cfg.getInt("SARPROC_NUM_LOOKS_REQUIRED");
  enable_max_looks=false;
  if(cfg.keywordExists("SARPROC_MAX_LOOKS_ALLOWED")){
    max_looks_allowed=cfg.getInt("SARPROC_MAX_LOOKS_ALLOWED");
    num_looks_to_skip=cfg.getInt("SARPROC_NUM_LOOKS_SKIPPED");
    cerr << "Warning: only the first " << max_looks_allowed
	 << " looks are used for each pixel in the BIDR!!!" << endl;
    enable_max_looks=true;
  }
  // Determine image resolution
  pixelsperdegree=cfg.getInt("BIDR_PIXELS_PER_DEGREE");
  

  // The IGNORE_GAPS mode is needed if you want to process a swath
  // from a LBDR file with large gaps. It omits some error checking which
  // fails if there are gaps.
  if(cfg.keywordExists("BIDR_IGNORE_GAPS")){
    ignore_gaps_mode=(bool)cfg.getInt("BIDR_IGNORE_GAPS");
  }
  else ignore_gaps_mode=false;

  if(cfg.keywordExists("BIDR_ALLOW_DECREASING_LON")){
    allow_decreasing_lon_mode=(bool)cfg.getInt("BIDR_ALLOW_DECREASING_LON");
  }
  else allow_decreasing_lon_mode=false;

  // The REPLACE_X Mode is used to output arbitrary parameters to the X
  // backplane. As a side effect the new X paarmeter is the tie breaker for
  // choosing which beam to process.
  if(cfg.keywordExists("BIDR_REPLACE_X")){
    replace_X_mode=(bool)cfg.getInt("BIDR_REPLACE_X");
    if(replace_X_mode){
      replace_X_param=cfg.str("BIDR_REPLACE_X_PARAM");
    }
  }
  else replace_X_mode=false;

  if(cfg.keywordExists("SAR_PROC_USE_CONFIG_DATA_TAKE_NUMBER")){
    use_config_data_take_number_  =
      cfg.getInt("SAR_PROC_USE_CONFIG_DATA_TAKE_NUMBER");
  }
  if (use_config_data_take_number_)
  {
    data_take_id_ = cfg.getInt("data_take_number");
    if (dbg.level) {
      dbg.file << "Data take ID from config file = " << data_take_id_ << endl;
    }
  }

  // Read or set other configuration parameters that are used to set values
  // in the PDS labels

  producer_full_name_ = PDSLabel::replaceSpaces(cfg.str("PRODUCER_FULL_NAME"));
  target_name_ = cfg.str("target");
  source_product_id_ = l1b.productID();
  if (source_product_id_.length() == 0)
  {
    // Should never occur, but set value to a non-empty string to avoid
    // writing an illegal PDS statement
    cerr << "Warning:  cannot parse PRODUCT_ID from " <<
      L1B_ACTIVE_MODE_FILENAME << " in BIDR::BIDR()" << endl;
    source_product_id_ = "UNKNOWN";
  }
  product_version_ = cfg.getInt("BIDR_PRODUCT_VERSION_ID");
  mission_phase_name_ = cfg.str("MISSION_PHASE_NAME");
  data_set_map_projection_catalog_ = cfg.str("DATA_SET_MAP_PROJECTION_CATALOG");
  software_version_ = cfg.str("BIDR_SOFTWARE_VERSION_ID");

  // A valid FLYBY_ID is the character "T" followed by either "a", "A", or
  // a one- to three-digit number.  The number may not consist of all zeroes,
  // but leading zeroes are allowed (although not strictly legal).
  string flyby_id = cfg.str("FLYBY_ID");
  flyby_id_pass_ = "";
  int flyby_id_len = flyby_id.length();
  if (flyby_id_len < 2 || flyby_id_len > 4)
  {
    ErrorMessage e("BIDR::config:  illegal value " + flyby_id +
      " for FLYBY_ID");
    e.throwMe();
  }
  else if (flyby_id_len == 2)
  {
    flyby_id_pass_ = flyby_id.substr(1);
    if (flyby_id[1] != 'a' && flyby_id[1] != 'A' &&
        !(flyby_id[1] >= '1' && flyby_id[1] <= '9'))
    {
      ErrorMessage e("BIDR::config:  illegal value " + flyby_id +
	" for FLYBY_ID");
      e.throwMe();
    }
  }
  else
  {
    flyby_id_pass_ = flyby_id.substr(1);
    // Strip leading zeroes
    while (flyby_id_pass_[0] == '0' && flyby_id_pass_.length() > 1)
    {
      flyby_id_pass_ = flyby_id_pass_.substr(1);
    }
    if (flyby_id_pass_[0] == '0')
    {
      ErrorMessage e("BIDR::config:  illegal value " + flyby_id +
	" for FLYBY_ID");
      e.throwMe();
    }
    for (unsigned int i = 0; i < flyby_id_pass_.length(); i++)
    {
      if (flyby_id_pass_[i] < '0' || flyby_id_pass_[i] > '9')
      {
	ErrorMessage e("BIDR::config:  illegal value " + flyby_id +
	  " for FLYBY_ID");
	e.throwMe();
      }
    }
  }
  flyby_id_pass_ = PDSLabel::toUpper(flyby_id_pass_);

  if(flyby_id_pass_ == "A" && procMode!=TOPO && procMode!=TOPO_PLUS){

      char resp[5];
      cerr << "Should the Flyby ID be TA?  ([nN]/yY)  ";
      cout << "Should the Flyby ID be TA?  ([nN]/yY)  ";
      cin.getline(resp, 5);
      if (resp[0] != 'y' && resp[0] != 'Y')
      {
        ErrorMessage e("Fix the FLYBY_ID in the config file.  Exiting.");
	e.throwMe();
      }
  }
  // Determine segment number
  segment_id_=cfg.getInt("BIDR_SEGMENT_ID");
  if(segment_id_<0 || segment_id_>49){
    fprintf(stderr,"Fatal Error: Bad BIDR_SEGMENT_ID. Should be between 0 and 49\n");
    fprintf(stdout,"Fatal Error: Bad BIDR_SEGMENT_ID. Should be between 0 and 49\n");
    exit(1);
  }
  // Determine Latitude Bounds from L1B file
  // May need to move this functionality to a more suitable place later
  // Preprocessor could compute this and put in header.
  // or SARProcParams could initialize a number of things like this at once.
  // Right now I am placing it here so that the Projection object is available

  double maxlat=-pi/2;
  double minlat=pi/2;
  double firstlon=0;
  double minlon=4*pi;
  double maxlon=-4*pi;
  int valid_burst_no=0;
  int minlon_idx=0;
  int data_take_id_lbdr=0;
  while(!l1b.eof()){

    // Reading individual parameters instead of whole record 
    // for faster (but more fragile) access 
    // 
    l1b.readParameter("sclk"); // needed for PDS label start/stop timestamps
    l1b.readParameter("brst"); // needed for PDS label start/stop timestamps
    l1b.readParameter("t"); // needed for PDS label start/stop timestamps
    l1b.readParameter("beam_number"); // needed for single beam mode check
    l1b.readParameter("adc"); // needed for isSAR check
    l1b.readParameter("csr"); // needed for isCAl check in badSARData
    l1b.readParameter("engineer_qual_flag"); // needed for badSARData check
    l1b.readParameter("science_qual_flag"); // needed for badSARData check
    l1b.readParameter("act_centroid_lon"); // needed to compute lat bounds
    l1b.readParameter("act_centroid_lat"); // needed to compute lat bounds
    l1b.readParameter("record_id"); // needed for spp::inRegion
    l1b.readParameter("time_from_closest_approach"); // needed for spp.inRegion
    l1b.readParameter("time_from_epoch"); // needed for spp.inRegion
    l1b.readParameter("sab_counter"); // needed for quality override check
    l1b.readParameter("dtn"); // needed for PDS label product ID
    l1b.readParameter("num_pulses_received"); // needed for check in skipBurst
    l1b.skipRecord(); // go to beginning of next next record


    // only consider records which will be used by the SAR processor

    // perform checks to see if burst should be processed or skipped
    if(spp.skipBurst(proj,l1b)) continue;


    valid_burst_no++; // update count of bursts to be processed

    // set the start time only for the first valid burst; set the stop time
    // to be that of the current burst
    if (!sclk_start_.valid())
    {
      sclk_start_ = Time(l1b.t);
    }
    sclk_stop_ = Time(l1b.t);

    // read the data take ID, which should be the same for all bursts
    data_take_id_lbdr = l1b.dtn;

    // update lat,lon bounds if necessary
    double slat,slon;
    double lat,lon;
    slat=l1b.act_centroid_lat.getInUnits("rad");
    slon=l1b.act_centroid_lon.getInUnits("rad");
    lat=proj.latInRad(slon,slat);
    lon=proj.lonInRad(slon,slat);
    if(lat < minlat) minlat=lat;
    if(lat > maxlat) maxlat=lat;

    if(valid_burst_no==1){
      firstlon=lon;
    }
  
    // make certain longitude is in range (first_lon-180, first_lon+180)
    while(lon>firstlon+pi) lon-=2*pi;
    while(lon<firstlon-pi) lon+=2*pi;

    // This now works w/o worrying about wraparound so long as the
    // whole pass has less than 180 degrees of SAR coverage
    // which should always be the case for Cassini

    if(lon < minlon){
      minlon=lon;
      minlon_idx=valid_burst_no;
    }
    if(lon > maxlon) maxlon=lon;
  }

  // Check for no valid bursts
  if(valid_burst_no==0){
    cerr << "Fatal Error:BIDR constructor found no processable bursts." << endl;
    cerr << " Check to make sure all necesary parameters are being read for the various checks." << endl;
    cerr << " Check to make sure the region, radar mode," << endl;
    cerr << "and single beam SAR processor options are configured correctly." << endl;
    exit(1);
  }
  l1b.gotoFirstRecord(); // goes back to beginning of the file.

  // convert minlat and maxlat to integer degrees add desired pad

  // NEED TO MAKE THIS READ FROM THE CONFIG FILE
  Uvar latitude_pad=cfg["BIDR_LATLON_PAD_EACH_SIDE"];
  double lpf=latitude_pad.getInUnits("deg");
  minlat=floor(minlat*radtodeg-lpf);
  maxlat=ceil(maxlat*radtodeg+lpf);
  if(minlat<-90) minlat=-90; 
  if(maxlat>90) maxlat=90; 

  // compute the starting longitude for the grid and the direction of
  // motion 
  bool lon_inc; // true if longitude is increasing

  // This equation fails if a short time region of multiple beam bursts
  // is being processed but for that case direction of motion is irrelevant
  // because the whole region being processed will fit in one longitude
  // window.
  lon_inc=(minlon_idx<valid_burst_no/2);
  if(!lon_inc && !allow_decreasing_lon_mode){
    cerr << "Warning: BIDR::BIDR computed that oblcyl longitude is NOT"
	 << endl << "Warning (cont): increasing with time. " << endl;
    cerr << "Warning: lon_inc is being explicitly set to true. " << endl;
    cerr << "Warning( cont): OK for short time runs ...." << endl;
    lon_inc=true;
  }
  double start_lon;


  if(lon_inc) start_lon=floor(minlon*radtodeg-lpf);
  else start_lon=ceil(maxlon*radtodeg+lpf);

  double lon_width=ceil((maxlon-minlon)*radtodeg+2*lpf);

  bool lon360=false;
  int lon_buffer_in_deg=9;
  Uvar lon_buffer_sz;
  if(cfg.keywordExists("BIDR_360_LON_BUFFER")){
    lon360=cfg.getInt("BIDR_360_LON_BUFFER");
  }

  if(cfg.keywordExists("BIDR_LON_BUFFER_SIZE")){
    lon_buffer_sz=cfg["BIDR_LON_BUFFER_SIZE"];
   // this must be a integer factor of 360 (degrees) or longitude wraps around
   // incorrectly
    lon_buffer_in_deg=(int)(lon_buffer_sz.getInUnits("deg")+0.5);
    int chk=360%lon_buffer_in_deg;
    if(chk!=0){
      fprintf(stderr,"BAD BIDR_LON_BUFFER_SIZE should be factor of 360\n");
      exit(1);
    } 
  }


  // set up LatLonGrid
  grid=LatLonGrid(pixelsperdegree,minlat,maxlat,start_lon,lon_width,lon_inc,lon360,lon_buffer_in_deg);

  // Allocate arrays based on resolution and latitude bounds
  allocateArrays(spp);

  // Initialize standard lat/lon buffers
  int num_lats=grid.numLats();
  int num_lons=grid.numLons();
  int abs_start_i=grid.absStartValidIndex();
  int start_i=grid.relLonIndex(abs_start_i);

  // The following two double for loops fail if lon_inc is false
  // So we have made sure it never can be. The OblCyl definition
  // precludes it anyway.

  if(lon_inc){
    for(int i=start_i;i<num_lons;i++){
      int abs_i=abs_start_i+i-start_i;
      for(int j=0;j<num_lats;j++){
	double lon=grid.lonInRad(abs_i);
	double lat=grid.latInRad(j);
	slat_bfr[i][j]=proj.standardLatInRad(lon,lat);
	slon_bfr[i][j]=proj.standardLonInRad(lon,lat);
      }
    }
  
    for(int i=0;i<start_i;i++){
      int abs_i=abs_start_i+num_lons-start_i+i;
      for(int j=0;j<num_lats;j++){
	double lon=grid.lonInRad(abs_i);
	double lat=grid.latInRad(j);
	slat_bfr[i][j]=proj.standardLatInRad(lon,lat);
	slon_bfr[i][j]=proj.standardLonInRad(lon,lat);
      }
    }
  }
  else{
    int abs_end_i=abs_start_i -num_lons +1;
    for(int i=start_i;i<num_lons;i++){
      int abs_i=abs_end_i+i-start_i;
      for(int j=0;j<num_lats;j++){
	double lon=grid.lonInRad(abs_i);
	double lat=grid.latInRad(j);
	slat_bfr[i][j]=proj.standardLatInRad(lon,lat);
	slon_bfr[i][j]=proj.standardLonInRad(lon,lat);
      }
    }
  
    for(int i=0;i<start_i;i++){
      int abs_i=abs_end_i+num_lons-start_i+i;
      for(int j=0;j<num_lats;j++){
	double lon=grid.lonInRad(abs_i);
	double lat=grid.latInRad(j);
	slat_bfr[i][j]=proj.standardLatInRad(lon,lat);
	slon_bfr[i][j]=proj.standardLonInRad(lon,lat);
      }
    }
  }
  if (dbg.level) {
    double latinrad = grid.latInRad(0);
    double loninrad = grid.lineNumberToLonInRad(0);
    dbg.file << "OC lat, lon at (0, 0) (rad): " << latinrad << " " << loninrad << endl;
    double latindeg = radtodeg * proj.standardLatInRad(loninrad, latinrad);
    double lonindeg = radtodeg * proj.standardLonInRad(loninrad, latinrad);
    dbg.file << "St lat, lon at (0, 0) (deg): " << latindeg << " " << lonindeg << endl;
  }

  if (!use_config_data_take_number_)
  {
    // Set data take ID to the value read from the last processable burst
    data_take_id_ = data_take_id_lbdr;
    if (dbg.level) {
      dbg.file << "Data take ID from LBDR file = " << data_take_id_ << endl;
    }
  }
}

void BIDR::openFiles(const string& prefix, const string& mode, int enable_noise_sub){

  if(mode!="w"){
    ErrorMessage e("BIDR::openFiles mode must be w.");
    e.throwMe();
  }

  string inc_file_name = prefix + ".inc";
  string beammask_file_name = prefix + ".beammask";
  string s0_uncorr_file_name = prefix + ".s0_uncorr";
  string s0_corr_file_name;
  if(enable_noise_sub) s0_corr_file_name = prefix + ".s0nsqt_corr";
  else s0_corr_file_name = prefix + ".s0_corr";
  string X_file_name;
  if(!enable_noise_sub) X_file_name=prefix + ".X";
  else X_file_name=prefix+".s0nsqt_uncorr";
  string std_file_name = prefix + ".s0nsqt_uncorr_std";
  string s0ne_file_name = prefix + ".s0nsqt_noise_equiv";
  string lat_file_name = prefix + ".lat";
  string lon_file_name = prefix + ".lon";
  string start_burst_num_file_name = prefix + ".start_burst_num";
  string end_burst_num_file_name = prefix + ".end_burst_num";
  string num_looks_file_name = prefix + ".num_looks";
  string dop_file_name = prefix + ".doppler";
  string range_file_name = prefix + ".range";
  string nom_burst_num_file_name = prefix + ".nominal_burst_num";
  string topoplus_file_name = prefix + ".topoplus";
  if (!auto_overwrite_bidr_)
  {
    // If any of the BIDR output files exist, ask the user whether
    // to overwrite all the files.  Any response other than "y"
    // or "Y" is interpreted as "no."
    if (fileExists(inc_file_name) ||
        fileExists(beammask_file_name) ||
        fileExists(s0_uncorr_file_name) ||
        fileExists(s0_corr_file_name) ||
        fileExists(lat_file_name) ||
        fileExists(lon_file_name) ||
        fileExists(dop_file_name) ||
        fileExists(range_file_name) ||
        fileExists(nom_burst_num_file_name) ||
        fileExists(start_burst_num_file_name) ||
        fileExists(end_burst_num_file_name) ||
        fileExists(X_file_name) ||
        fileExists(std_file_name) ||
	fileExists(s0ne_file_name) ||
        fileExists(num_looks_file_name))
    {
      char resp[5];
      cerr << "Overwrite existing BIDR output files?  ([nN]/yY)  ";
      cout << "Overwrite existing BIDR output files?  ([nN]/yY)  ";
      cin.getline(resp, 5);
      if (resp[0] != 'y' && resp[0] != 'Y')
      {
        ErrorMessage e("BIDR output files exist and will not be overwritten.  Exiting.");
	e.throwMe();
      }
    }
  }

  X_file=fopen(X_file_name, mode);
  if(output_enable_[INC])inc_file=fopen(inc_file_name, mode);
  if(output_enable_[BEAMMASK])beammask_file=fopen(beammask_file_name, mode);
  if(output_enable_[S0_UNCORR])s0_uncorr_file=fopen(s0_uncorr_file_name, mode);
  if(output_enable_[S0_CORR])s0_corr_file=fopen(s0_corr_file_name, mode);
  if(output_enable_[S0_STD])std_file=fopen(std_file_name, mode);
  if(output_enable_[S0_NOISE_EQUIV])s0ne_file=fopen(s0ne_file_name, mode);
  if(output_enable_[LAT])lat_file=fopen(lat_file_name, mode);
  if(output_enable_[LON])lon_file=fopen(lon_file_name, mode);
  if(output_enable_[START_BURST_NUM])
    start_burst_num_file=fopen(start_burst_num_file_name, mode);
  if(output_enable_[END_BURST_NUM])
    end_burst_num_file=fopen(end_burst_num_file_name, mode);
  if(output_enable_[NUM_LOOKS])
    num_looks_file=fopen(num_looks_file_name, mode);
  if(output_enable_[NOMINAL_BURST_NUM])
    nom_burst_num_file=fopen(nom_burst_num_file_name, mode);
  if(output_enable_[DOPPLER])
    dop_file=fopen(dop_file_name, mode);
  if(output_enable_[RANGE])
    range_file=fopen(range_file_name, mode);
  // fopen overload in Utils.cpp takes string inputs and throws errors on failure
  if(procMode==TOPO_PLUS) topoplus_file=fopen(topoplus_file_name, mode);
}



void
BIDR::allocateArrays(const SARProcParams& spp)
{

  DebugInfo dbg("BIDR::allocateArrays");
  int num_lats=grid.numLats();
  int num_lons=grid.numLons();
  num_beams=5; // in single beam mode memory is conserved
  if(spp.single_beam_mode) num_beams=1;
  unsigned int bytes=0;
  int szarrayc3=size_array(sizeof(complex<float>),3, num_lons,num_beams,num_lats);
  int szarrayf3=size_array(sizeof(float),3, num_lons,num_beams,num_lats);
  int szarrayi3=size_array(sizeof(int),3, num_lons,num_beams,num_lats);
  int szarrayf2=size_array(sizeof(float),2, num_lons,num_lats);
  area_pixel=(float*)malloc(sizeof(float)*num_lats);
  float radii[3];
  radii[0] = default_target_radii[PositionVector::X].km();
  radii[1] = default_target_radii[PositionVector::Y].km();
  radii[2] = default_target_radii[PositionVector::Z].km();
  float rmag=sqrt(radii[0]*radii[0]+radii[1]*radii[1]+radii[2]*radii[2]);
  float pixresnom=(1.0/pixelsperdegree)*degtorad*rmag;
  pixresnom*=pixresnom;
  for(int j=0;j<num_lats;j++){
    float lat_in_rad=grid.latInRad(j);
    area_pixel[j]=cos(lat_in_rad)*pixresnom;
  }
  if(output_enable_[S0_UNCORR]){
    s0_uncorr_bfr=
      (float***)make_array(sizeof(complex<float>),3, num_lons,num_beams,num_lats);
    bytes+=szarrayc3;
  }
  s0_nsub_uncorr_bfr=
    (float***)make_array(sizeof(complex<float>),3, num_lons,num_beams,num_lats);
  bytes+=szarrayc3;

  s0ne_bfr=
    (float***)make_array(sizeof(complex<float>),3, num_lons,num_beams,num_lats);
  bytes+=szarrayc3;

  if(output_enable_[INC]){
    inc_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
    bytes+=szarrayf3;
  }

  if(output_enable_[DOPPLER]){
    dop_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
    bytes+=szarrayf3;
  }

  if(output_enable_[RANGE]){
    range_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
    bytes+=szarrayf3;
  }

  if(output_enable_[RANGE]||output_enable_[DOPPLER]
     ||output_enable_[NOMINAL_BURST_NUM]){
    max_X_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
    bytes+=szarrayf3;
  }

  X_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
  bytes+=szarrayf3;

  if(output_enable_[START_BURST_NUM]){
    start_burst_num=
      (int***)make_array(sizeof(int),3,num_lons,num_beams,num_lats); 
    bytes+=szarrayi3; 
  }
  if(output_enable_[END_BURST_NUM]){
    end_burst_num=(int***)make_array(sizeof(int),3,num_lons,num_beams,num_lats);    bytes+=szarrayi3; 
  }

  if(output_enable_[NOMINAL_BURST_NUM]){
    nom_burst_num=(int***)make_array(sizeof(int),3,num_lons,num_beams,num_lats);    bytes+=szarrayi3; 
  }


  num_looks=(int***)make_array(sizeof(int),3,num_lons,num_beams,num_lats);
  bytes+=szarrayi3; 

  slat_bfr = (double**)make_array(sizeof(double),2,num_lons,num_lats);
  slon_bfr = (double**)make_array(sizeof(double),2,num_lons,num_lats);
  bytes+=2*szarrayf2;

  bursts_since_update=(int*)make_array(sizeof(int),1,num_lons);
  beams_found=(bool**)make_array(sizeof(bool),2,num_lons,num_beams);
  checksum_ = (int *) make_array(sizeof(int), 1, BIDRTYPECOUNT);

  beammask_io_bfr = (char*) make_array(sizeof(char), 1, num_lats);
  s0_corr_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  lat_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  lon_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  range_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  dop_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  inc_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  X_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  std_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  s0ne_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  s0_uncorr_io_bfr = (float*) make_array(sizeof(float), 1, num_lats);
  start_burst_num_io_bfr = (int*) make_array(sizeof(int), 1, num_lats);
  nom_burst_num_io_bfr = (int*) make_array(sizeof(int), 1, num_lats);
  end_burst_num_io_bfr = (int*) make_array(sizeof(int), 1, num_lats);
  num_looks_io_bfr = (unsigned char*) make_array(sizeof(unsigned char), 1, num_lats);
  bytes+=num_lats*(sizeof(char)+7*sizeof(float)+4*sizeof(int));
  
  if(procMode==TOPO || procMode==TOPO_PLUS){
    dsigma0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    topo_inc1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    topo_inc2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    topo_lat=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    topo_lon=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    topo_wnl=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    mid_ds0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    slope_ds0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    rms_ds0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    noisefloor1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    noisefloor2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    power1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    power2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    ambrat1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    ambrat2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    dsigma0std=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    dsigma0bias=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    sumsigma0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    overlap_col=(int**)make_array(sizeof(int),2,4,MAX_BIDR_LINES);
    overlap_width=(int**)make_array(sizeof(int),2,4,MAX_BIDR_LINES);
    noise_equiv_s0=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats); 
    Xambig_bfr=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats); 
    resfactor=(float***)make_array(sizeof(float),3,num_lons,num_beams,num_lats);
    bytes+=szarrayf2*6 +szarrayf3*2;
  }
  if(procMode==TOPO_PLUS){
    topoplus_s01=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s02=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0ne1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0ne2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0baqscale1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0baqscale2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0amb1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    topoplus_s0amb2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
  }
  if(!s0_uncorr_bfr || !s0_nsub_uncorr_bfr || !s0ne_bfr || !X_bfr || !bursts_since_update
     || !beams_found || !num_looks || !slat_bfr || !slon_bfr || !checksum_
     || !beammask_io_bfr || !s0_corr_io_bfr || !lat_io_bfr || !lon_io_bfr
     || !s0_uncorr_io_bfr || !start_burst_num_io_bfr 
     || !end_burst_num_io_bfr || !num_looks_io_bfr || !X_io_bfr
     || !std_io_bfr || !s0ne_io_bfr
     || !range_io_bfr || !dop_io_bfr || !nom_burst_num_io_bfr || !area_pixel){
    cerr << "BIDR allocation failing for num_lons = " << num_lons
	 <<", num_lats = " << num_lats << ", num_beams = " << num_beams
         << endl;
    ErrorMessage e("BIDR object failed to allocate buffers ...\nTry reducing PIXELS_PER_DEGREE or using SINGLE_BEAM_MODE");
    e.throwMe();
  }
  if((procMode==NORMAL)&& (!start_burst_num || !end_burst_num || !inc_bfr)){
    cerr << "BIDR allocation failing for num_lons = " << num_lons
	 <<", num_lats = " << num_lats << ", num_beams = " << num_beams
         << endl;
    ErrorMessage e("BIDR object failed to allocate buffers ...\nTry reducing PIXELS_PER_DEGREE or using SINGLE_BEAM_MODE");
    e.throwMe();
  }

  if((procMode==STEREO)&& (!range_bfr || !dop_bfr || !nom_burst_num)){
    cerr << "BIDR allocation failing for num_lons = " << num_lons
	 <<", num_lats = " << num_lats << ", num_beams = " << num_beams
         << endl;
    ErrorMessage e("BIDR object failed to allocate buffers ...\nTry reducing PIXELS_PER_DEGREE or using SINGLE_BEAM_MODE");
    e.throwMe();
  }

  if((procMode==TOPO || procMode==TOPO_PLUS)&& (!dsigma0 || !dsigma0std || !dsigma0bias ||
			 !sumsigma0 || !overlap_col || !overlap_width ||
			 !noise_equiv_s0 || !Xambig_bfr || !resfactor || !topo_inc1 ||
			 !topo_inc2 || !topo_lat || !topo_lon || !topo_wnl ||
			 !mid_ds0 || !slope_ds0 || !rms_ds0 || !inc_bfr ||
			 !noisefloor1 || !noisefloor2 || !ambrat1 || !ambrat2 
			 || !power1 || !power2)){
    cerr << "BIDR allocation failing for num_lons = " << num_lons
	 <<", num_lats = " << num_lats << ", num_beams = " << num_beams
         << endl;
    ErrorMessage e("BIDR object failed to allocate buffers ...\nTry reducing PIXELS_PER_DEGREE or using SINGLE_BEAM_MODE");
    e.throwMe();
  }

  if((procMode==TOPO_PLUS)&& ( !topoplus_s01 || !topoplus_s02 ||
			       !topoplus_s0ne1 || !topoplus_s0ne2 ||
			       !topoplus_s0baqscale1 || !topoplus_s0baqscale2 ||
			       !topoplus_s0amb1 || !topoplus_s0amb2)
     ){
    cerr << "BIDR TOPO_PLUS allocation failing for num_lons = " << num_lons
	 <<" MAX_TOPO_OVERLAPWIDTH=" << MAX_TOPO_OVERLAP_WIDTH
         << endl;
    ErrorMessage e("BIDR object failed to allocate buffers ...\nTry reducing PIXELS_PER_DEGREE or using SINGLE_BEAM_MODE");
    e.throwMe();
  }

  // initialize arrays
  for(int i=0;i<num_lons;i++){
    for(int b=0;b<num_beams;b++){
      for(int j=0;j<num_lats;j++){
	if(output_enable_[START_BURST_NUM])start_burst_num[i][b][j]=-1;
        if(output_enable_[END_BURST_NUM])end_burst_num[i][b][j]=-1;
        if(output_enable_[INC])inc_bfr[i][b][j]=0;
	num_looks[i][b][j]=0;
        if(procMode==STEREO) max_X_bfr[i][b][j]=0;
        if(procMode==TOPO || procMode==TOPO_PLUS){
	  resfactor[i][b][j]=0;
	  noise_equiv_s0[i][b][j]=0;
	  Xambig_bfr[i][b][j]=0;
	}
        X_bfr[i][b][j]=0;
        s0_uncorr_bfr[i][b][j]=0;
        s0_nsub_uncorr_bfr[i][b][j]=0;
        s0ne_bfr[i][b][j]=0;
      }      
    }
  }
  if(procMode==TOPO || procMode==TOPO_PLUS){
    for(int i=0;i<4;i++){
      for(int j=0;j<MAX_BIDR_LINES;j++){
	dsigma0[i][j]=0;
	sumsigma0[i][j]=0;
	dsigma0std[i][j]=0;
	dsigma0bias[i][j]=0;
	overlap_col[i][j]=0;
	overlap_width[i][j]=0;
	topo_inc1[i][j]=0;
	topo_inc2[i][j]=0;
	topo_lat[i][j]=0;
	topo_lon[i][j]=0;
	topo_wnl[i][j]=0;
	mid_ds0[i][j]=0;
	slope_ds0[i][j]=0;
	rms_ds0[i][j]=0;
	power1[i][j]=0;
	power2[i][j]=0;
	ambrat1[i][j]=0;
	ambrat2[i][j]=0;
	noisefloor1[i][j]=0;
	noisefloor2[i][j]=0;
      }
    }
  }

 if(procMode==TOPO_PLUS){
    for(int i=0;i<4;i++){
      for(int j=0;j<MAX_BIDR_LINES;j++){
	for(int k=0;k<MAX_TOPO_OVERLAP_WIDTH;k++){
	  topoplus_s01[i][j][k]=BAD_VALUE;
	  topoplus_s02[i][j][k]=BAD_VALUE;
	  topoplus_s0ne1[i][j][k]=BAD_VALUE;
	  topoplus_s0ne2[i][j][k]=BAD_VALUE;
	  topoplus_s0baqscale1[i][j][k]=BAD_VALUE;
	  topoplus_s0baqscale2[i][j][k]=BAD_VALUE;
	  topoplus_s0amb1[i][j][k]=BAD_VALUE;
	  topoplus_s0amb2[i][j][k]=BAD_VALUE;
	}
      }
    }
  }

  for(int i=0;i<num_lons;i++){
      bursts_since_update[i]=0;
      for(int b=0;b<num_beams;b++){
	beams_found[i][b]=false;
      }
      for(int j=0;j<num_lats;j++){
        slat_bfr[i][j]=0;
        slon_bfr[i][j]=0;
      }
  }

  for(int i=0;i<BIDRTYPECOUNT;i++){
    checksum_[i] = 0;
  }


  if(dbg.level){
    dbg.file << "BIDR allocates " << bytes 
	     << " bytes of array + " << sizeof(BIDR) << " bytes of object = ";
    bytes+=sizeof(L1I);
    dbg.file << bytes << " bytes Total" << endl;
    dbg.file << "  " << "num_lons=" << num_lons <<", num_lats=" << num_lats
	     << ", num_beams=" << num_beams << endl;
  }
}    

void
BIDR::smoothOverlap(int**col, int** width){
  // set up arrays for fits

#define DEBUG_SMOOTH
#ifdef DEBUG_SMOOTH
  FILE* dfp = fopen("dbgsmooth.txt","w");
#endif  
  for(int i=0;i<4;i++){

    // compute means and stds to normalize data for fit
    double line_mean=0, line_std=0;
    double col_mean=0, col_std=0;
    double width_mean=0, width_std=0;
    int line_min=0, line_max=0;
    int numsamp=0;
    bool min_found=false;
    for(int j=0;j<MAX_BIDR_LINES;j++){
      // check valid overlap
      if(width[i][j]!=-1){
         if(!min_found){
	   line_min=j;
           min_found=true;
	 }
	 line_max=j;
         numsamp++;
         line_mean+=j;
         line_std+=j*j;
         col_mean+=col[i][j];
         col_std+=col[i][j]*col[i][j];
         width_mean+=width[i][j];
         width_std+=width[i][j]*width[i][j];
      }
    }
    line_mean=line_mean/numsamp;
    col_mean=col_mean/numsamp;
    width_mean=width_mean/numsamp;
    line_std=line_std-numsamp*line_mean*line_mean;
    col_std=col_std-numsamp*col_mean*col_mean;
    width_std=width_std-numsamp*width_mean*width_mean;
    line_std=sqrt(line_std/(numsamp-1));
    col_std=sqrt(col_std/(numsamp-1));
    width_std=sqrt(width_std/(numsamp-1));

  // compute normalize data for fits
    int N=line_max-line_min+1;


    Dvec x("x",N);
    Dvec yc("yc",N);
    Dvec yw("yw",N);
    for(int j=0;j<N;j++){
      x(j)=((line_min+j)-line_mean)/line_std;
      yc(j)=(col[i][line_min+j]-col_mean)/col_std;
      yw(j)=(width[i][line_min+j]-width_mean)/width_std;
    }
  // perform 15th order polynomial fit on overlap and col 
  // do fits from first valid (nonzero col, width>-1) to last 
  // exclude zero values from fit.
    int porder=12;
    Dvec pc("pc",porder+1);
    Dvec pw("pw",porder+1);
    polyfit(pc,yc,x);
    polyfit(pw,yw,x);

#ifdef DEBUG_SMOOTH
    fprintf(dfp,"pw=[");
    for(int j=0;j<porder+1;j++) fprintf(dfp,"%g ",pw(j));
    fprintf(dfp,"]\n");
    fprintf(dfp,"pc=[");
    for(int j=0;j<porder+1;j++) fprintf(dfp,"%g ",pc(j));
    fprintf(dfp,"]\n");
    fprintf(dfp,"width_mean=%g width_std=%g col_mean=%g col_std=%g line_mean=%g line_std=%g line_min=%d line_max=%d\n",width_mean,width_std,col_mean,col_std,line_mean,line_std,line_min,line_max);
#endif
    double dcol[N];
    double dwidth[N];


    for(int j=0;j<N;j++){
      dcol[j]=0;
      dwidth[j]=0;
      for(int k=porder;k>0;k--){
	dcol[j]+=pc(k);
        dcol[j]*=x(j);
	dwidth[j]+=pw(k);
        dwidth[j]*=x(j);
      }
      dcol[j]+=pc(0);
      dwidth[j]+=pw(0);
#ifdef DEBUG_SMOOTH
      fprintf(dfp,"FIT %d %d %g %g %g %g %g\n",i,j,x(j),yc(j),dcol[j],yw(j),dwidth[j]);
#endif
      dcol[j]*=col_std;
      dcol[j]+=col_mean;

      dwidth[j]*=width_std;
      dwidth[j]+=width_mean;
    }


    // replace values with fit rounded to nearest nonnegative integer 

    for(int j=0;j<MAX_BIDR_LINES;j++){

 
      if(j<line_min){
	col[i][j]=(int)(dcol[0]+0.5);
	width[i][j]=-1;
      }
      else if(j>line_max){
	col[i][j]=(int)(dcol[N-1]+0.5);
	width[i][j]=-1; 
      }

      else{
#ifdef DEBUG_SMOOTH
        fprintf(dfp,"%d %d %d %d ",i,j,col[i][j],width[i][j]);
#endif
	col[i][j]=(int)(dcol[j-line_min]+0.5);
	width[i][j]=(int)(dwidth[j-line_min]*1+0.5);

#ifdef DEBUG_SMOOTH
        fprintf(dfp,"%d %d\n",col[i][j],width[i][j]);
#endif
      }
  
     
    } // end loop over j
  } // end loop over i
}

void 
BIDR::setOverlap(int** col, int** width){
  overlap_set=true;
  
  if(procMode!=TOPO && procMode!=TOPO_PLUS){
    ErrorMessage e("BIDR::setOverlap only works in TOPO and TOPO_PLUS modes");
    e.throwMe();
  }
 

  // set fixed overlap for candidate height SAR processor runs
  for(int i=0;i<4;i++){
    for(int j=0;j<MAX_BIDR_LINES;j++){
      overlap_col[i][j]=col[i][j];
      overlap_width[i][j]=width[i][j];
      if( procMode==TOPO_PLUS && 2*overlap_width[i][j]+1 > MAX_TOPO_OVERLAP_WIDTH){
	cerr << "In BIDR.h MAX_TOPO_OVERlAP_WIDTH = " << MAX_TOPO_OVERLAP_WIDTH << endl;
        fprintf(stderr,"This is smaller than width for profile # %d line %d which is %d\n",i,j,2*width[i][j]+1);
        exit(1);
      }
    }
  }
}
void
BIDR::updateTopoArrays(int first_i){
  float noise_floor1[4]={0,0,0,0}, noise_floor2[4]={0,0,0,0};
  float s01[4]={0,0,0,0},s02[4]={0,0,0,0};
  float SNR1[4]={0,0,0,0},SNR2[4]={0,0,0,0};
  float EY[4]={0,0,0,0}, EXY[4]={0,0,0,0}, EXX[4]={0,0,0,0};
  float MSS[4]={0,0,0,0}; // mean sum square s01-s02;
  float P1[4]={0,0,0,0}, P2[4]={0,0,0,0}; // mean powers
  float sumX1[4]={0,0,0,0}, sumX2[4]={0,0,0,0}; // sums of X factor
  float sumXamb1[4]={0,0,0,0}, sumXamb2[4]={0,0,0,0}; // sums of Xambiguity
  int N[4]={0,0,0,0};
  float truenumlooks1[4]={0,0,0,0},truenumlooks2[4]={0,0,0,0};
  
  // estimate intermediate arrays dsigma0, overlap_col and overlap_width
  for(int b=0;b<4;b++){
    bool found_overlap=false;
    int j0=0,j1=0;
    float nl1,nl2;
    int num_lats=grid.numLats();
//#define DEBUGTOPO
    for(int j=0;j<num_lats;j++){
      float inc1, inc2,s01corr,s02corr, s0ne1, s0ne2;

      if(num_looks[first_i][b][j]>=num_looks_required){
	nl1=num_looks[first_i][b][j]*resfactor[first_i][b][j];
      }
      else{
	nl1=0;
      }
      if(num_looks[first_i][b+1][j]>=num_looks_required){
	nl2=num_looks[first_i][b+1][j]*resfactor[first_i][b+1][j];
      }
      else{
	nl2=0;
      }

      // only use data at least 3 dB above the noise floor
      if(!overlap_set){
	if(s0_uncorr_bfr[first_i][b][j]<2*noise_equiv_s0[first_i][b][j])
	   nl1=0;
	if(s0_uncorr_bfr[first_i][b+1][j]<2*noise_equiv_s0[first_i][b+1][j])
	   nl2=0;
      }

      // if overlap is set only use pixels in range
      else{
        int col=overlap_col[b][num_lines_of_output];
        int width=overlap_width[b][num_lines_of_output];
	if(j<col-width || j>col+width){
	  nl1=0;
	  nl2=0;
	}
      }
      if(nl1>0 && nl2>0){
	if(!found_overlap){
	  found_overlap=true;
	  j0=j;
	}
	j1=j;

        inc1=inc_bfr[first_i][b][j];
        inc2=inc_bfr[first_i][b+1][j];
	s01corr=incAngleCorrect(s0_nsub_uncorr_bfr[first_i][b][j],inc1,0,0,inc_correction_model);
	s02corr=incAngleCorrect(s0_nsub_uncorr_bfr[first_i][b+1][j],inc2,0,0,inc_correction_model);
	s01[b]+=s01corr*nl1;
        if(isnan(s01[b])){
	  printf("BIDR::UpdateTopoArrays Found NaN in s01[%d], nl1=%g s01corr=%g s0_nsub_uncorr[i][b][j]=%g inc1=%g i=%d b=%d j=%d\n",b,nl1,s01corr,s0_nsub_uncorr_bfr[first_i][b][j],inc1,first_i,b,j);
	  exit(1);
	}
    
        s02[b]+=s02corr*nl2;

        P1[b]+=X_bfr[first_i][b][j]*s0_uncorr_bfr[first_i][b][j]*nl1;
        P2[b]+=X_bfr[first_i][b+1][j]*s0_uncorr_bfr[first_i][b+1][j]*nl2;

        sumX1[b]+=X_bfr[first_i][b][j]*nl1;
        sumX2[b]+=X_bfr[first_i][b+1][j]*nl2;

        sumXamb1[b]+=Xambig_bfr[first_i][b][j]*nl1;
        sumXamb2[b]+=Xambig_bfr[first_i][b+1][j]*nl2;

        float x=j-overlap_col[b][num_lines_of_output];
        float y= s01corr-s02corr;
	EY[b]+=y;
        MSS[b]+=y*y;
        EXY[b]+=x*y;
	EXX[b]+=x*x;
	N[b]++;
	truenumlooks1[b]+=nl1;
	truenumlooks2[b]+=nl2;
	s0ne1=incAngleCorrect(noise_equiv_s0[first_i][b][j],inc1,0,0,inc_correction_model);
	s0ne2=incAngleCorrect(noise_equiv_s0[first_i][b+1][j],inc2,0,0,inc_correction_model);
        noise_floor1[b]+=s0ne1*nl1;
        noise_floor2[b]+=s0ne2*nl2;
#ifdef DEBUGTOPO
    if(num_lines_of_output>880 && num_lines_of_output<920 && b==0)
      printf("line %d overlap %d j0=%d j1=%d s01=%g s02=%g tnl1=%g tnl2=%g nf1=%g nf2=%g\n",num_lines_of_output,b,
	     j0,j1,s01[b],s02[b],truenumlooks1[b],truenumlooks2[b],
	     noise_floor1[b],noise_floor2[b]);
#endif
      }
    } // end column (latitude)loop
    s01[b]/=truenumlooks1[b];
    s02[b]/=truenumlooks2[b];
    noise_floor1[b]/=truenumlooks1[b];
    noise_floor2[b]/=truenumlooks2[b];
    P1[b]/=truenumlooks1[b];
    P2[b]/=truenumlooks1[b];
    SNR1[b]=s01[b]/noise_floor1[b];
    SNR2[b]=s02[b]/noise_floor2[b];
    MSS[b]=MSS[b]/N[b];
    EY[b]=EY[b]/N[b];
    EXY[b]=EXY[b]/N[b];
    EXX[b]=EXX[b]/N[b];

    if(! overlap_set){
      overlap_col[b][num_lines_of_output]=(j1+j0)/2;
      overlap_width[b][num_lines_of_output]=(j1-j0)/2;
    }
    if(truenumlooks1[b]==0) overlap_width[b][num_lines_of_output]=0;
#ifdef DEBUGTOPO
    if(num_lines_of_output>880 && num_lines_of_output<920 && b==0)
      printf("Final normed line=%d s01=%g s02=%g nf1=%g nf2=%g SNR1=%g SNR2=%g col=%d width=%d\n",num_lines_of_output,
	     s01[b],s02[b],noise_floor1[b],noise_floor2[b],SNR1[b],SNR2[b],
	     overlap_col[b][num_lines_of_output],
	     overlap_width[b][num_lines_of_output]);
#endif
  } // end beam loop
  // compute remaining topo arrays
  for(int i=0;i<4;i++){
    if(truenumlooks1[i]>0){
      dsigma0[i][num_lines_of_output]=s01[i]-s02[i];
      sumsigma0[i][num_lines_of_output]=s01[i]+s02[i];
      dsigma0std[i][num_lines_of_output]=sqrt((1+2/SNR1[i]+1/(SNR1[i]*SNR1[i]))/truenumlooks1[i]+(1+2/SNR2[i]+1/(SNR2[i]*SNR2[i]))/truenumlooks2[i]);
      dsigma0bias[i][num_lines_of_output]=0.2*noise_floor1[i]+0.2*noise_floor2[i];
      topo_inc1[i][num_lines_of_output]=inc_bfr[first_i][i][overlap_col[i][num_lines_of_output]];
      topo_inc2[i][num_lines_of_output]=inc_bfr[first_i][i+1][overlap_col[i][num_lines_of_output]];
      topo_lat[i][num_lines_of_output]=slat_bfr[first_i][overlap_col[i][num_lines_of_output]];
      topo_lon[i][num_lines_of_output]=slon_bfr[first_i][overlap_col[i][num_lines_of_output]];
      topo_wnl[i][num_lines_of_output]=(truenumlooks1[i]+truenumlooks2[i])/2.0;
      mid_ds0[i][num_lines_of_output]=EY[i];
      slope_ds0[i][num_lines_of_output]=EXY[i]/EXX[i];
      rms_ds0[i][num_lines_of_output]=sqrt(MSS[i]);
      noisefloor1[i][num_lines_of_output]=noise_floor1[i];
      noisefloor2[i][num_lines_of_output]=noise_floor2[i];
      ambrat1[i][num_lines_of_output]=sumXamb1[i]/sumX1[i];
      ambrat2[i][num_lines_of_output]=sumXamb2[i]/sumX2[i];
      power1[i][num_lines_of_output]=P1[i];
      power2[i][num_lines_of_output]=P2[i];
    }
    else{
      dsigma0[i][num_lines_of_output]=BAD_VALUE;
      sumsigma0[i][num_lines_of_output]=BAD_VALUE;
      dsigma0std[i][num_lines_of_output]=BAD_VALUE;
      dsigma0bias[i][num_lines_of_output]=BAD_VALUE;
      topo_inc1[i][num_lines_of_output]=BAD_VALUE;
      topo_inc2[i][num_lines_of_output]=BAD_VALUE;
      topo_lat[i][num_lines_of_output]=BAD_VALUE;
      topo_lon[i][num_lines_of_output]=BAD_VALUE;
      topo_wnl[i][num_lines_of_output]=BAD_VALUE;
      mid_ds0[i][num_lines_of_output]=BAD_VALUE;
      slope_ds0[i][num_lines_of_output]=BAD_VALUE;
      rms_ds0[i][num_lines_of_output]=BAD_VALUE;
      noisefloor1[i][num_lines_of_output]=BAD_VALUE;
      noisefloor2[i][num_lines_of_output]=BAD_VALUE;
      power1[i][num_lines_of_output]=BAD_VALUE;
      power2[i][num_lines_of_output]=BAD_VALUE;
      ambrat1[i][num_lines_of_output]=BAD_VALUE;
      ambrat2[i][num_lines_of_output]=BAD_VALUE;
    }
#ifdef DEBUGTOPO
    if(num_lines_of_output>880 && num_lines_of_output<920 && i==0)
      printf("Final dsigma0=%g sums0=%g dsigma0std=%g dsigma0bias=%g\n",
	     dsigma0[i][num_lines_of_output],
	     sumsigma0[i][num_lines_of_output],
	     dsigma0std[i][num_lines_of_output],
	     dsigma0bias[i][num_lines_of_output]);
#endif
  }
}
void
BIDR::topoProcessAndWrite(L1I& l1i, SARProcParams& spp, 
			  bool single_burst_mode, float trial_height){
  // for now topo outputs a subset of NOMINAL backplanes so I
  // can merely call stripProcessAnd Write
  // this will need to change if I wind up outputting more
  stripProcessAndWrite(l1i,spp,single_burst_mode);
}
      
void
BIDR::stripProcessAndWrite(L1I& l1i, SARProcParams& spp, bool single_burst_mode){

  // set the direction to use for the X test when choosing a beam
  switch(l1i.replace_X_param){
  case L1I::NONE:
    Xcheckdir=1;
    break;
  case L1I::RANGE_RES:
    Xcheckdir=-1;
    break;
  case L1I::AZIMUTH_RES:
    Xcheckdir=-1;
    break;
  case L1I::LOOKVEL_ANGLE:
    Xcheckdir=1;
    break;
  case L1I::RES_ELEMENT_ANGLE:
    Xcheckdir=1;
    break;
  case L1I::AMBIG_RATIO:
    Xcheckdir=-1;
    break;
  default:
    Xcheckdir=1;
    break;
  }
  static bool first_call=true;
  static DebugInfo dbg("BIDR::stripProcessAndWrite");


  bool bad_margin=false; // flag to determine if usable area is clipped.
 
  // get burst and beam numbers
  int burst_no=l1i.sab_counter;
  int bn=l1i.beam_number-1;

  // necessary to save memory alloaction for single beam mode
  if(spp.single_beam_mode) bn=0;


  // determine value to replace X with if desired
  if(replace_X_mode){
    replace_X_value=l1i.getParam(replace_X_param);
  }

  // determine loose boundaries on grid for burst image
  // Geometry conversion from (t,range, doppler) Standard(slat,slon) 
  // performed by L1i  from Standard (slat,slon) to OblCyl(lat,lon) is
  // performed by the OblCylProj object
  l1i.computeLatLonBounds(proj); 

  // sets up simple model for converting lat/lon to/from OblCyl lat/lon
  // that only applies to burst
  proj.initBurstTransformation(); 

  int bounds_ij[4];

  // Get absolute boundary indices (w/o performing modulus to restrict to
  // grid)
  // If the boundaries extend beyond the edge of the current
  grid.getBoundary(l1i.lonlat_bounds, bounds_ij);

  int i_min=bounds_ij[0];
  int i_max=bounds_ij[1];
  int j_min=bounds_ij[2];
  int j_max=bounds_ij[3];

  if(dbg.level){
    dbg.file << "BIDR::stripProcessAndWrite:Burst #" << burst_no << " Beam #" << l1i.beam_number << endl;
  }

  if(dbg.level){
    dbg.file << "BIDR::stripProcessAndWrite:Lon Bounds (" << l1i.lonlat_bounds[0]*radtodeg <<","
	     << l1i.lonlat_bounds[1]*radtodeg  << ")" << endl;
    dbg.file << "BIDR::stripProcessAndWrite:LonGrid  Absolute Bounds (" << i_min <<","
	     << i_max  << ") and Relative Bounds (" << grid.relLonIndex(i_min) <<","
	     << grid.relLonIndex(i_max)  << ")" << endl;

  }

  // check to see if burst longitude bounds are in valid range
  bool minvalid =grid.checkValidLongitude(i_min);
  bool maxvalid =grid.checkValidLongitude(i_max);
  
  // Make sure both min and max longitudes are valid first time through
  // if not throw error
  if((!minvalid || !maxvalid) && first_call && !single_burst_mode && !ignore_gaps_mode){
    float bufl1 = grid.absStartValidIndex()/(float)pixelsperdegree;
    float bufl2 = grid.absEndValidIndex()/(float)pixelsperdegree;
    float bl1 = i_min/(float)pixelsperdegree;
    float bl2 = i_max/(float)pixelsperdegree;
    ErrorMessage e("BIDR::StripProcessAndWrite  Out of Bounds Longitude Line Needed for first burst (Burst#"+toStr(burst_no)+") \nTry Increasing BIDR_LATLON_PAD_EACH_SIDE OR BIDR_LON_BUFFER_SIZE (default 9 deg)\n Burst Lon Range=("+toStr(bl1)+","+toStr(bl2)+") Buffer lon range=("+toStr(bufl1)+","+toStr(bufl2)+") lon_inc is "+toStr(grid.lonIncreasing()));
    first_call=false;
    e.throwMe();
  } 
  first_call=false;



  // Check to see if the buffer is full and we need to output line(s)
  // output lines as necessary
  while(!minvalid || !maxvalid){


  // Check to see if we need a line that has already been written
  //  (throws error)
    if( !single_burst_mode && !ignore_gaps_mode && ((grid.lonIncreasing() && !minvalid) || 
	(!grid.lonIncreasing() && !maxvalid))){
      ErrorMessage e("BIDR::stripProcessAndWrite  Previously Written Longitude Line Needed for Burst#"+toStr(burst_no)+"\nTry increasing BIDR_LON_BUFFER_SIZE default 9.0 degrees, needs to be a factor of 360");
      e.throwMe();
    }
    if(dbg.level){
      dbg.file << "Minvalid " << minvalid << " Maxvalid " << maxvalid << " outputting line and shifting buffer ..." << endl;
    }
    outputLine(spp);
    minvalid =grid.checkValidLongitude(i_min);
    maxvalid =grid.checkValidLongitude(i_max);
  }

  // output bounding box to standard info line
  spp.cout_line+=toStr(j_min);
  spp.cout_line+=" ";
  spp.cout_line+=toStr(j_max);
  spp.cout_line+=" ";
  spp.cout_line+=toStr(grid.absIndexToLineNumber(i_min));
  spp.cout_line+=" ";
  spp.cout_line+=toStr(grid.absIndexToLineNumber(i_max));
  spp.cout_line+=" ";

  int valid_pixels_found=0;
  // loop over grid points within boundaries

  double burst_resolution=l1i.resolutionArea(spp);
  for(int abs_i=i_min;abs_i<=i_max;abs_i++){


    int i=grid.relLonIndex(abs_i);
 
    // set bursts_since_update to 0 
    // A longitude Line is still needed if it is within the LatLon boundary of a burst
    // even if it is not in the useable area
    bursts_since_update[i]=0;

    // indicate beam bn+1 has been found for this longitude line.
    // for single beam mode this is unnecessary but harmless.
    beams_found[i][bn]=true;
    

    for(int j=j_min;j<=j_max;j++){

      // compute doppler and range at point
      double fdop_in_Hz, range_in_km;

      float slat = slat_bfr[i][j];
      float slon = slon_bfr[i][j];

      // get doppler and range for pixel
      // may change this to operate on projected lat and lon for speed
      l1i.latLonToDopplerAndRange(slon,slat, fdop_in_Hz,
			      range_in_km);


      // check to see if pixel is on same side of nadir as boresight.
      bool goodlook = l1i.goodSideOfNadir(slon,slat);
      bool inbounds = false;


//#define DEBUG_FIX_REROUTE
#ifdef  DEBUG_FIX_REROUTE
      if(fabs(range_in_km-57105.4)< 1 && fabs(fdop_in_Hz-764348)<2 && goodlook){
	cout << range_in_km << " " << fdop_in_Hz << " " << goodlook << " " << inbounds << endl;
      }
#endif      

      // check to see if pixel is in useable area and initialize Interpolation
      if(goodlook) inbounds=l1i.initInterpolate(fdop_in_Hz, range_in_km);



      bool max_looks_exceeded = enable_max_looks && 
	(num_looks[i][bn][j]>=max_looks_allowed+num_looks_to_skip);

      if(inbounds) valid_pixels_found++;
      // if doppler/range within processor window
      if(inbounds && !max_looks_exceeded){
	
        // if edges are useable we do not have enough margin
        // throw an error
	if(j==j_min || j==j_max || abs_i==i_min || abs_i== i_max){
          bad_margin=true;
 	}

        // this code allows the first N looks of a pixel to be skipped


	// add to num_looks
	num_looks[i][bn][j]++;

	if(num_looks[i][bn][j]<num_looks_to_skip){
	  // this code allows the first N looks of a pixel to be skipped
	  continue;
	}



              // interpolate s0 to point, add to s0_uncorr_bfr[beamnum-1]
	float s0=l1i.sigma0Interpolate(fdop_in_Hz, range_in_km);
        if(mark_nadir_amb){
	  if(l1i.nadirAmb(fdop_in_Hz,range_in_km)){
	    s0=-1000;
	  }
	}
	s0_uncorr_bfr[i][bn][j]+=s0;
               // interpolate X to point, add to X_bfr[beamnum-1]
	float X = l1i.XInterpolate();
        float s0ns= (s0*X-spp.quant_noise_offset-spp.thermal_noise_offset*spp.quant_noise_scale)/spp.quant_noise_scale;
        
        s0ns/=X;
        float s0ne= spp.quant_noise_offset/(X*spp.quant_noise_scale)+spp.thermal_noise_offset/X;
        s0_nsub_uncorr_bfr[i][bn][j]+=s0ns;
        s0ne_bfr[i][bn][j]+=s0ne;
        

        if(procMode==STEREO){
	  if(max_X_bfr[i][bn][j]<X){
            max_X_bfr[i][bn][j]=X;
	    nom_burst_num[i][bn][j]=burst_no;
            range_bfr[i][bn][j]=l1i.rawRange(fdop_in_Hz,range_in_km);
            dop_bfr[i][bn][j]=l1i.rawDoppler(fdop_in_Hz);
	  }
	}
	else{
	  // compute incidence angle at point, add to inc_bfr[beam_num-1]
	  float inc = l1i.getIncidenceAngle(slat,slon); 
	  inc_bfr[i][bn][j] += inc;
	}
        if(procMode==TOPO || procMode==TOPO_PLUS){
          
	  resfactor[i][bn][j]+=area_pixel[j]/burst_resolution;
	  noise_equiv_s0[i][bn][j]+=(spp.thermal_noise_offset*spp.quant_noise_scale+spp.quant_noise_offset)/X/spp.quant_noise_scale; 
	  Xambig_bfr[i][bn][j]+=l1i.XAmbigInterpolate();
	}

        if(procMode==NORMAL){
	  // update start and end_burst_num
	  if(start_burst_num[i][bn][j]<0) start_burst_num[i][bn][j]=burst_no;
	  end_burst_num[i][bn][j]=burst_no;
	
	}

        if(replace_X_mode) X=replace_X_value;
	X_bfr[i][bn][j]+=X;


        if(dbg.level>1 && j==dbg.level){
	  dbg.file << "Pixel abs_i=" << abs_i << " i=" << i
		   << " j=" << j << " bn=" << bn 
		   << " num_looks=" << num_looks[i][bn][j] << endl; 
	}
      }
    }
  }

  // warn if no valid pixels were found
  if(valid_pixels_found==0){
    cerr << endl;
    cerr << "Warning BIDR::stripProcessAndWrite found no valid pixels for burst# " << burst_no << endl;
  }


  if(dbg.level){
    dbg.file << "BIDR::stripProcessAndWrite:Burst #" << burst_no <<  " " << valid_pixels_found << " Valid pixels were found."<< endl;
  }

  // warn if usable area is clipped in processing
  if(bad_margin){
    cerr << endl;
    cerr << "Warning: BIDR::stripProcessAndWrite: Not enough margin for burst# "
      +toStr(burst_no) << endl; 
    cerr << 
      "Warning(cont.): Examine L1I::computeLatLonBounds and increase the BURST_PAD_FACTOR" << endl;
    cerr << "Warning(cont) L1I::boundaryReport follows ...." << endl;
    l1i.reportBoundary(proj);
  }


  // add 1 to all valid bursts_since_update data.
  incrementBurstsSinceUpdate();
    
} 

// output oldest Longitude Line to BIDR files
void BIDR::outputLine(const SARProcParams& spp){

  static DebugInfo dbg("BIDR::outputLine");

  // get oldest index
  int first_i=grid.firstValidLongitudeIndex();

  unsigned int numlats=grid.numLats();


  // No output in the event that BIDR region forcing is employed 
  // and we are outside the forced region
  if(!spp.force_bidr_region || (num_lines_of_output >= spp.forced_first_line
		       && num_lines_of_output <= spp.forced_last_line)){

   if(output_beam_deltas) outputBeamDeltas(spp);
   if(dbg.level){
      dbg.file << "BIDR::outputLine" 
	       << num_lines_of_output+1
	       << " out of " << grid.totalNumLons() 
               << " Relative index =" << first_i 
               << " Absolute index =" << grid.absStartValidIndex() << endl;
      float lon_in_rad=grid.lineNumberToLonInRad(num_lines_of_output);
      float lat_in_rad=grid.latInRad(grid.numLats()/2);
      float slat = proj.standardLatInRad(lon_in_rad,lat_in_rad);
      float slon = proj.standardLonInRad(lon_in_rad,lat_in_rad);
      dbg.file << "Oblique Lon,Lat of center pixel is: ( "
	       << lon_in_rad*radtodeg << "," << lat_in_rad*radtodeg << endl;
      dbg.file << "PE Lon,Lat of center pixel is: ( "
	       << slon*radtodeg << "," << slat*radtodeg << endl; 
      
   }
      // for each latitude pixel in line
  for(unsigned int j=0;j<numlats;j++){
    // foreach beam (NUM_REQUIRED_LOOKS must be obtained for each beam
    //               regardless of feathering)
    for(int b=0;b<num_beams;b++){
      num_looks[first_i][b][j]-=num_looks_to_skip;
      if(dbg.level>1 && j==(unsigned int)dbg.level){
	dbg.file << "Pixel first_i=" << first_i << " j=" << j
		 << " b=" << b << " num_looks=" << num_looks[first_i][b][j]
		 << endl;
      }
      
      if(num_looks[first_i][b][j]<num_looks_required 
	 || isnan(s0_uncorr_bfr[first_i][b][j])){
	s0_uncorr_bfr[first_i][b][j]=BAD_VALUE;
	s0_nsub_uncorr_bfr[first_i][b][j]=BAD_VALUE;
	s0ne_bfr[first_i][b][j]=BAD_VALUE;
	if(output_enable_[INC])inc_bfr[first_i][b][j]=BAD_VALUE;
	X_bfr[first_i][b][j]=BAD_VALUE;
	if(output_enable_[RANGE])range_bfr[first_i][b][j]=BAD_VALUE;
	if(output_enable_[DOPPLER])dop_bfr[first_i][b][j]=BAD_VALUE;
	if(output_enable_[NOMINAL_BURST_NUM])nom_burst_num[first_i][b][j]=-1;
	if(output_enable_[START_BURST_NUM])start_burst_num[first_i][b][j]=-1;
	if(output_enable_[END_BURST_NUM])end_burst_num[first_i][b][j]=-1;
	if(procMode==TOPO || procMode==TOPO_PLUS){
	  resfactor[first_i][b][j]=0;
	  noise_equiv_s0[first_i][b][j]=BAD_VALUE;
	  Xambig_bfr[first_i][b][j]=BAD_VALUE;
	}
      }
      // normalize pixel for each beam (divide by num_looks)
      else{
	X_bfr[first_i][b][j]/=num_looks[first_i][b][j];
	s0_uncorr_bfr[first_i][b][j]/=num_looks[first_i][b][j];
	s0_nsub_uncorr_bfr[first_i][b][j]/=num_looks[first_i][b][j];
	s0ne_bfr[first_i][b][j]/=num_looks[first_i][b][j];
	if(output_enable_[INC])inc_bfr[first_i][b][j]/=num_looks[first_i][b][j];        if(procMode==TOPO || procMode==TOPO_PLUS){
	  resfactor[first_i][b][j]/=num_looks[first_i][b][j];
	  noise_equiv_s0[first_i][b][j]/=num_looks[first_i][b][j];
	  Xambig_bfr[first_i][b][j]/=num_looks[first_i][b][j];
	}
      }
    }  // End of beam loop
    
    // merge multiple beams and computed beam mask
    // may handle multiple strategies using up to
    // all of longitude line X values for each beam 
    // and parameters in spp object
    // Appropriately handles BAD_VALUES;
    // perform beam to beam weighting and
    // compute beam mask
	  
    unsigned char beam_mask=0;
    float s0_uncorr=BAD_VALUE, X=BAD_VALUE, s0ns=BAD_VALUE, s0ne=BAD_VALUE;
    float inc=BAD_VALUE;
    int burst_num_1=0,burst_num_2=0;
    int numlooks=0;
    float range=BAD_VALUE,doppler=BAD_VALUE;
    multipleBeamMerge(spp,first_i,j,beam_mask,s0_uncorr,inc,burst_num_1,
		      burst_num_2,numlooks,X,s0ns,s0ne,doppler,range);

    // update topoplus arrays
    if(procMode==TOPO_PLUS & overlap_set){
      // determine relevant indices into topoplus arrays
      // (place in overlap and profile number)
      for(int pn=0;pn<4;pn++){
	if((int)j<=overlap_col[pn][num_lines_of_output]+overlap_width[pn][num_lines_of_output] &&
	   (int)j>=overlap_col[pn][num_lines_of_output]-overlap_width[pn][num_lines_of_output] && 
	   overlap_width[pn][num_lines_of_output] > 0){
	  int j2=(int)j-overlap_col[pn][num_lines_of_output]+overlap_width[pn][num_lines_of_output];
	  updateTopoPlusArrays(spp,first_i,j,pn,num_lines_of_output,j2);
	}
      }
    }

    // fetch latitude and longitude in standard projection
    // in radians
    float slat = slat_bfr[first_i][j];
    float slon = slon_bfr[first_i][j];
    
    // compute incidence angle corrected sigma0
    // appropriately handle BAD_VALUES
    float s0_corr=BAD_VALUE;
    if(output_enable_[S0_CORR]){
      if(spp.noise_subtraction_on)
	s0_corr = incAngleCorrect(s0ns,inc,slat,slon,inc_correction_model);
      else
	s0_corr = incAngleCorrect(s0_uncorr,inc,slat,slon,inc_correction_model);
    }
    // update checksum
    checksum_[BEAMMASK] += beam_mask;

    // All lats and lon in radians until output in degrees by multiplying
    // by rtd
    float slat_in_degrees=slat*radtodeg;
    float slon_in_degrees=slon*radtodeg;
    float inc_in_degrees=inc*radtodeg;
    if(numlooks<num_looks_required || !output_enable_[INC]) inc_in_degrees=BAD_VALUE;
    
    // write pixel to appropriate io buffers
    beammask_io_bfr[j]=beam_mask;
    s0_corr_io_bfr[j]=s0_corr;
    lat_io_bfr[j]=slat_in_degrees; 
    lon_io_bfr[j]=positiveWestLon(slon_in_degrees);
    inc_io_bfr[j]=inc_in_degrees;
    range_io_bfr[j]=range;
    dop_io_bfr[j]=doppler;
    s0_uncorr_io_bfr[j]=s0_uncorr;
    start_burst_num_io_bfr[j]=burst_num_1;
    nom_burst_num_io_bfr[j]=burst_num_1;
    end_burst_num_io_bfr[j]=burst_num_2;
    if(numlooks>=255){
      num_looks_io_bfr[j]=(unsigned char)255;
    }
    else{
      num_looks_io_bfr[j]=(unsigned char)numlooks;
    }
    X_io_bfr[j]=X;
    // HACK to output s0_nsub_uncorr to X file until we add the new file
    if(spp.noise_subtraction_on) X_io_bfr[j]=s0ns;
    s0ne_io_bfr[j]=s0ne;

    float s0std;
    if(s0ne==BAD_VALUE || s0ns==BAD_VALUE) s0std=BAD_VALUE;
    else s0std=sqrt((fabs(s0ns)*fabs(s0ns)+2*fabs(s0ns)*s0ne+s0ne*s0ne)/numlooks);

    std_io_bfr[j]=s0std;
    // Update values to be written to PDS label
    if (minimum_latitude_ > slat_in_degrees)
      minimum_latitude_ = slat_in_degrees;
    if (maximum_latitude_ < slat_in_degrees)
      maximum_latitude_ = slat_in_degrees;
    // Shift slon_in_degrees to [-180, 180) to avoid discontinuity at 0 degrees
    if (dbg.level==2)
    {
      dbg.file << "BIDR::outputLine:  ";
      dbg.file << "min, max lon and slon_in_degrees before shift:" <<
        minimum_longitude_ << " " << maximum_longitude_ << " " <<
	slon_in_degrees << endl;
    }
    // set first_longitude for first pixel
    if(first_longitude_==-1000){
      first_longitude_=slon_in_degrees;
      minimum_longitude_=slon_in_degrees+180.0;
      maximum_longitude_=slon_in_degrees-180.0;      
    }
    while (slon_in_degrees >= first_longitude_+180.0) { slon_in_degrees -= 360.0; }
    while (slon_in_degrees < first_longitude_-180.0) { slon_in_degrees += 360.0; }
    if (minimum_longitude_ > slon_in_degrees)
      minimum_longitude_ = slon_in_degrees;
    if (maximum_longitude_ < slon_in_degrees)
      maximum_longitude_ = slon_in_degrees;
    if (dbg.level==2)
    {
      dbg.file << "BIDR::outputLine:  ";
      dbg.file << "min, max lon and slon_in_degrees after shift:" <<
        minimum_longitude_ << " " << maximum_longitude_ << " " <<
	slon_in_degrees << endl;
    }
  } // end of latitude loop

  
  // output i/o buffers to the BIDR files  
  switch (procMode){
  case NORMAL:
    if(     (fwrite((void*)beammask_io_bfr,sizeof(char),numlats,beammask_file)!=numlats) 
	    ||  (fwrite((void*)s0_corr_io_bfr,sizeof(float),numlats,s0_corr_file)!=numlats) 
	    ||  (fwrite((void*)lat_io_bfr,sizeof(float),numlats,lat_file)!=numlats) 
	    ||  (fwrite((void*)lon_io_bfr,sizeof(float),numlats,lon_file)!=numlats) 
	    ||  (fwrite((void*)inc_io_bfr,sizeof(float),numlats,inc_file)!=numlats) 
	    ||  (fwrite((void*)s0_uncorr_io_bfr,sizeof(float),numlats,s0_uncorr_file)!=numlats) 
	    ||  (fwrite((void*)start_burst_num_io_bfr,sizeof(int),numlats,start_burst_num_file)!=numlats) 
	    ||  (fwrite((void*)end_burst_num_io_bfr,sizeof(int),numlats,end_burst_num_file)!=numlats) 
	    ||  (fwrite((void*)num_looks_io_bfr,sizeof(unsigned char),numlats,num_looks_file)!=numlats)
	    ||  (fwrite((void*)X_io_bfr,sizeof(float),numlats,X_file)!=numlats) 
	    ||  (fwrite((void*)std_io_bfr,sizeof(float),numlats,std_file)!=numlats) 
	    ||  (fwrite((void*)s0ne_io_bfr,sizeof(float),numlats,s0ne_file)!=numlats) 
	    ){
      
      cerr << "BIDR::outputLine error writing BIDR files." << endl;
      exit(1);
    }
    break;

  case TOPO:
  case TOPO_PLUS:
    updateTopoArrays(first_i);
    if(     (fwrite((void*)beammask_io_bfr,sizeof(char),numlats,beammask_file)!=numlats) 

	    ||  (fwrite((void*)s0_uncorr_io_bfr,sizeof(float),numlats,s0_uncorr_file)!=numlats) 
	    ||  (fwrite((void*)num_looks_io_bfr,sizeof(unsigned char),numlats,num_looks_file)!=numlats)
	    ||  (fwrite((void*)X_io_bfr,sizeof(float),numlats,X_file)!=numlats) 
	    ){
      
      cerr << "BIDR::outputLine error writing BIDR files." << endl;
      exit(1);
    }
    break;
  case STEREO:
    if(     (fwrite((void*)beammask_io_bfr,sizeof(char),numlats,beammask_file)!=numlats) 

	    ||  (fwrite((void*)range_io_bfr,sizeof(float),numlats,range_file)!=numlats) 
	    ||  (fwrite((void*)dop_io_bfr,sizeof(float),numlats,dop_file)!=numlats) 
	    ||  (fwrite((void*)s0_uncorr_io_bfr,sizeof(float),numlats,s0_uncorr_file)!=numlats) 
	    ||  (fwrite((void*)nom_burst_num_io_bfr,sizeof(int),numlats,nom_burst_num_file)!=numlats) 
	    ||  (fwrite((void*)X_io_bfr,sizeof(float),numlats,X_file)!=numlats) 
	    ){
      
      cerr << "BIDR::outputLine error writing BIDR files." << endl;
      exit(1);
    }
    break;
  default:
    cerr << "Bad BIDR::procMode in BIDR::outputLine" << endl;
    exit(1);
  } // end of procMode switch
  } // end of forced region conditional
  // advance Grid valid region by one longitude index
  grid.advance();
  num_lines_of_output++;
  // reinitialize the line which was just written
  initLongitudeLine(first_i);
  // calculate a new line of standard lat/lon values
  int numlons=grid.numLons();
  for(unsigned int j=0;j<numlats;j++){
    double lon_in_rad=grid.lineNumberToLonInRad(num_lines_of_output+numlons-1);
    double lat_in_rad=grid.latInRad(j);
    slat_bfr[first_i][j]=proj.standardLatInRad(lon_in_rad,lat_in_rad);
    slon_bfr[first_i][j]=proj.standardLonInRad(lon_in_rad,lat_in_rad);
  }
  if (dbg.level)
  {
    dbg.file << "BIDR::outputLine:" << endl;
    dbg.file << "Maximum latitude (deg):      " << maximum_latitude_  << endl;
    dbg.file << "Minimum latitude (deg):      " << minimum_latitude_  << endl;
    dbg.file << "PE Maximum longitude (deg):  " << maximum_longitude_ << endl;
    dbg.file << "PE Minimum longitude (deg):  " << minimum_longitude_ << endl;
  }
}

void BIDR::outputBeamDeltas(const SARProcParams& spp){
  if(procMode==STEREO || procMode==TOPO || procMode==TOPO_PLUS){
    cerr << "BIDR::outputBeamDeltas fails in STEREO or TOPO mode" << endl;
    exit(1);
  }
  int i=grid.firstValidLongitudeIndex();
  int rowno=num_lines_of_output+1;

  // compute SABs
  int sab[5];
  int n=grid.numLats();
  if(spp.single_beam_mode){
    ErrorMessage e("Cannot Output Beam Deltas in single_beam_mode");
    e.throwMe();
  }
  for(int b=0;b<5;b++){
    int max_sab=-1;
    int min_sab=1000000;
    for(int j=0;j<n;j++){
      if(start_burst_num[i][b][j]>=0 && start_burst_num[i][b][j]<min_sab)
	min_sab=start_burst_num[i][b][j];
      if(end_burst_num[i][b][j]>max_sab)
	max_sab=end_burst_num[i][b][j];
       
    }
    int k=(max_sab-min_sab)/5;
    k=k/2;
    sab[b]=min_sab+k*5;
    if(sab[b]>200000) return;
  }
  // compute range_offsets (needs spp.set_burst_range_dop = true )
  float roff[5];
  if(!spp.set_burst_range_dop){
    ErrorMessage e("Cannot Output Beam Deltas unless spp.set_burst_range_dop is true");
    e.throwMe();
  }
  for(int b=0;b<5;b++) roff[b]=spp.range_offset[sab[b]];
  // compute deltas and js
  float dbeam[4],dbeam_top[4],dbeam_bot[4],inc[4];
  float joverlap[4];
  int jno[4];
  for(int k=0;k<4;k++){
    float dval1=0,dval2=0;
    int dnum=0;
    float jval=0;
    int jnum=0;
    float incval=0;
    for(int j=0;j<n;j++){
      if(num_looks[i][k][j]>=num_looks_required &&
	 num_looks[i][k+1][j]>=num_looks_required){
	jnum++;
        jval+=j;
        dnum++;
        float s01=s0_uncorr_bfr[i][k][j]/num_looks[i][k][j];
        float s02=s0_uncorr_bfr[i][k+1][j]/num_looks[i][k+1][j];
        dval1+=s01;
        dval2+=s02;
        incval+=(inc_bfr[i][k][j]/num_looks[i][k][j])*radtodeg;
#define DEBUGBD
#ifdef DEBUGBD
int n1=num_looks[i][k][j];
int n2=num_looks[i][k+1][j];
float inc1=(inc_bfr[i][k][j]/num_looks[i][k][j])*radtodeg;
float inc2=(inc_bfr[i][k+1][j]/num_looks[i][k+1][j])*radtodeg;
 float x1=X_bfr[i][k][j]/n1;
 float x2=X_bfr[i][k+1][j]/n2;
 int iinc1=(int)(floor(inc1*1000+0.5));
 int iinc2=(int)(floor(inc2*1000+0.5));
 int is01=(int)(floor(s01*1000+0.5));
 int is02=(int)(floor(s02*1000+0.5));
 int ix1=(int)(floor(x1+0.5));
 int ix2=(int)(floor(x2+0.5));
if(fwrite(&rowno,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&j,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&k,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&is01,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&is02,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&ix1,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&ix2,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&iinc1,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&iinc2,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&n1,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&n2,sizeof(int),1,beam_deltas_file)!=1)
  {
    ErrorMessage e("Cannot write to beam_deltas_file");
    e.throwMe();
  }
#endif        
      }
    }
    joverlap[k]=jval/jnum;
    dbeam[k]=dval1/dval2;
    jno[k]=jnum;
    dbeam_top[k]=dval1;
    dbeam_bot[k]=dval2;
    inc[k]=incval/jnum;
  }
  // write to file
  int iroff[5],ij[4],idbeam[4],itop[4],ibot[4],iinc[4];
  for(int k=0;k<5;k++) iroff[k]=int(floor(roff[k]*1000 + 0.5));
  for(int k=0;k<4;k++){
    ij[k]=int(floor(joverlap[k]*1000 + 0.5));
    idbeam[k]=int(floor(dbeam[k]*1000 + 0.5));
    iinc[k]=int(floor(inc[k]*1000 + 0.5));
    itop[k]=int(floor(dbeam_top[k]*10000 + 0.5));
    ibot[k]=int(floor(dbeam_bot[k]*10000 + 0.5));
  }
#ifndef DEBUGBD
  if(fwrite(&rowno,sizeof(int),1,beam_deltas_file)!=1
     ||fwrite(&(sab[0]),sizeof(int),5,beam_deltas_file)!=5
     || fwrite(&(iroff[0]),sizeof(int),5,beam_deltas_file)!=5
     || fwrite(&(ij[0]),sizeof(int),4,beam_deltas_file)!=4
       || fwrite(&(idbeam[0]),sizeof(int),4,beam_deltas_file)!=4
       || fwrite(&(itop[0]),sizeof(int),4,beam_deltas_file)!=4
       || fwrite(&(ibot[0]),sizeof(int),4,beam_deltas_file)!=4
       || fwrite(&(iinc[0]),sizeof(int),4,beam_deltas_file)!=4
       || fwrite(&(jno[0]),sizeof(int),4,beam_deltas_file)!=4){
    ErrorMessage e("Cannot write to beam_deltas_file");
    e.throwMe();
  }
#endif
  printf("ROW %d SABS %d %d %d %d %d ROFF %g %g %g %g %g\n",
	 rowno,sab[0],sab[1],sab[2],sab[3],sab[4],roff[0],roff[1],
	 roff[2],roff[3],roff[4]);
  printf("JS %g %g %g %g DBEAMS %g %g %g %g INC %g %g %g %g \n",
	 joverlap[0],joverlap[1],joverlap[2],joverlap[3],
	 dbeam[0],dbeam[1],dbeam[2],dbeam[3],inc[0],inc[1],inc[2],inc[3]);
}
void BIDR::initLongitudeLine(int i)
{
  int num_lats=grid.numLats();
  for(int b=0;b<num_beams;b++){
    for(int j=0;j<num_lats;j++){
      if(output_enable_[START_BURST_NUM])start_burst_num[i][b][j]=-1;
      if(output_enable_[END_BURST_NUM])end_burst_num[i][b][j]=-1;
      num_looks[i][b][j]=0;
      if(output_enable_[INC])inc_bfr[i][b][j]=0;
      if(procMode==STEREO) max_X_bfr[i][b][j]=0;
      X_bfr[i][b][j]=0;
      s0_uncorr_bfr[i][b][j]=0;
      s0_nsub_uncorr_bfr[i][b][j]=0;
      s0ne_bfr[i][b][j]=0;
      if(procMode==TOPO || procMode==TOPO_PLUS){
	noise_equiv_s0[i][b][j]=0;
	resfactor[i][b][j]=0;
	Xambig_bfr[i][b][j]=0;
      }
    }      
  }
  bursts_since_update[i]=0;
  for(int b=0;b<num_beams;b++){
    beams_found[i][b]=false;
  }
  for(int j=0;j<num_lats;j++){
    slat_bfr[i][j]=0;
    slon_bfr[i][j]=0;
  }
}

BIDR::~BIDR()
{
  int num_lons=grid.numLons();
  int num_lats=grid.numLats();
  free(area_pixel);
  free_array((void*)s0_uncorr_bfr,3,num_lons,num_beams,num_lats);
  free_array((void*)s0_nsub_uncorr_bfr,3,num_lons,num_beams,num_lats);
  free_array((void*)s0ne_bfr,3,num_lons,num_beams,num_lats);
  free_array((void*)X_bfr,3,num_lons,num_beams,num_lats);
  free_array((void*)bursts_since_update,1,num_lons);
  free_array((void*)beams_found,2,num_lons,num_beams);
  free_array((void*)num_looks,3,num_lons,num_beams,num_lats);
  free_array((void*)slat_bfr,2,num_lons,num_lats);
  free_array((void*)slon_bfr,2,num_lons,num_lats);
  free_array((void*)checksum_,1,BIDRTYPECOUNT);
  free_array((void*)beammask_io_bfr,1,num_lats);
  free_array((void*)s0_corr_io_bfr,1,num_lats);
  free_array((void*)lat_io_bfr,1,num_lats);
  free_array((void*)lon_io_bfr,1,num_lats);
  free_array((void*)inc_io_bfr,1,num_lats);
  free_array((void*)range_io_bfr,1,num_lats);
  free_array((void*)dop_io_bfr,1,num_lats);
  free_array((void*)s0_uncorr_io_bfr,1,num_lats);
  free_array((void*)start_burst_num_io_bfr,1,num_lats);
  free_array((void*)end_burst_num_io_bfr,1,num_lats);
  free_array((void*)nom_burst_num_io_bfr,1,num_lats);
  free_array((void*)num_looks_io_bfr,1,num_lats);

  if(procMode==NORMAL){
    free_array((void*)inc_bfr,3,num_lons,num_beams,num_lats);
    free_array((void*)start_burst_num,3,num_lons,num_beams,num_lats);
    free_array((void*)end_burst_num,3,num_lons,num_beams,num_lats);
  }
  if(procMode==STEREO){
    free_array((void*)range_bfr,3,num_lons,num_beams,num_lats);
    free_array((void*)dop_bfr,3,num_lons,num_beams,num_lats);
    free_array((void*)max_X_bfr,3,num_lons,num_beams,num_lats);
    free_array((void*)nom_burst_num,3,num_lons,num_beams,num_lats);    
  }
  if(procMode==TOPO || procMode==TOPO_PLUS){
    free_array((void*)inc_bfr,3,num_lons,num_beams,num_lats);
    free_array(dsigma0,2,4,MAX_BIDR_LINES);
    free_array(topo_inc1,2,4,MAX_BIDR_LINES);
    free_array(topo_inc2,2,4,MAX_BIDR_LINES);
    free_array(topo_lat,2,4,MAX_BIDR_LINES);
    free_array(topo_lon,2,4,MAX_BIDR_LINES);
    free_array(topo_wnl,2,4,MAX_BIDR_LINES);
    free_array(mid_ds0,2,4,MAX_BIDR_LINES);
    free_array(slope_ds0,2,4,MAX_BIDR_LINES);
    free_array(rms_ds0,2,4,MAX_BIDR_LINES);
    free_array(power1,2,4,MAX_BIDR_LINES);
    free_array(power2,2,4,MAX_BIDR_LINES);
    free_array(noisefloor1,2,4,MAX_BIDR_LINES);
    free_array(noisefloor2,2,4,MAX_BIDR_LINES);
    free_array(ambrat1,2,4,MAX_BIDR_LINES);
    free_array(ambrat2,2,4,MAX_BIDR_LINES);
    free_array(dsigma0std,2,4,MAX_BIDR_LINES);
    free_array(dsigma0bias,2,4,MAX_BIDR_LINES);
    free_array(sumsigma0,2,4,MAX_BIDR_LINES);
    free_array(overlap_col,2,4,MAX_BIDR_LINES);
    free_array(overlap_width,2,4,MAX_BIDR_LINES);
    free_array((void*)resfactor,3,num_lons,num_beams,num_lats);
    free_array((void*)noise_equiv_s0,3,num_lons,num_beams,num_lats);
    free_array((void*)Xambig_bfr,3,num_lons,num_beams,num_lats);
  }
  if(procMode==TOPO_PLUS){
    free_array(topoplus_s01,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s02,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0ne1,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0ne2,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0baqscale1,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0baqscale2,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0amb1,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
    free_array(topoplus_s0amb2,3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH);
  }
  closeFiles();
}

void BIDR::closeFiles(){
  if(beam_deltas_file)fclose(beam_deltas_file);
  if(beammask_file) fclose(beammask_file);
  if(X_file) fclose(X_file);
  if(std_file) fclose(std_file);
  if(s0ne_file) fclose(s0ne_file);
  if(s0_corr_file) fclose(s0_corr_file);
  if(lat_file) fclose(lat_file);
  if(lon_file) fclose(lon_file);
  if(range_file) fclose(range_file);
  if(dop_file) fclose(dop_file);
  if(inc_file) fclose(inc_file);
  if(s0_uncorr_file) fclose(s0_uncorr_file);
  if(start_burst_num_file) fclose(start_burst_num_file);
  if(end_burst_num_file) fclose(end_burst_num_file);
  if(nom_burst_num_file) fclose(nom_burst_num_file);
  if(num_looks_file) fclose(num_looks_file);
  if(topoplus_file) fclose(topoplus_file);
  
}

bool BIDR::longitudeLineComplete(int i, const SARProcParams& spp){
  bool retval=false;
  DebugInfo dbg("BIDR::longitudeLineComplete");
  if(dbg.level){
    dbg.file << "BIDR::LongitudeLineComplete called for index="<<i << endl;
    dbg.file << " burst_since_update="<< bursts_since_update[i] 
	     << " beams_found=[";
    for(int c=0;c<num_beams;c++) dbg.file << " " << beams_found[i][c];
    dbg.file << " ]" << endl;
  }
  // need to not have been updated for last 10 bursts to be complete
  if(bursts_since_update[i]>10){
    retval=true;
  }
  return(retval);
}


void  BIDR::multipleBeamMerge(const SARProcParams& spp, int i, int j, 
			  unsigned char& beam_mask, float& s0_uncorr, 
			float& inc, int& burstno1, int& burstno2, 
			      int& numlooks, float& X, float& s0ns, float& s0ne,
			      float& doppler, float& range){
  int b=0; // selected beam index 0-4

  // Single beam mode := only 1 beam available select that beam
  if (spp.single_beam_mode){
    b=0; // allows less memory usage for single beam mode
  }


  // Pick the best beam, does not do a weighted sum
  // best beam definition is currently
  // must looks, and if looks are equal highest Xfactor
  // this should handle BAD_VALUES correctly since it just passes the values
  // for the selected beam
  else if(!feather_beams){
    b=-1;
    int most_looks=0;
    float bestX=-1e99;
    if(Xcheckdir==-1) bestX=1e99;


    if((procMode==TOPO || procMode==TOPO_PLUS) && overlap_set){
      int k=num_lines_of_output;
      bool beam1_left=(overlap_col[0][k]<overlap_col[1][k]);
      if(beam1_left==true){
	if(j<overlap_col[0][k]-overlap_width[0][k]) b=0;
	else if(j<overlap_col[0][k]) b=1;
	else if(j<=overlap_col[0][k]+overlap_width[0][k]) b=0;
	

	else if(j<overlap_col[1][k]-overlap_width[1][k]) b=1;
	else if(j<overlap_col[1][k]) b=2;
	else if(j<=overlap_col[1][k]+overlap_width[1][k]) b=1;

	else if(j<overlap_col[2][k]-overlap_width[2][k]) b=2;
	else if(j<overlap_col[2][k]) b=3;
	else if(j<=overlap_col[2][k]+overlap_width[2][k]) b=2;
	
	else if(j<overlap_col[3][k]-overlap_width[3][k] ) b=3;
	else if(j<overlap_col[3][k]) b=4;
	else if(j<=overlap_col[3][k]+overlap_width[3][k] ) b=3;
	else b=4;
      }
      else{
	if(j>overlap_col[0][k]+overlap_width[0][k]) b=0;
	else if(j>overlap_col[0][k]) b=1;
	else if(j>=overlap_col[0][k]-overlap_width[0][k]) b=0;

	else if(j>overlap_col[1][k]+overlap_width[1][k]) b=1;
	else if(j>overlap_col[1][k]) b=2;
	else if(j>=overlap_col[1][k]-overlap_width[1][k]) b=1;

	else if(j>overlap_col[2][k]+overlap_width[2][k]) b=2;
	else if(j>overlap_col[2][k]) b=3;
	else if(j>=overlap_col[2][k]-overlap_width[2][k]) b=2;
	
	else if(j>overlap_col[3][k]+overlap_width[3][k]) b=3;
	else if(j>overlap_col[3][k]) b=4;
	else if(j>=overlap_col[3][k]-overlap_width[3][k]) b=3;
	else b=4;
      }
      if(num_looks[i][b][j]==0) b=-1;
    }
    if(dominant_beam3 && num_looks[i][2][j]>num_looks_required){
      b=2;
    }
 
    // nominal beam selection is used if neither special cases applies
    // or TOPO beam selection fails due to nonoverlapping case.
    if(b==-1){
      b=0;
      for(int t=0;t<num_beams;t++){
	if(num_looks[i][t][j]>most_looks){
	  most_looks=num_looks[i][t][j];
	  b=t;
	  bestX=X_bfr[i][t][j];
	}
	else if(num_looks[i][t][j]==most_looks){
	  if(Xcheckdir>0 && X_bfr[i][t][j]>bestX){
	    b=t;
	    bestX=X_bfr[i][t][j];
	  }
	  else if(Xcheckdir<0 && X_bfr[i][t][j]<bestX){
	    b=t;
	    bestX=X_bfr[i][t][j];
	  }
	}
      }
    }
  }
  
  if(!feather_beams){ 
    // Compute output quantities from selected beam
    if(spp.single_beam_mode)
      beam_mask=0x01 << spp.beam_number-1; 
    else beam_mask=0x01 << b;
    
    s0_uncorr=s0_uncorr_bfr[i][b][j];
    if(procMode==NORMAL){
      inc=inc_bfr[i][b][j];
      burstno1=start_burst_num[i][b][j];
      burstno2=end_burst_num[i][b][j];
    }

    if(procMode==TOPO || procMode==TOPO_PLUS){
      inc=inc_bfr[i][b][j];
    }

    if(procMode==STEREO){
      range=range_bfr[i][b][j];
      doppler=dop_bfr[i][b][j];
      burstno1=nom_burst_num[i][b][j];
    }
    numlooks=num_looks[i][b][j];  
    X=X_bfr[i][b][j];
    s0ns=s0_nsub_uncorr_bfr[i][b][j];
    s0ne=s0ne_bfr[i][b][j];
    // if numlooks is zero set beam mask to zero as well 
    // (all other bad values should automatically be handled
    // Notice that for nonzero numlooks some beam is selected even if 
    // there are not enough looks to produce a valid s0
    if(numlooks==0) beam_mask=0;
  
  }
  // beam feathering code
  if(feather_beams){
    if(procMode==STEREO || procMode==TOPO || procMode==TOPO_PLUS){
      cerr << "BIDR::multipleBeamMerge Beam Feathering incompatible with STEREO or TOPO mode" << endl;
      exit(1);
    }
    // loop through beams
    bool isvalid=false; // boolean to check if any beam is valid
    float sumw=0.0; // sum of weights normailiztion constant
    float w;
    int bmsk;

    // intialize pixel value quantities
    beam_mask=0;
    s0_uncorr=0;
    inc=0;
    int burstno1=100000000;
    int burstno2=-1;
    X=0;
    s0ns=0;
    s0ne=0;
    numlooks=0;
    for(int b=0;b<5;b++){
      if(num_looks[i][b][j]>=num_looks_required){
	isvalid=true;
	w=1.0/sqrt((float)num_looks[i][b][j]); // per beam weight
        sumw+=w;  // keep track of sum of weight
        bmsk=0x01 << b; // compute bit for beam mask

        // accumulate backplane and image pixel values over beams
        beam_mask= beam_mask | bmsk;
        s0_uncorr+=s0_uncorr_bfr[i][b][j]*w; 
	inc+=inc_bfr[i][b][j]*w; // for completeness (could just take first
                                 // valid value)
        if(start_burst_num[i][b][j] < burstno1) 
	  burstno1=start_burst_num[i][b][j];
        if(end_burst_num[i][b][j] > burstno2) 
	  burstno2=end_burst_num[i][b][j];
        numlooks+=num_looks[i][b][j];
	X+=X_bfr[i][b][j]*w;  
	s0ns+=s0_nsub_uncorr_bfr[i][b][j]*w;  
	s0ne+=s0ne_bfr[i][b][j]*w;  

      }
    }
    
   
    // handle no looks case
    if(!isvalid){
      beam_mask=0;
      s0_uncorr=BAD_VALUE;
      s0ns=BAD_VALUE;
      s0ne=BAD_VALUE;
      inc=BAD_VALUE;
      X=BAD_VALUE;
      numlooks=0; // looks for beams with less than the required number of
                   // looks are not counted when feathering is applied
      burstno1=-1;
      burstno2=-1;
    }

    else{
      // normalize pixel values
      inc=inc/sumw;
      s0_uncorr=s0_uncorr/sumw;
      s0ns=s0ns/sumw;
      s0ne=s0ne/sumw;
      X=X/sumw;
    }
  } // end beam feathering case
}

void  
BIDR::updateTopoPlusArrays(const SARProcParams& spp, int i, int j, int pn, int abs_i, int col){
  topoplus_s01[pn][abs_i][col]=s0_nsub_uncorr_bfr[i][pn][j];
  topoplus_s02[pn][abs_i][col]=s0_nsub_uncorr_bfr[i][pn+1][j];
  topoplus_s0ne1[pn][abs_i][col]=noise_equiv_s0[i][pn][j];
  topoplus_s0ne2[pn][abs_i][col]=noise_equiv_s0[i][pn+1][j];
  topoplus_s0baqscale1[pn][abs_i][col]=s0_nsub_uncorr_bfr[i][pn][j]/(s0_uncorr_bfr[i][pn][j]-noise_equiv_s0[i][pn][j]);
  topoplus_s0baqscale2[pn][abs_i][col]=s0_nsub_uncorr_bfr[i][pn+1][j]/(s0_uncorr_bfr[i][pn+1][j]-noise_equiv_s0[i][pn+1][j]);

  // need to fix these
  topoplus_s0amb1[pn][abs_i][col]=BAD_VALUE;
  topoplus_s0amb2[pn][abs_i][col]=BAD_VALUE;
  
}

float 
incAngleCorrect(float s0_uncorr, float inc, float slat,
		      float slon, corrTypeE inc_correction_model){
  float retval, s,c,s2,c4,f1,f2,f3;
  static int val=BAD_VALUE_HEX;
  static float* ptr=(float*) &val;
  static float BAD_VALUE=*ptr;
  if(s0_uncorr==BAD_VALUE || inc == BAD_VALUE){
    retval=BAD_VALUE;
  }
  else if(inc_correction_model==RKIRK01TA){
    retval=s0_uncorr*sqrt(2)*sin(inc);
  }
  else if(inc_correction_model==ENCELADUS_WYE1){
    retval=s0_uncorr*2.9165/(3.71*pow(cos(inc),1.46));
  }
  else if(inc_correction_model==RHEA_WYE1){
    retval=s0_uncorr*1.6930/(2.15*pow(cos(inc),1.45));
  }
  else if(inc_correction_model==DIONE_WYE1){
    retval=s0_uncorr*1.4253/(1.825*pow(cos(inc),1.5));
  }
  else if(inc_correction_model==HAGHAG){
    /** old version
    s=sin(inc);
    s2=s*s;
    c=cos(inc);
    c4=c*c*c*c;
    f1=(c4+771.3241*s2);
    f1=3.3167/sqrt(f1*f1*f1);
    f2=(c4+18.6535*s2);
    f2=0.3749/sqrt(f2*f2*f2);
    f3=0.3456*pow(c,1.7531);
    retval=s0_uncorr*0.2013/(f1+f2+f3);
    ***/

    s=sin(inc);
    s2=s*s;
    c=cos(inc);
    c4=c*c*c*c;
    f1=(c4+893.9677*s2);
    f1=2.8126/sqrt(f1*f1*f1);
    f2=(c4+34.1366*s2);
    f2=0.5824/sqrt(f2*f2*f2);
    f3=0.3767*pow(c,1.9782);
    retval=s0_uncorr*0.2907/(f1+f2+f3);

  }
  /**
  else if(inc_correction_model==DBLINE041104){
    retval=s0_uncorr*dbline_s045/
      (pow(10,0.1*(dbline_slope*inc+dbline_offset)));
  }
  else if(inc_correction_model==VENUS){
    retval=s0_uncorr*venus_s045/muhleman_backscatter(venus_k1,venus_k2,inc);
  }
  **/
  else retval=BAD_VALUE;
  return(retval);
}



float 
theorBackscatter(float inc, corrTypeE inc_correction_model){
  float retval, s,c,s2,c4,f1,f2,f3;
  static int val=BAD_VALUE_HEX;
  static float* ptr=(float*) &val;
  static float BAD_VALUE=*ptr;

  // avoid numerical problems for inc angles near 90
  if(inc>= 0.99999*pi/2) return(0.0);

  if(inc == BAD_VALUE){
    retval=BAD_VALUE;
  }
  else if(inc_correction_model==RKIRK01TA){
    retval=1/sin(inc);
  }
  else if(inc_correction_model==ENCELADUS_WYE1){
    retval=(3.71*pow(cos(inc),1.46));
  }
  else if(inc_correction_model==RHEA_WYE1){
    retval=(2.15*pow(cos(inc),1.45));
  }
  else if(inc_correction_model==DIONE_WYE1){
    retval=(1.825*pow(cos(inc),1.5));
  }
  else if(inc_correction_model==HAGHAG){
    s=sin(inc);
    s2=s*s;
    c=cos(inc);
    c4=c*c*c*c;
    f1=(c4+893.9677*s2);
    f1=2.8126/sqrt(f1*f1*f1);
    f2=(c4+34.1366*s2);
    f2=0.5824/sqrt(f2*f2*f2);
    f3=0.3767*pow(c,1.9782);
    retval=f1+f2+f3;

  }
  else retval=BAD_VALUE;
  return(retval);
}

void 
BIDR::flush(const SARProcParams& spp){

  static DebugInfo dbg("BIDR::flush");

  // output remaining valid longitude lines and extra lines so that longitude boundary is in integer degrees.
  while(num_lines_of_output<grid.totalNumLons()){
    outputLine(spp);
  }
  if(procMode==TOPO_PLUS && overlap_set){
    int n1=4, n2=MAX_BIDR_LINES, n3=MAX_TOPO_OVERLAP_WIDTH;
    if(!fwrite(&n1,1,sizeof(int),topoplus_file) ||
       !fwrite(&n2,1,sizeof(int),topoplus_file) || 
       !fwrite(&n3,1,sizeof(int),topoplus_file) ||
       !write_array(topoplus_file,topoplus_s01,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s02,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0ne1,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0ne2,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0baqscale1,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0baqscale2,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0amb1,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH) ||
       !write_array(topoplus_file,topoplus_s0amb2,sizeof(float),3,4,MAX_BIDR_LINES,MAX_TOPO_OVERLAP_WIDTH)){
      cerr << "Error writing TOPOPLUS file\n" << endl;
      exit(1);
    }
    // debugging tools
    //cout << "I MADE IT" << endl;
    //fclose(topoplus_file);
    //exit(1);
  }
  if (dbg.level)
  {
    dbg.file << "BIDR::flush:" << endl;
    dbg.file << "Maximum latitude (deg):      " << maximum_latitude_  << endl;
    dbg.file << "Minimum latitude (deg):      " << minimum_latitude_  << endl;
    dbg.file << "PE Maximum longitude (deg):  " << maximum_longitude_ << endl;
    dbg.file << "PE Minimum longitude (deg):  " << minimum_longitude_ << endl;
  }
  
}

void
BIDR::incrementBurstsSinceUpdate(){
  int numlons=grid.numLons();;
  for(int i=0;i<=numlons;i++) bursts_since_update[i]++;
}

//=====================================================================
// writeHeaders()
//
// Write PDS labels to each output BIDR product.
//=====================================================================
void BIDR::writeHeaders()
{
  writeHeader(inc_file, INC);
  writeHeader(beammask_file, BEAMMASK);
  writeHeader(s0_uncorr_file, S0_UNCORR);
  writeHeader(X_file, S0NSQT_UNCORR);
  writeHeader(std_file, S0_STD);
  writeHeader(s0ne_file, S0_NOISE_EQUIV);
  writeHeader(s0_corr_file, S0_CORR);
  writeHeader(lat_file, LAT);
  writeHeader(lon_file, LON);
  writeHeader(dop_file, DOPPLER);
  writeHeader(range_file, RANGE);
  writeHeader(start_burst_num_file, START_BURST_NUM);
  writeHeader(end_burst_num_file, END_BURST_NUM);
  writeHeader(nom_burst_num_file, NOMINAL_BURST_NUM);
  writeHeader(num_looks_file, NUM_LOOKS);
}

//=====================================================================
// rewriteHeaders()
//
// Rewrite the PDS labels in each output BIDR product with updated
// values.
//=====================================================================
void BIDR::rewriteHeaders()
{
  rewriteHeader(inc_file, INC);
  rewriteHeader(beammask_file, BEAMMASK);
  rewriteHeader(s0_uncorr_file, S0_UNCORR);
  rewriteHeader(X_file, S0NSQT_UNCORR);
  rewriteHeader(std_file, S0_STD);
  rewriteHeader(s0ne_file, S0_NOISE_EQUIV);
  rewriteHeader(s0_corr_file, S0_CORR);
  rewriteHeader(lat_file, LAT);
  rewriteHeader(lon_file, LON);
  rewriteHeader(range_file, RANGE);
  rewriteHeader(dop_file, DOPPLER);
  rewriteHeader(start_burst_num_file, START_BURST_NUM);
  rewriteHeader(end_burst_num_file, END_BURST_NUM);
  rewriteHeader(nom_burst_num_file, NOMINAL_BURST_NUM);
  rewriteHeader(num_looks_file, NUM_LOOKS);
}

//=====================================================================
// writeHeader()
//
// Write a PDS label to a BIDR output file.  The file pointer is
// assumed to point to the beginning of the file.
//=====================================================================
void BIDR::writeHeader(FILE* f, BIDRTypeE bidr_type)
{
  if(!output_enable_[bidr_type]) return;
  string bidr_label = constructLabel(bidr_type);
  if (bidr_label == EMPTY_STRING)
  {
    ErrorMessage e("BIDR::writeHeader:  PDS label creation failed");
    e.throwMe();
  }
  if(fwrite((void*)bidr_label.c_str(),sizeof(char),bidr_label.length(),f)!=bidr_label.length()){
    ErrorMessage e("Error rewriting BIDR label");
    e.throwMe();
  }
}

//=====================================================================
// rewriteHeader()
//
// Update the PDS label in a BIDR output file.  This method is intended
// for use when some of the values in the PDS label need to be updated
// and after some (if not all) of the data has been written to the output
// file after the original label.  rewriteHeader() does not update the
// label on a per-value basis; instead, it regenerates the label and
// overwrites all of the existing label.  The actual recalculation of
// new values is contained within constructLabel().  rewriteHeader() also
// assumes that constructLabel() will return a label of exactly the same
// length (in bytes) as the original label.
//=====================================================================
void BIDR::rewriteHeader(FILE* f, BIDRTypeE bidr_type)
{
  if(!output_enable_[bidr_type]) return;
  // Save the current position in the output file and reset the file
  // pointer to the beginning of the file.  After the new label has
  // been written, restore the file pointer to its original position.
  int pos = ftell(f);
  rewind(f);
  string bidr_label = constructLabel(bidr_type);
  if (bidr_label == EMPTY_STRING)
  {
    ErrorMessage e("BIDR::rewriteHeader:  PDS label creation failed");
    e.throwMe();
  }

  if(fwrite((void*)bidr_label.c_str(),sizeof(char),bidr_label.length(),f)!=bidr_label.length()){
    ErrorMessage e("Error rewriting BIDR label");
    e.throwMe();
  }

  if(fseek(f,pos,0)!=0){
    ErrorMessage e("BIDR::rewriteHeader:  unable to seek back to original file pointer");
    e.throwMe();    
  }
}

//=====================================================================
// constructLabel()
//
// This method defines the actual keywords and values that appear in
// the PDS labels of the BIDR products.  The keywords and values are
// specific to the BIDR products.  The string returned by this method
// contains a properly formatted PDS label for the specified type of
// BIDR data and is to be prepended to the output data.
//=====================================================================
string BIDR::constructLabel(BIDRTypeE bidr_type)
{
  //
  // Construct list of keyword-value pairs for BIDR PDS label
  //
  KeyValSeq list;
  int data_records = grid.totalNumLons();
  int record_length = grid.numLats();
  string product_creation_time = Time::getCurrentUtcDOYTime();
  stringstream buf(ios::out);
  double line_projection_offset;
  double sample_projection_offset;
  double radii[3];
  double rotation_matrix[3][3];

  DebugInfo dbg("BIDR::constructLabel");

  if (check_pds_string_lengths_)
  {
    // Populate the table that holds the limits on the lengths of PDS
    // character values.  (The actual checking is done within the
    // constructor for KeyValPair().)
    PDSLabel::setPDSCharValLimits();
  }

  //
  // Calculate values before beginning to construct the label.
  //

  // A_AXIS_RADIUS, B_AXIS_RADIUS, C_AXIS_RADIUS
  // For the general case, the A, B, and C radii are ordered from longest
  // to shortest
  radii[0] = default_target_radii[PositionVector::X].km();
  radii[1] = default_target_radii[PositionVector::Y].km();
  radii[2] = default_target_radii[PositionVector::Z].km();
  sort(radii, radii+3);

  // LINE_PROJECTION_OFFSET, SAMPLE_PROJECTION_OFFSET
  // start_lon and min_lat are values at the edges of the grid, while
  // the line and sample projection offsets are calculated relative to
  // the center of the first pixel in the grid.  As a result, calculations
  // with start_lon and min_lat must apply a half-pixel offset.
  double pixel_offset = 0.5/pixelsperdegree;
  double start_lon = grid.startLon() * radtodeg;
  double min_lat = grid.minLat() * radtodeg;
  double lon_firstpix = start_lon + pixel_offset;
  double lat_firstpix = min_lat + pixel_offset;

  // Wrap longitude to range [-180, 180) degrees
  while (lon_firstpix < -180.0) lon_firstpix += 360.0;
  while (lon_firstpix >= 180.0) lon_firstpix -= 360.0;

  if (dbg.level)
  {
    dbg.file << "BIDR::constructLabel:" << endl;
    dbg.file << "Start lon of grid (deg):       " << setprecision(9) <<
      std::fixed << start_lon << endl;
    dbg.file << "Minimum lat of grid (deg):     " << setprecision(9) <<
      std::fixed << min_lat << endl;
    dbg.file << "Pixel (1,1) center lat (deg):  " << setprecision(9) <<
      std::fixed << lat_firstpix << endl;
    dbg.file << "Pixel (1,1) center lon (deg):  " << setprecision(9) <<
      std::fixed << lon_firstpix << endl;
  }
  line_projection_offset = -lon_firstpix * pixelsperdegree;
  sample_projection_offset = -lat_firstpix * pixelsperdegree;

  // OBLIQUE_PROJ_X_AXIS_VECTOR, OBLIQUE_PROJ_Y_AXIS_VECTOR,
  // OBLIQUE_PROJ_Z_AXIS_VECTOR
  proj.rotationMatrix(rotation_matrix);

  //
  // Create the list of keyword-value pairs in the order in which they
  // are to be written to the label.
  //

  list.push_back(KeyValPair("PDS_VERSION_ID", pds_version_id_, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("/*       FILE FORMAT AND LENGTH */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("RECORD_TYPE", record_type_, !IS_ODL_STRING));
  list.push_back(KeyValPair("RECORD_BYTES",
    record_length*sample_bits_[bidr_type]/BITS_PER_BYTE));
  //
  // Initialize the next two values to their smallest possible values.
  //
  list.push_back(KeyValPair("FILE_RECORDS", data_records+1));
  list.push_back(KeyValPair("LABEL_RECORDS", 1));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair(
    "/*       POINTERS TO START RECORDS OF OBJECTS IN FILE */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  //
  // Initialize the next value to its smallest possible value.
  //
  list.push_back(KeyValPair("^IMAGE", 2));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("/*       IMAGE DESCRIPTION */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("DATA_SET_ID", data_set_id_, IS_ODL_STRING));
  list.push_back(KeyValPair("DATA_SET_NAME", data_set_name_, IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_INSTITUTION_NAME",
    producer_institution_name_, IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_ID", producer_id_, !IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_FULL_NAME", producer_full_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCT_ID", productID(bidr_type),
    !IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCT_VERSION_ID", product_version_,
    2, NO_UNITS));
  list.push_back(KeyValPair("INSTRUMENT_HOST_NAME", instrument_host_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_HOST_ID", instrument_host_id_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_NAME", instrument_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_ID", instrument_id_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("TARGET_NAME", PDSLabel::toUpper(target_name_),
    !IS_ODL_STRING));
  //
  // Leave fractional seconds in for now
  //
  list.push_back(KeyValPair("START_TIME", sclk_start_.utc("ISOD"),
    !IS_ODL_STRING));
  //
  // Leave fractional seconds in for now
  //
  list.push_back(KeyValPair("STOP_TIME", sclk_stop_.utc("ISOD"),
    !IS_ODL_STRING));
  //
  // Since PDS defines SPACECRAFT_CLOCK_START_COUNT and
  // SPACECRAFT_CLOCK_STOP_COUNT as character-valued, format both values as
  // strings so that their lengths can be limit-checked.
  //
  buf << std::setw(9) << std::setfill('0') <<
    sclk_start_.sclk(cassini_str);
  list.push_back(KeyValPair("SPACECRAFT_CLOCK_START_COUNT", buf.str(),
    !IS_ODL_STRING));
  buf.seekp(0, ios::beg);
  buf << std::setw(9) << std::setfill('0') <<
    sclk_stop_.sclk(cassini_str);
  list.push_back(KeyValPair("SPACECRAFT_CLOCK_STOP_COUNT", buf.str(),
    !IS_ODL_STRING));
  //
  // Leave fractional seconds in for now
  //
  list.push_back(KeyValPair("PRODUCT_CREATION_TIME", product_creation_time,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("SOURCE_PRODUCT_ID", source_product_id_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("MISSION_PHASE_NAME", mission_phase_name_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("MISSION_NAME", mission_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("SOFTWARE_VERSION_ID", software_version_,
    IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair(
    "/*       DESCRIPTION OF OBJECTS CONTAINED IN FILE */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("OBJECT", "IMAGE", !IS_ODL_STRING));
  list.push_back(KeyValPair("LINES", data_records));
  list.push_back(KeyValPair("LINE_SAMPLES", record_length));
  list.push_back(KeyValPair("SAMPLE_TYPE", sample_type_[bidr_type],
    IS_ODL_STRING));
  list.push_back(KeyValPair("SAMPLE_BITS", sample_bits_[bidr_type]));
  list.push_back(KeyValPair("CHECKSUM", checksum(bidr_type), 9, NO_UNITS));
  list.push_back(KeyValPair("SCALING_FACTOR", scaling_factor_, 8));
  list.push_back(KeyValPair("OFFSET", offset_, 8));
  list.push_back(KeyValPair("MISSING_CONSTANT", missing_constant_[bidr_type],
    IS_ODL_STRING));
  list.push_back(KeyValPair("NOTE", note_[bidr_type], IS_ODL_STRING));
  list.push_back(KeyValPair("END_OBJECT", "IMAGE", !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("OBJECT", "IMAGE_MAP_PROJECTION",
    !IS_ODL_STRING));
  list.push_back(KeyValPair("^DATA_SET_MAP_PROJECTION",
    data_set_map_projection_catalog_, IS_ODL_STRING));
  list.push_back(KeyValPair("MAP_PROJECTION_TYPE", map_projection_type_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("A_AXIS_RADIUS", radii[0], 6, "km"));
  list.push_back(KeyValPair("B_AXIS_RADIUS", radii[1], 6, "km"));
  list.push_back(KeyValPair("C_AXIS_RADIUS", radii[2], 6, "km"));
  list.push_back(KeyValPair("FIRST_STANDARD_PARALLEL",
    first_standard_parallel_, IS_ODL_STRING));
  list.push_back(KeyValPair("SECOND_STANDARD_PARALLEL",
    second_standard_parallel_, IS_ODL_STRING));
  list.push_back(KeyValPair("POSITIVE_LONGITUDE_DIRECTION",
    positive_longitude_direction_, !IS_ODL_STRING));
  list.push_back(KeyValPair("CENTER_LATITUDE", proj_center_latitude_,
    6, "deg"));
  list.push_back(KeyValPair("CENTER_LONGITUDE", proj_center_longitude_,
    6, "deg"));
  list.push_back(KeyValPair("REFERENCE_LATITUDE", proj.referenceLatitude(),
    6, "deg"));
  list.push_back(KeyValPair("REFERENCE_LONGITUDE", proj.referenceLongitude(),
    6, "deg"));
  list.push_back(KeyValPair("LINE_FIRST_PIXEL", line_first_pixel_));
  list.push_back(KeyValPair("LINE_LAST_PIXEL", data_records));
  list.push_back(KeyValPair("SAMPLE_FIRST_PIXEL", sample_first_pixel_));
  list.push_back(KeyValPair("SAMPLE_LAST_PIXEL", record_length));
  list.push_back(KeyValPair("MAP_PROJECTION_ROTATION",
    map_projection_rotation_, 1, NO_UNITS));
  list.push_back(KeyValPair("MAP_RESOLUTION", (double) pixelsperdegree, 1,
    "pix/deg"));
  list.push_back(KeyValPair("MAP_SCALE", mapScale(), 8, "km/pix"));
  list.push_back(KeyValPair("MAXIMUM_LATITUDE", maximum_latitude_,
    8, 12, "deg"));
  list.push_back(KeyValPair("MINIMUM_LATITUDE", minimum_latitude_,
    8, 12, "deg"));
  // Shift longitudes to [0, 360) before positive-west coordinate conversion.
  // Define extra variables for the shifted longitudes so that we don't
  // change the initial value of maximum_longitude_ before BIDR::outputLine
  // is called for the first time.
  float easternmost_longitude = maximum_longitude_;
  float westernmost_longitude = minimum_longitude_;
  while (easternmost_longitude >= 360.0) { easternmost_longitude -= 360.0; }
  while (easternmost_longitude <    0.0) { easternmost_longitude += 360.0; }
  while (westernmost_longitude >= 360.0) { westernmost_longitude -= 360.0; }
  while (westernmost_longitude <    0.0) { westernmost_longitude += 360.0; }
  list.push_back(KeyValPair("EASTERNMOST_LONGITUDE",
    positiveWestLon(easternmost_longitude), 8, 12, "deg"));
  list.push_back(KeyValPair("WESTERNMOST_LONGITUDE",
    positiveWestLon(westernmost_longitude), 8, 12, "deg"));
  list.push_back(KeyValPair("LINE_PROJECTION_OFFSET",
    line_projection_offset, 3, NO_UNITS));
  list.push_back(KeyValPair("SAMPLE_PROJECTION_OFFSET",
    sample_projection_offset, 3, NO_UNITS));
  list.push_back(KeyValPair("OBLIQUE_PROJ_POLE_LATITUDE", proj.poleLatitude(),
    6, "deg"));
  list.push_back(KeyValPair("OBLIQUE_PROJ_POLE_LONGITUDE", proj.poleLongitude(),
    6, "deg"));
  list.push_back(KeyValPair("OBLIQUE_PROJ_POLE_ROTATION", proj.poleRotation(),
    6, "deg"));
  list.push_back(KeyValPair("OBLIQUE_PROJ_X_AXIS_VECTOR", rotation_matrix[0],
    3, 8, NO_UNITS));
  list.push_back(KeyValPair("OBLIQUE_PROJ_Y_AXIS_VECTOR", rotation_matrix[1],
    3, 8, NO_UNITS));
  list.push_back(KeyValPair("OBLIQUE_PROJ_Z_AXIS_VECTOR", rotation_matrix[2],
    3, 8, NO_UNITS));
  list.push_back(KeyValPair("LOOK_DIRECTION",
    PDSLabel::toUpper(look_direction_), !IS_ODL_STRING));
  list.push_back(KeyValPair("COORDINATE_SYSTEM_TYPE", coordinate_system_type_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("COORDINATE_SYSTEM_NAME", coordinate_system_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("END_OBJECT", "IMAGE_MAP_PROJECTION",
    !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("END", EMPTY_STRING, !IS_ODL_STRING));

  //
  // Create the PDS label in the proper format.
  //
  PDSLabel label(list, record_length*sample_bits_[bidr_type]/BITS_PER_BYTE,
    data_records, "^IMAGE");

  if (dbg.level)
  {
    if (bidr_type == LAT)
    {
      dbg.file << "BIDR::constructLabel:  PDS label for LAT file is" << endl;
      dbg.file << label.label() << endl;
    }
  }
  return(label.label());
}

//=====================================================================
// mapScale()
//
// Return the value of the MAP_SCALE keyword in the PDS label.
// The calculation assumes that Titan is a sphere; i.e., its x-, y-,
// and z-axes have the same lengths.
//=====================================================================
double
BIDR::mapScale()
{
  double target_radius =
    (default_target_radii.magnitude()).km()/sqrt(3.0);

  return(target_radius/(pixelsperdegree * radtodeg));
}

//=====================================================================
// productID()
//
// Determine the value of the PRODUCT_ID keyword in the BIDR PDS label.
// This value is dependent upon the BIDR data type, the map resolution,
// the (standard) latitude, longitude, and hemisphere at the center of
// the file, the data take ID, the (Titan) flyby ID, and the product
// version.
//=====================================================================
string
BIDR::productID(const BIDRTypeE bidr_type)
{
  string product_id = "BI";
  double lat_centerrad, lon_centerrad, slat_center, slon_center;
  int clat, clon;
  stringstream buf_stream(ios::app|ios::out);

  DebugInfo dbg("BIDR::productID");

  // BIDR data type
  // Naming conventions are included for the internal (nondeliverable)
  // products to avoid problems with the output filenames.
  switch(bidr_type)
  {
    case S0_CORR:
      product_id += "F";
      break;
    case S0_CORR_DB:
      product_id += "B";
      break;
    case S0_UNCORR:
      product_id += "U";
      break;
    case S0NSQT_UNCORR:
      product_id += "S";
      break;
    case S0_STD:
      product_id += "D";
      break;
    case S0_NOISE_EQUIV:
      product_id += "X";
      break;
    case INC:
      product_id += "E";
      break;
    case LAT:
      product_id += "T";
      break;
    case LON:
      product_id += "N";
      break;
    case BEAMMASK:
      product_id += "M";
      break;
    case START_BURST_NUM:
      product_id += "J";
      break;
    case END_BURST_NUM:
      product_id += "K";
      break;
    case NUM_LOOKS:
      product_id += "L";
      break;
    case NOMINAL_BURST_NUM:
    case RANGE:
    case DOPPLER:
      product_id+= "Z";
      break;
    default:
      ErrorMessage e("BIDR::productID:  unknown BIDR type " + toStr(bidr_type));
      e.throwMe();
      break;
  }
  product_id+= "Q";
  // Map resolution
  switch(pixelsperdegree)
  {
  case 1:
    product_id += "A";
    break;
  case 2:
    product_id += "B";
    break;
  case 4:
    product_id += "C";
    break;
  case 8:
    product_id += "D";
    break;
  case 16:
    product_id += "E";
    break;
  case 32:
    product_id += "F";
    break;
  case 64:
    product_id += "G";
    break;
  case 128:
    product_id += "H";
    break;
  case 256:
    product_id += "I";
    break;
  case 512:
    product_id += "J";
    break;
  default:
    product_id += "X";
  }

  // Determine the latitude, hemisphere (N/S), and longitude at the center
  // of the file.  The oblique cylindrical center latitude is the average of
  // the minimum and maximum oblique cylindrical latitudes in the file.  The
  // center longitude is calculated similarly.  (The calculation does not need
  // to account for wraparound, as the grid longitudes are either monotonically
  // increasing or monotonically decreasing.)  Convert the center latitude and
  // longitude to standard coordinates to obtain the characters needed for
  // the product ID.
  double maxlatrad = grid.latInRad(grid.numLats()-1);
  lat_centerrad = (grid.minLat() + maxlatrad)/2;
  double endlon;
  if (grid.lonIncreasing())
  {
    endlon = grid.startLon() + grid.totalNumLons()*degtorad/pixelsperdegree;
  }
  else
  {
    endlon = grid.startLon() - grid.totalNumLons()*degtorad/pixelsperdegree;
  }
  lon_centerrad = (grid.startLon() + endlon)/2;
  slat_center = radtodeg * proj.standardLatInRad(lon_centerrad, lat_centerrad);
  slon_center = radtodeg * proj.standardLonInRad(lon_centerrad, lat_centerrad);
  clat = round_double(slat_center);
  clon = round_double(positiveWestLon(slon_center));
  if (dbg.level) {
    dbg.file << "BIDR::productID:" << endl;
    dbg.file << "Center lat (deg):    " << slat_center << endl;
    dbg.file << "PW center lon (deg): " << positiveWestLon(slon_center) << endl;
    dbg.file << "Min lat (rad):       " << grid.minLat() << endl;
    dbg.file << "Max lat (rad):       " << maxlatrad << endl;
    dbg.file << "Start lon (rad):     " << grid.startLon() << endl;
    dbg.file << "End lon (rad):       " << endlon << endl;
  }
  buf_stream << std::setw(2) << std::setfill('0') << fabs(clat);
  if (clat > 0)
  {
    buf_stream << "N";
  } else
  {
    buf_stream << "S";
  }
  // clon is already in the range [0, 360) from proj.standardLonInRad()
  buf_stream << std::setw(3) << std::setfill('0') << clon;

  // Data take ID
  buf_stream << "_D" << std::setw(3) << std::setfill('0') << data_take_id_;

  // Flyby ID
  // flyby_id_pass_ is the string following the "T" in the flyby ID.
  // flyby_id_pass_ is set by BIDR::config, which checks the string for
  // validity, strips leading zeroes, and converts it to upper-case when
  // needed.
  buf_stream << "_T";
  for (unsigned int i = 0; i < 3-flyby_id_pass_.length(); i++)
  {
    buf_stream << "0";
  }
  buf_stream << flyby_id_pass_;
  buf_stream << "S" << std::setw(2) << std::setfill('0') << segment_id_;
  // Product version
  buf_stream << "_V" << std::setw(2) << std::setfill('0') << product_version_;

  product_id += buf_stream.str();
  return product_id;
}

//=====================================================================
// checksum()
//
// Return the value of the CHECKSUM keyword in the PDS label for the
// specified type of BIDR data.
//=====================================================================
int
BIDR::checksum(const BIDRTypeE bidr_type)
{
  return checksum_[bidr_type];
}

//----------------------------------------------------------------------
// filePrefix()
//
// Assumes UNIX-style directory specifications
//----------------------------------------------------------------------
string
BIDR::filePrefix(const string filename)
{
  size_t len = filename.length();
  size_t i, j;

  // Find the last occurrence of the directory separator in the input
  // filename.  If the string doesn't contain a separator, point to the
  // beginning of the string.
  if ((i = filename.rfind("/")) == string::npos)
  {
    i = 0;
  }
  else
  {
    i++;
  }
  string basename = filename.substr(i, len);
  len = basename.length();

  // Find the "." that terminates the file prefix.
  if ((j = basename.rfind(".")) == string::npos)
  {
    j = len;
  }

  return(basename.substr(0, j));
}

void BIDR::initProcMode(){
  for(int c=0;c<BIDRTYPECOUNT;c++) output_enable_[c]=false;
  switch(procMode){
  case NORMAL: 
    output_enable_[INC]=true;
    output_enable_[BEAMMASK]=true;
    output_enable_[S0_UNCORR]=true;
    output_enable_[S0NSQT_UNCORR]=true;
    output_enable_[S0_STD]=true;
    output_enable_[S0_NOISE_EQUIV]=true;
    output_enable_[S0_CORR]=true;
    output_enable_[LAT]=true;
    output_enable_[LON]=true;
    output_enable_[START_BURST_NUM]=true;
    output_enable_[END_BURST_NUM]=true;
    output_enable_[NUM_LOOKS]=true;
    break;
  case TOPO: 
  case TOPO_PLUS:
    output_enable_[INC]=true;
    output_enable_[BEAMMASK]=true;
    output_enable_[S0_UNCORR]=true;
    output_enable_[S0NSQT_UNCORR]=true;
    output_enable_[NUM_LOOKS]=true;
    break;
  case STEREO:
    output_enable_[BEAMMASK]=true;
    output_enable_[S0_UNCORR]=true;
    output_enable_[S0NSQT_UNCORR]=true;
    output_enable_[DOPPLER]=true;
    output_enable_[RANGE]=true;
    output_enable_[NOMINAL_BURST_NUM]=true;
    break;
  default:
    cerr << "Fatal error: BIDR has bad procMode" << endl;
    exit(1);
    break;
  }
}
