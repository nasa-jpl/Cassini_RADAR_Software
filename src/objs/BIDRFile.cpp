static const char rcs_id_bidrfile_c[] =
  "@(#) $Id: BIDRFile.cpp,v 11.7 2017/04/18 18:30:35 cveerama Exp $";

#include <iomanip>
#include <sstream>
#include <string>
#include "BIDRFile.h"
#include "DebugInfo.h"
#include "Error.h"
#include "PDSLabel.h"
#include "SimpleArray.h"

using std::cerr;
using std::cout;
using std::endl;
using std::setfill;
using std::setw;
using std::stringstream;

BIDRFile::BIDRFile()
: FileMgr(), header_handled_(false), bte_(INVALID),num_lats_(0),
  num_lons_(0), x_("x_",1), y_("y_",1), y2_("y2_",1), current_lon_idx_(-1),
  skip_invalid_(false),skip_low_(false), use_dB_(false), data_valid_(false),
  min_lon_(0), max_lon_(360)
{
}

BIDRFile::BIDRFile(const string& filename, const string& mode)
  : FileMgr(filename,mode), header_handled_(false), bte_(INVALID),num_lats_(0),
    num_lons_(0), x_("x_",1), y_("y_",1), y2_("y2_",1), current_lon_idx_(-1),
    skip_invalid_(false),skip_low_(false), use_dB_(false), data_valid_(false),
    min_lon_(0), max_lon_(360)
{
  // set up bad value float by specific bits.
  int val=BAD_VALUE_HEX;
  float* ptr=(float*) &val;
  BAD_VALUE=*ptr;
}

BIDRFile::~BIDRFile()
{
}

// must be called after readHeader and readLine
void BIDRFile::setDataValid(){

  data_valid_=true;

  // check range
  float lon=start_lon_;
  if(lon_inc_) lon+=(current_lon_idx_+0.5)*res_;
  else lon-=(current_lon_idx_+0.5)*res_;
  while(lon>max_lon_) lon-=360;
  while(lon<min_lon_) lon+=360;
  if(lon>max_lon_) data_valid_=false;

  // check for invalid or low valued cases and ignore if desired
  if(skip_invalid_){
    int num_bad=0;
    for(int c=0;c<num_lats_;c++){
      if(y_(c)==invalid_val_) num_bad++;
      else if(skip_low_ && y_(c)<=low_val_) num_bad++;
    }

    if(num_bad==num_lats_) data_valid_=false;
  }
}

void BIDRFile::readHeader()
{
  DebugInfo dbg("BIDRFile::readHeader");

  if(header_handled_){
    ErrorMessage e("BIDRFile attempt to read header twice.");
    e.throwMe();
  }
  string filename = name();
  int val;
  double map_resolution = 0;
  double line_projection_offset = 0;
  double sample_projection_offset = 0;

  product_id_ = EMPTY_STRING;
  PDSLabel header(*this);
  if (header.label() == EMPTY_STRING)
  {
    // PDS label not found in input file.  Assume that the file contains
    // the old-style binary header instead.
    cerr << "BIDRFile::readHeader:  " <<
      "Warning:  PDS label not found in " << filename << endl;

    read(min_lat_);
    read(start_lon_);
    read(pixelsperdegree);
    read(val);
    lon_inc_=(bool)val;
    read(num_lats_);
    read(num_lons_);
    read(val);
    bte_=(BIDRTypeE) val;
    // compute max values, res and x vector
    res_=1/(float)pixelsperdegree;
    d_res_=1/(double)pixelsperdegree;
  }
  else
  {
    label_length_ = header.labelLength();

    if (!header.getNumeric("LINE_SAMPLES", &num_lats_))
    {
      ErrorMessage e("Keyword LINE_SAMPLES not found in " + filename);
      e.throwMe();
    }
    if (!header.getNumeric("LINES", &num_lons_))
    {
      ErrorMessage e("Keyword LINES not found in " + filename);
      e.throwMe();
    }
    if (!header.getNumeric("MAP_RESOLUTION", &map_resolution))
    {
      ErrorMessage e("Keyword MAP_RESOLUTION not found in " + filename);
      e.throwMe();
    }
    if (!header.getString("PRODUCT_ID", &product_id_))
    {
      ErrorMessage e("Keyword PRODUCT_ID not found in " + filename);
      e.throwMe();
    }
    if (!header.getNumeric("LINE_PROJECTION_OFFSET", &line_projection_offset))
    {
      ErrorMessage e("Keyword LINE_PROJECTION_OFFSET not found in " +
	filename);
      e.throwMe();
    }
    if (!header.getNumeric("SAMPLE_PROJECTION_OFFSET",
      &sample_projection_offset))
    {
      ErrorMessage e("Keyword SAMPLE_PROJECTION_OFFSET not found in " +
	filename);
      e.throwMe();
    }
    if (!header.getNumeric("OBLIQUE_PROJ_POLE_LATITUDE", &pole_latitude_))
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_POLE_LATITUDE not found in " +
	filename);
      e.throwMe();
    }
    if (!header.getNumeric("OBLIQUE_PROJ_POLE_LONGITUDE", &pole_longitude_))
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_POLE_LONGITUDE not found in " +
	filename);
      e.throwMe();
    }
    if (!header.getNumeric("OBLIQUE_PROJ_POLE_ROTATION", &pole_rotation_))
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_POLE_ROTATION not found in " +
	filename);
      e.throwMe();
    }
    // Convert oblique cylindrical projection parameters from degrees
    // to radians and from positive-west to positive-east coordinates
    pole_latitude_ *= degtorad;
    pole_longitude_ = (360.0 - pole_longitude_) * degtorad;
    pole_rotation_ *= degtorad;

    // Fetch elements of the rotation matrix for use in conversions from
    // standard to oblique cylindrical coordinates.

    if (header.getVector("OBLIQUE_PROJ_X_AXIS_VECTOR", 3,
        rotation_matrix_[0]) != 3)
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_X_AXIS_VECTOR not found in " +
	filename);
      e.throwMe();
    }

    if (header.getVector("OBLIQUE_PROJ_Y_AXIS_VECTOR", 3,
        rotation_matrix_[1]) != 3)
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_X_AXIS_VECTOR not found in " +
	filename);
      e.throwMe();
    }

    if (header.getVector("OBLIQUE_PROJ_Z_AXIS_VECTOR", 3,
        rotation_matrix_[2]) != 3)
    {
      ErrorMessage e("Keyword OBLIQUE_PROJ_X_AXIS_VECTOR not found in " +
	filename);
      e.throwMe();
    }

    bte_ = getBIDRType();
    // compute max values, res and x vector
    if (map_resolution != 0)
    {
      res_ = 1 / map_resolution;
      d_res_ = 1 / map_resolution;
    }
    else
    {
      ErrorMessage e("Value of MAP_RESOLUTION in " + filename + " is 0.0");
      e.throwMe();
    }
    label_ = header.label();
    header.createListFromLabel(label_, label_list_);
    if (dbg.level)
    {
      dbg.file << "Label:" << endl;
      dbg.file << label_;
    }
  }

  pixelsperdegree = int(1/res_);
  double pixel_offset = 0.5/pixelsperdegree;
  min_lat_ = (-sample_projection_offset / map_resolution) - pixel_offset;
  start_lon_ = (-line_projection_offset / map_resolution) - pixel_offset;
  d_min_lat_ = (-sample_projection_offset / map_resolution) - pixel_offset;
  d_start_lon_ = (-line_projection_offset / map_resolution) - pixel_offset;
  while (start_lon_ < 0.0) start_lon_ += 360.0;
  while (start_lon_ >= 360.0) start_lon_ -= 360.0;

  max_lat_=min_lat_+res_*num_lats_;
  lon_inc_ = true; // Always true as per definition of map projection

  if (dbg.level)
  {
    dbg.file << "BIDRFile::readHeader:" << endl;
    dbg.file << "Min lat:    " << min_lat_ << endl;
    dbg.file << "Start lon:  " << start_lon_ << endl;
    dbg.file << "Pixels/deg: " << 1/res_ << endl;
    dbg.file << "Lon incr:   " << lon_inc_ << endl;
    dbg.file << "Num lats:   " << num_lats_ << endl;
    dbg.file << "Num lons:   " << num_lons_ << endl;
    dbg.file << "BIDR type:  " << bte_ << endl;
    dbg.file << "Resolution: " << res_ << endl;
    dbg.file << "Max lat:    " << max_lat_ << endl;
  }

  // due type specific handling
  switch(bte_){
  case INC:
    y_units_="degrees";
    type_string_="Incidence Angle";
    invalid_val_=-1;
    low_val_=0;
    break;
  case BEAMMASK:
    y_units_="";
    type_string_="Beam Mask";
    invalid_val_=0;
    low_val_=0;
    break;
  case S0_UNCORR:
    y_units_="dB";
    type_string_="Sigma0 Uncorrected";
    if(use_dB_){
      invalid_val_=-99;
      low_val_=-50;
    }
    else{
      invalid_val_=BAD_VALUE;
      low_val_=0;
    }
    break;
  case S0NSQT_UNCORR:
    y_units_="dB";
    type_string_="Sigma0 Uncorrected";
    if(use_dB_){
      invalid_val_=-99;
      low_val_=-50;
    }
    else{
      invalid_val_=BAD_VALUE;
      low_val_=0;
    }
    break;
  case S0_CORR:
    y_units_="dB";
    type_string_="Sigma0 Corrected";
    if(use_dB_){
      invalid_val_=0;
      low_val_=0;
    }
    else{
      invalid_val_=BAD_VALUE;
      low_val_=0;
    }
    break;
  case S0_NOISE_EQUIV:
    y_units_="";
    type_string_="Noise Equivalent Sigma0";
    if(use_dB_){
      invalid_val_=-99;
      low_val_=-50;
    }
    else{
      invalid_val_=BAD_VALUE;
      low_val_=0;
    }
    break;
  case S0_STD:
    y_units_="";
    type_string_="Sigma0 Standard Deviation";
    if(use_dB_){
      invalid_val_=-99;
      low_val_=-50;
    }
    else{
      invalid_val_=BAD_VALUE;
      low_val_=0;
    }
    break;
  case LAT:
    y_units_="degrees";
    type_string_="Latitude";
    invalid_val_=-91;
    low_val_=-91;
    break;
  case LON:
    y_units_="degrees";
    type_string_="Longitude";
    invalid_val_=-10000;
    low_val_=-10000;
    break;
  case START_BURST_NUM:
    y_units_="";
    type_string_="Starting Burst Number";
    invalid_val_=-1;
    low_val_=-1;
    break;
  case END_BURST_NUM:
    y_units_="";
    type_string_="Ending Burst Number";
    invalid_val_=-1;
    low_val_=-1;
    break;
  case NUM_LOOKS:
    y_units_="";
    type_string_="Number of Looks";
    invalid_val_=0;
    low_val_=0;
    break;
  default:
    ErrorMessage e("BIDRFile::readHeader:  " + toStr(bte_) +
      " is not a valid BIDR type");
    e.throwMe();
    break;
  }
  header_handled_=true;

  // allocate arrays
  x_.resize(num_lats_);
  y_.resize(num_lats_);
  y2_.resize(num_lats_);
  for(int c=0;c< num_lats_;c++){
    x_(c)=min_lat_+(c+0.5)*res_;
  }
  // Move the file pointer to the end of the header.
  setPosition(label_length_);
}

void BIDRFile::fixPole(double tca){
  int convert_from_iau=1;
  if(convert_from_iau){
  // compute MIAU
  double w=22.5769768; // IAU TITAN
  double theta0=189.64;
  double t0=0;
  double dtr=pi/180;

  double rtd=1/dtr;
  double second_to_day=1.0/(24.0*3600);
  double second_to_cent=1.0/(100.0*365.25*24*3600);
  double t=tca-t0;
  double T=t*second_to_cent;
  double S1=(29.80-52.1*T)*dtr;
  double ra= 36.41 -0.036*T+2.66*sin(S1);
  double dec= 83.94-0.004*T - 0.30*cos(S1);
  double d=t*second_to_day;
  double theta= 189.64 + w*d-2.64*sin(S1);
  theta=theta*dtr;

  double alpha=pi/2-dec*dtr;
  double beta=pi/2+ra*dtr;



  // rotate desired pole location to prime meridian
  double M3[3][3]={{cos(beta),sin(beta),0},{-sin(beta),cos(beta),0},{0,0,1}};
  // rotate pole to z
  double M2[3][3]={{1,0,0},{0,cos(alpha),sin(alpha)},{0,-sin(alpha),cos(alpha)}};
  // rotate prime meridian to its correct location
  double M1[3][3]={{cos(theta),sin(theta),0},{-sin(theta),cos(theta),0},{0,0,1}};
  mxm_c(M1,M2,M1);
  mxm_c(M1,M3,rotation_matrix_iau);
  //#define DEBUG
#ifdef DEBUG
  fprintf(stderr,"IAU ra=%g dec=%g theta=%g w=%g d=%g\n",ra,dec,theta,w,d);
#endif  
// compute M6para
  t0=207731850.598693073;
  t=tca-t0;
  ra=39.4826606;
  dec=83.4279446;
  w=22.57809919;
  double dradt=-30.1046278;
  double ddecdt=0;
  double dwdt=0.05231705;
  double cent=t/(365.25*100*24*3600);
  double cent0=t0/(365.25*100*24*3600);
  alpha=pi/2-dec*dtr-ddecdt*cent*dtr;
  beta=pi/2+ra*dtr+dradt*cent*dtr;
  double raout=ra+dradt*cent;
  double decout=dec+ddecdt*cent;

  d=t*second_to_day;

  S1=(29.80-52.1*cent0)*dtr;
  double raiau=36.41 -0.036*cent0+2.66*sin(S1);
  double fudgelon=0.01166991978; // forces longitude to be IAU(0,0) to be 0 at t0.
  double dra=raiau-ra+fudgelon;
  theta0=189.64+dra+22.5769768*t0*second_to_day -2.64*sin(S1);
  theta=theta0*dtr + (w +dwdt/2*cent)*d*dtr;
  double thetaout=fmod(theta/dtr,360);
  
  // rotate desired pole location to prime meridian
  double M3N[3][3]={{cos(beta),sin(beta),0},{-sin(beta),cos(beta),0},{0,0,1}};
  // rotate pole to z
  double M2N[3][3]={{1,0,0},{0,cos(alpha),sin(alpha)},{0,-sin(alpha),cos(alpha)}};
  // rotate prime meridian to its correct location
  double M1N[3][3]={{cos(theta),sin(theta),0},{-sin(theta),cos(theta),0},{0,0,1}};
  mxm_c(M1N,M2N,M1N);
  mxm_c(M1N,M3N,rotation_matrix_6para);

  // modify rot_mat
  #ifdef DEBUG
  fprintf(stderr,"6para ra=%g dec=%g theta=%g w=%g\n",ra,dec,theta,w);
#endif  
  double tmp[3][3], tmp2[3][3];
  invert_c(rotation_matrix_6para,tmp);
  mxm_c(rotation_matrix_,rotation_matrix_iau,tmp2);

#ifdef DEBUG
fprintf(stderr,"DEBUG OUTPUT FOR FIXPOLE for t=%23.9f\n",tca);
fprintf(stderr,"OLD ROTATION MATRIX\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n\n",  rotation_matrix_[0][0],rotation_matrix_[0][1],rotation_matrix_[0][2],rotation_matrix_[1][0],rotation_matrix_[1][1],rotation_matrix_[1][2],rotation_matrix_[2][0],rotation_matrix_[2][1],rotation_matrix_[2][2]);
fprintf(stderr,"IAU ROTATION MATRIX\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n\n",  rotation_matrix_iau[0][0],rotation_matrix_iau[0][1],rotation_matrix_iau[0][2],rotation_matrix_iau[1][0],rotation_matrix_iau[1][1],rotation_matrix_iau[1][2],rotation_matrix_iau[2][0],rotation_matrix_iau[2][1],rotation_matrix_iau[2][2]);
fprintf(stderr,"6 para ROTATION MATRIX\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n\n",  rotation_matrix_6para[0][0],rotation_matrix_6para[0][1],rotation_matrix_6para[0][2],rotation_matrix_6para[1][0],rotation_matrix_6para[1][1],rotation_matrix_6para[1][2],rotation_matrix_6para[2][0],rotation_matrix_6para[2][1],rotation_matrix_6para[2][2]);

#endif
  mxm_c(tmp2,tmp,rotation_matrix_); 
#ifdef DEBUG
fprintf(stderr,"NEW ROTATION MATRIX\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n%12.9f %12.9f %12.9f\n\n",  rotation_matrix_[0][0],rotation_matrix_[0][1],rotation_matrix_[0][2],rotation_matrix_[1][0],rotation_matrix_[1][1],rotation_matrix_[1][2],rotation_matrix_[2][0],rotation_matrix_[2][1],rotation_matrix_[2][2]);
 exit(1);
#endif

  /********  modify projection angles ********/
  // The Z basis vector gives the latitude (gamma_p) and longitude
  // (lambda_p) of the north pole in terms of the standard system.
  pole_latitude_ = atan2(rotation_matrix_[2][2],
    sqrt(pow(rotation_matrix_[2][0],2) + pow(rotation_matrix_[2][1],2)));
  pole_longitude_ = atan2(rotation_matrix_[2][1], rotation_matrix_[2][0]);
  
  // The X basis vector gives the reference latitude and longitude.
  double reference_latitude_ = atan2(rotation_matrix_[0][2],
    sqrt(pow(rotation_matrix_[0][0],2) + pow(rotation_matrix_[0][1],2)));
  double reference_longitude_ = atan2(rotation_matrix_[0][1], rotation_matrix_[0][0]);


  // Determine the pole rotation theta_p.  Use the equation for
  // tan(lambda_a + theta_p).  Set (gamma, lambda) to the reference
  // latitude and longitude; this point lies on the X basis vector
  // and thus lambda_a = 0.
  pole_rotation_ = atan2(cos(reference_latitude_)*sin(reference_longitude_ -
    pole_longitude_), sin(pole_latitude_)*cos(reference_latitude_)*
    cos(reference_longitude_ - pole_longitude_) -
    cos(pole_latitude_)*sin(reference_latitude_));
  }
}
void BIDRFile::readLine(){
  if(!header_handled_){
    ErrorMessage e("BIDRFile::readLine header not handled");
    e.throwMe();
  }
  current_lon_idx_++;
  char val;
  float fval;
  int ival;
  for(int c=0;c< num_lats_;c++){
    switch(bte_){
    case INC:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_CORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case LAT:
    case LON:
      read(fval);
      // convert s0 value to dB if necessary
      if(use_dB_ && (bte_ == S0_UNCORR || bte_ == S0_CORR || bte_ == S0NSQT_UNCORR) ){
	if(fval==BAD_VALUE) fval=-99;  // BAD VALUE  goes to -99 dB
        else if(fval==0) fval=-50;     // zero sigma0 goes to -50 dB
	else {
	  fval=10*log10(fval);
	  if(fval<-50)  fval=-50;        // low sigma0 goes to -50 dB  
	}

      }
      y_(c)=fval;
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      read(ival);
      y_(c)=(double)ival;
      break;
    case NUM_LOOKS:
    case BEAMMASK:
    case S0_CORR_DB:
      read(val);
      y_(c)=(double)val;
      break;
    default:
      ErrorMessage e("BIDRFile::readLine Invalid BUDR Type");
      e.throwMe();
      break;
    }
  }
  setDataValid();
}

void BIDRFile::slideShow(int lons_per_plot, int lon_skip, bool plot_average){

  if(!header_handled_){
    ErrorMessage e("BIDRFile::slideShow Need to read header first.");
    e.throwMe();
  }

  // Color string
  if(lons_per_plot>10){
    ErrorMessage e("BIDRFile::slideShow Greater than 10 lines per plot");
    e.throwMe();
  }
  string color_list[10]={"black","red","green","blue","brown","gray","violet","cyan","magenta","yellow"};

  // only plots the plots which it can complete
  double num_plots_f=floor((double)num_lons_/(double)lons_per_plot/(double)lon_skip);
  int num_plots=(int)num_plots_f;
  if(num_plots!=num_plots_f){
    cerr << "BIDRFile::slideShow Warning skipping incomplete plot ...";
  }

  // Main loop of slide show routine
  for(int c=0;c<num_plots;c++){
    plot_.clear();
//    plot_.setTool("xmgr");
    plot_.setTool("xmgrace");
    plot_.setXlabel("Latitude(degrees)");
    plot_.setYlabel(type_string_+"("+y_units_+")");
    plot_.setLegendTextSize(0.7);
    y2_=0;
    int num_valid_lines=0;
    for(int i=0;i<lons_per_plot;i++){
      int num_valid=0;
      for(int j=0;j<lon_skip;j++){
	readLine();
	if(data_valid_){
	  num_valid++;
	  if(plot_average) y2_+=y_;
	}
      }
      if(num_valid==0) continue; // skip plots of invalid data;
      if(plot_average) y2_/=num_valid;
      else y2_=y_;
      plot_.addXY(x_,y2_,line("solid",color_list[i]),sym("circle",color_list[i],0.5));
      float lon=start_lon_;
      if(lon_inc_) lon+=(current_lon_idx_+0.5)*res_;
      else lon-=(current_lon_idx_+0.5)*res_;
      while(lon<0) lon+=360;
      while(lon>360)lon-=360;
      plot_.addLegend("Lon="+toStr(lon));
      num_valid_lines++;
    }
    if(num_valid_lines) plot_.show("x");
  }

}

//=====================================================================
// getBIDRType()
//
// Decode the BIDR data type in a product ID.
//=====================================================================
BIDRTypeE
BIDRFile::getBIDRType()
{
  // Just examine the first three characters -- don't attempt to validate
  // the entire string
  if (product_id_.length() < 3)
  {
    ErrorMessage e("BIDRFile::getBIDRType:  invalid product ID " +
      product_id_);
    e.throwMe();
  }
  if (product_id_[0] != 'B' || product_id_[1] != 'I') {
    ErrorMessage e("BIDRFile::getBIDRType:  " + product_id_ +
      " is not a valid BIDR product ID");
    e.throwMe();
  }

  char bidr_type = product_id_[2];

  switch (bidr_type)
  {
    case 'F':
      return S0_CORR;
      break;
    case 'B':
      return S0_CORR_DB;
      break;
    case 'D':
      return S0_STD;
      break;
    case 'U':
      return S0_UNCORR;
      break;
    case 'S':
      return S0NSQT_UNCORR;
      break;
    case 'E':
      return INC;
      break;
    case 'T':
      return LAT;
      break;
    case 'N':
      return LON;
      break;
    case 'M':
      return BEAMMASK;
      break;
    case 'J':
      return START_BURST_NUM;
      break;
    case 'K':
      return END_BURST_NUM;
      break;
    case 'L':
      return NUM_LOOKS;
      break;
    case 'X':
      return S0_NOISE_EQUIV;
      break;
    default:
      ErrorMessage e("BIDRFile::getBIDRType:  " + toStr(bidr_type) +
         " is not a valid BIDR type");
      e.throwMe();
      break;
  }
  // Should never be reached
  return INVALID;
}

//=====================================================================
// getScaleFactor()
//
// Compare the new resolution (in pixels/km) with the resolution of the
// existing file and return the scale factor to use when averaging or
// subsampling the data.  Return -1 if the new resolution is not one of
// the allowed values for BIDR data or if the new resolution is greater
// than the existing resolution.
//=====================================================================
int BIDRFile::getScaleFactor(int new_resolution)
{
  if (new_resolution != 1 &&
      new_resolution != 2 &&
      new_resolution != 4 &&
      new_resolution != 8 &&
      new_resolution != 16 &&
      new_resolution != 32 &&
      new_resolution != 64 &&
      new_resolution != 128 &&
      new_resolution != 256 &&
      new_resolution != 512) return ((int) RESOLUTION_INVALID);

  if (new_resolution > pixelsperdegree) return ((int) RESOLUTION_TOO_HIGH);

  return(pixelsperdegree / new_resolution);
}

//=====================================================================
// standardLatInRad()
//
// This method duplicates the calculation in OblCylProj::standardLatInRad().
// standardLatInRad() has been placed in this module to avoid modifying
// Projections.{cpp,h} for the time being.
//=====================================================================
double BIDRFile::standardLatInRad(double lonrad, double latrad)
{
  double slatrad = asin( sin(pole_latitude_)*sin(latrad) -
    cos(pole_latitude_)*cos(latrad)*cos(lonrad + pole_rotation_) );
  return(slatrad);
}

//=====================================================================
// standardLonInRad()
//
// This method duplicates the calculation in OblCylProj::standardLonInRad().
// standardLonInRad() has been placed in this module to avoid modifying
// Projections.{cpp,h} for the time being.
//=====================================================================
double BIDRFile::standardLonInRad(double lonrad, double latrad)
{
  double slonrad = pole_longitude_ +
    atan2(cos(latrad)*sin(lonrad + pole_rotation_),
    sin(pole_latitude_)*cos(latrad)*cos(lonrad + pole_rotation_) +
    cos(pole_latitude_)*sin(latrad));
  while (slonrad < 0) { slonrad += 2*pi; }
  while (slonrad >= 2*pi) { slonrad -= 2*pi; }
  return(slonrad);
}

//=====================================================================
// OCLatFromGrid()
//=====================================================================
double BIDRFile::OCLatFromGrid(double sample)
{
  return((d_min_lat_ + 0.5*d_res_ + (sample-1)*d_res_ ) * degtorad);
}

//=====================================================================
// OCLonFromGrid()
//=====================================================================
double BIDRFile::OCLonFromGrid(double line)
{
  return((d_start_lon_ + 0.5*d_res_ + (line-1)*d_res_ ) * degtorad);
}

//=====================================================================
// latInRad()
//
// Given a latitude and longitude in the standard system, return the
// corresponding latitude in the oblique system.
// The input parameters and return value are in radians.
// This method duplicates the calculation in OblCylProj::latInRad().
//=====================================================================
double BIDRFile::latInRad(double slonrad, double slatrad)
{
  double latrad;
  double s[3], obl[3];

  s[0] = cos(slatrad)*cos(slonrad);
  s[1] = cos(slatrad)*sin(slonrad);
  s[2] = sin(slatrad);
  rotate_vector(s, const_cast<double(*)[3]>(rotation_matrix_), obl);
  latrad = atan2(obl[2], sqrt(pow(obl[0],2) + pow(obl[1],2)));
  return(latrad);
}

//=====================================================================
// lonInRad()
//
// Given a latitude and longitude in the standard system, return the
// corresponding longitude in the oblique system.
// The input parameters and return value are in radians.
// This method duplicates the calculation in OblCylProj::lonInRad().
//=====================================================================
double BIDRFile::lonInRad(double slonrad, double slatrad)
{
  double lonrad;
  double s[3], obl[3];

  s[0] = cos(slatrad)*cos(slonrad);
  s[1] = cos(slatrad)*sin(slonrad);
  s[2] = sin(slatrad);
  rotate_vector(s, const_cast<double(*)[3]>(rotation_matrix_), obl);
  lonrad = atan2(obl[1], obl[0]);
  while (lonrad < 0) { lonrad += 2*pi; }
  while (lonrad >= 2*pi) { lonrad -= 2*pi; }
  return(lonrad);
}

//=====================================================================
// OCLatToGrid()
//=====================================================================
double BIDRFile::OCLatToGrid(double slatrad)
{
  double sample;

  sample = 1 + ((slatrad * radtodeg) - d_min_lat_ - 0.5 * d_res_)/d_res_;
  return(sample);
}

//=====================================================================
// OCLonToGrid()
//=====================================================================
double BIDRFile::OCLonToGrid(double slonrad)
{
  double line;

  line = 1 + ((slonrad * radtodeg) - d_start_lon_ - 0.5 * d_res_)/d_res_;

  //** presume no BIDR file exceeds 270 degrees of arc **/
  //** this allows negatives for latlon just before start ***/
  //** of swath                                           ***/

  if(line>(270/d_res_)) line-=360/d_res_;
  
  return(line);
}

BIDRFileUpdates::BIDRFileUpdates(char *input_file, char *output_file,
  int flagged_ppd_new, int ppd_new)
: bf_in_(input_file, "r"),
  bf_out_(output_file, "w"),
  line_start_(1),
  line_end_(0),
  pixel_start_(1),
  pixel_end_(0),
  orig_line_start_(1),
  orig_line_end_(0),
  orig_pixel_start_(1),
  orig_pixel_end_(0),
  scale_factor_(1),
  truncated_(0)
{
  int remainder;

  bf_in_.readHeader();
  line_end_ = bf_in_.numLons();
  pixel_end_ = bf_in_.numLats();
  orig_line_end_ = bf_in_.numLons();
  orig_pixel_end_ = bf_in_.numLats();

  // no change in ppd case
  if (flagged_ppd_new==0 and ppd_new==0){
    ppd_=bf_in_.pixelsPerDegree();
    scale_factor_=1;
  }
  else{
    if (flagged_ppd_new)
      {
	scale_factor_ = bf_in_.getScaleFactor(ppd_new);
	if (scale_factor_ == (int) RESOLUTION_INVALID)
	  {
	    ErrorMessage e("Invalid value " + toStr(ppd_new) + " specified for " +
			   " new resolution.  Valid resolutions are 2, 8, 32, 128, 256.");
	    e.throwMe();
	  }
	if (scale_factor_ == (int) RESOLUTION_TOO_HIGH)
	  {
	    ErrorMessage e(string("New resolution must be lower than that of ") +
			   string("input file.  Valid resolutions are 2, 8, 32, 128, 256."));
	    e.throwMe();
	  }
      }
    ppd_ = ppd_new;


    remainder = (line_end_ - line_start_ + 1) % scale_factor_;
    if (remainder)
      {
	cerr << "Adjustments needed to line count" << endl;
	line_end_ -= remainder;
	cerr << "New end line = " << line_end_ << endl;
	truncated_++;
      }
    remainder = (pixel_end_ - pixel_start_ + 1) % scale_factor_;
    if (remainder)
      {
	cerr << "Adjustments needed to pixel count" << endl;
	pixel_end_ -= remainder;
	cerr << "New end pixel = " << pixel_end_ << endl;
	truncated_++;
      }
    
  }
  label_list_ = bf_in_.labelList();
  readInputToBuffer();
}

BIDRFileUpdates::~BIDRFileUpdates()
{
  bf_in_.close();
  bf_out_.close();
}

//=====================================================================
// average()
//=====================================================================
int BIDRFileUpdates::average()
{
  unsigned char **cbuf = (unsigned char **) NULL;
  int **ibuf = (int **) NULL;
  float **fbuf = (float **) NULL;

  switch(bf_in_.type())
  {
    case NUM_LOOKS:
    case BEAMMASK:
    case S0_CORR_DB:
      averageBuffer(cbuf);
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:    
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      averageBuffer(fbuf);
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      averageBuffer(ibuf);
      break;
    default:
      break;
  }
  return 1;
}


//=====================================================================
// flip()
//=====================================================================
int BIDRFileUpdates::flip()
{
  char **cbuf = (char **) NULL;
  int **ibuf = (int **) NULL;
  float **fbuf = (float **) NULL;

  switch(bf_in_.type())
  {
    case NUM_LOOKS:
    case BEAMMASK:
    case S0_CORR_DB:
      flipBuffer(cbuf);
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      flipBuffer(fbuf);
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      flipBuffer(ibuf);
      break;
    default:
      break;
  }
  return 1;
}

int BIDRFileUpdates::readInputToBuffer()
{
  int i, j;
  int numlats = bf_in_.numLats();
  int numlons = bf_in_.numLons();
  char **cbuf;
  int **ibuf;
  float **fbuf;

  switch(bf_in_.type())
  {
    case BEAMMASK:
    case NUM_LOOKS:
    case S0_CORR_DB:
      cbuf = (char **) make_array(sizeof(char), 2, numlons, numlats);
      bf_in_.setPosition(bf_in_.labelLength());
      for (i = 0; i < numlons; i++)
      {
        for (j = 0; j < numlats; j++) bf_in_.read(cbuf[i][j]);
      }
      buf_ = cbuf;
      buf_numlats_ = numlats;
      buf_numlons_ = numlons;
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      fbuf = (float **) make_array(sizeof(float), 2, numlons, numlats);
      bf_in_.setPosition(bf_in_.labelLength());
      for (i = 0; i < numlons; i++)
      {
        for (j = 0; j < numlats; j++) bf_in_.read(fbuf[i][j]);
      }
      buf_ = (char **) fbuf;
      buf_numlats_ = numlats;
      buf_numlons_ = numlons;
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      ibuf = (int **) make_array(sizeof(int), 2, numlons, numlats);
      bf_in_.setPosition(bf_in_.labelLength());
      for (i = 0; i < numlons; i++)
      {
        for (j = 0; j < numlats; j++) bf_in_.read(ibuf[i][j]);
      }
      buf_ = (char **) ibuf;
      buf_numlats_ = numlats;
      buf_numlons_ = numlons;
      break;
    default:
      break;
  }
  return(bf_in_.getPosition());
}

//=====================================================================
// averageBuffer()
//=====================================================================
void BIDRFileUpdates::averageBuffer(unsigned char **cbuf)
{
  char cval = 0;

  int i, j;
  int ioff, joff;
  int pixelsum = 0;
  int nvalid = 0;
  int *count;
  int numlats;
  int numlons;

  DebugInfo dbg("BIDRFileUpdates::averageBuffer(unsigned char **)");

  checksum_ = 0;

  if (cbuf != (unsigned char **) NULL) return;

  numlons = (line_end_-line_start_+1)/scale_factor_;
  numlats = (pixel_end_-pixel_start_+1)/scale_factor_;

  if (dbg.level)
  {
    dbg.file << "Line start, end = " << line_start_ << " " << line_end_ << endl;
    dbg.file << "Pixel start, end = " << pixel_start_ << " " << pixel_end_ << endl;
    dbg.file << "New longitude dimension = " << numlons << endl;
    dbg.file << "New latitude dimension  = " << numlats << endl;
  }

  cbuf = (unsigned char **) make_array(sizeof(char), 2, numlons, numlats);
  count = (int *) make_array(sizeof(int), 1, numlats);
  for (i = 0; i < numlons; i++) {for (j = 0; j < numlats; j++) cbuf[i][j] = 0;}
  for (j = 0; j < numlats; j++) count[j] = 0;

  for (i = line_start_-1; i < line_end_; i+=scale_factor_)
  {
    cerr << "Completed output row " << i << " out of " << line_end_-line_start_-1 << endl;
      for (j = pixel_start_-1; j < pixel_end_; j+=scale_factor_)
      {
	for (ioff = 0; ioff < scale_factor_; ioff++)
	  {

	    for (joff = 0; joff < scale_factor_; joff++)
	      {
		cval = buf_[i+ioff][j+joff];
		if (cval != (char) (*(bf_in_.invalidValue())))
		  {
		    nvalid++;
		    if(bf_in_.type()==BEAMMASK) pixelsum |= cval;
		    else pixelsum+= cval;
		  }
	      }
	  }
	count[(j-pixel_start_+1)/scale_factor_] += nvalid;
	if(bf_in_.type() == NUM_LOOKS && pixelsum > 255) pixelsum=255;
	if(bf_in_.type() == S0_CORR_DB) pixelsum/=nvalid;
	cbuf[(i-line_start_+1)/scale_factor_][(j-pixel_start_+1)/scale_factor_] = (unsigned char) pixelsum;
 	nvalid = 0;
	pixelsum = 0;
    }
    for (j = 0; j < numlats; j++) { count[j] = 0; }
  }
  for (i = 0; i < numlons; i++)
  {
    for (j = 0; j < numlats; j++)
    {
      checksum_ += cbuf[i][j];
    }
  }

  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_);
  buf_ = (char **) cbuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateLabel();
}

void BIDRFileUpdates::averageBuffer(int **ibuf)
{
  int ival = 0;
  int ** ibuf_ = (int**) buf_;
  int i, j;
  int ioff, joff;
  int pixelsum = 0;
  int nvalid = 0;
  int *count;
  int numlats;
  int numlons;


  DebugInfo dbg("BIDRFileUpdates::averageBuffer(int **)");
  if(bf_in_.type()==START_BURST_NUM) pixelsum=100000000;
  checksum_ = 0;

  if (ibuf != (int **) NULL) return;

  numlons = (line_end_-line_start_+1)/scale_factor_;
  numlats = (pixel_end_-pixel_start_+1)/scale_factor_;

  if (dbg.level)
  {
    dbg.file << "Line start, end = " << line_start_ << " " << line_end_ << endl;
    dbg.file << "Pixel start, end = " << pixel_start_ << " " << pixel_end_ << endl;
    dbg.file << "New longitude dimension = " << numlons << endl;
    dbg.file << "New latitude dimension  = " << numlats << endl;
  }

  ibuf = (int **) make_array(sizeof(int), 2, numlons, numlats);
  count = (int *) make_array(sizeof(int), 1, numlats);
  for (i = 0; i < numlons; i++) {for (j = 0; j < numlats; j++) ibuf[i][j] = 0;}
  for (j = 0; j < numlats; j++) count[j] = 0;

  for (i = line_start_-1; i < line_end_; i+=scale_factor_)
  {
    cerr << "Completed output row " << i << " out of " << line_end_-line_start_-1 << endl;
    for (j = pixel_start_-1; j < pixel_end_; j+=scale_factor_)
      {
	for (ioff = 0; ioff < scale_factor_; ioff++)
	  {

	    for (joff = 0; joff < scale_factor_; joff++)
	      {
		ival = ibuf_[i+ioff][j+joff];
		if (ival != (*(bf_in_.invalidValue())) )
		  {
		    nvalid++;
		    if(bf_in_.type()==START_BURST_NUM && pixelsum > ival)
		      pixelsum=ival;
		    else if(bf_in_.type()==END_BURST_NUM && pixelsum < ival)
		      pixelsum=ival;
		  }
	      }
	  }
	ibuf[(i-line_start_+1)/scale_factor_][(j-pixel_start_+1)/scale_factor_] = pixelsum;
	nvalid = 0;
	pixelsum = 0;
        if(bf_in_.type()==START_BURST_NUM) pixelsum=100000000;
      }
    for (j = 0; j < numlats; j++) { count[j] = 0; }
  }
  for (i = 0; i < numlons; i++)
  {
    for (j = 0; j < numlats; j++)
    {
      checksum_ += ibuf[i][j];
    }
  }

  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_*sizeof(int));
  buf_ = (char **) ibuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateLabel(sizeof(int));
}

void BIDRFileUpdates::averageBuffer(float **fbuf)
{
  int ibad=BAD_VALUE_HEX;
  float* ptr=(float*)&ibad;
  float BAD_VALUE=*ptr;

  float fval = 0;
  float ** fbuf_ = (float**) buf_;
  int i, j;
  int ioff, joff;
  double pixelsum = 0;
  int nvalid = 0;
  int *count;
  int numlats;
  int numlons;


  DebugInfo dbg("BIDRFileUpdates::averageBuffer(float **)");
  checksum_ = 0;

  if (fbuf != (float **) NULL) return;

  numlons = (line_end_-line_start_+1)/scale_factor_;
  numlats = (pixel_end_-pixel_start_+1)/scale_factor_;

  if (dbg.level)
  {
    dbg.file << "Line start, end = " << line_start_ << " " << line_end_ << endl;
    dbg.file << "Pixel start, end = " << pixel_start_ << " " << pixel_end_ << endl;
    dbg.file << "New longitude dimension = " << numlons << endl;
    dbg.file << "New latitude dimension  = " << numlats << endl;
  }

  fbuf = (float **) make_array(sizeof(float), 2, numlons, numlats);
  count = (int *) make_array(sizeof(int), 1, numlats);
  for (i = 0; i < numlons; i++) {for (j = 0; j < numlats; j++) fbuf[i][j] = 0;}
  for (j = 0; j < numlats; j++) count[j] = 0;

  for (i = line_start_-1; i < line_end_; i+=scale_factor_)
  {
    cerr << "Completed output row " << i << " out of " << line_end_-line_start_-1 << endl;
      for (j = pixel_start_-1; j < pixel_end_; j+=scale_factor_)
      {
	for (ioff = 0; ioff < scale_factor_; ioff++)
	  {

	    for (joff = 0; joff < scale_factor_; joff++)
	      {

		fval = fbuf_[i+ioff][j+joff];
 		if (fval != BAD_VALUE)
		  {
		    nvalid++;
		    pixelsum+=fval;
		  }
	      }
	  }
        if(nvalid!=0){
	  fbuf[(i-line_start_+1)/scale_factor_][(j-pixel_start_+1)/scale_factor_] = (float)pixelsum/nvalid;
	}
	else{
	  fbuf[(i-line_start_+1)/scale_factor_][(j-pixel_start_+1)/scale_factor_] = BAD_VALUE;
	}
	nvalid = 0;
	pixelsum = 0;
      }
    for (j = 0; j < numlats; j++) { count[j] = 0; }
  }

  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_*sizeof(float));
  buf_ = (char **) fbuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateLabel(sizeof(float));
}


//=====================================================================
// flipBuffer()
//=====================================================================
void BIDRFileUpdates::flipBuffer(char **cbuf)
{
  char ** cbuf_ = (char**) buf_;
  int i, j;
  int numlats;
  int numlons;



  if (cbuf != (char **) NULL) return;

  numlons = line_end_-line_start_+1;
  numlats = pixel_end_-pixel_start_+1;
  cbuf = (char **) make_array(sizeof(char), 2, numlons, numlats);
  for (i = 0; i < numlons; i++)
    {
      cerr << "Completed row " << i  << " out of " << numlons << endl;
      int i2 = numlons-i-1;
      for (j = 0; j < numlats; j++){
	cbuf[i][j]=cbuf_[i2][j];
      }
    }
  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_);
  buf_ = (char **) cbuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateFlippedLabel(sizeof(char));
}


//=====================================================================
// flipBuffer()
//=====================================================================
void BIDRFileUpdates::flipBuffer(float **fbuf)
{
  float ** fbuf_ = (float**) buf_;
  int i, j;
  int numlats;
  int numlons;



  if (fbuf != (float **) NULL) return;

  numlons = line_end_-line_start_+1;
  numlats = pixel_end_-pixel_start_+1;
  fbuf = (float **) make_array(sizeof(float), 2, numlons, numlats);
  for (i = 0; i < numlons; i++)
    {
      cerr << "Completed row " << i << " out of " << numlons << endl;
      int i2 = numlons-i-1;
      for (j = 0; j < numlats; j++){
	fbuf[i][j]=fbuf_[i2][j];
      }
    }
  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_*sizeof(float));
  buf_ = (char **) fbuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateFlippedLabel(sizeof(float));
}


//=====================================================================
// flipBuffer()
//=====================================================================
void BIDRFileUpdates::flipBuffer(int **ibuf)
{
  int ** ibuf_ = (int**) buf_;
  int i, j;
  int numlats;
  int numlons;


  if (ibuf != (int **) NULL) return;

  numlons = line_end_-line_start_+1;
  numlats = pixel_end_-pixel_start_+1;
  ibuf = (int **) make_array(sizeof(int), 2, numlons, numlats);
  for (i = 0; i < numlons; i++)
    {
      cerr << "Completed row " << i << " out of " << numlons << endl;
      int i2 = numlons-i-1;
      for (j = 0; j < numlats; j++){
	ibuf[i][j]=ibuf_[i2][j];
      }
    }
  // Replace the old buffer with the new data
  free_array((void*)buf_, 2, buf_numlons_, buf_numlats_*sizeof(int));
  buf_ = (char **) ibuf;
  buf_numlats_ = numlats;
  buf_numlons_ = numlons;
  updateFlippedLabel(sizeof(int));
}

//=====================================================================
// setProductID()
//
// Update a BIDR product ID for a new resolution and possible change
// to the center location.
//=====================================================================
void BIDRFileUpdates::setProductID()
{
  char new_char;

  if (product_id_.length() < 12)
  {
    ErrorMessage e("BIDRFileUpdates::setProductID:  invalid product ID " +
      product_id_);
    e.throwMe();
  }
  // Just examine the first two characters -- don't attempt to validate
  // the entire string
  if (product_id_[0] != 'B' || product_id_[1] != 'I') {
    ErrorMessage e("BIDRFileUpdates::setProductID:  " + product_id_ +
      " is not a valid BIDR product ID");
    e.throwMe();
  }

  switch (ppd_)
  {
    case 1:
      new_char = 'A';
      break;
    case 2:
      new_char = 'B';
      break;
    case 4:
      new_char = 'C';
      break;
    case 8:
      new_char = 'D';
      break;
    case 16:
      new_char = 'E';
      break;
    case 32:
      new_char = 'F';
      break;
    case 64:
      new_char = 'G';
      break;
    case 128:
      new_char = 'H';
      break;
    case 256:
      new_char = 'I';
      break;
    case 512:
      new_char = 'J';
      break;
    default:
      ErrorMessage e("BIDRFileUpdates::setProductID:  invalid resolution " +
         toStr(ppd_));
      e.throwMe();
      break;
  }
  product_id_[4] = new_char;

  // Update the center latitude and/or longitude if the area covered by the
  // new product differs from that of the original.  This segment of code
  // addresses the general case in which the new product covers only a
  // subset of the area covered by the original product and has edges that
  // do not fall on integer latitude or longitude boundaries (in the oblique
  // cylindrical grid).  For most (all?) cases in which this method will be
  // used, however, no changes will be needed to the center latitude and
  // longitude.
  if ( (orig_pixel_start_ != pixel_start_) ||
       (orig_pixel_end_ != pixel_end_) ||
       (orig_line_start_ != line_start_) ||
       (orig_line_end_ != line_end_) ||
       truncated_ ) {
    int clat, clon; // in degrees
    double lat_centerrad, lon_centerrad;
    double slat, slon; // in degrees
    double start_lon, end_lon; // in radians
    lat_centerrad = (bf_in_.OCLatFromGrid(pixel_end_) +
                     bf_in_.OCLatFromGrid(pixel_start_))/2;
    end_lon = bf_in_.OCLonFromGrid(line_end_);
    start_lon = bf_in_.OCLonFromGrid(line_start_);
    // Translate longitudes to [-pi, pi) to avoid discontinuity at 0 rad
    while (end_lon >= pi) { end_lon -= 2*pi; }
    while (end_lon < -pi) { end_lon += 2*pi; }
    while (start_lon >= pi) { start_lon -= 2*pi; }
    while (start_lon < -pi) { start_lon += 2*pi; }
    lon_centerrad = (start_lon + end_lon)/2;
    slat = radtodeg * bf_in_.standardLatInRad(lon_centerrad, lat_centerrad);
    slon = radtodeg * bf_in_.standardLonInRad(lon_centerrad, lat_centerrad);
    clat = round_double(slat);
    clon = round_double(BIDR::positiveWestLon(slon));

    // Construct the new product ID.
    string prefix = product_id_.substr(0, 5);
    string suffix = product_id_.substr(11);
    stringstream buf_stream(ios::app|ios::out);
    buf_stream << prefix << std::setw(2) << std::setfill('0') << fabs(clat);
    if (clat > 0)
    {
      buf_stream << "N";
    }
    else
    {
      buf_stream << "S";
    }
    buf_stream << std::setw(3) << std::setfill('0') << clon << suffix;
    product_id_ = buf_stream.str();
  }
}

//=====================================================================
// updateLabel()
//
// Update the PDS keywords that are affected by the new record size,
// image resolution, etc.
//=====================================================================
void BIDRFileUpdates::updateLabel(int pixelsize, bool flipped)
{
  int ival;
  double map_scale = 0.0;
  KeyValSeq header_keyvals;
  KeyValSeq new_keyvals;

  DebugInfo dbg("BIDRFileUpdates::updateLabel");

  // Create a new label by making a copy of the original product label.

  PDSLabel header(label_list_, buf_numlats_*pixelsize, buf_numlons_, "^IMAGE");
  string header_text = header.label();
  (void) header.createListFromLabel(header_text, header_keyvals);

  // Update the record counts and offset to the first data record.
  if (!header.getNumeric("FILE_RECORDS", &ival))
  {
    ErrorMessage e("Keyword FILE_RECORDS not found in " + bf_in_.name());
    e.throwMe();
  }
  new_keyvals.push_back(KeyValPair("FILE_RECORDS", ival));
  if (!header.getNumeric("LABEL_RECORDS", &ival))
  {
    ErrorMessage e("Keyword LABEL_RECORDS not found in " + bf_in_.name());
    e.throwMe();
  }
  new_keyvals.push_back(KeyValPair("LABEL_RECORDS", ival));
  if (!header.getNumeric("^IMAGE", &ival))
  {
    ErrorMessage e("Keyword ^IMAGE not found in " + bf_in_.name());
    e.throwMe();
  }
  new_keyvals.push_back(KeyValPair("^IMAGE", ival));

  // PDSLabel constructor doesn't automatically update RECORD_BYTES; must
  // set explicitly instead
  new_keyvals.push_back(KeyValPair("RECORD_BYTES", buf_numlats_*pixelsize));

  // Update the file dimensions, etc.
  new_keyvals.push_back(KeyValPair("LINES", buf_numlons_));
  new_keyvals.push_back(KeyValPair("LINE_SAMPLES", buf_numlats_));
  new_keyvals.push_back(KeyValPair("LINE_LAST_PIXEL", buf_numlons_));
  new_keyvals.push_back(KeyValPair("SAMPLE_LAST_PIXEL", buf_numlats_));
  new_keyvals.push_back(KeyValPair("CHECKSUM", checksum_));
  new_keyvals.push_back(KeyValPair("MAP_RESOLUTION",
    (double) bf_in_.pixelsPerDegree()/scale_factor_, 1, "pix/deg"));
  string product_creation_time = Time::getCurrentUtcDOYTime();
  new_keyvals.push_back(KeyValPair("PRODUCT_CREATION_TIME", product_creation_time,
    !IS_ODL_STRING));
  if (!header.getNumeric("MAP_SCALE", &map_scale))
  {
    ErrorMessage e("Keyword MAP_SCALE not found in " + bf_in_.name());
    e.throwMe();
  }
  new_keyvals.push_back(KeyValPair("MAP_SCALE", map_scale*scale_factor_,
    8, "km/pix"));

  // Update PRODUCT_ID with the new resolution and, if needed, a new
  // center latitude and/or longitude
  product_id_ = bf_in_.productID();
  setProductID();
  new_keyvals.push_back(KeyValPair("PRODUCT_ID", product_id_, !IS_ODL_STRING));

  // Update the line and sample projection offsets by deriving the starting
  // longitude and latitude from the original offsets, then recalculating the
  // offsets with the new resolution
  double pixel_offset_orig = 0.5/(double)bf_in_.pixelsPerDegree();
  double pixel_offset_new  = 0.5/((double)bf_in_.pixelsPerDegree()/scale_factor_);
  double line_projection_offset;
  double sample_projection_offset;
  double firstpix_orig;
  double firstpix_new;
  double start_lon;
  double min_lat;

  if (!header.getNumeric("LINE_PROJECTION_OFFSET", &line_projection_offset))
  {
    ErrorMessage e("Keyword LINE_PROJECTION_OFFSET not found in " + bf_in_.name());
    e.throwMe();
  }
  firstpix_orig = -line_projection_offset/bf_in_.pixelsPerDegree();
  start_lon = firstpix_orig - pixel_offset_orig;
  firstpix_new  = start_lon + pixel_offset_new;
  // Wrap longitude to range [-180, 180) degrees

  // modify the LINE_PROJECTION_OFFSET if buffer was flipped
  if(flipped){
    firstpix_new-=(float)(buf_numlons_-1)/(float)(bf_in_.pixelsPerDegree()/scale_factor_);
  }
  while (firstpix_new < -180.0) firstpix_new += 360.0;
  while (firstpix_new >= 180.0) firstpix_new -= 360.0;
  line_projection_offset = -firstpix_new * (bf_in_.pixelsPerDegree()/scale_factor_);
  new_keyvals.push_back(KeyValPair("LINE_PROJECTION_OFFSET", line_projection_offset,
    3, NO_UNITS));
  if (!header.getNumeric("SAMPLE_PROJECTION_OFFSET", &sample_projection_offset))
  {
    ErrorMessage e("Keyword SAMPLE_PROJECTION_OFFSET not found in " + bf_in_.name());
    e.throwMe();
  }
  firstpix_orig = -sample_projection_offset/bf_in_.pixelsPerDegree();
  min_lat = firstpix_orig - pixel_offset_orig;
  firstpix_new  = min_lat + pixel_offset_new;
  sample_projection_offset = -firstpix_new * (bf_in_.pixelsPerDegree()/scale_factor_);
  new_keyvals.push_back(KeyValPair("SAMPLE_PROJECTION_OFFSET", sample_projection_offset,
    3, NO_UNITS));

  // Recalculate the maximum/minimum latitude and easternmost/westernmost
  // longitude, using the original data as a basis for the calculations.
  // Do the calculations in radians using the positive-east convention for
  // longitude; then convert to degrees and positive-west at the end.
  int i, j;
  double di, dj;
  double maximum_latitude = -pi/2;
  double minimum_latitude = pi/2;
  double easternmost_longitude = -pi;
  double westernmost_longitude = pi;
  double first_longitude = -1000;

  double *lat;
  double lon;
  double slat, slon;
  lat = (double *) make_array(sizeof(double), 1, pixel_end_+1);
  for (j = 0; j < pixel_start_; j++) { lat[j] = 0.0; }
  for (j = pixel_start_; j <= pixel_end_; j++) {
    dj = j;
    lat[j] = bf_in_.OCLatFromGrid(dj);
  }
  for (i = line_start_; i <= line_end_; i++)
  {
    di = i;
    lon = bf_in_.OCLonFromGrid(di);
    for (j = pixel_start_; j <= pixel_end_; j++)
    {
      slat = bf_in_.standardLatInRad(lon, lat[j]);
      slon = bf_in_.standardLonInRad(lon, lat[j]);
      // Move longitude to [-pi, pi) to avoid discontinuity at 0 rad
      if(first_longitude==-1000){
	first_longitude=slon;
        easternmost_longitude=slon-pi;
        westernmost_longitude=slon+pi;      
      }
      while (slon >= first_longitude+pi) { slon -= 2*pi; }
      while (slon < first_longitude-pi) { slon += 2*pi; }
      if (slat > maximum_latitude) { maximum_latitude = slat; }
      if (slat < minimum_latitude) { minimum_latitude = slat; }
      if (slon > easternmost_longitude) { easternmost_longitude = slon; }
      if (slon < westernmost_longitude) { westernmost_longitude = slon; }
    }
  }
  // Move longitudes to [0, 2*pi)
  while (easternmost_longitude >= 2*pi) { easternmost_longitude -= 2*pi; }
  while (easternmost_longitude < 0)     { easternmost_longitude += 2*pi; }
  while (westernmost_longitude >= 2*pi) { westernmost_longitude -= 2*pi; }
  while (westernmost_longitude < 0)     { westernmost_longitude += 2*pi; }
  // Convert to degrees and positive-west
  maximum_latitude *= radtodeg;
  minimum_latitude *= radtodeg;
  easternmost_longitude *= radtodeg;
  westernmost_longitude *= radtodeg;
  easternmost_longitude = BIDR::positiveWestLon(easternmost_longitude);
  westernmost_longitude = BIDR::positiveWestLon(westernmost_longitude);
  new_keyvals.push_back(KeyValPair("MAXIMUM_LATITUDE", maximum_latitude,
    8, 12, "deg"));
  new_keyvals.push_back(KeyValPair("MINIMUM_LATITUDE", minimum_latitude,
    8, 12, "deg"));
  new_keyvals.push_back(KeyValPair("EASTERNMOST_LONGITUDE", easternmost_longitude,
    8, 12, "deg"));
  new_keyvals.push_back(KeyValPair("WESTERNMOST_LONGITUDE", westernmost_longitude,
    8, 12, "deg"));

  // Update the label with the new values
  label_list_ = updateValues(label_list_, new_keyvals);

  // Print new label for debugging
  if (dbg.level)
  {
    PDSLabel scratch(label_list_, buf_numlats_*pixelsize, buf_numlons_, "^IMAGE");
    if (scratch.label() != EMPTY_STRING) {
      dbg.file << "New PDS label:" << endl;
      dbg.file << scratch.label() << endl;
    }
  }

}


void BIDRFileUpdates::updateFlippedLabel(int pixelsize){
  bool flipped=true;
  unsigned char** cbuf;
  switch(bf_in_.type())
  {
    case BEAMMASK:
    case NUM_LOOKS:
    case S0_CORR_DB:
      cbuf = (unsigned char **) buf_;
      checksum_=0;
      for (int i = 0; i < buf_numlons_; i++)
      {
        for (int j = 0; j < buf_numlats_; j++) checksum_+=(int)cbuf[i][j];
      }
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      checksum_=0;
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      checksum_=0;
      break;
    default:
      break;
  }
  updateLabel(pixelsize,flipped);
}

int BIDRFileUpdates::writeBufferToOutput()
{
  int i, j;
  char **cbuf;
  int **ibuf;
  float **fbuf;
  int pixelsize=1;
  switch(bf_in_.type())
  {
    case BEAMMASK:
    case NUM_LOOKS:
    case S0_CORR_DB:
      pixelsize=sizeof(char);
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      pixelsize=sizeof(float);
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      pixelsize=sizeof(int);
      break;
    default:
      break;
  }
  PDSLabel header(label_list_, buf_numlats_*pixelsize, buf_numlons_, "^IMAGE");
  bf_out_.setPosition(0);
  bf_out_.write(header.label());
  switch(bf_in_.type())
  {
    case BEAMMASK:
    case NUM_LOOKS:
    case S0_CORR_DB:
      cbuf = (char **) buf_;
      for (i = 0; i < buf_numlons_; i++)
      {
        for (j = 0; j < buf_numlats_; j++) bf_out_.write(cbuf[i][j]);
      }
      break;
    case INC:
    case LAT:
    case LON:
    case S0_UNCORR:
    case S0NSQT_UNCORR:
    case S0_NOISE_EQUIV:
    case S0_STD:
    case S0_CORR:
      fbuf = (float **) buf_;
      for (i = 0; i < buf_numlons_; i++)
      {
        for (j = 0; j < buf_numlats_; j++) bf_out_.write(fbuf[i][j]);
      }
      break;
    case START_BURST_NUM:
    case END_BURST_NUM:
      ibuf = (int **) buf_;
      for (i = 0; i < buf_numlons_; i++)
      {
        for (j = 0; j < buf_numlats_; j++) bf_out_.write(ibuf[i][j]);
      }
      break;
    default:
      break;
  }
  return(bf_out_.getPosition());
}

