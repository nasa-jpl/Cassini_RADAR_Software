//----------------------------------------------------------------------------
//
// Beam.cpp
//
// This file contains method definitions for the Beam handling classes
// This class provides an interface to return antenna gain for a given
// direction vector.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_beam_c[] =
  "@(#) $Id: Beam.cpp,v 11.6 2012/09/27 20:01:58 richw Exp $";

//---------------
// Other includes
//---------------

#include <strings.h>
#include <string>
#include <math.h>
#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include "Error.h"
#include "Units.h"
#include "Beam.h"
#include "Array.h"
#include "Constants.h"
#include "DebugInfo.h"
#include "Utils.h"

using std::cout;
using std::cerr;
using std::endl;
using std::greater;

//---------------------------------------------------
// Static const member initialization for class Beam
//---------------------------------------------------

const unsigned int Beam::Nazim_ = 400;
const unsigned int Beam::Nelev_ = 1200;
const double Beam::beam_broadening_factor_ = 0.88;  

//--------------------------------------------
// Standard interface methods for class Beam
//--------------------------------------------

//--------------
// Constructors
//--------------

//--------------------------------------------------------------------------
// Beam()
//
// Dummy constructor sets some default values for the size and spacing
// of the sampled beam pattern.  No specific beam number, frame etc is
// specified.  These must be specified before this object can be used.
//---------------------------------------------------------------------------

Beam::Beam() 
  : azim_step_antenna_(0.01*degtorad,"rad"),
    elev_step_antenna_(0.01*degtorad,"rad"),
    azim_start_antenna_(-2.0*degtorad,"rad"),
    elev_start_antenna_(-6.0*degtorad,"rad"),
    azim_step_antenna_in_rad_(0.01*degtorad),
    elev_step_antenna_in_rad_(0.01*degtorad),
    azim_start_antenna_in_rad_(-2.0*degtorad),
    elev_start_antenna_in_rad_(-6.0*degtorad),
    compute_beam_pattern_(false),
    read_beam_pattern_(false),
    beam_parameter_set_(false),
    beam_widths_computed_(false),
    beam_num_set_(false),
    boresight_set_(false),
    pattern_loaded_(false),
    beam_number_(0),
    beam_width_azim_oneway_("beam_width_azim_oneway_"),
    beam_width_elev_oneway_("beam_width_elev_oneway_"),
    beam_width_azim_twoway_("beam_width_azim_twoway_"),
    beam_width_elev_twoway_("beam_width_elev_twoway_"),
    beam_gain_("beam_gain_"), 
    beam_max_gain_(0.0),
    pattern_filename_("none")
  {
    beam_gain_.resize(Nazim_,Nelev_);
  }

//--------------------------------------------------------------------------
// Beam(unsigned int beamnum, Config& cfg)
//
// Construct from config file for specified beam number (1-5).
// Beam pattern is either loaded from an external file,
// or computed from a sinc pattern (controlled by beam_pattern_source keyword)
// with widths specified in the config file.  See config method below.
//---------------------------------------------------------------------------

Beam::Beam(unsigned int beamnum, Config& cfg) 
  : azim_step_antenna_(0.01*degtorad,"rad"),
      elev_step_antenna_(0.01*degtorad,"rad"),
      azim_start_antenna_(-2.0*degtorad,"rad"),
      elev_start_antenna_(-6.0*degtorad,"rad"),
      azim_step_antenna_in_rad_(0.01*degtorad),
      elev_step_antenna_in_rad_(0.01*degtorad),
      azim_start_antenna_in_rad_(-2.0*degtorad),
      elev_start_antenna_in_rad_(-6.0*degtorad),
      compute_beam_pattern_(false),
      read_beam_pattern_(false),
      beam_parameter_set_(false),
      beam_widths_computed_(false),
      beam_num_set_(true),
      boresight_set_(false),
      pattern_loaded_(false),
      beam_number_(beamnum),
      beam_width_azim_oneway_("beam_width_azim_oneway_"),
      beam_width_elev_oneway_("beam_width_elev_oneway_"),
      beam_width_azim_twoway_("beam_width_azim_twoway_"),
      beam_width_elev_twoway_("beam_width_elev_twoway_"),
      beam_gain_("beam_gain_"), 
      beam_max_gain_(0.0),
      pattern_filename_("none")
  {
  beam_gain_.resize(Nazim_,Nelev_);
  config(cfg,beamnum);
  }

//--------------------------------------------------------------------------
// Beam(unsigned int beamnum, Config& cfg, read_pattern)
//
// Construct from config file for specified beam number (1-5).
// Beam pattern is either loaded from an external file (read_pattern true),
// or computed from a sinc pattern (read_pattern false) with widths specified
// in the config file.  See config method below.
// This method should only be used if the user does not want to follow
// the behavior specified in the config file by the beam_pattern_source keyword.
//---------------------------------------------------------------------------

Beam::Beam(unsigned int beamnum, Config& cfg, bool read_pattern) 
  : azim_step_antenna_(0.01*degtorad,"rad"),
      elev_step_antenna_(0.01*degtorad,"rad"),
      azim_start_antenna_(-2.0*degtorad,"rad"),
      elev_start_antenna_(-6.0*degtorad,"rad"),
      azim_step_antenna_in_rad_(0.01*degtorad),
      elev_step_antenna_in_rad_(0.01*degtorad),
      azim_start_antenna_in_rad_(-2.0*degtorad),
      elev_start_antenna_in_rad_(-6.0*degtorad),
      compute_beam_pattern_(!read_pattern),
      read_beam_pattern_(read_pattern),
      beam_parameter_set_(false),
      beam_widths_computed_(false),
      beam_num_set_(true),
      boresight_set_(false),
      pattern_loaded_(false),
      beam_number_(beamnum),
      beam_width_azim_oneway_("beam_width_azim_oneway_"),
      beam_width_elev_oneway_("beam_width_elev_oneway_"),
      beam_width_azim_twoway_("beam_width_azim_twoway_"),
      beam_width_elev_twoway_("beam_width_elev_twoway_"),
      beam_gain_("beam_gain_"), 
      beam_max_gain_(0.0),
      pattern_filename_("none")
  {
  beam_gain_.resize(Nazim_,Nelev_);
  config(cfg,beamnum,read_pattern);
  }

//-------------------------------------------------------------------
// config(cfg,beamnum)
//
// Read beam specification parameters from config file and setup the
// specified beam number.
// Pattern loading behavior is controlled by config keywords.
// For sinc patterns, the widths are loaded.
// For external patterns, the filenames are obtained from the config
// file.  The patterns are not actually loaded here.  Instead, they
// are loaded when they are first required (load on demand).
//------------------------------------------------------------------

void Beam::config(Config& cfg, unsigned int beamnum)
  {
  string beam_pattern_source = cfg.str("beam_pattern_source");
  if (beam_pattern_source == "file")
    {
    config(cfg,beamnum,true);
    }
  else if (beam_pattern_source == "sinc_model")
    {
    config(cfg,beamnum,false);
    }
  else
    {
    ErrorMessage e("beam_pattern_source is invalid");
    e.throwMe();
    }
  }

//-------------------------------------------------------------------
// config(cfg,beamnum,read_pattern)
//
// Read beam specification parameters from config file and setup the
// specified beam number.
// For sinc patterns (read_pattern false), the widths are loaded.
// For external patterns, the filenames are obtained from the config
// file.  The patterns are not actually loaded here.  Instead, they
// are loaded when they are first required (load on demand).
//------------------------------------------------------------------

void Beam::config(Config& cfg, unsigned int beamnum, bool read_pattern)
  {
  read_beam_pattern_ = read_pattern;
  compute_beam_pattern_ = !read_pattern;

  beam_number_ = beamnum;

  // Assign frame accordingly
  frame_ = Frame(beam_frame_spice_id[beam_number_-1], cassini_spice_id);

  // Extra frame used by special interface routines
  nominal_frame_ = frame_;

  char c[100];

  // Since pattern files are scaled to beam3 peak gain = unity 
  // if read_pattern is TRUE we read beam 3 peak gain in initially
  // it will be reset to the maximum gain in the pattern by loadBeamPattern.

  Uvar mindb=cfg["beam_mindB_gain"];
  beam_min_gain_= pow(10,mindb.getInUnits("")/10.0);
  if(read_pattern)
    {
    sprintf(c,"beam_maxdB_gain3");
    }
  else
    {
    sprintf(c,"beam_maxdB_gain%d",beam_number_);
    }

  Uvar x=cfg[c];
  beam_max_gain_= pow(10,x.getInUnits("")/10.0);
  beam_parameter_set_= true;

  read_beam_pattern_=read_pattern;
  compute_beam_pattern_=!read_pattern;
  if(read_beam_pattern_)
    {
    // Get beam pattern filename
    char c[20];
    sprintf(c,"beam_%d_pattern_file",beam_number_);
    pattern_filename_=cfg.str(c);
    }
  else
    {
    // Assign widths corresponding to beam number
    sprintf(c,"beam_azi_%d",beam_number_);
    beam_width_azim_oneway_=cfg[c];
    sprintf(c,"beam_elev_%d",beam_number_);
    beam_width_elev_oneway_=cfg[c];
    }

  spacecraft_ = cfg.str("spacecraft");
  }

//-----------------
// Get/Set methods
//-----------------

//-----------------------------------------------------------
// g = getGainOneWay(u)
//
// Return the one-way gain in the direction u for this beam.
// Linear value (not dB)
//-----------------------------------------------------------

double Beam::getGainOneWay(DirectionVector u)
  {
  if (beam_parameter_set_== false)
    {
    ErrorMessage e("Beam::getGainOneWay: Beam not configured yet");
    e.throwMe();
    }

  Uvar azim,elev;
  u.representIn(frame_);
  u.getAzimuthElevation(azim,elev);
  return(bilinear(azim,elev));
  }

//------------------------------------------------
// g = getMaxGain()
//
// Return the one-way max gain for this beam.
// Linear value (not dB)
//------------------------------------------------

double Beam::getMaxGain()
  {
  if (beam_parameter_set_== false)
    {
    ErrorMessage e("Beam::getMaxGain: No beam parameter is read in");
    e.throwMe();
    }

  if (!pattern_loaded_)
    {  // load on demand
    if (read_beam_pattern_)
      {
      loadBeamPatternFromFile();
      }
    else
      {
      computeSincBeamPattern();
      }
    }
  return(beam_max_gain_);
  }

//---------------------------------------------------------------------
// getBoresight(t)
//
// This method returns a DirectionVector constructed for time t with
// this beam boresight vector.  The vector will be constructed in the
// corresponding beam frame.  The time t governs how the vector will
// transform into other time dependent frames.
//---------------------------------------------------------------------

DirectionVector Beam::getBoresight(const Time& t) const
  {
  Frame fbeam("CASSINI_RADAR_" + toStr(beam_number_), spacecraft_);
  DirectionVector boresight("boresight",fbeam,t,0,0,1);
  return(boresight);
  }

//--------------------------------
// i = getBeamNumber()
//
// Get this Beam Number (1-5)
//-------------------------------

unsigned int Beam::getBeamNumber() const throw(ErrorMessage)
  {
  if (beam_num_set_ == false)
    {
    ErrorMessage e("Beam number has not been set");
    e.throwMe();
    }
  return (beam_number_);
  }

//-----------------------------------------
//Get beam width in the Azi angle direction
//-----------------------------------------

Uvar Beam::getAzimuthWidthOneWay()
  {
  
  if( beam_parameter_set_== false)
    {
    ErrorMessage e("No beam parameter is read in");
    e.throwMe();
    }

  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }

  return(beam_width_azim_oneway_);
  }

Uvar Beam::getAzimuthWidthTwoWay()
  {
  if( beam_parameter_set_== false)
    {
    ErrorMessage e("No beam parameter is read in");
    e.throwMe();
    }

  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }

  return(beam_width_azim_twoway_);
  }


//------------------------------------------
//Get beam width in the elev angle direction
//------------------------------------------

Uvar Beam::getElevationWidthOneWay()
  { 
  if( beam_parameter_set_== false)
    {
    ErrorMessage e("no beam parameter is read in");
    e.throwMe();
    }

  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }

  return( beam_width_elev_oneway_);
  }

Uvar Beam::getElevationWidthTwoWay()
  { 
  if( beam_parameter_set_== false)
    {
    ErrorMessage e("no beam parameter is read in");
    e.throwMe();
    }

  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }

  return( beam_width_elev_twoway_);
  }


//------------------------------------------
// Get beam frame
//------------------------------------------

const Frame& Beam::getFrame(Time& t) const throw(ErrorMessage)
  { 
  if( beam_parameter_set_== false)
    {
    throw ErrorMessage("Beam not configured, no frame available");
    }

  t=t_;
  return( frame_);
  }

void Beam::estimateBeamWidthsFromPattern() throw (ErrorMessage)
{
  // One-Way Azimuth 3 dB width
  Uvec azim=elevAndGainToAzim(0,-3);
  if(azim.size()==2){
    beam_width_azim_oneway_=azim(1)-azim(0);
  }
  else{
    throw ErrorMessage("Azimuth One-way Beam Width Undefined");
  }

  // Two-Way Azimuth 3 dB width
  azim=elevAndGainToAzim(0,-1.5);
  if(azim.size()==2){
    beam_width_azim_twoway_=azim(1)-azim(0);
  }
  else{
    throw ErrorMessage("Azimuth Two-way Beam Width Undefined");
  }
  
  // One-Way Elevation 3 dB width
  Uvec elev=azimAndGainToElev(0,-3);
  if(elev.size()==2){
    beam_width_elev_oneway_=elev(1)-elev(0);
  }
  else{
    throw ErrorMessage("Elevation One-way Beam Width Undefined");
  }

  // Two-Way Elevation 3 dB width  
  elev=azimAndGainToElev(0,-1.5);
  if(elev.size()==2){
    beam_width_elev_twoway_=elev(1)-elev(0);
  }
  else{
    throw ErrorMessage("Elevation Two-way Beam Width Undefined");
  }
  beam_widths_computed_=true;
}

//-----------------------------------------------------------
// loadBeamPatternFromFile
//
// This method reads an antenna pattern from an existing file
//  not been set up yet (supposedly provided by M. J. from 
// radiometetry measurements)
//-----------------------------------------------------------

void Beam::loadBeamPatternFromFile()
  {
    // If pattern file name is "none"
    if(strcasecmp(pattern_filename_.c_str(),"none")==0){
      throw ErrorMessage("Beam::loadBeamPatternFromFile: no file available");
    }

    
    FileMgr patfile(pattern_filename_.c_str(),"r");
    for(unsigned int i=0;i<Nazim_;i++){
      for(unsigned int j=0;j<Nelev_;j++){
	patfile.read(beam_gain_(i,j));
      }
    }

    // Beam patterns in files are relative to peak gain of Beam 3
    // beam_max_gain_ is initially set to peak gain of beam 3 and then
    // reset.
    beam_gain_*=beam_max_gain_;
    beam_max_gain_=0;
    for(unsigned int i=0;i<Nazim_;i++){
      for(unsigned int j=0;j<Nelev_;j++){
	if(beam_max_gain_ < beam_gain_(i,j)){
	  beam_max_gain_=beam_gain_(i,j);
	}
      }
    }
    read_beam_pattern_=true;
    compute_beam_pattern_=false;
    pattern_loaded_ = true;
    computeMaxGainLine();//requires to have set pattern_loaded_ flag=true
   
 }

Uvar Beam::maxGainLineToAzim(const Uvar& elev){
  if (!pattern_loaded_)
    {  // load on demand
      if (read_beam_pattern_)
	{
	  loadBeamPatternFromFile();
	}
      else
	{
	  computeSincBeamPattern();
	}
    }

  if (beam_num_set_ == false)
    {
      throw ErrorMessage("Beam number is not set");
    }
  if( beam_parameter_set_== false)
    {
      throw ErrorMessage("no beam parameter is read in");
    }
  static Uvar azim;
  azim=max_gain_line_azim_offset_+ max_gain_line_slope_*elev;
  return(azim);
}

void Beam::computeMaxGainLine(){
  Uvar delev(0.01*degtorad,"rad");  // elevation step size for samples
  Uvec elev_limits=azimAndGainToElev(0.0,-10.0); // get -10 dB elevations
                                                 // for 0 azimuth

  // make sure limits are in ascending order
  // set elev to minimum limit
  Uvar elev=elev_limits(0);
  Uvar azim;
  if(elev_limits(0)>elev_limits(1)){
    elev_limits(0)=elev_limits(1);
    elev_limits(1)=elev;
    elev=elev_limits(0);
  }

  // compute sampled (azim,elev) points on max gain curve
  // fit to line
  int N=0;
  Uvar sumAzim=0;
  Uvar sumElev=0;
  Uvar sumAzimElev=0;
  Uvar sumElevElev=0;
  while( elev< elev_limits(1)){
    azim=elevAndMaxGainToAzim(elev);
    N++;
    sumElev+=elev;
    sumAzim+=azim;
    sumAzimElev+=azim*elev;
    sumElevElev+=elev*elev;
    elev+=delev;
  }
  max_gain_line_azim_offset_=(sumAzimElev-sumAzim)/(sumElev-N);
  max_gain_line_slope_=(sumAzim- max_gain_line_azim_offset_*N)/(sumElev);
}

//----------------------------------
//fit beam contour line using best fit ellipse
//-----------------------------------
void Beam::computeBestFitEllipse(Uvec& azimuth, Uvec& elevation, const double& gain_db)
  {
    azimuth.resize(4);
    elevation.resize(4);
    //start with 200 points along the gain line
    unsigned int Np = 400;
    Uvec azim_input(" ",Np), elev_input("",Np);
    get1wayBeamGain(azim_input, elev_input, gain_db);
    //add more points at zero elevation and zero azimuth
    unsigned int extra_pts=100;//should be even, multiple of 2
    Uvec azim_fit(" ",Np+extra_pts), elev_fit(" ",Np+extra_pts);
    for (unsigned int i=0; i<azim_input.size();++i)
      {
	azim_fit(i)= azim_input(i);
	elev_fit(i)=elev_input(i);
      }
    Uvec dumm("",2);
    for (unsigned int i=0;i<extra_pts/2;++i)
      {
	dumm=azimAndGainToElev(Uvar(0,"rad"),gain_db);
	azim_fit(azim_input.size()+i)= Uvar(0,"rad");
	azim_fit(azim_input.size()+extra_pts/2+i)=Uvar(0,"rad");
	elev_fit(elev_input.size()+i)= dumm(0);
	elev_fit(elev_input.size()+extra_pts/2+i)=dumm(1);
      }
    //for (unsigned int i=0;i<extra_pts/4;++i)
    //{
    //dumm = elevAndGainToAzim(Uvar(0,"rad"),gain_db);
    //elev_fit(elev_input.size()+extra_pts/2+i)= Uvar(0,"rad");
    //elev_fit(elev_input.size()+extra_pts/2+extra_pts/4+i)=Uvar(0,"rad");
    //azim_fit(azim_input.size()+extra_pts/2+i)=dumm(0);
    //azim_fit(azim_input.size()+extra_pts/2+extra_pts/4+i)=dumm(1);
    //}
    Dvec x("x", azim_fit.size()),y("y",elev_fit.size());
    Dvec x_out("",4), y_out("",4);
    for(unsigned int i=0;i<azim_fit.size();++i)
      {
	x(i)=azim_fit(i).getInUnits("deg");
	y(i)= elev_fit(i).getInUnits("deg");
      }
    EllipseFit(x,y,x_out,y_out);
    for(unsigned int i=0;i<4;++i)
      {
	azimuth(i)= Uvar(x_out(i)*degtorad,"rad");
	elevation(i)=Uvar(y_out(i)*degtorad,"rad");
      }
    
  }
//-------------------------------------------------------------------------
// computeSincBeamPattern
//
// This method constructs an antenna beam pattern.
//-------------------------------------------------------------------------

void Beam::computeSincBeamPattern()
  { 
  if( beam_parameter_set_== false)
  {
    throw ErrorMessage("no beam parameter is read in");
  }
     
  //--------------------------------------------------------------------
  //Assign azi and elev angles from 
  //  azim_start to -azim_start:
  //  elev_start to -elev_start
  //--------------------------------------------------------------------
  
 
  
  //cout<<"Generating Antenna Beam Pattern"<<endl;

  for (unsigned int i = 0; i < Nazim_ ; ++i)
  {
    Uvar azi;
    azi = azim_start_antenna_ + i * azim_step_antenna_;

    for (unsigned int j = 0; j < Nelev_; ++j)
      {
	Uvar elev;
	elev = elev_start_antenna_ + j*elev_step_antenna_;
	Uvar au,bv;
	double sau,sbv;

	au = beam_broadening_factor_ * azi;
	au /= beam_width_azim_oneway_;

	bv = beam_broadening_factor_* elev;
	bv /= beam_width_elev_oneway_;

	sau = sinc(au.getValue()); 
	sbv = sinc(bv.getValue());

	beam_gain_(i,j) = sau * sau * sbv * sbv *beam_max_gain_; 
	if (beam_gain_(i,j) < 0.0) 
	  {
	    throw ErrorMessage("Error: beam gain is less than 0.0");
	  }
      }
  }   
  compute_beam_pattern_ = true;
  read_beam_pattern_ = false;
  pattern_loaded_ = true; 
  computeMaxGainLine();//require to set pattern_loaded_  = true
 

  //------------------------------------------------------
  // Display beam parameters on the screen
  //------------------------------------------------------
  //cout <<"beam width azi "<<beam_width_azim_<<endl;
  //cout << "beamwidth elev "<< beam_width_elev_<<endl;
  //cout << "maximum gain "<<beam_max_gain_<<endl;
  //cout <<" beam broadening factor "<<beam_broadening_factor_<<endl;  

  }

Uvec Beam::azimAndGainToElev(const Uvar& azim, double relative_gain_in_dB)
  {
  if (!pattern_loaded_)
    {  // load on demand
    if (read_beam_pattern_)
      {
      loadBeamPatternFromFile();
      }
    else
      {
      computeSincBeamPattern();
      }
    }

  if (beam_num_set_ == false)
    {
      throw ErrorMessage("Beam number is not set");
    }
  if( beam_parameter_set_== false)
    {
      throw ErrorMessage("no beam parameter is read in");
    }
  double target_gain=beam_max_gain_*(pow(10,0.1*relative_gain_in_dB));
  double azim_idx=getAzimuthIndex(azim);
  Dvec gains=beam_gain_.getInterpRow(azim_idx);
  Dvec elev_idx=gains.findInterpolatedIndices(target_gain);
  if(elev_idx.size()!=2 && DebugInfo::allWarnings)
    {
      cerr << "Warning: " << elev_idx.size() 
	   << " (Not 2) Different Elevations "
	   << "with relative gain " << relative_gain_in_dB << " dB "
	   << "and Azimuth " << azim << endl;
    }
  Uvec elev("elev",elev_idx.size());
  for(unsigned int c=0;c<elev_idx.size();c++){
    elev(c)=elev_start_antenna_+elev_step_antenna_*elev_idx(c);
  }

  return(elev);
  }


Uvec Beam::elevAndGainToAzim(const Uvar& elev, double relative_gain_in_dB)
{
  if (!pattern_loaded_)
    {  // load on demand
    if (read_beam_pattern_)
      {
      loadBeamPatternFromFile();
      }
    else
      {
      computeSincBeamPattern();
      }
    }

  if (beam_num_set_ == false)
    {
      throw ErrorMessage("Beam number is not set");
    }
  if( beam_parameter_set_== false)
    {
      throw ErrorMessage("no beam parameter is read in");
    }
  double target_gain=beam_max_gain_*(pow(10,0.1*relative_gain_in_dB));
  double elev_idx=getElevationIndex(elev);
  Dvec gains=beam_gain_.getInterpCol(elev_idx);
  Dvec azim_idx=gains.findInterpolatedIndices(target_gain);
  if(azim_idx.size()!=2 && DebugInfo::allWarnings)
    {
      cerr << "Warning: " << azim_idx.size() << " (Not 2) Different Azimuths "
	   << "with relative gain " << relative_gain_in_dB << " dB "
	   << "and Elevation " << elev << endl;
    }
  Uvec azim("azim",azim_idx.size());
  for(unsigned int c=0;c<azim_idx.size();c++){
    azim(c)=azim_start_antenna_+azim_step_antenna_*azim_idx(c);
  }
  return(azim);
}
    


Uvar Beam::elevAndMaxGainToAzim(const Uvar& elev)
  {
  if (!pattern_loaded_)
    {  // load on demand
      if (read_beam_pattern_)
	{
	loadBeamPatternFromFile();
	}
      else
	{
	computeSincBeamPattern();
	}
    }
    
  if (beam_num_set_ == false)
    {
      throw ErrorMessage("Beam number is not set");
    }
  if( beam_parameter_set_== false)
    {
      throw ErrorMessage("no beam parameter is read in");
    }   
  double gain_db=10.0*log(bilinear(Uvar(0,"rad"),elev)/getMaxGain())/log(10.0);
  if (gain_db <= -15 ) ErrorMessage("Beam.cpp::elevAndMaxGainToAzim:for a given elevation angle, beam gain is lower than -15dB ").throwMe();

  //do not use //
  //Dvec gains=beam_gain_.getInterpCol(elev_idx);
  //unsigned int azi_idx;
  //gains.max(azi_idx);//find maximum
  
  Uvec azim("azimuth ",2);
    if (gain_db > -3.0)
      {
	azim=elevAndGainToAzim(elev,-3.0);
      }
    else if (gain_db <= -3.0 && gain_db > -5.0) 
      {
	azim=elevAndGainToAzim(elev,-5.0);
      }
    else if (gain_db <=-5.0 && gain_db > -10.0)
      {
	azim=elevAndGainToAzim(elev,-10.0);
      }
    else
      {
	azim=elevAndGainToAzim(elev,-15.0);
      }
   
    //need special care if there are more than 2 values returned
    Uvar return_value=Uvar(0,"rad");
    if(azim.size() == 2)
      {
	return_value =(azim(0)+azim(1))/2 ;
      }
    else if(azim.size() >2)
      {
	for(unsigned int i=0;i<azim.size();++i)
	  {
	    return_value+= azim(i);
	  }
	return_value/=double(azim.size());
      }
    else
      {
	ErrorMessage("Beam.cpp::elevAndMaxGainToAzim:no azimuth angles are returned ").throwMe();
      }
    return(return_value);
  

  }		
		


//-------------------------------------------
//Get One-way  beam gain points (azimuth, elevation)
//-------------------------------------------
void Beam::get1wayBeamGain(Uvec& azimuth, Uvec& elevation, const double& gain_db)
  {
    if (beam_num_set_ == false)
      {
	throw ErrorMessage("Beam number is not set");
      }
    if( beam_parameter_set_== false)
      {
	throw ErrorMessage("no beam parameter is read in");
      }
    if (!beam_widths_computed_)
      {
      estimateBeamWidthsFromPattern();
      }
    if (gain_db < - 15)
      {
	throw ErrorMessage("beam.cpp::get1wayBeamGain: get1wayBeamGain works for gain values larger than -15db");
      }

    unsigned int azim_size = azimuth.size();
    unsigned int elevation_size = elevation.size();
    if (azim_size != elevation_size)
      {
	throw ErrorMessage("Two arrays have different size");
      }
    if (azim_size%4 != 0)
      {
	cout<<"azimuth and elevation array sizes will be adjusted to multiple of 8 "<<endl;
	unsigned int n = (unsigned int) azim_size/8;
	n *= (n+1)*8;
	azimuth.resize(n);
	elevation.resize(n);
      }
    unsigned int half_size = elevation.size()/2;
    unsigned int quarter_size = elevation.size()/4;
    unsigned int eighth_size = elevation.size()/8;

    Uvar zero_degree = Uvar(0*degtorad,"rad");
    Uvar elev_top, elev_bottom;
    Uvar azim_top,azim_bottom;

    
    Uvec two_elements("two elements array", 2);//dummy array  
    two_elements=azimAndGainToElev(zero_degree,gain_db);
    min_max(elev_bottom, elev_top, two_elements);//fix min and max
    elev_bottom += Uvar(elev_step_antenna_in_rad_,"rad");
    elev_top -=Uvar(elev_step_antenna_in_rad_,"rad");
    two_elements=elevAndGainToAzim(zero_degree,gain_db);
    min_max(azim_bottom,azim_top, two_elements);//fix min and max
    azim_bottom += Uvar(azim_step_antenna_in_rad_,"rad");
    elev_top -=Uvar(azim_step_antenna_in_rad_,"rad");

    //debug
    //cout<<elev_top.getInUnits("deg")<< " "<<elev_bottom.getInUnits("deg")<<endl;
    //cout<<azim_top.getInUnits("deg")<<" "<<azim_bottom.getInUnits("deg")<<endl;
    Uvar elev_start = zero_degree;
    for (unsigned int i = 0; i< eighth_size;++i)
      {//get two azimuth angles for each elevation angle between 0 and elev_top
	Uvar x = elev_start + (elev_top - elev_start)*double(i)/double(eighth_size-1);
	two_elements = zero_degree;//reset the container
	two_elements = elevAndGainToAzim(x,gain_db);
	elevation(i) = x;
	elevation(i+eighth_size)= x;
	min_max(azimuth(i+eighth_size),azimuth(i), two_elements);
      }
    elev_start = -Uvar(elev_step_antenna_in_rad_,"rad");
    for (unsigned int i = 0; i< eighth_size;++i)
      {//get two azimuth angles for each elevation angle between 0 and elev_bottom
	Uvar x = elev_start +(elev_bottom-elev_start)*double(i)/double(eighth_size-1);
	two_elements = zero_degree;//reset the container
	two_elements = elevAndGainToAzim(x,gain_db);
	elevation(i+quarter_size) = x;
	elevation(i+quarter_size+eighth_size)= x;
	min_max(azimuth(i+quarter_size),azimuth(i+quarter_size+eighth_size),two_elements);
      }
    
    Uvar azim_start = zero_degree;
    for(unsigned int i=0; i<eighth_size;++i)
      {//get two elevation angles for each azimuth betwen 0 and azi_top
	Uvar x = azim_start + (azim_top - azim_start)*double(i)/double(eighth_size-1);
	two_elements = zero_degree;//reset the container
	two_elements = azimAndGainToElev(x,gain_db);
	azimuth(half_size+i) = x;
	azimuth(half_size+eighth_size+i) = x;
	min_max(elevation(half_size+i),elevation(half_size+eighth_size+i),two_elements);
      }
    azim_start= -Uvar(azim_step_antenna_in_rad_,"rad");
    for(unsigned int i=0; i<eighth_size;++i)
      {//get two elevation angles for each azimth between 0 and azi_bottom
	Uvar x = azim_start+(azim_bottom-azim_start)*double(i)/double(eighth_size-1);
	two_elements = zero_degree;//reset the container
	two_elements = azimAndGainToElev(x,gain_db);
	azimuth(half_size+quarter_size+i) = x;
	azimuth(half_size+quarter_size+eighth_size+i) = x;
	min_max(elevation(half_size+quarter_size+i),elevation(half_size+quarter_size+eighth_size+i),two_elements);
      }

    //orderly arranged container
    vector<Uvar> elevation_output, azimuth_output;   
    elevation_output.clear();
    azimuth_output.clear();
   
    for(unsigned int i = 0; i < 2 ;++i)
    for(unsigned int j = 0; j < 2 ;++j)
    {
      vector<Uvar> dummy;
      dummy.clear();
      vector<Uvar> elev_q, azim_q;
      elev_q.clear();
      azim_q.clear();
      map<Uvar,Uvar> elev_azim_map;
      elev_azim_map.clear();
      map<Uvar,Uvar>::const_iterator map_pos;

      if(i ==0 && j == 0)
	{//first quadrant
	  for (unsigned int i=0; i<elevation.size();++i)
	    {
	      if(elevation(i)>= zero_degree && azimuth(i)>=zero_degree)
		{
		  dummy.push_back(elevation(i));
		  elev_q.push_back(elevation(i));
		  azim_q.push_back(azimuth(i));
		  elev_azim_map[elevation(i)]= azimuth(i);
		}
	    }
	  sort(dummy.begin(),dummy.end());//increasing elevation angle
	  for(unsigned int k=0;k<dummy.size();++k)
	    {
	      map_pos = elev_azim_map.find(dummy[k]);
	      if(map_pos == elev_azim_map.end()) ErrorMessage("Beam.cpp::get1wayBeamgain, map finding error").throwMe();
	      elevation_output.push_back(map_pos->first);
	      azimuth_output.push_back(map_pos->second);
	      //cout<<"elev and azi "<<map_pos->first<<" "<< map_pos->second<<endl;
	    }
	}
      
      else if(i ==0 && j == 1)
	{//second quadrant
	  for (unsigned int i=0; i<elevation.size();++i)
	    {
	      if(elevation(i)>= zero_degree && azimuth(i)<zero_degree)
		{
		  dummy.push_back(elevation(i));
		  elev_q.push_back(elevation(i));
		  azim_q.push_back(azimuth(i));
		  elev_azim_map[elevation(i)]= azimuth(i);
		}
	    }
	  sort(dummy.begin(),dummy.end(), greater<Uvar>());//decreasing elevation angle
	  for(unsigned int k=0;k<dummy.size();++k)
	    {
	      map_pos = elev_azim_map.find(dummy[k]);
	      if(map_pos == elev_azim_map.end()) ErrorMessage("Beam.cpp::get1wayBeamgain, map finding error").throwMe();
	      elevation_output.push_back(map_pos->first);
	      azimuth_output.push_back(map_pos->second);
	      //cout<<"elev and azi "<<map_pos->first<<" "<< map_pos->second<<endl;
	    }
	}
      
      else if(i ==1 && j == 0)
	{//third  quadrant
	  for (unsigned int i=0; i<elevation.size();++i)
	    {
	      if(elevation(i)< zero_degree && azimuth(i)<zero_degree)
		{
		  dummy.push_back(elevation(i));
		  elev_q.push_back(elevation(i));
		  azim_q.push_back(azimuth(i));
		  elev_azim_map[elevation(i)]= azimuth(i);
		}
	    }
	  sort(dummy.begin(),dummy.end(), greater<Uvar>());//decrasing elevation angle
	  for(unsigned int k=0;k<dummy.size();++k)
	    {
	      map_pos = elev_azim_map.find(dummy[k]);
	      if(map_pos == elev_azim_map.end()) ErrorMessage("Beam.cpp::get1wayBeamgain, map finding error").throwMe();
	      elevation_output.push_back(map_pos->first);
	      azimuth_output.push_back(map_pos->second);
	      //cout<<"elev and azi "<<map_pos->first<<" "<< map_pos->second<<endl;
	    }
	}
      else if (i==1 && j==1)
	{
	  for (unsigned int i=0; i<elevation.size();++i)
	    {
	      if(elevation(i)< zero_degree && azimuth(i) >= zero_degree)
		{
		  dummy.push_back(elevation(i));
		  elev_q.push_back(elevation(i));
		  azim_q.push_back(azimuth(i));
		  elev_azim_map[elevation(i)]= azimuth(i);
		}
	    }
	  sort(dummy.begin(),dummy.end());//increasing  elevation angle
	  for(unsigned int k=0;k<dummy.size();++k)
	    {
	      map_pos = elev_azim_map.find(dummy[k]);
	      if(map_pos == elev_azim_map.end()) ErrorMessage("Beam.cpp::get1wayBeamgain, map finding error").throwMe();
	      elevation_output.push_back(map_pos->first);
	      azimuth_output.push_back(map_pos->second);
	      //cout<<"elev and azi "<<map_pos->first<<" "<< map_pos->second<<endl;
	    }
	}
      else
	{
	  ErrorMessage("No quadrant assigned").throwMe();
	}
    }
    if( elevation_output.size()!= elevation.size()) ErrorMessage("Beam.cpp::get1wayBeamGain(): output size mismatch").throwMe();
    for(unsigned int i=0;i<elevation_output.size();++i)
      {//rearrange output arrays
	elevation(i) = elevation_output[i];
	azimuth(i)=azimuth_output[i];
      } 
  } 


double Beam::getAzimuthIndex(const Uvar& azim) const
{
  Uvar tmp= (azim - azim_start_antenna_)/ azim_step_antenna_;
  return(tmp.getValue());
}

double Beam::getElevationIndex(const Uvar& elev) const
{
  Uvar tmp=(elev - elev_start_antenna_)/elev_step_antenna_;
  return(tmp.getValue());
} 

void Beam::integrateBeamPattern(double& solid_angle, double& directivity_dB, double& beam_integral){
  
  if (!pattern_loaded_)
    {  // load on demand
      if (read_beam_pattern_)
	{
	  loadBeamPatternFromFile();
	}
      else
	{
	  computeSincBeamPattern();
	}
    }

  double d2r=degtorad;
 
  
  double phi,theta,x,y,z;
  double azi, elev,gain;

  //theta span, step
  double theta_start=0.0;//deg
  double theta_end= -elev_start_antenna_.getInUnits("deg");//deg
  unsigned int Num_theta_step=1000;
  double theta_step= (theta_end - theta_start)/float(Num_theta_step);
	
  //phi span and setp
  double phi_start=-180;
  double phi_end= 180;
  unsigned int  Num_phi_step=360;
  double phi_step=(phi_end - phi_start)/float(Num_phi_step);
	

  beam_integral=0.0;//reset
  solid_angle=0.0;
  directivity_dB=-100;//reset

  for(unsigned int j=0;j<Num_theta_step;++j){
    for(unsigned int k=0;k<Num_phi_step;++k){
      theta= theta_start +float(j)*theta_step;//0 to 6 degree
      phi=  phi_start+ float(k)*phi_step;//180 to 180
      
      theta= theta * d2r;
      phi= phi*d2r;
      
      x= sin(theta)*cos(phi);
      y= sin(theta)*sin(phi);
      z= cos(theta);
      
      elev= asin(y);
      azi=asin(x/cos(elev));
      
      gain=bilinear(azi,elev);
	    
      solid_angle=solid_angle+sin(theta)*( theta_step*d2r)*(phi_step*d2r);
      beam_integral= beam_integral + gain*sin(theta)*( theta_step*d2r)*(phi_step*d2r);
	   
    }
  }
  directivity_dB= 10.0*log( 4.0*pi/(beam_integral/getMaxGain()))/log(10.0);
}





//-----------------------------------------------------
//bilinear: return beam gain for a given set of azi and 
// elev angles
//-----------------------------------------------------

double Beam::bilinear(const Uvar& azim, const Uvar& elev)
  {
  if (!pattern_loaded_)
    {  // load on demand
    if (read_beam_pattern_)
      {
      loadBeamPatternFromFile();
      }
    else
      {
      computeSincBeamPattern();
      }
    }

  double a, b;   
  a=getAzimuthIndex(azim);
  b=getElevationIndex(elev);

  double gain=0.0;
  // attempt generic bilinear interpolation   
  try
    {
    gain = ::bilinear(a,b,beam_gain_); 
    }
  // on out of range failure set gain to zero
  catch(ErrorMessage& e)
    {
    gain=0;
    }

  if(gain==0) gain=beam_min_gain_;
  return(gain);
  }


//-----------------------------------------------------
//bilinear: return beam gain for a given set of azi and 
// elev angles
//-----------------------------------------------------

double Beam::bilinear(double azim_in_rad, double elev_in_rad)
  {
  if (!pattern_loaded_)
    {  // load on demand
    if (read_beam_pattern_)
      {
      loadBeamPatternFromFile();
      }
    else
      {
      computeSincBeamPattern();
      }
    }

  double a, b;   
  a=(azim_in_rad-azim_start_antenna_in_rad_)/ azim_step_antenna_in_rad_;
  b=(elev_in_rad-elev_start_antenna_in_rad_)/ elev_step_antenna_in_rad_;


  double gain=0.0;
  // attempt generic bilinear interpolation   
  try
    {
    gain = ::bilinear(a,b,beam_gain_); 
    }
  // on out of range failure set gain to zero
  catch(ErrorMessage& e)
    {
    gain=0;
    }

  if(gain==0) gain=beam_min_gain_;
  return(gain);
  }


//--------------------------------------------
// Special interface methods for class Beam
//--------------------------------------------

//---------------
// Constructors
//---------------

//--------------------------------------------------------------------------
// Beam(unsigned int beamnum, Config& cfg,Time t)
//
//---------------------------------------------------------------------------

Beam::Beam(unsigned int beamnum, Config& cfg, bool read_pattern, Time t) 
    :azim_step_antenna_(0.01*degtorad,"rad"),
     elev_step_antenna_(0.01*degtorad,"rad"),
     azim_start_antenna_(-2.0*degtorad,"rad"),
     elev_start_antenna_(-6.0*degtorad,"rad"),
     azim_step_antenna_in_rad_(0.01*degtorad),
     elev_step_antenna_in_rad_(0.01*degtorad),
     azim_start_antenna_in_rad_(-2.0*degtorad),
     elev_start_antenna_in_rad_(-6.0*degtorad),
     compute_beam_pattern_(!read_pattern),
     read_beam_pattern_(read_pattern),
     beam_parameter_set_(false),
     beam_widths_computed_(false),
     beam_num_set_(true),
     boresight_set_(false),
     pattern_loaded_(false),
     beam_number_(beamnum),
     beam_width_azim_oneway_("beam_width_azim_oneway_"),
     beam_width_elev_oneway_("beam_width_elev_oneway_"),
     beam_width_azim_twoway_("beam_width_azim_twoway_"),
     beam_width_elev_twoway_("beam_width_elev_twoway_"),
     beam_gain_("beam_gain_"), 
     beam_max_gain_(0.0),
     pattern_filename_("none")
 {
  beam_gain_.resize(Nazim_,Nelev_);
  config(cfg,beamnum,read_pattern);
  setBoresight(t);
 }

// get boresight
DirectionVector Beam::getBoresight() const
  { 
  if (!boresight_set_)
    {
    ErrorMessage e("Boresight not set.");
    e.throwMe();
    }
  return(boresight_);
  }

// Set Boresight 
DirectionVector Beam::setBoresight(const Time& t, const Uvar& terr) 
  throw(ErrorMessage)
{ 
 if (beam_num_set_ == false)
   {
     throw ErrorMessage("Beam number has not been set");
    
   }

 if( beam_parameter_set_== false)
   {
     throw ErrorMessage("Beam not configured, no frame available");
   }
 // terr is a time error in the rotation of the beam frame with respect to
 DirectionVector x("x",nominal_frame_,t+Time(terr),1,0,0);
 DirectionVector y("y",nominal_frame_,t+Time(terr),0,1,0);
 DirectionVector z("z",nominal_frame_,t+Time(terr),0,0,1);
 PositionVector origin("origin",nominal_frame_,t,0,0,0);
 Frame inertial_frame("J2000","Earth");
 x.representIn(inertial_frame);
 y.representIn(inertial_frame);
 z.representIn(inertial_frame);
 origin.representIn(inertial_frame);
 x.setTime(t);
 y.setTime(t);
 z.setTime(t);

 frame_=Frame(origin,x,y,z);
 boresight_=DirectionVector("boresight_",frame_,t,0,0,1);
 t_=t;
 boresight_set_=true;
 return(boresight_);
}

// Set Boresight 
DirectionVector Beam::setBoresight(const Time& t) throw(ErrorMessage)
{ 
 if (beam_num_set_ == false)
   {
     throw ErrorMessage("Beam number has not been set");
    
   }

 if( beam_parameter_set_== false)
   {
     throw ErrorMessage("Beam not configured, no frame available");
   }
 frame_=nominal_frame_;
 boresight_=DirectionVector("boresight_",frame_,t,0,0,1);
 t_=t;
 boresight_set_=true;
 return(boresight_);
}

// Beam::GetExtrema
// Returns vectors for various points in the gain pattern
// CENTER = boresight
// TOP = + 1/2 3dB (2-way) elevation width
// BOTTOM = - 1/2 3dB (2-way) elevation width
// LEFT = - 1/2 3dB (2-way) azimuth width
// RIGHT = + 1/2 3dB (2-way) azimuth width

DirectionVector Beam::getExtrema(BeamExtremaType etype)
{
  if(!boresight_set_){
    throw ErrorMessage("Boresight not set.");
  }
  DirectionVector tmp=boresight_; // sets time and frame of look vector
  // 1.5 dB one-way Gain = 3 dB two-way gain
  Uvec azim=elevAndGainToAzim(0,-1.5);
  Uvec elev=azimAndGainToElev(0,-1.5);
  
  if(azim.size()!=2 || elev.size()!=2){
    throw ErrorMessage("Beam extrema undefined");
  }

  
  switch (etype){
  case CENTER:
    tmp.setAzimuthElevation(0,0);
    break;
  case TOP:
    tmp.setAzimuthElevation(0,elev(1));
    break;
  case BOTTOM:
    tmp.setAzimuthElevation(0,elev(0));
    break;
  case LEFT:
    tmp.setAzimuthElevation(azim(0),0);
    break;
  case RIGHT:
    tmp.setAzimuthElevation(azim(1),0);
    break;
  default:
    throw ErrorMessage("Beam2::getExtrema: Bad input enum_type value");
  }
  return(tmp);
}

//------------------------------------------
//Get beam frame
//------------------------------------------
const Frame& Beam::getNominalFrame() const throw(ErrorMessage)
  { 

  if( beam_parameter_set_== false)
  {
    throw ErrorMessage("Beam not configured, no frame available");
  }

    return( nominal_frame_);
  }

//-----------------------------------------------------------
// beamFractionToDirection
//-----------------------------------------------------------

DirectionVector Beam::beamFractionToDirection(Uvar azim_frac, Uvar elev_frac)
{
  if(! boresight_set_){
    throw ErrorMessage("Boresight not set.");
  }
  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }
  DirectionVector tmp=boresight_;
  Uvar azim, elev;
  azim=azim_frac*beam_width_azim_twoway_;
  elev=elev_frac*beam_width_elev_twoway_;
  tmp.setAzimuthElevation(azim,elev);
  return(tmp);
}

//------------------------
// directionToBeamFraction
//------------------------

Uvec Beam::directionToBeamFraction(DirectionVector v)
{
  if(! boresight_set_){
    throw ErrorMessage("Boresight not set.");
  }
  if (!beam_widths_computed_)
    {
    estimateBeamWidthsFromPattern();
    }
  Uvec beam_frac("beam_frac",2);
  Uvar azim, elev;
  v.getAzimuthElevation(azim,elev);
  beam_frac(0)=azim/beam_width_azim_twoway_;
  beam_frac(1)=elev/beam_width_elev_twoway_;
  return(beam_frac);
}

