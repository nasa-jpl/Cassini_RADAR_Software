//----------------------------------------------------------------------------
// Frame.cpp
//
// This file contains method definitions for the Frame handling classes
// These classes use the NAIF Spice toolkit to provide automatic handling
// of position and state vectors and the coordinate transformations between
// them.
//
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_frame_c[] =
  "@(#) $Id: Frame.cpp,v 11.9 2017/03/13 22:22:25 richw Exp $";

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>
// Below is just for ifrnum_
#include <SpiceZfc.h>

//---------------
// Other includes
//---------------

#include <strings.h>
#include <iostream>
#include <string>
#include <sstream>
#include "Frame.h"
#include "Time.h"
#include "Constants.h"
#include "Profiling.h"
#include "Interpolate.h"
#include "TargetGeom.h"
#include "Flyby.h"
#include "SARProcParams.h"
#include "DebugInfo.h"

using std::cout;
using std::endl;
using std::cerr;

//---------------------
// DEFINITIONS
//--------------------

#define SPICEMSGLEN 1841
#define SPICEMODNAMELEN 11
#define ROTMAT_SIZE 45
#define QUAT_SIZE 20
#define INTERPOLATE_WINDOW_SIZE 5

//-----------------------------
// Global variables set by Frame::config
//-----------------------------
SpiceInt cassini_spice_id=0;
SpiceInt titan_spice_id=0;
SpiceInt default_target_spice_id=0;
SpiceInt sun_spice_id;
SpiceInt earth_spice_id;
SpiceInt iau_titan_frame_spice_id=0;
SpiceInt j2000_frame_spice_id=0;
SpiceInt beam_frame_spice_id[5]={0,0,0,0,0};
SpiceInt default_target_frame_spice_id=0;
PositionVector default_target_radii;
string default_target_name;
bool default_ring_target = false;

//--------------------------------
// Static Methods for class Frame
//--------------------------------

//----------
// selfTest
//----------

bool Frame::selfTest()
  {
  cout << endl << "No specific tests setup for Frame" << endl;
  /**
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/000728_SK_SOI_T45_82.bsp");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas00049.tsc");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/naif0007.tls");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas_radar_v11.ti");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/cas_v31.tf");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/pck00006.tpc");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/981005_PLTEPH-DE405S.bsp");
  spiceLoad("/u/svejk0/richw/dat/cassini/naif/sat083.bsp");
  spiceLoad("/u/svejk0/richw/dat/cassini/tour_sim/naif/pdt_c_kernel_p8e.bc");

  Time t("2005-301T04:30:37.840");
  Frame ftitan("IAU_TITAN",cassini_str);
  PositionVector r("r",ftitan,t,
    Uvar(0,km_str),Uvar(3575,km_str),Uvar(0,km_str));
  FloatVector v("v",ftitan,t,
    Uvar(6,km_per_s_str),Uvar(0,km_per_s_str),Uvar(0,km_per_s_str));
  DirectionVector x("x",ftitan,t,1,0,0);
  DirectionVector z("z",-r);
  DirectionVector y("y",cross(z,x));
  Frame fbeam("Beam3",r,x,y,z);
  DirectionVector bore("bore",fbeam,t,0,0,1);
  DirectionVector b2("b2",bore);
  b2.representIn(ftitan);
  b2.representIn(fbeam);
  if (b2 != bore) return(false);
  StateVector s1(r,v);
  Uvar a(2.0);
  StateVector s2 = s1;
  s2 *= a;
  s2 = a*s1;
  StateVector s3 = s2 - s1;
  */

  return(true);
  }

//-------
// Setup
//-------

//----------------------------------------------------------
// spiceLoad(filename)
//
// This static method loads data from a Spice kernel file.
// It must be called to load all the relevant kernel files
// before any Frame or geometry work can be done.
//----------------------------------------------------------

void Frame::spiceLoad(const string& filename)
  // No Exceptions
  {
  furnsh_c(filename.c_str());
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  }


//----------------------------------------------------------
// spiceUnLoad(filename)
//
// This static method unloads data from a Spice kernel file.
// 
//----------------------------------------------------------

void Frame::spiceUnLoad(const string& filename)
  // No Exceptions
  {
  unload_c(filename.c_str());
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  }


//-------------------------------------------------------
// config(cfg)
//
// Perform all necessary spiceloads using config object
//-------------------------------------------------------

void Frame::config(Config & cfg, bool need_ckernel)
  {
  //----------------------------
  // set up spice error reporting
  //----------------------------

  erract_c("SET",0,"RETURN");
  errdev_c("SET",0,"NULL");
  reset_c();

  //------------------------------
  // Load configuration parameters
  //------------------------------

  string prefix;
  if (cfg.keywordExists("spice_directory"))
    {
    prefix = cfg.str("spice_directory");
    }
  string tour_ephemeris = prefix + cfg.str("tour_ephemeris");
  string sclk_kernel = prefix + cfg.str("sclk_kernel");
  string time_kernel = prefix + cfg.str("time_kernel");
  string instrument_kernel = prefix + cfg.str("instrument_kernel");
  string frame_kernel = prefix + cfg.str("frame_kernel");
  string planetary_constants_kernel =
    prefix + cfg.str("planetary_constants_kernel");
  string solarsystem_ephemeris_kernel =
    prefix + cfg.str("solarsystem_ephemeris_kernel");
  string planet_ephemeris_kernel = prefix + cfg.str("planet_ephemeris_kernel");
  
  //-------------------------------------------
  // Load spice kernel files (except ckernels)
  //-------------------------------------------

  spiceLoad(tour_ephemeris);
  spiceLoad(sclk_kernel);
  spiceLoad(time_kernel);
  spiceLoad(instrument_kernel);
  spiceLoad(frame_kernel);
  spiceLoad(planetary_constants_kernel);
  spiceLoad(solarsystem_ephemeris_kernel);
  spiceLoad(planet_ephemeris_kernel);

  //-----------------------
  // Handle ckernel files
  //-----------------------

  loadCkernelSet(cfg,need_ckernel);

  // compute Frame specific commonly used frame and origin IDs
  target_=cfg.str("target");
  if(target_ =="NONE" || target_ =="None" || target_ =="none" ||target_=="Source" || target_=="source" || target_=="SOURCE")
    {
      //source : microwave source 
    }
  else
    {
      target_frame_="IAU_" + target_;
      spice_target_id(target_,target_origin_id_);
      spice_frame_id(target_frame_,target_frame_id_);
    }
  spice_target_id(cassini_str,cassini_origin_id_);
  spice_frame_id("CASSINI_RADAR_1",beam_frame_ids_[0]);
  spice_frame_id("CASSINI_RADAR_2",beam_frame_ids_[1]);
  spice_frame_id("CASSINI_RADAR_3",beam_frame_ids_[2]);
  spice_frame_id("CASSINI_RADAR_4",beam_frame_ids_[3]);
  spice_frame_id("CASSINI_RADAR_5",beam_frame_ids_[4]);

  // set up globals for commonly used frame and origin IDs and target radii
  spice_target_id("Cassini",cassini_spice_id);
  spice_target_id("Titan",titan_spice_id);
  spice_target_id("Sun",sun_spice_id);
  spice_target_id("Earth",earth_spice_id);
  spice_frame_id("IAU_TITAN",iau_titan_frame_spice_id);
  spice_frame_id("J2000",j2000_frame_spice_id);
  default_target_spice_id=target_origin_id_;
  default_target_frame_spice_id=target_frame_id_;
  default_target_name=target_;

  if (cfg.str("tracking_option") == "target_ring_radius_sweep_closest")
    {
    default_ring_target = true;
    }

  if(target_=="none" || target_=="NONE"||target_=="None" || target_=="Source" || target_=="source" || target_=="SOURCE")
    {
      default_target_radii=PositionVector();
    }
  else
    {
      spice_target_radii(target_,Frame(target_frame_id_,target_origin_id_),
		     Uvar(0,"s"),target_origin_id_,default_target_radii);
    }
  for(int i=0;i<5;i++) beam_frame_spice_id[i]=beam_frame_ids_[i];
  
  //---------------------------------------------------
  // load geometry speed up configuration parameters
  //---------------------------------------------------

  if(geom_mode_==PRECOMPUTE)
    {
    target_geometry_filename_=cfg.str("intermediate_geometry_file");
    target_=cfg.str("target");
    target_frame_="IAU_" + target_;
    Flyby tmp(cfg);
    Uvar t1,t2;
    t1=tmp.startTime()+tmp.epochTime();
    t2=tmp.endTime()+tmp.epochTime();
    Uvar dt=cfg["geometry_time_step"];
    time_step_=dt.getInUnits(seconds_str);
    num_time_steps_=(int)floor((t2-t1).getInUnits(seconds_str)/time_step_)+
      INTERPOLATE_WINDOW_SIZE+1;
    start_time_ = t1.getInUnits(seconds_str)-time_step_ *
      (INTERPOLATE_WINDOW_SIZE+1)/2;
    allocateSpecialGeometryArrays();
    computeSpecialGeometryArraysFromSpice();
    }
  if(geom_mode_==READ_FILE)
    {
    target_geometry_filename_=cfg.str("intermediate_geometry_file");
    target_=cfg.str("target");
    target_frame_="IAU_" + target_;
    readGeometryFromFile();
    computeDerivativeArrays();
    }
  
  
  }

//-------------------------------------------------------
// loadCkernelSet(cfg)
//
// Load all the ckernel files listed in the config object.
//-------------------------------------------------------

void Frame::loadCkernelSet(Config& cfg, bool must_find)
  {
  // Support the use of a directory prefix for ckernel files
  string ck_prefix;
  if (cfg.keywordExists("ckernel_directory"))
    {
    ck_prefix = cfg.str("ckernel_directory");
    }

  if (cfg.keywordExists("ckernel"))
    {  // If ckernel is present, load it first
    string ckernel = ck_prefix + cfg.str("ckernel");
    spiceLoad(ckernel);
    loaded_ckernel_filenames_.push_back(ckernel);
    }

  // Try to find ckernel anyway if required (allows interactive mode)
  else if (must_find){
    string ckernel = ck_prefix + cfg.str("ckernel");
    spiceLoad(ckernel);
    loaded_ckernel_filenames_.push_back(ckernel);
  }

  if (cfg.keywordExists("ckernel_1"))
    {  // Load set of ckernels if present
    unsigned int i = 1;
    while (1)
      {
      if (cfg.keywordExists("ckernel_" + toStr(i)))
        {
        string ckernel = ck_prefix + cfg.str("ckernel_" + toStr(i));
        spiceLoad(ckernel);
        loaded_ckernel_filenames_.push_back(ckernel);
        }
      else
        {
        break;  // only load a sequential set of names
        }
      ++i;
      }
    }

  }

//-------------------------------------------------------
// unloadCkernelSet(cfg)
//
// Unload all the currently loaded ckernel files listed
// in the config object.
//-------------------------------------------------------

void Frame::unloadCkernelSet()
  {
  for (list<string>::const_iterator cur_ckernel =
         loaded_ckernel_filenames_.begin();
       cur_ckernel != loaded_ckernel_filenames_.end();
       ++cur_ckernel)
    {
    spiceUnLoad(*cur_ckernel);
    }
  loaded_ckernel_filenames_.clear();
  }
    
//----------------------------------------------------------
// getCkernelList(ckernel_list)
//
// Fill the supplied list with the loaded ckernel filenames.
//----------------------------------------------------------

void Frame::getCkernelList(list<string>& loaded_ckernel_filenames)
  {
  loaded_ckernel_filenames = loaded_ckernel_filenames_;
  }

//----------------------------------------------------------
// numCkernelsLoaded()
//
// Returns the number of ckernel files currently loaded from
// the config file.
//----------------------------------------------------------

unsigned int Frame::numCkernelsLoaded()
  {
  return(loaded_ckernel_filenames_.size());
  }

//-------------------------
// Static initializations
//-------------------------

list<string> Frame::loaded_ckernel_filenames_;
Frame::GeometryModeE Frame::geom_mode_=Frame::USE_SPICE;
string Frame::target_geometry_filename_;
string Frame::target_frame_;
string Frame::target_;
double* Frame::quat_=NULL; 
double*  Frame::sc_vel_;   
double*  Frame::sc_pos_;   
double*  Frame::sc_acc_;   
double*  Frame::sc_jerk_;   
double*  Frame::d2qdt2_;   
double Frame::current_time_=-3e13; // initialized to 1000000 BC
bool Frame::rotmat_available_[5];
bool Frame::sc_pos_available_;
bool Frame::sc_vel_available_;
double Frame::start_time_;
double Frame::time_step_;
int Frame::num_time_steps_;
double* Frame::rotmat_interp_;
double* Frame::sc_vel_interp_;
double* Frame::sc_pos_interp_;

SpiceInt Frame::next_frame_id_=1;  // 0 reserved for invalid_frame

SpiceInt Frame::target_frame_id_;
SpiceInt Frame::target_origin_id_;
SpiceInt Frame::cassini_origin_id_;
SpiceInt Frame::beam_frame_ids_[5];

//--------------------------------
// Object Methods for class Frame
//--------------------------------

//--------------
// Constructors
//--------------

//-----------------------------------------------------------------------
// Frame()
//
// Default constructor makes an illegal frame which will cause an error
// exception if it is actually used.  It is intended for situations where
// a Frame needs to be built before it is loaded with valid data.
// Loading is done by assignment.
//-----------------------------------------------------------------------

Frame::Frame()
  : spice_frame_(true), spice_frame_id_(0), frame_id_(0),
    spice_origin_id_(0),  
    special_frame_number_(-1), beam_number_(0)
  {  }

Frame::Frame(const Frame& f)
 : spice_frame_(f.spice_frame_), 
    spice_frame_id_(f.spice_frame_id_), 
    frame_id_(f.frame_id_),
    spice_origin_id_(f.spice_origin_id_),  
    special_frame_number_(f.special_frame_number_), 
    beam_number_(f.beam_number_)
{
  if(!spice_frame_){
    for(int i=0;i<3;i++){
      rel_origin_[i]=f.rel_origin_[i];
      for(int j=0;j<3;j++){
	spice_xform_[i][j]=f.spice_xform_[i][j];
	spice_invxform_[i][j]=f.spice_invxform_[i][j];
      }
    }
  }
}

// Simplest Spice frame constructor
Frame::Frame(int frame_id, int origin_id)
  : spice_frame_(true), 
    spice_frame_id_(frame_id), frame_id_(frame_id),
    spice_origin_id_(origin_id)  
{ 
  setSpecialFrameNumber();
}

//-----------------------------------------
// Streamlined Non-Spice frame constructors
//-----------------------------------------

//-----------------------------------------------------------------------------
// Frame(frame,rel_origin)
//
// This constructor sets up a new Frame with the same
// orientation as the input frame, and with its origin translated to the
// indicated position.
//-----------------------------------------------------------------------------

Frame::Frame(const Frame& frame,
  PositionVector rel_origin)
  : spice_frame_(false), 
    spice_frame_id_(frame.spice_frame_id_), 
    spice_origin_id_(frame.spice_origin_id_),  
    special_frame_number_(frame.special_frame_number_), 
    beam_number_(frame.beam_number_)
  {
  spice_xform_[0][0] = frame.spice_xform_[0][0];
  spice_xform_[0][1] = frame.spice_xform_[0][1];
  spice_xform_[0][2] = frame.spice_xform_[0][2];
  spice_xform_[1][0] = frame.spice_xform_[1][0];
  spice_xform_[1][1] = frame.spice_xform_[1][1];
  spice_xform_[1][2] = frame.spice_xform_[1][2];
  spice_xform_[2][0] = frame.spice_xform_[2][0];
  spice_xform_[2][1] = frame.spice_xform_[2][1];
  spice_xform_[2][2] = frame.spice_xform_[2][2];

  spice_invxform_[0][0] = frame.spice_invxform_[0][0];
  spice_invxform_[0][1] = frame.spice_invxform_[0][1];
  spice_invxform_[0][2] = frame.spice_invxform_[0][2];
  spice_invxform_[1][0] = frame.spice_invxform_[1][0];
  spice_invxform_[1][1] = frame.spice_invxform_[1][1];
  spice_invxform_[1][2] = frame.spice_invxform_[1][2];
  spice_invxform_[2][0] = frame.spice_invxform_[2][0];
  spice_invxform_[2][1] = frame.spice_invxform_[2][1];
  spice_invxform_[2][2] = frame.spice_invxform_[2][2];


  setSpiceRelativeOrigin(rel_origin);
  getNewFrameID();
  }

//-----------------------------------------------------------------------------
// Frame(origin,x,y,z)
//
// Construct a new frame.
// The new frame has an origin relative to some defining frame,
// and cardinal axes x,y,z relative to some defining frames.
// The new frame will be expressed relative to the frame used to define the
// origin.
// Note that the new frame is fixed relative to the defining frame because
// no data is available on the relative motion of the origin or axes.
// An exception is thrown if the times of the input vectors don't match.
//-----------------------------------------------------------------------------

Frame::Frame(const PositionVector& origin,
  DirectionVector x, DirectionVector y, DirectionVector z)
  : spice_frame_(false), spice_frame_id_(origin.frame_.spice_frame_id_), 
    spice_origin_id_(origin.frame_.spice_origin_id_),  
    special_frame_number_(origin.frame_.special_frame_number_), 
    beam_number_(origin.frame_.beam_number_)
  {
  if (origin.timeInSeconds() != x.timeInSeconds() || origin.timeInSeconds() 
      != y.timeInSeconds() ||
      origin.timeInSeconds() != z.timeInSeconds())
    {
      ErrorMessage e("Time mismatch while constructing new frame");
      e.throwMe();
    }

  setSpiceRelativeRotation(x,y,z,origin.frame_);
  setSpiceRelativeOrigin(origin);
  getNewFrameID();
  }

//-----------------------------------------------------------------------------
// Frame(x,y,z)
//
// Construct a new frame.
// The new frame has the same origin as direction vector x.
// As a result, this frame is suitable for the construction of rotated frame
// with its origin unchanged w.r.t. the origin of x.
// If the frame of x is in spice frame, this new frame will be recogized 
// by spice routines
// Note: Although this frame is constructed by simply  rotating a spice frame
// , it is not any more spice frame.
//-----------------------------------------------------------------------------

Frame::Frame(DirectionVector x,
	     DirectionVector y, DirectionVector z)
  : spice_frame_(false), spice_frame_id_(x.frame_.spice_frame_id_), 
    spice_origin_id_(x.frame_.spice_origin_id_),  
    special_frame_number_(x.frame_.special_frame_number_), 
    beam_number_(x.frame_.beam_number_)
  {
  if (x.timeInSeconds() != y.timeInSeconds() 
      || y.timeInSeconds() != z.timeInSeconds())
    {
      ErrorMessage e("Time mismatch while constructing new frame:");
    e.throwMe();
    }

  setSpiceRelativeRotation(x,y,z,x.frame_);
  rel_origin_[0]=x.frame_.rel_origin_[0];
  rel_origin_[1]=x.frame_.rel_origin_[1];
  rel_origin_[2]=x.frame_.rel_origin_[2];
  getNewFrameID();
  }


//----------------------------------------------------------------------------
// Frame(frame_name,origin_name)
//
// This constructor sets up a Spice frame with the origin at the indicated
// Spice object.  Unrecognized object names cause an error exception.
// Unrecognized frame names will cause a Spice error abort.
//-----------------------------------------------------------------------------

Frame::Frame(const string& frame_name, const string& origin_name)
  : spice_frame_(true) 
  {
    spice_target_id(origin_name,spice_origin_id_);
    spice_frame_id(frame_name,spice_frame_id_);
    *this=Frame(spice_frame_id_,spice_origin_id_);
  }

//-----------------------------------------------------------------------------
// Frame(frame_name,frame,rel_origin)
//
// This constructor sets up a new Frame called frame_name with the same
// orientation as the input frame, and with its origin translated to the
// indicated position.
// frame_name is ignored!!!
//-----------------------------------------------------------------------------

Frame::Frame(const string& frame_name, const Frame& frame,
  PositionVector rel_origin)
  {
    *this=Frame(frame,rel_origin);
  }

//-----------------------------------------------------------------------------
// Frame(frame_name,origin,x,y,z)
//
// Construct a new frame with the indicated name. (Name is ignored)
// The new frame has an origin relative to some defining frame,
// and cardinal axes x,y,z relative to some defining frames.
// The new frame will be expressed relative to the frame used to define the
// origin.
// Note that the new frame is fixed relative to the defining frame because
// no data is available on the relative motion of the origin or axes.
// An exception is thrown if the times of the input vectors don't match.
//-----------------------------------------------------------------------------

Frame::Frame(const string& frame_name, const PositionVector& origin,
  DirectionVector x, DirectionVector y, DirectionVector z)
{
  *this=Frame(origin,x,y,z);
}

//-----------------------------------------------------------------------------
// Frame(frame_name,x,y,z)
//
// Construct a new frame with the indicated name. (Ignores name.)
// The new frame has the same origin as direction vector x.
// As a result, this frame is suitable for the construction of rotated frame
// with its origin unchanged w.r.t. the origin of x.
// If the frame of x is in spice frame, this new frame will be recogized 
// by spice routines
// Note: Although this frame is constructed by simply  rotating a spice frame
// , it is not any more spice frame.
//-----------------------------------------------------------------------------

Frame::Frame(const string& frame_name, DirectionVector x,
	     DirectionVector y, DirectionVector z)
{
  *this=Frame(x,y,z);
}

//-----------------------------------------------------------------------------
// Frame(frame_name,frame,rot)
//
// Construct a new frame with the indicated name relative to the indicated
// defining frame. (Ignores name.)
// The new frame has the same origin as the defining frame,
// and orientation specified by the ordered rotation in rot.
// Note that the new frame is fixed relative to the defining frame because
// no data is available on the relative motion of the origin or axes.
//-----------------------------------------------------------------------------

Frame::Frame(const string& frame_name, const Frame& frame, const Rotation& rot)  : spice_frame_(false), spice_frame_id_(frame.spice_frame_id_), 
    spice_origin_id_(frame.spice_origin_id_),  
    special_frame_number_(frame.special_frame_number_), 
    beam_number_(frame.beam_number_)
  {
  setSpiceRelativeRotation(rot,frame);
  rel_origin_[0] = frame.rel_origin_[0];
  rel_origin_[1] = frame.rel_origin_[1];
  rel_origin_[2] = frame.rel_origin_[2];
  getNewFrameID();
  }


//----------------------
// Geometry speedup methods
//-----------------------


void Frame::allocateSpecialGeometryArrays(){
  cleanUp();
  quat_=(double*)malloc(sizeof(double)*QUAT_SIZE*num_time_steps_);
  d2qdt2_=(double*)malloc(sizeof(double)*QUAT_SIZE*num_time_steps_);
  rotmat_interp_=(double*)malloc(sizeof(double)*ROTMAT_SIZE);
  sc_pos_=(double*)malloc(sizeof(double)*3*num_time_steps_);
  sc_vel_=(double*)malloc(sizeof(double)*3*num_time_steps_);
  sc_acc_=(double*)malloc(sizeof(double)*3*num_time_steps_);
  sc_jerk_=(double*)malloc(sizeof(double)*3*num_time_steps_);
  sc_pos_interp_=(double*)malloc(sizeof(double)*3);
  sc_vel_interp_=(double*)malloc(sizeof(double)*3);  
  current_time_=-3e13; // initialized to 1000000 BC
  for(int i=0;i<5;i++) rotmat_available_[i]=false;
  sc_pos_available_=false;
  sc_vel_available_=false;
}

// This must be called prior to Config or else USE_SPICE will be assumed
void Frame::setGeometryMode(Frame::GeometryModeE gm){
  geom_mode_=gm;
}

void Frame::writeGeometryToFile(){
  FILE* fp=fopen(target_geometry_filename_.c_str(),"w");
  if(fp==NULL){
    ErrorMessage e("Cannot create geometry file " + 
		   target_geometry_filename_);
    e.throwMe();
  }
  if(fwrite((void*) &start_time_, sizeof(double),1 ,fp) != 1 ||
     fwrite((void*) &time_step_, sizeof(double),1 ,fp) != 1 ||
     fwrite((void*) &num_time_steps_, sizeof(int),1 ,fp) != 1 ||
     fwrite((void*) quat_, sizeof(double), num_time_steps_*QUAT_SIZE, fp)
     !=  (unsigned) num_time_steps_*QUAT_SIZE  || 
     fwrite((void*) sc_pos_, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3 ||
     fwrite((void*) sc_vel_, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3){
    
    ErrorMessage e("Error writing geometry file " 
		   + target_geometry_filename_);
    e.throwMe();
  }
  fclose(fp);
}
void Frame::compareGeometryArraysToFile(const string& file){
  FILE* fp=fopen(file.c_str(),"r");
  if(fp==NULL){
    ErrorMessage e("Cannot open comparison geometry file " + 
		   file);
    e.throwMe();
  }
  double start_time, time_step;
  int num_time_steps;
  if(fread((void*) &start_time, sizeof(double),1 ,fp) != 1 ||
     fread((void*) &time_step, sizeof(double),1 ,fp) != 1 ||
     fread((void*) &num_time_steps, sizeof(int),1 ,fp) != 1){
    ErrorMessage e("Error reading geometry file " 
		   + target_geometry_filename_);
    e.throwMe();
  }

    
  cout.precision(10);

  cout << " Comparing new geometry in file " << target_geometry_filename_
       << " with old geometry in file " << file << endl;

  cout << " Start time is (new,old,diff)=" 
       << "(" << start_time_ << " s," << start_time << " s,"
       << start_time_-start_time << " s)" << endl;
  cout << " Time step is (new,old,diff)=" 
       << "(" << time_step_ << " s," << time_step << " s,"
       << time_step_-time_step << " s)" << endl;
  cout << " Number of time steps is (new,old,diff)=" 
       << "(" << num_time_steps_ << "," << num_time_steps << ","
       << num_time_steps_-num_time_steps << ")" << endl;


  // check to make sure arrays in file are the same size as arrays in memory
  // if not report error and return normally
  if(num_time_steps!=num_time_steps_){
    cout << "Comparison geometry file has wrong number of time steps.";
    fclose(fp);
    return;
  }
  double* quat=(double*)malloc(sizeof(double)*QUAT_SIZE*num_time_steps_);
  double* sc_pos=(double*)malloc(sizeof(double)*3*num_time_steps_);
  double* sc_vel=(double*)malloc(sizeof(double)*3*num_time_steps_);

  if(fread((void*) quat, sizeof(double), num_time_steps_*QUAT_SIZE, fp)
     !=  (unsigned) num_time_steps_*QUAT_SIZE  || 
     fread((void*) sc_pos, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3 ||
     fread((void*) sc_vel, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3){
    
    ErrorMessage e("Error reading comparison geometry file " 
		   + file);
    e.throwMe();
  }
  fclose(fp);

  double max_quat_err[QUAT_SIZE];
  double max_pos_err[3]={0,0,0};
  double max_vel_err[3]={0,0,0};
  double abs_max_quat_err=0;
  for(int i=0;i<QUAT_SIZE;i++) max_quat_err[i]=0;



 
  for(int n=0;n<num_time_steps_;n++){
    int offset1=n*QUAT_SIZE;
    int offset2=3*n;
    for(int i=0;i<QUAT_SIZE;i++){
      double err = fabs(quat[offset1+i]-quat_[offset1+i]);
      if(err>max_quat_err[i])max_quat_err[i]=err;
      if(err>abs_max_quat_err)abs_max_quat_err=err;
    }			 
    for(int i=0;i<3;i++){
      double err = fabs(sc_pos[offset2+i]-sc_pos_[offset2+i]);
      if(err>max_pos_err[i])max_pos_err[i]=err;
      err = fabs(sc_vel[offset2+i]-sc_vel_[offset2+i]);
      if(err>max_vel_err[i])max_vel_err[i]=err;
    }
  }

  for(int b=0;b<5;b++){
    cout <<"   BEAM " << b+1 <<" Maximum Quaternion difference:" << endl;
    int start_idx= 4*b;
    for(int j=0;j<4;j++){
	    cout << max_quat_err[start_idx+j] << "\t";
    }
    cout << endl;
    
  }
  cout << " Absolute Maximum Quaternion Component Difference is "
       << abs_max_quat_err << endl;
  cout << " Maximum S/C position error is (x,y,z) = " 
       << "(" << max_pos_err[0] <<" km," << max_pos_err[1] << "km, "
       << max_pos_err[2] << "km)" << endl;
  cout << " Maximum S/C velocity error is (x,y,z) = " 
       << "(" << max_vel_err[0] <<" km/s," << max_vel_err[1] << "km/s, "
       << max_vel_err[2] << "km/s)" << endl;
  free(quat);
  free(sc_pos);
  free(sc_vel);
}
void Frame::readGeometryFromFile(){
  FILE* fp=fopen(target_geometry_filename_.c_str(),"r");
  if(fp==NULL){
    ErrorMessage e("Cannot open geometry file " + 
		   target_geometry_filename_);
    e.throwMe();
  }
  if(fread((void*) &start_time_, sizeof(double),1 ,fp) != 1 ||
     fread((void*) &time_step_, sizeof(double),1 ,fp) != 1 ||
     fread((void*) &num_time_steps_, sizeof(int),1 ,fp) != 1){
    ErrorMessage e("Error reading geometry file " 
		   + target_geometry_filename_);
    e.throwMe();
  }
  allocateSpecialGeometryArrays();
  if(fread((void*) quat_, sizeof(double), num_time_steps_*QUAT_SIZE, fp)
     !=  (unsigned) num_time_steps_*QUAT_SIZE  || 
     fread((void*) sc_pos_, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3 ||
     fread((void*) sc_vel_, sizeof(double), num_time_steps_*3, fp) 
     != (unsigned) num_time_steps_*3){
    
    ErrorMessage e("Error reading geometry file " 
		   + target_geometry_filename_);
    e.throwMe();
  }
  fclose(fp);
}


void Frame::cleanUp(){
  if(quat_!=NULL){
    free(quat_);
    free(d2qdt2_);
    free(rotmat_interp_);
    free(sc_pos_);
    free(sc_pos_interp_);
    free(sc_vel_);
    free(sc_vel_interp_);
    free(sc_acc_);
    free(sc_jerk_);
    quat_=NULL;
  }
}
void Frame::computeDerivativeArrays(){
  // setup time array (x)
  // and 1-D time varying quantity array (y)
  // and intermed 2nd derivative d2y
  double* y = (double*) malloc(sizeof(double)*num_time_steps_);
  double* d2y = (double*) malloc(sizeof(double)*num_time_steps_);
  double* x = (double*) malloc(sizeof(double)*num_time_steps_);
  for(int i=0;i<num_time_steps_;i++){
    x[i]=start_time_+time_step_*i;
  }

  
  // compute quaternion second derivative

  for(int i=0;i<QUAT_SIZE;i++){
    
    // copy 1-D array to y
    for(int j=0;j<num_time_steps_;j++){
      y[j]=quat_[i+QUAT_SIZE*j];
    }

    // compute second derivative
    cubic_spline(x,y,num_time_steps_,1e33,1e33,d2y);

    // copy second derivative to storage array
    for(int j=0;j<num_time_steps_;j++){
      d2qdt2_[i+QUAT_SIZE*j]=d2y[j];
    }
  }
  

  // compute s/c acceleration
  for(int i=0;i<3;i++){
    
    // copy 1-D array to y
    for(int j=0;j<num_time_steps_;j++){
      y[j]=sc_pos_[i+3*j];
    }

    // compute second derivative
    double dy0=sc_vel_[i];
    double dyn=sc_vel_[i+(num_time_steps_-1)*3];
    cubic_spline(x,y,num_time_steps_,dy0,dyn,d2y);

    // copy second derivative to storage array
    for(int j=0;j<num_time_steps_;j++){
      sc_acc_[i+3*j]=d2y[j];
    }
  }

  // compute s/c jerk
  for(int i=0;i<3;i++){
    
    // copy 1-D array to y
    for(int j=0;j<num_time_steps_;j++){
      y[j]=sc_vel_[i+3*j];
    }

    // compute second derivative
    double dy0=sc_acc_[i];
    double dyn=sc_acc_[i+(num_time_steps_-1)*3];
    cubic_spline(x,y,num_time_steps_,dy0,dyn,d2y);

    // copy second derivative to storage array
    for(int j=0;j<num_time_steps_;j++){
      sc_jerk_[i+3*j]=d2y[j];
    }
  }
  // free intermediate arrays
  free(y);
  free(x);
  free(d2y);
}
void Frame::interpolateRotMat(int beam_number){
  if(!rotmat_available_[beam_number-1]){
    if(beam_number <1 || beam_number > 5){
      ErrorMessage e("Frame::interpolateRotMat beam_number outside range");
      e.throwMe();
    }
    int b=beam_number-1;
    int rmstart_idx=9*b;
    int qstart_idx=4*b;
    int qskip=QUAT_SIZE;
    double q[4];

    // interpolate quaternion



    
    double normq=0;
    for(int i=0;i<4;i++){
      q[i]=interpolate1D(quat_,d2qdt2_,qstart_idx+i,qskip,current_time_);
      normq+=q[i]*q[i];
    }

    normq=sqrt(normq);

    // normalize quaternion and flip sign if necessary

   
    bool flip_sign = (q[0]<0);
    for(int i=0;i<4;i++){
      q[i]/=normq;
      if(flip_sign) q[i]=-q[i];
    }

   
    /*** commented out because cubic spline was found to be more accurate 
    // spherical linear interpolation
    // q_interp= ((q2*(-q1)))^alpha*q1
    // normalizaton and sign flip (if necessary) are in inter...Quat..()
    interpolateQuaternion(quat_,qstart_idx,qskip,current_time_,q);
    ***/
    
    // convert quaternion to rotation matrix
    // (don't use spice)

    double* m=&(rotmat_interp_[rmstart_idx]);
    double w=q[0], x=q[1], y=q[2], z=q[3];
    m[0]=1 - 2*y*y - 2*z*z;
    m[1]=2*x*y - 2*w*z;
    m[2]=2*x*z + 2*w*y;

    m[3]=2*x*y + 2*w*z;
    m[4]=1 - 2*x*x - 2*z*z;
    m[5]=2*y*z - 2*w*x;

    m[6]=2*x*z - 2*w*y;
    m[7]=2*y*z + 2*w*x;
    m[8]=1 - 2*x*x - 2*y*y;

    rotmat_available_[b]=true;
  }
}

void Frame::interpolateQuaternion(double* qs, 
				    int start_idx, int skip, 
				    double t, double* q, int debug){

  // interp equation is q=((q2*q1^-1)^alpha)*q1;
  // compute bounding times
  if(t<=start_time_){
    ErrorMessage e("Frame: attempt to interpolate quaternion before start_time");
    e.throwMe();
  }

  int ti=(int)floor((t-start_time_)/time_step_);
  if(ti>=num_time_steps_-1){
    ErrorMessage e("Frame: attempt to interpolate quaternion after end_time");
    e.throwMe();
  }
  double t1=start_time_+ti*time_step_;

  // compute interpolation coefficient
  double alpha=(t-t1)/time_step_;

  if(debug){
    cout<< "InterpolateQuaternion debug ...." << endl;
    cout << "alpha=" << alpha << " t1(rel)= " << t1-start_time_  << " t(rel)= "
       << t-start_time_  << endl;
  }
  // get boundary quaternions
  int i=ti*QUAT_SIZE+start_idx;
  double* q1=&(qs[i]);
  double* q2=&(qs[i+skip]);

  if(q1[0]<0){
    q1[0]=-q1[0];
    q1[1]=-q1[1];
    q1[2]=-q1[2];
    q1[3]=-q1[3];
  }

  if(q2[0]<0){
    q2[0]=-q2[0];
    q2[1]=-q2[1];
    q2[2]=-q2[2];
    q2[3]=-q2[3];
  }
  if(debug){
  cout << " q1:";

  printQuaternion(q1);
  cout << endl;

  cout << " q2:";
  printQuaternion(q2);
  cout << endl;
  }
  // compute q1^-1 
  double q1i[4];
  q1i[0]=q1[0];
  q1i[1]=-q1[1];
  q1i[2]=-q1[2];
  q1i[3]=-q1[3];

  if(debug){
    cout << " q1^-1:";
    printQuaternion(q1i);
    cout << endl;
  }

  // multiply q2 by q1^-1 
  double qtrans[4];
  multiplyQuaternions(q2,q1i,qtrans);

  if(debug){
    cout << " (q2*q1^-1):";
    printQuaternion(qtrans);
    cout << endl;
  }

  // raise qtrans to the alpha power (flip sign if necessary to go short way)
  if(qtrans[0]<0){
    qtrans[0]=-qtrans[0];
    qtrans[1]=-qtrans[1];
    qtrans[2]=-qtrans[2];
    qtrans[3]=-qtrans[3];
  }
  raiseQuaternionToPower(qtrans,alpha);

  if(debug){
    cout << " (q2*q1^-1)^alpha:";
    printQuaternion(qtrans);
    cout << endl;
  }

  // perform final multiply
  multiplyQuaternions(qtrans,q1,q);

  if(debug){
    cout << " ((q2*q1^-1)^alpha)*q1):";
    printQuaternion(q);
    cout << endl;
  }

  // flip sign of [q] so that q[0] is positive
  if(q[0]<0){
    q[0]=-q[0];
    q[1]=-q[1];
    q[2]=-q[2];
    q[3]=-q[3];
  } 

 if(debug){
   cout << " Final output quaternion:";
   printQuaternion(q);
   cout << endl;
 }
}
void Frame::printQuaternion( double q[4]){
  double cosold=q[0];
  double sinold=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  double angold;
  if(cosold==0) angold=pi/2;
  else angold=atan(sinold/cosold);
  cout.precision(20);
  if(sinold==0) sinold=1;
  cout << "Rotation angle :" << angold*2*180/pi << "degrees, V=["
       << q[1]/sinold <<"," << q[2]/sinold <<"," << q[3]/sinold << "]";
}
void Frame::raiseQuaternionToPower(double q[4], double a){
  // zero rotation quat special case
  if(q[0]==1) return;

  double cosold=q[0];
  double sinold=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  double angold;
  if(cosold==0) angold=pi/2;
  else angold=atan(sinold/cosold);
  double cosnew=cos(angold*a);
  double sinnew=sin(angold*a);
  q[0]=cosnew;
  q[1]=q[1]*sinnew/sinold;
  q[2]=q[2]*sinnew/sinold;
  q[3]=q[3]*sinnew/sinold;   
}

  void Frame::multiplyQuaternions(double q1[4], double q2[4], double q[4]){

    // quaternion definition
    // q=[w,x,y,z]=[w,V]

    // quaternion multiply def
    //q = q1q2=[w1w2-dot(V1,V2),(w1V2+w2V1+cross(V1,V2))]

    double w1=q1[0], x1=q1[1], y1=q1[2], z1=q1[3];

    double w2=q2[0], x2=q2[1], y2=q2[2], z2=q2[3];

    q[0]=w1*w2-x1*x2-y1*y2-z1*z2;

    q[1]=w1*x2+w2*x1+y1*z2-z1*y2;
    q[2]=w1*y2+w2*y1+z1*x2-x1*z2;
    q[3]=w1*z2+w2*z1+x1*y2-y1*x2;

  }
  double Frame::interpolate1D(double* array, double* d2arraydt,
			    int start_idx, int skip, 
			      double t)
{
  double y;
  if(! interpolate_cubic_spline(start_time_, time_step_, array, d2arraydt, 
				num_time_steps_, t, &y, start_idx, skip)){
    ErrorMessage 
      e("Attempt to interpolate outside time range of geometry arrays");
    e.throwMe();

  }
  return(y);
}

void Frame::interpolateSCVel(){
  if(!sc_vel_available_){
    int skip=3;
    for(int i=0;i<3;i++){
      sc_vel_interp_[i]=interpolate1D(sc_vel_, sc_jerk_,i,skip,current_time_);
    }
    sc_vel_available_=true;
  }
}

void Frame::interpolateSCPos(){
  if(!sc_pos_available_){
    int skip=3;
    for(int i=0;i<3;i++){
      sc_pos_interp_[i]=interpolate1D(sc_pos_, sc_acc_,i,skip,current_time_);
    }
    sc_pos_available_=true;
  }
}

void Frame::updateTime(double t_in_s){
  if(t_in_s!=current_time_){
    current_time_=t_in_s;
    sc_pos_available_=false;
    sc_vel_available_=false;
    for(int i=0;i<5;i++){
      rotmat_available_[i]=false;
    }
  }
}


void Frame::updateTime(const Time& t){
  double t_in_s=get_in_base_units(t.et());
  updateTime(t_in_s);
}

void Frame::computeSpecialGeometryArraysFromSpice(){

  GeometryModeE old = geom_mode_;
  geom_mode_=USE_SPICE; //set this so that computations will use SPICE

  Frame tf(target_frame_id_,target_origin_id_);

    for(int i=0;i<num_time_steps_;i++){
      double t=start_time_+time_step_*i;

      // set temporary ephemeris array
      double vel[3], pos[3];

      // Compute spacecraft state in target body fixed coordinates
      // using spice
      StateVector s;
      tf.ephemeris(s, cassini_origin_id_, t);
      pos[X]=get_in_base_units(s.r_[X]);
      pos[Y]=get_in_base_units(s.r_[Y]);
      pos[Z]=get_in_base_units(s.r_[Z]);
      vel[X]=get_in_base_units(s.v_[X]);
      vel[Y]=get_in_base_units(s.v_[Y]);
      vel[Z]=get_in_base_units(s.v_[Z]);

      // copy to ephemeris arrays
      for (int j=0;j<3;j++){
	sc_vel_[i*3+j]=vel[j];
	sc_pos_[i*3+j]=pos[j];
      }
 
      for(int b=0;b<5;b++){
        // set up beam frame and temporary arrays;
	Frame bf(beam_frame_ids_[b],cassini_origin_id_);
	double b2t[3][3];

        // Compute Rotation Matrix for each beam
	rotationMatrix(bf,tf,t,b2t);

        // convert into quat (using SPICE);
	int start_idx= 4*b;        
        SpiceDouble q[4];
        m2q_c(b2t,q);


        //--------------------------------------------------------------
        // to keep the quaternions varying in a nice continuous way over
        // time we make sure that the xyz vector (q[1],q[2],q[3])
        // doesn't flip sign from one time sample to the next.
        // realize this allows the rotation to be "the long way around"
        // but this effect is corrected after the quaternion is interpolated
        //---------------------------------------------------------------


        double vmax=0.0;
        int max_j=1;

        // find maximal rotation axis component
        for(int j=0;j<4;j++){
           if(j>0 && fabs(q[j])>vmax) {
	    vmax=fabs(q[j]);
            max_j=j;          
	  }
	}

        // check to see if maximal component changed sign from last time
        // sample
	bool flipped=false;
	if(i!=0){
	  double testval=q[max_j]/quat_[(i-1)*QUAT_SIZE+start_idx+max_j];
	  if(testval<0) flipped=true;
	}

        // if so negate quaternion to preserve continuity
	if(flipped){
	  for(int j=0;j<4;j++){
	    q[j]=-q[j];
	  }
	}

	// copy quaternion into storage array
	for(int j=0;j<4;j++){
	  quat_[i*QUAT_SIZE+start_idx+j]=q[j];
	}
      }
    }

    // do pre-computation of derivative necessary for cublic spline interp
    computeDerivativeArrays();
    geom_mode_=old; // reset geometry mode to previous value (should be
                        // PRECOMPUTE)
}

// Routine rotates each quaternion in geometry array so that the
// boresight points to a particular range and doppler value as determined
// by the correction polynomial in the SARProcParams object
// a single minimum angle rotation is performed
// When debug level is 1 the angle and the axes of rotation are output
void Frame::adjustSpecialGeometryArrays(Config& cfg){
  SARProcParams spp(cfg,1);
  if(!spp.set_doppler_centroid) return;
  DebugInfo dbg("Frame::adjustSpecialGeometryArrays");
  GeometryModeE old = geom_mode_;
  geom_mode_=USE_SPICE; //set this so that computations will use SPICE
  Frame tf(target_frame_id_,target_origin_id_);
  // for now hack lambda in
  Uvar lambda=speed_light/Uvar(13.78e9,"Hz");
  for(int i=0;i<num_time_steps_;i++){
    double t=start_time_+time_step_*i;
    Time time;
    time.setEt(t);
    Frame bf(beam_frame_ids_[3],cassini_origin_id_);
    DirectionVector b3bore_old(bf,t,0,0,1);

    // compute old beam 3 boresight
    TargetGeom tg(time);
    tg.setTarget();
    b3bore_old.representIn(tf);
    StateVector sc_state;
    tf.ephemeris(sc_state,cassini_origin_id_,t);
    tg.setState(sc_state);
    tg.setLookDirection(b3bore_old);
    dbg.file << "A" << endl;
    Uvar range=tg.range();
    Uvar doppler=tg.doppler(lambda);
   
    // compute new beam 3 boresight
    TargetGeom tg2(time);
    tg2.setTarget();
    tg2.setState(sc_state);
    Uvar drange,ddopc;
    Uvar dummy1,dummy2,dummy3;
    spp.setDopplerProfile(time,3,dummy1,ddopc,drange,dummy2,dummy3);
    
     
    tg2.setRangeDopplerInTargetFrame(range+drange,doppler+ddopc,lambda);
    dbg.file << "B:" << range << " " << drange << " " << doppler << " "
	     << ddopc << " " <<lambda << endl;
    DirectionVector b3bore_new;
    try{
    b3bore_new=tg2.lookDirection();
    }
    catch(ErrorMessage e){
      continue;
    }
    dbg.file << "C" << endl;
    Uvar ang=b3bore_old.angle(b3bore_new);
    double a=ang.getInUnits("rads");
    DirectionVector axis=cross(b3bore_old,b3bore_new);
    double q1[4],q2[4],q[4];
    double cos_a=cos(a/2);
    double sin_a=sin(a/2);
    q2[0]= cos_a;
    q2[1]= sin_a*axis[DirectionVector::X];
    q2[2]= sin_a*axis[DirectionVector::Y];
    q2[3]= sin_a*axis[DirectionVector::Z];
    double norm=sqrt(q2[0]*q2[0]+q2[1]*q2[1]+q2[2]*q2[2]+q2[3]*q2[3]);
    for(int b=0;b<1;b++){
      int start_idx=4*b;
      for(int d=0;d<4;d++){
	q2[d]=q2[d]/norm;
	q1[d]=quat_[i*QUAT_SIZE+start_idx+d];
      }
      tf.multiplyQuaternions(q1,q2,q);

      //--------------------------------------------------------------
      // to keep the quaternions varying in a nice continuous way over
      // time we make sure that the xyz vector (q[1],q[2],q[3])
      // doesn't flip sign from one time sample to the next.
      // realize this allows the rotation to be "the long way around"
      // but this effect is corrected after the quaternion is interpolated
      //---------------------------------------------------------------


      double vmax=0.0;
      int max_j=1;

      // find maximal rotation axis component
      for(int j=0;j<4;j++){
	if(j>0 && fabs(q[j])>vmax) {
	  vmax=fabs(q[j]);
	  max_j=j;          
	}
      }

      // check to see if maximal component changed sign from last time
      // sample
      bool flipped=false;
      if(i!=0){
	double testval=q[max_j]/quat_[(i-1)*QUAT_SIZE+start_idx+max_j];
	if(testval<0) flipped=true;
      }

      // if so negate quaternion to preserve continuity
      if(flipped){
	for(int j=0;j<4;j++){
	  q[j]=-q[j];
	}
      }

      // copy quaternion into storage array
      for(int j=0;j<4;j++){
	quat_[i*QUAT_SIZE+start_idx+j]=q[j];
      } 

    } // end beam loop
    if(dbg.level){
      dbg.file << "t= " << t-start_time_ << " angle=" << a*180/pi << endl;
    }
  } // end time step loop

    // do pre-computation of derivative necessary for cublic spline interp
  computeDerivativeArrays();
  geom_mode_=old; // reset geometry mode to previous value (should be
  // PRECOMPUTE)
  dbg.file << "Z" << endl;
}

// Determines accuracy of interpolation of S/C position and velocity in
// Target body fixed coordinates and interpolations of Rotation Matrices
// Also tests PositionVector::representIn with and without spice

void Frame::testGeometryArrays(){

  
  Frame tf(target_frame_id_,target_origin_id_);

  GeometryModeE old = geom_mode_;
  geom_mode_=USE_SPICE; //set this so that checking computations will work 
                        // correctly


  // Initialize test state vector in target body fixed coordinates
  
  double max_vel_err=0.0, max_pos_err=0.0, max_rotmat_err=0.0;
  double max_pv_err=0.0;
  double time_max_pv_err, time_max_pos_err, time_max_vel_err,
    time_max_rotmat_err;
  SpiceDouble true_rotmat[3][3];
  SpiceDouble approx_rotmat[3][3];
  int beam_max_pv_err;
    for(int i=0;i<(num_time_steps_-1);i++){
      double t=start_time_+time_step_*i+ time_step_/2.0;
      double end_time=start_time_+time_step_*(num_time_steps_) -
	time_step_/2.0;
      Time t_time;

      Uvar x(0,km_str), y(0,km_str), z(0,km_str);
      PositionVector v1,v2;
      t_time.setEt(Uvar(t,seconds_str));
      PositionVector ontarget("ontarget",tf,t_time,x,y,z);
      updateTime(t_time);

      // set temporary ephemeris array
      double vel[3], pos[3];

      // Compute spacecraft state in target body fixed coordinates
      // using spice
      StateVector s;
      tf.ephemeris(s, cassini_origin_id_, t);
      pos[X]=get_in_base_units(s.r_[X]);
      pos[Y]=get_in_base_units(s.r_[Y]);
      pos[Z]=get_in_base_units(s.r_[Z]);
      vel[X]=get_in_base_units(s.v_[X]);
      vel[Y]=get_in_base_units(s.v_[Y]);
      vel[Z]=get_in_base_units(s.v_[Z]);


      // interpolate sc_vel_ and sc_pos_;
      interpolateSCPos();
      interpolateSCVel();

      double vel_err[3], pos_err[3];
      // compute ephemeris interpolation errors
      for (int j=0;j<3;j++){
	vel_err[j]=fabs(vel[j]-sc_vel_interp_[j]);
	if(vel_err[j]>max_vel_err){
	  max_vel_err=vel_err[j];
	  time_max_vel_err=t;
	}
	pos_err[j]=fabs(pos[j]-sc_pos_interp_[j]);
	if(pos_err[j]>max_pos_err){
	  max_pos_err=pos_err[j];
          time_max_pos_err=t;
	}
      }
      cout << "TIME " << t-start_time_ <<"(from start):" << endl;
      // loop over beams
      for(int b=0;b<5;b++){
	cout <<"   BEAM " << b+1 <<":" << endl;
        // set up beam frame and temporary rotmat array
	Frame bf(beam_frame_ids_[b],cassini_origin_id_);
	double b2t[3][3];

        // Compute Rotation Matrix for each beam from spice
	rotationMatrix(bf,tf,t,b2t);
        
	// Interpolate rotmat from rotmat_ array
	interpolateRotMat(b+1);


        // Compute rotmat interpolation errors
	int start_idx= 9*b;
 
	for(int j=0;j<3;j++){
	  for(int k=0;k<3;k++){
	    double err;

	    err=fabs(rotmat_interp_[start_idx+3*j+k]-b2t[j][k]);
	    cout << "      (j,k)=(" << j <<"," << k << 
	      ") (True,Interp,Error)=";
	    cout << "(" << b2t[j][k] << "," << rotmat_interp_[start_idx+3*j+k]
		 << "," << err << ")" << endl;
	    if(err> max_rotmat_err){ 
	      max_rotmat_err=err;
	      time_max_rotmat_err=t;
	    }
	  }
	}
	// compute error in PositionVector::representIn()
	geom_mode_=USE_SPICE;
	v1=ontarget;
	v1.representIn(bf);
        geom_mode_=PRECOMPUTE;
	v2=ontarget;
	v2.representIn(bf);
	geom_mode_=USE_SPICE;
	double pv_err=get_in_base_units((v1-v2).magnitude());
	if(pv_err>max_pv_err){
	  max_pv_err=pv_err;
          time_max_pv_err=t;
          beam_max_pv_err=b;
          for(int j=0;j<3;j++)
	    for(int k=0;k<3;k++){
	      true_rotmat[j][k]=b2t[j][k];
	      approx_rotmat[j][k]=rotmat_interp_[9*b+3*j+k];
	    }
	}
        cout << "Error in RepresentIn from TBF to Beam Frame is " << pv_err
	     << " km." << endl;
	cout << "   END BEAM " << b+1 << endl;
      }
      cout << "END TIME " << t << " Test " 
	   << 100*(t-start_time_)/(end_time-start_time_)
	   << "  Percent Complete " << endl;
    }

    // output error results
    cout << " Maximum 1-D position error is " << max_pos_err << " km " 
	 << time_max_pos_err -start_time_ << " s after start" << endl;
    cout << " Maximum 1-D velocity error is " << max_vel_err << " km/s " 
	 << time_max_vel_err -start_time_ << " s after start" << endl;
	 
    cout << " Maximum ROTMAT component error is " << max_rotmat_err 
	 << time_max_rotmat_err -start_time_ << " s after start" << endl;
    cout << " Maximum PositionVectorRepresentIn error is " << max_pv_err
	 << " km." << time_max_pv_err-start_time_ << " s after start" << endl;

    cout << "Breakdown of RepresentIn worst case......" << endl;
    cout << "True quaternion :";
    SpiceDouble q[4];
    m2q_c(true_rotmat,q);
    printQuaternion(q); 
    cout << endl;
    cout << "Quaternion computed from rotmat_interp_:";
    m2q_c(approx_rotmat,q);
    printQuaternion(q);
    //    interpolateQuaternion(quat_,4*beam_max_pv_err,QUAT_SIZE,
    //			  time_max_pv_err,q,1);
    geom_mode_=old;
}
  
void
Frame::setSpecialFrameNumber(){
  if(spice_frame_id_==target_frame_id_ 
     && spice_origin_id_==target_origin_id_) 
    {
      special_frame_number_=0;
      return;
    }
  
  else if(spice_origin_id_==cassini_origin_id_){
    for (int c=1;c<=5;c++){
      if (spice_frame_id_==beam_frame_ids_[c-1]){
	special_frame_number_=c;
	return;
      }
    }
  }

  special_frame_number_=-1;
  return;
}

//---------------
// Predicates
//---------------

bool Frame::operator==(const Frame& b) const // No Exceptions
  {
  return(frame_id_==b.frame_id_ && spice_origin_id_ == b.spice_origin_id_);
  }

bool Frame::operator!=(const Frame& b) const // No Exceptions
  {
  return(!(Frame::operator==(b)));
  }

bool Frame::isSpice() const // No Exceptions
  {
  return(spice_frame_);
  }

//---------------
// Other methods
//---------------

//-----------------------------
// name_str = name();
//
// Returns name of this Frame
//-----------------------------

string Frame::name() const // No Exceptions
  {
    string s;
    if(spice_frame_){
      spice_frame_name(frame_id_,s);
    }
    else{
      s="User defined frame #" + toStr(frame_id_);
    }
    return(s);
  }

//-----------------------------
// name_str = spiceFrameName();
//
// Returns name of the defining spice frame of this Frame
//-----------------------------

string Frame::spiceFrameName() const // No Exceptions
  {
    string s;
    spice_frame_name(spice_frame_id_,s);
    return(s);
  }
//-------------------------------------------------------------------------
// ephemeris(state,target_name,time,correction);
//
// Returns state (position,velocity) of target in this Frame at specified
// time, using the indicated correction ("LT+S","LT","NONE").
// Note, currently it is not possible to define a new frame that is moving
// with respect to a spice frame.  Only fixed rotations and translations
// are allowed.
//-------------------------------------------------------------------------

void Frame::ephemeris(StateVector& s, const string& target, const Time& t,
  const string& corr) const
  
  {
    if(strcasecmp(corr.c_str(),"none")!=0){
      ErrorMessage e("Frame::ephemeris No automatic light speed corrections");
      e.throwMe();
    }
    double t_in_s=get_in_base_units(t.et());
    SpiceInt target_id;
    spice_target_id(target,target_id);
    ephemeris(s,target_id,t_in_s);
  }

//-------------------------------------------------------------------------
// ephemeris(state,target_origin_id,time_in_s);
//
// Returns state (position,velocity) of target in this Frame at specified
// time, using no relativistic correction
// Note, currently it is not possible to define a new frame that is moving
// with respect to a spice frame.  Only fixed rotations and translations
// are allowed.
//-------------------------------------------------------------------------

void Frame::ephemeris(StateVector& s, SpiceInt target_id, 
		      double time_in_s) const
  
  {
  double cs=get_time();
  double ss[6];

  //----------------------------------------------------
  // First get ephemeris in defining spice frame
  //----------------------------------------------------

  // Fast geometry case Not using Spice, This is a Beam or Target Body
  // Fixed Frame,
  // and we are obtaining ephemeris for Cassini
  if(geom_mode_!=USE_SPICE && special_frame_number_>=0 && 
     target_id == Frame::cassini_origin_id_){

    // interpolate s/c ephemeris in Target Body Fixed frame
    if(special_frame_number_==0){
      Frame::updateTime(time_in_s);
      Frame::interpolateSCPos();
      Frame::interpolateSCVel();
      
      // Put interpolated positions and velocities into StateVector.
      
      
      ss[0]=sc_pos_interp_[0];
      ss[1]=sc_pos_interp_[1];
      ss[2]=sc_pos_interp_[2];
      ss[3]=sc_vel_interp_[0];
      ss[4]=sc_vel_interp_[1];
      ss[5]=sc_vel_interp_[2];
    }

    // S/C ephemeris in Beam Frame is all zeros
    else{
      ss[0]=0;
      ss[1]=0;
      ss[2]=0;
      ss[3]=0;
      ss[4]=0;
      ss[5]=0;
    }
  }
  
  // slow USE_SPICE geometry computation
  else{

    if(Frame::geom_mode_!=USE_SPICE && 
       special_frame_number_>0 && Frame::target_origin_id_==target_id){

      cerr << "Warning: Trying to find ephemeris of target body in beam frame"
	   << endl;
      cerr << "This is not supported with fast geometry calculations ..." 
	   << endl;
      cerr << "Spice will be called and computational efficiency sacrificed."
	   << endl;
    }
    SpiceDouble lt;

    //-------------------------
    // get state in a Spice frame
    //-------------------------
    char frname[33];
    spice_frame_name(spice_frame_id_,frname,33);
    spkgeo_c(target_id,time_in_s,frname,spice_origin_id_,ss,&lt);
    if(failed_c()){GeomError e(GeomError::spice_misc,time_in_s); e.throwMe();}
  }     

  //---------------------------------------------------
  // if not in a spice frame do additional computations
  //----------------------------------------------------
  if(!spice_frame_){
    SpiceDouble ro[3];
    ro[PositionVector::X] = rel_origin_[PositionVector::X];
    ro[PositionVector::Y] = rel_origin_[PositionVector::Y];
    ro[PositionVector::Z] = rel_origin_[PositionVector::Z];
	
    // Apply origin offset
    ss[0] -= ro[PositionVector::X];
    ss[1] -= ro[PositionVector::Y];
    ss[2] -= ro[PositionVector::Z];
	
    // Transform to this frame (which is fixed relative to the spice frame)
    mxv_c(spice_invxform_, ss, ss);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    mxv_c(spice_invxform_, &ss[3], &ss[3]);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  }

  // move position and velocity into StateVector.
  // (only fixed frame definitions allowed)
  
  s.r_[X] = Uvar(ss[0],km_str);
  s.r_[Y] = Uvar(ss[1],km_str);
  s.r_[Z] = Uvar(ss[2],km_str);
  s.v_[X] = Uvar(ss[3],km_per_s_str);
  s.v_[Y] = Uvar(ss[4],km_per_s_str);
  s.v_[Z] = Uvar(ss[5],km_per_s_str);
  
  // Set the remaining StateVector fields (except the name)
  s.setTimeInSeconds(time_in_s);
  s.frame_ = *this;
  double ce=get_time();
  frame_ephemeris_time+=ce-cs;
  }



//--------------------------------------------------------------------------
// q = quaternion(x,y,z);
//
// Returns an array of 4 doubles that represent the attitude given by
// the vectors x,y,z relative to *this frame as a unit quaternion.
//--------------------------------------------------------------------------

Array1D<double> Frame::quaternion(const DirectionVector& x,
  const DirectionVector& y, const DirectionVector& z) const
  
  {
  // Setup origin using *this so newframe also uses *this.
  PositionVector origin("origin",*this,x.time(),0,0,0);
  Frame newframe(origin,x,y,z);
  return(quaternion(newframe,x.time()));
  }

//--------------------------------------------------------------------------
// q = quaternion(f,t);
//
// Returns an array of 4 doubles that represent the attitude give by
// the frame f relative to *this frame as a unit quaternion.
//--------------------------------------------------------------------------

Array1D<double> Frame::quaternion(const Frame& f, const Time& t) const
  {
  SpiceDouble m[3][3];
  SpiceDouble et;
  t.getEt(et);
  rotationMatrix(*this,f,et,m);
  SpiceDouble quats[4];
  m2q_c(m,quats);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  Array1D<double> q("quaternion",4);
  q(0) = quats[0];
  q(1) = quats[1];
  q(2) = quats[2];
  q(3) = quats[3];
  return(q);
  }

//--------------------------------------------------------------------------
// axialVectors(f,t,x,y,z);
//
// Returns the axial vectors (x,y,z) of frame f at time t represented in
// *this frame.
//--------------------------------------------------------------------------

void Frame::axialVectors(const Frame& f, const Time& t,
  DirectionVector& x, DirectionVector& y, DirectionVector& z) const
  {
  SpiceDouble m[3][3];
  SpiceDouble et;
  t.getEt(et);
  rotationMatrix(*this,f,et,m);
  t.getEt(x.t_);
  y.t_ = x.t_;
  z.t_ = x.t_;
  x.frame_ = *this;
  y.frame_ = *this;
  z.frame_ = *this;
  x[DirectionVector::X] = m[0][0];
  x[DirectionVector::Y] = m[0][1];
  x[DirectionVector::Z] = m[0][2];
  y[DirectionVector::X] = m[1][0];
  y[DirectionVector::Y] = m[1][1];
  y[DirectionVector::Z] = m[1][2];
  z[DirectionVector::X] = m[2][0];
  z[DirectionVector::Y] = m[2][1];
  z[DirectionVector::Z] = m[2][2];
  }


//-------------------------------------------
// SLOW WARNING
// In order to avoid changing the spice calls
// this operation is inefficient when spice calls are
// performed
//----------------------------------------------------------------------------
// rotationMatrix(frame_from,frame_to,et,m)
//
// Return the rotation matrix that will
// transform a 3-element vector assumed to be in the 'from' frame
// to the indicated 'to' frame at the indicated ephemeris time.
// This rotation matrix does not account for any origin changes.
//----------------------------------------------------------------------------
 
void Frame::rotationMatrix(const Frame& frame_from, 
  const Frame& frame_to, SpiceDouble et, SpiceDouble m[3][3]) 
  // No Exceptions
  {

 // Fast geometry case: not using Spice,
 // One frame is target body fixed and the other is a beam frame
  
  if(geom_mode_!=USE_SPICE && frame_from.special_frame_number_>=0 
     && frame_to.special_frame_number_>=0 && 
     (frame_to.special_frame_number_==0 || 
      frame_from.special_frame_number_==0)){
    
    // identity transform special case
    if(frame_from.special_frame_number_==frame_to.special_frame_number_){
      m[0][0]=1; m[0][1]=0; m[0][2]=0;
      m[1][0]=0; m[1][1]=1; m[1][2]=0;
      m[2][0]=0; m[2][1]=0; m[2][2]=1;
    }
    else{
      bool beam_to_target;
      int bn;
      if(frame_from.special_frame_number_==0){
	beam_to_target=false;
	bn=frame_to.special_frame_number_;
      }
      else{
	beam_to_target=true;
	bn=frame_from.special_frame_number_;
      }
      Frame::updateTime(et);
      Frame::interpolateRotMat(bn);

      int offset1=9*(bn-1);
      for(int i=0;i<3;i++){ 
	for(int j=0;j<3;j++){
          // transpose intepolated matrix if target to beam desired
	  int offset2;
	  if(beam_to_target) offset2 = offset1+3*i+j;
	  else offset2 = offset1+3*j+i;
	  m[i][j]=rotmat_interp_[offset2];
	}
      }
    }
  }
  else if (frame_from.isSpice() && frame_to.isSpice())
    {  // both frames are in NAIF/Spice system
    pxform_c(frame_from.name().c_str(), frame_to.name().c_str(), et, m);
    if(failed_c()){GeomError e(GeomError::spice_misc,et); e.throwMe();}
    }
  else if (frame_from.isSpice())
    {  // only the from frame is a Spice frame
    // Transformation from spice frame
    pxform_c(frame_from.name().c_str(),frame_to.spiceFrameName().c_str(),et,m);
    if(failed_c())
      {
      GeomError e(GeomError::spice_misc,et); 
      e.throwMe();
      }
    // Combined transformation to the output (to) frame 
    mxm_c(frame_to.spice_invxform_, m, m);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    }
  else if (frame_to.isSpice())
    {  // only the destination (to) frame is a Spice frame
    // Transform to the output (to) frame 
    pxform_c(frame_from.spiceFrameName().c_str(),frame_to.name().c_str(),et,m);
    if(failed_c())
      {
      GeomError e(GeomError::spice_misc,et); 
      e.throwMe();
      }
    // Combined transformation to the output frame
    mxm_c(m, frame_from.spice_xform_, m);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    }
  else
    {  // Neither the from or to frames are Spice frames
    // Transform to the spice frame used to define the output (to) frame 
    pxform_c(frame_from.spiceFrameName().c_str(), 
	     frame_to.spiceFrameName().c_str(),
	     et, m);
    if(failed_c())
      {
      GeomError e(GeomError::spice_misc,et); 
      e.throwMe();
      }
    // Combined transformations 
    mxm_c(frame_to.spice_invxform_, m, m);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    mxm_c(m, frame_from.spice_xform_, m);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    }
  }


// Various routines for getting values directly from intermediate geometry 
// (and other) global values

// assumes sphere outputs value in km
double Frame::directAccessTargetRadius(const Time t){
  double r=get_in_base_units(default_target_radii.magnitude())/sqrt(3.0);
  return(r);
}

void Frame::directAccessRotMatB2T(const Time t,int bn,SpiceDouble m[3][3]){
  updateTime(t.et());
  interpolateRotMat(bn);
  int off=(bn-1)*9;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m[i][j]=rotmat_interp_[off+i*3+j];
    }
  }
}
void Frame::directAccessRotMatT2B(const Time t,int bn,SpiceDouble m[3][3]){
  updateTime(t.et());
  interpolateRotMat(bn);
  int off = (bn-1)*9;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m[j][i]=rotmat_interp_[off+i*3+j];
    }
  }
}
void Frame::directAccessSCPos(const Time t, double v[3]){
  updateTime(t.et());
  interpolateSCPos();
  for(int i=0;i<3;i++){
    v[i]=sc_pos_interp_[i];
  }
}
 void Frame::directAccessSCVel(const Time t, double v[3]){
  updateTime(t.et());
  interpolateSCVel();
  for(int i=0;i<3;i++){
    v[i]=sc_vel_interp_[i];
  }
}
 void Frame::directAccessSCAcc(const Time t, double v[3]){
  updateTime(t.et());
  int j=(int)floor(current_time_-start_time_);
  if(j<0 | j>= num_time_steps_){
    ErrorMessage 
      e("Attempt to interpolate outside time range of geometry arrays");
    e.throwMe();
  }
  double dj=current_time_-j*time_step_-start_time_;
  for(int i=0;i<3;i++){
    v[i]=sc_acc_[j*3+i]+sc_jerk_[j*3+i]*dj;
  }
}
//-----------------------------------------------------------------------------
// setSpiceRelativeRotation(rot,base_frame)
//
// Set this frame object to have the indicated Rotation relative to the
// base frame (ie., set spice_xform_, spice_invxform).
// The base frame may or may not be a spice frame.
//-----------------------------------------------------------------------------

void Frame::setSpiceRelativeRotation(const Rotation& rot,
  const Frame& base_frame)
  {
  DirectionVector x("",base_frame,Uvar(0),1,0,0);
  DirectionVector y("",base_frame,Uvar(0),0,1,0);
  DirectionVector z("",base_frame,Uvar(0),0,0,1);
  rot.rotatedAxes(x,y,z);
  setSpiceRelativeRotation(x,y,z,base_frame);
  }

//-----------------------------------------------------------------------------
// setSpiceRelativeRotation(x,y,z,base_frame)
//
// Set this frame object to have the indicated axes relative to the
// base frame (ie., set spice_xform_, spice_invxform).
// The base frame may or may not be a spice frame.
// An exception is thrown if the times of the input vectors don't match.
//-----------------------------------------------------------------------------

void Frame::setSpiceRelativeRotation(DirectionVector x, DirectionVector y,
  DirectionVector z, const Frame& base_frame)
  {
  if (x.timeInSeconds() != y.timeInSeconds() 
      || x.timeInSeconds() != z.timeInSeconds())
    {
      ErrorMessage e("Frame::setRelativeRotation: Time mismatch");
      e.throwMe();
    }

  x.representIn(base_frame);
  y.representIn(base_frame);
  z.representIn(base_frame);

  //----------------------------------------------------------------------
  // Set up the transformation to take a vector
  // back to the base frame (spice frame or not)
  //----------------------------------------------------------------------

  spice_xform_[0][0] = x[DirectionVector::X];
  spice_xform_[1][0] = x[DirectionVector::Y];
  spice_xform_[2][0] = x[DirectionVector::Z];

  spice_xform_[0][1] = y[DirectionVector::X];
  spice_xform_[1][1] = y[DirectionVector::Y];
  spice_xform_[2][1] = y[DirectionVector::Z];

  spice_xform_[0][2] = z[DirectionVector::X];
  spice_xform_[1][2] = z[DirectionVector::Y];
  spice_xform_[2][2] = z[DirectionVector::Z];

  //------------------------------------------------------------------------
  // If the base frame is itself defined relative to a Spice frame,
  // then we need to merge the two transformations.
  //------------------------------------------------------------------------

  if (!base_frame.isSpice())
    {  // base_frame is not in spice frame, so merge the transformations
    mxm_c(base_frame.spice_xform_, spice_xform_, spice_xform_);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    }

  //-------------------------------------------------------
  // Finally, setup the inverse transform (the transpose)
  //-------------------------------------------------------

  spice_invxform_[0][0] = spice_xform_[0][0];
  spice_invxform_[0][1] = spice_xform_[1][0];
  spice_invxform_[0][2] = spice_xform_[2][0];

  spice_invxform_[1][0] = spice_xform_[0][1];
  spice_invxform_[1][1] = spice_xform_[1][1];
  spice_invxform_[1][2] = spice_xform_[2][1];

  spice_invxform_[2][0] = spice_xform_[0][2];
  spice_invxform_[2][1] = spice_xform_[1][2];
  spice_invxform_[2][2] = spice_xform_[2][2];
  
  }

//-----------------------------------------------------------------------------
// setSpiceRelativeOrigin(origin)
//
// Set this frame object to have the indicated origin position.
// (ie., set rel_origin_).
// The origin position frame may or may not be a spice frame.
//-----------------------------------------------------------------------------

void Frame::setSpiceRelativeOrigin(const PositionVector& origin)
  {

    if(spice_frame_id_!=origin.frame_.spice_frame_id_ ||
       spice_origin_id_!=origin.frame_.spice_origin_id_){
      ErrorMessage e("setSpiceRelativeOrigin: Spice frame mismatch");
      e.throwMe();
    }

  //----------------------------------------------------------------------
  // Set up the translation to take a vector
  // back to the defining frame for the origin position (spice frame or not)
  //----------------------------------------------------------------------

  rel_origin_[X] = origin.km(PositionVector::X);
  rel_origin_[Y] = origin.km(PositionVector::Y);
  rel_origin_[Z] = origin.km(PositionVector::Z);

  //------------------------------------------------------------------------
  // If the defining frame is itself defined relative to a Spice frame,
  // then we need to merge the two translations.
  //------------------------------------------------------------------------

  if (!origin.frame_.isSpice())
    {
    // Move origin representation back to defining spice frame
    SpiceDouble rel_origin[3];
    string ustr(km_str);
    rel_origin[0] = rel_origin_[X];                
    rel_origin[1] = rel_origin_[Y];                
    rel_origin[2] = rel_origin_[Z];
    mxv_c(origin.frame_.spice_xform_, rel_origin, rel_origin);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    rel_origin_[X] = rel_origin[0];
    rel_origin_[Y] = rel_origin[1];
    rel_origin_[Z] = rel_origin[2];
    }

  }

void
Frame::getNewFrameID(){
  frame_id_=Frame::next_frame_id_++;
  if(frame_id_<0){
    ErrorMessage e("Too many user defined frames");
    e.throwMe();
  }
}

//----------------------------------
// Methods for class PositionVector
//----------------------------------

//--------------
// Constructors
//--------------

// Default constructor sets a default name, frame and time(set by Time default).
PositionVector::PositionVector()
  // No Exceptions
  : name_(NULL), frame_(), t_(0.0)
  {
  v_[X].setkm(0.0);
  v_[Y].setkm(0.0);
  v_[Z].setkm(0.0);
  }

// Default constructor sets a default frame and time(set by Time default).
PositionVector::PositionVector(const string& name)
  // No Exceptions
  : name_(NULL), frame_(), t_(0.0)
  {
    name_=new string;
    *name_=name;
    v_[X].setkm(0.0);
    v_[Y].setkm(0.0);
    v_[Z].setkm(0.0);
  }

// Fully specified constructor
PositionVector::PositionVector(const string& name, const Frame& frame,
  const Time& t, const Uvar& x, const Uvar& y, const Uvar& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(get_in_base_units(t))
  {
    name_=new string;
    *name_=name;
    v_[X] = x;
    v_[Y] = y;
    v_[Z] = z;
  }

PositionVector::PositionVector(const Frame& frame,
  const Time& t, const Uvar& x, const Uvar& y, const Uvar& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(get_in_base_units(t))
  {
    v_[X] = x;
    v_[Y] = y;
    v_[Z] = z;
  }

PositionVector::PositionVector(const Frame& frame,
  const double& t, const double& x, const double& y, const double& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(t)
  {
    v_[X].setkm(x);
    v_[Y].setkm(y);
    v_[Z].setkm(z);
  }

// Constructor that takes the frame, and time from another PositionVector
//-- No longer copies name 
PositionVector::PositionVector(const PositionVector& r,
  const Uvar& x, const Uvar& y, const Uvar& z)
  // No Exceptions
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  v_[X] = x;
  v_[Y] = y;
  v_[Z] = z;
  }

// copy constructor
// no longer copies name
PositionVector::PositionVector(const PositionVector& r)
  // No Exceptions
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  v_[X] = r.v_[X];
  v_[Y] = r.v_[Y];
  v_[Z] = r.v_[Z];
  }

// Type conversion DirectionVector -> PositionVector
// no longer copies name
PositionVector::PositionVector(const DirectionVector& u)
  // No Exceptions
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  }

// Type conversion FloatVector -> PositionVector
//---------- no longer copies name
PositionVector::PositionVector(const FloatVector& v)
  // No Exceptions
  : name_(NULL), frame_(v.frame_), t_(v.t_)
  {
  v_[X] = v.v_[X];
  v_[Y] = v.v_[Y];
  v_[Z] = v.v_[Z];
  }

// Destructor
PositionVector::~PositionVector()
  {
    if(name_) delete name_;
  }

//---------------
// Operators
//---------------

// unary -
PositionVector PositionVector::operator-() const // No Exceptions
  {
  PositionVector rv(*this);
  rv.v_[X] = -v_[X];
  rv.v_[Y] = -v_[Y];
  rv.v_[Z] = -v_[Z];
  return(rv);
  }

// +=
PositionVector& PositionVector::operator+=(const PositionVector& r)
  
  {
  v_[X] += r.v_[X];
  v_[Y] += r.v_[Y];
  v_[Z] += r.v_[Z];
  return(*this);
  }

// +=
PositionVector& PositionVector::operator+=(const Uvar& a)
  
  {
  v_[X] += a;
  v_[Y] += a;
  v_[Z] += a;
  return(*this);
  }

// -=
PositionVector& PositionVector::operator-=(const PositionVector& r)
  
  {
  v_[X] -= r.v_[X];
  v_[Y] -= r.v_[Y];
  v_[Z] -= r.v_[Z];
  return(*this);
  }

// -=
PositionVector& PositionVector::operator-=(const Uvar& a)
  
  {
  v_[X] -= a;
  v_[Y] -= a;
  v_[Z] -= a;
  return(*this);
  }


// *=
PositionVector& PositionVector::operator*=(const DirectionVector& r)
  
  {
    v_[X] *= r.v_[X];
    v_[Y] *= r.v_[Y];
    v_[Z] *= r.v_[Z];
    return(*this);
  }

// *=
PositionVector& PositionVector::operator*=(const double& a)
  
  {
  v_[X] *= a;
  v_[Y] *= a;
  v_[Z] *= a;
  return(*this);
  }

// *=
PositionVector& PositionVector::operator*=(const Uvar& a)
  
  {
  v_[X] *= a;
  v_[Y] *= a;
  v_[Z] *= a;
  return(*this);
  }

// Commented out because result is not a position
// /=
//PositionVector& PositionVector::operator/=(const PositionVector& r)
//  
//  {
//  v_[X] /= r.v_[X];
//  v_[Y] /= r.v_[Y];
//  v_[Z] /= r.v_[Z];
//  return(*this);
//  }

// /=
PositionVector& PositionVector::operator/=(const double& a)
  
  {
  v_[X] /= a;
  v_[Y] /= a;
  v_[Z] /= a;
  return(*this);
  }

// Assignment
PositionVector& PositionVector::operator=(const PositionVector& r) // No Exceptions
  {
  if (v_ == r.v_) return(*this);  // cover self assignment
  // Assignment leaves *this's name unchanged, but copies the data.
  frame_ = r.frame_;
  t_ = r.t_;
  v_[X] = r.v_[X];
  v_[Y] = r.v_[Y];
  v_[Z] = r.v_[Z];
  return(*this);
  }

// Indexing

Uvar& PositionVector::operator[](const indexE& i) // No Exceptions
  {
  return(v_[i]);
  }

const Uvar & PositionVector::operator[](const indexE& i)
  const // No Exceptions
  {
  return(v_[i]);

  }

double PositionVector::km(const indexE& i) const
{
  return(v_[i].km());
}

void PositionVector::setkm(const indexE& i, double value) 
{
  v_[i].setkm(value);
}

//---------------
// Predicates
//---------------

bool PositionVector::operator==(PositionVector b) const
  
  {
  b.representIn(frame_);
  return(t_==b.t_ && v_[X] == b.v_[X] &&
    v_[Y] == b.v_[Y] && v_[Z] == b.v_[Z]);
  }

bool PositionVector::operator!=(PositionVector b) const
  
  {
  return(!(PositionVector::operator==(b)));
  }

bool PositionVector::timeMatch(const PositionVector& b) const
  
  {
  return(t_ == b.t_);
  }

bool PositionVector::timeMatch(const DirectionVector& b) const
  
  {
  return(t_ == b.t_);
  }

bool PositionVector::timeMatch(const FloatVector& b) const
  
  {
  return(t_ == b.t_);
  }

bool PositionVector::frameMatch(const PositionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool PositionVector::frameMatch(const DirectionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool PositionVector::frameMatch(const FloatVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

//-----------------------------------------------------------------
// isSpice()
//
// Returns true if the frame this PositionVector is represented in
// is a Spice frame.
//-----------------------------------------------------------------

bool PositionVector::isSpice() const // No Exceptions
  {
  return(frame_.isSpice());
  }

//---------------
// Other methods
//---------------

//------------------------------------
// name_str = name();
//
// Returns name of this PositionVector
//------------------------------------

string PositionVector::name() const // No Exceptions
  {
    if(name_) return(*name_);
    else return("");
  }

//------------------------------------
// t = time();
//
// Returns time of this PositionVector
//------------------------------------

Time PositionVector::time() const // No Exceptions
  {
    Time t;
    t.setEt(t_);
    return(t);
  }

//------------------------------------
// setTime();
//
// Sets time of this PositionVector
//------------------------------------

void PositionVector::setTime(const Time& t) // No Exceptions
  {
   t.getEt(t_);
  }

//------------------------------------
// t = timeInSeconds();
//
// Returns time of this PositionVector in seconds
//------------------------------------

double PositionVector::timeInSeconds() const // No Exceptions
  {
     return(t_);
  }


void PositionVector::setTimeInSeconds(double t) // No Exceptions
  {
  t_ = t;
  }

//------------------------------------------------------------------
// name_str = frameName();
//
// Returns name of the frame this PositionVector is represented in.
//------------------------------------------------------------------

string PositionVector::frameName() const // No Exceptions
  {
  return(frame_.name());
  }


//----------------------------------------------------------------------------
// representInSpiceFrame(frame)
//
// Transform *this PositionVector object to the indicated frame defined at
// the same time.
//----------------------------------------------------------------------------

PositionVector& PositionVector::representInSpiceFrame()
  {
  if (frame_.isSpice()) return(*this);  // nothing to do


  SpiceDouble v[3];

  v[X]=km(X);
  v[Y]=km(Y);
  v[Z]=km(Z);

  mxv_c(frame_.spice_xform_, v, v);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}

  v[X]-=frame_.rel_origin_[X];
  v[Y]-=frame_.rel_origin_[Y];
  v[Z]-=frame_.rel_origin_[Z];

  setkm(X,v[X]);
  setkm(Y,v[Y]);
  setkm(Z,v[Z]);

  frame_=Frame(frame_.spice_frame_id_,frame_.spice_origin_id_);
  
  return(*this);  
  }
//----------------------------------------------------------------------------
// representIn(frame)
//
// Transform *this PositionVector object to the indicated frame defined at
// the same time.
//----------------------------------------------------------------------------

PositionVector& PositionVector::representIn(const Frame& frame)
  {
  if (frame_ == frame) return(*this);  // nothing to do

  SpiceDouble et=t_;



  SpiceDouble v[3];

  //-----------------------------------------------
  // If PositionVector is not in spice frame
  // first transform it to its defining spice frame
  //------------------------------------------------


  if(!frame_.spice_frame_){
    representInSpiceFrame();
  } 


   
  // Handles special case of target body fixed to beam frame
  // when speedup geometry is enabled
  if(Frame::geom_mode_!=Frame::USE_SPICE &&
     frame_.special_frame_number_==0 &&
     frame.special_frame_number_>0){

    StateVector s;
    frame_.ephemeris(s,frame.spice_origin_id_,et);
    // Add in origin offset
    v_[X] -= s.r_[X];
    v_[Y] -= s.r_[Y];
    v_[Z] -= s.r_[Z];

    v[X] = km(X);
    v[Y] = km(Y);
    v[Z] = km(Z);
    frameRotate(frame_,frame,et,v);

    setkm(X,v[X]);
    setkm(Y,v[Y]);
    setkm(Z,v[Z]);
  }
  else{

    v[X] = km(X);
    v[Y] = km(Y);
    v[Z] = km(Z);

    frameRotate(frame_,frame,et,v);

    setkm(X,v[X]);
    setkm(Y,v[Y]);
    setkm(Z,v[Z]);


    if (!frame.isSpice() || frame_.spice_origin_id_ != frame.spice_origin_id_)
      {  // need to translate to the new origin
	// Get old frame's origin in the new frame.
	StateVector s;
	frame.ephemeris(s,frame_.spice_origin_id_,et);
	// Add in origin offset
	v_[X] += s.r_[X];
	v_[Y] += s.r_[Y];
	v_[Z] += s.r_[Z];
      }
  }

  //------------------------------------------------------------------------
  // Set this's frame to the new frame (since this is now represented in the
  // new frame).
  //------------------------------------------------------------------------

  frame_ = frame;

  return(*this);
  }

//----------------------------------------------------------------------------
// representIn(r)
//
// Transform *this PositionVector object to the same frame that r is
// represented in.
// If the times of *this and r do not match, an exception is thrown.
//----------------------------------------------------------------------------

PositionVector& PositionVector::representIn(const PositionVector& r)
  {
  if (t_ != r.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + r.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(r.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this PositionVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

PositionVector& PositionVector::representIn(const DirectionVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this PositionVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

PositionVector& PositionVector::representIn(const FloatVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//-------------------------------------------------------
// magnitude()
//
// Return the magnitude of *this vector (including units)
//-------------------------------------------------------

Uvar PositionVector::magnitude() const 
  {
  return(sqrt(v_[X]*v_[X] + v_[Y]*v_[Y] + v_[Z]*v_[Z]));
  }

//-------------------------------------------------
// scale(mag)
//
// Scale *this to have the indicated magnitude.
// Throw an exception if attempting to scale zero.
//-------------------------------------------------

void PositionVector::scale(Uvar mag)
  {
  Uvar cur_mag = magnitude();
  if (cur_mag.getValue() == 0.0)
    {
    GeomError e("Tried to scale zero PositionVector: " + name(),
      GeomError::scale_zero);
    e.throwMe();
    }
  v_[X] *= mag/cur_mag;
  v_[Y] *= mag/cur_mag;
  v_[Z] *= mag/cur_mag;
  }

//-------------------------------------------------------------------------
// angle(v)
//
// Compute the subtended angle between *this and v.
// This method does NOT alter *this.
// Result is Uvar with units "rads".
//-------------------------------------------------------------------------

Uvar PositionVector::angle(const PositionVector& v) const
  {
  if (frame_ != v.frame_)
    {
    GeomError e("Tried to find angle between PositionVectors: " + name() +
      " and " + v.name() + " defined in different Frames (ambiguous)",
      GeomError::ambiguity_error);
    e.throwMe();
    }
  DirectionVector v1(*this);
  DirectionVector v2(v);
  return(Uvar(acos(dot(v1,v2)),"rad"));
  }

// Routine for rotating a Position Vector about a given direction vector
// Modifies the components of the calling vector and outputs a reference to it
PositionVector& PositionVector::rotateAbout(const DirectionVector& v, const Uvar& angle){
  SpiceDouble m[3][3];
  SpiceDouble vec[3];
  double t = 1-cos(angle);
  double s = sin(angle);
  double c = cos(angle);
  double x = v.v_[X];
  double y = v.v_[Y];
  double z = v.v_[Z];
  m[0][0]=t*x*x+c;
  m[0][1]=t*x*y+s*z;
  m[0][2]=t*x*z-s*y;
  m[1][0]=t*x*y-s*z;
  m[1][1]=t*y*y+c;
  m[1][2]=t*y*z+s*x;
  m[2][0]=t*x*y+s*y;
  m[2][1]=t*y*z-s*x;
  m[2][2]=t*z*z+c;
  vec[0]=get_in_base_units(v_[X]);
  vec[1]=get_in_base_units(v_[Y]);
  vec[0]=get_in_base_units(v_[Z]);
  mxv_c(m,vec,vec);
  v_[X]=Uvar(vec[0],"km");
  v_[Y]=Uvar(vec[1],"km");
  v_[Z]=Uvar(vec[2],"km");
  return(*this);
}

//-------------------------------------------------------------------------
// invert()
//
// Take the reciprocal of each element.
// If any element is zero, an exception is thrown.
//-------------------------------------------------------------------------

// Commenetd out because result is not a position
//void PositionVector::invert()
//  {
// if (v_[X].getValue() == 0 || v_[Y].getValue() == 0 || v_[Z].getValue() == 0)
//    {
//    GeomError e("Attempted to divide by zero element of PositionVector: " +
//     name());
//    e.throwMe();
//    }
//  v_[X] = 1.0/v_[X];
//  v_[Y] = 1.0/v_[Y];
//  v_[Z] = 1.0/v_[Z];
//  }

//----------------------------------------------------------------------------
// frameRotate(frame_from,frame_to,et,v)
//
// Transform a 3-element vector assumed to be in the 'from' frame
// to the indicated 'to' frame all defined at the indicated ephemeris time.
// This transformation represents the vector in the orientation of the new
// frame, but does not account for any origin changes.
//----------------------------------------------------------------------------

void PositionVector::frameRotate(const Frame& frame_from, 
  const Frame& frame_to, SpiceDouble et, SpiceDouble v[3])
  {

  SpiceDouble m[3][3];

  // get rotation matrix
  Frame::rotationMatrix(frame_from,frame_to,et,m);

  // multiply
  mxv_c(m,v,v);
 
  //cout << "frameRotate:" << endl;
  //cout << "  from to: " << frame_from.name() << " " << frame_to.name() << endl; 
  //cout << "  m = " << m[0][0] << " " << m[0][1] << " " << m[0][2] << endl;
  //cout << "      " << m[1][0] << " " << m[1][1] << " " << m[1][2] << endl;
  //cout << "      " << m[2][0] << " " << m[2][1] << " " << m[2][2] << endl;
  }

//------------------------------------
// Binary operators for PositionVector
//------------------------------------

// +
PositionVector operator+(PositionVector a1, const PositionVector& a2)
  {
  return(a1 += a2);
  }

// +
PositionVector operator+(PositionVector a1, const Uvar& a2)
  {
  return(a1 += a2);
  }

// +
PositionVector operator+(const Uvar& a1, PositionVector a2)
  {
  return(a2 += a1);
  }

// -
PositionVector operator-(PositionVector a1, const PositionVector& a2)
  {
  return(a1 -= a2);
  }

// -
PositionVector operator-(PositionVector a1, const Uvar& a2)
  {
  return(a1 -= a2);
  }

// -
PositionVector operator-(const Uvar& a1, PositionVector a2)
  {
  a2 = -a2;
  return(a2 += a1);
  }



// *
PositionVector operator*(PositionVector a1, const DirectionVector& a2)
  {
    return(a1 *= a2);
  }

// *
PositionVector operator*(PositionVector a1, const Uvar& a2)
  {
  return(a1 *= a2);
  }

// *
PositionVector operator*(const Uvar& a1, PositionVector a2)
  {
  return(a2 *= a1);
  }

// /
DirectionVector operator/(const PositionVector& a1, const PositionVector& a2)
  {
    double x= a1.km(PositionVector::X)/a2.km(PositionVector::X);
    double y= a1.km(PositionVector::Y)/a2.km(PositionVector::Y);
    double z= a1.km(PositionVector::Z)/a2.km(PositionVector::Z);
    DirectionVector d(a1.frame_,a1.t_,x,y,z);
    return(d);
  }

// /
PositionVector operator/(PositionVector a1, const Uvar& a2)
  {
  return(a1 /= a2.km());
  }

// commented out because result is not a position
// /
//PositionVector operator/(const Uvar& a1, PositionVector a2)
//  {
//  a2.invert();
//  return(a2 *= a1.km());
//  }

// <<
ostream& operator<<(ostream& s, const PositionVector& r)
  {
  s << "[" << r[PositionVector::X] << ","
    << r[PositionVector::Y] << ","
    << r[PositionVector::Z] << "]";
  return(s);
  }

//----------------------
// Supporting functions
//----------------------

//--------------------------------------------------------------------------
// val = dot(v1,v2)
//
// Form the vector dot product of v1 and v2 at time t, returned in val.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the dot product.
// Combinations of Position and Direction and Float Vectors are also supported
// with separate functions.
//--------------------------------------------------------------------------

Uvar dot(const PositionVector& v1, PositionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[PositionVector::X]*v2[PositionVector::X] +
         v1[PositionVector::Y]*v2[PositionVector::Y] +
         v1[PositionVector::Z]*v2[PositionVector::Z]);
  }

Uvar dot(const PositionVector& v1, DirectionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[PositionVector::X]*v2[DirectionVector::X] +
         v1[PositionVector::Y]*v2[DirectionVector::Y] +
         v1[PositionVector::Z]*v2[DirectionVector::Z]);
  }

Uvar dot(const DirectionVector& v1, PositionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[DirectionVector::X]*v2[PositionVector::X] +
         v1[DirectionVector::Y]*v2[PositionVector::Y] +
         v1[DirectionVector::Z]*v2[PositionVector::Z]);
  }

Uvar dot(const PositionVector& v1, FloatVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[PositionVector::X]*v2[FloatVector::X] +
         v1[PositionVector::Y]*v2[FloatVector::Y] +
         v1[PositionVector::Z]*v2[FloatVector::Z]);
  }

Uvar dot(const FloatVector& v1, PositionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[FloatVector::X]*v2[PositionVector::X] +
         v1[FloatVector::Y]*v2[PositionVector::Y] +
         v1[FloatVector::Z]*v2[PositionVector::Z]);
  }

Uvar dot(const FloatVector& v1, DirectionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[FloatVector::X]*v2[DirectionVector::X] +
         v1[FloatVector::Y]*v2[DirectionVector::Y] +
         v1[FloatVector::Z]*v2[DirectionVector::Z]);
  }

Uvar dot(const DirectionVector& v1, FloatVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[DirectionVector::X]*v2[FloatVector::X] +
         v1[DirectionVector::Y]*v2[FloatVector::Y] +
         v1[DirectionVector::Z]*v2[FloatVector::Z]);
  }

//--------------------------------------------------------------------------
// v = cross(v1,v2)
//
// Form the vector cross product of v1 and v2, returned in v.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the cross product.  Thus, the result will be in v1's frame.
// Combinations of Position and Direction and Float Vectors are also supported
// with separate functions.
// The result of a cross product is always a FloatVector except for the
// special case of crossing two DirectionVectors which yields another
// DirectionVector.
// The time of the result matches the times of the inputs which are required
// to match.
//--------------------------------------------------------------------------

FloatVector cross(const PositionVector& v1, PositionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[PositionVector::Y] * v2[PositionVector::Z] -
    v1[PositionVector::Z] * v2[PositionVector::Y];
  result[FloatVector::Y] = v1[PositionVector::Z] * v2[PositionVector::X] -
    v1[PositionVector::X] * v2[PositionVector::Z];
  result[FloatVector::Z] = v1[PositionVector::X] * v2[PositionVector::Y] -
    v1[PositionVector::Y] * v2[PositionVector::X];
  return(result);
  }

FloatVector cross(const PositionVector& v1, DirectionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[PositionVector::Y] * v2[DirectionVector::Z] -
    v1[PositionVector::Z] * v2[DirectionVector::Y];
  result[FloatVector::Y] = v1[PositionVector::Z] * v2[DirectionVector::X] -
    v1[PositionVector::X] * v2[DirectionVector::Z];
  result[FloatVector::Z] = v1[PositionVector::X] * v2[DirectionVector::Y] -
    v1[PositionVector::Y] * v2[DirectionVector::X];
  return(result);
  }

FloatVector cross(const DirectionVector& v1, PositionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[DirectionVector::Y] * v2[PositionVector::Z] -
    v1[DirectionVector::Z] * v2[PositionVector::Y];
  result[FloatVector::Y] = v1[DirectionVector::Z] * v2[PositionVector::X] -
    v1[DirectionVector::X] * v2[PositionVector::Z];
  result[FloatVector::Z] = v1[DirectionVector::X] * v2[PositionVector::Y] -
    v1[DirectionVector::Y] * v2[PositionVector::X];
  return(result);
  }

FloatVector cross(const PositionVector& v1, FloatVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[PositionVector::Y] * v2[FloatVector::Z] -
    v1[PositionVector::Z] * v2[FloatVector::Y];
  result[FloatVector::Y] = v1[PositionVector::Z] * v2[FloatVector::X] -
    v1[PositionVector::X] * v2[FloatVector::Z];
  result[FloatVector::Z] = v1[PositionVector::X] * v2[FloatVector::Y] -
    v1[PositionVector::Y] * v2[FloatVector::X];
  return(result);
  }

FloatVector cross(const FloatVector& v1, PositionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }
  
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[FloatVector::Y] * v2[PositionVector::Z] -
    v1[FloatVector::Z] * v2[PositionVector::Y];
  result[FloatVector::Y] = v1[FloatVector::Z] * v2[PositionVector::X] -
    v1[FloatVector::X] * v2[PositionVector::Z];
  result[FloatVector::Z] = v1[FloatVector::X] * v2[PositionVector::Y] -
    v1[FloatVector::Y] * v2[PositionVector::X];
  return(result);
  }

FloatVector cross(const FloatVector& v1, DirectionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[FloatVector::Y] * v2[DirectionVector::Z] -
    v1[FloatVector::Z] * v2[DirectionVector::Y];
  result[FloatVector::Y] = v1[FloatVector::Z] * v2[DirectionVector::X] -
    v1[FloatVector::X] * v2[DirectionVector::Z];
  result[FloatVector::Z] = v1[FloatVector::X] * v2[DirectionVector::Y] -
    v1[FloatVector::Y] * v2[DirectionVector::X];
  return(result);
  }

FloatVector cross(const DirectionVector& v1, FloatVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[DirectionVector::Y] * v2[FloatVector::Z] -
    v1[DirectionVector::Z] * v2[FloatVector::Y];
  result[FloatVector::Y] = v1[DirectionVector::Z] * v2[FloatVector::X] -
    v1[DirectionVector::X] * v2[FloatVector::Z];
  result[FloatVector::Z] = v1[DirectionVector::X] * v2[FloatVector::Y] -
    v1[DirectionVector::Y] * v2[FloatVector::X];
  return(result);
  }

//----------------------------------
// Methods for class DirectionVector
//----------------------------------

//--------------
// Constructors
//--------------

// Default constructor sets a default frame and time.
DirectionVector::DirectionVector()
  // No Exceptions
  : name_(NULL), frame_(), t_(0.0)
  {
  v_[X] = 0.0;
  v_[Y] = 0.0;
  v_[Z] = 0.0;
  }

// Default constructor sets a default frame and time.
DirectionVector::DirectionVector(const string& name)
  // No Exceptions
  : name_(NULL), frame_(), t_(0.0)
  {
    name_=new string;
    *name_=name;
    v_[X] = 0.0;
    v_[Y] = 0.0;
    v_[Z] = 0.0;
  }


// Fully specified constructor
DirectionVector::DirectionVector(const string& name, const Frame& frame,
  const double& t, const double& x, const double& y, const double& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(t)
  {
    name_=new string;
    *name_=name;
    v_[X]=x;
    v_[Y]=y;
    v_[Z]=z;
    unity_scale(true);
  }

DirectionVector::DirectionVector(const string& name, const Frame& frame,
  const Time& t, const double& x, const double& y, const double& z)
  // No Exceptions
  : name_(NULL), frame_(frame)
  {
    name_=new string;
    *name_=name;
    t.getEt(t_);
    v_[X]=x;
    v_[Y]=y;
    v_[Z]=z;
    unity_scale(true);
  }

// Fully specified except name
DirectionVector::DirectionVector(const Frame& frame,
  const double& t, const double& x, const double& y, const double& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(t)
  {
    v_[X]=x;
    v_[Y]=y;
    v_[Z]=z;
    unity_scale(true);
  }

DirectionVector::DirectionVector(const Frame& frame,
  const Time& t, const double& x, const double& y, const double& z)
  // No Exceptions
  : name_(NULL), frame_(frame)
  {
    t.getEt(t_);
    v_[X]=x;
    v_[Y]=y;
    v_[Z]=z;
    unity_scale(true);
  }

// Construct from PositionVector by scaling to unit magnitude
DirectionVector::DirectionVector(const string& name, const PositionVector& r)
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
    name_=new string;
    *name_=name;
    v_[X] = r.v_[X].km();
    v_[Y] = r.v_[Y].km();
    v_[Z] = r.v_[Z].km();
    unity_scale();
  }

// copy constructor with new name
DirectionVector::DirectionVector(const string& name, const DirectionVector& u)
  // No Exceptions
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  name_=new string;
  *name_=name;
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  }

// copy from FloatVector with new name
DirectionVector::DirectionVector(const string& name, const FloatVector& u)
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
    name_=new string;
    *name_=name;
    v_[X] = get_in_base_units(u.v_[X]);
    v_[Y] = get_in_base_units(u.v_[Y]);
    v_[Z] = get_in_base_units(u.v_[Z]);
    unity_scale();
  }

// copy constructor
DirectionVector::DirectionVector(const DirectionVector& u)
  // No Exceptions
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  }

// copy from PositionVector constructor
DirectionVector::DirectionVector(const PositionVector& r)
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  v_[X] = r.v_[X].km();
  v_[Y] = r.v_[Y].km();
  v_[Z] = r.v_[Z].km();
  unity_scale();
  }

// copy from FloatVector
DirectionVector::DirectionVector(const FloatVector& u)
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  v_[X] = get_in_base_units(u.v_[X]);
  v_[Y] = get_in_base_units(u.v_[Y]);
  v_[Z] = get_in_base_units(u.v_[Z]);
  unity_scale();
  }

// Destructor
DirectionVector::~DirectionVector()
  {
    if(name_) delete name_;
  }

//---------------
// Operators
//---------------

// unary -
DirectionVector DirectionVector::operator-() const // No Exceptions
  {
  DirectionVector rv(*this);
  rv.v_[X] = -v_[X];
  rv.v_[Y] = -v_[Y];
  rv.v_[Z] = -v_[Z];
  return(rv);
  }

// +=
DirectionVector& DirectionVector::operator+=(const DirectionVector& r)
  // No Exceptions
  {
  v_[X] += r.v_[X];
  v_[Y] += r.v_[Y];
  v_[Z] += r.v_[Z];
  unity_scale();
  return(*this);
  }

// +=
//DirectionVector& DirectionVector::operator+=(const double& a)
// No Exceptions
//{
//v_[X] += a;
//v_[Y] += a;
//v_[Z] += a;
//  return(*this);
//}

// -=
DirectionVector& DirectionVector::operator-=(const DirectionVector& r)
  // No Exceptions
  {
  v_[X] -= r.v_[X];
  v_[Y] -= r.v_[Y];
  v_[Z] -= r.v_[Z];
  unity_scale();
  return(*this);
  }

// -=
//DirectionVector& DirectionVector::operator-=(const double& a)
// No Exceptions
//  {
//v_[X] -= a;
//v_[Y] -= a;
//v_[Z] -= a;
//return(*this);
//}

// *=
DirectionVector& DirectionVector::operator*=(const DirectionVector& r)
  // No Exceptions
  {
  v_[X] *= r.v_[X];
  v_[Y] *= r.v_[Y];
  v_[Z] *= r.v_[Z];
  unity_scale();
  return(*this);
  }

// *=
//DirectionVector& DirectionVector::operator*=(const double& a)
// No Exceptions
//  {
//v_[X] *= a;
//v_[Y] *= a;
//v_[Z] *= a;
//return(*this);
//}

// /=
DirectionVector& DirectionVector::operator/=(const DirectionVector& r)
  // No Exceptions
  {
  v_[X] /= r.v_[X];
  v_[Y] /= r.v_[Y];
  v_[Z] /= r.v_[Z];
  unity_scale();
  return(*this);
  }

// /=
//DirectionVector& DirectionVector::operator/=(const double& a)
// No Exceptions
//  {
//v_[X] /= a;
//v_[Y] /= a;
//v_[Z] /= a;
//return(*this);
//}

// Assignment
DirectionVector& DirectionVector::operator=(const DirectionVector& u) // No Exceptions
  {
  if (v_ == u.v_) return(*this);  // cover self assignment
  // Assignment leaves *this's name unchanged, but copies the data.
  frame_ = u.frame_;
  t_ = u.t_;
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  return(*this);
  }

// Indexing

double& DirectionVector::operator[](const indexE& i) // No Exceptions
  {
  return(v_[i]);
  }

const double& DirectionVector::operator[](const indexE& i)
  const // No Exceptions
  {
  return(v_[i]);
  }

//---------------
// Predicates
//---------------

bool DirectionVector::operator==(DirectionVector b) const
  {
  b.representIn(frame_);
  return(timeMatch(b) && v_[X] == b.v_[X] &&
    v_[Y] == b.v_[Y] && v_[Z] == b.v_[Z]);
  }

bool DirectionVector::operator!=(DirectionVector b) const
  {
  return(!(DirectionVector::operator==(b)));
  }

bool DirectionVector::timeMatch(const DirectionVector& b) const
  {
  return(t_ == b.t_);
  }

bool DirectionVector::timeMatch(const PositionVector& b) const
  {
  return(t_ == b.t_);
  }

bool DirectionVector::timeMatch(const FloatVector& b) const
  {
  return(t_ == b.t_);
  }

bool DirectionVector::timeMatch(const double& t) const
  {
  return(t_ == t);
  }

bool DirectionVector::timeMatch(const Time& t) const
  {
    double et;
    t.getEt(et);
    return(t_ == et );
  }

bool DirectionVector::frameMatch(const DirectionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool DirectionVector::frameMatch(const PositionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool DirectionVector::frameMatch(const FloatVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

//-----------------------------------------------------------------
// isSpice()
//
// Returns true if the frame this DirectionVector is represented in
// is a Spice frame.
//-----------------------------------------------------------------

bool DirectionVector::isSpice() const // No Exceptions
  {
  return(frame_.isSpice());
  }

//-----------------
// Get/Set methods
//-----------------

//------------------------------------------------------------------------
// {get,set}Planetodetic(lat,lon)
//
// Get and set this DirectionVector to point at the indicated
// coordinates.  Latitude is planetodetic and varies from -90 to +90 deg
// with the positive direction along the z-axis.
// Longitude varies from -180 deg to +180 deg with zero along the x-axis
// and positive following the right hand rule (ie., towards +y).
//------------------------------------------------------------------------

void DirectionVector::getPlanetodetic(Uvar& lat, Uvar& lon) const
  {
  Uvar theta;
  getSpherical(theta,lon);
  lat = Uvar(pi/2,"rad") - theta;
  }

void DirectionVector::setPlanetodetic(const Uvar& lat, const Uvar& lon)
  {
  Uvar theta = Uvar(pi/2,"rad") - lat;
  setSpherical(theta,lon);
  }

//------------------------------------------------------------------------
// {get,set}AzimuthElevation(azi,elev)
//
// Get and set this DirectionVector to point at the indicated
// azimuth and elevation direction in *this's frame.
// Azimuth and elevation are defined for the Cassini RADAR as follows:
// Elevation is the complement of the angle between the *this and the
// y-axis of *this's frame.   It varies from -90 (along -y) to +90 (along +y).
// ie., elevation is the angle from the x-z plane.
// Azimuth is the angle between the projection of *this onto the x-z plane
// of this's frame and the z-axis of *this's frame.  Positive azimuth
// follows the right hand rule with respect to the y-axis.  ie., positive
// azimuth means a deflection from the +z-axis towards the +x-axis.
// Azimuth is defined over 360 degrees, but in practise, values will be
// restricted to -90 (along -x) to 90 (along +x).
//------------------------------------------------------------------------

void DirectionVector::getAzimuthElevation(Uvar& azi, Uvar& elev) const
  {
  elev = Uvar(asin(v_[Y]),"rad");
  azi = Uvar(asin(v_[X]/cos(elev)),"rad");  
  }

void DirectionVector::setAzimuthElevation(const Uvar& azi, const Uvar& elev)
  {
  v_[X] = sin(azi)*cos(elev);
  v_[Y] = sin(elev);
  v_[Z] = cos(azi)*cos(elev);
  }

//------------------------------------------------------------------------
// {get,set}RADEC(ra,dec)
//
// Get and set this DirectionVector to point at the indicated
// right ascension (RA) and declination (DEC).
// RA and DEC are closely related to standard spherical coordinates.
// RA is the standard spherical phi expressed in the range [0,360] deg.
// DEC is the complement of the standard sherical theta.
// Technically, RA,DEC coordinates are only valid in the earth centered
// J2000 frame, but they can be used with any coordinate frame if desired.
//------------------------------------------------------------------------

void DirectionVector::getRADEC(Uvar& ra,
  Uvar& dec) const 
  {
  getSpherical(dec,ra);
  if (ra < Uvar(0,"rad")) ra += Uvar(2*pi,"rad");
  dec = Uvar(pi/2,"rad") - dec;
  }

void DirectionVector::setRADEC(const Uvar& ra,
  const Uvar& dec) 
  {
  setSpherical(Uvar(pi/2,"rad") - dec, ra);
  }

//------------------------------------------------------------------------
// {get,set}Spherical(theta,phi)
//
// Get and set this DirectionVector to point at the indicated
// standard spherical coordinates.  theta is the angle from the +z-axis
// and varies from 0 to 180 deg.
// phi varies from -180 deg to +180 deg with zero along the x-axis
// and positive following the right hand rule (ie., towards +y).
//------------------------------------------------------------------------

void DirectionVector::getSpherical(Uvar& theta, Uvar& phi) const
  {
  theta = Uvar(acos(v_[Z]),"rad");
  phi = Uvar(atan2(v_[Y],v_[X]),"rad");
  }

void DirectionVector::setSpherical(const Uvar& theta, const Uvar& phi)
  {
  v_[X] = sin(theta)*cos(phi);
  v_[Y] = sin(theta)*sin(phi);
  v_[Z] = cos(theta);
  }

//---------------
// Other methods
//---------------

//-------------------------------------
// name_str = name();
//
// Returns name of this DirectionVector
//-------------------------------------

string DirectionVector::name() const // No Exceptions
  {
    if(name_) return(*name_);    
    else return("");
  }

//------------------------------------
// t = time();
//
// Returns time of this DirectionVector
//------------------------------------

Time DirectionVector::time() const // No Exceptions
  {
    Time t;
    t.setEt(t_);
    return(t);
  }

//------------------------------------
// t = timeInSeconds();
//
// Returns time in seconds of this DirectionVector
//------------------------------------


//------------------------------------
// setTime();
//
// Sets time of this DirectionVector
//------------------------------------

void DirectionVector::setTime(const Time& t) // No Exceptions
  {
   t.getEt(t_);
  }

//------------------------------------
// t = timeInSeconds();
//
// Returns time of this PositionVector in seconds
//------------------------------------

double DirectionVector::timeInSeconds() const // No Exceptions
  {
     return(t_);
  }


void DirectionVector::setTimeInSeconds(double t) // No Exceptions
  {
  t_ = t;
  }

//------------------------------------------------------------------
// name_str = frameName();
//
// Returns name of the frame this DirectionVector is represented in.
//------------------------------------------------------------------

string DirectionVector::frameName() const // No Exceptions
  {
  return(frame_.name());
  }


// Routine for rotating a Direction Vector about a given direction vector
// Modifies the components of the calling vector and outputs a reference to it
DirectionVector& DirectionVector::rotateAbout(const DirectionVector& v, const Uvar& angle){
  SpiceDouble m[3][3];
  double t = 1-cos(angle);
  double s = sin(angle);
  double c = cos(angle);
  double x = v.v_[X];
  double y = v.v_[Y];
  double z = v.v_[Z];
  m[0][0]=t*x*x+c;
  m[0][1]=t*x*y+s*z;
  m[0][2]=t*x*z-s*y;
  m[1][0]=t*x*y-s*z;
  m[1][1]=t*y*y+c;
  m[1][2]=t*y*z+s*x;
  m[2][0]=t*x*y+s*y;
  m[2][1]=t*y*z-s*x;
  m[2][2]=t*z*z+c;
  mxv_c(m,v_,v_);
  return(*this);
}
//----------------------------------------------------------------------------
// representIn(frame)
//
// Transform *this DirectionVector object to the indicated frame defined at
// the same time.
//----------------------------------------------------------------------------

DirectionVector& DirectionVector::representIn(const Frame& frame)
  {
  if (frame_ == frame) return(*this);  // nothing to do

  SpiceDouble et;
  et=t_;
  PositionVector::frameRotate(frame_,frame,et,v_);
  unity_scale();

  //------------------------------------------------------------------------
  // Set this's frame to the new frame (since this is now represented in the
  // new frame).
  //------------------------------------------------------------------------

  frame_ = frame;

  return(*this);
  }

//----------------------------------------------------------------------------
// representIn(r)
//
// Transform *this DirectionVector object to the same frame that r is
// represented in.
// If the times of *this and r do not match, an exception is thrown.
//----------------------------------------------------------------------------

DirectionVector& DirectionVector::representIn(const PositionVector& r)
  {
  if (t_ != r.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + r.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(r.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this DirectionVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

DirectionVector& DirectionVector::representIn(const DirectionVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this DirectionVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

DirectionVector& DirectionVector::representIn(const FloatVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//-------------------------------------------------------------------------
// angle(v)
//
// Compute the subtended angle between *this and v.
// Result is Uvar with units "rads".
// If the times of *this and v do not match, an exception is thrown.
//-------------------------------------------------------------------------

Uvar DirectionVector::angle(DirectionVector v) const
  {
  v.representIn(frame_);
  return(Uvar(acos(dot(*this,v)),"rad"));
  }

//-------------------------------------------------------------------------
// invert()
//
// Take the reciprocal of each element.
// If any element is zero, an exception is thrown.
//-------------------------------------------------------------------------

void DirectionVector::invert()
  {
  if (v_[X] == 0 || v_[Y] == 0 || v_[Z] == 0)
    {
    GeomError e("Attempted to divide by zero element of DirectionVector: " +
      name());
    e.throwMe();
    }
  v_[X] = 1.0/v_[X];
  v_[Y] = 1.0/v_[Y];
  v_[Z] = 1.0/v_[Z];
  unity_scale();
  }

//----------------------------------------------------------------------------
// dptr = data()
//
// Return a pointer to the internal vector of doubles.
// This violates encapsulation, but is convenient when using these
// vectors with spice library routines.
//----------------------------------------------------------------------------

double* DirectionVector::data()
  {
  return(v_);
  }

//----------------------------------------------------
// unity_scale(allow_zero)
//
// Scale *this to have unit magnitude.
// DirectionVectors should always have unit magnitude.
// allow_zero if true lets zero vectors get by unchanged.
//   By default, allow_zero is false.
//----------------------------------------------------

void DirectionVector::unity_scale(bool allow_zero)
  {
  double cur_mag2 = v_[0]*v_[0] + v_[1]*v_[1] + v_[2]*v_[2];

  // check if already unity magnitude if so do noting
  if(cur_mag2==1.0) return;


  if (cur_mag2 == 0.0)
    {
      if(!allow_zero){
	GeomError e("Tried to scale zero DirectionVector: " + name(),
		    GeomError::scale_zero);
        //double xx = 1/0;
	e.throwMe();
      }
      else return;
    }

  double cur_mag=sqrt(cur_mag2);
  v_[X] /= cur_mag;
  v_[Y] /= cur_mag;
  v_[Z] /= cur_mag;
  }

//------------------------------------
// Binary operators for DirectionVector
//------------------------------------

// +
DirectionVector operator+(DirectionVector a1, const DirectionVector& a2)
  // No Exceptions
  {
  return(a1 += a2);
  }

// +
//DirectionVector operator+(DirectionVector a1,
//const double& a2)
// No Exceptions
//  {
//return(a1 += a2);
//}

// +
//DirectionVector operator+(const double& a1,
//DirectionVector a2)
// No Exceptions
//  {
//return(a2 += a1);
//}

// -
DirectionVector operator-(DirectionVector a1, const DirectionVector& a2)
  // No Exceptions
  {
  return(a1 -= a2);
  }

// -
//DirectionVector operator-(DirectionVector a1,
//const double& a2)
// No Exceptions
//{
//return(a1 -= a2);
//}

// -
//DirectionVector operator-(const double& a1,
//DirectionVector a2)
// No Exceptions
//{
//a2 = -a2;
//return(a2 += a1);
//}

// *
DirectionVector operator*(DirectionVector a1, const DirectionVector& a2)
  // No Exceptions
  {
  return(a1 *= a2);
  }

// *
//DirectionVector operator*(DirectionVector a1,
//const double& a2)
// No Exceptions
//{
//return(a1 *= a2);
//}

// *
//DirectionVector operator*(const double& a1,
//DirectionVector a2)
// No Exceptions
//{
//return(a2 *= a1);
//}

// /
DirectionVector operator/(DirectionVector a1, const DirectionVector& a2)
  // No Exceptions
  {
  return(a1 /= a2);
  }

// /
//DirectionVector operator/(DirectionVector a1,
//const double& a2)
// No Exceptions
//{
//return(a1 /= a2);
//}

// /
//DirectionVector operator/(const double& a1,
//DirectionVector a2)
// No Exceptions
//{
//a2.invert();
//return(a2 *= a1);
//}

// <<
ostream& operator<<(ostream& s, const DirectionVector& u)
  {
  s << "[" << u[DirectionVector::X] << ","
    << u[DirectionVector::Y] << ","
    << u[DirectionVector::Z] << "]";
  return(s);
  }

//----------------------
// Supporting functions
//----------------------

//--------------------------------------------------------------------------
// val = dot(v1,v2)
//
// Form the vector dot product of v1 and v2, returned in val.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the dot product.
//--------------------------------------------------------------------------

double dot(const DirectionVector& v1, DirectionVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  double retval=v1[DirectionVector::X]*v2[DirectionVector::X] +
         v1[DirectionVector::Y]*v2[DirectionVector::Y] +
         v1[DirectionVector::Z]*v2[DirectionVector::Z];

  // force dot product of direction vectors to be in range [-1 1]
  // eliminates pernicious bug which leads to acos(1.000000000002)=nan
  if(retval>1.0) retval=1.0;
  else if(retval<-1.0) retval=-1.0;
  return(retval);
  }

//--------------------------------------------------------------------------
// v = cross(v1,v2)
//
// Form the vector cross product of v1 and v2, returned in v.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the cross product.  Thus, the result will be in v1's frame.
//--------------------------------------------------------------------------

DirectionVector cross(const DirectionVector& v1, DirectionVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  DirectionVector result(v1);

  result[DirectionVector::X] = v1[DirectionVector::Y] * v2[DirectionVector::Z] -
    v1[DirectionVector::Z] * v2[DirectionVector::Y];
  result[DirectionVector::Y] = v1[DirectionVector::Z] * v2[DirectionVector::X] -
    v1[DirectionVector::X] * v2[DirectionVector::Z];
  result[DirectionVector::Z] = v1[DirectionVector::X] * v2[DirectionVector::Y] -
    v1[DirectionVector::Y] * v2[DirectionVector::X];
  double mag = sqrt(result[DirectionVector::X]*result[DirectionVector::X]+
		    result[DirectionVector::Y]*result[DirectionVector::Y]+
		    result[DirectionVector::Z]*result[DirectionVector::Z]);
  //need to make the magnitude to be 1
  result[DirectionVector::X]/=mag;
  result[DirectionVector::Y]/=mag;
  result[DirectionVector::Z]/=mag;
  return(result);
  }

//----------------------------------
// Methods for class FloatVector
//----------------------------------

//--------------
// Constructors
//--------------

// Default constructor sets a default frame and time.
FloatVector::FloatVector()
  // No Exceptions
  : name_(NULL), frame_()
  {
  v_[X] = 0.0;
  v_[Y] = 0.0;
  v_[Z] = 0.0;
  }

// Default constructor sets a default frame and time.
FloatVector::FloatVector(const string& name)
  // No Exceptions
  : name_(NULL), frame_("J2000","Earth")
  {
  name_=new string;
  *name_=name;
  v_[X] = 0.0;
  v_[Y] = 0.0;
  v_[Z] = 0.0;
  }

// Fully specified constructor
FloatVector::FloatVector(const string& name, const Frame& frame,
  const Time& t, const Uvar& x, const Uvar& y, const Uvar& z)
  // No Exceptions
  : name_(NULL), frame_(frame)
  {
  name_=new string;
  *name_=name;
  t.getEt(t_);
  v_[X] = x;
  v_[Y] = y;
  v_[Z] = z;
  }

// Fully specified constructor (except name)
FloatVector::FloatVector(const Frame& frame,
  const double& t, const Uvar& x, const Uvar& y, const Uvar& z)
  // No Exceptions
  : name_(NULL), frame_(frame), t_(t)
  {
  v_[X] = x;
  v_[Y] = y;
  v_[Z] = z;
  }

// Construct from PositionVector
FloatVector::FloatVector(const string& name, const PositionVector& r)
  // No Exceptions
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  name_=new string;
  *name_=name;
  v_[X] = r.v_[X];
  v_[Y] = r.v_[Y];
  v_[Z] = r.v_[Z];
  }

// copy constructor with new name
FloatVector::FloatVector(const string& name, const FloatVector& u)
  // No Exceptions
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  name_=new string;
  *name_=name;
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  }

// copy constructor
FloatVector::FloatVector(const FloatVector& u)
  // No Exceptions
  : name_(NULL), frame_(u.frame_), t_(u.t_)
  {
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  }

// copy from PositionVector constructor
FloatVector::FloatVector(const PositionVector& r)
  // No Exceptions
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  v_[X] = r.v_[X];
  v_[Y] = r.v_[Y];
  v_[Z] = r.v_[Z];
  }

// Destructor
FloatVector::~FloatVector()
  {
    if(name_) delete name_;
  }

//---------------
// Operators
//---------------

// unary -
FloatVector FloatVector::operator-() const // No Exceptions
  {
  FloatVector rv(*this);
  rv.v_[X] = -v_[X];
  rv.v_[Y] = -v_[Y];
  rv.v_[Z] = -v_[Z];
  return(rv);
  }

// +=
FloatVector& FloatVector::operator+=(const FloatVector& r)
  {
  v_[X] += r.v_[X];
  v_[Y] += r.v_[Y];
  v_[Z] += r.v_[Z];
  return(*this);
  }

// +=
FloatVector& FloatVector::operator+=(const Uvar& a)
  {
  v_[X] += a;
  v_[Y] += a;
  v_[Z] += a;
  return(*this);
  }

// -=
FloatVector& FloatVector::operator-=(const FloatVector& r)
  {
  v_[X] -= r.v_[X];
  v_[Y] -= r.v_[Y];
  v_[Z] -= r.v_[Z];
  return(*this);
  }

// -=
FloatVector& FloatVector::operator-=(const Uvar& a)
  {
  v_[X] -= a;
  v_[Y] -= a;
  v_[Z] -= a;
  return(*this);
  }

// *=
FloatVector& FloatVector::operator*=(const FloatVector& r)
  {
  v_[X] *= r.v_[X];
  v_[Y] *= r.v_[Y];
  v_[Z] *= r.v_[Z];
  return(*this);
  }

// *=
FloatVector& FloatVector::operator*=(const Uvar& a)
  {
  v_[X] *= a;
  v_[Y] *= a;
  v_[Z] *= a;
  return(*this);
  }

// /=
FloatVector& FloatVector::operator/=(const FloatVector& r)
  {
  v_[X] /= r.v_[X];
  v_[Y] /= r.v_[Y];
  v_[Z] /= r.v_[Z];
  return(*this);
  }

// /=
FloatVector& FloatVector::operator/=(const Uvar& a)
  {
  v_[X] /= a;
  v_[Y] /= a;
  v_[Z] /= a;
  return(*this);
  }

// Assignment
FloatVector& FloatVector::operator=(const FloatVector& u) // No Exceptions
  {
  if (v_ == u.v_) return(*this);  // cover self assignment
  // Assignment leaves *this's name unchanged, but copies the data.
  frame_ = u.frame_;
  t_ = u.t_;
  v_[X] = u.v_[X];
  v_[Y] = u.v_[Y];
  v_[Z] = u.v_[Z];
  return(*this);
  }

// Indexing

Uvar& FloatVector::operator[](const indexE& i) // No Exceptions
  {
  return(v_[i]);
  }

const Uvar& FloatVector::operator[](const indexE& i)
  const // No Exceptions
  {
  return(v_[i]);
  }

//---------------
// Predicates
//---------------

bool FloatVector::operator==(FloatVector b) const
  {
  b.representIn(frame_);
  return(timeMatch(b) && v_[X] == b.v_[X] &&
    v_[Y] == b.v_[Y] && v_[Z] == b.v_[Z]);
  }

bool FloatVector::operator!=(FloatVector b) const
  {
  return(!(FloatVector::operator==(b)));
  }

bool FloatVector::timeMatch(const FloatVector& b) const
  {
  return(t_ == b.t_);
  }

bool FloatVector::timeMatch(const PositionVector& b) const
  {
  return(t_ == b.t_);
  }

bool FloatVector::timeMatch(const DirectionVector& b) const
  {
  return(t_ == b.t_);
  }

bool FloatVector::timeMatch(const double& t) const
  {
  return(t_ == t);
  }
bool FloatVector::timeMatch(const Time& t) const
  {
    double et;
    t.getEt(et);
    return(t_ == et);
  }

bool FloatVector::frameMatch(const FloatVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool FloatVector::frameMatch(const PositionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

bool FloatVector::frameMatch(const DirectionVector& b) const
  // No Exceptions
  {
  return(frame_ == b.frame_);
  }

//-----------------------------------------------------------------
// isSpice()
//
// Returns true if the frame this FloatVector is represented in
// is a Spice frame.
//-----------------------------------------------------------------

bool FloatVector::isSpice() const // No Exceptions
  {
  return(frame_.isSpice());
  }

//-----------------
// Get/Set methods
//-----------------

//---------------
// Other methods
//---------------

//-------------------------------------
// name_str = name();
//
// Returns name of this FloatVector
//-------------------------------------

string FloatVector::name() const // No Exceptions
  {
   if (name_) return(*name_);
   else return("");
  }

//------------------------------------
// t = time();
//
// Returns time of this FloatVector
//------------------------------------

Time FloatVector::time() const // No Exceptions
  {
    Time t;
    t.setEt(t_);
    return(t);
  }
  

//------------------------------------
// setTime();
//
// Sets time of this FloatVector
//------------------------------------

void FloatVector::setTime(const Time& t) // No Exceptions
  {
    t.getEt(t_);
  }

//------------------------------------
// t = timeInSeconds();
//
// Returns time of this PositionVector in seconds
//------------------------------------

double FloatVector::timeInSeconds() const // No Exceptions
  {
     return(t_);
  }


void FloatVector::setTimeInSeconds(double t) // No Exceptions
  {
  t_ = t;
  }
//------------------------------------------------------------------
// name_str = frameName();
//
// Returns name of the frame this FloatVector is represented in.
//------------------------------------------------------------------

string FloatVector::frameName() const // No Exceptions
  {
  return(frame_.name());
  }

//----------------------------------------------------------------------------
// representIn(frame)
//
// Transform *this FloatVector object to the indicated frame defined at
// the same time.  Note that just like DirectionVectors, FloatVectors do
// not account for origin translation.
//----------------------------------------------------------------------------

FloatVector& FloatVector::representIn(const Frame& frame)
  {
  if (frame_ == frame) return(*this);  // nothing to do

  SpiceDouble et;
  et=t_;

  // xfer to SPICE compatible storage
  SpiceDouble v[3];
  v[X] = v_[X].getValue();
  v[Y] = v_[Y].coerce(v_[X]);
  v[Z] = v_[Z].coerce(v_[X]);

  PositionVector::frameRotate(frame_,frame,et,v);

  v_[X].setValue(v[X]);
  v_[Y].setValue(v[Y]);
  v_[Z].setValue(v[Z]);

  //------------------------------------------------------------------------
  // Set this's frame to the new frame (since this is now represented in the
  // new frame).
  //------------------------------------------------------------------------

  frame_ = frame;

  return(*this);
  }

//----------------------------------------------------------------------------
// representIn(r)
//
// Transform *this FloatVector object to the same frame that r is
// represented in.
// If the times of *this and r do not match, an exception is thrown.
//----------------------------------------------------------------------------

FloatVector& FloatVector::representIn(const PositionVector& r)
  {
  if (t_ != r.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + r.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(r.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this FloatVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

FloatVector& FloatVector::representIn(const DirectionVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//----------------------------------------------------------------------------
// representIn(u)
//
// Transform *this FloatVector object to the same frame that u is
// represented in.
// If the times of *this and u do not match, an exception is thrown.
//----------------------------------------------------------------------------

FloatVector& FloatVector::representIn(const FloatVector& u)
  {
  if (t_ != u.t_)
    {
    GeomError e("Tried to represent " + name() + " in the frame of "
      + u.name() + " defined at a different time",GeomError::time_mismatch);
    e.throwMe();
    }
  return(representIn(u.frame_));
  }

//-------------------------------------------------------
// magnitude()
//
// Return the magnitude of *this vector (including units)
//-------------------------------------------------------

Uvar FloatVector::magnitude() const 
  {
  return(sqrt(v_[X]*v_[X] + v_[Y]*v_[Y] + v_[Z]*v_[Z]));
  }

//-------------------------------------------------
// scale(mag)
//
// Scale *this to have the indicated magnitude.
// Throw an exception if attempting to scale zero.
//-------------------------------------------------

void FloatVector::scale(Uvar mag)
  {
  Uvar cur_mag = magnitude();
  if (cur_mag.getValue() == 0.0)
    {
    GeomError e("Tried to scale zero PositionVector: " + name(),
      GeomError::scale_zero);
    e.throwMe();
    }
  v_[X] *= mag/cur_mag;
  v_[Y] *= mag/cur_mag;
  v_[Z] *= mag/cur_mag;
  }

//-------------------------------------------------------------------------
// angle(v)
//
// Compute the subtended angle between *this and v.
// Result is Uvar with units "rads".
//-------------------------------------------------------------------------

Uvar FloatVector::angle(const FloatVector& v) const
  {
  DirectionVector v1(*this);
  DirectionVector v2(v);
  return(Uvar(acos(dot(v1,v2)),"rad"));
  }

//-------------------------------------------------------------------------
// invert()
//
// Take the reciprocal of each element.
// If any element is zero, an exception is thrown.
//-------------------------------------------------------------------------

void FloatVector::invert()
  {
  if (v_[X] == 0 || v_[Y] == 0 || v_[Z] == 0)
    {
    GeomError e("Attempted to divide by zero element of FloatVector: " +
      name());
    e.throwMe();
    }
  v_[X] = 1.0/v_[X];
  v_[Y] = 1.0/v_[Y];
  v_[Z] = 1.0/v_[Z];
  }

//------------------------------------
// Binary operators for FloatVector
//------------------------------------

// +
FloatVector operator+(FloatVector a1, const FloatVector& a2)
  {
  return(a1 += a2);
  }

// +
FloatVector operator+(FloatVector a1, const Uvar& a2)
  {
  return(a1 += a2);
  }

// +
FloatVector operator+(const Uvar& a1, FloatVector a2)
  {
  return(a2 += a1);
  }

// -
FloatVector operator-(FloatVector a1, const FloatVector& a2)
  {
  return(a1 -= a2);
  }

// -
FloatVector operator-(FloatVector a1, const Uvar& a2)
  {
  return(a1 -= a2);
  }

// -
FloatVector operator-(const Uvar& a1, FloatVector a2)
  {
  a2 = -a2;
  return(a2 += a1);
  }

// *
FloatVector operator*(FloatVector a1, const FloatVector& a2)
  {
  return(a1 *= a2);
  }

// *
FloatVector operator*(FloatVector a1, const Uvar& a2)
  {
  return(a1 *= a2);
  }

// *
FloatVector operator*(const Uvar& a1, FloatVector a2)
  {
  return(a2 *= a1);
  }

// /
FloatVector operator/(FloatVector a1, const FloatVector& a2)
  {
  return(a1 /= a2);
  }

// /
FloatVector operator/(FloatVector a1, const Uvar& a2)
  {
  return(a1 /= a2);
  }

// /
FloatVector operator/(const Uvar& a1, FloatVector a2)
  {
  a2.invert();
  return(a2 *= a1);
  }

// <<
ostream& operator<<(ostream& s, const FloatVector& u)
  {
  s << "[" << u[FloatVector::X] << ","
    << u[FloatVector::Y] << ","
    << u[FloatVector::Z] << "]";
  return(s);
  }

//----------------------
// Supporting functions
//----------------------

//--------------------------------------------------------------------------
// val = dot(v1,v2)
//
// Form the vector dot product of v1 and v2, returned in val.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the dot product.
//--------------------------------------------------------------------------

Uvar dot(const FloatVector& v1, FloatVector v2)
  {
  if (!v1.frameMatch(v2)) v2.representIn(v1);
  return(v1[FloatVector::X]*v2[FloatVector::X] +
         v1[FloatVector::Y]*v2[FloatVector::Y] +
         v1[FloatVector::Z]*v2[FloatVector::Z]);
  }

//--------------------------------------------------------------------------
// v = cross(v1,v2)
//
// Form the vector cross product of v1 and v2, returned in v.
// If v1 and v2 are in different frames, v2 is represented in v1's frame
// before doing the cross product.  Thus, the result will be in v1's frame.
//--------------------------------------------------------------------------

FloatVector cross(const FloatVector& v1, FloatVector v2)
  {
  if (!v1.timeMatch(v2))
    {
    GeomError e("Tried to cross vectors: " + v1.name() + " and "
      + v2.name() + " defined at different times",GeomError::time_mismatch);
    e.throwMe();
    }

  if (!v1.frameMatch(v2)) v2.representIn(v1);
  FloatVector result(v1);

  result[FloatVector::X] = v1[FloatVector::Y] * v2[FloatVector::Z] -
    v1[FloatVector::Z] * v2[FloatVector::Y];
  result[FloatVector::Y] = v1[FloatVector::Z] * v2[FloatVector::X] -
    v1[FloatVector::X] * v2[FloatVector::Z];
  result[FloatVector::Z] = v1[FloatVector::X] * v2[FloatVector::Y] -
    v1[FloatVector::Y] * v2[FloatVector::X];
  return(result);
  }

//----------------------------------
// Methods for class StateVector
//----------------------------------

//--------------
// Constructors
//--------------

// Default constructor sets default frame and time.
StateVector::StateVector()
  // No Exceptions
  : name_(NULL), frame_()
  {
  }

// Default constructor sets default frame and time.
StateVector::StateVector(const string& name)
  // No Exceptions
  : name_(NULL), frame_()
  {
    name_=new string;
    *name_=name;
  }

// Construct from position and velocity vectors
StateVector::StateVector(const string& name, const PositionVector& r,
  const FloatVector& v)
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
    name_=new string;
    *name_=name;
    setState(r,v);
  }

// Construct from position and velocity vectors
StateVector::StateVector(const PositionVector& r, const FloatVector& v)
  : name_(NULL), frame_(r.frame_), t_(r.t_)
  {
  setState(r,v);
  }

// copy constructor
StateVector::StateVector(const StateVector& s)
  // No Exceptions
  : name_(NULL), frame_(s.frame_), t_(s.t_)
  {
  r_[X] = s.r_[X];
  r_[Y] = s.r_[Y];
  r_[Z] = s.r_[Z];
  v_[X] = s.v_[X];
  v_[Y] = s.v_[Y];
  v_[Z] = s.v_[Z];
  }

// Destructor
StateVector::~StateVector()
  {
    if(name_) delete name_;
  }

//-------------
// Operators
//-------------

// unary -
StateVector StateVector::operator-() const // No Exceptions
  {
  StateVector rv(*this);
  rv.r_[X] = -r_[X];
  rv.r_[Y] = -r_[Y];
  rv.r_[Z] = -r_[Z];
  rv.v_[X] = -v_[X];
  rv.v_[Y] = -v_[Y];
  rv.v_[Z] = -v_[Z];
  return(rv);
  }

// +=
StateVector& StateVector::operator+=(const StateVector& s)
  {
  r_[X] += s.r_[X];
  r_[Y] += s.r_[Y];
  r_[Z] += s.r_[Z];
  v_[X] += s.v_[X];
  v_[Y] += s.v_[Y];
  v_[Z] += s.v_[Z];
  return(*this);
  }

// -=
StateVector& StateVector::operator-=(const StateVector& s)
  {
  r_[X] -= s.r_[X];
  r_[Y] -= s.r_[Y];
  r_[Z] -= s.r_[Z];
  v_[X] -= s.v_[X];
  v_[Y] -= s.v_[Y];
  v_[Z] -= s.v_[Z];
  return(*this);
  }

// *=
StateVector& StateVector::operator*=(const StateVector& s)
  {
  r_[X] *= s.r_[X];
  r_[Y] *= s.r_[Y];
  r_[Z] *= s.r_[Z];
  v_[X] *= s.v_[X];
  v_[Y] *= s.v_[Y];
  v_[Z] *= s.v_[Z];
  return(*this);
  }

// *=
StateVector& StateVector::operator*=(const Uvar& a)
  {
  r_[X] *= a;
  r_[Y] *= a;
  r_[Z] *= a;
  v_[X] *= a;
  v_[Y] *= a;
  v_[Z] *= a;
  return(*this);
  }

// /=
StateVector& StateVector::operator/=(const StateVector& s)
  {
  r_[X] /= s.r_[X];
  r_[Y] /= s.r_[Y];
  r_[Z] /= s.r_[Z];
  v_[X] /= s.v_[X];
  v_[Y] /= s.v_[Y];
  v_[Z] /= s.v_[Z];
  return(*this);
  }

// /=
StateVector& StateVector::operator/=(const Uvar& a)
  {
  r_[X] /= a;
  r_[Y] /= a;
  r_[Z] /= a;
  v_[X] /= a;
  v_[Y] /= a;
  v_[Z] /= a;
  return(*this);
  }

// Assignment
StateVector& StateVector::operator=(const StateVector& s)
  {
  if (v_ == s.v_ && r_ == s.r_) return(*this);  // cover self assignment
  // Assignment leaves *this's name unchanged, but copies the data.
  frame_ = s.frame_;
  t_ = s.t_;
  r_[X] = s.r_[X];
  r_[Y] = s.r_[Y];
  r_[Z] = s.r_[Z];
  v_[X] = s.v_[X];
  v_[Y] = s.v_[Y];
  v_[Z] = s.v_[Z];
  return(*this);
  }

//------------
// Predicates
//------------

bool StateVector::timeMatch(const StateVector& b) const
  {
  return(t_ == b.t_);
  }

bool StateVector::timeMatch(const double& t) const
  {
  return(t_ == t);
  }

bool StateVector::timeMatch(const Time& t) const
  {
    double et;
    t.getEt(et);
    return(t_ == et);
  }

//---------------------
// Get and Set methods
//---------------------

PositionVector StateVector::position() const // No Exceptions
  {
  return(PositionVector(frame_,t_,r_[X].km(),r_[Y].km(),r_[Z].km()));
  }

FloatVector StateVector::velocity() const // No Exceptions
  {
  return(FloatVector(frame_,t_,v_[X],v_[Y],v_[Z]));
  }

//----------------------------------------------------------------------------
// setState(r,v)
//
// Set the position part of this StateVector to match r and the velocity to v.
// The frame and time of r will be taken for this StateVector,
// and the associated velocity vector will be transformed to the
// same frame, and the time will be checked for a match.
// Mismatched times or transformation errors will cause an error
// exception.
//----------------------------------------------------------------------------

void StateVector::setState(const PositionVector& r, FloatVector v)
  {
  if (!r.timeMatch(v))
    {
    ErrorMessage e(
      "Time mismatch between position: " + r.name() +
      " and velocity: " + v.name() + " for state vector: " + name());
    e.throwMe();
    }
  if(!r.frameMatch(v))
    {
    ErrorMessage e(
      "Frame mismatch between position: " + r.name() +
      " and velocity: " + v.name() + " for state vector: " + name());
    e.throwMe();
    }
  r_[X] = r.v_[X];
  r_[Y] = r.v_[Y];
  r_[Z] = r.v_[Z];
  v_[X] = v.v_[X];
  v_[Y] = v.v_[Y];
  v_[Z] = v.v_[Z];
  frame_ = r.frame_;
  t_ = r.t_;
  }

//----------------------------------------------------------------------------
// representIn(frame)
//
// Transform *this StateVector object to the indicated frame defined at
// the same time.  This means the the position and velocity vectors
// are each represnted in the specified frame.
//----------------------------------------------------------------------------

/******** THis is commented out because it is WRONG needs a sxform SPICE call
********* since frames may rotate with respect to each other
StateVector& StateVector::representIn(const Frame& frame)
  {

  // --- Erroneous stuff removed to avoid confusion later

  }

  ********* Bad routine commented out
*********************************/
//---------------
// Other methods
//---------------

//-------------------------------------
// name_str = name();
//
// Returns name of this StateVector
//-------------------------------------

string StateVector::name() const // No Exceptions
  {
  if (name_) return(*name_);
  else return("");
  }

//------------------------------------
// t = time();
//
// Returns time of this StateVector
//------------------------------------

Time StateVector::time() const // No Exceptions
  {
   Time t;
   t.setEt(t_);
   return(t);
  }

//------------------------------------
// setTime();
//
// Sets time of this StateVector
//------------------------------------

void StateVector::setTime(const Time& t) // No Exceptions
  {
  t.getEt(t_);;
  }

//------------------------------------
// t = timeInSeconds();
//
// Returns time of this PositionVector in seconds
//------------------------------------

double StateVector::timeInSeconds() const // No Exceptions
  {
     return(t_);
  }


void StateVector::setTimeInSeconds(double t) // No Exceptions
  {
  t_ = t;
  }

//-------------------------------------------------------------------------
// invert()
//
// Take the reciprocal of each element.
// If any element is zero, an exception is thrown.
//-------------------------------------------------------------------------

void StateVector::invert()
  {
  if (r_[X].getValue() == 0 || r_[Y].getValue() == 0 || r_[Z].getValue() == 0
   || v_[X].getValue() == 0 || v_[Y].getValue() == 0 || v_[Z].getValue() == 0)
    {
    GeomError e("Attempted to divide by zero element of StateVector: " +
      name());
    e.throwMe();
    }
  r_[X] = 1.0/r_[X];
  r_[Y] = 1.0/r_[Y];
  r_[Z] = 1.0/r_[Z];
  v_[X] = 1.0/v_[X];
  v_[Y] = 1.0/v_[Y];
  v_[Z] = 1.0/v_[Z];
  }

//------------------------------------
// Binary operators for StateVector
//------------------------------------

// +
StateVector operator+(StateVector a1, const StateVector& a2)
  {
  return(a1 += a2);
  }

// -
StateVector operator-(StateVector a1, const StateVector& a2)
  {
  return(a1 -= a2);
  }

// *
StateVector operator*(StateVector a1, const StateVector& a2)
  {
  return(a1 *= a2);
  }

// *
StateVector operator*(StateVector a1, const Uvar& a2)
  {
  return(a1 *= a2);
  }

// *
StateVector operator*(const Uvar& a1, StateVector a2)
  {
  return(a2 *= a1);
  }

// /
StateVector operator/(StateVector a1, const StateVector& a2)
  {
  return(a1 /= a2);
  }

// /
StateVector operator/(StateVector a1, const Uvar& a2)
  {
  return(a1 /= a2);
  }

// /
StateVector operator/(const Uvar& a1, StateVector a2)
  {
  a2.invert();
  return(a2 *= a1);
  }

// <<
ostream& operator<<(ostream& s, const StateVector& state)
  {
  s << "[" << state.position() << ","
    << state.velocity() << "]";
  return(s);
  }

//-------------------------
// Class Rotation Methods
//-------------------------

//--------------
// Constructors
//--------------

//-----------------------------------------------------------------------------
// Rotation(r1,r2,r3,o1,o2,o3)
//
// Construct a rotation from an ordered set of rotations (r1-3) around the axes
// of a coordinate frame.  The order of rotations is given by o1-3.
// Thus, Rotation(0,10deg,0,2,0,1) performs a rotation of 10 degrees around
// the x axis.
// Rotation signs follow the right hand rule, and o1-3 are zero offset with
// 0 -> x, 1 -> y, 2 -> z.
//-----------------------------------------------------------------------------

Rotation::Rotation(const Uvar& r1, const Uvar& r2, const Uvar& r3,
  unsigned int o1, unsigned int o2, unsigned int o3)
  // No Exceptions
  : r1_(r1), r2_(r2), r3_(r3), o1_(o1), o2_(o2), o3_(o3)
  {
  }

//----------------------------------------------------------------------------
// rotatedAxes(x,y,z)
//
// Set the direction vectors x,y,z that result when *this Rotation is applied
// to the coordinate axes of the frame currently used by x.
// The DirectionVectors will all be expressed in the x frame.
// Only the components are set - the name and time are left unchanged.
//----------------------------------------------------------------------------

void Rotation::rotatedAxes(DirectionVector& x, DirectionVector& y,
  DirectionVector& z) const
  {
  SpiceDouble cr = cos(r1_);
  SpiceDouble sr = sin(r1_);
  SpiceDouble cp = cos(r2_);
  SpiceDouble sp = sin(r2_);
  SpiceDouble cy = cos(r3_);
  SpiceDouble sy = sin(r3_);

  SpiceDouble trans_matrix[3][3];
  SpiceDouble pmatrix[3][3];

  trans_matrix[0][0] = 1;
  trans_matrix[0][1] = 0;
  trans_matrix[0][2] = 0;
  trans_matrix[1][0] = 0;
  trans_matrix[1][1] = 1;
  trans_matrix[1][2] = 0;
  trans_matrix[2][0] = 0;
  trans_matrix[2][1] = 0;
  trans_matrix[2][2] = 1;

  //-------------------------------------------------------
  // Apply the individual rotation transformations in turn
  //-------------------------------------------------------

  for (unsigned int j=0; j <= 2; ++j)
    {
    if (o1_ == j)
      {
      pmatrix[0][0] = 1;
      pmatrix[0][1] = 0;
      pmatrix[0][2] = 0;
      pmatrix[1][0] = 0;
      pmatrix[1][1] = cr;
      pmatrix[1][2] = sr;
      pmatrix[2][0] = 0;
      pmatrix[2][1] = -sr;
      pmatrix[2][2] = cr;
      }
    if (o2_ == j)
      {
      pmatrix[0][0] = cp;
      pmatrix[0][1] = 0;
      pmatrix[0][2] = -sp;
      pmatrix[1][0] = 0;
      pmatrix[1][1] = 1;
      pmatrix[1][2] = 0;
      pmatrix[2][0] = sp;
      pmatrix[2][1] = 0;
      pmatrix[2][2] = cp;
      }
    if (o3_ == j)
      {
      pmatrix[0][0] = cy;
      pmatrix[0][1] = sy;
      pmatrix[0][2] = 0;
      pmatrix[1][0] = -sy;
      pmatrix[1][1] = cy;
      pmatrix[1][2] = 0;
      pmatrix[2][0] = 0;
      pmatrix[2][1] = 0;
      pmatrix[2][2] = 1;
      }

    mxm_c(pmatrix, trans_matrix, trans_matrix);
    if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
    }

  //----------------------------------------------
  // Make sure all the vectors are in the x-frame
  //----------------------------------------------

  y.representIn(x);
  z.representIn(x);

  //---------------------------------------------------
  // Load components from final transformation matrix
  //---------------------------------------------------

  x[DirectionVector::X] = trans_matrix[0][0];
  x[DirectionVector::Y] = trans_matrix[0][1];
  x[DirectionVector::Z] = trans_matrix[0][2];

  y[DirectionVector::X] = trans_matrix[1][0];
  y[DirectionVector::Y] = trans_matrix[1][1];
  y[DirectionVector::Z] = trans_matrix[1][2];

  z[DirectionVector::X] = trans_matrix[2][0];
  z[DirectionVector::Y] = trans_matrix[2][1];
  z[DirectionVector::Z] = trans_matrix[2][2];
  }




//----------------------------------
// Class RotationalVelocity
//----------------------------------
 RotationalVelocity::RotationalVelocity(const Frame& fsc, 
					const Frame& ftarget,
					const Time& t) 
{
  //first need angular velocity vector
  //
  Uvar delta_t = Uvar(10,"s");//10 ticks
  double et;

  //rotation matrix delta_t/2 later
  SpiceDouble ftarget_to_fsc2[3][3];
  Time t2= t+delta_t/2.0;
  t2.getEt(et);
  Frame::rotationMatrix(ftarget,fsc,et,ftarget_to_fsc2);

   

  //rotation matrix delta_t/2 earlier
  SpiceDouble ftarget_to_fsc1[3][3];
  Time t1 = t -delta_t/2.0;
  t1.getEt(et);
  Frame::rotationMatrix(ftarget,fsc, et, ftarget_to_fsc1);
   

  //rotation matrix from ftarget_to_fsc1 to ftarget_to_fsc2
  SpiceDouble rot_matrix[3][3];
  mxmt_c(ftarget_to_fsc2,ftarget_to_fsc1, rot_matrix);
 

  //compute rotation axis and angle: vector rotation
  SpiceDouble* rot_axis;
  rot_axis = new SpiceDouble[3];
  SpiceDouble rot_angle;
  raxisa_c(rot_matrix, rot_axis, &rot_angle);

  //------------------
  //It turns out that raxisa describes
  // a vector transformation.  For
  // cooridnate transformation, rot_axis
  // is needed to be reversed
  // In addition, those rotation axes are those of
  // fsc!!!  As a result, we need to represent
  // rotation axis in target frame
  //------------------
  for(unsigned int i=0;i<3;++i) {
    rot_axis[i]*=-1.0;
    //cout<<"rot axis "<< i << " "<<rot_axis[i]<<endl;
  }


 //--------------------
  //to describe rotational velocity in target frame
  //--------------------
  FloatVector rotation("",fsc,t,  
		       Uvar(rot_angle,"rad")/delta_t * rot_axis[0],
		       Uvar(rot_angle,"rad")/delta_t * rot_axis[1],
		       Uvar(rot_angle,"rad")/delta_t * rot_axis[2]);
  rotation.representIn(ftarget);

  v_[0] = rotation[FloatVector::X];
  v_[1] = rotation[FloatVector::Y];
  v_[2] = rotation[FloatVector::Z];
  t_=t;
  delete[] rot_axis;
}
//-----------------------
//--------------------------------
// Supporting external functions
//--------------------------------

//------------------------------------------------------------
// spice_surfpt(rsurf,found,observer,ulook,radii)
//
// Computes surface intercept point using spice routine.
// The frame of the result will be set to match the target frame
// taken from the input PositionVector target_radii.
// If no intercept is found, found is set to false.
//------------------------------------------------------------

void spice_surfpt(PositionVector& rsurf, bool& found, PositionVector observer,
  DirectionVector ulook, const PositionVector& target_radii)
  {

  // Put observer position and look direction in target frame
  observer.representIn(target_radii);
  ulook.representIn(target_radii);

  //---------------------------------------------------------------
  // Load spice compatible storage and call spice routine surfpt_c
  //---------------------------------------------------------------

  SpiceDouble r[3],u[3],radii_x,radii_y,radii_z,rs[3];
  SpiceBoolean ifound;
  r[PositionVector::X] = observer[PositionVector::X].getInUnits(km_str);
  r[PositionVector::Y] = observer[PositionVector::Y].getInUnits(km_str);
  r[PositionVector::Z] = observer[PositionVector::Z].getInUnits(km_str);
  u[DirectionVector::X] = ulook[DirectionVector::X];
  u[DirectionVector::Y] = ulook[DirectionVector::Y];
  u[DirectionVector::Z] = ulook[DirectionVector::Z];
  radii_x = target_radii[PositionVector::X].getInUnits(km_str);
  radii_y = target_radii[PositionVector::Y].getInUnits(km_str);
  radii_z = target_radii[PositionVector::Z].getInUnits(km_str);
  surfpt_c(r,u,radii_x,radii_y,radii_z,rs,&ifound);
  if(failed_c())
    {
    found=false;
    reset_c();
    }
  if (ifound) found = true; else found = false;

  if (ifound)
    {
    rsurf = PositionVector(target_radii,
      Uvar(rs[PositionVector::X],km_str),
      Uvar(rs[PositionVector::Y],km_str),
      Uvar(rs[PositionVector::Z],km_str));
    }

  }

//------------------------------------------------------------
// spice_ringplane_intercept(rint,found,observer,ulook)
//
// Computes ringplane surface intercept point using spice routine.
// The frame of the result will be set to match the input look vector frame.
// If no intercept is found, found is set to false.
//------------------------------------------------------------

void spice_ringplane_intercept(PositionVector& rint, bool& found,
  PositionVector observer, DirectionVector ulook,
  const PositionVector& target_radii)
  {

  // Put observer position in look frame
  observer.representIn(ulook);

  //---------------------------------------------------------------
  // Load spice compatible storage and call spice routine inrypl_c
  //---------------------------------------------------------------

  int nxpts;
  SpiceDouble zvec[3] = {0,0,1};
  SpiceDouble rsc[3],rs[3],ul[3];
  SpicePlane ringplane;
  cout << "zvec = " << zvec[0] << " " << zvec[1] << " " << zvec[2]  << endl;
  nvc2pl_c(zvec,0.0,&ringplane);
  SpiceBoolean ifound;
  rsc[PositionVector::X] = observer[PositionVector::X].getInUnits(km_str);
  rsc[PositionVector::Y] = observer[PositionVector::Y].getInUnits(km_str);
  rsc[PositionVector::Z] = observer[PositionVector::Z].getInUnits(km_str);
  cout << "rsc = " << rsc[0] << " " << rsc[1] << " " << rsc[2]  << endl;
  cout << "ulook = " << ulook << endl;
  inrypl_c(rsc,ulook.data(),&ringplane,&nxpts,rs);
  cout << "rs = " << rs[0] << " " << rs[1] << " " << rs[2]  << endl;
  if (nxpts != 1)
    {
    found = false;
    }
  else
    {
    found = true;
    }
  if(failed_c())
    {
    found=false;
    reset_c();
    }

  if (ifound)
    {
    if (vnorm_c(rs) < target_radii[PositionVector::X].getInUnits("km"))
      {  // Inside body of target (Saturn), so use regular body intercept.
      spice_surfpt(rint,found,observer,ulook,target_radii);
      }
    else
      {
      rint = PositionVector(observer,
        Uvar(rs[PositionVector::X],km_str),
        Uvar(rs[PositionVector::Y],km_str),
        Uvar(rs[PositionVector::Z],km_str));
      }
    }

  }

//------------------------------------------------------------
// spice_nearpt(rnadir,altitude,observer,radii)
//
// Computes nadir intercept point and altitude using spice routine.
// The frame of the result will be set to match the target frame
// taken from the input PositionVector target_radii.
//------------------------------------------------------------

void spice_nearpt(PositionVector& rnadir, Uvar& alt,
  PositionVector observer, const PositionVector& target_radii)
  {
  // Put observer position in target frame
  observer.representIn(target_radii);

  //---------------------------------------------------------------
  // Load spice compatible storage and call spice routine nearpt_c
  //---------------------------------------------------------------

  SpiceDouble r[3],radii_x,radii_y,radii_z,rs[3],h;
  r[PositionVector::X] = observer[PositionVector::X].getInUnits(km_str);
  r[PositionVector::Y] = observer[PositionVector::Y].getInUnits(km_str);
  r[PositionVector::Z] = observer[PositionVector::Z].getInUnits(km_str);
  radii_x = target_radii[PositionVector::X].getInUnits(km_str);
  radii_y = target_radii[PositionVector::Y].getInUnits(km_str);
  radii_z = target_radii[PositionVector::Z].getInUnits(km_str);
  nearpt_c(r,radii_x,radii_y,radii_z,rs,&h);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}

  rnadir = PositionVector(target_radii,
    Uvar(rs[PositionVector::X],km_str),
    Uvar(rs[PositionVector::Y],km_str),
    Uvar(rs[PositionVector::Z],km_str));
  alt = Uvar(h,km_str);
  }

//------------------------------------------------------------
// spice_npedln(radii,rpoint,udir,rnear,dist)
//
// Computes nearest point on a line to an ellipsoid using the spice routine
// npedln.
// Inputs:
//   rx,ry,rz are the body axes of the ellipsoid
//   rpoint is a point on the line
//   udir points along the line
// Outputs:
//   rnear is the nearest point on the line.
//   dist is the distance from the line to the nearest point.
//------------------------------------------------------------

void spice_npedln(const Uvar& rx, const Uvar& ry, const Uvar& rz,
  const PositionVector& rpoint, DirectionVector udir,
  PositionVector& rnear, Uvar& dist)
  {
  // Put point and direction in same frame
  udir.representIn(rpoint);

  //---------------------------------------------------------------
  // Load spice compatible storage and call spice routine nearpt_c
  //---------------------------------------------------------------

  SpiceDouble r[3],radii_x,radii_y,radii_z,u[3],pnear[3],d;
  r[0] = rpoint[PositionVector::X].getInUnits(km_str);
  r[1] = rpoint[PositionVector::Y].getInUnits(km_str);
  r[2] = rpoint[PositionVector::Z].getInUnits(km_str);
  u[0] = udir[DirectionVector::X];
  u[1] = udir[DirectionVector::Y];
  u[2] = udir[DirectionVector::Z];
  radii_x = rx.getInUnits(km_str);
  radii_y = ry.getInUnits(km_str);
  radii_z = rz.getInUnits(km_str);
  npedln_c(radii_x,radii_y,radii_z,r,u,pnear,&d);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}

  rnear = PositionVector(rpoint,
    Uvar(pnear[0],km_str),
    Uvar(pnear[1],km_str),
    Uvar(pnear[2],km_str));
  dist = Uvar(d,km_str);
  }

//------------------------------------------------------------------
// spice_target_value(target_name,value_name,value)
//
// Fetch a value (identified by value_name) for a specific target
// (identified by target_name), and return that single value.
// If the value_name corresponds to a multi-component return value,
// then an error exception is generated.  Unrecognized names will
// also cause error exceptions.  Some specially checked value names
// will be given appropriate units.  All others are assumed to have
// no units.
//------------------------------------------------------------------

void spice_target_value(const string& name, const string& value_name,
  Uvar& value)
  {
  SpiceInt id;
  spice_target_id(name,id);
  SpiceInt dim;
  SpiceDouble values[10];
  bodvar_c(id,value_name.c_str(),&dim,values);
  if (failed_c())
    {
    GeomError e(GeomError::spice_misc);
    e.throwMe();
    }
  else if (dim != 1)
    {
    GeomError e("spice_target_value(" + name + "," + value_name +
      "): Multiple values for single return", GeomError::spice_misc);
    e.throwMe();
    }

  if (value_name == "GM")
    {
    value = Uvar(values[0],"km km km/(s s)");
    }
  else
    {
    value = Uvar(values[0]);
    }
  }

void spice_frame_name(SpiceInt id, SpiceChar* frname, int length){
  frmnam_c(id,length,frname);
  if(frname[0]==0){
    ErrorMessage e("No Spice Frame with ID "+toStr(id));
    e.throwMe();
  }
}

void spice_frame_name(SpiceInt id, string& name){
  SpiceChar frname[33];
  frmnam_c(id,33,frname);
  name=frname;
  if(frname[0]==0){
    ErrorMessage e("No Spice Frame with ID "+toStr(id));
    e.throwMe();
  }
}

void spice_frame_id(const string& name, SpiceInt& id)
{
  namfrm_c(name.c_str(),&id);
  if (id==0){
    ErrorMessage e("No Spice frame called: " + name);
    e.throwMe();
  } 
}

void spice_target_id(const string& name, SpiceInt& id)
{
  SpiceBoolean ifound;
  bodn2c_c(name.c_str(), &id, &ifound);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  if (ifound == SPICEFALSE)
    {
    GeomError e("spice_target_id: No Spice target called: " + name,
      GeomError::spice_bad_target);
    e.throwMe();
    }
}

void spice_target_radii(const string& name, const Frame& frame, const Time& t,
			SpiceInt id, PositionVector& radii)
  {
  SpiceInt Nvals;
  SpiceDouble vals[3];
  bodvar_c(id, "RADII", &Nvals, vals);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  if (Nvals != 3)
    {
    GeomError e("spice_target_radii: Can't get radii for: " + name,
      GeomError::spice_bad_target_radii);
    e.throwMe();
    }
    radii = PositionVector("target_radii", frame, t,
    Uvar(vals[0],km_str), Uvar(vals[1],km_str), Uvar(vals[2],km_str) );

}

void spice_target_radii(const string& name, Uvec& radii)
  {
  SpiceInt id;
  spice_target_id(name,id);
  SpiceInt Nvals;
  SpiceDouble vals[3];
  bodvar_c(id, "RADII", &Nvals, vals);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  if (Nvals != 3)
    {
    GeomError e("spice_target_radii: Can't get radii for: " + name,
      GeomError::spice_bad_target_radii);
    e.throwMe();
    }
  radii.resize(3);
  radii(0) = Uvar(vals[0],km_str);
  radii(1) = Uvar(vals[1],km_str);
  radii(2) = Uvar(vals[2],km_str);
  }

Uvar spice_utc_to_et(const string& utc_str){
  SpiceDouble ephemtime;
  utc2et_c(utc_str.c_str(),&ephemtime);
  if(failed_c()){
    GeomError e(GeomError::spice_misc); 
    e.throwMe();
  }
  return(Uvar(ephemtime,seconds_str));
}

bool isUtc(const string& utc_str)
  {
  SpiceDouble et;
  utc2et_c(utc_str.c_str(),&et);
  return(!failed_c());
  }

Uvar spice_sclk_to_et(SpiceInt sc_id,unsigned int sclk){
  typedef std::ostringstream OSTRINGSTREAM;
  OSTRINGSTREAM os;
  os << sclk;
  SpiceDouble ephemtime;
  scs2e_c(sc_id,toStr(os).c_str(),&ephemtime);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  return(Uvar(ephemtime,seconds_str));
}

string spice_et_to_utc(const Uvar& et,const string& format)
{
  char str[27];
  et2utc_c(et.getInUnits(seconds_str),format.c_str(),3,27,str);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  string utc_str(str);
  return(utc_str);
}

double spice_et_to_encoded_sclk(SpiceInt sc_id, const Uvar& et)
  {
  double esclk;
  sce2c_c(sc_id,et.getInUnits(seconds_str),&esclk);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  return(esclk);
  }

Time spice_encoded_sclk_to_et(SpiceInt sc_id, double esclk)
  {
  double et;
  sct2e_c(sc_id,esclk,&et);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  return(Time(Uvar(et,seconds_str)));
  }

unsigned int spice_decode_sclk(SpiceInt sc_id, double esclk){
  SpiceChar sclk_string[20];
  scdecd_c(sc_id,esclk,30,sclk_string);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}

  // remove partition number and convert to int
  char string2[40];
  int p;
  sscanf(sclk_string,"%d/%s\n",&p,string2);
  unsigned int sclk=atol(string2);

  return(sclk);
}

void spice_write_ckernel(const string& filename, unsigned int Ncomchar, 
			 SpiceDouble begtime,
			 SpiceDouble endtime, SpiceInt inst_id, 
			 string ref_string, SpiceBoolean ang_vel_flag, 
			 string seg_id_string, int Nrecords, 
			 SpiceDouble spice_sclkdp[], 
			 SpiceDouble spice_quats[][4],
			 SpiceDouble spice_avvs[][3], 
			 SpiceDouble spice_start_x_time[]){

  SpiceInt nint=1; // by default
  SpiceInt ckhan;
  ckopn_c(filename.c_str(),filename.c_str(),Ncomchar,&ckhan);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  ckw03_c(ckhan,begtime,endtime,inst_id,ref_string.c_str(),ang_vel_flag,
	  seg_id_string.c_str(),Nrecords,spice_sclkdp,spice_quats,spice_avvs,
	  nint,spice_start_x_time);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
  ckcls_c(ckhan);
  if(failed_c()){GeomError e(GeomError::spice_misc); e.throwMe();}
}

//----------------------------------------------//
// GeomError routines                           //
//----------------------------------------------//

bool GeomError::spice_handling_enabled_=true;

 GeomError::GeomError(GeomError::errorE err_type, double et) 
   // No exceptions
   : error_type(err_type),
     error_time(et)
    {
    if (error_type == unspecified)
      msg = "Unspecified Geometry Error";
    else if (error_type == scale_zero)
      msg = "Geometry Error: Tried to scale a zero vector";
    else if (error_type == time_mismatch)
      msg = "Geometry Error: Times do not match";
    else if (error_type == spice_nosolution)
      msg = "Geometry Error: Spice routine failed to find a solution";
    else if (error_type == read_error)
      msg = "Geometry Error: read error";
    else if (error_type == write_error)
      msg = "Geometry Error: write error";
    else if (error_type == ambiguity_error)
      msg = "Geometry Error: Ambiguous solution";
    else if (error_type == spice_misc){
      msg = "Spice Error: ";
      initSpiceError();
    }
    else if (error_type == spice_bad_ckernel)
      msg = "Spice Error: bad Ckernel";
    else if (error_type == spice_bad_target)
      msg = "Spice Error: bad Target";
    else if (error_type == spice_bad_target_radii)
      msg = "Spice Error: bad Target Radii";

    }

  GeomError::GeomError(const string& emsg, 
    GeomError::errorE err_type, double et) 
    //No exceptions
    : error_type(err_type),
      error_time(et)
    {
    msg = emsg;
    if(error_type == spice_misc)
      initSpiceError();
    }

  void GeomError::initSpiceError()
    {
      if(!spice_handling_enabled_) return; // avoid infinite recursion
      // Append Spice Error Info to Message String
      SpiceChar longmsg[SPICEMSGLEN];
      SpiceChar modname[SPICEMODNAMELEN];
      getmsg_c("long",SPICEMSGLEN,longmsg);
      msg=msg+" "+longmsg+"\n"+"Trace:";
     
      SpiceInt depth;
      trcdep_(&depth);
      for (SpiceInt i=1; i<=depth; i++)
	{
	  trcnam_( &i , modname , SPICEMODNAMELEN);
          msg=msg+modname;
          if(i<depth) msg=msg+"-->";
	}



      // reset spice error flag
      reset_c();
      spice_handling_enabled_=false; // avoid infinite recursion
      // check for bad ckernel
      if(error_time!=0){
        try{
	  Time t=Time(Uvar(error_time,seconds_str));
	  Frame j2000("J2000","Earth");
	  Frame cassini_inertial("J2000",cassini_str);
	  // code to make sure error is localized to ckernel
	  PositionVector sc_pos("sc_pos",cassini_inertial,t,0,0,0);
	  sc_pos.representIn(j2000);
	}
	catch(GeomError e){
          // Not a ckernel error so return
	  spice_handling_enabled_=true;
	  reset_c();
	  return;
	}
	SpiceInt j2000_id=1, fsc_id=-82000;
        SpiceDouble rotmat[9];        
        refchg_(&j2000_id,&fsc_id,&error_time,rotmat);
       
        if(failed_c()){
	  error_type=spice_bad_ckernel;
	}
        
      }
      spice_handling_enabled_=true;
      reset_c();
    }
  void GeomError::throwMe()
    {
      throw *this;
    }
