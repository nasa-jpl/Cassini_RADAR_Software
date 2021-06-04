//----------------------------------------------------------------------------
// Utils.cpp
//
// This file contains utility function definitions.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_utils_c[] =
  "@(#) $Id: Utils.cpp,v 11.7 2017/03/22 21:45:42 richw Exp $";

//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//---------------
// Other includes
//---------------

#include <stdlib.h>
#include <strings.h>
#include <string>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "Frame.h"
#include "Units.h"
#include "Time.h"
#include "Constants.h"
#include "TemplateUtils.h"
#include "math.h"
#include "SARFunctions.h"
#include "TargetGeom.h"
#include "Plot.h"
#include "SimpleArray.h"
#include "Utils.h"
using std::cout;
using std::endl;
using std::cerr;

//---------------------------------------------------------------------------
// set_unit_option(unit_option)
//
// Set the behavior of UnitVar's based on the string unit_option.
// unit_option is a command line option argument and can be one of the
// following: auto, base, none.
//   auto - set fully automatic mode
//   base - set enforce_base_units mode
//   none - set no_unit_support mode.
// If unit_option is NULL, then base is selected by default.
// An invalid unit option will cause an error exception.
//---------------------------------------------------------------------------

void set_unit_option(char* unit_option)
  {
  if (unit_option == NULL)
    {
    UnitVar::setMode("enforce_base_units");  // default
    }
  else
    {
    string unit_mode(unit_option);
    if (unit_mode == "auto")
      UnitVar::setMode("automatic");
    else if (unit_mode == "base")
      UnitVar::setMode("enforce_base_units");
    else if (unit_mode == "none")
      UnitVar::setMode("no_unit_support");
    else
      {
      ErrorMessage e("Invalid unit option: " + unit_mode);
      e.throwMe();
      }
    }
  }

//---------------------------------------------------------------------------
// str = toStr(arg)
//
// Convert various types to a std::string.
//---------------------------------------------------------------------------

string toStr(const unsigned int& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(const int& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(const unsigned long& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(const long& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(const double& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(const float& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

string toStr(string& arg)
  {
  ostringstream ost;
  ost << arg;
  return(ost.str());
  }

//---------------------------------------------------------------------------
// str = toStr(ost)
//
// Convert the string contents of an output strstream to a std::string.
// This is needed because the .str() method of strstream does not attach the
// 0 (end of string marker) in the string buffer at the end of the string.
//---------------------------------------------------------------------------

string toStr(ostringstream& ost)
  {
#ifdef __STRSTREAM__
  return(string(ost.str(),ost.pcount()));
#else
  return(ost.str());
#endif
  }

//---------------------------------------------------------------------------
// str = toStr(arg,min_width)
//
// Convert the integer argument to a string with specified minimum width
// and zero fill when needed.
//---------------------------------------------------------------------------

string toStr(unsigned int arg, unsigned int min_width)
  {    
  ostringstream ost;
  ost.width(min_width);
  ost.fill('0');
  ost << arg;
  // Following conditional compilation is needed to allow gcc-2.95 compilers
  // to successfully replace stringstream with strstream.
#ifdef __STRSTREAM__
  return(string(ost.str(),ost.pcount()));
#else
  return(ost.str());
#endif
  } 

//---------------------------------------------------------------------
// string t =get_next_token(string& s, const string& w)
// String parsing routine. Returns first encountered sequence of characters
// in s separated by characters in w. (I.e. w=" \t\n" for whitespace) 
// resets s to start of beginning of next sequence
//----------------------------------------------------------------------

string get_next_token(string& s, const string& w){
  string t="";
  bool w_char_found=true;
  unsigned int s_idx=0;
  //skip initial separator characters
  while(w_char_found && s_idx<s.size())
    {
      w_char_found=false;
      for(unsigned int c=0;c<w.size();c++){
	if (w[c]==s[s_idx]){
	  w_char_found = true;
	  break;
	}
      }
      if(w_char_found) s_idx++;
    }
  //append sequence to t
  while(!w_char_found && s_idx<s.size())
    {
      w_char_found=false;
      for(unsigned int c=0;c<w.size();c++){
	if (w[c]==s[s_idx]){
	  w_char_found = true;
	  break;
	}
      }
      if(!w_char_found){
	t.append(1,s[s_idx]);
	s_idx++;
      }
    }
  //find index of start of next sequence
  while(w_char_found && s_idx < s.size())
    {
     w_char_found=false;
      for(unsigned int c=0;c<w.size();c++){
	if (w[c]==s[s_idx]){
	  w_char_found = true;
	  break;
	}
      }
      if(w_char_found) s_idx++;
    }
  // reset s and return t
  if(s_idx>s.size()) s="";
  else s=&s[s_idx];
  return(t);
}

FILE* fopen(const string& name, const string& mode){
  FILE* retval= fopen(name.c_str(),mode.c_str());
  if(retval==NULL){
    ErrorMessage e("Error: Cannot open file "+name+ " with mode "+mode);
    e.throwMe();
  }
  return(retval);
}

//---------------------------------------------------------------------------
// fileExists(const string& filename)
//
// Check for file existence using stat(2).  This method returns true if
// stat() succeeds.  If stat() fails with error ENOENT and the filename is
// not an empty string, then the file is presumed not to exist and the
// method returns false.  If stat() fails under any other circumstance,
// fileExists() throws an error.  Refer to the man page for stat() for the
// specific failure modes; in general, the same failures will also prevent the
// file from later being opened for reading or writing.
//---------------------------------------------------------------------------
bool fileExists(const string& filename)
{
  struct stat buf;
  if (filename.length() == 0)
  {
    // Filename is an empty string
    ErrorMessage e("Utils::fileExists:  specified filename is an empty string");
    e.throwMe();
  }
  if (stat(filename.c_str(), &buf) == -1)
  {
    if (errno == ENOENT)
    {
      // File does not exist
      return false;
    }
    else
    {
      // Some other error occurred while trying to get the file's status
      string msg = "Utils::fileExists:  stat() failed for " + filename;
      ErrorMessage e(msg);
      e.throwMe();
    }
  }
  // stat() succeeded; file exists
  return true;
}

//----------------------------------------------------------------------
// void getlookDirection(DirectionVector& look, 
//		      unsigned int& solution,
//		      const Frame& ftitan,
//		      const StateVector& sc_state,		      
//		      const Uvar& range, 
//		      const Uvar& doppler,
//		      const Uvar& lambda,
//		      const Uvar& radius) throw (Unit::UnitError,ErrorMessage)
//
//    DirectionVector look
//    unsigned int solution: 0 - no solution,  1- solution
//    ftitan
//    sc_state
//    range : range
//    doppler: doppler
//    lambda: 
//    radius: target radius
//  For given range and doppler, it will determine the look direction vector 
//  based on sc_position and velocity direction vector represented in beam 
//  frame
//----------------------------------------------------------------------------

void getlookDirection(DirectionVector& look,
		      unsigned int& solution,
		      const PositionVector& sc_pos,
		      const FloatVector& sc_vel,      
		      const Uvar& range, 
		      const Uvar& doppler,
		      const Uvar& lambda,
		      const Uvar& radius) throw (Unit::UnitError,ErrorMessage)
  {  
    

     
    double p1,p2, p3, v1,v2,v3;
    DirectionVector dir_pos,dir_vel;
    Time t =sc_pos.time();
    solution = 0; 
    int sign = 1;
    //______________________________
    //Determine direction vectors of
    //sc's position and velocity
    // in beam frame
    //______________________________    
    Frame fbeam("CASSINI_RADAR_" + toStr(3),"Cassini");  
    dir_pos = sc_pos;
    dir_pos.representIn(fbeam);
     
    dir_vel = sc_vel;
    dir_vel.representIn(fbeam);

    p1 = dir_pos[DirectionVector::X];
    p2 = dir_pos[DirectionVector::Y];
    p3 = dir_pos[DirectionVector::Z];
       
    v1 = dir_vel[DirectionVector::X];
    v2 = dir_vel[DirectionVector::Y];
    v3 = dir_vel[DirectionVector::Z];

    //_______________________________________
    //Determine roll direction from direction
    //of sc_state in beam frame
    //_______________________________________
    if (p2 >= 0.0) 
    {
    //roll_direction is "right";
    sign = 1;
    }
    else  
    {
    //roll_direction is "left";
    sign = -1;
    }
  
   
    //__________________________________
    //determine unitless dist1 and freq1
    //__________________________________
    Uvar  dist1 =(radius*radius - sc_pos.magnitude()*sc_pos.magnitude()-range*range)/(2.0 * range *sc_pos.magnitude());
    Uvar freq1 =  doppler * lambda / (2.0 * sc_vel.magnitude());	 	  
    //_______________________________________
    //get magnitude of each values
    //_______________________________________
    double dist = dist1.getValue();
    double freq = freq1.getValue();

    
    //-------------------------------------------------------------------
    //Equations to find a look direction (b1,b2,b3)
    //freq = vx * sin(azi) + vy*cos(azi)*sin(elev) + vz*cos(azi)*cos(elev)
    //dist=  px * sin(azi) + py*cos(azi)*sin(elev) + pz*cos(azi)*cos(elev)
    //--------------------------------------------------------------------
    //simplified solution
    // v1 * b1 + v2 * b2 + v3 * b3 = v_b = freq
    // p1 * b1 + p2 * b2 + p3 * b3 = p_b = dist
    // b1^2 + b2^2 + b3^2 = 1.0
    //usually v1 and p3 are close to +-1
    // b1 = (freq - v2 * b2 - v3*b3)/v1
    // p1 * (freq - v2 * b2 - v3*b3)/v1 + p2 * b2 + p3 * b3 = p_b = dist
    // b2* (p2 - p1 * v2/v1) + b3*(p3 -p1 *v3/v1)= dist - p1*freq/v1
    //*************** Danger: p1 is very small number
    // A * b2 + B * b3 = C
    //************* C is virtually independent of freq because p1 is small
    // b3 = (C - A *b2)/B = D + E*b2
    // b1 = F * b2 + G 
    // v1 * b1 + v2 * b2 + v3 *(D + E*b2) = freq
    // v1* b1 = -(v2 + v3*E) b2 +freq - v3 *D
    // b1 = -(v2/v1 + v3*E/v1) + (freq - v3*D)/v1
    // b1 = F * b2 + G 
    //************* F is very small 
    // b1^2 + b2^2 + b3^2 = 1
    // (F*b2 + G)^2 + b2^2 + (D+E*b2)^2 = 1.0
    //(F^2 + 1 + E^2)* b2^2 + 2*b2 *(FG+DE) + G^2+D^2-1 =0
    //H*b^2 + 2I *b2 + J = 0
    // b2 =[ -I +(-)sqrt(I*I - H*J)]/H	
	
    double D,E,F,G,H,I,J;
    double sol,ans;
    double b1, b2, b3;

    //A = p2 - p1 * v2/v1;
    //B = p3 -p1 *v3/v1 ;
    //C= dist - p1*freq/v1;


    D = (dist - p1*freq/v1) /(p3 -p1 *v3/v1) ;// C/B 
    E = -(p2 - p1 * v2/v1)/ (p3 -p1 *v3/v1)  ;//-A/B;
    F = -v2/v1 - v3*E/v1;
    G = freq/v1 - v3*D/v1;    
    H = F*F + 1.0 + E*E;
    I = F*G + D*E;
    J = G*G + D*D - 1.0;    
    sol = I*I - H * J;
    if (sol >= 0.0)
      {      
      ans = (-I + double(sign)*sqrt(sol))/H;
      b2 = ans;
      b1 = F*b2 + G;
      b3 = D + E*b2;    
     
      //___________________________________________
      //Set the look direction vector
      //Need better way of defining look vector 
      //___________________________________________
      DirectionVector look_("look",fbeam,t,b1,b2,b3);
      look = look_;
      solution = 1;
      }    
    else
      {
	look=DirectionVector();
      }
  }

//---------------------------
//simple calculation of number of looks
//----------------------------
float number_of_looks(const StateVector& sc_state,
		      const Uvar& range,
		      const Uvar& doppler,
		      const Uvar& bpd,
		      const Uvar& lambda,
		      const Uvar& pri,
		      const Uvar& process_bandwidth,
		      const unsigned int& number_of_pulses )
 {
   float looks;
   double cos_theta = 
     (doppler * lambda /(2.0 * sc_state.velocity().magnitude())).getInUnits("");
   if (cos_theta > 1.0)
     {
     cout << "number_of_looks:  Error, cos_theta = " << cos_theta
         << " greater than 1.  Looks cannot be computed, set to 1" << endl;
     return(1);
     //ErrorMessage e("Utils.cpp:: doppler/(2v/lambda) is larger than 1");
     //e.throwMe();
     }
   Uvar ML_freq_bp = 2.0 * sc_state.velocity().magnitude() 
     * sc_state.velocity().magnitude()*bpd/(lambda*range);
   ML_freq_bp *= (1.0-cos_theta * cos_theta);
   Uvar ML_freq_on = 1/(pri *double(number_of_pulses));

  //frequency resolution=1/(integration time(=pri * N_p))
   //cout<<"frequency shift "<<ML_freq_bp<<endl;
  //-------------- SH's simplied model produces the same result----------
  //double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  //Uvar ML_freq_bp = 
  //2.0 * velocity.magnitude()* Vst_bore * bpd*sin_theta/(lambda*range_bore);
  //note: Vst_bore = V * (1-cos(theta)^2)^0.5 --> same result!!!!!
  //for detailed information, see c:/ygim/radar_work/iso_doppler_profile
  //-------------------------------------------------------------------------

  //unsigned int ML_i_fon = 
  //(unsigned int) (ML_freq_on / dop_interval).getInUnits("");
  
  looks = float(  (process_bandwidth/ML_freq_bp).getInUnits(""));
  //debug
  //cout<<"number of looks "<< looks<<endl;
  return(looks);
 }
//----------------------
//based on ambiguity speed up algoritm
//-----------------------

void compute_all_geometry(const Time& t,
			  const Frame& ftarget,
			  const Frame& fbeam,
			  const Frame& ftrack,
			  const StateVector& state,
			  const Uvar& radius,
			  const Uvec& range_grid,
			  const Uvec& doppler_grid,
			  Umat& alongtracks,
			  Umat& crosstracks,
			  Array2D<DirectionVector>& look_directions,
			  Umat& incidence_angles,
			  Array2D<PositionVector>& positions,
			  Array2D<unsigned int>& no_of_solution_grid)
  {
    unsigned int Nrange = range_grid.size();
    unsigned int Ndoppler = doppler_grid.size();
    //check time
  
    //if(t!=Uvar(state.time(),"s")) ErrorMessage("Utils.cpp:compute_all_geometry: time mismatch").throwMe();
    if(t!=state.time()) ErrorMessage("Utils.cpp:compute_all_geometry: time mismatch").throwMe();
    //recycled variables
    double a,b,c;
    DirectionVector target_to_beam_x,target_to_beam_y,target_to_beam_z;
    DirectionVector target_to_track_x,target_to_track_y,target_to_track_z;
   
    //transformation 
    fbeam.axialVectors(ftarget,
		       t,
		       target_to_beam_x,
		       target_to_beam_y,
		       target_to_beam_z);
    ftrack.axialVectors(ftarget,
			t,
			target_to_track_x,
			target_to_track_y,
			target_to_track_z);

    //recycled look vectors
    DirectionVector lookInBeamFrame("beam frame direction",fbeam,t,0,0,1);
    DirectionVector dirToSurface("track frame direction",ftarget,t,0,0,1);
   
    //radius in km
    double radius_in_km = radius.getInUnits("km");

    //-------------- spacecraft state vector------
    DirectionVector pos = state.position();
    DirectionVector vel = state.velocity();
    DirectionVector dir_pos = state.position();
    DirectionVector dir_vel = state.velocity();
    DirectionVector dir_cross = cross(dir_vel,dir_pos); 
    DirectionVector ulook(" from sc to surface",ftarget,t,0,0,1);

    //-----------------------------------------------
    //determine roll direction 
    //------------------------------------------------
    DirectionVector bore("boresight",fbeam,t,0,0,1);
    double roll_direction=dot(bore,dir_cross);

    //---------------dummy variables to be used for solving RangeDoppler
    double p1,p2,p3, v1,v2,v3;
    double sol;
    int sign=1;
    
    p1 = pos[DirectionVector::X];
    p2 = pos[DirectionVector::Y];
    p3 = pos[DirectionVector::Z];
    
    v1 = vel[DirectionVector::X];
    v2 = vel[DirectionVector::Y];
    v3 = vel[DirectionVector::Z];
    
 

    double u1,u2,u3;
    double A,B,C,D,E,F;
    double range_in_radius;
    double position_in_radius = (state.position().magnitude()/radius).getInUnits("");
    double compute_roll;
    Uvar lambda = speed_light/carrier_frequency;
 
    double* d_range= new double[range_grid.size()];
    double* d_doppler = new double[doppler_grid.size()];
    for (unsigned int i_range = 0; i_range<Nrange; ++i_range)
      d_range[i_range]=range_grid(i_range).getInUnits("km");

    for (unsigned int j_dop =0; j_dop<Ndoppler; ++j_dop)
      d_doppler[j_dop]=doppler_grid(j_dop).getInUnits("Hz");
    
    double d_radius = radius.getInUnits("km");
    double flatform_speed = state.velocity().magnitude().getInUnits("km/s");
    double d_lambda = lambda.getInUnits("km");


    for (unsigned int i_range = 0; i_range<Nrange; ++i_range)
    for (unsigned int j_dop =0; j_dop<Ndoppler; ++j_dop)
      {	
	//------------- Instead of relying on calling targetgeom, 
	//let's solve everyting here
	//-----------------------------------------------
	range_in_radius = d_range[i_range]/d_radius;
	A =  (1 - position_in_radius*position_in_radius 
	      - range_in_radius*range_in_radius)
	  /(2.0 * range_in_radius * position_in_radius);
	B =  d_lambda * d_doppler[j_dop]/(2.0*flatform_speed);
	
	C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
	D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
	E = (A - p1 *C)/p3;
	F = -(p2 + p1*D)/p3;
	
	sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);
	if (sol >= 0)
	  {
	    no_of_solution_grid(i_range,j_dop) = 1;   	 
	    u2 = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
	    u1 = C + D*u2;
	    u3 = E + F*u2;
	    
	    //set look vector from sc to surface in target frame
	    ulook[DirectionVector::X] = u1;
	    ulook[DirectionVector::Y] = u2;
	    ulook[DirectionVector::Z] = u3;
	    
	    //--------------------------------
	    //set surface intercept point
	    //--------------------------------
	    positions(i_range,j_dop) = state.position();
	    positions(i_range,j_dop) += range_grid(i_range)*ulook;
     
	    //-------------------------------------------
	    //Find directly the along/cross distances
	    //by using axialVectors between ftrack and ftitan
	    // surface_intercept_grid is defined in target frame
	    // in order to calculate along/cross, we need to transform 
	    //the surface vector
	    // in ftrack frame
	    //-------------------------------------------   
	    
	    dirToSurface = positions(i_range,j_dop);//dir vector in target frame
	    compute_roll=dir_cross[DirectionVector::X]*dirToSurface[DirectionVector::X]
	      +dir_cross[DirectionVector::Y]*dirToSurface[DirectionVector::Y]
	      +dir_cross[DirectionVector::Z]*dirToSurface[DirectionVector::Z];
	    
	    if (compute_roll*roll_direction > 0)
	      {
		sign = 1;
	      }
	    else
	      {
		sign=-1;
		u2 = ( -(C*D + E*F)+sqrt( sol))/(D*D + F*F + 1.0);
		u1 = C + D*u2;
		u3 = E + F*u2;
		
		//set look vector from sc to surface in target frame
		ulook[DirectionVector::X] = u1;
		ulook[DirectionVector::Y] = u2;
		ulook[DirectionVector::Z] = u3;
		
		//--------------------------------
		//set surface intercept point
		//--------------------------------
		positions(i_range,j_dop) = state.position();
		positions(i_range,j_dop) += range_grid(i_range)*ulook;
		dirToSurface = positions(i_range,j_dop);//dir vector in target frame
	      }
      	    
	    a = dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::X]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::X]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::X];
	    
	    b= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Y]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Y];
	    
	    c= dirToSurface[DirectionVector::X]*target_to_track_x[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Y]*target_to_track_y[DirectionVector::Z]
	      +dirToSurface[DirectionVector::Z]*target_to_track_z[DirectionVector::Z];
	    
	    //---------------------------------------------
	    //theta = -acos(c)+pi/2;//change to latitude
	    //phi = atan2(b,a);
	    //-------------------------------------------------
	    crosstracks(i_range,j_dop).setValue(radius_in_km * (pi/2.0-acos(c)));
	    alongtracks(i_range,j_dop).setValue(-radius_in_km *atan2(b,a));
	    
	    //----------------------------------------------------------------
	    //if (radius_in_km * (pi/2.0-acos(c)) == 0 ||
	    //  radius_in_km*atan2(b,a) == 0) 
	    //ErrorMessage("surface solution, but along/cross are 0").throwMe();
	    //---------------------------------------------------
	    //direction vector look is defined in target frame
	    // need to transform into beam frame
	    //reuse double variables a,b, and c
	    //ulook is defined in target frame
	    // so need a rotation to beam frame
	    //---------------------------------------------------
	    a = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::X]
	      +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::X]
	      +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::X];
	    
	    b = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Y]
	      +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Y]
	      +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Y];
	    
	    c = ulook[DirectionVector::X] * target_to_beam_x[DirectionVector::Z]
	      +ulook[DirectionVector::Y] * target_to_beam_y[DirectionVector::Z]
	      +ulook[DirectionVector::Z] * target_to_beam_z[DirectionVector::Z];
	    lookInBeamFrame[DirectionVector::X] = a;//in beam frame
	    lookInBeamFrame[DirectionVector::Y] = b;//in beam frame
	    lookInBeamFrame[DirectionVector::Z] = c;//in beam frame
	    
	    look_directions(i_range,j_dop)=lookInBeamFrame;
	    //----------------------
	    //set incidence angle
	    //------------------------
	    a = dot(-ulook,dirToSurface);//reuse a
	    incidence_angles(i_range,j_dop)=  Uvar(acos(a),"rad");
	  }
	else
	  {
	    no_of_solution_grid(i_range,j_dop) = 0;//no surface intercept point
	  }
      }
    delete[] d_range;
    delete[] d_doppler;
  }


void compute_range_doppler(const Time& t,
			   const Frame& ftarget,
			   const string& target_name,
			   const Uvar& radius,
			   const Frame& fbeam,
			   const Uvar& lambda,
			   const Uvec& azim_grid,
			   const Uvec& elev_grid,
			   Array2D<PositionVector>& position_grid,
			   Umat& area_pixel,
			   Umat& incidence_angle_pixel,
			   Umat& range_pixel,
			   Umat& doppler_pixel,
			   Imat& valid_pixel)
{
 
  double** range_grid;
  double** doppler_grid;
  double** incidence_angle_grid;
  int** valid_grid;

  range_grid= (double**)make_array(sizeof(double),2,azim_grid.size(),elev_grid.size());
  doppler_grid=(double**)make_array(sizeof(double),2,azim_grid.size(),elev_grid.size());
  incidence_angle_grid=(double**)make_array(sizeof(double),2,azim_grid.size(),elev_grid.size());
  valid_grid=(int**)make_array(sizeof(int),2,azim_grid.size(),elev_grid.size());

  StateVector sc;
  ftarget.ephemeris(sc,"Cassini",t,"NONE");
  DirectionVector look("",fbeam,t,0,0,1);
  TargetGeom tg(t);
  for(unsigned int i=0;i<azim_grid.size();++i)
    for(unsigned int j=0;j<elev_grid.size();++j){
      look.setAzimuthElevation(azim_grid(i),elev_grid(j));
      tg.setState(sc);
      tg.setTarget(target_name, ftarget);
      tg.setLookDirection(look);
      if(tg.foundSurfaceIntercept()){
	valid_grid[i][j]=1;
	range_grid[i][j]=tg.range().getInUnits("km");
	doppler_grid[i][j]=tg.doppler(lambda).getInUnits("Hz");
	incidence_angle_grid[i][j]=tg.incidenceAngle().getInUnits("rad");
	position_grid(i,j)=tg.surfaceIntercept();
      }
      else{
	valid_grid[i][j]=0;
      }
      tg.reset(t);
    }

  for(unsigned int i=0;i<azim_grid.size()-1;++i)
    for(unsigned int j=0;j<elev_grid.size()-1;++j){
      if(valid_grid[i][j]*valid_grid[i+1][j]
	 *valid_grid[i][j+1]*valid_grid[i+1][j+1]==1){
	valid_pixel(i,j)=1;
	double dummy= range_grid[i][j]+range_grid[i][j+1];
	dummy += range_grid[i+1][j]+range_grid[i+1][j+1];
	range_pixel(i,j)=Uvar(dummy,"km")/4.0;

	dummy=doppler_grid[i][j]+doppler_grid[i][j+1];
	dummy+= doppler_grid[i+1][j] + doppler_grid[i+1][j+1];
	doppler_pixel(i,j)=Uvar(dummy,"Hz")/4.0;

	dummy=incidence_angle_grid[i][j]+incidence_angle_grid[i][j+1];
	dummy+=incidence_angle_grid[i+1][j]+incidence_angle_grid[i+1][j+1];
	incidence_angle_pixel(i,j)=Uvar(dummy,"rad")/4.0;

	area_pixel(i,j)=surfacearea(position_grid(i,j),
				     position_grid(i,j+1),
				     position_grid(i+1,j),
				     radius);
	area_pixel(i,j)+=surfacearea(position_grid(i+1,j+1),
				      position_grid(i,j+1),
				      position_grid(i+1,j),
				      radius);
	
      }
      else{
	valid_pixel(i,j)=0;
      }
    }

  free_array((void**)range_grid,2,azim_grid.size(),elev_grid.size());
  free_array((void**)incidence_angle_grid,2,azim_grid.size(),elev_grid.size());
  free_array((void**)doppler_grid,2,azim_grid.size(),elev_grid.size());
  free_array((void**)valid_grid,2,azim_grid.size(),elev_grid.size());
  /*
  double vx,vy,vz;
  vx=vel[FloatVector::X].getInUnits("m");
  vy=vel[FloatVector::Y].getInUnits("m");
  vz=vel[FloatVector::Z].getInUnits("m");
  double speed= vel.magnitude().getInUnits("m/s");
  double lambda_d=lambda.getInUnits("m");

  PositionVector target_center("",ftarget,t,0,0,0);
DirectionVector
  double dx,dy,dz;//direction x,y,z
  for(unsigned int i=0;i<azim_grid.size();++i)
    for(unsigned int j=0;j<elev_grid.size();++j){
      dx = sin(azim_grid(i))*cos(elev_grid(i));
      dy = sin(elev_grid(j));
      dz = cos(azim_grid(i))*cos(elev_grid(j));
      doppler_grid(i,j)=2.0*speed/lambda_d*(dx*vx + dy*vy + dz*vz);
    }
  */

  
}

void compute_min_max_range(const Time& t,
			   const StateVector& sc_state,
			   const string& target,
			   const Frame& fbeam,
			   const Uvec& azi_angles,
			   const Uvec& elev_angles,
			   Uvar& min_range,
			   Uvar& max_range)
  {
    if(azi_angles.size()!=elev_angles.size()) ErrorMessage("size mismatch");
    unsigned int N=azi_angles.size();
    Uvec ranges("",N);
    TargetGeom tg(t);
    for(unsigned int i=0;i<N;++i){
      DirectionVector look("",fbeam,t,0,0,1);
      look.setAzimuthElevation(azi_angles(i),elev_angles(i));
      tg.setState(sc_state);
      tg.setTarget(target);
      tg.setLookDirection(look);
      ranges(i)=tg.range();
      tg.reset(t);
    }
    min_max(min_range,max_range,ranges);
  }
void compute_bore_range_doppler(const Time& t,
				const StateVector& sc_state,
				const string& target,
				const Frame& fbeam,
				const Uvar& lambda,
				Uvar& bore_range,
				Uvar& bore_doppler)
  {
    DirectionVector bore("",fbeam,t,0,0,1);
    TargetGeom tg(t);
    tg.setState(sc_state);
    tg.setTarget(target);
    tg.setLookDirection(bore);
    bore_range=tg.range();
    bore_doppler = tg.doppler(lambda);
  }
//-----------------------------------------------
// min_max(min,max,data)
//
// return min and max values of Uvec data set
//-----------------------------------------------

void min_max(Uvar& min, Uvar& max, const Uvec& data )
  throw (Unit::UnitError,ErrorMessage)
  {
  unsigned int i_size = data.size();
  min = data(0);
  max = data(0);
  for (unsigned int i = 0; i < i_size; i++)
    {
    if (min >= data(i)) min = data(i);
    if (max <= data(i)) max = data(i);
    }
  
  }

void min_max(Uvar& min, Uvar& max, const Umat& data )
  throw (Unit::UnitError,ErrorMessage)
{
  unsigned int i_size = data.rows();
  unsigned int j_size = data.cols();
  min = data(0,0);
  max = data(0,0);
  for (unsigned int i = 0; i < i_size; i++)
  for (unsigned int j = 0; j < j_size; j++)
    {
      if (min >= data(i,j)) min = data(i,j);
      if (max <= data(i,j)) max = data(i,j);
    }
}

void min_max(Uvar& min, Uvar& max, const vector<Uvar>& data)
  
  {
    unsigned int N= data.size();
    if(N==0) 
      ErrorMessage("Utils.cpp:: min_max: container size is 0").throwMe();
    min = data[0];
    max = data[0];
    for(unsigned int i=0;i<N;++i){
      if(min>= data[i]) min=data[i];
      if(max<= data[i]) max=data[i];
    }
  }



void min_max( double& min, double& max, const Dvec& value)
{
  min= value(0);
  max= value(0);
  for(unsigned int i=0;i<value.size();++i){
    if(min>value(i)) min= value(i);
    if(max<value(i)) max=value(i);
  }
}

//----------------------------------------------------------------------------
// polyfit(poly_coeff, y, x)
//
// return coefficients of polynomial fit
// Solve the following matrix
// Y(data_number) = poly_coeff(1 x i_poly_order) 
//                        * X (i_poly_order x data_number)
//  Y * X_transpose = poly_coeff * X * X_transpose
//  poly_coeff = Y * X_transpose * inv(X * X_transpose)
//-----------------------------------------------------------------------------

void  polyfit(Dvec& poly_coeff, const Dvec&  y,  const Dvec&  x)
     throw (ErrorMessage)
  {
    unsigned int order = poly_coeff.size();
    unsigned int data_number = y.size();
    poly_coeff = 0.0;
    if (data_number != x.size())
      {
	throw ErrorMessage("where y = ax, number of y data is not same as number of x data");
      }
    if (order > data_number)
      {
	throw ErrorMessage("Not enough data to make a polynormial fit");
      }
    
  
 


  //-----------------------------------------
  //(1)struct x matrix
  // 1     1       1
  // x     x       x
  // x^2   x^2     x^2
  // ..
  // ..
  // x^n   x^n     x^n
  //
  //
  // (2) construct x_transposed 
  // x_trans(i,j) = x(j,i)
  // (3) construct square matrix (order x order)
  //  by multiplying x(order x data_number) * x_trans(data_numberx order)
  //----------------------------------------

  Dmat x_mat("x_matrix",order,data_number);
  Dmat x_trans("altitdue_matrix_transposed",data_number,order);
  for (unsigned int i = 0;i <order; i++)
  for (unsigned int j = 0; j < data_number; j++)
    {
      x_mat(i,j) = pow(x(j),i);
    }
   

  for (unsigned int i = 0;i <data_number; i++)
  for (unsigned int j = 0; j < order; j++)
    {
      x_trans(i,j) = x_mat(j,i);
    }
  

  Dmat  x_square("x_square_matrix",order,order);
  x_square = 0.0;

  
  for (unsigned int i = 0;i <order; i++)
  for (unsigned int j = 0; j < order; j++)
    {    
      for (unsigned int k = 0 ; k < data_number; k++)
	{
	  x_square(i,j)+= x_mat(i,k) *  x_trans(k,j);
	}
    }
 
  

  //-----------------------------------------------------
  //(4) Find an inverse matrix of x_square -> 
  //                         x_inverse_square
  //
  //-----------------------------------------------------
  

  Dmat x_inverse_square("x_inverse_square_matrix",order,order);
  x_inverse_square = x_square;
  x_inverse_square.inv();
  
  Dmat I_diag("1_matrix",order,order);
  I_diag = 0.0;
  for (unsigned int i = 0;i <order; i++)
  for (unsigned int j = 0; j < order; j++)
  {    
    for (unsigned int k = 0 ; k < order; k++)
      {
	I_diag(i,j) += x_square(i,k) * x_inverse_square(k,j);      
      }
    if (i == j)
      {
	if(fabs(I_diag(i,i)-1.0) > 1e-2)
	  { 
	    //cout<<"diagonal "<< I_diag(i,i)<<endl;
	    throw ErrorMessage("Matrix multiplication: A x A^-1 != 1");
	    
	  } 
      }
    else
      {
	if(fabs(I_diag(i,j)) > 1e-2)
	  { 
	    //cout<<"off diag "<< I_diag(i,j)<<endl;
	    throw ErrorMessage("Matrix multiplication: A x A^-1 != 1");
	   
	  } 
      }
  }
  


  //--------------------------------------------------------
  //(5) multiply alt_transpose to y: y(1xdata number)*
  // alt_trans(datanumber x order) = y_alt_trans (1 x order)
  //
  //-------------------------------------------------------
  Dvec y_x_trans("y_x_x_transposed",order);
  y_x_trans = 0.0;
  for (unsigned i = 0; i < order; i++)
  for (unsigned j = 0 ; j < data_number; j++)
    {
      y_x_trans(i) += y(j) * x_trans(j,i) ;
    }
  

  //----------------------------------------
  //(6) Finally, y_alt_trans = poly_coeff * alt_square
  // poly_coeff (1 x order) = y_alt_trans(1 x order) * alt_inverse_square(order x order)
  //---------------------------------------

  for (unsigned int i = 0; i < order; i++)
  for (unsigned int j = 0; j< order; j++) 
    {
      poly_coeff(i) += y_x_trans(j) * x_inverse_square(j,i);
    }
  

  }

//--------------------------- end of poly_fit-----------------


//-------------------------------------------------
//fitting data to a straight line: using chi2
// y = a + bx
//-----------------------------------------------------
void linearfit(const vector<Uvar>& x, const vector<Uvar>& y,
	       const vector<Uvar>& sig, const int& mwt,
	       Uvar& a, Uvar& b,
	       Uvar& sig_a, Uvar& sig_b,
	       double& chi2, double& q)
  {
    unsigned int N=x.size();
    Dvec x_data("",N), y_data("",N),sig_data("",N);
    if(x.size()!=y.size()
       || x.size()!= sig.size())
      ErrorMessage("utils.cpp::linearfit: input data size mismatch").throwMe();
  
    Uvar s_x,s_y;
    s_x = 0.0; s_y = 0.0;
    for(unsigned  int i=0;i<N;++i)
      {
	x_data(i) = x[i].getInUnits(x[0]);
	y_data(i)=y[i].getInUnits(y[0]);
	sig_data(i)=sig[i].getInUnits(sig[0]);
	s_x += fabs(x_data(i));
	s_y += fabs(y_data(i));
      }
    double tmp_a,tmp_b;//a+bx
    double siga,sigb;
    linearfit(x_data,y_data,sig_data,mwt,tmp_a,tmp_b,siga,sigb,chi2,q);
    //get unit correction
    s_x /= s_x;//
    s_y /= s_y;

    a = tmp_a* s_y;//add unit
    sig_a = siga * s_y;

    b =tmp_b* (s_y/s_x);//add unit
    sig_b = sigb * (s_y/s_x);
  
  }

 void linearfit(const Uvec& x, const Uvec& y,
		const Uvec& sig, const int& mwt,
		Uvar& a, Uvar& b,
		Uvar& sig_a, Uvar& sig_b,
		double& chi2, double& q)
  {
    unsigned int N = x.size();
    Dvec x_data("",N);
    Dvec y_data("",N);
    Dvec sig_data("",N);
    
    if(x.size()!=y.size()
       || x.size()!= sig.size())
      ErrorMessage("utils.cpp::linearfit: input data size mismatch").throwMe();
    
    Uvar s_x,s_y;
    s_x = 0.0;s_y=0.0;
    for(unsigned int i=0;i<N;++i)
      {
	x_data(i)=x(i).getInUnits(x(0));
	y_data(i)=y(i).getInUnits(y(0));
	sig_data(i)=sig(i).getInUnits(sig(0));
	s_x += fabs(x_data(i));
	s_y += fabs(y_data(i));
      }
    double tmp_a,tmp_b;//a+bx
    double siga,sigb;
    linearfit(x_data,y_data,sig_data,mwt,tmp_a,tmp_b,siga,sigb,chi2,q);
    //get unit correction
    s_x /= s_x;//
    s_y /= s_y;

    a = tmp_a* s_y;//add unit
    sig_a = siga * s_y;

    b =tmp_b *(s_y/s_x);//add unit
    sig_b = sigb * (s_y/s_x);
    
  }



void linearfit(const Dvec& x, const Dvec& y,
	       const Dvec& sig, const int& mwt,
	       double& a, double& b,
	       double& sig_a, double& sig_b,
	       double& chi2, double& q)
  {
    int N=x.size();
    if(x.size()!=y.size()
       || x.size()!= sig.size())
      ErrorMessage("utils.cpp::linearfit: input data size mismatch").throwMe();
    
    int i;
    double wt=0.0,t=0.0,sxoss=0.0;
    double sx=0.0, sy=0.0,st2=0.0,ss,sigdat=0.0;
    
    b=0.0;
    if(mwt){//weighted dat
      ss=0.0;
      for(i=0;i<N;++i){
	wt = 1.0/sig(i)/sig(i);
	ss += wt;
	sx += x(i)*wt;
	sy += y(i)*wt;}
    }
    else{//no weight
      for(i=0;i<N;++i){
	sx +=x(i);
	sy +=y(i);    }
      ss=N;
    }
    sxoss = sx/ss;
    
    if(mwt){
      for(i=0;i<N;++i)  {
	t=(x(i)-sxoss)/sig(i);
	st2 += t * t;
	b += t*y(i)/sig(i);}
    }
    else{
      for(i=0;i<N;++i){
	  t = x(i)-sxoss;
	  st2 +=t*t;
	  b +=t*y(i); }
    }

    b /= st2;
    a=(sy-sx*b)/ss;
    sig_a = sqrt((1.0+sx*sx/(ss*st2))/ss);
    sig_b=sqrt(1.0/st2);

    chi2=0.0;
    if(mwt==0){//no weight
      for(i=0;i<N;++i)
	chi2 += (y(i)-a - b*x(i)) * (y(i)-a-b*x(i));
      q=1.0;
      sigdat=sqrt(chi2/(N-2));
      sig_a *=sigdat;
      sig_b *=sigdat;
    }
    else{
      for(i=0;i<N;++i)
	chi2 +=(y(i)-a-b*x(i))*(y(i)-a-b*x(i))/sig(i)/sig(i);
      //  q=gammainc(0.5*(N-2),0.5*chi2);does not work
      q=1.0;
      cout<<"Incomplte Gamma function is not implemented yet"<<endl;
    }
    
  }

//-------------------------------------
//void bitset(unsigned int& target, 
//const unsigend int& b1, const unsigned int& b2,const unsigned int& setvalue)
// bitset utility for 16 bit unsigned int variable
//-------------------------------------
void bitset(unsigned short& target, const unsigned int& b1, 
	    const unsigned int& b2,const unsigned int& setvalue)
  {
  if (b2 < b1) 
    {
    ErrorMessage("bad position range").throwMe();
    }
  if (b1 > (8*sizeof(target)-1) || b2 > (8*sizeof(target)-1)) 
    {
    ErrorMessage("bad bit position").throwMe();
    }
  unsigned int max_setvalue = 1;//for given range bits
  for (unsigned short i = b1; i < b2+1;++i)
    {
      max_setvalue *=2;
      //for example,for 4 bits, max set values should be less than 2^4
    }  
  if (setvalue >= max_setvalue) 
    {
      ErrorMessage("Not enough ranges to set the value").throwMe();
    }
  unsigned short  pos_mask = 0;
  for (unsigned int i=b1; i <=b2;++i)
    {
    unsigned short submask = 1;
    if (i != 0)
      {
      for (unsigned int j =0 ; j <i;++j)
	{
	submask *= 2;
	}
      }
    pos_mask += submask;
    } 
  unsigned int neg_mask = ~pos_mask;
  unsigned short clear_bits_in_range= target & neg_mask;
  unsigned short bit_shifted_setvalue = setvalue <<b1;
  target= clear_bits_in_range + bit_shifted_setvalue;
  }

//------------------------------
//bitset of unsigned int
//-----------------------------
void bitset(unsigned int& target,
	    const unsigned int& bit,
	    const unsigned int& setvalue)
{
  if (bit > (8*sizeof(target)-1))
    {
      ErrorMessage("bad bit position").throwMe();
    }
  if(setvalue>1)
    {
      ErrorMessage("bitset should be either 0 or 1 ").throwMe();
    }
  

  unsigned int  pos_mask = 0;
  for (unsigned int i=0; i <=bit;++i)
    {
      unsigned int submask = 1;
      if(i!=0)
	{
	  for (unsigned int j =1 ; j <i;++j)
	    {
	      submask *= 2;
	    }
	}
      pos_mask += submask;
    }   
  unsigned int neg_mask = ~pos_mask;
  unsigned int clear_bits_in_range= target & neg_mask;
  unsigned int bit_shifted_setvalue = setvalue <<bit;
  target= clear_bits_in_range + bit_shifted_setvalue;
}

//----------------------------------------------
//void bitdisplay(const unsigned short& target)
//----------------------------------------------
void bitdisplay(const unsigned short& target)
  {
  unsigned int i_size = 8 *sizeof(target);
  for (unsigned int i = 0; i< i_size;++i)
    {
    cout<<bitget(target,i_size-i-1,i_size-i-1);
    }  
  cout<<endl;
  }

//-----------------------------------------
//Return Muhleman backscatter coefficient
//----------------------------------------
// Modified to use k1/sin(inc) if k2=0
double muhleman_backscatter(const double& k1, const double& k2, const Uvar& incidenceangle)
  {
    if (k1 <=0 || k2 < 0)
      {
      ErrorMessage e("Backscattering coefficient should be positive");
      e.throwMe();
      }

     if (incidenceangle < Uvar(0,"rad"))
      {
      ErrorMessage e("Incidence angle should be positive");
      e.throwMe();
      }

     if(k2==0){
       if(incidenceangle==Uvar(0,"rad")){
	 return 10.0;
       }
       else{
	 return MIN(10.0, k1/sin(incidenceangle));
       }
     }
     else{
       return (k1* cos(incidenceangle)
	    / pow(sin(incidenceangle)+k2*cos(incidenceangle),3.0));
     }
  }

//----------------------------------------
// Modified to use k1/sin(inc) if k2=0
double muhleman_backscatter(const double& k1, const double& k2, const double& incidenceangle_in_rad)
  {
    if (k1 <=0 || k2 < 0)
      {
      ErrorMessage e("Backscattering coefficient should be positive");
      e.throwMe();
      }
    if (incidenceangle_in_rad < 0.0)
      {
      ErrorMessage e("Incidence angle should be positive");
      e.throwMe();
      }

    if(k2==0){
      if(incidenceangle_in_rad==0){
	 return 10.0;
       }
       else{
	 return MIN(10.0, k1/sin(incidenceangle_in_rad));
       }
     }
    return (k1* cos(incidenceangle_in_rad)
	    / pow(sin(incidenceangle_in_rad)+k2*cos(incidenceangle_in_rad),3.0));
  }

double Lambertian(const double& x, const Uvar& angle){
  return( pow(cos(angle),x));
}

//----------------------------------------------
//calculate area bounded by three position vectors on the surface of a sphere
//---------------------------------------------
Uvar surfacearea(const PositionVector& p1, const PositionVector& p2, const PositionVector& p3, const Uvar& sphere_radius)
  {
    Uvar side_a = sphere_radius * fabs(p1.angle(p2).getInUnits("rad"));
    Uvar side_b = sphere_radius * fabs(p2.angle(p3).getInUnits("rad"));
    Uvar side_c = sphere_radius * fabs(p3.angle(p1).getInUnits("rad"));
    Uvar cos_c = (side_a*side_a + side_b*side_b - side_c*side_c)
      /(2.0 * side_a*side_b);
    Uvar angle_c = Uvar(acos(cos_c.getInUnits("")),"rad");
    return(0.5 * sin(angle_c) * side_a * side_b);
  }

//-----------------------
// Math functions
//-----------------------

//------------------------
//int round(x): round double
//---------------------------
int round_double(const double& x)
  {
    double y;
    if (x >= 0.0) y = x + 0.5;
    else y = x - 0.5;
    return(int(y));
  }


//---------------------------------
// double sinc(double d)
//  return sin(pi * d)/ pi * d
//---------------------------------

double  sinc(double d)
  {
  if (d == 0.0) 
    {
    return 1.0;
    }
  else
    { 
    return ( sin(pi * d)/(pi * d) );
    }
  }

//-----------------------------------------------------
// bilinear(a,b,d)
//
// bilinear: return value in 2-d double array d
// corresponding to doubel indices a and b
// This version of bilinear works in index space, so
// arguments a and b are floating point values in the
// range [0,N-1] where N is the number of rows,cols in d. 
// For a more general purpose version, see the routine following
// this one.
//-----------------------------------------------------

double bilinear(double a, double b, const Dmat& d) throw(ErrorMessage)
  {
  unsigned int Nrows, Ncols;
  d.size(Nrows,Ncols);
  if( a < 0 || b < 0 || a >= (double)Nrows || b >= (double)Ncols)
    {
    throw ErrorMessage("Bilinear interpolation out of range.");
    }

  unsigned int i1=int(a);
  unsigned int i2=i1+1;
  unsigned int j1=int(b);
  unsigned int j2=j1+1;

  // take in account edge values
  if(i1==Nrows-1) i2=i1;
  if(j1==Ncols-1) j2=j1;

  double cr2=a-i1;
  double cr1=1-cr2;
  double cc2=b-j1;
  double cc1=1-cc2;
    
  double retval=cr1*cc1*d(i1,j1)+cr1*cc2*d(i1,j2)+cr2*cc1*d(i2,j1)+
                cr2*cc2*d(i2,j2);
 
  return (retval);
  }


//-----------------------------------------------------
// Ubilinear(double a, double b,const Umat& d)
//: return value is Uvar
// bilinear: return value in 2-d double array d
// corresponding to doubel indices a and b
// This version of bilinear works in index space, so
// arguments a and b are floating point values in the
// range [0,N-1] where N is the number of rows,cols in d. 
// For a more general purpose version, see the routine following
// this one.
//-----------------------------------------------------

Uvar bilinear(double a, double b, const Umat& d) throw(ErrorMessage)
  {
  unsigned int Nrows, Ncols;
  d.size(Nrows,Ncols);
  if( a < 0 || b < 0 || a >= (double)Nrows || b >= (double)Ncols)
    {
    throw ErrorMessage("Bilinear interpolation out of range.");
    }

  unsigned int i1=int(a);
  unsigned int i2=i1+1;
  unsigned int j1=int(b);
  unsigned int j2=j1+1;

  // take in account edge values
  if(i1==Nrows-1) i2=i1;
  if(j1==Ncols-1) j2=j1;

  double cr2=a-i1;
  double cr1=1-cr2;
  double cc2=b-j1;
  double cc1=1-cc2;
    
  Uvar retval=cr1*cc1*d(i1,j1)+cr1*cc2*d(i1,j2)+cr2*cc1*d(i2,j1)+
                cr2*cc2*d(i2,j2);
 
  return (retval);
  }

//-----------------------------------------------------------------------------
// bilinear(x1,x2,y,out_of_bounds_value,newx1,newx2,newy)
//
// General purpose bilinear interpolation of a matrix of data.
// y is an array of data for the corresponding independent variables x1 and x2.
// Thus, y(i,j) = f(x1(i),x2(j)).
// A new matrix is constructed for the new vectors of independent variables.
// Thus, newy(i,j) = f(newx1(i),newx2(j)) using bilinear interpolation from
// x1,x2,y.  If newx1(i) or newx2(j) is out of bounds (not inside the range
// of values in x1,x2) then newy(i,j) will be set to the supplied
// out_of_bounds_value.
// The input newy will be resized to match newx1 and newx2.  Any original
// contents will be destroyed.
// Note, x1 and x2 are assumed to be sorted in ascending order.
//-----------------------------------------------------------------------------

void bilinear(const Uvec& x1, const Uvec& x2, const Umat& y,
  const Uvar& out_of_bounds_value,
  const Uvec& newx1, const Uvec& newx2, Umat& newy)
  {
  if (x1.size() != y.rows() || x2.size() != y.cols())
    {
    ErrorMessage e("bilinear: mismatched x1,x2,y");
    e.throwMe();
    }

  newy.resize(newx1.size(),newx2.size());

  int ix1 = 0;
  for (unsigned int i=0; i < newx1.size(); i++)
    {
    unsigned int k1;
    for (k1=ix1; k1 < x1.size(); k1++)
      {  // locate each newx1 within x1
      if (newx1(i) < x1(k1)) break;
      }
    ix1 = k1 - 1;
    int ix2 = 0;
    for (unsigned int j=1; j < newx2.size(); j++)
      {
      unsigned int k2;
      for (k2=ix2; k2 < x2.size(); k2++)
        {  // locate each newx2 within x2
        if (newx2(j) < x2(k2)) break;
        }
      ix2 = k2 - 1;
      if (ix1 == -1 || unsigned(ix1+1) == x1.size() ||
          ix2 == -1 || unsigned(ix2+1) == x2.size())
        {  // out of bounds
        ix1 = 0;  // reset to avoid violating array boundaries
        ix2 = 0;
        newy(i,j) = out_of_bounds_value;
        }
      else
        {  // bilinear interpolation
        Uvar ax1 = x1(ix1);
        Uvar ax2 = x2(ix2);
        Uvar ax1a = x1(ix1+1);
        Uvar ax2a = x2(ix2+1);
        Uvar p1 = y(ix1,ix2);
        Uvar p2 = y(ix1+1,ix2);
        Uvar p3 = y(ix1+1,ix2+1);
        Uvar p4 = y(ix1,ix2+1);
        double t = ((newx1(i) - ax1) / (ax1a - ax1)).getInUnits("");
        double u = ((newx2(j) - ax2) / (ax2a - ax2)).getInUnits("");
        newy(i,j) = (1-t)*(1-u)*p1 + t*(1-u)*p2 + t*u*p3 + (1-t)*u*p4;
        }
      }
    }
  }








//-------------------------------------------------------------
//decode Cassini data
// baq_decode
// baq mode : unsigned short
// baq threshold: unsigned short baq_threshold[0..23] for 24 BAQ blocks
//  BAQ look up table represents baq values from 0 to 127 to 0 to 254
//   to have a precision of 1/2
//  As resutl, real BAQ value from BAQ threshold array is
//     *************  float(BAQ)/2.0 *********************
// data: rdat_ (8-bit unsigned char vector : cdata) like footer_ and header_
// total number of data: Nrdat_ per burst 
// pul: number of pulses per burst
// number of data per pulse: Ns_ Nb = Nrad_/pulse
// number of data per each block: Ns= Ns_Nb/Nb_  (Nb_ = 24 for cassini radar system)
// decoding table
// 8_to_8 straight ==> 256 levels 127.5....-127.5
// 8_to_4 MSB      ==> 16 levels  120, 104, 88, 72, 56, 40,24,8,
//                                -8,-24,-40, -56, -72, -88, -104, -120
// 8_to_2 MSB      ==> 4 levels  96, 32, -32, -96
// 8_to_1 MSB      ==> 2 levels  64, -64
// 8_to_2 BAQ      ==> 4 levels  d[]* float(baq_threshold)/2.0
//                                 d[]={0.49,1.64,-0.49,-.164}
//                                  bit={ 00  01     10     11}
// 8_to_4 BAQ      ==> 16 levels d[]* float(baq_threshold)/2.0
//                     d[]={0.117, 0.355, 0.60, 0.861, 1.148, 1.479,1.891,2.498
//                        -0.117, -0.355, -0.60, -0.861, -1.148, -1.479, -1.891,-2.491}
//                   bit={0000,0001,0010,0011,0100,0101,0110,0111,
//                        1000,1001,1010,1011,1100,1101,1110,1111}
// things to be resolved: factor of 2 difference d[j]/2  
//Why 2 here in em_arc_ver2????????
// As of Nov 21, 2002, I understand why!  Because BAQ threshold value represent a value
// between 0 and 127 at an interval of 0.5
//----------------------------------------------------------------
void decodeCassiniData(const vector<unsigned char>& encoded_data,
		       const Array1D<unsigned int>& baq_threshold,
		       Ieb& ieb,
		       Array1D<float>& decoded_data,
		       unsigned int& Ndecoded,
		       float& offset)
  {
    Ndecoded = 0;//reset Ndecoded
    decoded_data = 0.0;
    offset = 0.0;
    //--------------------------------------------------------------
    //When r_mode is one of  4, 5,6,7 ( no BAQ values stored) or larger
    // than 12 ( r_mode >=12: spare),
    //we do not expect to have any data.
    //For other modes, number of pulses in the receiving window and
    // Nrdat_ must be non-zero
    // Return when rw == 0 or Nrdat_ == 0
    //-------------------------------------------------------------------
    unsigned int Nrdat = encoded_data.size(); 
    if( (ieb.getMod()>3  && ieb.getMod() <8) || ieb.getMod()>11) 
      {
	if (Nrdat != 0) 
	  {//passive data
	    throw ErrorMessage("Error: non-zero radar data in passive mode"); 
	  }
	return;
      }   
    else
      {//active mode
	if(Nrdat == 0)
	  {
	    cout<<"Warning::Nrdat is zero while radar is in active mode"<<endl;
	    return;
	  }
      }
      
    unsigned int baq_mode = ieb.getBaq();
    unsigned int Ndat_max = 32*1024;
    unsigned int Nb=24;
    unsigned int Nword = ieb.getNumberOfWords()-ieb.Nheader-ieb.Nfooter;
    unsigned int number_data_pri = ieb.getPointsInsidePRI();//points inside pri
    Ndecoded = ieb.getRadarDataPoints();//number of data to be decoded
    decoded_data.resize(Ndecoded);
    int rw_tmp = (int) (ieb.getPul());
    rw_tmp += ieb.getTroInPriUnits();  
    unsigned int rw = (unsigned int) rw_tmp;
    if(2*Nword != Nrdat)
      {
	cout<<"ieb mode "<<ieb.getMod()<<endl;
	cout<<"ieb baq "<<ieb.getBaq()<<endl;
	cout<<"adc  "<<ieb.getAdc()<<endl;
	cout<<"Number of words "<< Nword<<endl;
	cout<<"number data pri "<<number_data_pri<<endl;
	cout<<"number of pulse "<< ieb.getPul()<<endl;
	cout<<"number of bytes "<<Nrdat<<endl;
	//ErrorMessage("Utils.cpp: encode data size mismatch").throwMe();
      }
    //dummy variables
    int Ns_Nb, Ns, j;
    unsigned short iblock;
    double x;
    string input;
   
   
    //----------------------------------
    //how radar echo data are  divided?
    //assume 15 k word of data block
    // 2* 15 k = 30 k 8-bit data
    //for 8-2 bit BAQ
    //total number of data 4 *30 k = 120 k
    //number of rw : 10 (number of pulses in the receive window)
    //per each pulse: 120 k/ 10 = 12 k (= Ns_Nb)
    //for each block, 12 k / 24 (= Nb_) = 50 (=Ns)
    // for i_th data
    //  i % Ns_Nb : 0.. 12 k
    // (i % Ns_Nb) % Ns --> jth baq block (0..23)
    // if (i%Ns_Nb) % Ns > Nb_ --> insert into the last baq block(23)
    //-------------------------------------------------------
    if (Ndecoded > Ndat_max) 
      {//larger than 32K bytes
	cout<<"rw "<< rw<<endl;
	cout<<"number data pri "<< number_data_pri<<endl;
	throw ErrorMessage("Ndecoded > Ndat_max");
	return;
      }

    //----------------------------
    //Before decoding, take care of byte swapping
    //----------------------------   
    Array1D<unsigned char> radar_buffer("radar buffer",Nrdat);
    vector<unsigned char>::const_iterator  p= encoded_data.begin();
    unsigned int i_buffer_count =0;
    for (p = encoded_data.begin(); p != encoded_data.end();p++)
      {
	radar_buffer(i_buffer_count) = *p;
	i_buffer_count++;
      }
    for (unsigned int i = 0; i < Nrdat-1;i = i+2)
      {
	unsigned char a = radar_buffer(i);
	unsigned char b = radar_buffer(i+1);
	radar_buffer(i) = b;
	radar_buffer(i+1) = a;//byte swapping
      }
    
    
    //-------------------------------------------------------
    //cout <<"BAQ mode started"<<"baq_mode "<<baq_mode<<endl;
    //-------------------------------------------------------
    if(baq_mode == 0) //8-to- 2 bit BAQ (SAR mode)
      {      
	unsigned int i_data_count = 0;    
	Ns_Nb = number_data_pri;
	Ns = int(Ns_Nb / Nb);   
	double  d[4]={0.49,1.64,-0.49,-1.64};//00, 01, 10, 11 = 0.49 1.64 -0.49  -1.64    
	
	for (unsigned int i_buffer = 0; i_buffer < Nrdat;++i_buffer)
	for (int i = 0; i < 8; i = i + 2)
	  {
	    j = chrbitget(radar_buffer(i_buffer),i,i+1);
	    iblock = i_data_count % Ns_Nb;
	    iblock = int(iblock/Ns);
	    if (iblock >= Nb) {iblock = Nb -1 ;}
	    x =d[j]*(float)(baq_threshold(iblock))/2.0;	

	    if(i_data_count <Ndecoded)
	    decoded_data(i_data_count)= x;
	    else
	      {
		cout<<"ieb mode "<<ieb.getMod()<<endl;
		cout<<"ieb baq "<<ieb.getBaq()<<endl;
		cout<<"adc  "<<ieb.getAdc()<<endl;
		cout<<"Number of words "<< Nword<<endl;
		cout<<"number data pri "<<number_data_pri<<endl;
		cout<<"number of pulse "<< ieb.getPul()<<endl;
		cout<<"tro in pri unit "<<ieb.getTroInPriUnits()<<endl;
		cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
		cout<<endl;
	      }
	    i_data_count++;
	  }
	
	if (i_data_count < Ndecoded )
	  {
	    throw ErrorMessage(" number of decodede< collected data");
	  }        
      }     	
    
    else if (baq_mode == 1) //8-to-1 bit BAQ: sign bit and threshold only 
      {
	unsigned int i_data_count=0;
	Ns_Nb = number_data_pri;
	Ns = int(Ns_Nb/Nb);	  
	double  d[2]= {1.0, -1.0};//0 1 = 1.0 -1.0
	for (unsigned int i_buffer = 0; i_buffer < Nrdat;++i_buffer)
	for (int i = 0; i < 8; i++)
	{
	  j = chrbitget(radar_buffer(i_buffer),i,i);
	  iblock = i_data_count % Ns_Nb;
	  iblock = int(iblock/Ns);
	  if (iblock >= Nb) {iblock = Nb -1 ;}
	  x = d[j]*(float)(baq_threshold(iblock))/2.0;  

	  if(i_data_count<Ndecoded)
	  decoded_data(i_data_count) = x;
	  else
	    {
	      cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
	    }

	  i_data_count ++;
	}		
	    
	if (i_data_count < Ndecoded )
	  {
	    throw ErrorMessage(" number of decodede< collected data");
	  }  
      }  
    
    else if (baq_mode == 2)
      {
	cout<<"8-to 0 (standard sci data but no data)"<< baq_mode<<endl;
	Ndecoded = 0;
	return;
      }
    else if (baq_mode == 3) //ALTL DATA COMPRESSION
      {	   
	//ALTL compression method produces 16 bit data
	//number of data:number_data_pri
	//data type: sum of absolute value
	//However, encoded_data  is 8 bit.  So, need to add even and odd set to
	//construct new data
	
	//--------------------------
	//temporary holder
	//size = Nrdat :same as 8 bit size
	//real data size : pri_number (Nrdat /2)
	// Caution: first 32 bits corresponds to DC offset
	//First check data size
	//-------------------------------------------------
	
	if ( (Nrdat- 4) !=  2* number_data_pri)
	  ErrorMessage("Sab.cpp::ALTL compressed: number_data_pri does not match to 8-bit buffer size").throwMe();
	
	//---------------------------
	//Second, check if twos compliment was used
	//---------------------------
	offset = 0.0;
	unsigned int sign = chrbitget(radar_buffer(0),7,7);
	if (sign == 0)
	  {//normal representation
	    offset +=  chrbitget(radar_buffer(0),0,7)*pow(2,24);
	    offset +=  chrbitget(radar_buffer(1),0,7)*pow(2,16);
	    offset +=  chrbitget(radar_buffer(2),0,7)*pow(2,8);
	    offset +=  chrbitget(radar_buffer(3),0,7);    
	  }
	else 
	  {//twos compliment was used
	    // (2^B - N) : B = 32, N= normal way of reading
	    // change sign when the first bit is 1
	    offset = pow(2,32);
	    offset -=  chrbitget(radar_buffer(0),0,7)*pow(2,24);
	    offset -=  chrbitget(radar_buffer(1),0,7)*pow(2,16);
	    offset -=  chrbitget(radar_buffer(2),0,7)*pow(2,8);
	    offset -=  chrbitget(radar_buffer(3),0,7);  
	    offset *= -1.0;
	  }
	
	cout<<"total number of data"<<rw*number_data_pri<<endl;
	offset /=(double) (rw*number_data_pri);
	cout<<"altl dc offset"<<offset<<endl;
	
	for (unsigned int i_byte = 0; i_byte < 4;++i_byte)
	  {
	    cout<<i_byte<<"th byte "<<endl;
	    for (unsigned int i_bit =0; i_bit < 8;++i_bit)
	      {
		cout<<chrbitget(radar_buffer(i_byte),i_bit,i_bit);
	      }
	    cout<<endl;
	  }
	
	Ndecoded = number_data_pri;
	decoded_data.resize(Ndecoded);
	for (unsigned int i = 0; i < number_data_pri;++i)
	  {
	    unsigned int value1 = radar_buffer(2*i+4);
	    unsigned int value2 = radar_buffer(2*i+1+4);
	    decoded_data(i) = (float(value1) + float(value2)*256.0)/float(rw);
	  }
      }	
    else if (baq_mode == 4) // 8-4 MSB
      { 
	unsigned int i_data_count = 0;
	double d[8]={8.0,24.0,40.0,56.0,72.0,88.0,104.0,120.0};
	// 0000 0001 0001 ....1000 1001..
	for (unsigned int i_buffer = 0; i_buffer < Nrdat;++i_buffer)
	  {
	    j = chrbitget(radar_buffer(i_buffer),0,2);  
	    if(i_data_count<Ndecoded)
	      {
		decoded_data(i_data_count) = d[j]; 
		if (chrbitget(radar_buffer(i_buffer),3,3)==1) 
		  decoded_data(i_data_count)= -decoded_data(i_data_count);
	      }
	    else
	      {
		if(j!=0)
		  cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
	      }

	    i_data_count ++;  	    

	    j = chrbitget(radar_buffer(i_buffer),4,6);  
	    if(i_data_count<Ndecoded)
	      {
		decoded_data(i_data_count)=d[j] ;
		if (chrbitget(radar_buffer(i_buffer),7,7) == 1)  
		  decoded_data(i_data_count)= -decoded_data(i_data_count);
	      }
	    else
	      {
		if(j!=0)
		  cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
	      }
	   
	    i_data_count ++;
	  }	        

	if (i_data_count < Ndecoded )
	  {
	    throw ErrorMessage("baq_mode=4: decoded_data < actual data");
	  }

      }	  
    else if(baq_mode ==5) //8-to-8 straight
      {
	unsigned int i_data_count = 0;
	double  d[2]={1.0, -1.0};	   
	for (unsigned int i_buffer = 0; i_buffer < Nrdat;++i_buffer)
	  {	     
	    j = chrbitget(radar_buffer(i_buffer),7,7);
	    x = ( (float)(chrbitget(radar_buffer(i_buffer),0,6)) + 0.5 ) *d[j];  
	    if(i_data_count<Ndecoded)
	    decoded_data(i_data_count) = x;
	    else
	      {
		ErrorMessage("Utils.cpp: this should not happen").throwMe();
	      }
	    i_data_count ++;
	  }

	if (i_data_count < Ndecoded )
	  {	
	  throw ErrorMessage("baq_mode=5: decoded_data < actual data");
	  }     

      }
    else if (baq_mode ==6 || baq_mode ==7 )//low and high ALT mode 8-to-4 BAQ
      {
	unsigned int i_data_count = 0;
	Ns_Nb= number_data_pri;
	Ns = int(Ns_Nb/Nb);// number of data per each baq_threshold	     
	double d[8]={0.117,0.355,0.600,0.861,1.148,1.479,1.891,2.498};
	for (unsigned int i_buffer = 0; i_buffer < Nrdat;++i_buffer)
	  { 		
	    j = chrbitget(radar_buffer(i_buffer),0,2);
	    iblock = i_data_count % Ns_Nb;
	    iblock = int(iblock/Ns);
	    if (iblock >= Nb) {iblock = Nb -1 ;}
	    x =  d[j]*(float)(baq_threshold(iblock))/2.0; 
	    if(i_data_count<Ndecoded)
	      {
		decoded_data(i_data_count) = x;
		if(chrbitget(radar_buffer(i_buffer),3,3)==1) 
		  decoded_data(i_data_count) = -decoded_data(i_data_count);
	      }
	    else
	      {
		if(j!=0)
		  cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
	      }

	    i_data_count++;

	    j = chrbitget(radar_buffer(i_buffer),4,6);
	    iblock = i_data_count % Ns_Nb;
	    iblock = int(iblock/Ns);
	    if (iblock >= Nb) {iblock = Nb -1 ;}
	    x =  d[j]*(float)(baq_threshold(iblock))/2.0;  
	    if (i_data_count > Ndat_max) {throw ErrorMessage("i_data_count > Ndat_max");}		
	    
	    if(i_data_count<Ndecoded)
	      {
		decoded_data(i_data_count)=x;
		if(chrbitget(radar_buffer(i_buffer),7,7) == 1) 
		  decoded_data(i_data_count) = -decoded_data(i_data_count);
	      }
	    else
	      {
		if(j!=0)
		  cout<<"i_data_count "<<i_data_count<<" Ndecoded "<<Ndecoded<<" value "<<j<<endl;
	      }

	    i_data_count++;		
	  }
	if (i_data_count < Ndecoded )
	  {
	    throw ErrorMessage("baq_mode=6 or 7: decoded_data < actual data");
	  }
      }
    else 
      {
	throw ErrorMessage("no baq_mode found");
	Ndecoded = 0;
	return;
      }    
    if(ieb.getBaq()!=3) offset = decoded_data.mean();

  }


//**************************************************************************
// This java code is  an interactive demo of the first ellipse-specific
// direct fitting method presented in the papers:

//    M. Pilu, A. Fitzgibbon, R.Fisher ``Ellipse-specific Direct
//    least-square Fitting '' , IEEE International Conference on Image
//    Processing, Lausanne, September 1996. (poscript) (HTML) 

//    A. Fitzgibbon, M. Pilu , R.Fisher ``Direct least-square fitting of
//    Ellipses '' , International Conference on Pattern Recognition, Vienna,
//    August 1996. (poscript) - Extended version available as DAI Research
//    Paper #794

// The demo can be tried out at   
//   http://www.dai.ed.ac.uk/students/maurizp/ElliFitDemo/demo.html

// The code was written by Maurizio Pilu , University of Edinburgh.  The
// applet's graphic interface was much inspired by the Curve Applet
// written by Michael Heinrichs at SFU, Vancouver:
// http://fas.sfu.ca:80/1/cs/people/GradStudents/heinrica/personal/curve.html 
//**************************************************************************
//------------------
// solve: pvec[1]*x^2 + pvec[2]*y^2+pvec[3]*x*y+pvec[4]*x+pvec[5]*y+pvec[6]=0
//------------------

void EllipseFit(const Dvec& x_input, 
		const Dvec& y_input, 
		Dvec& ellipse_x_out,
		Dvec& ellipse_y_out)
  {
    if(x_input.size() != y_input.size()) ErrorMessage("input x,y data size mismatch").throwMe();
    if(x_input.size() < 6) ErrorMessage("Need more than 6 data points").throwMe();
    int np = x_input.size();
  
    // method doesn't work for fewer than 6 points
    if (np<6) return;
    
    double x[np+1],y[np+1];
    // Array1D<double>  x("x",np),y("y",np);//unit less
    for(int i =0;i<np;++i)
      {
	x[i+1]= x_input(i);
	y[i+1]= y_input(i);
      }

    //-------------------------------------------------------
    //solve for 
    // p[1]*x*x + p[2]*x*y + p[3]*y*y + p[4]*x +p[5]*y + p[6] =0
    //-----------------------------------------------------------

    double **D = new double*[np+1];
    double pvec[7];
    int c;
    for( c=0;c<np+1;c++) D[c]=new double[7];
    double *S[7];
    for( c=0;c<7;c++) S[c]=new double[7];
    double *Const[7];
    double *temp[7];
    double *L[7],*invL[7]; 
    double *C[7]; 
    double *V[7]; 
    double *sol[7];
    for( c=0;c<7;c++) 
      {
	L[c]=new double[7];
	invL[c] = new double[7];
	Const[c]=new double[7];
	temp[c]=new double[7];
	C[c]=new double[7];
	V[c]=new double[7];
	sol[c]=new double[7];
      }
    
    double d[7];
    double tx,ty;
    int nrot=0;    
    int i,j;
    
  
    // initialize the constraint matrix
    for(int i=0;i<7;i++)
    for(int j=0;j<7;j++)
      {
	Const[i][j]=0;
      }
    Const[1][3]=-2;
    Const[2][2]=1;
    Const[3][1]=-2;   
	
   
		
    // Now first fill design matrix
    for ( i=1; i <= np; i++)
      { 
	tx = x[i];
	ty = y[i];
	D[i][1] = tx*tx;
	D[i][2] = tx*ty;
	D[i][3] = ty*ty;
	D[i][4] = tx;
	D[i][5] = ty;
	D[i][6] = 1.0;
      }
		
    //pm(Const,"Constraint");
    // Now compute scatter matrix  S
    A_TperB(D,D,(double**)S,np,6,np,6);
    //pm(S,"Scatter");
		
    choldc((double**)S,6,(double**)L);    
    //pm(L,"Cholesky");
    
    inverse((double**)L,(double**)invL,6);
    //pm(invL,"inverse");
    
    AperB_T((double**)Const,(double**)invL,(double**)temp,6,6,6,6);
    AperB((double**)invL,(double**)temp,(double**)C,6,6,6,6);
    //pm(C,"The C matrix");
    
    jacobi((double**)C,6,d,(double**)V,nrot);
    //pm(V,"The Eigenvectors");  /* OK */
    //pv(d,"The eigevalues");
    
    A_TperB((double**)invL,(double**)V,(double**)sol,6,6,6,6);
    //pm(sol,"The GEV solution unnormalized");  /* SOl */
    
    // Now normalize them 
    for ( j=1;j<=6;j++)  /* Scan columns */
      {
	double mod = 0.0;
	for (i=1;i<=6;i++)
	  mod += sol[i][j]*sol[i][j];
	for (i=1;i<=6;i++)
	  sol[i][j] /=  sqrt(mod); 
      }
    
    //pm(sol,"The GEV solution");  /* SOl */
    
    double zero=10e-20;
    int  solind=0;
    for (i=1; i<=6; i++)
      {
	if (d[i]<0 && fabs(d[i])>zero)     
	  solind = i;
      }
    // Now fetch the right solution
    for (j=1;j<=6;j++)
      pvec[j] = sol[j][solind];
    //pv(pvec,"the solution");
    
    // DEBUG: display ellipse rotation
    float theta = 0.5*atan(pvec[2]/(pvec[3]-pvec[1]));

    //----------------------
    //new coordinate system: 
    //AA[1]*x*x + 2*AA[2]*x*y +AA[3]*y*y+ 2*AA[4]*x + 2*AA[5]*y+ AA[6] = 0
    // transformation
    // x_old = cos(theta)*x + sin(theta)*y
    // y_old = -sin(theta)*x + cos(theta)*y
    //-----------------------
    double AA[7];
    AA[1]=pvec[1]*cos(theta)*cos(theta)- pvec[2]*cos(theta)*sin(theta)+pvec[3]*sin(theta)*sin(theta);
    AA[2]=0.5*pvec[2]*(cos(theta)*cos(theta)-sin(theta)*sin(theta))+(pvec[1]-pvec[3])*sin(theta)*cos(theta);//=0.0
    AA[3]=pvec[1]*sin(theta)*sin(theta)+pvec[2]*sin(theta)*cos(theta)+pvec[3]*cos(theta)*cos(theta);
    AA[4]= 0.5*pvec[4]*cos(theta)- 0.5*pvec[5]*sin(theta);
    AA[5]= 0.5*pvec[4]*sin(theta) + 0.5*pvec[5]*cos(theta);
    AA[6] = pvec[6];
    if(fabs(AA[2]) > 0.001) ErrorMessage("transformation error").throwMe();
    
    //------------------------
    // finally equation looks like this
    // (x + AA[4]/AA[1]) ^2/(1.0/AA[1]) + (y + AA[5]/AA[3])^2/ (1.0/ccc) =
    //  -ggg + ddd^2/aaa + fff^2/ccc
    //----------------------------
    //cout<<aaa<<" "<<bbb<<" "<<ccc<<" "<<ddd<<" "<<fff<<" "<<ggg<<endl;
    //debug
    //for(c=1;c<=6;++c) cout<<AA[c]<<" ";
    //cout<<endl;
    
    double x_axis_pts[3],y_axis_pts[3];
    AA[6] = -AA[6] + AA[4]*AA[4]/AA[1] + AA[5]*AA[5]/AA[3];
    if(AA[6]/AA[1] < 0 || AA[6]/AA[3]<0 ) ErrorMessage("EllipseFit fails").throwMe();
    x_axis_pts[1] =-AA[4]/AA[1] + sqrt(AA[6]/AA[1]);
    x_axis_pts[2] = -AA[4]/AA[1] - sqrt(AA[6]/AA[1]);
    y_axis_pts[1] = -AA[5]/AA[3] + sqrt(AA[6]/AA[3]);
    y_axis_pts[2] = -AA[5]/AA[3] -sqrt(AA[6]/AA[3]);
 
    
    // x_old = cos(theta)*x + sin(theta)*y
    // y_old = -sin(theta)*x + cos(theta)*y
    ellipse_x_out.resize(4);
    ellipse_y_out.resize(4);
    
    if(AA[1] >= AA[3])
      {//x_axis becomes major
	// x_old = cos(theta)* x_major
	// y_old = -sin(theta)* x_major
	ellipse_x_out(0)= x_axis_pts[1] * cos(theta);
	ellipse_x_out(1) = x_axis_pts[2] * cos(theta);
	ellipse_y_out(0)=  x_axis_pts[1] * sin(theta);
	ellipse_y_out(1) = x_axis_pts[2]* sin(theta);
	// x_old = sin(theta)* y_minor
	// y_old = cos(theta)* y_minor
	ellipse_x_out(2)= -y_axis_pts[1] * sin(theta);
	ellipse_x_out(3) = -y_axis_pts[2] * sin(theta);
	ellipse_y_out(2)=  y_axis_pts[1] * cos(theta);
	ellipse_y_out(3) = y_axis_pts[2]* cos(theta);

      }
    else
      {//y becomes major
	// x_old = sin(theta)* y_major
	// y_old = cos(theta)* y_major
	ellipse_x_out(0)= -y_axis_pts[1] * sin(theta);
	ellipse_x_out(1) = -y_axis_pts[2] * sin(theta);
	ellipse_y_out(0)=  y_axis_pts[1] * cos(theta);
	ellipse_y_out(1) = y_axis_pts[2]* cos(theta);
	// x_old = cos(theta)* x_minor
	// y_old = -sin(theta)* x_minor
	ellipse_x_out(2)= x_axis_pts[1] * cos(theta);
	ellipse_x_out(3) = x_axis_pts[2] * cos(theta);
	ellipse_y_out(2)=  x_axis_pts[1] * sin(theta);
	ellipse_y_out(3) = x_axis_pts[2]* sin(theta);
      }
    

  }

void jacobi(double **a, int n, double *d , double **v, int nrot)      
  {
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c;
    
    double *b = new double[n+1];
    double *z = new double[n+1];
    
    for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
    }
    nrot=0;
    for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++)
	  sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
	/*    free_vector(z,1,n);
	      free_vector(b,1,n);  */
	return;
      }
      if (i < 4)
	tresh=0.2*sm/(n*n);
      else
	tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++) {
	  g=100.0*fabs(a[ip][iq]);
	  if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
	      && fabs(d[iq])+g == fabs(d[iq]))
	    a[ip][iq]=0.0;
	  else if (fabs(a[ip][iq]) > tresh) {
	    h=d[iq]-d[ip];
	    if (fabs(h)+g == fabs(h))
	      t=(a[ip][iq])/h;
	    else {
	      theta=0.5*h/(a[ip][iq]);
	      t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	      if (theta < 0.0) t = -t;
	    }
	    c=1.0/sqrt(1+t*t);
	    s=t*c;
	    tau=s/(1.0+c);
	    h=t*a[ip][iq];
	    z[ip] -= h;
	    z[iq] += h;
	    d[ip] -= h;
	    d[iq] += h;
	    a[ip][iq]=0.0;
	    for (j=1;j<=ip-1;j++) {
	      ROTATE(a,j,ip,j,iq,tau,s);
	    }
	    for (j=ip+1;j<=iq-1;j++) {
	      ROTATE(a,ip,j,j,iq,tau,s);
	    }
	    for (j=iq+1;j<=n;j++) {
	      ROTATE(a,ip,j,iq,j,tau,s);
	    }
	    for (j=1;j<=n;j++) {
	      ROTATE(v,j,ip,j,iq,tau,s);
	    }
	    ++nrot;
	  }
	}
      }
      for (ip=1;ip<=n;ip++) {
	b[ip] += z[ip];
	d[ip]=b[ip];
	z[ip]=0.0;
      }
    }
    //printf("Too many iterations in routine JACOBI");
  }
    
	
void ROTATE(double **a, int i, int j, int k, int l, double tau, double s) 
  {
    double g,h;
    g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);
    a[k][l]=h+s*(g-h*tau);
  }
	
//  Perform the Cholesky decomposition    
// Return the lower triangular L  such that L*L'=A  
void choldc(double **a, int n, double **l)
  {
    int i,j,k;
    double sum;
    double *p = new double[n+1];
    
    for (i=1; i<=n; i++)  {
      for (j=i; j<=n; j++)  {
	for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
	if (i == j) {
	  if (sum<=0.0)  
	    // printf("\nA is not poitive definite!");
	    {}
	  else 
	    p[i]=sqrt(sum); }
	else 
	  {
	    a[j][i]=sum/p[i];
	  }
      }
    }       
    for (i=1; i<=n; i++)  
      for (j=i; j<=n; j++)  
	if (i==j)
	  l[i][i] = p[i];
	else
	  {
	    l[j][i]=a[j][i];  
	    l[i][j]=0.0;
	  }    
  }


/********************************************************************/
/**    Calcola la inversa della matrice  B mettendo il risultato   **/
/**    in InvB . Il metodo usato per l'inversione e' quello di     **/
/**    Gauss-Jordan.   N e' l'ordine della matrice .               **/
/**    ritorna 0 se l'inversione  corretta altrimenti ritorna     **/
/**    SINGULAR .                                                  **/
/********************************************************************/
int inverse(double **TB, double **InvB, int N) 
  {  
    int k,i,j,p,q;
    double mult;
    double D,temp;
    double maxpivot;
    int npivot;
    
    double **B = new double*[N+1];
    double **A = new double*[N+1];
    double **C = new double*[N+1];
    for(i=0;i<N+1;i++)
      {
	B[i]=new double[N+2];
	A[i]=new double[2*N+2];
	C[i]=new double[N+1];
      }
    double eps = 10e-20;
    
    
    for(k=1;k<=N;k++)
      for(j=1;j<=N;j++)
	B[k][j]=TB[k][j];
    
    for (k=1;k<=N;k++)
      {
	for (j=1;j<=N+1;j++)
	  A[k][j]=B[k][j];
	for (j=N+2;j<=2*N+1;j++)
	  A[k][j]=(float)0;
	A[k][k-1+N+2]=(float)1;
      }
    for (k=1;k<=N;k++)
      {
	maxpivot=fabs((double)A[k][k]);
	npivot=k;
	for (i=k;i<=N;i++)
	  if (maxpivot<fabs((double)A[i][k]))
	    {
	      maxpivot=fabs((double)A[i][k]);
	      npivot=i;
	    }
	if (maxpivot>=eps)
	  {      if (npivot!=k)
	    for (j=k;j<=2*N+1;j++)
	      {
		temp=A[npivot][j];
		A[npivot][j]=A[k][j];
		A[k][j]=temp;
	      } ;
	  D=A[k][k];
	  for (j=2*N+1;j>=k;j--)
	    A[k][j]=A[k][j]/D;
	  for (i=1;i<=N;i++)
	    {
	      if (i!=k)
		{
		  mult=A[i][k];
		  for (j=2*N+1;j>=k;j--)
		    A[i][j]=A[i][j]-mult*A[k][j] ;
		}
	    }
	  }
	else
	  {  // printf("\n The matrix may be singular !!") ;
	    return(-1);
	  };
      }
    /**   Copia il risultato nella matrice InvB  ***/
    for (k=1,p=1;k<=N;k++,p++)
      for (j=N+2,q=1;j<=2*N+1;j++,q++)
	InvB[p][q]=A[k][j];
    return(0);
  }            /*  End of INVERSE   */



void AperB(double **A, double **B, double **res, 
	   int righA, int colA, int righB, int colB) 
  {
    int p,q,l;                                      
    for (p=1;p<=righA;p++)    
      {
	for (q=1;q<=colB;q++)                        
	  { 
	    res[p][q]=0.0;                            
	    for (l=1;l<=colA;l++)                     
	      {
		res[p][q]=res[p][q]+A[p][l]*B[l][q];  
	      }
	  }    
      }
  }                                                 

void A_TperB(double **A, double **B, double **res,
	     int righA, int colA, int righB, int colB) 
  {
    int p,q,l;                                      
    for (p=1;p<=colA;p++)                        
      {
	for (q=1;q<=colB;q++)                        
	  { 
	    res[p][q]=0.0;                            
	    for (l=1;l<=righA;l++)                    
	      {
		res[p][q]=res[p][q]+A[l][p]*B[l][q];  
	      }
	  }    
      }
  }

void AperB_T(double **A, double **B, double **res,
	     int righA, int colA, int righB, int colB) 
  {
    int p,q,l;                                      
    for (p=1;p<=colA;p++)                         
      {
	for (q=1;q<=colB;q++)                        
	  { 
	    res[p][q]=0.0;                            
	    for (l=1;l<=righA;l++)                    
	      {
		res[p][q]=res[p][q]+A[p][l]*B[q][l];  
	      }
	  }    
      }
  }




//---------------------------------------
//compute target frame's rotation angles
// and rotation rate in J2000
//---------------------------------------

void getTBC_Frame_RA_DEC_ROTATION_in_J2000(const Time& t,
						const string& target_name,
						Uvar& ra, 
						Uvar& dec, 
						Uvar& rotation_angle)
  {
    double et;
    t.getEt(et);
    if(strcasecmp(target_name.c_str(),"none")!=0){
      Frame ftarget("IAU_"+target_name,target_name);
      Frame fj2000("J2000",target_name);
      
      SpiceDouble FJ2000_TARGET[3][3];
      SpiceDouble* angles;
      angles = new SpiceDouble[3];
      Frame::rotationMatrix(fj2000,ftarget,et,FJ2000_TARGET);
      //frame rotation matrix from J2000 to target frame
      m2eul_c(FJ2000_TARGET,3,2,3,angles,angles+1,angles+2);
      ra =       Uvar(angles[2],"rad");
      dec=       Uvar(pi/2-angles[1],"rad");//90-pole_declination
      rotation_angle= Uvar(angles[0],"rad");
      delete [] angles;
      DirectionVector z_body(ftarget,et,0,0,1);
      z_body.representIn(fj2000);
      z_body.getRADEC(ra,dec);
    }
    else{//no target
      cout<<"---------------- Warning ---------------"<<endl;
      cout<<"No target is specified, zero_degree angle will be returned"<<endl;
      cout<<"-----------------------------------------"<<endl;
      ra =       Uvar(0,"rad");
      dec=       Uvar(0,"rad");
      rotation_angle= Uvar(0,"rad");
    }

  }


Uvar getTBC_Frame_Rotation_Rate_in_J2000(const Time& t,
					 const string& target_name,
					 const Uvar& interval)
  {
    if(interval <Uvar(0,"s")) ErrorMessage("Flyby.cpp:: getTarget_Rotation_Rate: time interval should be larger than 0 s").throwMe();
    Uvar rate=Uvar(0,"rad/s");
    double et;
    if(strcasecmp(target_name.c_str(),"none")!=0){
      Frame ftarget("IAU_"+target_name,target_name);
      Frame fj2000("J2000",target_name);
      SpiceDouble fj2000_to_ftarget2[3][3], fj2000_to_ftarget1[3][3];
      
      //at time t2    
      Time t2= t+interval/2.0;
      t2.getEt(et);
      Frame::rotationMatrix(fj2000,ftarget,et,fj2000_to_ftarget2);
      
      //at time t1
      Time t1 = t -interval/2.0;
      t1.getEt(et);
      Frame::rotationMatrix(fj2000,ftarget,et,fj2000_to_ftarget1);
    
      //rotation matrix
      SpiceDouble rot[3][3];
      mxmt_c(fj2000_to_ftarget2,fj2000_to_ftarget1,rot);
      
      //rotation angle and axis
      SpiceDouble rot_angle;
      SpiceDouble* rot_axis;
      rot_axis = new SpiceDouble[3];
      raxisa_c(rot, rot_axis, &rot_angle);
    
      rate= Uvar(rot_angle,"rad")/interval ;
      delete[] rot_axis;
    }
    return(rate);
  }




//--------------------------------
//compute spacecraft yaw, pitch roll
//-------------------------------
void computeYaw_Pitch_Roll(const StateVector& sc_state,
			   const Frame& ftarget,
			   const Frame& fsc,
			   Uvar& yaw,
			   Uvar& pitch, 
			   Uvar& roll)
  {
    StateVector sc= sc_state;
    Time t= sc.time();
    double et;
    t.getEt(et);
    sc.position().representIn(ftarget);
    sc.velocity().representIn(ftarget);
    //---------------
    //SCH(or TCN) cooridnate
    //---------------
    DirectionVector sc_z(" ",sc.position());
    DirectionVector sc_y(" ",cross(sc_z,sc.velocity()));
    DirectionVector sc_x(" ",cross(sc_y,sc_z));
    Frame ftcn("tnc",sc.position(),sc_x,sc_y,sc_z);//ftcn coordinate 

    //SPACE CRAFT FRAME
    DirectionVector x("",fsc,t,1,0,0);
    DirectionVector y("",fsc,t,0,1,0);
    DirectionVector z("",fsc,t,0,0,1);
    
    double z_rotation = dot(sc_x, x);
    //cout<<"z rotation "<<z_rotation<<endl;
    if(z_rotation < 0)
      {
	x = -x;
	y= -y;
      }
    
    Frame FSC("new sc frame",sc.position(),x,y,z);
    
    //transformation matrix
    SpiceDouble fsc_to_ftcn[3][3];
    Frame::rotationMatrix(FSC,ftcn,et,fsc_to_ftcn);//from FSC to ftcn
    
    SpiceDouble* Spice_angle;
    Spice_angle = new SpiceDouble[3];

    // from FTCN to FSC: 1-2-3 (roll pich yaw)
    //from FSC to FTCN: 3-2-1 (yaw pich roll)
    m2eul_c(fsc_to_ftcn,3,2,1,Spice_angle, Spice_angle+1,Spice_angle+2);
    yaw= Uvar(*Spice_angle,"rad");//roll
    pitch=Uvar(*(Spice_angle+1),"rad");//ptch
    roll=Uvar(*(Spice_angle+2),"rad");//yaw
    
    delete[] Spice_angle;

  }

void  fit_pitch_polynomial( const string& before_or_after_epoch,
			    const Uvar& start_time, 
			    const Uvar& end_time,
			    const Uvar& pitch_rate,
			    Dvec& poly_fit,
			    Uvar& mid_time){
  Uvar peak_pitch;
  Uvar t1,t2,t3;
  if(before_or_after_epoch=="inbound"){
    if(!(end_time<start_time)) ErrorMessage("Utils.cpp::fit_pitch_polynomial: end time should be earlier than start time").throwMe();
    peak_pitch= Uvar( (pitch_rate *(start_time - end_time)/2.0).getInUnits("rad"),"rad");
    t1 =  end_time;
    t3=  start_time;
    t2= (t1+t3)/2.0;
  }
  else if (before_or_after_epoch=="outbound"){
    if(!(end_time>start_time))ErrorMessage("Utils.cpp::fit_pitch_polynomial: end time should be later than start time").throwMe();
    peak_pitch= Uvar( (pitch_rate *(end_time - start_time)/2.0).getInUnits("rad"),"rad");
    t1= start_time;
    t3= end_time;
    t2= (t1+t3)/2.0;
  }
  else{
    ErrorMessage("Utils.cpp:: fit_pitch_polynomialstring option should be either inbound or outbound").throwMe();
  }
  mid_time= t2;


  cout<<"start, mid, and end time "<< t1<<" "<<t2<<" "<<t3<<endl;
  cout<<"peak pitch bias "<< peak_pitch.getInUnits("deg")<<" deg"<<endl;
  //condition
  // pitch at start time and end time=0
  // pitch rate at start and end time=0
  // pitch value at (start+middle)/2= peak value
  // fit a + b(t-t2)+c(t-t2)^2 + d(t-t2)^3+e(t-t2)^4
  Dvec B("",5);
  B(0)=0.0;
  B(1)=peak_pitch.getInUnits("deg");
  B(2)=0.0;
  B(3)=0.0;
  B(4)=0.0;


  Dmat A("",5,5);
  //at t=t1
  double x = (t1-t2).getInUnits("min");
  A(0,0)=1; A(0,1)= x ; A(0,2)= pow(x,2); A(0,3)= pow(x,3); A(0,4)=pow(x,4);
  //at t=t2
  A(1,0)=1; A(1,1)=0.0 ; A(1,2)=0.0; A(1,3)=0.0; A(1,4)=0.0;
 
  //at t=t3
  x=(t3-t2).getInUnits("min");
  A(2,0)=1; A(2,1)=x; A(2,2)= pow(x,2); A(2,3)=pow(x,3); A(2,4)=pow(x,4);

  //time derivative at t1
  x = (t1-t2).getInUnits("min");
  A(3,0)=0.0; A(3,1)= 1.0;A(3,2)= 2.0*x; A(3,3)=3.0*pow(x,2); A(3,4)=4.0*pow(x,3);

  //time derivative at t3
   x=(t3-t2).getInUnits("min");
   A(4,0)=0.0; A(4,1)=1.0; A(4,2)=2.0*x; A(4,3)=3.0*pow(x,2); A(4,4)=4.0*pow(x,3);

   A.inv();

   poly_fit=0.0;
   for (unsigned int i=0; i<5 ;++i)
     for(unsigned int j=0;j<5;++j)
       poly_fit(i)+= A(i,j)*B(j);
}



//---------------------------------
//Utility function for doppler centroid tracking
//--------------------------------

//compute doppler and doppler derivatives
void dop_derivatives(const Uvar& radius, 
		     const Uvar& altitude,
		     const FloatVector& v,
		     const Uvar& yaw,    
		     const Uvar& pitch, 
		     const Uvar& azimuth,
		     const Uvar& height_ref, 
		     const Uvar& range,
		     const Uvar& lambda, 
		     const string& look_direction,
		     const double& r_scale,
		     double&  r_f,
		     double& r_dfdy,
		     double& r_dfdp,
		     double& r_dfdvs,
		     double& r_dfdvc,
		     double& r_dfdhref,
		     double& r_dfdvh,
		     double& r_dfdwvl,
		     double& r_dfdasa)
  {
    //convert everything in double
    double r_radcur = radius.getInUnits("m");
    double  r_sch[3];
   
    r_sch[0]=0;
    r_sch[1]= 0;
    r_sch[2]= altitude.getInUnits("m");
       
    double r_schvel[3];
    r_schvel[0]=v[FloatVector::X].getInUnits("m/s");
    r_schvel[1]=v[FloatVector::Y].getInUnits("m/s");
    r_schvel[2]=v[FloatVector::Z].getInUnits("m/s");
    
    //vector rotation: opposite of frame rotation
    double r_yaw = yaw.getInUnits("rad");
    double r_pitch = pitch.getInUnits("rad");
    
    //azimuth bias
    double r_azesa = azimuth.getInUnits("rad");

    double r_href = height_ref.getInUnits("m");
    double r_wvl = lambda.getInUnits("m");
    double r_range = range.getInUnits("m");
    double i_lrl=-1.0;//default 
    if(look_direction=="left" || look_direction=="Left") i_lrl=1.0;
    else if (look_direction=="right" || look_direction=="Right") i_lrl = -1.0;
    else ErrorMessage("SchDopFunctions.cpp: invalid look direction").throwMe();
    //main work
    dop_derivatives( r_radcur,
		     r_sch,
		     r_schvel,
		     r_yaw,
		     r_pitch,
		     r_azesa,
		     r_href,
		     r_range,
		     r_wvl,
		     i_lrl,
		     r_scale,
		     r_f,
		     r_dfdy,
		     r_dfdp,
		    r_dfdvs,
		     r_dfdvc,
		     r_dfdvh,
		     r_dfdhref,
		     r_dfdwvl,
		     r_dfdasa);
  }

//call dop derivative using doubles: no use of uvar
void  dop_derivatives(const double& r_radcur,
		      double r_sch[3],
		      double r_schvel[3],
		      const double& r_yaw,
		      const double& r_pitch,
		      const double& r_azesa,
		      const double& r_href,
		      const double& r_range,
		      const double& r_wvl,
		      const double& i_lrl,
		      const double& r_scale,
		      double& r_f,
		      double& r_dfdy,
		      double& r_dfdp,
		      double&r_dfdvs,
		      double& r_dfdvc,
		      double& r_dfdvh,
		      double& r_dfdhref,
		      double& r_dfdwvl,
		      double& r_dfdasa)
{
    double  r_x1, r_x2, r_l3, r_look;
    int path=1;//if path 1, scott code
    //r_l3: angle between nadir and look direction: 
    // roll angle is dependent on r_l3, pitch and yaw angle
    r_x1 = (r_radcur + r_sch[2]);
    r_x2 = (r_radcur + r_href);
    r_l3 = (r_x1*r_x1 + r_range*r_range - r_x2*r_x2)/(2.0*r_x1*r_range);
    r_look = acos((r_l3 -sin(r_azesa)*sin(r_pitch))/(cos(r_pitch)*cos(r_azesa)));
    
    //debug
    //r_look *= -1.0;//vector rotation : opposite angle of frame rotation
    //cout<<"yaw pitch azi "<< r_yaw*radtodeg<<" "<<r_pit ch*radtodeg<<" "
    //<<r_azesa*radtodeg<<endl;
    //cout<<r_sch[2]<<endl;    
    //cout<<"r_x1 r_x2 r_l3 r_l3: " <<r_x1<<" "<<r_x2<<" "
    //<<acos(r_l3)*radtodeg<<endl;
    //cout<<"r_look "<<r_look*radtodeg<<endl;
    //cout<<"velocity :"<< r_schvel[0]<<"  "<<r_schvel[1]<<" "
    //<< r_schvel[2]<<endl;
    
       	
    //  some useful derivatives
    double  r_dlookdhref, r_dcoslookdp, r_dsinlookdp;
    double  r_dcoslookdasa,r_dsinlookdasa;
    r_dlookdhref = (1.0/(sin(r_look)*cos(r_pitch)*cos(r_azesa)))
      *(r_x2/(r_x1*r_range));


    double r_lvsch[3], r_dlvdy[3], r_dlvdp[3], r_dlvdasa[3];
    r_lvsch[0] = (cos(r_look)*sin(r_pitch)*cos(r_yaw) 
		  + sin(r_look)*sin(r_yaw)*i_lrl)* cos(r_azesa) 
      - sin(r_azesa)*cos(r_pitch)*cos(r_yaw);//OK: confirmed
   
    r_lvsch[1] = (-cos(r_look)*sin(r_pitch)*sin(r_yaw) 
		  + sin(r_look)*cos(r_yaw)*i_lrl)* cos(r_azesa) 
      + sin(r_azesa)*cos(r_pitch)*sin(r_yaw);//OK:confirmed
    r_lvsch[2] = -cos(r_look)*cos(r_pitch)*cos(r_azesa) 
      - sin(r_azesa)*sin(r_pitch);//OK:confirmed

    //cout<<"look vector in sch frame "<< r_lvsch[0]<<" "<<r_lvsch[1]<<" "<<r_lvsch[2]<<endl;
    if(path)
      {
	//cout<<"SH code "<<endl;
	r_dcoslookdp = (1.0/(cos(r_pitch)*cos(r_azesa)))
	  *((sin(r_look)*sin(r_pitch) +sin(r_azesa))/cos(r_pitch));
	r_dcoslookdasa = (1.0/cos(r_azesa))
	  *((cos(r_look)*tan(r_azesa))
	    /cos(r_pitch) +   tan(r_pitch)/cos(r_azesa));
	
	r_dlvdy[0] = (-cos(r_look)*sin(r_pitch)*sin(r_yaw) 
		      + sin(r_look)*cos(r_yaw)*i_lrl)*cos(r_azesa) 
	  +  sin(r_azesa)*cos(r_pitch)*sin(r_yaw);
	r_dlvdy[1] = (-cos(r_look)*sin(r_pitch)*cos(r_yaw) 
		      - sin(r_look)*sin(r_yaw)*i_lrl)*cos(r_azesa) 
	  +  sin(r_azesa)*cos(r_pitch)*cos(r_yaw);
	r_dlvdy[2] = 0.      ;  
	                        
    
	r_dlvdp[0] = (cos(r_look)*cos(r_pitch)*cos(r_yaw) 
		      + sin(r_pitch)*cos(r_yaw)*r_dcoslookdp 
		      + r_dsinlookdp*sin(r_yaw)*i_lrl)*cos(r_azesa) 
	  - sin(r_azesa)*cos(r_yaw)*r_dcoslookdp;
	r_dlvdp[1] = (-cos(r_look)*cos(r_pitch)*sin(r_yaw)
		      -sin(r_pitch)*sin(r_yaw)*r_dcoslookdp+
		      r_dsinlookdp*cos(r_yaw)*i_lrl)*cos(r_azesa) 
	  + sin(r_azesa)*sin(r_yaw)*r_dcoslookdp ;
	r_dlvdp[2] = (cos(r_look)*sin(r_pitch) 
		      - cos(r_pitch)*r_dcoslookdp)*cos(r_azesa) 
	  -  sin(r_azesa)*r_dsinlookdp ;
	r_dlvdasa[0] = -cos(r_azesa)*cos(r_pitch)*cos(r_yaw) 
	  - sin(r_azesa)*(cos(r_look)*sin(r_pitch)*cos(r_yaw) 
			  + sin(r_look)*sin(r_yaw)*i_lrl) 
	  + cos(r_azesa)*(sin(r_pitch)*cos(r_yaw)*r_dcoslookdasa
			  + i_lrl*sin(r_yaw)*r_dsinlookdasa);
	r_dlvdasa[1] = cos(r_azesa)*sin(r_yaw)*cos(r_pitch) 
	  -sin(r_azesa)*(-cos(r_look)*sin(r_pitch)*sin(r_yaw) 
			 + sin(r_look)*cos(r_yaw)*i_lrl) 
	  + cos(r_azesa)*(-sin(r_pitch)*sin(r_yaw)*r_dcoslookdasa 
			  + i_lrl*cos(r_yaw)*r_dsinlookdasa);
	r_dlvdasa[2] = -cos(r_azesa)*sin(r_pitch) 
	  + sin(r_azesa)*cos(r_pitch)*cos(r_look) 
	  - cos(r_azesa)*cos(r_pitch)*r_dcoslookdasa;
      }
    else
      {
	//cout<<"my code "<<endl;
	r_dcoslookdp = (1.0/(cos(r_pitch)*cos(r_pitch)*cos(r_azesa)))
	  *(r_l3*sin(r_pitch) +sin(r_azesa));
	r_dcoslookdasa = (1.0/(cos(r_azesa)*cos(r_azesa)*cos(r_pitch)))
	  *(r_l3*sin(r_azesa)+sin(r_pitch));
  	r_dlvdy[0] = (-cos(r_look)*sin(r_pitch)*sin(r_yaw) 
		      + sin(r_look)*cos(r_yaw)*i_lrl)*cos(r_azesa) 
	  +  sin(r_azesa)*cos(r_pitch)*sin(r_yaw);
	r_dlvdy[1] = (-cos(r_look)*sin(r_pitch)*cos(r_yaw) 
		      - sin(r_look)*sin(r_yaw)*i_lrl)*cos(r_azesa) 
	  +  sin(r_azesa)*cos(r_pitch)*cos(r_yaw);
	r_dlvdy[2] = 0.0      ;  
	                        
    
	r_dlvdp[0] = (cos(r_look)*cos(r_pitch)*cos(r_yaw) 
		      + sin(r_pitch)*cos(r_yaw)*r_dcoslookdp 
		      + r_dsinlookdp*sin(r_yaw)*i_lrl)*cos(r_azesa) 
	  +sin(r_azesa)*sin(r_pitch)*cos(r_yaw);
	  //- sin(r_azesa)*cos(r_yaw)*r_dcoslookdp;
	r_dlvdp[1] = (-cos(r_look)*cos(r_pitch)*sin(r_yaw)
		      -sin(r_pitch)*sin(r_yaw)*r_dcoslookdp+
		      r_dsinlookdp*cos(r_yaw)*i_lrl)*cos(r_azesa) 
	  -sin(r_azesa)*sin(r_pitch)*sin(r_yaw);
	  //+ sin(r_azesa)*sin(r_yaw)*r_dcoslookdp ;
	r_dlvdp[2] = (cos(r_look)*sin(r_pitch) 
		      - cos(r_pitch)*r_dcoslookdp)*cos(r_azesa) 
	  -sin(r_azesa)*cos(r_pitch);
	  // -  sin(r_azesa)*r_dsinlookdp ;

	r_dlvdasa[0] = -cos(r_azesa)*cos(r_pitch)*cos(r_yaw) 
	  - sin(r_azesa)*(cos(r_look)*sin(r_pitch)*cos(r_yaw) 
			  + sin(r_look)*sin(r_yaw)*i_lrl) 
	  + cos(r_azesa)*(sin(r_pitch)*cos(r_yaw)*r_dcoslookdasa
			  + i_lrl*sin(r_yaw)*r_dsinlookdasa);
	r_dlvdasa[1] = cos(r_azesa)*sin(r_yaw)*cos(r_pitch) 
	  -sin(r_azesa)*(-cos(r_look)*sin(r_pitch)*sin(r_yaw) 
			 + sin(r_look)*cos(r_yaw)*i_lrl) 
	  + cos(r_azesa)*(-sin(r_pitch)*sin(r_yaw)*r_dcoslookdasa 
			  + i_lrl*cos(r_yaw)*r_dsinlookdasa);
	r_dlvdasa[2] = -cos(r_azesa)*sin(r_pitch) 
	  + sin(r_azesa)*cos(r_pitch)*cos(r_look) 
	  - cos(r_azesa)*cos(r_pitch)*r_dcoslookdasa;
      }
    r_dsinlookdp = -(1.0/tan(r_look))*r_dcoslookdp;
    r_dsinlookdasa = -(1.0/tan(r_look))*r_dcoslookdasa;
    
  
    
    //double r_lvsch[3];
    
    //double M11, M13, M21, M23, M31, M33;//matrix element
    //double M12, M22, M32;
    //M11 = cos(r_yaw)*cos(r_pitch);
    //M12 = cos(r_yaw)*sin(r_pitch)*sin(r_look)*i_lrl + sin(r_yaw)*cos(r_look);
    //M13 = -cos(r_yaw)*sin(r_pitch)*cos(r_look) + sin(r_yaw)*sin(r_look);
    //M21 = -sin(r_yaw)*cos(r_pitch);
    //M22 = -sin(r_yaw)*sin(r_pitch)*sin(r_look)*i_lrl 
    //+ cos(r_yaw)*cos(r_look);
    //M23 = sin(r_yaw)*sin(r_pitch)*cos(r_look) + cos(r_yaw)*sin(r_look);
    //M31 = sin(r_pitch);
    //M32 = -cos(r_pitch)*cos(r_look);
    //M33 = cos(r_pitch)*cos(r_look);
    //r_lvsch[0]= -M11*sin(r_azesa) - M13*cos(r_azesa); 
    //r_lvsch[1]= -M21*sin(r_azesa) - M23*cos(r_azesa);
    //r_lvsch[2]= -M31*sin(r_azesa) - M33*cos(r_azesa);
    //cout<<"M11 and M13 "<<M11<<" "<<M12<<" "<<M13<<endl;
    //cout<<"M21 and M23 "<<M21<<" "<<M22<<" "<<M23<<endl;
    //cout<<"M31 and M33 "<<M31<<" "<<M32<<" "<<M33<<endl;
    //cout<<"look vector "<< r_lvsch[0]<<" "<<r_lvsch[1]<<" "
    //<< r_lvsch[2]<<endl;
    //cout<<"norm "<< r_lvsch[0]* r_lvsch[0] + r_lvsch[1]* r_lvsch[1] 
    //+ r_lvsch[2]* r_lvsch[2]<<endl;
    
    //---------------------------------------------------------
    //    compute the derivatives of the look vector wrt href
    //---------------------------------------------------------
    double r_dlvdhref[3];
    r_dlvdhref[0] = (-sin(r_look)*sin(r_pitch)*cos(r_yaw) 
		     + cos(r_look)*sin(r_yaw)*i_lrl)*cos(r_azesa);
    r_dlvdhref[1] = (sin(r_look)*sin(r_pitch)*sin(r_yaw) 
		     + cos(r_look)*cos(r_yaw)*i_lrl)*cos(r_azesa);
    r_dlvdhref[2] = sin(r_look)*cos(r_pitch)*cos(r_azesa) ;                                     
    r_dlvdhref[0] = r_dlvdhref[0]*r_dlookdhref ;
    r_dlvdhref[1] = r_dlvdhref[1]*r_dlookdhref ;
    r_dlvdhref[2] = r_dlvdhref[2]*r_dlookdhref ;
    
    //---------------------------------------------------------
    //    compute the velocity vector derivatives
    //---------------------------------------------------------
    double r_dvdvs[3], r_dvdvc[3], r_dvdvh[3];
    r_dvdvs[0] = 1.0    ;
    r_dvdvs[1] = 0.;
    r_dvdvs[2] = 0.;
    
    r_dvdvc[0] = 0.;
    r_dvdvc[1] = 1.0  ;
    r_dvdvc[2] = 0.;
    
    r_dvdvh[0] = 0.;
    r_dvdvh[1] = 0.;
    r_dvdvh[2] = 1.0;
    
       
    r_f = (2.0/r_wvl)*dot(r_schvel,r_lvsch);
    r_dfdwvl = -r_f/(r_wvl*r_scale);
    r_dfdy = (2.0/r_wvl)*dot(r_schvel,r_dlvdy) /r_scale;
    r_dfdp = (2.0/r_wvl)*dot(r_schvel,r_dlvdp)/r_scale;
    r_dfdhref = (2.0/r_wvl)*dot(r_schvel,r_dlvdhref);
    r_dfdasa = (2.0/r_wvl)*dot(r_schvel,r_dlvdasa)/r_scale;
    r_dfdvs = (2.0/r_wvl)* dot(r_dvdvs,r_lvsch);
    r_dfdvc = (2.0/r_wvl)* dot(r_dvdvc,r_lvsch);
    r_dfdvh = (2.0/r_wvl)* dot(r_dvdvh,r_lvsch);
    
}


//-------------------------------------------------------
//   This program is copied version of S.H.'s svdvecfit.f
//   FILE NAME: svdvecfit.f
//   
//   DATE WRITTEN: 01/02/95 
//   
//   PROGRAMMER: Scott Hensley
//   
//   FUNCTIONAL DESCRIPTION: This routine does a least squares fit 
//   to a vector valued observation least squares problem.
//   
//   ROUTINES CALLED: gaussj,svbksb,svdcmp,funcs
//   
//   NOTES: funcs is a user supplied function giving the jacobian
//   of the observation parameters wrt to fit parameters. This routine
//   is a generalization of Numerical Recipes svdfit. Note that this
//   routine can also be used in a nonlinear least squares procedure
//   by iterating properly.
//
//   Solves the least problem 
//
//             T   -1     -1     T   -1 
//    A = (AMAT COV  AMAT)  (AMAT COV  )VOBS 
//
//    where AMAT is the jacobain of the observations vs parameters,
//    COV is the covriance matrix of observations
//    and VOBS is the vector of observations. 
//
//    r_a should be passed in with current best estimate of values
//   
//   UPDATE LOG: 
//         
//  4/17/95 - Reversed order of r_vecin, r_vobs, and r_cov    SJS
//            revmoved r_vt, cleaned up parameter list
//   
//----------------------------------------------------------

//--------------------------
// i_mp: number of input points (range, doppler)
// i_rd: number of measurements for each point : 1
// i_fp: number of fit parameters: 18
// r_vecin(i_fp, i_mp): values of 18 parameters
// r_vobs(i_rd,i_mp): difference between
//                                       doppler data and estimation
// r_covarr: (i_rd,i_rd,i_mp): all set to 1.0 (some kind of significance)
// i_np: 2 + number of beams (means yaw, pictch, and each beam's azi bais)
// r_a(i_np): inial estimate  least squares for each point
// r_at2(i_np): delta to add to previous solution
// r_u(i_np,i_np) svd matrix, orthogonal
// r_v(i_np,i_np) svd matrix, orthogonal
// r_w(i_np): svd matrix , orthognal
// l_chiq: boolian
//-----------------------------
void svdvecfit(const unsigned int& i_mp,
	       const unsigned int& i_rd,
	       const unsigned int& i_fp,
	       const Dmat&  r_vecin,
	       const Dmat&  r_vobs,
	       const D3D& r_cov,
	       const unsigned int& i_np,
	       Dvec& r_a,
	       Dvec& r_at2,
	       Dmat& r_u,
	       Dmat& r_v,
	       Dvec& r_w,
	       Dmat& r_chisq,
	       const bool& l_chisq,
	       const Array1D<int>& i_paramest,
	       const Array1D<int>& i_usedata)
{
  //check size
  if(r_vecin.rows()!= i_fp ) 
    ErrorMessage("SchDopFunctions.cpp: r_vecin.rows()!=i_fp ").throwMe();
  if(r_vobs.rows()!=i_rd )
    ErrorMessage("SchDopFunctions.cpp: r_vobs.rows()!=i_rd ").throwMe();
  if(r_cov.sizeDim1()!=i_rd || r_cov.sizeDim2()!=i_rd )
    ErrorMessage("SchDopFunctions.cpp: r_cov is not a 3D matrix of i_rd x i_rd").throwMe();
  if(r_a.size()!=i_np|| r_at2.size()!=i_np) ErrorMessage("SchDopFunctions.cpp: r_a size is not i_np").throwMe();
  if(r_u.rows()!= i_np || r_u.cols()!=i_np||
     r_v.rows()!=i_np || r_v.cols()!=i_np ||
     r_w.size()!=i_np )
    ErrorMessage("SchDopFunctions.cpp: r_u, r_v, and r_w  size is not i_npxi_np").throwMe();
  if(r_chisq.rows()!=i_rd )
    ErrorMessage("SchDopFunctions.cpp: r_u, r_v, and r_w  size is not i_rdxi_mp").throwMe();

  //--------------------------- some initial settings------------
  unsigned int I_RDE=1;
  unsigned int I_NPE=7;// 2(yaw and pitch) + beam numbers
  if(I_RDE != i_rd) 
    ErrorMessage("SchDopFunctions.cpp: I_RDE != i_rd").throwMe();
  if(I_NPE != i_np) 
    ErrorMessage("SchDopFunctions.cpp: I_NPE != i_np").throwMe();
  

  unsigned int i,j,k,i_pts;
  double R_TOL, R_LAMBDA;
  Dmat r_covtemp("",I_RDE,I_RDE);
  Dmat r_am("",I_NPE,I_RDE);
  Dmat r_amat("",I_NPE,I_NPE);
  Dvec r_ptot("",I_NPE);
  double r_wmax,r_thres;
  Dmat r_b("",I_RDE,1);
  Dvec r_chird("",I_RDE);
  R_TOL = 1.0e-20;
  R_LAMBDA=1.0;
  
  //reset r_u =0
  //reset r_ptot=0;
  for(i=0;i<i_np;++i)
    {
      for(j=0;j<i_np;++j)
	{
	  r_u(i,j)=0.0;
	}
      r_ptot(i)=0.0;
    }


  //debugging parameters
  //vector<Uvar> dfdy,dfdp;
  //vector<Uvar> index;
  //loop over input points
  
  for (i_pts=0;i_pts<i_mp;++i_pts)
    {
      //invert covariant matrix of the observation
      for(i=0;i<i_rd;++i)
	for(j=0;j<i_rd;++j)
	  r_covtemp(i,j) = r_cov(i,j,i_pts);//was set to be 1
      
      if(i_rd==1) 
	r_covtemp(i_rd-1,i_rd-1) = 1.0/r_covtemp(i_rd-1,i_rd-1);//our case
      else gaussj(r_covtemp,r_b);//general case -need to implement later

      //--------------------------------------------
      //for(i=0;i<i_fp;++i)
      //cout<<" i r_vec "<< i<<" "<<r_vecin(i,i_pts)<<endl;
      //get required Jacobian matrix: use only differentials
      // such as df/dy, df/dp,df/da
      // single data points
      //r_amat(0,0) = df/dy 
      //r_amat(0,1) = df/dp
      //r_amat(0,2) =df/dbeam_azi
      //----------------------------------------------
    
      funcs(i_pts, i_rd,i_fp, r_vecin.getCol(i_pts),i_np,r_a,r_amat,
	    i_paramest,i_usedata);

      //------------------------------------
      //dfdy.push_back(r_amat(0,0));
      //dfdp.push_back(r_amat(0,1));
      //Uvar  range= r_vecin(10,i_pts) + r_vecin(17,i_pts)*r_vecin(11,i_pts);//
      //index.push_back(range/1000.0);
      //if(i_pts==0){
      //cout<<"first element in schdop "<< r_vecin.getCol(0)<<endl;
      //cout<<"return "<<r_amat(0,0)<<" "<<r_amat(0,1)<<endl;
      //}
      //debug
      //for(i=0;i<i_rd;++i)
      //for(j=0;j<i_np;++j)
      //cout<<" svdfecfit i j r_amat "<<i<< " "<<j<< " "<< r_amat(i,j)<<endl;
      //   multiply amat transpose by the inverse cov matrix
      // r_am = r_amat * r_covtemp:: basically keep df/dy,df/dp
      //--------------------------------------------------------

      for( i=0;i<i_np;++i){
	for(  j=0;j<i_rd;++j){
	  r_am(i,j) = 0.0;
	  for( k=0;k<i_rd;++k)
	    r_am(i,j) = r_am(i,j) + r_amat(k,i)*r_covtemp(k,j);//differential/variance
	}
      }

      //--------------------------------------      
      //debug
      //for( i=0; i< i_np;++i)
      //for (j=0;j<i_rd;++j)
      //cout<<"i,j,r_am = "<< i<<" "<<j<<" "<<r_am(i,j)<<endl;    
      //for( k=0;k<i_rd;++k)
      //for( j=1;j<i_np;++j)
      //  cout<<"k j r_amat "<< k<<" "<<j<<" "<< r_amat(k,j)<<endl;
      //    multiply am by amat
      //------------------------
      //Sum up all r_u = [r_amat *r_covtemp]  x [r_amt ] 
      // form a matrix [ (df/dy)**2, (df/dy)(df/dp) ....]
      //------------------------
      
      for( i=0;i<i_np;++i)
      for( j=0;j<i_np;++j)
	{
	  for( k=0;k<i_rd;++k)
	    r_u(i,j) = r_u(i,j) + r_am(i,k)*r_amat(k,j);
	}
      //---------------------------------------------------------
      //   multilpy am by vobs
      //cout<<"i_pts r_vobs=  "<<i_pts<< " "<<r_vobs(0,i_pts)<<endl;
      // form a matrix
      //  delta_f * df/dy
      //-------------------------------------------------
      for(i=0; i<i_np;++i)
	for(k=0;k<i_rd;++k)
	  r_ptot(i) = r_ptot(i) + r_am(i,k)*r_vobs(k,i_pts);
    }//for each points  
    
  //-----------------------------------------------------
  //Plot a;
  //  a.addXY(index,"",dfdy," ",line("solid","red",1),sym("none"));
  //  a.addXY(index,"",dfdp," ",line("solid","black",1),sym("none"));
  //  a.show("x");  
  //   find the SVD of the r_u matrix
  //debug      
  //cout<<"r_u matrix before sdvcmp"<<endl;
  //for( i=0;i<i_np;++i){
  //cout<<endl;
  //for(j=0;j<i_np;++j)
  //  cout<<r_u(i,j)<<" ";
  //}
  //cout<<endl;
  //
  //cout<<"r_ptot "<<r_ptot<<endl;
  //debug  
  //Dmat temp_u("",r_u.rows(),r_u.cols());
  //temp_u=r_u;
  //Dvec temp_w("",r_w.size());
  //temp_w=r_w;
  //Dmat temp_v("",r_v.rows(),r_v.cols());
  //temp_v =r_v;
  //svdcmp_cpp(temp_u, temp_w, temp_v);
  //cout<<"After avdcmp tmp "<<endl;
  //for( i=0;i<i_np;++i){
  //cout<<endl;
  //for(j=0;j<i_np;++j)
  // cout<<temp_u(i,j)<<" ";
  //}
  //-------------------------------------------

  svdcmp_cpp(r_u,r_w,r_v);

  //------------------------------------------
  //cout<<endl;
  //cout<<"After svdcmp r_u "<<endl;
  //for( i=0;i<i_np;++i){
  //cout<<endl;
  //for(j=0;j<i_np;++j)
  // cout<<r_u(i,j)<<" ";
  //}
  //cout<<endl; 
  //cout<<"r_w "<<r_w<<endl;
  //
  //cout<<endl;
  //cout<<"r_v "<<endl;
  //for( i=0;i<i_np;++i){
  //cout<<endl;
  //for(j=0;j<i_np;++j)
  //  cout<<r_v(i,j)<<" ";
  //}
  //---------------------------------------------
   
  r_wmax = 0.0;
  for( i=0;i<i_np;++i){
    if(r_w(i) > r_wmax)
      r_wmax = r_w(i);}
  r_thres = r_wmax*R_TOL;
  for( i=0;i<i_np;++i){
    if(r_w(i) < r_thres)
      r_w(i) = 0.0;}

  //--------------------------------
  //cout<<"r_thres = "<<r_thres<<endl;
  //for( i=0;i<i_np;++i)
  //cout<<"w = "<<i<<" "<<r_w(i)<<endl;
  // use the svbksb routine to solve for the desired parameters
  //  cout<<"Entering svbksb..."<<endl;
  //---------------------------------
  //SVBKSB
  //--------------------------------
  
  svbksb(r_u,r_w,r_v,r_ptot,r_at2);
 
      
  //    update the r_a vector
  //cout<<"r_at2 "<<r_at2<<endl;
  //cout<<"r_at2 /number "<< r_at2/double(i_mp)<<endl;
  
  for( i=0;i<i_np;++i){
    r_at2(i) = r_at2(i)*i_paramest(i);
    r_a(i) = r_at2(i)/R_LAMBDA + r_a(i);
  }
  
  //    evaluate the chisq array (linearized version)
  if(l_chisq)
    {
      //    loop over data points
      for(i=0;i<i_rd;++i) 
	r_chird(i) = 0.0;
      double chisq=0.0;
      for( i=0;i<i_mp;++i)
	{
	  funcs(i,i_rd,i_fp,r_vecin.getCol(i),i_np,r_a,r_amat,
		i_paramest,i_usedata);
	  
	  for( j=0;j<i_rd;++j){
	    r_chisq(j,i) = 0.0;
	    for( k=0; k<i_np;++k)
	      r_chisq(j,i) = r_chisq(j,i) + r_amat(j,k)*r_at2(k);
	    //cout<<"r_chisq = "<<i<<" "<<j<<" "<<r_chisq(j,i)<<" "<<r_vobs(j,i)<<endl;
	    r_chisq(j,i) = r_covtemp(j,j)*(r_chisq(j,i)-r_vobs(j,i))
	      *(r_chisq(j,i) -  r_vobs(j,i));
	    chisq += r_chisq(j,i);
	    r_chird(j) = r_chird(j) + r_chisq(j,i);
	  }
	}
      
      chisq = sqrt(chisq/(2.*double(i_mp)));
      //  cout<<"r_chisq =  "<<chisq<<endl;
    }   
}

//---------------------
//svdvar
//----------------------




void  funcs(const unsigned int& i_q,
	    const unsigned int& i_rd,
	    const unsigned int& i_fp,
	    const Dvec& r_vecin,
	    const unsigned int& i_ma,
	    const Dvec& r_a,
	    Dmat& r_amat,
	    const Array1D<int>& i_paramest,
	    const Array1D<int>& i_usedata)
{
  double r_radcur                     ;//radius of curvature
  double r_sch[3]                     ;//SCH vector                 
  double r_schvel[3]                  ;//SCH velocity vector
  double r_yaw                        ;//yaw angle
  double r_pitch                      ;//pitch angle
  double r_href                       ;//reference height
  double r_range                      ;//range 
  double r_wvl                        ;//wavelength
  double r_scale                      ;//scale for some derivatives
  double i_lrl                       ;//left/right looking indicator
  double r_azesa                      ;//azimuth steering angle
  
  double r_f                          ;//Doppler
  double r_dfdy                       ;//derivative wrt to yaw
  double r_dfdp                       ;//derivative wrt to pitch
  double r_dfdvs                      ;//derivative wrt to s velocity
  double r_dfdvc                      ;//derivative wrt to c velocity
  double r_dfdvh                      ;//derivative wrt to h velocity
  double r_dfdhref                    ;//derivative wrt to height reference
  double r_dfdwvl                     ;//derivative wrt to wavelength
  double r_dfdasa                     ;//derivative wrt to azimuth steering angle
  
  unsigned int i,j,i_beam;
  unsigned int MAX_BEAMS=5;
  //  unsigned int  i_numset = 2 + MAX_BEAMS;
  unsigned int i_yaw=0;//first parameter
  unsigned int i_pitch=1;//second parameter
  // unsigned int i_rdd=1;
  Array1D<unsigned int> i_asaa("",MAX_BEAMS);

  


  r_radcur = r_vecin(0);
  r_sch[0] = r_vecin(1) ;
  r_sch[1] = r_vecin(2);
  r_sch[2] = r_vecin(3);
  r_schvel[0] = r_vecin(4);
  r_schvel[1] = r_vecin(5);
  r_schvel[2] = r_vecin(6);
  r_yaw = r_vecin(7);
  r_pitch = r_vecin(8);
  r_href = r_vecin(9);
  r_range = r_vecin(10) + r_vecin(17)*r_vecin(11);//bin number * delta_r

     
  r_wvl = r_vecin(12);
  r_scale = r_vecin(13);
  i_beam = (unsigned int) round_double(r_vecin(14));
  i_lrl = r_vecin(15);
  r_azesa = r_vecin(16);

  for(i=0;i<MAX_BEAMS;++i)
    i_asaa(i) = 2 + i;//help index beam1's parameter is 2
  

  //  compute derivatives

  dop_derivatives(r_radcur,r_sch,r_schvel,r_yaw,r_pitch,r_azesa,
		  r_href,r_range,r_wvl,i_lrl,r_scale,r_f,
		  r_dfdy,r_dfdp,r_dfdvs,r_dfdvc,
		  r_dfdvh,r_dfdhref,r_dfdwvl,r_dfdasa);

  //debug
  //if(fabs(r_dfdy)<0.001)
  //{
  //  cout<<"beam number "<<i_beam<<endl;
  //  cout<<"beam azi "<<r_azesa*radtodeg<<endl;
  //  cout<<"df/dy r_amat(0,0) "<<r_dfdy<<" "<<r_amat(0,0)<<endl;
  //  cout<<"range "<<r_range/1000.0<<endl;      
  //  cout<<"r_radcur "<< r_radcur/1000.0<<endl;
  //  cout<<"r_sch "<<r_sch[0]/1000.0<<" "<<r_sch[1]/1000.0<<" "<<r_sch[2]/1000.0<<endl;
  //  cout<<"r_scv "<<r_schvel[0]/1000.0<<" "<<r_schvel[1]/1000.0<<" "<<r_schvel[2]/1000.0<<endl;
  //  cout<<"r_yaw "<< r_yaw*radtodeg<<endl;
  //  cout<<"r_pitch "<<r_pitch*radtodeg<<endl;
  //  cout<<"r_azesa "<<r_azesa*radtodeg<<endl;
  //  cout<<"r_dfdy "<<r_dfdy<<endl;
  //  cout<<"r_dfdp "<<r_dfdy<<endl;
  //  cout<<"r_dfda "<<r_dfdasa<<endl;
  //}
  //     populate matrix
  r_amat(0,i_yaw) = r_dfdy;//first parameter: i_yaw=0
  r_amat(0,i_pitch) = r_dfdp;//second parameter: i_pitch=1
  for( i=0; i<MAX_BEAMS;++i)
    {
      if(i ==  i_beam)
	r_amat(0,i_asaa(i)) = r_dfdasa;
      else
	r_amat(0,i_asaa(i)) = 0.0;
    }
  //     set to zero the derivatives all parameters not being estimated

  //cout<<"fitting parameters "<<i_paramest<<endl;
  //cout<<"used data "<<i_usedata<<endl;

  for( i=0;i<i_ma;++i)
    for(j=0;j<i_rd;++j)
      r_amat(j,i) = r_amat(j,i)*i_paramest(i)*i_usedata(j);
 
  //debug
  //for( i=0;i<i_ma;++i)
  //for(j=0;j<i_rd;++j)
  //  cout<<"functs: j, i, r_mat "<<j<<" "<<i<<" "<< r_amat(j,i) <<endl;   
}






void svdvar(Dmat& r_v, Dvec& r_w, Dmat& r_cvm)
{
  unsigned int m;
  m = r_v.rows();
  if(r_v.rows() != r_v.cols()) 
    ErrorMessage("r_v is not a square matrix ").throwMe();
  if(r_w.size()!=m) ErrorMessage("r_w has a wrong size").throwMe();
  if(r_cvm.rows()!= m || r_cvm.cols()!=m) r_cvm.resize(m,m);
  unsigned int i,j,k;
  double sum;

  //variance of fitting coefficient
  Dvec wti("",m);
  for(i=0;i<m;++i){
    wti(i)=0.0;
    if(r_w(i)!=0.0) wti(i)=1.0/(r_w(i)*r_w(i));
  }
   
  for(i=0;i<m;++i){
    for(j=0; j<=i;j++){
      sum=0.0;
      for(k=0;k<m;++k)
	sum=sum+r_v(i,k)*r_v(j,k)*wti(k);
      r_cvm(i,j)=sum;
      r_cvm(j,i)=sum;
    }
  }
}


//----------------------------------
//Start of gaussj
//----------------------------------
void gaussj(Dmat& a,   Dmat& b)   
{
  unsigned int n,m;
  n=a.rows();// n x n square matrix
  if(a.cols()!=m) ErrorMessage("gaussj: a is not a square matrix").throwMe();
  
  m = b.cols();//n x m matrix
  if(b.rows()!=n) ErrorMessage("gaussj: b is not a  matrix of n x m").throwMe();
  
  unsigned int  i,icol,irow,j,k,l,ll;
  Array1D<unsigned int> indxc("",n),indxr("",n),ipiv("",n);
  double big,dum,pivinv;
  
  for(j=0;j<n;++j) ipiv(j)=0;
  for( i=0;i<n;++i){
    big=0.0;
    for(j=0;j<n;++j){
      if(ipiv(j)!=1){
	for( k=0;k<n;++k){
	  if (ipiv(k)==0){
	    if (fabs(a(j,k))>=big){
	      big=fabs(a(j,k));
	      irow=j;
	      icol=k;
	    }
	  }
	  else if (ipiv(k)>1) 
	    ErrorMessage("singular matrix in gaussj").throwMe();
	}
      }
    }      
    ipiv(icol)=ipiv(icol)+1;
    if (irow!=icol){
      for( l=0; l<n;++l){
	dum=a(irow,l);
	a(irow,l)=a(icol,l);
	a(icol,l)=dum;}
      for( l=0;l<m;++l){
	dum=b(irow,l);
	b(irow,l)=b(icol,l);
	b(icol,l)=dum;}
    }
    indxr(i)=irow;
    indxc(i)=icol;
    if (a(icol,icol)==0) ErrorMessage("singular matrix in gaussj").throwMe();
    pivinv=1.0/a(icol,icol);
    a(icol,icol)=1.0;
    for(l=0;l<n;++l)
      a(icol,l)=a(icol,l)*pivinv;
    for(l=0;l<m;++l)
      b(icol,l)=b(icol,l)*pivinv;

    for( ll=0; ll<n;++ll){
      if(ll!=icol){
	dum=a(ll,icol);
	a(ll,icol)=0.0;
	for(l=0;l<n;++l)
	  a(ll,l)=a(ll,l)-a(icol,l)*dum;
	for( l=0;l<m;++l)
	  b(ll,l)=b(ll,l)-b(icol,l)*dum;
      }
    }
  }

  for(  l=n-1; l >=0;--l){
    if(indxr(l)!=indxc(l)){
      for( k=0;k<n;++k){
	dum=a(k,indxr(l));
	a(k,indxr(l))=a(k,indxc(l));
	a(k,indxc(l))=dum;}
      if(l==0) break;//out of loop , do not subtract 1
    }
  }
}

//----------------------------------------------------------
// End of Gaussj
//--------------------------------------------------------
//wrapper for the use of c-style singular value decomposition
void  svdcmp(Dmat& a,	     Dvec& w,	     Dmat& v)
  {
    unsigned int m = a.rows();
    unsigned int n = a.cols();
    if(v.rows()!=n && v.cols()!=n) v.resize(n,n);
    if(w.size()!=n)  w.resize(n);
    
    int M, N;
    M = int(m)+1;
    N = int(n)+1;
    double* a_matrix[M];
    for(int i=0;i<M;++i)
      a_matrix[i]=new double[N];//MxN

    for(unsigned int i=0;i<m;++i)
      for(unsigned int j=0;j<n;++j)
	a_matrix[i+1][j+1]= a(i,j);

    //debug 
    //cout<<"Into svdcmp "<<endl;
    //for(unsigned int i=1;i<=m;++i){
    //cout<<endl;
    //for(unsigned int j=1;j<=n;++j)
    //cout<<a_matrix[i][j]<<"";}
    //cout<<endl;

    double* v_matrix[N];
    for(int i=0;i<N;++i)
      v_matrix[i]= new double[N];//NxN
    double w_array[N];
    svdcmp_c( (double**)a_matrix,m, n, w_array, (double**) v_matrix); 

    for(unsigned int i=0;i<m;++i)
      for(unsigned int j=0;j<n;++j)
	a(i,j)= a_matrix[i+1][j+1];

    for(unsigned int i=0;i<n;++i)  w(i)=w_array[i+1];

    for(unsigned int i=0;i<n;++i)
      for(unsigned int j=0;j<n;++j)
	v(i,j)=v_matrix[i+1][j+1];
  }

// C stype singular value decomposition
void svdcmp_c(double **a, int m, int n, double w[], double **v)
  {
    
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    vector<double> rv1;
    rv1.resize(n+1);//from 1 to n 
    g=scale=anorm=0.0;

    //cout<<"In svdcmp_nrc "<<endl;
    //for( i=1;i<=m;++i){
    //cout<<endl;
    //for(j=1;j<=n;++j)
    //cout<<a[i][j]<<"";}
    //cout<<endl;


    for (i=1;i<=n;i++) {
      l=i+1;
      rv1[i]=scale*g;
      g=s=scale=0.0;
      if (i <= m) {
	for (k=i;k<=m;k++) scale += fabs(a[k][i]);
	//for (k=i;k<=m;k++) 
	//{scale += fabs(a[k][i]);}
	if (scale) {
	  for (k=i;k<=m;k++) {
	    a[k][i] /= scale;
	    s += a[k][i]*a[k][i];
	  }
	  f=a[i][i];
	  g = -SIGN_C(sqrt(s),f);
	  h=f*g-s;
	  a[i][i]=f-g;
	  for (j=l;j<=n;j++) {
	    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	    f=s/h;
	    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	  }
	  for (k=i;k<=m;k++) a[k][i] *= scale;
	}
      }
      w[i]=scale *g;//cout<<"w[i] "<<w[i]<<endl;
      g=s=scale=0.0;
      if (i <= m && i != n) {
	for (k=l;k<=n;k++) scale += fabs(a[i][k]);
	if (scale) {
	  for (k=l;k<=n;k++) {
	    a[i][k] /= scale;
	    s += a[i][k]*a[i][k];
	  }
	  f=a[i][l];
	  g = -SIGN_C(sqrt(s),f);
	  h=f*g-s;
	  a[i][l]=f-g;
	  for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	  for (j=l;j<=m;j++) {
	    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	  }
	  for (k=l;k<=n;k++) a[i][k] *= scale;
	}
      }
      //debug
      //cout<<"anorm and w+rv1 "<< anorm<<" "<< fabs(w[i])+fabs(rv1[i])<<endl;
      anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
     
    }
    for (i=n;i>=1;i--) {
      if (i < n) {
	if (g) {
	  for (j=l;j<=n;j++)
	    v[j][i]=(a[i][j]/a[i][l])/g;
	  for (j=l;j<=n;j++) {
	    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	  }
	}
	for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
      }
      v[i][i]=1.0;
      g=rv1[i];
      l=i;
    }
    for (i=MIN(m,n);i>=1;i--) {
      l=i+1;
      g=w[i];
      for (j=l;j<=n;j++) a[i][j]=0.0;
      if (g) {
	g=1.0/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	  f=(s/a[i][i])*g;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (j=i;j<=m;j++) a[j][i] *= g;
      } else for (j=i;j<=m;j++) a[j][i]=0.0;
      ++a[i][i];
    }
    for (k=n;k>=1;k--) {
      for (its=1;its<=30;its++) {
	flag=1;

	
	for (l=k;l>=1;l--) {
	  nm=l-1;
	  //debug
	  //cout<<"fabs(rv1) and anorm"<<fabs(rv1[l])+(float)anorm
	  //<<" "<<(float) anorm<<endl;
	  //cout<<"fabs(w) + anorm "<< fabs(w[nm])+(float)anorm<<" "
	  //<<(float)anorm<<endl;
	  //cout<<"k l "<<k<<" "<<l<<endl;
	  if ((float)(fabs(rv1[l])+anorm) == (float)anorm) {
	    flag=0;
	    break;
	  }
	  if ((float)(fabs(w[nm])+anorm) == (float)anorm) break;
	}
	if (flag) {
	  c=0.0;
	  s=1.0;
	  for (i=l;i<=k;i++) {
	    f=s*rv1[i];
	    rv1[i]=c*rv1[i];
	    if ((float)(fabs(f)+anorm) == (float)anorm) break;
	    g=w[i];
	    h=pythag(f,g);
	    w[i]=h;
	    h=1.0/h;
	    c=g*h;
	    s = -f*h;
	    for (j=1;j<=m;j++) {
	      y=a[j][nm];
	      z=a[j][i];
	      a[j][nm]=y*c+z*s;
	      a[j][i]=z*c-y*s;
	    }
	  }
	}
	z=w[k];
	if (l == k) {
	  if (z < 0.0) {
	    w[k] = -z;
	    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	  }
	  break;
	}
	if (its == 30) ErrorMessage("no convergenece").throwMe();
	 
	
	x=w[l];
	nm=k-1;
	y=w[nm];
	g=rv1[nm];
	h=rv1[k];
	f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	g=pythag(f,1.0);
	f=((x-z)*(x+z)+h*((y/(f+SIGN_C(g,f)))-h))/x;
	c=s=1.0;
	for (j=l;j<=nm;j++) {
	  i=j+1;
	  g=rv1[i];
	  y=w[i];
	  h=s*g;
	  g=c*g;
	  z=pythag(f,h);
	  rv1[j]=z;
	  c=f/z;
	  s=h/z;
	  f=x*c+g*s;
	  g = g*c-x*s;
	  h=y*s;
	  y *= c;
	  for (jj=1;jj<=n;jj++) {
	    x=v[jj][j];
	    z=v[jj][i];
	    v[jj][j]=x*c+z*s;
	    v[jj][i]=z*c-x*s;
	  }
	  z=pythag(f,h);
	  w[j]=z;//cout<<"w[j] "<<w[j]<<endl;
	  if (z) {
	    z=1.0/z;
	    c=f*z;
	    s=h*z;
	  }
	  f=c*g+s*y;
	  x=c*y-s*g;
	  for (jj=1;jj<=m;jj++) {
	    y=a[jj][j];
	    z=a[jj][i];
	    a[jj][j]=y*c+z*s;
	    a[jj][i]=z*c-y*s;
	  }
	}
	rv1[l]=0.0;//cout<<"rv1[l] "<<rv1[l]<<endl;
	rv1[k]=f;//cout<<"rv1[k] "<<rv1[k]<<endl;
	w[k]=x;//cout<<"w(k) "<<w[k]<<endl;
      }
    }
  }

// Singular value decompostion for C++
void svdcmp_cpp(Dmat& a, Dvec& w, Dmat& v)
  {
    unsigned int m,n;
    m = a.rows();
    n = a.cols();
    
    if(v.rows()!=n && v.cols()!=n) v.resize(n,n);
    if(w.size()!=n) w.resize(n);

    unsigned int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    Dvec rv1("",n);
    rv1=0.0;

    //cout<<"In svdcmp_cpp: a_matrix "<<endl;
    //for(unsigned int i=0;i<m;++i){
    //cout<<endl;
    //for(unsigned int j=0;j<n;++j)
    //cout<<a(i,j)<<"";}
    //cout<<endl;

    
    g=scale=anorm=0.0;
    for (i=0;i<n;i++) {
      l=i+1;
      rv1(i)=scale*g;
      g=s=scale=0.0;
      if (i <= m) {
	for (k=i;k<m;k++) scale += fabs(a(k,i));
	if (scale) {
	  for (k=i;k<m;k++) {
	    a(k,i) /= scale;
	    s += a(k,i)*a(k,i);
	  }
	  f=a(i,i);
	  g = -SIGN_C(sqrt(s),f);
	  h=f*g-s;
	  a(i,i)=f-g;
	  for (j=l;j<n;j++) {
	    for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
	    f=s/h;
	    for (k=i;k<m;k++) a(k,j) += f*a(k,i);
	  }
	  for (k=i;k<m;k++) a(k,i) *= scale;
	}
      }
      w(i)=scale *g;//cout<<"w(i) "<<w(i)<<endl;
      g=s=scale=0.0;
      if (i <= m && i != n) {
	for (k=l;k<n;k++) scale += fabs(a(i,k));
	if (scale) {
	  for (k=l;k<n;k++) {
	    a(i,k) /= scale;
	    s += a(i,k)*a(i,k);
	  }
	  f=a(i,l);
	  g = -SIGN_C(sqrt(s),f);
	  h=f*g-s;
	  a(i,l)=f-g;
	  for (k=l;k<n;k++) rv1(k)=a(i,k)/h;
	  for (j=l;j<m;j++) {
	    for (s=0.0,k=l;k<n;k++) s += a(j,k)*a(i,k);
	    for (k=l;k<n;k++) a(j,k) += s*rv1(k);
	  }
	  for (k=l;k<n;k++) a(i,k) *= scale;
	}
      }
      anorm=MAX(anorm,(fabs(w(i))+fabs(rv1(i))));
    }
    for (i=n-1;i>=0;i--) {
      if (i < n) {
	if (g) {
	  for (j=l;j<n;j++)
	    v(j,i)=(a(i,j)/a(i,l))/g;
	  for (j=l;j<n;j++) {
	    for (s=0.0,k=l;k<n;k++) s += a(i,k)*v(k,j);
	    for (k=l;k<n;k++) v(k,j) += s*v(k,i);
	  }
	}
	for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
      }
      v(i,i)=1.0;
      g=rv1(i);
      l=i;
      if(i==0) break;//do not subtract 1 from 0
    }
    for (i=MIN(m,n)-1;i>=0;i--) {
      l=i+1;
      g=w(i);
      for (j=l;j<n;j++) a(i,j)=0.0;
      if (g) {
	g=1.0/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
	  f=(s/a(i,i))*g;
	  for (k=i;k<m;k++) a(k,j) += f*a(k,i);
	}
	for (j=i;j<m;j++) a(j,i) *= g;
      } else for (j=i;j<m;j++) a(j,i)=0.0;
      ++a(i,i);
      if(i==0) break;//do not subtract 1 from 0
    }
    for (k=n-1;k>=0;k--) {
      for (its=0;its<30;its++) {
	flag=1;
	for (l=k;l>=0;l--) {//this index should start from k 
	  nm=l-1;
	  //debug
	  //cout<<"fabs(rv1) and anorm"<<fabs(rv1(l))+(float)anorm
	  //<<" "<<(float) anorm<<endl;
	  //if(nm<n) cout<<"fabs(w) + anorm "<< fabs(w(nm))+(float)anorm<<" "
	  //<<(float)anorm<<endl;
	  //cout<<"k and l "<<k<<" "<<l<<endl;
	  //if(l==0) cout<<"rv(0) "<<rv1(l)<<endl;//when l==0, rv(0)=0
	  if ((float)(fabs(rv1(l))+anorm) == (float)anorm) {
	    flag=0;
	    break;
	  }
	  if ((float)(fabs(w(nm))+anorm) == (float)anorm) break;
	}
	if (flag) {
	  c=0.0;
	  s=1.0;
	  for (i=l;i<k;i++) {
	    f=s*rv1(i);
	    rv1(i)=c*rv1(i);
	    if ((float)(fabs(f)+anorm) == (float)anorm) break;
	    g=w(i);
	    h=pythag(f,g);
	    w(i)=h;//cout<<"wi "<< w(i)<<endl;
	    h=1.0/h;
	    c=g*h;
	    s = -f*h;
	    for (j=0;j<m;j++) {
	      y=a(j,nm);
	      z=a(j,i);
	      a(j,nm)=y*c+z*s;
	      a(j,i)=z*c-y*s;
	    }
	  }
	}
	z=w(k);
	if (l == k) {
	  if (z < 0.0) {
	    w(k) = -z;//cout<<"w(k) "<<w(k)<<endl;
	    for (j=0;j<n;j++) v(j,k) = -v(j,k);
	  }
	  break;
	}
	if (its == 30) ErrorMessage("no convergenece").throwMe();
	 
	x=w(l);
	nm=k-1;
	y=w(nm);
	g=rv1(nm);
	h=rv1(k);
	f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	g=pythag(f,1.0);
	f=((x-z)*(x+z)+h*((y/(f+SIGN_C(g,f)))-h))/x;
	c=s=1.0;
	for (j=l;j<=nm;j++) {
	  i=j+1;
	  g=rv1(i);
	  y=w(i);
	  h=s*g;
	  g=c*g;
	  z=pythag(f,h);
	  rv1(j)=z;
	  c=f/z;
	  s=h/z;
	  f=x*c+g*s;
	  g = g*c-x*s;
	  h=y*s;
	  y *= c;
	  for (jj=0;jj<n;jj++) {
	    x=v(jj,j);
	    z=v(jj,i);
	    v(jj,j)=x*c+z*s;
	    v(jj,i)=z*c-x*s;
	  }
	  z=pythag(f,h);
	  w(j)=z;//cout<<"w(j) "<<w(j)<<endl;
	  if (z) {
	    z=1.0/z;
	    c=f*z;
	    s=h*z;
	  }
	  f=c*g+s*y;
	  x=c*y-s*g;
	  for (jj=0;jj<m;jj++) {
	    y=a(jj,j);
	    z=a(jj,i);
	    a(jj,j)=y*c+z*s;
	    a(jj,i)=z*c-y*s;
	  }
	}
	rv1(l)=0.0;//cout<<"rv1(l) "<<rv1(l)<<endl;
	rv1(k)=f;//cout<<"rv1(k) "<<rv1(k)<<endl;
	w(k)=x;//cout<<"w(k) in the end "<<w(k)<<endl;
      }
      if(k==0) break;//do not subtract 1 from 0
    }


    /*
    //debug
    cout<<"Out svdcmp_cpp: a_matrix "<<endl;
    for(unsigned int i=0;i<m;++i){
      cout<<endl;
      for(unsigned int j=0;j<n;++j)
	cout<<a(i,j)<<"";}
    cout<<endl;

    //debug
    cout<<"Out svdcmp_cpp w-array "<< w<<endl;
    
    //debug
    cout<<"Out svdcmp: v_matrix "<<endl;
    for(unsigned int i=0;i<m;++i){
      cout<<endl;
      for(unsigned int j=0;j<n;++j)
	cout<<v(i,j)<<"";}
    cout<<endl;
    */
  }


//PYTHAG
 double  pythag(const double& a,
		const double& b)
{
 
  double  absa,absb;
  absa = fabs(a);
  absb = fabs(b);
  double return_value;

      if(absa > absb)
        return_value = absa*sqrt(1.0+(absb/absa)*(absb/absa));
      else{
        if(absb == 0.0)
	  return_value = 0.0;
        else
          return_value = absb*sqrt(1.0+(absa/absb)*(absa/absb));
      }
      
      return(return_value);
}


//Matrix equation solve (U w V^T) x = b, answer stored into x
void  svbksb(Dmat& u,
	     Dvec& w,
	     Dmat& v,
	     Dvec& b,
	     Dvec& x)
{//solve u * x = b: u is of m x n, b: size of m, x size of n
  unsigned int m,n;
  m = u.rows();
  n= u.cols();
  if(w.size()!=n) ErrorMessage("svbksb:: size of w should be n").throwMe();
  if(v.rows()!=n && v.cols()!=n) ErrorMessage("svbksb:: size of v should be nxn").throwMe();
  if(b.size()!=m) ErrorMessage("svbksb:: size of b should be m").throwMe(); 
  if(x.size()!=n) x.resize(n);//output size n

  //solving equation
  // Ax = b , where A(MxN) = U(MxN) w(N) V^T(NxN)
  // x = V (1/w)  (U^T b)
  unsigned int i,j,jj;
  
  double s;
  Dvec tmp("",n);

  for (j=0;j<n;j++) {
    s=0.0;
    if (w(j)) {
      for (i=0;i<m;i++) 
	s += u(i,j)*b(i);
      s /= w(j);
    }
    tmp(j)=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v(j,jj)*tmp(jj);
    x(j)=s;
  }
}

//SIGN_C: fron fortran SIGN function
double SIGN_C(const double& a, const double& b) 
 {
   double return_value;
   if(b >=0.0) return_value=fabs(a);
   else return_value = -fabs(a);
   return(return_value);
 }





//------------------------------------
//Fix Cassini Sab counter and S/C SCLK counter
//------------------------------------
void fixCassiniSabSclk(map<unsigned int, unsigned int>& full_path_sab_counter,
		       map<unsigned int, unsigned int>& full_path_sclk,
		       unsigned int& first_good_record,
		       unsigned int& num_sab_counter_roll_over,
		       const bool& plot_option)
  {
    unsigned int Nmax_sab_in_16bits = 65536;//sab counter is a 16 bit word


    //------------------------
    //check the container size
    //------------------------
    if(full_path_sab_counter.size() != full_path_sclk.size())
      ErrorMessage("sab_counter record size is different from sclk record size: "+toStr(full_path_sab_counter.size()) +" "+toStr(full_path_sclk.size())).throwMe();
    unsigned int complete_sab_count=full_path_sab_counter.size();

    //-----------------------------------------
    //debugging plot: if plot option is true
    //------------------------------------------
    vector<Uvar> y1,y2,y3;
    y1.clear();y2.clear();y3.clear();
    //store data for plot
    if(plot_option){
      for(unsigned int i=0;i<complete_sab_count;++i)
	y1.push_back(full_path_sab_counter[i]);
    }

    //--------------------------
    //(0)  Find the first good sab counter
    //---------------------------
    Array1D<unsigned int> r_buffer("",5);
    unsigned int start_sab_index=0;//default start value
    bool start_burst_determined=false;
    while(!start_burst_determined){
      for(unsigned int ii=0;ii<5;++ii){
	if( (start_sab_index+ii)< complete_sab_count)
	  r_buffer(ii)=full_path_sab_counter[start_sab_index+ii];
	else{
	  ErrorMessage("can not find a good starting sab index");
	}
      }
      bool good_buffer= true;
      for(unsigned int ii=0;ii<4;++ii){
	if(r_buffer(ii) >=  r_buffer(ii+1)) good_buffer=false; 
      }
      if(good_buffer){
	start_burst_determined=true;
      }
      else{
	start_burst_determined=false;
	start_sab_index++;//test the next one
      }
    }
    cout<<"starting sab index after checking sab counter "<< start_sab_index<<endl;
  
    //--------------------------
    //(1)  If there is a bad sab before the first good sab
    // fix it: reduce burst number one by one
    //---------------------------
    if(start_sab_index!=0){ 
      int  n_ref = full_path_sab_counter[start_sab_index];
      for( int ii=start_sab_index-1; ii>= 0;ii--){
	int  n = n_ref - ((int)start_sab_index - ii);
	if(n<1) {
	  cout<<"sab index will be set to 0, did not fix "<< endl;
	  full_path_sab_counter[(unsigned int) ii]=0;
	}
	else{
	  cout<<"New sab record will be written"<<endl;
	  full_path_sab_counter[(unsigned int)ii]=n;
	}
	cout<<"adjusted sab counter "<<ii<<" "<< n<<endl;
      }
    }


    //----------------------------------------------------
    //(2): fix any SAB bit error after start_sab_index: 
    // not fixing sab roll over
    //----------------------------------------------------
    r_buffer.resize(3);//change to ring buffer size of 3
    for(unsigned int ii=start_sab_index+1;ii<(complete_sab_count-1);++ii){
      r_buffer(0)= full_path_sab_counter[ii-1];
      r_buffer(1)= full_path_sab_counter[ii];
      r_buffer(2)= full_path_sab_counter[ii+1];
      if((r_buffer(0) < r_buffer(1)) &&
	(r_buffer(1) < r_buffer(2))) {}
      else if(r_buffer(1)> r_buffer(0) &&
	      r_buffer(1)>= r_buffer(2) &&
	      r_buffer(0) <r_buffer(2)){
	// 1-3-2 fix
	cout<<"(Up) bit error occurs at comple sab index of "<<ii<<endl;
	cout<<"ring buffer "<< r_buffer<<endl;
	cout<<"old sab counter "<< r_buffer(1)<<endl;
	full_path_sab_counter[ii]=(r_buffer(0)+r_buffer(2))/2;
	cout<<"new sab counter "<<full_path_sab_counter[ii]<<endl;
      }
      else if(r_buffer(1) <= r_buffer(0) &&
	      r_buffer(1) < r_buffer(2) &&
	      r_buffer(0) < r_buffer(2)){
	//3-1-2 fix
	cout<<"(Down)bit error occurs at comple sab index of "<<ii<<endl;
	cout<<"ring buffer "<< r_buffer<<endl;
	cout<<"old sab counter "<< r_buffer(1)<<endl;
	full_path_sab_counter[ii]=(r_buffer(0)+r_buffer(2))/2;
	cout<<"new sab counter "<<full_path_sab_counter[ii]<<endl;
      }
      else if((r_buffer(1) < r_buffer(0)) &&
	      (r_buffer(1) < r_buffer(2)) &&
	      (r_buffer(2) < r_buffer(0))){
	cout<<"Detect roll over, but handled later"<<endl;
      }
      else{
	cout<<"Unusual sab numbers, possibly right before roll over "<<r_buffer<<endl;
      }
    }


    //--------------------------------------------
    //(2-1) Last  record fix only if the last record's
    // sab counter decreases
    //--------------------------------------------
    if(full_path_sab_counter[complete_sab_count-1] 
       <= full_path_sab_counter[complete_sab_count-2]){
      cout<<"Last sab counter has a bit error "<<endl;
      full_path_sab_counter[complete_sab_count-1]=
	(full_path_sab_counter[complete_sab_count-1]+1)%Nmax_sab_in_16bits;
    }
    
    //store data for plot
    if(plot_option){
      for(unsigned int i=0;i<complete_sab_count;++i)
	y2.push_back(full_path_sab_counter[i]);
    }



    //--------------------------
    //(3) SAB counter roll over 
    //--------------------------
    num_sab_counter_roll_over= 0;//reset to 0
    for(unsigned int ii=1;ii<complete_sab_count-1;++ii){
      r_buffer(0)= full_path_sab_counter[ii-1];
      r_buffer(1)= full_path_sab_counter[ii];
      r_buffer(2)= full_path_sab_counter[ii+1];
      if((r_buffer(1) < r_buffer(0)) &&
	 (r_buffer(1) < r_buffer(2)) &&
	 (r_buffer(2) < r_buffer(0))){
	cout<<"sab counter roll over"<<endl;
	num_sab_counter_roll_over++;
	for(unsigned int jj=ii;jj<complete_sab_count;++jj)
	  full_path_sab_counter[jj]+= Nmax_sab_in_16bits;
      }
    }

    //-----------------------------
    //(3-1)take care of the last elements
    //------------------------------  
    if(full_path_sab_counter[complete_sab_count-1]
       <full_path_sab_counter[complete_sab_count-2]){
      cout<<"sab count roll over at the last element or bit error"<<endl;
      num_sab_counter_roll_over++;
      full_path_sab_counter[complete_sab_count-1] += Nmax_sab_in_16bits;
    }


    //----------------------
    //(3-2) Finally, check the whole record
    // after bit correction and unwrapp roll over
    //----------------------
    for(unsigned int ii=1;ii<complete_sab_count-1;++ii){
      r_buffer(0)= full_path_sab_counter[ii-1];
      r_buffer(1)= full_path_sab_counter[ii];
      r_buffer(2)= full_path_sab_counter[ii+1];
      if(!( (r_buffer(0) < r_buffer(1)) &&
	    (r_buffer(1) < r_buffer(2)))){
	//fix r_buffer(1)
	cout<<"Fix sab counter after unwrapping sab counter "<<endl;
	cout<<"Problem series "<< r_buffer<<endl;
	full_path_sab_counter[ii]=(r_buffer(0)+r_buffer(2))/2;
      }
    }

    //store data for plot
    if(plot_option){
      for(unsigned int i=0;i<complete_sab_count;++i)
	y3.push_back(full_path_sab_counter[i]);
    }


    //--------------------------
    //plot
    //------------------------
    if(plot_option){
      Plot a;
      a.addY(y1,"",line("solid","black",1),sym("none"));
      a.addLegend("original sab counter");
      a.addY(y2,"",line("solid","red",1),sym("none"));
      a.addLegend("bit correction");
      a.addY(y3,"",line("solid","green",1),sym("none"));
      a.addLegend("unwrapped");
      a.show("x");
      //clear the container
      y1.clear();y2.clear();y3.clear();
    }
    

    //store sclk for plot
    if(plot_option){
      for(unsigned int i=0;i<complete_sab_count;++i)
	y1.push_back((int)full_path_sclk[i]-(int)full_path_sclk[0]);
    }

    //-----------------------
    //(4) Find the first good record
    //----------------------------
    start_burst_determined=false;
    start_sab_index=0;
    r_buffer.resize(5);
    while(!start_burst_determined){
      for(unsigned int ii=0;ii<5;++ii){
	if( (start_sab_index+ii)< complete_sab_count)
	  r_buffer(ii)=full_path_sclk[start_sab_index+ii];
	else{
	  ErrorMessage("can not find a good starting sab index");
	}
      }
      bool good_buffer= true;
      for(unsigned int jj=0;jj<4;++jj){
	if(r_buffer(jj)>r_buffer(jj+1)) good_buffer=false; 
      }
      if(good_buffer){
	start_burst_determined=true;
      }
      else{
	start_burst_determined=false;
	start_sab_index++;//test the next one
      }
    }
    cout<<"First good record  after checking sclk "<< start_sab_index<<endl;
    first_good_record= start_sab_index;

    //----------------------------------------
    //(5) Correct: before the first good record 
    // Not implemented, Unlike SAB, I do not want
    // to assign possible wrong time there
    //------------------------------------------

    //---------------------------------
    //(6) Fix bit  error corruption 
    //-------------------------------
    r_buffer.resize(3);
    for(unsigned ii=start_sab_index;ii<(complete_sab_count-1);++ii){
      r_buffer(0)= full_path_sclk[ii-1];
      r_buffer(1)= full_path_sclk[ii];
      r_buffer(2)= full_path_sclk[ii+1];
      if( (r_buffer(0) < r_buffer(1)) &&
	  (r_buffer(2) < r_buffer(1)) &&
	  (r_buffer(0) <= r_buffer(2))){
	cout<<"(Up) bit corruption  at "<< ii <<endl;
	cout<<"original sequence(middle one is a problem "<<r_buffer<<endl;
	full_path_sclk[ii]=(r_buffer(0)+r_buffer(2))/2;
	cout<<"fixed one "<< full_path_sclk[ii]<<endl;
      }
      else if( (r_buffer(0) > r_buffer(1))&&
	       (r_buffer(2) > r_buffer(1))&&
	       (r_buffer(2) >= r_buffer(0))){
	cout<<"(Down) bit corruption at"<<ii<<endl;
	cout<<"original sequence(middle one is a problem "<<r_buffer<<endl;
	full_path_sclk[ii]=(r_buffer(0)+r_buffer(2))/2;
	cout<<"fixed one "<< full_path_sclk[ii]<<endl;
      }
      else{}
    }
    
    
    
    //------------------------
    //(7) Take care of the last element
    //------------------------
    if(full_path_sclk[complete_sab_count-1] 
       < full_path_sclk[complete_sab_count-2]){
      cout<<"the last element's sclk is corrupted "<<endl;
      cout<<"addition time correspond to the time difference of previous"<<endl;
      cout<<"two good records will be added "<<endl;
      full_path_sclk[complete_sab_count-1]=full_path_sclk[complete_sab_count-2];
      int i_diff=full_path_sclk[complete_sab_count-2]-
	full_path_sclk[complete_sab_count-3];
      if(i_diff<0) {
	cout<<"last element time correction failed "<<endl;
	i_diff=0;
      }
      full_path_sclk[complete_sab_count-1]+=i_diff;
    }

    //store sclk for plot
    if(plot_option){
      for(unsigned int i=0;i<complete_sab_count;++i)
	y2.push_back((int)full_path_sclk[i]-(int)full_path_sclk[0]);
    }
    
    //plot sclk
    if(plot_option){
      Plot b;
      b.addY(y1,"",line("solid","black",1),sym("none"));
      b.addLegend("original sclk");
      b.addY(y2,"",line("solid","red",1),sym("none"));
      b.addLegend("bit correction");
      b.show("x");
    }
  }


//-------------------------
//fix fin number roll over
//-----------------------
void fixFinNumberRollOver(map<unsigned int, unsigned int>& full_path_fin,
			  const unsigned int& number){
  unsigned int Nsize=full_path_fin.size();
  map<unsigned int, unsigned int>::const_iterator q;
  unsigned int number_of_roll_over=0;
 
  map<unsigned int, unsigned int> tmp_fin;


  for(map<unsigned int, unsigned int>::const_iterator p = full_path_fin.begin();
      p!=full_path_fin.end();++p){
    if(p==full_path_fin.begin()){
      tmp_fin[p->first]= p->second;
      continue;//do not bother with the first one
    }
    q = p;
    q--;//go to previous record
    if(q->second > p->second){
      number_of_roll_over++;
      cout<<"fin number roll over "<< q->second<<" "<<p->second<<endl;
    }
    tmp_fin[p->first]= p->second + number_of_roll_over * number;
  }
  cout<<"number of roll overs "<< number_of_roll_over<<endl; 
  if(tmp_fin.size()!= full_path_fin.size()) ErrorMessage("Utils.cpp::two fin containers have different sizes").throwMe();
  

  q=tmp_fin.begin();
  for(map<unsigned int, unsigned int>::const_iterator p = full_path_fin.begin();
     p!=full_path_fin.end() && q!=tmp_fin.end();++p, ++q){
    full_path_fin[p->first]= q->second;
  }
  
  if(Nsize!= full_path_fin.size()) ErrorMessage("Utils.cpp: original container size has changed").throwMe();
}


//-----------------
//fix clock error
//---------------
void  computeSCLKandBRST(map<unsigned int, unsigned int>& full_path_sab_counter,
			 map<unsigned int, unsigned int>& full_path_sclk,
			 map<unsigned int, unsigned int>& full_path_fin,
			 map<unsigned int, unsigned int>& full_path_bii,
			 map<unsigned int, unsigned int>& full_path_ctbc,
			 map<unsigned int,  Uvar>& full_path_brst,
			 map<unsigned int,  Uvar>& full_path_bpd,
			 unsigned int& sclk_correction,
			 unsigned int& brst_correction,
			 unsigned int& bpd_correction,
			 unsigned int& bii_correction){
 
  unsigned int Nmax=full_path_sab_counter.size();
  if(!(Nmax == full_path_fin.size() &&
       Nmax == full_path_bii.size() &&
       Nmax == full_path_ctbc.size() &&
       Nmax == full_path_brst.size() &&
       Nmax == full_path_bpd.size())){
    cout<<"Error: container size mismatch "<<endl;
    cout<<"expected container size "<< Nmax<<endl;
    cout<< full_path_fin.size()<<" "<< full_path_bii.size()<<endl;
    cout<< full_path_ctbc.size()<<" "<<full_path_brst.size()<<" "<< full_path_bpd.size()<<endl;
    ErrorMessage("Utils.cpp: computeSclkandBRST: container sizes do not match").throwMe();
  }


  sclk_correction=0;
  brst_correction=0;
  bpd_correction=0;

 


  //-----------------------
  //check fin number is correctly unwrapped
  //----------------------
  unsigned int Nfin_change=0;
  vector<unsigned int> fin_rec;
  map<unsigned int, unsigned int>::const_iterator p=full_path_fin.begin();
  fin_rec.push_back(p->second);
  
  for(map<unsigned int, unsigned int>::const_iterator q=full_path_fin.begin();
      q != full_path_fin.end();++q){
    if(q==full_path_fin.begin()) continue;
    p = q;
    p--;
    if(p->second > q->second) ErrorMessage("Fin number is not correctly unwrapped ").throwMe();
    if(p->second != q->second){
      Nfin_change++;
      fin_rec.push_back(q->second);
    }
  }
  cout<<"total number of fin number change "<< Nfin_change<<endl;
  unsigned int Nfin_max=fin_rec.size();
  cout<<"total number of  fins "<< Nfin_max<<endl;

  //--------------------------------
  //check sab number is increasing with burst
  //--------------------------------
  for(unsigned int ii=0;ii<Nmax-1;++ii){
    if(full_path_sab_counter[ii+1]< full_path_sab_counter[ii]){
      cout<<" current sab number "<< full_path_sab_counter[ii]<<endl;
      cout<<"  next  sab number "<< full_path_sab_counter[ii+1]<<endl;
      ErrorMessage("sab number does not increase with burst number ").throwMe();
    }
  }

  //-----------------------------
  //check sclk is increasing with sab number
  //------------------------------
  for(unsigned int ii=0;ii<Nmax-1;++ii){
    if(full_path_sclk[ii+1]< full_path_sclk[ii]){
      cout<<" current sclk "<< full_path_sclk[ii]<<endl;
      cout<<" next sclk "<< full_path_sclk[ii+1]<<endl;
      ErrorMessage("sclk does not increase with burst number ").throwMe();
    }
  }
  

  //-----------------------------
  //MAIN LOOP: make sclk and brst correction
  //------------------------------
  //Nfin: total number of different fin numbers
  // different elements are stored in fin_rec
  //-----------------------------------------

  for(unsigned int i=0;i<Nfin_max;++i){
    unsigned int fin_ref = fin_rec[i];
    vector<unsigned int> sab, sclk,bii, ctbc;
    vector<double> bpd, brst;
    sab.clear(); sclk.clear(); bii.clear(); ctbc.clear();
    bpd.clear();brst.clear();
   
    for(unsigned int j=0;j<Nmax;++j){
      if(full_path_fin[j]==fin_ref){
	sab.push_back(full_path_sab_counter[j]);
	sclk.push_back(full_path_sclk[j]);
	bii.push_back(full_path_bii[j]);
	ctbc.push_back(full_path_ctbc[j]);
	bpd.push_back(full_path_bpd[j].getInUnits("s"));
	brst.push_back(full_path_brst[j].getInUnits("s"));
      }//collect same fin record
    }//loop over all the records
    
    unsigned int NN=sab.size();
    if(NN != sclk.size() ||
       NN != bii.size()  ||
       NN != ctbc.size() ||
       NN != bpd.size()  ||
       NN != brst.size())
      ErrorMessage("Utils.cpp: container sizes for same fin number are different").throwMe();
    
    if(NN<=1){
      cout<<"Only one record for the fin number of "<< fin_ref<<endl;
      continue;//No need to do anything
    }

    
    //--------------------
    //ctbc start number: bii -1
    // ctbc start number: 255 when bii=0
    // time lapse from first record: [ctbc(0) - current_ctbc]*bpd 
    //--------------------
    bool fix_ctbc=false;
    for(unsigned int ii=0;ii<NN-1;++ii){
      if((ctbc[ii] > ctbc[ii+1]) && bii[ii]!=0) {
	fix_ctbc=true;
	cout<<"ctbc current and future "<< ctbc[ii]<<" "<< ctbc[ii+1]<<endl;
	if(ii!=0)	ctbc[ii]= ctbc[ii-1]-1;
	cout<<"after fix "<< ctbc[ii]<<" "<<ctbc[ii+1]<<endl;
      }
    }

    if(fix_ctbc==true){
      for(unsigned int ii=0;ii<NN-1;++ii){
	if((ctbc[ii] > ctbc[ii+1]) && bii[ii]!=0) {
	  ErrorMessage("could not fix ctbc ").throwMe();
	}
      }
    }
   

    //----------------
    //fix possible bpd error
    //If  any two nearby bursts have the same bpd
    // then take it as good (valid) bpd
    //-----------------
    bool found_good_bpd=false;
    double good_bpd;
   
    for(unsigned int ii=0;ii<NN-1;++ii){
      if(bpd[ii]==bpd[ii+1]){
	good_bpd = bpd[ii];
	found_good_bpd=true;
	break;
      }
    }
    if(!found_good_bpd) ErrorMessage("utils.cpp: for the same fin, any two nearby bpds are different").throwMe();

   
    for(unsigned int ii=0;ii<NN;++ii){
      if(bpd[ii]!= good_bpd){

	cout<<"bpd correction "<<endl;
	cout<<"original record "<< bpd[ii]<<endl;
	cout<<"new record "<< good_bpd<<endl;

	//replace container
	bpd[ii]=good_bpd;
	bpd_correction++;
	
	//replace original record
	for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	  if(p->second == sab[ii]) full_path_bpd[p->first] = Uvar(good_bpd,"s"); //replace
      }//if we found any
    }//check all bpd records


    //---------------------
    //fix possible bii error
    //-------------------------
    bool found_good_bii=false;
    unsigned int good_bii;
    for(unsigned int ii=0;ii<NN-1;++ii){
      if(bii[ii]==bii[ii+1]){
	good_bii = bii[ii];
	found_good_bii=true;
	break;
      }
    }
    if(!found_good_bii) ErrorMessage("Utils.cpp: for the same fin, could not find valid bii").throwMe();
    for(unsigned int ii=0;ii<NN;++ii){
      if(bii[ii]!= good_bii){

	cout<<"original bii "<<bii[ii]<<endl;
	cout<<"new bii record "<< good_bii<<endl;
	//replace container
	bii[ii]=good_bii;
	bii_correction++;

	//replace original bii
	for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	  if(p->second == sab[ii]) full_path_bii[p->first] = good_bii; //replace
      }
    }
      
    //-----------------------------------------------------------------------
    //find first good record by comparing time difference=bpd * sab_difference
    //---------------------------------------------------------------------
    unsigned int first_good_record;
    bool found_first_good_record=false;
    for(unsigned int ii=0;ii<NN-1;++ii){
      
      //-------------------------------
      //if the first record is not missing
      //-----------------------------
      if(bii[0]!=0){
	if(ii==0 &&  ctbc[0]==(bii[0]-1)){//this is the first record
	  //First record exists
	  first_good_record=0;
	  found_first_good_record=true;
	}
      }
      else if(bii[0]==0){
	if(ii==0 && ctbc[0]==255){//this is the first record
	  first_good_record=0;
	  found_first_good_record=true;
	}
      }
      else{}

      if(found_first_good_record) break;

      
      if(sclk[ii] > sclk[ii+1]){
	cout<<"current sclk "<< sclk[ii]<<endl;
	cout<<"next sclk "<< sclk[ii+1]<<endl;
	ErrorMessage("Utils.cpp::sclk is not increasing with burst number").throwMe();
      }
      double  diff_sclk = (double) (sclk[ii+1]- sclk[ii]);
      diff_sclk = diff_sclk* 4.0* 1000.0;//resolution is 0.25 ms
      double diff_brst=  brst[ii+1]- brst[ii];
      diff_brst = diff_brst*4.0*1000.0;
      
      unsigned int diff = (unsigned int) round_double(diff_sclk + diff_brst);
      unsigned int bpd_in_units_of_tick= (unsigned int) round_double(good_bpd*4.0*1000.0);
      
      if(diff==bpd_in_units_of_tick){
	first_good_record=ii;
	found_first_good_record=true;
	break;
      }
      else{
	cout<<"current sclk "<< sclk[ii]  <<" "<<brst[ii]<<endl;
	cout<<"next record"  << sclk[ii+1]<<" "<<brst[ii+1]<<endl;
	cout<<"diff in ticks  between two nearby bursts "<< diff<<endl;
	cout<<"valid bpd "<< good_bpd<<" "<<endl;
	cout<<"bpd in units of ticks "<< bpd_in_units_of_tick<<endl;
	cout<<"bii "<<bii[0]<<endl;
	cout<<endl;
      }
    }

    //-------------
    //if there is no good record
    //-----------------
    if(!found_first_good_record) ErrorMessage("none of two nearby are separated by bpd ").throwMe();

    //------------------
    //print record before good one
    //--------------------
    if(first_good_record!=0){
      for(unsigned int ii=0;ii<first_good_record;++ii){
	cout<<"first good record is not the first command after instruction upload "<<endl;
	cout<<"sclk brst  "<< sclk[ii]<<" "<< brst[ii]<<" "<<endl;
	cout<<"bpd in ms  "<< bpd[ii]<<" "<< good_bpd<<endl;
	cout<<endl;
      }
    }


    //--------------------------------------------------
    //correct sclk brst before good record  (sclk_correction, brst_correction)
    //---------------------------------------------------
    if(first_good_record != 0){
      //need correction before good record
      for(unsigned int jj=0;jj<first_good_record;++jj){

	//set bool
	bool correction_needed=false;

	//expected difference based on bpd
	unsigned int sab_diff= sab[first_good_record]- sab[jj];//positive 
	double burst_time_diff = good_bpd* sab_diff;//positive
	unsigned int  expected_tick_diff= (unsigned int) round_double(burst_time_diff*4000.0);
	
	//measured difference based on sclk and bpd
	unsigned int sclk_diff = sclk[first_good_record]- sclk[jj];
	int sclk_diff_tick= 4000 * sclk_diff;
	double brst_diff = brst[first_good_record] - brst[jj];
	int brst_diff_tick= (int) round_double( brst_diff*4000.0);
	if( (sclk_diff_tick +brst_diff_tick) < 0 ) ErrorMessage("Utils.cpp:should not happen").throwMe();
	unsigned int  sab_tick_diff =(unsigned int)( sclk_diff_tick + brst_diff_tick);

	unsigned int start_sclk=sclk[jj] ;
	double start_brst= brst[jj];

	if(expected_tick_diff > sab_tick_diff){
	  correction_needed=true;
	  //move jj burst backward
	  int tick_diff =  (int)expected_tick_diff - (int)sab_tick_diff;
	  for(int c=0;c<tick_diff;++c){
	    start_brst -= 0.00025;//1 tick is 0.25 ms
	    if(start_brst<0){
	      start_sclk--;
	      start_brst += 1.0;
	    }
	  }
	}
	else if(expected_tick_diff < sab_tick_diff){
	  correction_needed=true;
	  //move jj burst forward
	  int tick_diff = (int)sab_tick_diff - (int)expected_tick_diff;
	  for(int c=0;c<tick_diff;++c){
	    start_brst += 0.00025;
	    if(start_brst>1.0){
	      start_sclk++;
	      start_brst -=1.0;
	    }
	  }
	}
	else{}//do nothing
	
	if(correction_needed){
	  //cout<<"reverse correction before first good record "<<endl;
	  //cout<<"bpd in ticks "<< good_bpd*4000.0<<endl;
	  //cout<<"sab diff "<< sab_diff<<endl;
	  //cout<<"expected tick diff "<< expected_tick_diff<<endl;
	  //cout<<"sab based tick diff "<< sab_tick_diff<<endl;
	  //cout<<"sclk correction "<<endl;
	  //cout<<"old sclk and new sclk "<< sclk[jj]<<" "<<start_sclk<<endl;
	  //cout<<"brst correction "<<endl;
	  //cout<<"old brst and new brst  brst "<<brst[jj]<<" "<< start_brst<<endl;
	  //cout<<endl;
	  if(start_sclk != sclk[jj]) {
	    sclk_correction++;
	    for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	      if(p->second == sab[jj]) full_path_sclk[p->first] = start_sclk; //replace
	  }
	  if(start_brst != brst[jj]){
	    brst_correction++;
	    for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	      if(p->second == sab[jj]) full_path_brst[p->first] = Uvar(start_brst,"s"); //replace
	  }
	}//when correction is needed
      }//make correction before good record
    }//when first good record is not the first command
    
    
    //--------------------
    //fix sclk after good record
    //----------------------
    for(unsigned int jj=first_good_record+1;jj<NN;++jj){
      bool correction_needed=false;
 
      //expected difference based on bpd
      unsigned int sab_diff= sab[jj]- sab[first_good_record];//positive 
      double burst_time_diff = good_bpd* sab_diff;//positive
      unsigned int  expected_tick_diff= (unsigned int) round_double(burst_time_diff*4000.0);
      
      //measured difference based on sclk and bpd
      unsigned int sclk_diff = sclk[jj]- sclk[first_good_record];
      int sclk_diff_tick= 4000 * sclk_diff;
      double brst_diff = brst[jj] - brst[first_good_record];
      int brst_diff_tick= (int) round_double( brst_diff*4000.0);
      if( (sclk_diff_tick +brst_diff_tick) < 0 ) ErrorMessage("Utils.cpp:should not happen").throwMe();
      unsigned int  sab_tick_diff =(unsigned int)( sclk_diff_tick + brst_diff_tick);

      unsigned int start_sclk=sclk[jj] ;
      double start_brst= brst[jj];

      if(expected_tick_diff > sab_tick_diff){
	correction_needed=true;
	//move jj burst forward
	int tick_diff =  (int)expected_tick_diff - (int)sab_tick_diff;//supposed to be very small
	for(int c=0;c<tick_diff;++c){
	  start_brst +=0.00025;
	  if(start_brst>1.0){
	    start_sclk++;
	    start_brst -= 1.0;
	  }
	}
      }
      
      else if(expected_tick_diff < sab_tick_diff){
	correction_needed=true;
	//move jj burst backward
	int tick_diff = (int)sab_tick_diff - (int)expected_tick_diff;
	for(int c=0;c<tick_diff;++c){
	  start_brst -= 0.00025;
	  if(start_brst<0){
	    start_sclk--;
	    start_brst +=1.0;
	  }
	}
      }
      else{}//do nothing
      
      if(correction_needed){
	//cout<<"reference sclk brst "<< sclk[first_good_record]<<" "<< brst[first_good_record]<<endl;
	//cout<<"bpd in ticks "<< good_bpd*4000.0<<endl;
	//cout<<"sab diff "<< sab_diff<<endl;
	//cout<<"expected tick diff "<< expected_tick_diff<<endl;
	//cout<<"sab based tick diff "<< sab_tick_diff<<endl;
	//cout<<"sclk correction "<<endl;
	//cout<<"old sclk and new sclk "<< sclk[jj]<<" "<<start_sclk<<endl;
	//cout<<"brst correction "<<endl;
	//cout<<"old brst and new brst  brst "<<brst[jj]<<" "<< start_brst<<endl;
	//cout<<endl;
	if(start_sclk != sclk[jj]) {
	  sclk_correction++;
	  for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	    if(p->second == sab[jj]) full_path_sclk[p->first] = start_sclk; //replace
	}
	if(start_brst != brst[jj]){
	  brst_correction++;
	  for(p=full_path_sab_counter.begin();  p!= full_path_sab_counter.end();++p)
	    if(p->second == sab[jj]) full_path_brst[p->first] = Uvar(start_brst,"s"); //replace
	}
      }//when correction is needed
    }//make correction after first good record
  }//Loop over each fin number
}


//------------------------------------------
//compute burst time using
// trigger time and time from initial trigger
//
//-------------------------------------------
void  computeBurstTimeUsingTFI(map<unsigned int, unsigned int>& full_path_sab_counter,
			       map<unsigned int, unsigned int>& full_path_sclk,
			       map<unsigned int, unsigned int>& full_path_fin,
			       map<unsigned int, unsigned int>& full_path_bii,
			       map<unsigned int, unsigned int>& full_path_ctbc,
			       map<unsigned int,  Uvar>& full_path_brst,
			       map<unsigned int,  Uvar>& full_path_bpd,
			       map<unsigned int, Uvar>& full_path_tfi,
			       map<unsigned int, unsigned int>& full_path_trigger_time,
			       unsigned int& sclk_correction,
			       unsigned int& brst_correction,
			       unsigned int& bpd_correction,
			       unsigned int& bii_correction,
			       unsigned int& tfi_correction,
			       const int& trigger_time_adjust,
			       const unsigned int& show_text)
  {

    //check container size
    unsigned int Nmax=full_path_sab_counter.size();
    if(!(Nmax == full_path_fin.size() &&
	 Nmax == full_path_bii.size() &&
	 Nmax == full_path_ctbc.size() &&
	 Nmax == full_path_brst.size() &&
	 Nmax == full_path_bpd.size() &&
	 Nmax == full_path_tfi.size())){
      cout<<"Error: container size mismatch "<<endl;
      cout<<"expected container size "<< Nmax<<endl;
      cout<< full_path_fin.size()<<" "<< full_path_bii.size()<<endl;
      cout<< full_path_ctbc.size()<<" "<<full_path_brst.size()<<" "<< full_path_bpd.size()<<endl;
      cout<<full_path_tfi.size()<<endl;
      ErrorMessage("Utils.cpp: computeSclkandBRST: container sizes do not match").throwMe();
    }
    
    //reset correction counter
    sclk_correction=0;
    brst_correction=0;
    bpd_correction=0;
    bii_correction=0;
    tfi_correction=0;
    map<unsigned int, unsigned int>::const_iterator p;
    
    //--------------------------------------------
    //0: first possible trigger time bit corruption
    // remember: there is only one trigger time
    //---------------------------------------------
    unsigned int Ntrigger= full_path_trigger_time.size();
    cout<<"number of trigger times collected "<<Ntrigger<<endl;
    if(Ntrigger<10) ErrorMessage("number of trigger times collected less than 10: something suspicious").throwMe();
    
    //first element
    unsigned int trigger_time = full_path_trigger_time.begin()->second;;
    //compare with second element
    cout<<"first trigger time "<< trigger_time<<endl;
    if(trigger_time != (full_path_trigger_time.begin()++)->second){
      cout<<"first two trigger times do not match"<<endl;
      cout<<"trigger times "<< full_path_trigger_time[0]<<" "<< full_path_trigger_time[1]<<endl;
      cout<<"Possible bit corruption "<<endl;
      cout<<"Trigger time bit correction "<<endl;
      for(map<unsigned int, unsigned int>::const_iterator q=full_path_trigger_time.begin(); 
	  q!=full_path_trigger_time.end();++q){
	if(q==full_path_trigger_time.begin()) continue;
	p = q;
	p--;
	if(p->second ==q->second){
	  trigger_time = p->second;
	  cout<<"corrected trigger time "<< trigger_time<<endl;
	  break;
	}//if found, break
      }//search for good trigger
    }//if two nearby trigger time not correct
    else{
      cout<<"first two trigger time "<< full_path_trigger_time.begin()->second<<" "<<(full_path_trigger_time.begin()++)->second<<endl;
    }
    cout<<"trigger time from EGSE file "<< trigger_time<<endl;



    //-----------------------
    //1: check fin number is correctly unwrapped
    //----------------------
    cout<<"Fin number checking "<<endl;
    unsigned int Nfin_change=0;
    vector<unsigned int> fin_rec;
    vector<unsigned int> fin_transition;fin_transition.clear();
    p=full_path_fin.begin();
    fin_rec.push_back(p->second);
    fin_transition.push_back(0);
    

    for(map<unsigned int, unsigned int>::const_iterator q=full_path_fin.begin();
	q != full_path_fin.end();++q){
      if(q==full_path_fin.begin()) continue;
      p = q;
      p--;
      if(p->second > q->second) ErrorMessage("Fin number is not correctly unwrapped ").throwMe();
      if(p->second != q->second){
	Nfin_change++;
	fin_rec.push_back(q->second);
	fin_transition.push_back(q->first);//keep the last record's id number
      }
    }
    fin_transition.push_back(Nmax);//last record

    cout<<"total number of fin number change "<< Nfin_change<<endl;
    unsigned int Nfin_max=fin_rec.size();
    cout<<"total number of  fins "<< Nfin_max<<endl;
    
    //--------------------------------
    //2: check sab number is increasing with burst
    //--------------------------------
    cout<<"sab number checking "<<endl;
    for(map<unsigned int, unsigned int>::const_iterator q=full_path_sab_counter.begin();    q != full_path_sab_counter.end();++q){
       if(q== full_path_sab_counter.begin()) continue;
       p = q;
       p--;
       if(p->second>= q->second)
	 ErrorMessage("sab number does not increase with burst number ").throwMe();
    }
   

    cout<<"bii number checking"<<endl;
    //--------------------------
    //3: collect same fin number record with bii!=0
    //---------------------------
     for(unsigned int i=0;i<Nfin_max;++i){
       cout<<"precess fin number "<< fin_rec[i]<<endl;
       unsigned int fin_ref = fin_rec[i];
       vector<unsigned int> sab, sclk,bii, ctbc;
       vector<double> bpd, brst;
       vector<unsigned int> tfi;
       vector<unsigned  int> record_id;
       sab.clear(); sclk.clear(); bii.clear(); ctbc.clear();
       bpd.clear();brst.clear();record_id.clear();
       tfi.clear();
   
       for(unsigned int j=fin_transition[i];j<fin_transition[i+1];++j){
	 if(full_path_fin[j]==fin_ref){
	   record_id.push_back(j);
	   sab.push_back(full_path_sab_counter[j]);
	   sclk.push_back(full_path_sclk[j]);
	   bii.push_back(full_path_bii[j]);
	   ctbc.push_back(full_path_ctbc[j]);
	   bpd.push_back(full_path_bpd[j].getInUnits("s"));
	   brst.push_back(full_path_brst[j].getInUnits("s"));
	   tfi.push_back((unsigned int) round_double(full_path_tfi[j].getInUnits("s")));
	 }//collect same fin record
       }//loop over all the records
    
       unsigned int NN=sab.size();
       if(NN != sclk.size() ||
	  NN != bii.size()  ||
	  NN != ctbc.size() ||
	  NN != bpd.size()  ||
	  NN != brst.size() ||
	  NN != record_id.size())
	 ErrorMessage("Utils.cpp: container sizes for same fin number are different").throwMe();
       
       

       if(NN==1){
	 continue;//no fix
       }

       cout<<"sab length "<< sab[0]<<" "<<sab[sab.size()-1]<<endl;       
       //----------------------
       //3.0: tfi correction
       //------------------
       bool found_good_tfi=false;
       unsigned int good_tfi;
       for(unsigned int ii=0;ii<NN-1;++ii){
	 if(tfi[ii]==tfi[ii+1]){
	   good_tfi= tfi[ii];
	   found_good_tfi=true;
	   break;
	 }
       }
       if(!found_good_tfi) ErrorMessage("Utils.cpp:for the same fin, could not find valid tfi").throwMe();
       
       cout<<"good tfi "<< good_tfi<<endl;

       for(unsigned int ii=0;ii<NN;++ii){
	 if(tfi[ii]!=good_tfi){
	   cout<<"tfi error "<< sab[ii]<<endl;
	   cout<<"tfi old "<< tfi[ii]<<endl;
	   cout<<"good tfi "<< good_tfi<<endl;
	   //replace record
	   tfi[ii]=good_tfi;
	   tfi_correction++;
	   //replace original record
	   full_path_tfi[record_id[ii]] = good_tfi; //replace
	 }
       }



       //---------------------
       //3.1  fix possible bii error
       //-------------------------
       bool found_good_bii=false;
       unsigned int good_bii;
       for(unsigned int ii=0;ii<NN-1;++ii){
	 if(bii[ii]==bii[ii+1]){
	   good_bii = bii[ii];
	   found_good_bii=true;
	   break;
	 }
       }
       if(!found_good_bii) ErrorMessage("Utils.cpp: for the same fin, could not find valid bii").throwMe();
       for(unsigned int ii=0;ii<NN;++ii){
	 if(bii[ii]!= good_bii){
	   cout<<"original bii "<<bii[ii]<<endl;
	   cout<<"new bii record "<< good_bii<<endl;
	   //replace container
	   bii[ii]=good_bii;
	   bii_correction++;
	   
	   //replace original bii
	   full_path_bii[record_id[ii]] = good_bii; //replace
	 }
       }
       if(good_bii==0) continue;//fix this kind of stuff later
       


       //--------------------
       // Fix ctbc 
       //ctbc start number: bii -1
       // ctbc start number: 255 when bii=0
       // time lapse from first record: [ctbc(0) - current_ctbc]*bpd 
       //--------------------
       bool fix_ctbc=false;

       //---------------------------------
       //need to fix ctbc using ring buffer
       //--------------------------------
       if(bii[0]!=0 && ctbc[0]==0) {
	 //find the first non-zero
	 unsigned int first_valid_sab=0;
	 unsigned int first_valid_ctbc=0;
	 for(unsigned int ii=1; ii<NN;++ii){
	   if(ctbc[ii]!=0){
	     first_valid_sab=sab[ii];
	     first_valid_ctbc=ctbc[ii];
	     break;
	   }
	 }
	 if(first_valid_sab==0) ErrorMessage("need attention, error is too severe to fix").throwMe();
	 cout<<"fix first ctbc "<< ctbc[0]<< "first valid ctbc "<< first_valid_ctbc<<" "<< "valid sab "<< first_valid_sab<<endl;
	 
	 ctbc[0]= first_valid_ctbc + first_valid_sab -sab[0];
	 cout<<"after fixing first ctbc "<< ctbc[0]<< "first valid ctbc "<< first_valid_ctbc<<" "<< "valid sab "<< first_valid_sab<<endl; 
       }


       if(bii[0]!=0){
	 for(unsigned int ii=1;ii<NN-1;++ii){
	   //reset to 0
	   if(ctbc[ii] < ctbc[ii-1] && ctbc[ii]<= ctbc[ii+1]){
	     cout<<"near by ctbc "<< ctbc[ii-1]<<" "<< ctbc[ii]<<" "<< ctbc[ii-1]<<endl;
	     cout<<"before fix "<< ctbc[ii]<<endl;
	     int sab_diff= sab[ii]- sab[ii-1];
	     if(sab_diff<0) ErrorMessage("sab number is decreasing");
	     ctbc[ii] =  ctbc[ii-1] -(unsigned int) sab_diff;
	     cout<<"after fix "<< ctbc[ii]<<endl;
	   }
	   if( ctbc[ii] > ctbc[ii-1] && ctbc[ii] >= ctbc[ii+1]){
	     cout<<"near by ctbc "<< ctbc[ii-1]<<" "<< ctbc[ii]<<" "<< ctbc[ii-1]<<endl;
	     cout<<"before fix "<< ctbc[ii]<<endl;
	     int sab_diff= sab[ii]- sab[ii-1];
	     if(sab_diff<0) ErrorMessage("sab number is decreasing");
	     ctbc[ii] = ctbc[ii-1] - (unsigned int) sab_diff;
	     cout<<"after fix "<< ctbc[ii]<<endl;
	   }
	 }
       }

       for(unsigned int ii=1;ii<NN;++ii){
	 if((ctbc[ii] > ctbc[ii-1]) && bii[ii]!=0) {
	   fix_ctbc=true;
	   cout<<"before fix ctbc current and past bii sab "<< ctbc[ii]<<" "<< ctbc[ii-1]<<" " << ctbc[ii-2]<<" "<<bii[ii]<<" "<<sab[ii]<< endl;
	   if(ctbc[ii-1]!=0){
	     int sab_diff= sab[ii]- sab[ii-1];
	     if(sab_diff<0) ErrorMessage("sab number is decreasing");
	     ctbc[ii]= ctbc[ii-1]- (unsigned int) sab_diff;//take sab diff
	   }
	   else
	     ctbc[ii]=0;
	   cout<<"after fix [i-1] [i] "<< ctbc[ii-1]<<" "<<ctbc[ii]<<endl;
	 }
       }
       
       if(fix_ctbc==true){
	 for(unsigned int ii=1;ii<NN-1;++ii){
	   if((ctbc[ii] > ctbc[ii-1]) && bii[ii]!=0) {
	     ErrorMessage("could not fix ctbc ").throwMe();
	   }
	 }
       }
 
       //----------------
       //3.2: fix possible bpd error
       //If  any two nearby bursts have the same bpd
       // then take it as good (valid) bpd
       //-----------------
       bool found_good_bpd=false;
       double good_bpd;
       for(unsigned int ii=0;ii<NN-1;++ii){
	 if(bpd[ii]==bpd[ii+1]){
	   good_bpd = bpd[ii];
	   found_good_bpd=true;
	   break;
	 }
       }
       if(!found_good_bpd) ErrorMessage("utils.cpp: for the same fin, any two nearby bpds are different").throwMe();
       cout<<"fin record number "<<i<<endl;
       cout<<"good bpd "<< good_bpd<<endl;

       for(unsigned int ii=0;ii<NN;++ii){
	 if(bpd[ii]!= good_bpd){
	   cout<<"bpd error in sab "<< sab[ii]<<endl;
	   cout<<"bpd correction "<<endl;
	   cout<<"original record "<< bpd[ii]<<endl;
	   cout<<"new record "<< good_bpd<<endl;
	   //replace container
	   bpd[ii]=good_bpd;
	   bpd_correction++;
	
	   //replace original record
	   full_path_bpd[record_id[ii]] = Uvar(good_bpd,"s"); //replace
	 }//if we found any
       }//check all bpd records



      
      

       //-----------------------------------------------------------------------
       //3.3: find first brst of new instruction
       // otherwise first brst set to 0
       //-------------------------------------------
       double first_brst = 0.001;//default 1 ms
       if( (ctbc[0]+1)==bii[0]){//this is the first record
	 first_brst= brst[0];
       }

       if(show_text==1){
	 cout<<"first brst "<< first_brst<<endl;
       }
   

       //----------------
       //3.4: fix time 
       //compare two times: sclk + brst
       //                 : local_trigger_time+tfi+brst[0]+ctbc[i]
       //---------------
       unsigned int corr_sclk= abs(trigger_time_adjust);

       for(unsigned int ii=0;ii<NN;++ii){
	 unsigned int new_sclk= trigger_time + tfi[ii] ;
	 if( (ctbc[ii]+1)>bii[ii]) ErrorMessage("ctbc+1 larger than bii when bii!=0").throwMe();
	 double new_brst= first_brst + bpd[ii]*(bii[ii]-(ctbc[ii]+1));
	 unsigned int add_seconds = (unsigned int) new_brst;
	 new_sclk += add_seconds;
	 new_brst -= add_seconds;

	 if(corr_sclk!=0){
	   for(unsigned int k=0;k<corr_sclk;++k){
	     if(trigger_time_adjust>0) new_sclk++;
	     else new_sclk--;
	   }
	 }
     
	 if(show_text==1){
	   cout<<"sab and bpd "<< sab[ii]<<" "<<bpd[ii]<<endl;
	 }

	 //replace sclk
	 if(sclk[ii]!=new_sclk){ 
	   if(show_text==1)cout<<"old and new sclk         "<< sclk[ii]<<" "<<new_sclk<<endl;
	   sclk_correction++;
	   sclk[ii]= new_sclk;
	   full_path_sclk[record_id[ii]] = new_sclk; //replace
	 }
	 
	 
	 //replace brst
	  if(brst[ii]!=new_brst){
	    brst_correction++;
	    if(show_text==1)  cout<<"old and new brst         "<< brst[ii]<<" "<<new_brst<<endl;
	    brst[ii]= new_brst;
	    full_path_brst[record_id[ii]] = Uvar(new_brst,"s"); //replace
	  }
       }//for all the records with same fin number
     }//loop over different fin numbers with bii!=0


     //--------------------------
     //4: collect same fin number record with bii==0
     //---------------------------
     
     for(unsigned int i=0;i<Nfin_max;++i){
       unsigned int fin_ref = fin_rec[i];
       vector<unsigned int> sab, sclk,bii, ctbc;
       vector<double> bpd, brst;
       vector<unsigned int> tfi;
       vector<unsigned int> record_id;
       sab.clear(); sclk.clear(); bii.clear(); ctbc.clear();
       bpd.clear();brst.clear();record_id.clear();
       tfi.clear();
       
       for(unsigned int j=fin_transition[i];j<fin_transition[i+1];++j){
       	 if(full_path_fin[j]==fin_ref){
	   record_id.push_back(j);
	   sab.push_back(full_path_sab_counter[j]);
	   sclk.push_back(full_path_sclk[j]);
	   bii.push_back(full_path_bii[j]);
	   ctbc.push_back(full_path_ctbc[j]);
	   bpd.push_back(full_path_bpd[j].getInUnits("s"));
	   brst.push_back(full_path_brst[j].getInUnits("s"));
	   tfi.push_back((unsigned int) round_double(full_path_tfi[j].getInUnits("s")));
	 }//collect same fin record
       }//loop over all the records
       
       unsigned int NN=sab.size();
       if(NN != sclk.size() ||
	  NN != bii.size()  ||
	  NN != ctbc.size() ||
	  NN != bpd.size()  ||
	  NN != brst.size() ||
	  NN != record_id.size())
	 ErrorMessage("Utils.cpp: container sizes for same fin number are different").throwMe();
       
       
       
       if(NN==1){
	 continue;
       }
       if(bii[0]!=0) continue;//when bii==0
       
       
       //--------------------------------------------------------------------
       //4.2: find first sclk: compute missing sab between new instruction
       // and first sab record
       //-------------------------------------------
       //unsigned int record_first_sclk= sclk[0];
       //unsigned int theory_first_sclk= trigger_time + tfi[0];
       
       double first_brst=0.0;
       unsigned int first_sclk= sclk[0];
       
       
       //----------------
       //4.3: fix time 
       // since bii=0, ctbc becomes useless
       // use sab number difference
       // use first record as reference
       //---------------
       
       for(unsigned int ii=0;ii<NN;++ii){
	 unsigned int new_sclk=  first_sclk ;
	 double new_brst= first_brst+ bpd[ii]*(sab[ii]-sab[0]);
	 unsigned int add_seconds = (unsigned int) new_brst;
	 new_sclk += add_seconds;
	 new_brst -= add_seconds;
     

	 if(show_text==1)
	   cout<<"sab and bpd "<< sab[ii]<<" "<<bpd[ii]<<endl;
	

	 //replace sclk
	 if(sclk[ii]!=new_sclk){
	   if(show_text==1)  cout<<"old and new sclk         "<< sclk[ii]<<" "<<new_sclk<<endl;
	   sclk_correction++;
	   sclk[ii]= new_sclk;
	   full_path_sclk[record_id[ii]] = new_sclk; //replace
	 }


	 //replace brst
	 if(brst[ii]!=new_brst){
	   if(show_text)   cout<<"old and new brst         "<< brst[ii]<<" "<<new_brst<<endl;
	   brst_correction++;
	   brst[ii]= new_brst;
	   full_path_brst[record_id[ii]] = Uvar(new_brst,"s"); //replace
	 }
       }//for all the records with same fin number
     }//loop over different fin numbers
  }
     
     
  

//-------------------------------
//transformation matrix using doubles
//--------------------------------

void beam_xyz(const double& azi, const double& elev,
	      double& x, double& y, double& z){
  x=cos(elev)*sin(azi);
  y=sin(elev);
  z=sqrt(1.0 - x* x -y*y);
}

void beam_azim_elev(const double& x, const double& y, const double& z,
		    double& azim, double& elev)
{
  
  elev = asin(y/sqrt(x*x + y*y + z*z));
  azim= asin( x/cos(elev)/sqrt(x*x+y*y+z*z));
}


void compute_lat_lon(const double& px,
		     const double& py,
		     const double& pz,
		     double& lat,
		     double& lon){
  
  lon = atan2(py,px);
  lat = -acos(pz/sqrt(px*px +py*py + pz*pz))+pi/2.0;
}


void compute_TRMatrix_from_FSC_to_Ftarget(const double& sc_X_x,
					  const double& sc_X_y,
					  const double& sc_X_z,
					  const double& sc_Y_x,
					  const double& sc_Y_y,
					  const double& sc_Y_z,
					  const double& sc_Z_x,
					  const double& sc_Z_y,
					  const double& sc_Z_z,
					  const int& sc_x_to_vx,
					  Dmat& from_FSC_to_Ftarget)
{
  if(!(sc_x_to_vx==1 || sc_x_to_vx==-1)) ErrorMessage("sc x direction should be either 1 to -1").throwMe();
  if(!(from_FSC_to_Ftarget.cols()==3 || from_FSC_to_Ftarget.rows()==3)) ErrorMessage("matrix size is not 3x3").throwMe();
  from_FSC_to_Ftarget(0,0)=sc_X_x * sc_x_to_vx;
  from_FSC_to_Ftarget(0,1)=sc_Y_x * sc_x_to_vx;
  from_FSC_to_Ftarget(0,2)=sc_Z_x;

  from_FSC_to_Ftarget(1,0)=sc_X_y * sc_x_to_vx;
  from_FSC_to_Ftarget(1,1)=sc_Y_y * sc_x_to_vx;
  from_FSC_to_Ftarget(1,2)=sc_Z_y;

  from_FSC_to_Ftarget(2,0)=sc_X_z * sc_x_to_vx;;
  from_FSC_to_Ftarget(2,1)=sc_Y_z * sc_x_to_vx;;
  from_FSC_to_Ftarget(2,2)=sc_Z_z;


}


void compute_TRMatrix_product(const Dmat& from_B_to_A,
			      const Dmat& from_C_to_B,
			      Dmat& from_C_to_A)
{
 if(!(from_B_to_A.cols()==3 || from_B_to_A.rows()==3)) ErrorMessage("matrix size is not 3x3").throwMe();
 if(!(from_C_to_B.cols()==3 || from_C_to_B.rows()==3)) ErrorMessage("matrix size is not 3x3").throwMe();
 if(!(from_C_to_A.cols()==3 || from_C_to_A.rows()==3)) ErrorMessage("matrix size is not 3x3").throwMe();
 

 
 for(unsigned int ii=0;ii<3;++ii)
   for(unsigned int jj=0;jj<3;++jj){
     from_C_to_A(ii,jj)=0.0;
     for(unsigned int kk=0;kk<3;++kk)
       from_C_to_A(ii,jj) += from_B_to_A(ii,kk) * from_C_to_B(kk,jj);
   }

}



void compute_New_Vector(const Dmat& from_B_to_A,
			const double& x_in_B,
			const double& y_in_B,
			const double& z_in_B,
			double& x_in_A,
			double& y_in_A,
			double& z_in_A)
{
   if(!(from_B_to_A.cols()==3 || from_B_to_A.rows()==3)) ErrorMessage("matrix size is not 3x3").throwMe();

   x_in_A = from_B_to_A(0,0)* x_in_B + from_B_to_A(0,1)* y_in_B + from_B_to_A(0,2) * z_in_B;
   y_in_A = from_B_to_A(1,0)* x_in_B + from_B_to_A(1,1)* y_in_B + from_B_to_A(1,2) * z_in_B;
   z_in_A = from_B_to_A(2,0)* x_in_B + from_B_to_A(2,1)* y_in_B + from_B_to_A(2,2) * z_in_B;

}

void compute_range_doppler(const double& px,
			   const double& py,
			   const double& pz,
			   const double& vx,
			   const double& vy,
			   const double& vz,
			   const double& wavelength_in_m,
			   const double& look_x,
			   const double& look_y,
			   const double& look_z,
			   const double& radius,
			   double& range,
			   double& doppler,
			   bool& found_surface)
{
  
  doppler = 2.0/wavelength_in_m *1000.0 *(vx * look_x + vy*look_y + vz*look_z);

  double distance = sqrt(px*px + py*py +pz*pz);
  double dir_x = px/ distance;
  double dir_y = py/ distance;
  double dir_z = pz/distance;
  
  double dot_pos_look = dir_x * look_x + dir_y * look_y + dir_z * look_z;
  


  double root=radius*radius - distance*distance *(1-dot_pos_look*dot_pos_look);
  if(root<0.0) {
    range=0.0;
    found_surface = false;
    return;
  }
  range = -distance* dot_pos_look - sqrt(root);
  found_surface=true;
}


//------------
//fitting used by dlap
//------------
void dlap_min_fit(const Uvar& start_tracking,
	     const Uvar& lat_start,
	     const Uvar& dlat1_dt,
	     const Uvar& end_tracking,
	     const Uvar& lat_end ,
	     const Uvar& dlat2_dt,
	     const Uvar& target_lat,
	     const Uvar& target_CA_time,
	     Dvec& poly_coeff ){
  double t1= (start_tracking - target_CA_time).getInUnits("min");
  double lat1 = (lat_start - target_lat).getInUnits("deg");
  double lat1_dt= dlat1_dt.getInUnits("deg/min");
  double t2 = (end_tracking-target_CA_time).getInUnits("min");
  double lat2=(lat_end - target_lat).getInUnits("deg");
  double lat2_dt=dlat2_dt.getInUnits("deg/min");
  double_dlap_min_fit(t1,lat1,lat1_dt,
		      t2,lat2,lat2_dt,
		      poly_coeff);
}

void double_dlap_min_fit(const double& t1, 
		    const double& lat1, 
		    const double& dlat1_dt,
		    const double& t2,
		    const double& lat2,
		    const double&  dlat2_dt,
		    Dvec& poly_coeff )
{
  Dvec B("",4);
  Dmat A("",4,4);
  A(0,0)= pow(t1,2); A(0,1)=pow(t1,3); A(0,2)=pow(t1,4); A(0,3)=pow(t1,5); B(0)=lat1;
  A(1,0)= pow(t2,2); A(1,1)=pow(t2,3); A(1,2)=pow(t2,4); A(1,3)=pow(t2,5); B(1)=lat2;
  A(2,0)= 2.0*t1; A(2,1)= 3.0*pow(t1,2); A(2,2)= 4.0*pow(t1,3); A(2,3)= 5.0*pow(t1,4); B(2)=dlat1_dt;
  A(3,0)= 2.0*t2; A(3,1)= 3.0*pow(t2,2); A(3,2)= 4.0*pow(t2,3); A(3,3)= 5.0*pow(t2,4); B(3)=dlat2_dt;   
  A.inv();
  for(unsigned int c=0;c<4;++c){
    poly_coeff(c)=0.0;
    for(unsigned int d=0;d<4;++d)
      poly_coeff(c) += A(c,d)*B(d);
  }
}

void find(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list){

  index_list.clear();
  if(target_list.size()==0) return;

  for(unsigned int i=0;i<target_list.size();++i){
    if(target_list[i]== value){
      index_list.push_back(i);
    }
  }
}

void find_less_than_target_value(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list){
  index_list.clear();
  if(target_list.size()==0) return;
  
  for(unsigned int i=0;i<target_list.size();++i){
    if(target_list[i] <=  value){
      index_list.push_back(i);
    }
  }

}
void find_greater_than_target_value(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list){
index_list.clear();
  if(target_list.size()==0) return;
  
  for(unsigned int i=0;i<target_list.size();++i){
    if(target_list[i]>= value){
      index_list.push_back(i);
    }
  }
}


//-------------------
//Sphere Geometry: cppied from E.R.'s sphereical geometry computation method
//--------------------

double r_to_theta(double height, double radius, double range){
  
  double cost, scale , r_over_R;

  scale = 1. + height/radius;
  r_over_R = range/radius;

  cost = (scale*scale + r_over_R*r_over_R - 1)/(2.*r_over_R*scale);

  return acos(cost);

}

double r_to_theta(double height, double radius, double range, double dem_h){
  
  double cost, scale , r_over_R;
  double gamma;

  scale = 1. + height/radius;
  r_over_R = range/radius;
  gamma= dem_h/radius;
  cost = (scale*scale + r_over_R*r_over_R - (1.0+gamma)*(1.0+gamma))/(2.*r_over_R*scale);

  return acos(cost);

}



double theta_to_r(double height, double radius, double theta){
  double cost, scale;/* , r_over_R;*/
    cost = cos(theta);
    scale = 1. + height/radius;
    scale *= scale;

    return (height+radius)*cost*(1. - sqrt(1. - (scale - 1.)/
                                           (scale*cost*cost)));
}

double r_to_alpha(double height, double radius, double range){
  return asin(range*sin(r_to_theta(height,radius,range))/radius);
}

double  alpha_to_r(double height, double radius, double alpha){
   double cosa, scale;
   scale = 1. + height/radius;
   cosa = cos(alpha);

   return radius*sqrt(1. + scale*scale - 2.*scale*cosa);
}

double r_to_rho(double height, double radius, double range){
  return radius*r_to_alpha(height,radius,range);
}

double rho_to_r(double height, double radius, double rho){
   return alpha_to_r(height,radius,rho/radius);
}

double r_to_inc(double height, double radius, double range){
  return r_to_theta(height,radius,range) + r_to_alpha(height,radius,range);
}

double inc_to_r(double height, double radius, double inc){
  double cosi, scale;
  cosi = cos(inc);
  scale = 1. + height/radius;
  scale *= scale;
  return radius*cosi*(sqrt(1. + (scale - 1.)/(cosi*cosi)) - 1.);
}

double theta_to_alpha(double height, double radius, double theta){
  return r_to_alpha(height,radius,theta_to_r(height,radius,theta));
}

double alpha_to_theta(double height, double radius, double alpha){
  return r_to_theta(height,radius,alpha_to_r(height,radius,alpha));
}

double theta_to_rho(double height, double radius, double theta){
  return r_to_rho(height,radius,theta_to_r(height,radius,theta));
}

double rho_to_theta(double height, double radius, double rho){
   return r_to_theta(height, radius,rho_to_r(height, radius, rho));
}

double theta_to_inc(double height, double radius, double theta){
  return r_to_inc(height,radius,theta_to_r(height,radius,theta));
}

double inc_to_theta(double height, double radius, double inc){
  return r_to_theta(height,radius,inc_to_r(height,radius,inc));
}

double rho_to_inc(double height, double radius, double rho){
  return r_to_inc(height,radius,rho_to_r(height,radius,rho));
}


