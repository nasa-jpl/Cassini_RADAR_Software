//=============================================================================
// Utils.h
//
// This file contains utility function declarations and macros.
// 
// Macros
// 
// MAX(A,B)
// MIN(A,B)
//
// Interface summary:
//
// str = toStr(ostrstream);     Convert ostrstream.str() to std::string
// str = toStr(arg);            Convert standard type for string use
//
// getlookDirection(DirectionVector& look,
//		      unsigned int& solution,
//		      const PositionVector& sc_pos,
//		      const FloatVector& sc_vel,
//                  const Uvar& doppler,
//                  const Uvar& radius) 
//                 throw (Unit::UnitError,ErroMessage)
//
//    look : DirectionVector corresponding to a particular look direction
//    sc_pos_In_fbeam: DirectionVector of sc_pos in fbeam frame
//    sc_vel_In_fbeam: DirectionVector of sc_vel in fbeam frame
//    range : range
//    doppler: doppler
//  For given range and doppler, it will determine the look direction vector 
//  based on sc_position and velocity direction vector represented in beam 
//  frame
//
//
//
// polyfit(const ARRAY_DOUBLE1& poly_coeff, 
//        const ARRAY_DOUBLE1& angle,
//       const ARRAY_DOUBLE1& altitude)
//   throw (Unit::UnitError,ErrorMessage)
//
//
//  void bitset(unsigned int& target, 
//              const unsigned int& b1, 
//              const unsigned int& b2);
//         bitset manipulation utility for 16 bits unsigned int type
//
//
//  void bitdisplay(unsigned int& target)
//  display bits on the screen
// double muhleman_backscatter(const double& k1, const double& k2, const Uvar& incidenceangle);
//=============================================================================

#ifndef Utils_H
#define Utils_H

static const char rcs_id_utils_h[] =
  "@(#) $Id: Utils.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <stdio.h>
#include <string>
#include <sstream>

using std::string;
using std::ostringstream;

// Macros
#define MAX(A,B)      ((A)>(B)?(A):(B))
#define MIN(A,B)      ((A)<(B)?(A):(B))

//---------------------------------------------------------------------------
// str = toStr(arg)
//
// Functions to convert various types to strings.
//---------------------------------------------------------------------------

string toStr(const unsigned int& arg);
string toStr(const int& arg);
string toStr(const unsigned long& arg);
string toStr(const long& arg);
string toStr(const double& arg);
string toStr(const float& arg);
string toStr(ostringstream& arg);
string toStr(string& arg);

//---------------------------------------------------------------------------
// str = toStr(arg,min_width)
//
// Convert the integer argument to a string with specified minimum width
// and zero fill when needed.
//---------------------------------------------------------------------------

string toStr(unsigned int arg, unsigned int min_width);

// These inludes need to come after the toStr declarations to avoid a
// dependency on order for other include files.  Note that toStr is often
// used in template definitions in include files.
#include "Error.h"
#include "Units.h"
#include "Frame.h"
#include "Array.h"
#include "Time.h"
#include "Config.h"
#include "Ieb.h"

string get_next_token(string& s, const string& w);

void spice_setup(Config& cfg);

void getlookDirection(DirectionVector& look,
		      unsigned int& solution,
		      const PositionVector& sc_pos,
		      const FloatVector& sc_vel,
		      const Uvar& range, 
		      const Uvar& doppler,
		      const Uvar& lambda,
		      const Uvar& radius) 
     throw (Unit::UnitError,ErrorMessage);

float number_of_looks(const StateVector& sc_state,
		      const Uvar& range,
		      const Uvar& doppler,
		      const Uvar& bpd,
		      const Uvar& lambda,
		      const Uvar& pri,
		      const Uvar& process_bandwidth,
		      const unsigned int& number_of_pulses );

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
			  Array2D<unsigned int>& no_of_solution_grid);

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
			   Imat& valid_pixel);

void compute_min_max_range(const Time& t,
			   const StateVector& sc_state,
			   const string& target,
			   const Frame& fbeam,
			   const Uvec& azi_angles,
			   const Uvec& elev_angles,
			   Uvar& min_range,
			   Uvar& max_range);

void compute_bore_range_doppler(const Time& t,
				const StateVector& sc_state,
				const string& target,
				const Frame& fbeam,
				const Uvar& lambda,
				Uvar& bore_range,
				Uvar& bore_doppler);

// Nice fopen which returns FILE* and throws errors
FILE* fopen(const string& name, const string& mode);

// Check for file existence using stat(); throws errors
bool fileExists(const string& filename);

//-------------------
// find min and max
//---------------------
void min_max(Uvar& min, 
	     Uvar& max,
	     const Uvec& data )
     throw (Unit::UnitError,ErrorMessage);

void min_max(Uvar& min, 
	     Uvar& max,
	     const Umat& data )
     throw (Unit::UnitError,ErrorMessage);

void min_max(Uvar& min, 
	     Uvar& max,
	     const vector<Uvar>& data );
    
void min_max( double& min, double& max, const Dvec& value);
//--------------------------------
//polyfit
//---------------------------------
void  polyfit(Dvec& poly_coeff, 
	      const Dvec&  y,
	      const Dvec&  x)
     throw (ErrorMessage);
//--------------------------------
//linear fit
//---------------------------------
void linearfit(const vector<Uvar>& x, const vector<Uvar>& y,
	       const vector<Uvar>& sig, const int& mwt,
	       Uvar& a, Uvar& b,
	       Uvar& sig_a, Uvar& sig_b,
	       double& chi2, double& q);

 void linearfit(const Uvec& x, const Uvec& y,
		const Uvec& sig, const int& mwt,
		Uvar& a, Uvar& b,
		double& sig_a, double& sig_b,
		double& chi2, double& q);

 void linearfit(const Dvec& x, const Dvec& y,
		const Dvec& sig, const int& mwt,
		double& a, double& b,
		double& sig_a, double& sig_b,
		double& chi2, double& q);
//---------------------
//bit set and display
//-------------------
 void bitset(unsigned short& target, 
	     const unsigned int& b1, 
	     const unsigned int& b2,
	     const unsigned int& setvalue);

void bitset(unsigned int& target,
	       const unsigned int& bit,
	       const unsigned int& setvalue);

void bitdisplay(const unsigned short& target);

//-----------------------
//backscattering based on muhleman model
//------------------------
double muhleman_backscatter(const double& k1, 
			    const double& k2, 
			    const Uvar& incidenceangle);
double muhleman_backscatter(const double& k1, 
			    const double& k2, 
			    const double& incidenceangle_in_rad);
double Lambertian(const double& x,
		  const Uvar& angle);
//------------------------------
//surface area computataion
//------------------------------
Uvar surfacearea(const PositionVector& p1, 
		 const PositionVector& p2, 
		 const PositionVector& p3, 
		 const Uvar& sphere_radius);
//calculate  area bounded by three position vector on the surface of a sphere




//--------------------------
//round double to integer
//---------------------------
int round_double(const double& x);

//------------------------------------------
// Compute a two-dimensional sync function
// sinc
//------------------------------------------

double sinc(double d);

//-------------------------------------------
// Bilinear Interpolation
//-------------------------------------------

double bilinear(double a, double b, const Dmat& d) throw(ErrorMessage);
Uvar bilinear(double a, double b, const Umat& d) throw(ErrorMessage);
void bilinear(const Uvec& x1, const Uvec& x2, const Umat& y,
  const Uvar& out_of_bounds_value,
  const Uvec& newx1, const Uvec& newx2, Umat& newy);

//--------------------------
// Set unit handling option
//--------------------------
void set_unit_option(char* unit_option);


/**  Problems with this template definition on new compiler
//---------------------------------------------------------------------------
// str = toStr(arg)
//
// Convert the argument to a string using standard stream i/o.
//---------------------------------------------------------------------------

template<class T>
string toStr(const T& arg)
  {    
  ostringstream ost;
  ost << arg;
  // Following conditional compilation is needed to allow gcc-2.95 compilers
  // to successfully replace stringstream with strstream.
#ifdef __STRSTREAM__
  return(string(ost.str(),ost.pcount()));
#else
  return(ost.str());
#endif
  } 

//---------------------------------------------------------------------------
// str = toStr(ostringstream)
//
// Specialization of toStr for a stringstream input.
//---------------------------------------------------------------------------

template<> string toStr<ostringstream>(const ostringstream& arg)
  {
  return(arg.str());
  }

**/

//----------------------------------------------
//Decode Cassini Encoded Data
//----------------------------------------------
void decodeCassiniData(const vector<unsigned char>& encoded_data,
		       const Array1D<unsigned int>& baq_threshold,
		       Ieb& ieb,
		       Array1D<float>& decoded_data,
		       unsigned int& Ndecoded,
		       float& offset);
//-------------------------------------
//functions to be used for direct best fit ellipse
//-------------------------------------
void EllipseFit(const Dvec& x_input, 
		const Dvec& y_input, 
		Dvec& ellipse_x_out,
		Dvec& ellipse_y_out);

void jacobi(double **a, int n, double *d , double **v, int nrot)    ;
void choldc(double **a, int n, double **l);
int inverse(double **TB, double **InvB, int N) ;
void A_TperB(double **A, double **B, double **res,
	     int righA, int colA, int righB, int colB) ;
void AperB(double **A, double **B, double **res, 
	   int righA, int colA, int righB, int colB) ;

void AperB_T(double **A, double **B, double **res,
	     int righA, int colA, int righB, int colB) ; 
void ROTATE(double **a, int i, int j, int k, int l, double tau, double s) ;




//---------------------------------------
//compute target frame's rotation angles
// and rotation rate in J2000
//---------------------------------------

void getTBC_Frame_RA_DEC_ROTATION_in_J2000(const Time& t,
						const string& target_name,
						Uvar& ra, 
						Uvar& dec, 
						Uvar& rotation_angle);
Uvar getTBC_Frame_Rotation_Rate_in_J2000(const Time& t,
					 const string& target_name,
					 const Uvar& interval);

//-------------------------
//compute spacecraft yaw, pitch and roll
//--------------------------
void computeYaw_Pitch_Roll(const StateVector& sc_state,
			   const Frame& ftarget,
			   const Frame& fsc,
			   Uvar& yaw,
			   Uvar& pitch, 
			   Uvar& roll);
void  fit_pitch_polynomial( const string& before_or_after_epoch,
			    const Uvar& start_time, 
			    const Uvar& end_time,
			    const Uvar& pitch_rate,
			    Dvec& poly_fit,
			    Uvar& mid_time);

//------------------------------------------------------------------
// Utility functions for doppler centroid
// and related tools including  singular value decomposition(SVD) method
//--------------------------------------------------------------
//compute doppler and df/dr, df/dp, df/dy
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
		     double& r_dfdvh,
		     double& r_dfdhref,
		     double& r_dfdwvl,
		     double& r_dfdasa);

//doppler derivative using double array and variables
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
		      double& r_dfdasa);

//method of singular value decomposition
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
	       const Array1D<int>& i_usedata);
void  funcs(const unsigned int& i_q,
	    const unsigned int& i_rd,
	    const unsigned int& i_fp,
	    const Dvec& r_vecin,
	    const unsigned int& i_ma,
	    const Dvec& r_a,
	    Dmat& r_amat,
	    const Array1D<int>& i_paramest,
	    const Array1D<int>& i_usedata);
void svdvar(Dmat& r_v, Dvec& r_w, Dmat& r_cvm);
void gaussj(Dmat& a, 	    Dmat& b);
void  svdcmp(Dmat& a,	     Dvec& w,	     Dmat& v);
void svdcmp_c(double **a, int m, int n, double w[], double **v);
void svdcmp_cpp(Dmat& a, Dvec& w, Dmat& v);
double pythag(const double& a,
	      const double& b);
void  svbksb(Dmat& u,
	     Dvec& w,
	     Dmat& v,
	     Dvec& b,
	     Dvec& x);
double SIGN_C(const double& a, const double& b);




//------------------------------------
//Fix Cassini Sab counter and S/C SCLK counter
//------------------------------------
void fixCassiniSabSclk(map<unsigned int, unsigned int>& full_path_sab_counter,
		       map<unsigned int, unsigned int>& full_path_sclk,
		       unsigned int& first_good_record,
		       unsigned int& num_sab_counter_roll_over,
		       const bool& plot_option);

void fixFinNumberRollOver(map<unsigned int, unsigned int>& full_path_fin,
			  const unsigned int& number);

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
			 unsigned int& bii_corection);



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
			       const unsigned int& show_text);




//------------------------------------
//take care of transformation matrix
//----------------------------------
void beam_xyz(const double& azi, const double& elev,
	      double& x, double& y, double& z);

void beam_azim_elev(const double& x, const double& y, const double& z,
		    double& azim, double& elev);

void compute_lat_lon(const double& px,
		     const double& py,
		     const double& pz,
		     double& lat,
		     double& lon);

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
					  Dmat& from_FSC_to_Ftarget);

void compute_TRMatrix_product(const Dmat& from_B_to_A,
			      const Dmat& from_C_to_B,
			      Dmat& from_C_to_A);

void compute_New_Vector(const Dmat& from_B_to_A,
			const double& x_in_B,
			const double& y_in_B,
			const double& z_in_B,
			double& x_in_A,
			double& y_in_A,
			double& z_in_A);

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
		   bool& found_surface);



//----------------
//fit 5th poly
//  (t1, x1(t1), dx1/dt1,
//  t2, x2(t2), dx2/dt2,
//  t0, coeff)
// dx/dt0=0
//-----------------
void dlap_min_fit(const Uvar& start_tracking,
		  const Uvar& lat_start,
	     const Uvar& dlat1_dt,
	     const Uvar& end_tracking,
	     const Uvar& lat_end ,
	     const Uvar& dlat2_dt,
	     const Uvar& target_lat,
	     const Uvar& target_CA_time,
	     Dvec& poly_coeff );

void double_dlap_min_fit(const double& t1, 
		    const double& lat1, 
		    const double& dlat1_dt,
		    const double& t2,
		    const double& lat2,
		    const double&  dlat2_dt,
		    Dvec& poly_coeff );


void find(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list);
void find_less_than_target_value(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list);
void find_greater_than_target_value(const vector<Uvar>& target_list, const Uvar& value, vector<unsigned int>& index_list);



/***********************************************************************
Given the height and the planet radius and one of the following
parameters:

range: range from platform to surface
theta: look angle relative to the -vertical direction at the platform
alpha: interior angle from the center of the planet to the platform and
       scattering point
rho:   distance along the surface from the platform nadir point to the
       scattering point
inc:   the incidence angle measured relative to the local normal at the
       scattering point

the following set of functions allows for the conversion to any of the
other parameters in the list.

PROGRAMMER: E. Rodriguez
Copied from YGim's GLISTIN package June 18, 2007
***********************************************************************/

double r_to_theta(double height, double radius, double range);

double r_to_theta(double height, double radius, double range, double dem_h);

double theta_to_r(double height, double radius, double theta);

double r_to_alpha(double height, double radius, double range);

double alpha_to_r(double height, double radius, double alpha);

double r_to_rho(double height, double radius, double range);

double rho_to_r(double height, double radius, double rho);

double r_to_inc(double height, double radius, double range);

double inc_to_r(double height, double radius, double inc);

double theta_to_alpha(double height, double radius, double theta);

double alpha_to_theta(double height, double radius, double alpha);

double theta_to_rho(double height, double radius, double theta);

double rho_to_theta(double height, double radius, double rho);

double theta_to_inc(double height, double radius, double theta);

double inc_to_theta(double height, double radius, double inc);

double rho_to_inc(double height, double radius, double rho);



//SH's method of heading, reast,rnorth
void geo_hdg(const double& r_a, const double& r_e2,const double& r_lati, const double& r_loni,
	     const double& r_latf, const double& r_lonf, double&  r_geohdg);

void geo_dis(const double& r_a, const double& r_e2,const double& r_lati, const double& r_loni,
	     const double& r_latf, const double& r_lonf,double& r_geodis);
double  reast(const double& r_a, const double& r_e2,const double& r_lat);
double rnorth(const double& r_a, const double& r_e2,const double& r_lat);
double radius_along(const double& r_a, const double& r_e2,const double& r_hdg, const double& r_lat);
double radius_cross(const double& r_a, const double& r_e2,const double& r_hdg, const double& r_lat);


#endif 




