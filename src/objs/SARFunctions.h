#ifndef SARFUNCTIONS_H
#define SARFUNCTIONS_H

#include <string>
#include "Units.h"
#include "Error.h"
#include "Array.h"
#include "Ieb.h"


using std::complex;





// correlates x[i ... i+N-1] with f[0...N] and returns in d[0 ... N]
void range_compress(const fdata& x, unsigned int i, unsigned int N, 
		    const CFvec& f, CFvec& d);
void dig_chirp_hensley(double r_chirplen, double r_bw, double r_stepsize, 
		       double r_f0, double r_phase0, double r_sampfreq, 
		       int& i_ts, complex<float>* c_chirp);
float
get_point_target_echo(double tp, double t, 
		      double updown_freq_shift, 
		      double carrier_freq,  
		      double delay, double chirp_start_freq,
		      double chirp_rate, double  mag, double phi0 );

void dft(const CFvec& s, CFvec& S, unsigned int N);

float linear_fm_chirp(const Uvar& t, const Uvar& f0, const Uvar& cr, 
		      float mag, const Uvar& phi0);
float anti_aliased_linear_fm_chirp(const Uvar& t, const Uvar& f0, 
				   const Uvar& cr, float mag, 
				   const Uvar& phi0, 
				   const Uvar& pass_band_start,
				   const Uvar& pass_band_end);
complex<float> 
      complex_anti_aliased_linear_fm_chirp(const Uvar& t, const Uvar& f0, 
					   const Uvar& cr, float mag, 
					   const Uvar& phi0, 
					   const Uvar& pass_band_start,
					   const Uvar& pass_band_end);
void fft(const CFvec& f, CFvec& F);
void ifft(const CFvec& F, CFvec& f);
void fft_base(float* DATA, int NN, int ISIGN);
void fft(complex<float>*f, complex<float>*F,int N);
void ifft(complex<float>*F, complex<float>*f,int N);
void dieIfNotPowerofTwo(int N);


//--------------------------------------------
//Construction of chirp based on ieb paramters
//-------------------------------------------
float linear_chirp(const double& slow_time, 
		   const double& t, 
		   const double& fast_time,
		   const double& chirp_start_frequency_in_Hz, 
		   const double& chirp_rate,
		   const double& cfs_in_Hz,
		   const double& phase);
		   
complex<float> complex_linear_chirp(const double& slow_time, 
		    const double& t, 
		    const double& fast_time,
		    const double& chirp_start_frequency_in_Hz, 
		    const double& chirp_rate,
		    const double& cfs_in_Hz);

float digital_chirp(const double& slow_time,
		     const double& t, 
		     const double& fast_time,
		     const double& chirp_start_frequency_in_Hz, 
		     const int& nstep, 
		     const double& csd_in_s, 
		     const double& cfs_in_Hz);

complex<float> complex_digital_chirp(const double& slow_time,
				    const double& t, 
				    const double& fast_time,
				    const double& chirp_start_frequency_in_Hz, 
				    const int& nstep, 
				    const double& csd_in_s, 
				    const double& cfs_in_Hz);

//----------------------------------
//Single target echo: forward method:
//  10MHz sampling and adc subsampling
//-----------------------------------
void generate_rerouted_chirp(Ieb& ieb,
			     const double& phase,
			     fdata& echo);
void generate_echo_using_LFM(Ieb& ieb, 
			     const string& target,
			     const Frame& target_frame,
			     const Uvar& lat,
			     const Uvar& lon,
			     const Time& burst_time,
			     const double& phase,
			     fdata& echo);
			    //forward solution

void generate_echo_using_LFM_BIF(Ieb& ieb, 
				 const string& target,
				 const Frame& target_frame,
				 const Uvar& lat,
				 const Uvar& lon,
				 const Time& burst_time,
				 const double& phase,
				 fdata& echo,
				 const unsigned int& num_bursts_in_flight);


//---------------------------------
//Match filter based on Ieb
//---------------------------------
void compute_matched_filter(const double& adc,
			    const double& taup,
			    const double& chirp_rate,
			    const double& cfs,
			    const double& filter_start,
			    const double& filter_end,
			    const double& start_frequency,
			    const unsigned int& N,
			    unsigned int& i_npts,
			    CFvec& matched_filter);
void compute_matched_filter(Ieb& ieb,
			    const Uvar& doppler,
			    const Uvar& filter_start,
			    const Uvar& filter_end,
			    const unsigned int& N,
			    unsigned  int& i_npts,
			    CFvec&  matched_filter);

void compute_matched_filter_dig( double r_chirplen,
				 double r_bw,
				 double r_stepsize,
				 double r_phase0,
				 double r_sampfreq,
				 double r_f0,
				 int& i_ts,
				 const int& N,
				 CFvec& matched_filter);






// fast/simple vetor and matrix manipulation routines
double dot(double v1[3], double v2[3]);
double vector_dist(double x[3],double y[3]);
void rotate_vector(double vfrom[3], double rotmat[3][3], double vto[3]);

// optimized geometry routines used by SAR processor

int
fast_doppler_range_to_TBF_look(double doppler_in_hz, double range_in_km,
			       double target_rad,
			       const double p[3], const double v[3],
			       double lambda, double u1[3], double u2[3]);

int
get_surface_intercept_triaxial(double scpos[3], double look[3], double radii[3], double pos[3], double& r);

double getSphericalTriangleArea(double p1[3], double p2[3],
				double p3[3],double r,
				double& s12, double& s23,
				double& s31);

double getFlatTriangleArea(double p1[3], double p2[3],
				double p3[3]);

Uvar getSphericalTriangleArea(const PositionVector& p1, 
			      const PositionVector& p2, 
			      const PositionVector& p3, 
			      const Uvar r,
			      Uvar& s12, Uvar& s23,
			      Uvar& s31);

Uvar getLatLonPixelAreaOnSphere(const Uvar& lat, const Uvar& lon,
				const Uvar& latres, const Uvar& lonres);




//--------------------
//BAQ Encoding
//-------------------
void BAQEncode(Ieb&ieb, const Uvec& echo, SUvec& thresh, CharUvec& words);


//--------------------------------
// functions for Doppler centroid: pulse segmentation for doppler centroid
// computation
//-----------------------------
void ShiftEcho(Fvec& echo_data, const Uvar& dt, const Uvar& SR);
void Convert2Pow(unsigned int& n);
void Convert2DEcho(const Fvec& echo_data, Fmat& echo2D, const unsigned int& N_pulse, const unsigned int& N_in_pri);
void Convert2DEcho(const Fvec& echo_data, 
		   const unsigned int& start_index,
		   const unsigned int& end_index, 
		   const unsigned int& N_samples_per_window,
		   const unsigned int& N_pulses_inside_echo,
		   const unsigned int& N_samples_per_pri,
		   Fmat& echo2D);
void Convert2DReal2Complx(const Fmat& echo2D, CFmat& c_echo2D,
			  const int& band_selection);
void ConvertRealToComplex(const Fvec& echo, CFvec& c_echo,
			  const int& band_selection);
//void ComputPulseCorrelation(const CFmat& c_data2D, CFmat& c_corr2D);
//void AveragePulseCorrelation(const CFmat& c_corr2D, CFvec& c_corr1D, complex<float>& c_corr0D);
void ComputFractDopCentroid(const CFmat& echo2D_c, Dvec& dop_frac1D, double& dop_frac0D);
void ComputIntDopCentroid(const Uvar& bore_dop, const Uvar& pri, int& dop_int);
void CompareDopCentroid(const Uvar& exact_dop, const Fvec& dop_frac1D, const int& dop_int, const Uvar& pri);
void PRFHopping(const Uvec& PRFs, const Uvec& est_dop, const Uvec& fdop, Uvec& hop_dop);
void Plot1DReal(const Fvec& data1D);
void Plot2DReal(const Fmat& data2D);
void Plot1DComplx(const CFvec& c_data1D);
void Plot2DComplx(const CFmat& c_data2D);
unsigned int get_next_power_of_2(unsigned int n);
void simpleRangeDopplerCompression(const CFmat& echo2D_c,
				   const CFvec& fft_matched_filter,
				   const Uvar& prf,
				   const Uvar& bore_doppler,
				   const Uvar& adc,
				   const Uvar& range_start,
				   Uvec& range_axis,			   
				   Uvec& doppler_axis,
				   CFmat& echo2D_rc,
				   CFmat& echo2d_rdc);
				 

void simpleRangeDopplerCompression2(const CFmat& echo2D_c,
				    const CFvec& fft_matched_filter,
				    const Uvar& prf,
				    const Uvar& bore_doppler,
				    const Uvar& fdop_rate,
				    const Uvar& adc,
				    const Uvar& range_start,
				    Uvec& range_axis,			   
				    Uvec& doppler_axis,
				    CFmat& echo2D_rc,
				    CFmat& echo2d_rdc);

void simpleRangeDopplerCompression3(const CFmat& echo2D_c,
				   const CFvec& fft_matched_filter,
				   const Uvar& prf,
				   const Uvar& bore_doppler,
				   const Uvar& adc,
				   const Uvar& range_start,
				   const double& oversample_ratio,
				   Uvec& range_axis,			   
				   Uvec& doppler_axis,
				   CFmat& echo2D_rc,
				   CFmat& echo2d_rdc);
				
//copied from B.S. code
double getFastDfDt(double look[3], 
		   double sc_acc[3], 
		   double lambda_in_km,
		    double speed_in_km_s,
		   double range_in_km, 
		   double fdop);

//------------------------------------
// functions for echo simulation: (Y. Zhang)
//------------------------------------
Uvar echo_phase(const Uvar& tr, const Uvar& range,
		const Uvar& cfs, const Uvar& csd, 
                const Uvar& csf, const Uvar& pri); 
double solid_pulse(const Uvar& tr, const Uvar& range, 
		   const Uvar& pri, const Uvar& tau_p);
Uvar doppler4target(const Uvar& tp, const string& target_name, 
		    const Uvar& lat, const Uvar& lon); 
void SCState(const Uvar& tp, const string& target_name, 
	     Uvec& Rsc, Uvec& Vsc);
Uvar fitRange(const Uvec& Rsc, const Uvec& Vsc,
	      const Uvec& Rt, const Uvar& tp);
void SetTargets(Time& t, const string& target_name, 
		Uvec& latb, Uvec& lonb, Uvec& Xtb, Uvec& Ytb,
		Uvec& Ztb, const int& key);
void getBurstEcho(Ieb& ieb, const Time& epoch, const Uvar& t_x, 
		  const string& target_name, const Uvec& Rt, 
		  Uvec& rw_time, Uvec& echo);

//---------------------------------
//compute pulse segmentation
//--------------------------------
void computeStartEndIndices( Uvar& min_range,
			    const Uvar& max_range,
			    const Uvar& pri,
			    const Uvar& taup,
			    const Uvar& adc,
			    const unsigned int& pul,
			    const Uvar& rwd,
			    const unsigned int& Ndata_pri,
			    const unsigned int& Nradar_data,
			    unsigned int& index_start,
			    unsigned int& index_end,
			    unsigned int& Npulse);
#endif
