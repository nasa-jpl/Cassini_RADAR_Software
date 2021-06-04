//-------------------------------------------------------------------
// DopplerCentroid.h
// This class will not be used.
//-------------------------------------------------------------------
                 

#ifndef DopplerCentroid_H
#define DopplerCentroid_H

#include <string>
//----------------------------
//Forward declaration
//----------------------------
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Frame.h"



//-------------------------
// Class DopplerCentroid declaration
//-------------------------
class DopplerCentroid
  {
  public:
    //---------------------
    //constructor
    //--------------------
    DopplerCentroid();

    //set method: once for all
    void setLookDirection(const string& look_direction);
    void setTarget(const string& target);
    void setWavelength(const Uvar& lambda);
    void computeAzimuthSteeringAnglesAtEpoch(const Time& epoch_time);
    void setFittingParameters(const bool& fit_yaw,
			     const bool& fit_pitch);

    void reset();//reset container and wavelength  

    // load range and doppler data
    void loadCassini5BeamRecords(const vector<StateVector>& sc_state_record,
		    const vector<unsigned int>& beam_record_indicator,
		    const vector<unsigned int>& range_doppler_record_length,
		    const Umat& range_values,
		    const Umat& doppler_values);

    //void loadRangeDopplerData(const StateVector& sc_state,
    //		      const unsigned int& beam_id,
    //		      const Uvec& range,
    //		      const Uvec& doppler);
    
    //void loadRangeDopplerDataInkm_Hz(const StateVector& sc_state,
    //			     const unsigned int& beam_id,
    //			     const Dvec& range_in_km,
    //
    //	     const Dvec& doppler_in_KHz);

    //compute yaw, pitch angle
    void computeYawPitch(Uvar& yaw,
			 Uvar& delta_yaw,
			 Uvar& pitch,
			 Uvar& delta_pitch,
			 const bool& show_plot);
  
   

    //utility function
    void generateCkernelWithBias(const Time& start_time,
				 const Time& end_time,
				 const Uvar& time_step,
				 const Uvar& yaw_bias_from_fsc_to_ftcn,
				 const Uvar& pitch_bias_from_fsc_to_ftcn,
				 string& ck_filename);

    void generateCkernelTimeVaryingBias(const Time& start_time,
					const Time& end_time,
					const Time& epoch_time,
					const Uvar& time_step,
					const Dvec& yaw_coeff,
					const Dvec& pitch_coeff,
					const vector<Uvar>& ca_time,
					const Dvec& ca_diff_coeff,
					string& ck_filename);

    Uvar computeDopplerAtZeroAzimuth(const StateVector& sc_state_target_frame,
				     const unsigned int& beam_id,
				     const Uvar& range);
    void computeDopplerAtZeroAzimuth(const StateVector& sc_state_target_frame,
				     const unsigned int& beam_id,
				     const Uvec& ranges,
				     Uvec& dopplers);
    void getDopDerivativeYPAZ(Uvar& dfdy, Uvar& dfdp, Uvar& dfdasa);

  private:
    //things to be set for a flyby
    bool look_direction_set_;
    bool target_radius_set_;
    bool azimuth_steering_angle_set_;
    bool wavelength_set_;
    bool fitting_parameter_set_;
    bool l_chisq_;
   
    int x_direction_to_Px_;
    unsigned int Nbeam_;
    double  rscale_;
    Time t_ypaz_;
    Time t_beam3_;
    
    
    //flyby dependent parameters
    Frame ftarget_;
    Frame fj2000_;
    Frame fsc_;
    double lambda_;
    double look_direction_;
    double radius_;
    double href_;//terrain height=0

    //final result
    double yaw_bias_rad_;
    double pitch_bias_rad_;
    Array1D<int> i_paramest_;
    Array1D<int> i_usedata_;   
    Array1D<int> i_asaa_;

    //beam dependent information
    Array1D<double>  yaw_rad_;
    Array1D<double>  pitch_rad_;
    Array1D<double> roll_rad_;
    Dmat  position_;
    Dmat velocity_;
    Array1D<unsigned int> beam_record_indicator_;
    Array1D<double> beam_first_range_;
    Array1D<double> beam_delta_range_;
    Array1D<double> azimuth_angles_;
    Array1D<Uvar> azimuth_bias_angles_;
   
    //for each range and doppler information
    vector<unsigned int> record_beam_id_;
    vector<Time> record_time_;
    vector<double> range_;//data
    vector<double> doppler_;//data
    vector<double> doppler_ck_;//simulated doppler based  on ckernel
    vector<double> new_doppler_;
    vector<double> diff_doppler_;
 
    //need for svdvecfit
    double r_sch_[3];
    double r_schvel_[3];
    double r_f_, r_dfdy_,r_dfdp_,r_dfdvs_;
    double r_dfdvc_,r_dfdvh_,r_dfdhref_,r_dfdwvl_;
    double r_dfdasa_;
    unsigned int i_yaw_;
    unsigned int i_pitch_;
    unsigned int NPP_;
    unsigned int MP_;
    unsigned int i_mp_;
    unsigned int i_rd_;
    unsigned int i_fp_;
    unsigned int i_np_;
    Dmat r_vecin_;
    Dmat r_vobs_;
    D3D r_cov_;
    Dvec r_a_;
    Dvec r_at2_;
    Dmat r_u_;
    Dmat r_v_;
    Dvec r_w_;
    Dmat r_chisq_;
  
    Dmat r_cvm_;
  

    //finally, keep the record of time, yaw, pitch
    // yaw_bias, pitch_bias
    vector<Time> beam3_time_record_;
    vector<Uvar> yaw_record_ftcn_to_fsc_;
    vector<Uvar> pitch_record_ftcn_to_fsc_;
    vector<Uvar> roll_record_ftcn_to_fsc_;
   
  
    void computeTRMatrix(const StateVector& sc_state,
			 const unsigned int& beam_id);
    void displayFitting();
 
  };

//----------------------------------
//Relative integer tracking function
//-----------------------------------
void track_rel_int_doppler(const vector<unsigned int>& sab_record,
			   const vector<Uvar>& prf_record,
			   const vector<double>& fract_record,
			   const vector<Uvar>& geom_doppler,
			   const vector<int>& geom_int_doppler,
			   vector<int>& rel_int_doppler,
			   const bool& show_plot);
#endif
  
