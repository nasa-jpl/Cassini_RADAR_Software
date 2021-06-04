//-----------------------------------------------
//IebList.h
// This file contains the IebList class declaration.
// The IebList class provides an interface to handle IEB record file such as
// writing and reading
//
// This class will be shared with CASSINI_RADAR_RMSS group (GH). 
//
//
// class IebList;
//
//
//   Construction:
//  IebList(const string& filename,const string& filetype) throw(ErrorMessage)
//
//-----------------------
//supporting functions
//----------------------
// void load(const Ieb& ieb);

//--------------------------
//I/O
//--------------------------
//void readRecord(const Time& t, Ieb& ieb) throw(ErrorMessage);
//void writeRecord() throw(ErrorMessage);
//--------------------------------------------------------



#ifndef IebList_H
#define IebList_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class IebList;
#include <vector>
#include <iterator>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Ieb.h"
#include "Io.h"
#include "Ambiguity.h"
#include "Constants.h"
#include "IebProfile.h"
using std::string;
//-------------------------
// Class Ieb declaration
//-------------------------
class IebList
  {
  public: 
    
    //---------------------
    //typedef's and enums
    //---------------------
    typedef std::vector<unsigned short> sdata;//short integer data
    typedef std::vector<Ieb>::const_iterator Ieb_ptr;
  
    //--------------
    // Construction
    //-------------- 
    IebList(Config& cfg); 
    IebList(const Uvar& start_time, 
	    const Uvar& end_time,
	    Config& cfg);
   
    ~IebList();//destructor
    
    //-----------------------
    //config set up
    //----------------------
    void config();
   
    //------------------
    //ambiguity container:
    //-----------------
    vector<AmbiguityFNN> amb;
    vector<unsigned int> nlooks;
    vector<Uvar> ML_area;

    //-----------------------
    // Beam related....
    //{get} method to obtain 3dB points , beam objects
    // and prf and dutycycle polyfit coeff
    //-----------------------
    void computeBeam3dB_5dBpoints();
    void getBeam3dB_azimuth_points_at_0_elevation(Umat& azimuth_3dB);
    void getBeam3dB_elevation_points_at_0_azimuth(Umat& elevation_3dB);
    void computeRangeDoppler_3dB_5dB_Borepoints(const Time& t);
    

    Umat getBeam3dB_azimuth(); //Included by Chandini
    //------------------------
    //return the size of ieb records, start time, and end time
    // indexing
    //------------------------
    Time getStartTime();
    Time getEndTime();
    Time getFirstValidSlowFastFieldTime();
    unsigned int getNumberofIebRecords();
    Ieb& operator[] (const unsigned int& i);

    //--------------------------------
    //adding and  deleting ieb record at time t
    //-------------------------------
    void addIeb(const Ieb& ieb);
    void back_insert(const Ieb& ieb, bool& valid);//add one ieb at the end of ieblist
    void deleteIeb(const Ieb& ieb);
    Ieb getIeb_v0(const Time& t) ;//return valid ieb at time of t1 <= t
    Ieb getIeb(const Time& t) ;//return valid ieb at time of t1 <= t
    Ieb getLiveIvpUpdateIeb(const Time& t, const Uvar& delta_trigger);
    Ieb getIebOn(const Time& t);//retrun valid ieb at time of t1==t only
    //GH: add method to delete block of time.  see if this is faster.
    void deleteIebBlock(const Time& goodtime, const unsigned int& duration) ;
  
    //---------------------------
    //Generate Ideal IEB sequence
    //----------------------------
    void generateIdealIebSequence();

    //----------------------
    //Supporting function using Geometry/SC_attitude
    //-----------------------
    void computeSarRadarEcho(const unsigned int& beam_number,
			     const Time& t,
			     const Uvar& taup,
			     Uvar& echo); 
    //ieb is updated as needed
    void computeUsableCrosstrackExtentandGap(Uvar& extent,
					     Uvar& gap,
					     const Time& t, 
					     const Uvar& time_offset,
					     const Frame& track_frame);
    //what if ieb is fixed while geometry is being updated
    void computeUsableCrosstrackExtentandGap(Uvar& extent,
					     Uvar& gap,
					     const Time& t,
					     const Ieb& ieb135,
					     const Uvar& time_offset,
					     const Ieb& ieb24,
					     const Frame& track_frame);
    //return ambiguity output
    void getAmbiguity(vector<AmbiguityFNN>& amb_output);
    //compute azimuth and range resoltuion
    void computeAzimuthRangeResolution(Uvar& x_res,
				       Uvar& r_res, 
				       Uvar& rg_res,
				       const Time& t);

    //compute special chirp start frequency
    Uvar  findToneChirpStartFrequency(const Uvar& pri,const Uvar& adc);
    
    //compute looks at t
    //unsigned int computeLooks(const Time& t);
    //compute looks at t but with a fixed iebg
    //unsigned int computeLooks(const Time& t, const Ieb& ieb);
    //what is the minimum bpd to produce a desired number of looks
    //Uvar computeBpd(const Time& t, const unsigned int& looks);

    //-----------------------------------------
    //I/O : read/write IEBs from/into a file
    //-----------------------------------------
    void readAllIebRecords(const string& filename, const string& filetype);
    void writeAllIebRecords(const string& filename, const string& filetype);
    void getAllIebRecords(vector<Ieb>& target_ieb_list);

    //--------------------------
    //GH:I/O for RMSS (testbed)
    //--------------------------
    void writeAllwtkIebRecords(const string& filename, const int base);
    // will set BII, FIN, SIN
    unsigned int computeBii () ;
    unsigned int compressSlowfield(const int percentage) ;

    //-----------------------------
    //located PRF hopping
    //-------------------------------  
    void locatePrfHopping(Array1D<Time>& hopping_time);//transfer hopping time
    void locatePrfHoppingDivision(Array1D<Time>& prf_start_time,
				  Array1D<Time>& prf_end_time,
				  unsigned int& number_of_divisions);

    ///Included by CV///
    //----------------------------------
    // Calculate Range/Doppler Ambiguity
    //----------------------------------
    void computeDelta_Dop_Rng( Uvec& dop_delta, Uvec& g_azi, Uvec& range_delta, Uvec& g_elev, const Uvar& t);//Included by Chandini
    void gamb_ratio_dop_rng(Uvar& amb_ratio_dop, Uvar& amb_ratio_range, const Uvec& dop_delta, const Uvec& g_azi, const Uvec& range_delta, const Uvec& g_elev, const Uvar& t, const Uvar& fp_try);//Included by Chandini

    
    //--------------------------------------------
    //constant members
    //size of slow and field fields in units of word(16 bit)
    //Slow/Fast/Power/Tnc
    //---------------------------------------------
    static const unsigned int Nslowfield = 10;
    static const unsigned int Nfastfield = 7;
    static const unsigned int Npowerfield = 2;
    static const unsigned int Ntncfield = 5;
    
  private:
    //---------------------------
    //internal representation
    //---------------------------  
    bool allbeams_set_; 
    bool set_csf_times_pri_integer_;
    //-----------------
    //epoch time information
    //-----------------
    Time epoch_time_;

    //---------------------
    //slow and fast fields
    // this slow and fast field are used as container to store
    // them before/after writing/reading into file
    //---------------------
    sdata  slowfield_;   
    sdata  fastfield_;   
    sdata  powerfield_;
    sdata  tncfield_;
    
    //----------------
    //Internal ieb 
    //-----------------  
    vector<Ieb> ieb_list_;
    unsigned int data_take_number_;
    //-----------------
    //Internal frame and target name
    //-----------------
    string target_name_;
    Frame target_frame_;

    //--------------------
    //internal copy of config object
    //-------------------- 
    Config cfg_;//need to have an internal copy for ambiguity calculation using ieb in ieblist
    Uvar start_time_;
    Uvar  end_time_;
    
    
    //------------------------
    //Contain beam objects, four 3dB points, two in the azimuth axis at 0 elev
    // and the other two in the elevation axis at 0 azimuth 
    // pri and dutycycle polyfit
    //in order to pass objects and coefficients to ieb class
    //------------------------
    vector<Beam> Bvector_;
    vector<Frame> Fvector_;
    Umat azimuth_3dB_;
    Umat elevation_3dB_;
    Umat azimuth_5dB_;
    Umat elevation_5dB_;

    Umat range_3dB_bore_;
    Umat doppler_3dB_bore_;
    Umat range_5dB_bore_;
    Umat doppler_5dB_bore_;
 
    // ----------------------------
    // Ambiguity Calculation
    //  g_amb (Included by Chandini)
    //-----------------------------

    //Uvar amb_ratio_dop_;
    //Uvar amb_ratio_range_; 
    Uvar fp_try_;
    Uvar g_doplevel_;
    //Uvec dop_delta_;
    //Uvec g_azi_;
    Uvar g_rangelevel_;
    //Uvec range_delta_;
    //Uvec g_elev_;
    //int Nphi_;

    //---------------------
    //system temperature, radar power
    //squared_deviation_of_system_noise_cal
    // noise level when calibration was performed
    //---------------------
    Uvar Tsys_;
    Uvar Pt_;
    Uvar squared_deviation_of_system_noise_cal_altl_;
    Uvar squared_deviation_of_system_noise_cal_alth_;
    Uvar squared_deviation_of_system_noise_cal_sarl_;
    Uvar squared_deviation_of_system_noise_cal_sarh_;
    Uvar Pn_cal_altl_;
    double radar_buffer_std_deviation_altl_;
    double cal_gain_altl_;

    //------------------------------------
    //Radar mode dependent settings to be used for the 
    // generation of ideal IEBs
    //-------------------------------------      
  

    //------------------
    //parameter needed to calculate usable cross track extent
    //-------------------
    double pbw_ratio_;
    double muhleman_k1_;
    double muhleman_k2_;
    double amb_ratio_;
    double noise_equiv_sigma0_;
    double min_beam_gain_;

    //----------------
    //parameter used for tone data
    //---------------
    Uvar min_frequency_of_returned_cw_echo_;
    Uvar max_frequency_of_returned_cw_echo_;

  private:
    Uvar computePRF(const Uvar& t);//Included by CV
  };


#endif






