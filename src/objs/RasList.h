//-----------------------------------------------
//RasList
// This file contains the RasList  class declaration.
// handle intermediate file format for S. S.
//--------------------------------------------------


#ifndef RasList_H
#define RasList_H

//----------------------------
//Forward declaration
//----------------------------
class RasList;

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Io.h"
#include "Constants.h"
#include "Sab.h"

using std::string;
//-------------------------
// Class RasList declaration
//-------------------------
class RasList
  {
  public: 
    //--------------
    // Construction
    //-------------- 
    RasList(const string& filename, const string& mode); 

    ~RasList();
    //-----------------------------------------
    //Input method load sab from L0 or echo simulation program
    //-----------------------------------------
    void loadSab(const Sab& sab);
    void loadEncodedData(Ieb& ieb, 
			 const Time& transmit_time,
			 const unsigned int& burst_id,
			 const unsigned int& beam_id,
			 const Array1D<unsigned short>& baq,
			 const Array1D<unsigned char>& data);


    //--------------------------------
    //load ieb for ras data processing and checking
    //------------------------------    
    void decodeEcho( Ieb& ieb);
    void decodeBaqThresholds();

    //------------------------
    //find a  record
    //------------------------  
    bool findRecord(const unsigned int& n_sab);
    bool findRecord(const Time& t);

    //----------------------------
    //record manipulation
    //----------------------------
    void readFirstRecord();
    void readNextRecord();
    void skipOneRecord();
    unsigned int& operator[] (const unsigned int& n);//return sab number of nth record
    //--------------------------
  //compute record number
    //---------------------------
    unsigned int computeNumberOfRecords();
    void getFirstLastSabNumber(unsigned int& first_sab, unsigned int& last_sab);
   
    //---------------------
    //file position control
    //---------------------
    void rewind();
    bool eof();

  
    //-----------------------------
    //convert 1D data to 2D data
    //-------------------------------
    void getEchoData(Array1D<float>& echo1D);
    //void shiftEcho(const unsigned int& dn);//dt is supposed to be computed outside
   
    void get2Decho(const unsigned int& start_index,
		   const unsigned int& end_index,
		   const unsigned int& N_samples_per_window,
		   const unsigned int& N_pulses_inside_echo);
    void directFFT();
   
  
  
    //--------------------
    //plot echo, matched filter power spectrum
    //---------------------
    void displayEcho();
    void compare_baq_exponent(const Uvar& pri, const unsigned int& pul,
			      const unsigned int& tmp_baq,
			      const Uvar& adc,
			      const int& tro_int);
    void display2Decho(const unsigned int& n_pulse);
    void displayAll2Dechoes();
    void displayComplexPulse(const unsigned int& n_pulse);
   

    //--------------------------
    //perform bit analysis
    //--------------------------
    void getEchoStatistics(double& standard_deviation,
			   double& average);
			   
    void showSignBitStatistics(const unsigned int& n_sab);

    //------------------------------------------
    //public data folder    
    //These parameters are extracted from SAB
    //-------------------------------------------
    Time t;
    unsigned int sync;
    unsigned int sab_counter;
    unsigned int sclk;
    unsigned int brst;
    unsigned int beam_number;
    unsigned int baq_mode;
    unsigned int radar_mode;
   
    unsigned int valid_data_per_pri;
    unsigned int valid_data_per_burst;
    unsigned int radar_data_size_in_bytes;

    //public data folder
    // These are constructed after decoding data
    float offset;
    
    Array1D<float> echo_data;//decoded data (real, no imaginary part)
    CFvec fft_echo_data;
    Array1D<float> abs_echo_data;//absolute value of real data
    Array1D<unsigned int> baq_decoded;//decoded baq
    Array1D<Uvar> sign_bit;

    //2D formatted data
    Array2D<float> echo2D;
    CFmat echo2D_c;
  
    unsigned int pos;
    unsigned int neg;

    static const  unsigned int Nblock_ = 24;
    //-----------------------
    //Priviate variables/methods
    //-----------------------
  private:
    FileMgr file_;
    string filename_;
    string mode_;
   
    bool decoded_;
    bool echo2D_set_;
    vector<unsigned short>  baq_;//baq stored in 16 bit words
    vector<unsigned char>  rdata_;//compressed data
   
    //----------------------------------------------------------
    //useful record containter that holds sab, time
    // file position
    //these containers are set onece when rasfile is opened
    //-----------------------------------------------------
    //list of sab number
    vector<unsigned int> sabnumber_list_;
    vector<unsigned int>::const_iterator  sabnumber_list_ptr_;
    //sab number vs file position
    map<unsigned int, int> sabnumber_vs_fileposition_;
    map<unsigned int,int>::const_iterator sabnumber_vs_fileposition_ptr_;
    //list of time for each record
    vector<Time> time_list_;
    vector<Time>::const_iterator time_list_ptr_;
    //time vs sab record
    map<Time, unsigned int> time_vs_sabnumber_;
    map<Time, unsigned int>::const_iterator time_vs_sabnumber_ptr_;
   

    
    unsigned int ras_sab_size_before_buffer_;
    unsigned int ras_sab_size_after_sab_before_buffer_;
    unsigned int ras_up_to_sab_size_;
    unsigned int Ndop_bins_;
    unsigned int Nrange_bins_;
    
    void mapAllRecords();
    void readRecord();
    void writeRecord();
    void getComplex2Decho();   
    
  };


#endif




/*
//void convert2pow(unsigned int& n);
  void get2Decho2Pri();
    void avgPulseAmp();   
   
    void getComplex2Decho_d();
    void pulseCorrelation();
    void avgPulseCorrelation();
    void get1DPhaseDiff();
    void getPhaseDiff();
    void phaseUnwrapp();
    void averagePhaseDiff();
    void fractDopCentroid();
    void IntegDopCentroid(const Uvar& bore_dop);
    void TotalDopCentroid();
    void avgTotalDopCentroid();
    void FractDopCentroid();



  //doppler centroid data
    CFmat pulse_corr;
    CFvec avg_pulse_corr;
    Array2D<float> phase_diff;
    Array1D<float> avg_phase_diff;
    Array1D<float> fDOPrePRF;
    Array1D<float> pulse_amp;
    Uvec totalDop;
    complex<float> avgTotalCorr;
    Uvar avgTotalDop;
    Uvar fractAvgTotalDop ;
    int int_m;

    void displayComplexPulse();
 
    void plotDopCentroid();
    void plotPhaseDiff();
    void plotAvgPhaseDiff();
    void plotPulseCorrelation();
    void plotTotalDop();
    void plotAvgPulseAmp();
    void compareBoreSightDop(const Uvar& bore_dop);





*/
