//-------------------------------------------------------------------
// RasProc.h
// simple range and doppler compression for doppler centroid tracking
//-------------------------------------------------------------------
                 

#ifndef RasProc_H
#define RasProc_H

#include <string>
//----------------------------
//Forward declaration
//----------------------------
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "RasList.h"
#include "Ieb.h"



//-------------------------
// Class RasProc declaration
//-------------------------
class RasProc: public RasList
  {
  public:
    //---------------------
    //constructor
    //--------------------
    RasProc(const string& filename, const string& mode);

    ~RasProc();//destruction
    //------------------
    //load ieb
    //--------------------
    void loadIeb(Ieb& ieb);


    //--------------
    //get valid number of pulses
    //----------------
    unsigned int getValidNumberOfPulses();

    //------------------------
    //simple matched filter and range/azimuth image
    //----------------------
    void computeMatchedFilter(const Uvar& doppler_centroid);   
    void convert2Decho(const Uvar& range_start,  const Uvar& range_end);
    void convert2Decho1PRI(const Uvar& range_start,  const Uvar& range_end);
    void computeRangeAzimuthCompression(); 
    //compute fractional doppler centroid
    double computeFractionalDoppler();
    void computeRangeFractionalDoppler(Uvec& range, 
				       Dvec& fract_doppler, 
				       double&  mean_fract_doppler);

    //-------------------------
    //display
    //--------------------------
    void displayMatchedFilter();
    void displayRangeCompressedData();
    void displayAzimuthData();  
    void displayRangeDopplerImage();

  
    //matched filter
    CFvec matched_filter;//complex echo(t)
    CFvec fft_matched_filter;
  
    //range compressed data
    CFmat echo2D_rc;//range only
    CFmat echo2D_rdc;//range and doppler

    //----------------------------------
    //fractional doppler container
    //----------------------------------
    CFmat c_corr2D;
    CFvec c_corr1D;
    complex<float> c_corr0D;
    Dvec dop_frac1D;
    double dop_frac0D;

    //range and doppler axis  
    Uvec range_axis;
    Uvec doppler_axis;

  private:
    bool ieb_loaded_;
    bool matched_filter_set_;
    bool echo2D_set_;
    Uvar doppler_centroid_;
    Uvar range_start_;
    Uvar range_end_;
    Ieb ieb_;
    unsigned int N_samples_per_window_;
    unsigned int N_pulses_inside_echo_;
    unsigned int N_samples_inside_echo_;
  };
#endif

