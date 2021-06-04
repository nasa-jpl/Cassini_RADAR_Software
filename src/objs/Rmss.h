//-----------------------------------------------
//Rmss.h
// This file contains the Rmss class declaration.
// The Rmss class provides an interface to obtain information on Instrument
//    Expanded Block including SLOW/FAST fields from IEB data files.
// This program will create and read a file containing IEB values in increasing
// time
//
// This class will be shared with CASSINI_RADAR_RMSS group (GH). 
//
//
// class Rmss;
//
//
//   Construction:
//  Rmss(const string& filename,const string& mode,const string& filetype)
//
//
//   Predicates:
//
//   Get/set methods: 
//    
//   Supporting functions
//--------------------------------------------------------



#ifndef Rmss_H
#define Rmss_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class Rmss;
#include "Io.h"
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Ieb.h"
using std::string;


//-------------------------
// Class Ieb declaration
//-------------------------
class Rmss
  {
  public:
    

    //------------------------
    //typedef for vector of IEB
    //------------------------
    typedef Array1D<Ieb> Iebvec;

    //---------------
    //construction
    //---------------
    Rmss(const string& filename,const string& filetype) throw(ErrorMessage);
 

    //----------------------
    //I/O
    //---------------------
    void loadIeb(const Ieb& ieb)throw(ErrorMessage);
    void readRecordsfromFile() throw(ErrorMessage);
    void writeRecordstoFile() throw(ErrorMessage);  
    void close();
    //--------------------------
    //Supporting function
    //---------------------------
    void clear() throw(ErrorMessage);

    //------------------------------
    //RMSS public data fields
    //------------------------------
   
    Uvar prf;
    Uvar tau_p;
    Uvar BR;
    Iebvec ieb_sequence;
    //---------------------
    //slow field
    //---------------------
    Uvar s_tfi        ;//time from ieb trigger 
                       //16 bits,1sec LSB,range 0-18.20417 H
    unsigned int s_typ;//instruction type 2 bits 11_2 fixed
    unsigned int s_dtn;//data take number 8 bits 0-255
    unsigned int s_sin;//instruction number 8 bits 0 - 255
    unsigned int s_mod;//radar mode 4 bits
    unsigned int s_csr;//calibration source 4 bits
    Uvar s_adc;//adc sample rate 2 bits
    Uvar s_rcv;//receiver bandwidth 2 bits
    unsigned int s_tro;//transmitt/receive window offset 4 bit, 1 pri,-8 to 7
    unsigned int s_baq;//baq mode 3 bits
    unsigned int s_bem;//beam mask 5 bits
    unsigned int s_at1;//receiver attenuation for beam1 and 2, 12 bits
    unsigned int s_at3;//receiver attenuation for beam3,12 bits
    unsigned int s_at4;//receiver attenuation for beam4 and 5, 12 bits
    Uvar s_rip;        //radiometer integration time, 4 bits, 5 ms LSB 10-75 ms
    unsigned int s_rad;//radiometer window count 8 bits range: 1-255
    Uvar  s_csd; //chirp step duration (8 bit, 133.333 ns LSB)
    unsigned int s_csq;//chrip step quantity 12 bits (range: 2 -750)
    Uvar s_cfs;        //chirp frequency step size (16 bits, 1.788 Hz LSB)
   
    //------------------------------------------------
    //fast field
    //-------------------------------------------------
    Uvar f_tfi;        // time from ieb trigger 16 bits 1 sec LSB,
                       //0-18.20417 Hour
    unsigned int f_typ;//instruction type 2 bits 10_2(fixed)
    unsigned int f_fin;//instruction number 8 bits 0-255
    unsigned int f_bil;//busts in instruction 8 bits 1-255
    unsigned int f_pul;//pulses per burst 8 bits 0-255
    Uvar f_bpd;        //burst period 12 bits 1ms LSB, 10-4095
    Uvar f_rwd;        //receive window delay, 10 bits,
                       // one pri LSB, 0-1023 PRI
    Uvar f_pri;        //pulse repetition interval
                       // one ADC sample clock period range 0-1023
                       //only even values of pri_number is valid
                       //ADC sample rate
                       //250 kHz  pri_number * 4 micro second
                       //1 MHz    pri_number * 1 micro second
                       //2 MHz    pri_number * 0.5 micro second
                       //10 MHz   pri_number * 10 * 0.1 micro second
    Uvar f_csf;        //chrip start frequency 16 bits, 457.764 Hz LSB

    //--------------------
    //Max number of record
    //---------------------
    const static unsigned int Nmax = 7200;//number of seconds in 2 hr
  private:
   
  
    FileMgr fp_;
    string filename_;
    string ft_;
    Time t_;
    unsigned int ieb_counter_;
  };


#endif






