//-----------------------------------------------
//Ieb.h
// This file contains the Ieb class declaration.
// The Ieb class provides an interface to obtain information on Instrument
//    Expanded Block including SLOW/FAST fields from IEB data files.
// IEB sequence table in IEB data file will be created such that 
//    successive parameter sets increase with time and have a difference of at 
//    least one second (see Blue Book 3-18)
//
// This class will be shared with CASSINI_RADAR_RMSS group (GH). 
//
//
// class Ieb;
//
//
//   Construction:
//--------------
// Construction
//-------------- 
//    Ieb(const Time& t);
//    Ieb();
//   
// 
//
//---------------------------
//setup
//----------------------------
//    void config(Config& cfg) throw(ErrorMessage);
//      
//----------------
// Get/Set methods
//----------------
//    Time getTime() const  throw(ErrorMessage) ;
//    void setTime(const Time& t) throw(ErrorMessage);
//
//
//--------------------
//Primary inputs: 
//---------------------
//    void setTarget(const string& target_name, const Frame& ftitan)
//                   throw(ErrorMessage);   
//    void setTarget(const string& target_name) throw(ErrorMessage);
//    void setPrf() throw(ErrorMessage);
//    void setPulsewidth(const double& dutycycle) throw(ErrorMessage);
//    void setReceivewindow(const unsigned int& add_pri) throw(ErrorMessage);

//---------------------
//SupportingFunction
//---------------------
//    void generateSlowfield() throw(ErrorMessage);
//    void generateFastfield() throw(ErrorMessage);
//    void encode() throw(ErrorMessage);
//    void decode(const std::vector<unsigned short>& slow, 
//                const std::vector<unsigned short>& fast) throw(ErrorMessage);
//    void loadSlowfield(const std::vector<unsigned short> slow) 
//              throw (ErrorMessage);
//    void loadFastfield(const std::vector<unsigned short> fast)
//             throw(ErrorMessage);
//
//
//--------------------------------------------------------



#ifndef Ieb_H
#define Ieb_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class Ieb;
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "CAGTable.h"

using std::string;


//-------------------------
// Class Ieb declaration
//-------------------------
class Ieb
  {
  public:
    
    //---------------------
    //typedef's and enums
    //---------------------
    typedef std::vector<unsigned short> sdata;//short integer data

    //--------------
    // Construction
    //-------------- 
    Ieb(const Time& t);
    Ieb();
   
    
    //---------------------------------------
    // {get,set} time
    //--------------------------------------
    //{get,set} internal variable time
    Time getTime() const ;
    void setTime(const Time& t) ;
  
    //----------------------------
    //{get,set} ilxtype
    //4: every mode is valid (for ideal param generation)
    //3: slow/fast valid
    //2: only fast valid
    //1: only TNC valid
    //0: only power mode valid
    //reset method is OK at any time
    //-----------------------------
    unsigned int getIebMod() const;
    void setIebMod(const unsigned int& ilxtype);
    
    //--------------------------
    //Methods to access slowfield/fastfield parameters
    //-------------------------

    //{get,set} fastfield tfi
    Uvar getFastTfi() const ;
    void setFastTfi(const Uvar& tfi);
    
    //{get} fastfield typ
    unsigned int getFastTyp() const;

    //{get,set} fastfield fin
    unsigned int getFin() const ;
    void setFin(const unsigned int& fin);

    //{get,set} fastfield  bursts in instruction
    unsigned int getBii() const ;
    void setBii(const unsigned int& bii);

    //{get,set} fastfield pul
    unsigned int getPul() const ;
    void setNumberOfPulses(const unsigned int& Np);
    Uvar getPulseTransmitTime() ;
   

   
    void computeSAR_Pul_Rwd_Tro_Csf(const Uvar& bpd, const double& data_rate,
				    const Umat&  range,const Umat& doppler);
    void computeSAR_Pul_Rwd_Tro_Csf(const Umat& range,const Umat& doppler);
    void computeSAR_Pul_Rwd_Tro_Csf2(const Umat& range,const Umat& doppler);
    void computeSAR_Pul_Rwd_Tro_Csf3(const Umat& range,const Umat& doppler);
    void computeALT_Rwd_Tro_Csf(const Umat& range,const Umat& doppler);
    void computeALT_Rwd_Tro_Csf(const Umat& range,const Umat& doppler, const int& tro);
    void computeToneALT_Rwd_Tro(const Umat& range, const Uvar& chirp_start_freq);
    void computeToneALT_Rwd_Tro(const Umat& range, const Uvar& chirp_start_freq, const int& tro);
    //{get,set} fastfield burst period
    Uvar getBpd() const;
    unsigned int getBpdEncodedValue() const ;
    void setBpd(const Uvar& bpd);
    void increaseBpd (const Uvar& increase_bpd);
    //void computeBpd(const double& data_rate);

    //{get,set} fastfield  receive window delay
    Uvar getRwd();
    void setRwdInUnitsOfPri(const unsigned int& rwd);
    unsigned int getRwdInPriUnits() const;
    void setRwd(const Uvar& rwd);
    
    //{get,set,compute} fastfield PRI and dutycycle
    Uvar getPri() const;
    unsigned int getPriEncodedValue() const;
    void setPri(const Uvar& pri);
    void computePri(const Uvar& altitude,const Dvec& prf_alt_fit);
    void computeDutycycle(const Uvar& altitude,const Dvec& dutycycle_alt_fit);
    double getDutycycle() const;
    void setDutycycle(const double& dutycycle);

    //{get,set} fastfield csf: chirp start frequency
    Uvar getChirpStartFrequency() const ;
    unsigned int getChirpStartFrequencyEncodedValue() const;
    void setChirpStartFrequency(const Uvar& csf);
    void computeChirpStartFrequency(const Uvar& doppler);
    void computeChirpStartFrequency(const Uvar& min_doppler, 
				    const Uvar& max_doppler);
    
    //{get,set}slow field tfi
    Uvar getSlowTfi() const;
    void setSlowTfi(const Uvar& tfi);
   
    //slow field typ
    // no set method is allowed
    unsigned int getSlowTyp() const;

    //{get,set} slow field dtn
    unsigned int getDtn() const ;
    void setDtn(const unsigned int& dtn);
   
    //{get,set} slow field sin
    unsigned int getSin() const;
    void setSin(const unsigned int& sin);
   
    //{get,set} slow field radar  mod for normal operation
    string getModName() const;
    unsigned int getMod() const;
    void computeMod(const Uvar& adc, 
		    const unsigned int& csr, 
		    const unsigned int& baq, 
		    const unsigned int& bem,
		    const unsigned int& changeOfBW);
   
    //{get,set} slow field radar  mod for calibration
    void setWTKCalibration(const Uvar& adc,
			   const unsigned int& csr,
			   const unsigned int& bem);

    // a copy of WTK but some changes....
    //GH {get,set} slow field radar  mod for calibration
    void setGHCalibration(const Uvar& adc,
			  const unsigned int& csr,
	  		  const unsigned int& changeOfBW);

    //GH edit Csr to correct for Auto-gain FLTSW feature (if enable disable for 1 instruction)
    void disableAutogain();
    
    //{get,set}slow field calibration source
    unsigned int getCsr() const;
       
    //-------------------------------------------------------
    //only get method is avail for slow field adc and rcv 
    // because they are coupled to radar mode selection
    // only changes in radar mode can change these two parameters
    //--------------------------------------------------------
    Uvar getAdc() const;
    unsigned int getAdcEncodedValue() const;
    Uvar getRcv() const;
    unsigned int getRcvEncodedValue() const;

    //{get,set} slowfield TRO
    Uvar getTro() const ;
    int getTroInPriUnits();
    void setTroTime(const Uvar& tro_time);   
    void setTroInUnitsOfPri(const int& tro_number);
    void resetTroBasedOnCTRX(const unsigned int& ctrx);   
    void computeTro(const Uvar& min_range, const Uvar& max_range);
    

    //{get,set} slowfield BAQ
    unsigned int getBaq() const;
  
    //{get,set} slowfield beam mask
    unsigned int getBem()  const;
    unsigned int getNumberOfBeams();
   
    //set attenuator
    //void computeAttenuationSettings(const Uvar& Pn, 
    //	            const Uvec& Echo_power,
    //	            const Uvar& squared_deviation_for_known_source);
   
    void computeAttenuationSettings(const double& noise_bit,
				    const Uvar& Pn,
				    const Uvar& power_to_variance);

    //unsigned int  getSingleBeamAttenuationSetting(const Uvar& Pn, 	
    //		    const Uvar& Echo_power,
    //		    const Uvar& squared_deviation_for_known_source);

    void setAttenuation(const unsigned int& at1, 
			const unsigned int& at3, 
			const unsigned int& at4);
    
    //double radar_buffer_squared_deviation(const Uvar& Power_input, 
    //      const unsigned int& attenuation_in_dB,
    //       const Uvar& squared_deviation_for_known_input_noise_power);

  

    //{get,set} slowfield attenuator1
    unsigned int getAt1();
    unsigned int getAt1N1dB()  const;
    unsigned int getAt1N2dB()  const;
    unsigned int getAt1N3dB()  const;
   
    //{get,set} slowfield attenuator3
    unsigned int getAt3();
    unsigned int getAt3N1dB() const;
    unsigned int getAt3N2dB() const;
    unsigned int getAt3N3dB() const;
   

    //{get,set} slowfield attenuator4
    unsigned int getAt4();
    unsigned int getAt4N1dB() const;
    unsigned int getAt4N2dB() const;
    unsigned int getAt4N3dB() const;

    void getAttenuationBitMask(unsigned int& at1_bit_mask,
			       unsigned int& at3_bit_mask,
			       unsigned int& at4_bit_mask);

    //-------------------------------------
    //get the attenuation setting of a beam
    // two method:
    // beam number: 1 2 3 4 5
    // beam id    : 0 1 2 3 4
    // getAttenuationdBofBeamNumber(const unsigned int& beam_number)
    // getAttenuationdBofBeamId(const unsigned int& beam_id)
    //-------------------------------------
    unsigned int getAttenuationdB(const unsigned int& beam_number);
    double getAttenuation(const unsigned int& beam_number);
    

    //{get,set} slowfield radiometer integration time (RIP)
    Uvar getRip() const;
    unsigned int getRipEncodedValue() const;
    void setRip(const Uvar& rip);
    void computeRip();

    //{get,set} radiometer window count
    unsigned int getRad() const;
    void setRad(const unsigned int& rad);
    void computeRad();//on the assumption that bpd is set
    void computeRadandBpd(const unsigned int& bursts_in_flight);

    //------------------------------------------------
    //computeWaveform: set CSD, CSQ, and CFS altogether
    //-------------------------------------------------
    void computeWaveform();
    void constructWaveform(const Uvar& csd, const unsigned int& csq, const Uvar& cfs);
    Uvar getCsd()  const;
    unsigned int getCsdEncodedValue() const;
    unsigned int getCsq() const;
    Uvar getCfs() const;
    unsigned int getCfsEncodedValue() const;

    //------------------------------------------------
    //GH:resyncTimeTfi: method to reset TFIs after time has shifted
    //-------------------------------------------------
    void resyncTimeTfi (const Time& iebstarttime, const Time& ilxexecutiontime) ;
    //------------------------------------------------
    //GH:defaultRadiometer: method to generate radiometer IEB, 
    // Used: when at Active Range(ALT,SAR,..), but power not available
    //-------------------------------------------------
    void defaultRadiometer (const Time& iebstarttime, const Time& ilxexecutiontime);
    //------------------------------------------------
    //GH:interleave: method to generate default interleaved instruction, 
    // Used: ??? dfined by config file
    //-------------------------------------------------
    void interleave ( const int input_value ) ;
    //-------------------------------------------------
    //GH: to lessen the number of slowfields, average some of the waveforms
    //-------------------------------------------------
    void changeWaveform( const unsigned int& newcsq, const Uvar& newcfs );

    //-----------------------
    //get derived quanity taup = (csq+1)*csd
    //                    br = (csq+1)*cfs
    // get active time offset for doppler computation
    //-----------------------
    Uvar getTaup() const;
    Uvar getBr() const;
    Uvar getOffsetTimefromBurst() const;

    //--------------------------------
    //Method to access power variables
    //--------------------------------
    Uvar getPowerTfi() const;
    void setPowerTfi(const Uvar& tfi);

    unsigned int getPowerTyp() const;

    //power mode
    unsigned int getPmd() const;
    void setPmd(const unsigned int& pmd);

    //----------------------------------
    //Method to access t&c variables
    //----------------------------------
    Uvar getTncTfi() const;
    void setTncTfi(const Uvar& tfi);

    void setRL_ND_AutoRad(const Uvar& rl, const Uvar& max_rl,
			  const Uvar& nd, const Uvar& max_nd);
    void setRL_ND(const Uvar& rl, const Uvar& nd);
    Uvar getRL() const;
    Uvar getND() const;
   
    unsigned int getTncTyp() const;
    unsigned int getTnc () const;
    unsigned short getTca () const;
    unsigned short getTcb () const;
    unsigned short getTcc () const;

    //---------------------------------------------------------------------
    //SupportingFunction to encode/decode/load/get slow/fast/power/tnc fields
    //---------------------------------------------------------------------
    void encodeSlowFastfield() ;
    void decodeSlowFastfield();

    void encodePowerfield();// gh add to split power/tnc
    void decodePowerfield();// gh add to split power/tnc
    void encodeTncfield();// gh add to split power/tnc
    void decodeTncfield();// gh add to split power/tnc

    void encodePowerTncfield();
    void decodePowerTncfield();
    void loadSlowFastfield(const sdata& slow,
			   const sdata& fast ) ;
    void loadPowerfield(const sdata& powerfield);
    void loadTncfield(const sdata& tncfield ) ;

    sdata getSlowfield() const;
    sdata getFastfield() const;
    sdata getPowerfield() const;
    sdata getTncfield() const;

    

   
    unsigned int getPointsInsidePRI();
    unsigned int getRadarDataPoints();
    unsigned int getNumberOfWords();
    double getDataRate();
    void setMaxDataRate(const double& data_rate);
    void RwdAdjustForNBursts (const unsigned int& bursts_in_flight, Uvar& bpd);

    bool encodeSlowFastfieldComplete();
    bool decodeSlowFastfieldComplete();
    bool encodePowerfieldComplete(); // GH: added to split power/tnc
    bool decodePowerfieldComplete(); // GH: added to split power/tnc
    bool encodeTncfieldComplete(); // GH:added to split power tnc
    bool decodeTncfieldComplete(); // GH:added to split power tnc
    bool encodePowerTncfieldComplete();
    bool decodePowerTncfieldComplete();

    static bool  selfTest();

   
    

    //-----------------------
    //set all the flags false
    //-----------------------
    void clear();


    //---------------------------------------
    //constants: word size of each field(slow,fast,power,tnc)
    //           Nmax_buffer: 16K words, 32768 Bytes
    //           header and footer size
    //-----------------------------------------
    static const unsigned int Nslowfield = 10;
    static const unsigned int Nfastfield = 7;
    static const unsigned int Npowerfield = 2;
    static const unsigned int Ntncfield = 5;
    static const unsigned int Nmax_buffer=32768;//32K bytes
    static const unsigned int Nheader = 40;//40 words
    static const unsigned int Nfooter = 22;//22 words


  private:   
    //---------------------------
    //internal representation
    //---------------------------   
    bool time_set_; 
    bool dutycycle_set_;
    bool f_tfi_set_;
    //bool fin_set_; turned off (4/28/03)
    bool pul_set_;
    bool bpd_set_;
    bool rwd_set_;
    bool pri_set_;
    bool chirp_start_frequency_set_;
    bool s_tfi_set_;
    bool dtn_set_;
    //bool sin_set_; turned off (4/28/03)
    bool mod_set_;
    bool csr_set_;
    bool adc_set_;
    bool rcv_set_;
    bool tro_set_;
    bool baq_set_;
    bool bem_set_;
    bool at1N1_set_;
    bool at1N2_set_;
    bool at1N3_set_;
    bool at3N1_set_;
    bool at3N2_set_;
    bool at3N3_set_;
    bool at4N1_set_;
    bool at4N2_set_;
    bool at4N3_set_;
    bool rip_set_;
    bool rad_set_;
    bool csd_set_;
    bool csq_set_;
    bool chirp_frequency_step_set_;
    bool encode_slowfast_complete_;
    bool decode_slowfast_complete_;
    bool slowfastfield_loaded_;
    bool power_loaded_;
    bool tnc_loaded_;

    bool encode_power_complete_; // GH: add to split Pwr/tnc
    bool decode_power_complete_; // GH: add to split Pwr/tnc
    bool encode_tnc_complete_;// GH: add to split Pwr/tnc
    bool decode_tnc_complete_;// GH: add to split Pwr/tnc

    bool encode_powertnc_complete_;
    bool decode_powertnc_complete_;
    bool p_tfi_set_;
    bool t_tfi_set_;
    bool tnc_set_;
    bool tca_set_;
    bool tcb_set_;
    bool tcc_set_;
    bool data_rate_set_;
    //------------------
    //One ieb at time t_
    //-------------------
    Time t_;

    //----------------------------------------
    //string radar mode: SARH_AGC or SARH ...
    //------------------------------------------
    string radar_mode_;
  
    //------------------------------
    //derived quantities
    //IEB parameters not directly specified in slow/fast field
    //------------------------------
    Uvar prf_;
    double dutycycle_;
    Uvar BR_;
    Uvar tau_p_;

    //---------------------
    //ILXTYPE
    // 0 : PWR mode
    // 1 : TNC mode
    // 2 : Fast field only
    // 3 : Slow/Fast field (default)
    //----------------------
    unsigned int ilxtype_;

    //------------------------------------------------
    //fast field
    //-------------------------------------------------
    sdata fastfield_;   //array of 16bit words
    Uvar f_tfi_;      
    unsigned short f_tfi_e_;
    const static unsigned short f_typ_e_;//instruction type 2 bits 10_2(fixed) 
    unsigned short f_fin_e_;//instruction number 8 bits 0-255
    unsigned short f_bii_e_;//busts in instruction 8 bits 1-255
    unsigned short f_pul_e_;//pulses per burst 8 bits 0-255  
    Uvar f_bpd_;        //burst period 12 bits 1ms LSB, 10-4095
    unsigned short f_bpd_e_;
    Uvar f_rwd_;
    unsigned int   f_rwd_e_;        //receive window delay, 10 bits,
    Uvar f_pri_;        //pulse repetition interval
    unsigned short f_pri_e_;
    Uvar f_csf_;        //chirp start frequency 16 bits, 457.764 Hz LSB
    unsigned short f_csf_e_;

    //---------------------
    //slow field
    //---------------------
    sdata slowfield_;//array of 16bit words
    Uvar s_tfi_        ;//16 bits,1sec LSB,range 0-18.20417 H              
    unsigned short s_tfi_e_;//integer value corresponding unit valued s_tfi
    const static unsigned short s_typ_e_;//instruction type 2 bits 11_2 fixed
    unsigned short s_dtn_e_;//data take number 8 bits 0-255
    unsigned short s_sin_e_;//instruction number 8 bits 0 - 255
    unsigned short s_mod_e_;//radar mode 4 bits
    unsigned short s_csr_e_;//calibration source 4 bits
    Uvar s_adc_;
    unsigned short s_adc_e_;//adc sample rate 2 bits
    Uvar s_rcv_;//receiver bandwidth 2 bits
    unsigned short s_rcv_e_;
    Uvar s_tro_;
    unsigned short s_tro_e_;   
    unsigned short s_baq_e_;//baq mode 3 bits
    unsigned short s_bem_e_;//beam mask 5 bits
    unsigned int s_at1N1_dB_;
    unsigned short s_at1N1_e_;//receiver attenuation for beam1 and 2, 12 bits
    unsigned int  s_at1N2_dB_;
    unsigned short s_at1N2_e_;//receiver attenuation for beam1 and 2, 12 bits
    unsigned int s_at1N3_dB_;
    unsigned short s_at1N3_e_;//receiver attenuation for beam1 and 2, 12 bits 
    unsigned int s_at3N1_dB_;
    unsigned short s_at3N1_e_;//receiver attenuation for beam3,12 bits
    unsigned int s_at3N2_dB_;
    unsigned short s_at3N2_e_;//receiver attenuation for beam3,12 bits
    unsigned int  s_at3N3_dB_;
    unsigned short s_at3N3_e_;//receiver attenuation for beam3,12 bits
    unsigned int  s_at4N1_dB_;
    unsigned short s_at4N1_e_;//receiver attenuation for beam4 and 5, 12 bit
    unsigned int   s_at4N2_dB_;
    unsigned short s_at4N2_e_;//receiver attenuation for beam4 and 5, 12 bits
    unsigned int  s_at4N3_dB_;
    unsigned short s_at4N3_e_;//receiver attenuation for beam4 and 5, 12 bits
    Uvar s_rip_;    //radiometer integration time, 4 bits, 5 ms LSB 10-75 ms
    unsigned short s_rip_e_;
    unsigned short s_rad_e_;//radiometer window count 8 bits range: 1-255
    Uvar  s_csd_; //chirp step duration (8 bit, 133.333 ns LSB)
    unsigned short s_csd_e_;
    unsigned short s_csq_e_;//chirp step quantity 12 bits (range: 2 -750)   
    Uvar s_cfs_;        //chirp frequency step size (16 bits, 1.788 Hz LSB)
    unsigned short s_cfs_e_;

    //---------------------------------
    //Power  variables
    //----------------------------------
    sdata powerfield_;
    Uvar p_tfi_;
    unsigned int p_tfi_e_;
    const static unsigned short p_typ_e_;
    unsigned int pmd_e_;//power mode

    //---------------------------
    //TNC field
    //----------------------------
    sdata tncfield_;
    const static unsigned short t_typ_e_;
    Uvar t_tfi_;
    unsigned int t_tfi_e_;
    unsigned int tnc_e_;
    unsigned int rl_e_;
    Uvar rl_;
    unsigned int nd_e_;
    Uvar nd_;
    unsigned short tca_;//16bit word
    unsigned short tcb_;//16bit word
    unsigned short tcc_;//16bit word
    //------------------
    //data rate
    //------------------
    double data_rate_;

    //-------------
    //calibration table  variables
    //---------------
    CAGTable cag_table_;

  };


  //-----------------------------------
  //Utility to set three independent attenuators 
  // according to attenuation map and total attenuation value
  //-----------------------------
  void computeAttenuationMap(unsigned int& at1, unsigned int& at2, unsigned int& at3,const unsigned int& total_att);


  //convert beam mask (such as 11111) to encoded value ( 31)
  void setBeamMask(unsigned int& bem, const unsigned int& beam_mask);



#endif









