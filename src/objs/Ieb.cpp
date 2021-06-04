//----------------------------------------------------------------------------
//
// Ieb.cpp
//
// This file contains method definitions for the Ieb handling classes
//----------------------------------------------------------------------------


//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include "Error.h"
#include "Units.h"
#include "Ieb.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Utils.h"


using std::cout;
using std::cerr;
using std::endl;

//-------------------------------------
//SLOW FIELD PARAMETER SETTING
//------------------------------------


//--------------------------------
//Static const memeber initialization for class Ieb
//--------------------------------
   const unsigned short Ieb::s_typ_e_ = 3;
   const unsigned short Ieb::f_typ_e_ = 2;
   const unsigned short Ieb::t_typ_e_ = 1;
   const unsigned short Ieb::p_typ_e_ = 0;


//-------------------------
// Methods for class Ieb
//-------------------------

//--------------
// Constructors: contruct with a filename and file type
//  to write/read record
//--------------
Ieb::Ieb(const Time& t) 
  :time_set_(true), 
   dutycycle_set_(false),
   f_tfi_set_(false),
   //fin_set_(false),
   pul_set_(false),
   bpd_set_(false),
   rwd_set_(false),
   pri_set_(false),
   chirp_start_frequency_set_(false),
   s_tfi_set_(false),
   dtn_set_(false),
   //sin_set_(false),
   mod_set_(false),
   csr_set_(false),
   adc_set_(false),
   rcv_set_(false),
   tro_set_(false),
   baq_set_(false),
   bem_set_(false),
   at1N1_set_(false),
   at1N2_set_(false),
   at1N3_set_(false),
   at3N1_set_(false),
   at3N2_set_(false),
   at3N3_set_(false),
   at4N1_set_(false),
   at4N2_set_(false),
   at4N3_set_(false),
   rip_set_(false),
   rad_set_(false),
   csd_set_(false),
   csq_set_(false),
   chirp_frequency_step_set_(false),
   encode_slowfast_complete_(false),
   decode_slowfast_complete_(false),
   slowfastfield_loaded_(false),
   power_loaded_(false),
   tnc_loaded_(false),
   encode_power_complete_(false), // GH: add to split Pwr/tnc
   decode_power_complete_(false), 
   encode_tnc_complete_(false),
   decode_tnc_complete_(false),
   encode_powertnc_complete_(false),
   decode_powertnc_complete_(false),
   p_tfi_set_(false),
   t_tfi_set_(false),
   tnc_set_(false),
   tca_set_(false),
   tcb_set_(false),
   tcc_set_(false),
   data_rate_set_(false),
   t_(t)   
  {
  ilxtype_ = 3;//default value
  f_fin_e_ = 1;//default value
  s_sin_e_ = 1;//default value
  pmd_e_ = 17;//default value , radar modes
  f_bii_e_ = 1;//default value, no flag for bii  
  tnc_e_ = 0;//hamilton: default TNC value  
  tca_ = 0;//hamilton: default TNC value  
  tcb_ = 0;//hamilton: default TNC value  
  tcc_ = 0;//hamilton: default TNC value 
  slowfield_.resize(Nslowfield);
  fastfield_.resize(Nfastfield);
  powerfield_.resize(Npowerfield);
  tncfield_.resize(Ntncfield);
  for (unsigned int i = 0; i <Nslowfield;++i)
    {
    slowfield_[i] = 0;
    }
  for (unsigned int i = 0; i < Nfastfield;++i)
    {
    fastfield_[i] = 0;
    }
  for (unsigned int i = 0; i < Npowerfield;++i)
    {
    powerfield_[i]=0;
    }
  for (unsigned int i = 0 ; i < Ntncfield ;++i)
    {
    tncfield_[i]=0;
    }
  s_at1N1_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N1_dB_ = 0;
  s_at4N1_dB_ = 0;


  s_at1N2_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N2_dB_ = 0;
  s_at4N2_dB_ = 0;

  s_at1N3_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N3_dB_ = 0;
  s_at4N3_dB_ = 0;
  }

Ieb::Ieb() 
  :time_set_(false),   
   dutycycle_set_(false),
   f_tfi_set_(false),
   //fin_set_(false),
   pul_set_(false),
   bpd_set_(false),
   rwd_set_(false),
   pri_set_(false),
   chirp_start_frequency_set_(false),
   s_tfi_set_(false),
   dtn_set_(false),
   //sin_set_(false),
   mod_set_(false),
   csr_set_(false),
   rcv_set_(false),
   tro_set_(false),
   baq_set_(false),
   bem_set_(false),
  at1N1_set_(false),
  at1N2_set_(false),
  at1N3_set_(false),
  at3N1_set_(false),
  at3N2_set_(false),
  at3N3_set_(false),
  at4N1_set_(false),
  at4N2_set_(false),
  at4N3_set_(false),
  rip_set_(false),
  rad_set_(false),
  csd_set_(false),
  csq_set_(false),
  chirp_frequency_step_set_(false),
  encode_slowfast_complete_(false),
  decode_slowfast_complete_(false),
  slowfastfield_loaded_(false),
  power_loaded_(false),
  tnc_loaded_(false),
  encode_power_complete_(false), // GH: add to split Pwr/tnc
  decode_power_complete_(false),
  encode_tnc_complete_(false),
  decode_tnc_complete_(false),
  encode_powertnc_complete_(false),
  decode_powertnc_complete_(false),
  p_tfi_set_(false),
  t_tfi_set_(false),
  tnc_set_(false),
  tca_set_(false),
  tcb_set_(false),
   tcc_set_(false),
   data_rate_set_(false)
  {
  ilxtype_ = 3;//default value
  f_fin_e_ = 1;//default value
  s_sin_e_ = 1;//default value
  pmd_e_ = 17;//default value , radar modes
  f_bii_e_ = 1;//default value, no flag for bii  
  tnc_e_ = 0;//hamilton: default TNC value  
  tca_ = 0;//hamilton: default TNC value  
  tcb_ = 0;//hamilton: default TNC value  
  tcc_ = 0;//hamilton: default TNC value 
  slowfield_.resize(Nslowfield);
  fastfield_.resize(Nfastfield);
  powerfield_.resize(Npowerfield);
  tncfield_.resize(Ntncfield);
  for (unsigned int i = 0; i <Nslowfield;++i)
    {
    slowfield_[i] = 0;
    }
  for (unsigned int i = 0; i < Nfastfield;++i)
    {
    fastfield_[i] = 0;
    }
  for (unsigned int i = 0; i < Npowerfield;++i)
    {
    powerfield_[i]=0;
    }
  for (unsigned int i = 0 ; i < Ntncfield ;++i)
    {
    tncfield_[i]=0;
    }
  s_at1N1_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N1_dB_ = 0;
  s_at4N1_dB_ = 0;

  s_at1N2_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N2_dB_ = 0;
  s_at4N2_dB_ = 0;

  s_at1N3_dB_ = 0;//has not been determined yet, for TA, autogain is likely
  s_at3N3_dB_ = 0;
  s_at4N3_dB_ = 0;
  }

//----------------------
//{get,set} time 
//-----------------------
Time Ieb::getTime() const
  {
  if (!time_set_) 
    {
      ErrorMessage("Ieb::getTime: no time is set").throwMe();
    }
  return(t_);
  }

void Ieb::setTime(const Time& t) 
  {
  if (time_set_)
    {
      ErrorMessage("Ieb::setTime:Time has been set already").throwMe();
    }
  t_=t;
  time_set_ = true;  
  }

//----------------------------
//{get,set} ilxtype
//4: all mode is valid: need for ideal param generation
//3: slow/fast valid
//2: only fast valid
//1: only TNC valid
//0: only power mode valid
//reset method is OK
//-----------------------------
unsigned int Ieb::getIebMod() const
  {
    if (ilxtype_ > 4) 
      ErrorMessage("getIlxtype: ilxtype is larger than 4").throwMe();
    return(ilxtype_);
  }
void Ieb::setIebMod(const unsigned int& ilxtype)
  {
  if (ilxtype > 4) 
      ErrorMessage("setIlxtype: ilxtype is larger than 4: "+toStr(ilxtype)).throwMe();
  ilxtype_ = ilxtype; 
  }

//----------------------
//{get,set} fastfield TFI
//---------------------
Uvar Ieb::getFastTfi() const
  {
    if(!f_tfi_set_)
      {
	ErrorMessage("Ieb::getFastTfi:no tfi is set").throwMe();
      }
      return(f_tfi_);
  }

void Ieb::setFastTfi(const Uvar& tfi)
  {
    if(f_tfi_set_)
      {
	ErrorMessage("Ieb::getFastTfi: tfi is already set").throwMe();
      }
    if (tfi < Uvar(0,"s") || tfi > Uvar(18.20417 * 3600.0,"s"))
      {
	ErrorMessage("Ieb::setFastTfi: tfi is out of range").throwMe();
      }
  f_tfi_ = tfi;
  f_tfi_e_ = (unsigned short) f_tfi_.getInUnits("s");
  f_tfi_ = Uvar(double(f_tfi_e_),"s"); 
  f_tfi_set_=true;
  }

//--------------------------
//{get} fastfield typ
//--------------------------
unsigned int Ieb::getFastTyp() const
  {
  return (f_typ_e_);
  }

//---------------------
//{get,set} fastfield fin
//----------------------
unsigned int Ieb::getFin() const
  {
    return (f_fin_e_);
  }
void Ieb::setFin(const unsigned int& fin)
  {
    unsigned int fin_tmp= fin;
    if (fin_tmp>= 256)
      {
	fin_tmp= fin_tmp % 256;//roll over BB13-7
      }
    f_fin_e_ = fin_tmp;  
  }

//--------------------------------
//{get,set} bursts in instruction
//--------------------------------
unsigned int Ieb::getBii() const
  {
  return(f_bii_e_);
  }
void Ieb::setBii(const unsigned int& bii)
  {
   
    if(bii >= 256 )
      { 
	ErrorMessage e("Ieb::setBii: bii should be between 1 - 255");
	e.throwMe();
      }
    if (bii==0)
      {
	cout<<"Warning: attempt to use Perpetual Instruction"<<endl;
      }
    if(!bem_set_)
      {
	ErrorMessage("Ieb::setBii: beam mask is not set").throwMe();
      }
    
  
    f_bii_e_ = bii;
  }


//------------------------
//{get,set} fastfield pul
//------------------------
unsigned int Ieb::getPul() const
  {
    if(!pul_set_)
      {
	ErrorMessage("Ieb::getPul:no pulse number is set").throwMe();
      }
    return(f_pul_e_);
  }

void Ieb::setNumberOfPulses(const unsigned int& Np)
  {
    if (pul_set_) 
      {
	ErrorMessage("Ieb::setNumberofPulse: pulse number is already set").throwMe();
      }
    if (!pri_set_)
      {
	ErrorMessage("Ieb::setNumberofPul: Can't convert pulse_on_time into integer number of pulses without pri setting").throwMe();
      }
    if (Np >= 256) 
      {
	cout << "Ieb:setNumberOfPulses = " << Np << ".   WARNING: Will change to PULmax of 254." << endl; 
	// add 11/07/2007 to catch pulse greater than 256
	//ErrorMessage("Ieb::setNumberOfPulses: pulse number is out of range").throwMe();
      }
    // add 11/07/2007 to catch pulse greater than 256
    if (Np >= 256) f_pul_e_ = 254;
	f_pul_e_ = Np;

    pul_set_ = true;

  }
Uvar Ieb::getPulseTransmitTime() 
  {
    if(!pul_set_)
      {
	ErrorMessage("Ieb::getPul:no pulse number is set").throwMe();
      }
    return(double(f_pul_e_) * getPri());
  }





//----------------------------------------------
//compute Pul, Rwd, Tro, and Csf for a given
// bpd and data rate
//Usefull when keeping a certain number of looks 
// by choosing small bpd
//----------------------------------------------
void Ieb::computeSAR_Pul_Rwd_Tro_Csf(const Uvar& bpd, const double& data_rate, 
				     const Umat& range, const Umat& doppler)
  {
    if (!pri_set_) ErrorMessage("Ieb.cpp::computeSAR_pul_rwd_tro_csf_using bpd and data rate: pri is not set" ).throwMe();
    if(!bem_set_) ErrorMessage("Ieb.cpp: beam mask is not set").throwMe();
    //set bpd
    setBpd(bpd);  
  
    //from data rate, number of words * 16/bpd = data rate
    //let's not to overshoot
    unsigned int words = (unsigned int) (data_rate * getBpd().getInUnits("s") / 16.0) -1 ;
   
    //assume 8-2 baq
    if (s_baq_e_ != 0){ ErrorMessage("Ieb.cpp::computeSAR_pul_rwd_tro_csf_using bpd and data rate: baq mode is not 8-2").throwMe();}

    //total number of data points 8 * number of words
    unsigned int data_points = 8 * words;

    // how many pulses? data_points/f_pri_e
    unsigned int const_looks_receive_window = (unsigned int) (data_points/getPointsInsidePRI());

    //compare with normal method
    //use which beam is being used
    // range(beam_number, 4 (5dB) + 1 (bore))
    vector<Uvar> range_variation,doppler_variation;
    for(unsigned int i=0;i<5;++i){
      if(bitget(s_bem_e_,i,i)==1){
	for(unsigned int j=0;j<5;++j){
	  range_variation.push_back(range(i,j));
	  doppler_variation.push_back(doppler(i,j));
	}
      }
    }
    if( (range_variation.size()!=doppler_variation.size()) && range_variation.size()==0){
      cout<<"Size of range and doppler spread containers "<< range_variation.size()<<" "<< doppler_variation.size()<<endl;
      ErrorMessage("Ieb.cpp: 5dB range/doppler spread computation error: data not collected properly").throwMe();
    }
    Uvar min_range,max_range;
    min_max(min_range,max_range,range_variation);
    Uvar min_doppler,max_doppler;
    min_max(min_doppler,max_doppler,doppler_variation);
    if (s_bem_e_ == 0){  ErrorMessage("Ieb.cpp(computeSAR_Pul_Rwd_Tro_Csf1.1): Attempt to run SAR with no (zero) beams.").throwMe();}

    //calculate roundtrip time for inner most  beam's  5dB point
    Uvar roundtrip_time = 2.0 * min_range/speed_light;
    unsigned int  roundtrip_time_in_pri 
       = (unsigned int) (roundtrip_time/getPri()).getInUnits("");
    unsigned int pul=  roundtrip_time_in_pri  - 1;
    if (const_looks_receive_window > pul) ErrorMessage("Ieb.cpp::computeSAR_Pul_Rwd_Csf: using bpd and data rate: you can have enough pulses without limiting bpd").throwMe();
     setNumberOfPulses( const_looks_receive_window);
     setRwdInUnitsOfPri(roundtrip_time_in_pri);
     setTroInUnitsOfPri(0);//zero tro
     computeChirpStartFrequency(min_doppler,max_doppler);
  }
//-----------------------------------------------------------
//compute Pul, Rwd, Tro, and Csf for SAR
// Number of pulses is chosen to fill up the roundtrip time
// for beam 1
//-----------------------------------------------------------
// Algorithm #1
//------------------------------
// PUL = 2*Rmin(1)/c/PRI - 1
// RWD = 2*Rmin(1)/c/PRI
// TRO = 2*(Rmax(5)-Rmin(1))
//------------------------------

void Ieb::computeSAR_Pul_Rwd_Tro_Csf(const Umat& range, const Umat& doppler)
 {
   unsigned int sdb_buffer_size = 0, Max_buffer_size = 0; 
   unsigned int old_pul = 0;

   //cout << "GH:computeSAR_Pul_Rwd_Tro_Csf.BAQ = " << getBaq() << "  Max Buffer = " << Max_buffer_size << endl;
   // this need to become equations based in TB.  dependant on BPD, BUF, BAQ......
   switch ( getBaq() ){ // really 32767. however make it even and give a little pad
   case (0): // 8to2 BAQ
     Max_buffer_size = 32000;
     break;
   case (1): // 8to1 BAQ
     Max_buffer_size = 32000;
     break;
   case (2): // 8to0
     Max_buffer_size = 32000;
     break;
   case (3): // PRI summation
     Max_buffer_size = 32000;
     break;
   case (4): // 8to4 MSB
     Max_buffer_size = 32000;
     break;
   case (5): // 8 straight.  Test bed start to drop around 10050  (why??)  Just at 365 kbps??
     if (data_rate_ > 240.0) {
        Max_buffer_size = 19900 ; //  test bed data discovered during TA(C43 run)
        cout << "Warning(Ieb.cpp) Buffer limited to " << Max_buffer_size << " during high date rates and straight8."  << endl;
     }
     else Max_buffer_size = 32000;
     break;
   case (6): // 8to4 BAQ
     Max_buffer_size = 32000;
     break;
   case (7): // 8to4 BAQ
     Max_buffer_size = 32000;
     break;
   default: // should never get here
     Max_buffer_size = 32000;
     cout <<"ERROR:Ieb.cpp BAQ not defined??  Using default Max Buffer Sixe limit check: " <<Max_buffer_size<< endl;
   }

   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if (doppler.rows() != 5 || doppler.cols() != 5) {
     ErrorMessage("Data type mismatch for doppler").throwMe();}
   if( (pul_set_)||(rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
	cout << "MOD: " << s_mod_e_ << endl;
     ErrorMessage("Csf: At least one of pulrwd, tro, csf is already set").throwMe();
     }

   if(!bem_set_) ErrorMessage("Ieb.cpp: beam mask is not set").throwMe();

   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
//     if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeSAR_PUL_RWD_CSF: use set method for calibration mode").throwMe();
     
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeiSAR_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

	//sar mode

       //-----------------------------------------
       //collect range spread information
       // range(beam_number, 4 (5dB) + 1 (bore))
       //----------------------------------------------
       vector<Uvar> range_variation,doppler_variation;
       for(unsigned int i=0;i<5;++i){
	 if(bitget(s_bem_e_,i,i)==1){
	   for(unsigned int j=0;j<5;++j){
	     range_variation.push_back(range(i,j));
	     doppler_variation.push_back(doppler(i,j));
	   }
	 }
       }
       if( (range_variation.size()!=doppler_variation.size()) && range_variation.size()==0){
	 cout<<"Size of range and doppler spread containers "<< range_variation.size()<<" "<< doppler_variation.size()<<endl;
	 ErrorMessage("Ieb.cpp: 5dB range/doppler spread computation error: data not collected properly").throwMe();
       }
       Uvar min_range,max_range;
       min_max(min_range,max_range,range_variation);
       Uvar min_doppler,max_doppler;
       min_max(min_doppler,max_doppler,doppler_variation);
       

    if (s_bem_e_ == 0){  ErrorMessage("Ieb.cpp(computeSAR_Pul_Rwd_Tro_Csf1.2): Attempt to run SAR with no (zero) beams.").throwMe();}

     //calculate roundtrip time for inner most  beam's  5dB point
     Uvar roundtrip_time = 2.0 * min_range/speed_light;
     unsigned int  roundtrip_time_in_pri 
       = (unsigned int) (roundtrip_time/getPri()).getInUnits("");


     setNumberOfPulses( roundtrip_time_in_pri -1);
     setRwdInUnitsOfPri(roundtrip_time_in_pri); 
     
     //cout<<"min and max range "<<min_range<<" "<<max_range<<endl;
     //cout<<"get pul "<<getPul()<<endl;
     //cout<<"rwd "<<rwd<<endl;
     //cout<<"rwd "<<getRwd()<<endl;
     //-----------------------------------------------
     //computeTro using min and max range
     //------------------------------------------------    
     Uvar  tro_time = 2.0 * max_range/speed_light - getRwd();
     
     //if (s_csr_e_ == 1) tro_time = 0;

     int tro_number = int(ceil((tro_time/getPri()).getInUnits("")));   
     setTroInUnitsOfPri(tro_number);
     //2*min_range/speed_light should not be used, instead use RWD

     
     //gh 
     sdb_buffer_size = getRadarDataPoints() ;
     // Max_buffer  = 32767/
     cout << "data in buffer = " << sdb_buffer_size << ":max: " << Max_buffer_size << endl;
     if (sdb_buffer_size > Max_buffer_size ){
       old_pul = getPul();
       pul_set_ = false;
       setNumberOfPulses ( (unsigned int)(((old_pul+tro_number) * (Max_buffer_size)/sdb_buffer_size)-tro_number) );
       sdb_buffer_size = getRadarDataPoints() ;
       cout << "Ieb.cpp:Warning: changing PUL number:NEWdata in buffer = " << sdb_buffer_size << ":max:  " << Max_buffer_size << "  PULwas:" << old_pul << "PUL:" << getPul()<< endl;
     }


     //cout<<"getTro "<<getTro()<<" "<<getTroInPriUnits()<<endl;
     //---------------------------------------
     //computeCsf
     //--------------------------------------
     computeChirpStartFrequency(min_doppler,max_doppler);
     
     //cout<<"chirp start frequency "<<getChirpStartFrequency()<<endl; 
     }
   else
     {
       ErrorMessage("Invalid radar mode trying to run computeSAR_Pul_Rwd_Tro_Csf").throwMe();
     }
 }

//-----------------------------------------------------------
//compute Pul, Rwd, Tro, and Csf for SAR
// Number of pulses is chosen to fill up the roundtrip time
// for beam 5
//-----------------------------------------------------------
// Algorithm #2
//------------------------------
// PUL = 2*Rmin(5)/c/PRI - 1
// RWD = 2*Rmin(5)/c/PRI
// TRO = 2*(Rmax(5)-Rmin(5))
//------------------------------

void Ieb::computeSAR_Pul_Rwd_Tro_Csf2(const Umat& range, const Umat& doppler)
 {
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if (doppler.rows() != 5 || doppler.cols() != 5) {
     ErrorMessage("Data type mismatch for doppler").throwMe();}
   if( (pul_set_)||(rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
     ErrorMessage("Csf2:At least one of pulrwd, tro, csf is already set").throwMe();
     }

   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
     //if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeSAR_PUL_RWD_CSF: use set method for calibration mode").throwMe();
     
     if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeiSAR_PUL_RWD_CSF: incorrect setting of CSR").throwMe();


     for(unsigned int i=0;i<5;++i){
       if(bitget(s_bem_e_,i,i)==0)
	 ErrorMessage("This method requres all 5 beam use").throwMe();
     }
	   
     //sar mode 
     //-----find beam 1 (near) and beam 5 (far)-----
     unsigned int beam_1, beam_5;
     if(range(4,4) > range(0,4))
       {
         beam_5 = 4; beam_1 = 0;
       }
     else
       { 
	 beam_1 = 4; beam_5 = 0;
       } 
     Uvec range_beam_5("range_beam_5",5);
     for(unsigned int i = 0 ; i < 5 ; i++)
       range_beam_5(i) = range(beam_5,i);

     Uvar min_range,max_range;
     min_max(min_range,max_range,range_beam_5);
     //---------------------------------------------
     Uvar min_doppler,max_doppler;
     min_max(min_doppler,max_doppler,doppler);

     if (s_bem_e_ == 0){  ErrorMessage("Ieb.cpp(computeSAR_Pul_Rwd_Tro_Csf2): Attempt to run SAR with no (zero) beams.").throwMe();}

     //calculate roundtrip time for inner most  beam's  5dB point
     Uvar roundtrip_time = 2.0 * min_range/speed_light;
     unsigned int  roundtrip_time_in_pri 
       = (unsigned int) (roundtrip_time/getPri()).getInUnits("");
     setNumberOfPulses( roundtrip_time_in_pri -1);
     setRwdInUnitsOfPri(roundtrip_time_in_pri);

  
     //cout<<"min and max range "<<min_range<<" "<<max_range<<endl;
     //cout<<"get pul "<<getPul()<<endl;
     //cout<<"rwd "<<rwd<<endl;
     //cout<<"rwd "<<getRwd()<<endl;
     //-----------------------------------------------
     //computeTro using min and max range
     //------------------------------------------------    
     Uvar  tro_time = 2.0 * max_range/speed_light - getRwd();
     //2*min_range/speed_light should not be used, instead use RWD

     int tro_number = int(ceil((tro_time/getPri()).getInUnits("")));   
     setTroInUnitsOfPri(tro_number);
 
     //cout<<"getTro "<<getTro()<<" "<<getTroInPriUnits()<<endl;
     //---------------------------------------
     //computeCsf
     //--------------------------------------
     computeChirpStartFrequency(min_doppler,max_doppler);
     
     //cout<<"chirp start frequency "<<getChirpStartFrequency()<<endl; 
     }
   else
     {
       ErrorMessage("Invalid radar mode trying to run computeSAR_Pul_Rwd_Tro_Csf").throwMe();
     }
 }

//-----------------------------------------------------------
//compute Pul, Rwd, Tro, and Csf for SAR
// Number of pulses is chosen to fill up the roundtrip time
// for beam 5
//-----------------------------------------------------------
// Algorithm #3
//------------------------------
// PUL = 2*Rmin(5)/c/PRI - 1
// RWD = 2*Rmin(5)/c/PRI
// TRO = 2*(Rmax(1)-Rmin(5))
//------------------------------

void Ieb::computeSAR_Pul_Rwd_Tro_Csf3(const Umat& range, const Umat& doppler)
 {
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if (doppler.rows() != 5 || doppler.cols() != 5) {
     ErrorMessage("Data type mismatch for doppler").throwMe();}
   if( (pul_set_)||(rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
     ErrorMessage("Csf3:At least one of pulrwd, tro, csf is already set").throwMe();
     }

   for(unsigned int i=0;i<5;++i){
     if(bitget(s_bem_e_,i,i)==0)
       ErrorMessage("This method requres all 5 beam use").throwMe();
   }

   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
    // if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeSAR_PUL_RWD_CSF: use set method for calibration mode").throwMe();
     
     if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeiSAR_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

     //sar mode 
     //-----find beam 1 (near) and beam 5 (far)-----
     unsigned int beam_1, beam_5;
     if(range(0,4) > range(0,4))
       {
         beam_5 = 4; beam_1 = 0;
       }
     else
       { 
	 beam_1 = 4; beam_5 = 0;
       } 
     Uvec range_beam_1("range_beam_1",5);
     Uvec range_beam_5("range_beam_5",5);
     for(unsigned int i = 0 ; i < 5 ; i++)
       {
         range_beam_1(i) = range(beam_1,i);
         range_beam_5(i) = range(beam_5,i);
       }

     Uvar min_range_1,max_range_1;
     Uvar min_range_5,max_range_5;
     min_max(min_range_1,max_range_1,range_beam_1);
     min_max(min_range_5,max_range_5,range_beam_5);

     Uvar min_range,max_range;
     min_range = min_range_5;
     max_range = max_range_1;
     //---------------------------------------------
     Uvar min_doppler,max_doppler;
     min_max(min_doppler,max_doppler,doppler);

     if (s_bem_e_ == 0){  ErrorMessage("Ieb.cpp(computeSAR_Pul_Rwd_Tro_Csf3): Attempt to run SAR with no (zero) beams.").throwMe();}

     //calculate roundtrip time for nearest beam's 3dB or 5dB point
     Uvar roundtrip_time = 2.0 * min_range/speed_light;
     unsigned int  roundtrip_time_in_pri 
       = (unsigned int) (roundtrip_time/getPri()).getInUnits("");
     setNumberOfPulses( roundtrip_time_in_pri -1);
     setRwdInUnitsOfPri(roundtrip_time_in_pri);
  
     //cout<<"min and max range "<<min_range<<" "<<max_range<<endl;
     //cout<<"get pul "<<getPul()<<endl;
     //cout<<"rwd "<<rwd<<endl;
     //cout<<"rwd "<<getRwd()<<endl;
     //-----------------------------------------------
     //computeTro using min and max range
     //------------------------------------------------    
     Uvar  tro_time = 2.0 * max_range/speed_light - getRwd();
     //2*min_range/speed_light should not be used, instead use RWD

     int tro_number = int(ceil((tro_time/getPri()).getInUnits("")));   
     setTroInUnitsOfPri(tro_number);

     //cout<<"getTro "<<getTro()<<" "<<getTroInPriUnits()<<endl;
     //---------------------------------------
     //computeCsf
     //--------------------------------------
     computeChirpStartFrequency(min_doppler,max_doppler);
     
     //cout<<"chirp start frequency "<<getChirpStartFrequency()<<endl; 
     }
   else
     {
       ErrorMessage("Invalid radar mode trying to run computeSAR_Pul_Rwd_Tro_Csf").throwMe();
     }
 }


//------------------------
//Set Rwd, Tro, and Csf
// Number of pulses do not need to fill up the roundtrip time
// also Beam 3 only
//------------------------
void Ieb::computeALT_Rwd_Tro_Csf(const Umat& range, const Umat& doppler)
 {
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if (doppler.rows() != 5 || doppler.cols() != 5) {
     ErrorMessage("Data type mismatch for doppler").throwMe();}
   if(!pul_set_) ErrorMessage("computeRwd...-No pul set").throwMe();
   if((rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
       ErrorMessage("At least one of rwd, tro, csf is already set").throwMe();
     }
   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
       ErrorMessage("Sar radar mode but try to computeALT_RWD_Tro_Csf").throwMe();
     }
	 
   if ( s_mod_e_ == 1  ||s_mod_e_==9)
     {//ALTH mode
       // GH090104  if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

      //------------------
       //ALTH normal operation
       //-----------------
       //if (s_bem_e_ != 4) ErrorMessage("Attempt to run high alt using other than beam 3").throwMe(); 
       if (s_bem_e_ != 4) cout << "Ieb.cpp:Warning:Attempt to run altimeter using other than beam 3."  << endl;

       //set RWD
       if (s_csr_e_ !=1){
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       //assuming number of pulses 21
       //tro -6
       //receiver window 15
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri+3);
       setTroInUnitsOfPri(-6);

       //gh110624 set min on max if not defined in config file
	//POLICY on PUL during altimnetry:
	// if set by config.  check for roundtrip time and maxvalue
	// if set to 0.  set to roundtrip time and maxvalue
	// therefore connot set to 0... 
	//tro in rounttrip calculation.  so possible early return while transmitting
        unsigned int RoundtripPul = ( (int) (roundtrip_time_in_pri-getTroInPriUnits()/2) -1 );
        //unsigned int RoundtripPul = roundtrip_time_in_pri-1; //without TRO
//output some test messages: 20110726
	//unsigned int MaxPul = 21;
// test messages calculatenumber of pulses to fill buffer
	unsigned int maxaltbuffer = 30000; //kbytes
        //20110726 pul calculation added gh
	unsigned int buffFullPul = (unsigned int)( (maxaltbuffer / (f_pri_e_*10) ) - getTroInPriUnits() ) ;
        unsigned int MaxPul = buffFullPul ; // was 21...

//cout << "20110726:  buffer calc " << bufcalc << endl;


// ner code 20110726 gh


        // set by config however too large
        if (getPul() > RoundtripPul){ 
                cout << "Pulses violation roundtrip time.  Pul was:" << getPul() << " now: " << RoundtripPul << endl;
                pul_set_ = false;
                setNumberOfPulses ( RoundtripPul ) ;
        }
        //set to roundtrip time
        if (getPul() == 0){ //set to roundtrip time
                pul_set_ = false;
                setNumberOfPulses ( RoundtripPul ) ;
        }
        // insure roundtrip or use setting does not go over maximum allowed
        if ( getPul() > MaxPul ) {
                cout << "Warning:  Over ALTH max pulses. Was:" << getPul() << " now: " << MaxPul << endl;
		pul_set_ = false;
                setNumberOfPulses ( MaxPul );
        }
        // gh110624
cout << "20110726:  max values are (roundtrip, buffer, using) : " << RoundtripPul << ":"<< buffFullPul <<" = " << getPul()<< endl;

       //check rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip(PRI):" << roundtrip_time_in_pri << endl;
	   cout << "rwd(PRI):" << getRwdInPriUnits() << endl;
	   cout << "pul(PRI):" << getPul() << endl;
           cout << "rwd(s):" << rwd << endl;
	   cout << "tro(s):" << getTro() << endl;
	   cout << "rwdpul(PRI):" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTH): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //------------------------------
       //check data size
       //should be smaller than 32k
       //----------------------------------
       int receiver_window_in_pri_units= int(getPul())+getTroInPriUnits();
       unsigned int data_points = f_pri_e_*10 ;
       data_points *=(unsigned int) receiver_window_in_pri_units;
       if (data_points >=  Nmax_buffer)
	 {
	   ErrorMessage("computeALT_RWD_TRO_CSF: number of collected data points exceeds 32K bytes ").throwMe();
	 }
       
       computeChirpStartFrequency(doppler(2,4));//boresight doppler
     }
   
   else if (s_mod_e_ ==0 || s_mod_e_ == 8)
     {//ALTL (SCAT) mode normal operations
       // GH090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

  	if (s_bem_e_ != 4) ErrorMessage("Attempt to run low alt using other than beam3").throwMe();

       //set rwd
       if (s_csr_e_ !=1) {
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       // close it 6 pri later
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri-3);
       setTroInUnitsOfPri(6);

       // handle multiply burst_in_flight here....
      // cout << "GHbedub3: roundtrip(inpri): " <<roundtrip_time_in_pri<< " RWD(sec): " << rwd << " :RWD(pri)):"<< f_rwd_e_ << " BPD: " << getBpd()<<f_bpd_ <<  "nbursts: " << bursts_in_flight << endl;
       

       //check whether rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTL1): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }       
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //csf using beam3 boresight doppler
       computeChirpStartFrequency(doppler(2,4));//beam3's boresight only
     }
   else
     {
       ErrorMessage("No valid radar mode for computeALT_PUL_RWD_CSF").throwMe();
     }
 }
void Ieb::computeALT_Rwd_Tro_Csf(const Umat& range, const Umat& doppler, const int& tro)
{
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if (doppler.rows() != 5 || doppler.cols() != 5) {
     ErrorMessage("Data type mismatch for doppler").throwMe();}
   if(!pul_set_) ErrorMessage("computeRwd...-No pul set").throwMe();
   if((rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
       ErrorMessage("At least one of rwd, tro, csf is already set").throwMe();
     }
   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
       ErrorMessage("Sar radar mode but try to computeALT_RWD_Tro_Csf").throwMe();
     }
	 
   if ( s_mod_e_ == 1  ||s_mod_e_==9)
     {//ALTH mode
       // GH090104  if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

      //------------------
       //ALTH normal operation
       //-----------------
       //if (s_bem_e_ != 4) ErrorMessage("Attempt to run high alt using other than beam 3").throwMe(); 
       if (s_bem_e_ != 4) cout << "Ieb.cpp:Warning:Attempt to run altimeter using other than beam 3."  << endl;

       //set RWD
       if (s_csr_e_ !=1){
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       
       //------------------------------------
       //To increase valid time, we will open 
       // the receiver window -3 pri earlier and 
       //assuming number of pulses 21
       //tro -6
       //receiver window 15
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri+3);
       setTroInUnitsOfPri(tro);

       //check rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTH): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //------------------------------
       //check data size
       //should be smaller than 32k
       //----------------------------------
       int receiver_window_in_pri_units= int(getPul())+getTroInPriUnits();
       unsigned int data_points = f_pri_e_*10 ;
       data_points *=(unsigned int) receiver_window_in_pri_units;
       if (data_points >=  Nmax_buffer)
	 {
	   ErrorMessage("computeALT_RWD_TRO_CSF: number of collected data points exceeds 32K bytes ").throwMe();
	 }
       
       computeChirpStartFrequency(doppler(2,4));//boresight doppler
     }
   
   else if (s_mod_e_ ==0 || s_mod_e_ == 8)
     {//ALTL (SCAT) mode normal operations
       // GH090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

  	if (s_bem_e_ != 4) ErrorMessage("Attempt to run low alt using other than beam3").throwMe();

       //set rwd
       if (s_csr_e_ !=1) {
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       // close it 6 pri later
       //------------------------------------
       if( tro%2 !=0) ErrorMessage("Ieb.cpp:: scatt tro is not an even number").throwMe();
       /* Removed by RW 2017-08-16 to allow zero TRO in sa_290
       if (baq_set_ && s_baq_e_ == 5 && tro==0)
         {  // in normal 8-8 scatt, we expect a positive TRO
         ErrorMessage("Ieb.cpp:: TRO is zero").throwMe();
         }
       */
       if(tro<0) ErrorMessage("Ieb.cpp:: TRO is negative").throwMe();

       setRwdInUnitsOfPri(roundtrip_time_in_pri-(tro/2));
       setTroInUnitsOfPri(tro);

       // handle multiply burst_in_flight here....
      // cout << "GHbedub3: roundtrip(inpri): " <<roundtrip_time_in_pri<< " RWD(sec): " << rwd << " :RWD(pri)):"<< f_rwd_e_ << " BPD: " << getBpd()<<f_bpd_ <<  "nbursts: " << bursts_in_flight << endl;
       

       //check whether rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTL2): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }       
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //csf using beam3 boresight doppler
       computeChirpStartFrequency(doppler(2,4));//beam3's boresight only
     }
   else
     {
       ErrorMessage("No valid radar mode for computeALT_PUL_RWD_CSF").throwMe();
      }
}
//---------------------------
//speical alt waveform with a given chirp start freqeuncy
//--------------------------
void Ieb::computeToneALT_Rwd_Tro(const Umat& range, const Uvar& chirp_start_freq)
 {
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if(!pul_set_) ErrorMessage("computeRwd...-No pul set").throwMe();
   if((rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
       ErrorMessage("At least one of rwd, tro, csf is already set").throwMe();
     }
   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
       ErrorMessage("Sar radar mode but try to computeALT_RWD_Tro_Csf").throwMe();
     }
   
   if ( s_mod_e_ == 1  ||s_mod_e_==9)
     {//ALTH mode
       //GH 090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

	//------------------
       //ALTH normal operation
       //-----------------
       //if (s_bem_e_ != 4) ErrorMessage("Attempt to run high alt using other than beam 3").throwMe(); 
       if (s_bem_e_ != 4) cout << "Ieb.cpp:Warning:Attempt to run altimeter using other than beam 3."  << endl;

       //set RWD
       if (s_csr_e_ !=1) {
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       //assuming number of pulses 21
       //tro -6
       //receiver window 15
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri+3);
       setTroInUnitsOfPri(-6);

       //check rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTH): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       
       //------------------------------
       //check data size
       //should be smaller than 32k
       //----------------------------------
       int receiver_window_in_pri_units= int(getPul())+getTroInPriUnits();
       unsigned int data_points = f_pri_e_*10 ;
       data_points *=(unsigned int) receiver_window_in_pri_units;
       if (data_points >=  Nmax_buffer)
	 {
	   ErrorMessage("computeALT_RWD_TRO_CSF: number of collected data points exceeds 32K bytes ").throwMe();
	 }
        setChirpStartFrequency(chirp_start_freq);
     }
   
   else if (s_mod_e_ ==0 || s_mod_e_ == 8)
     {//ALTL (SCAT) mode normal operations
       //GH 090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
	
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();
	
       if (s_bem_e_ != 4) ErrorMessage("Attempt to run low alt using other than beam 3").throwMe();
       
       //set rwd
       if (s_csr_e_ != 1){
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       // close it 6 pri later
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri-3);
       setTroInUnitsOfPri(6);

       // handle multiply burst_in_flight here....
      // cout << "GHbedub3: roundtrip(inpri): " <<roundtrip_time_in_pri<< " RWD(sec): " << rwd << " :RWD(pri)):"<< f_rwd_e_ << " BPD: " << getBpd()<<f_bpd_ <<  "nbursts: " << bursts_in_flight << endl;
       

       //check whether rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTL3): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }       
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //set chirp start freq
       setChirpStartFrequency(chirp_start_freq);
     }
   else
     {
       ErrorMessage("No valid radar mode for computeALT_PUL_RWD_CSF").throwMe();
     }
 }


//------------------------
//Manual TRO setting
//------------------------
 void Ieb::computeToneALT_Rwd_Tro(const Umat& range, const Uvar& chirp_start_freq, const int& tro)
{
   if(!mod_set_) ErrorMessage("Radar mode is not set").throwMe();
   if(!pri_set_) ErrorMessage("Pri is not set").throwMe();
   if (range.rows() != 5 || range.cols() != 5) {
     ErrorMessage("Data type mismatch for range ").throwMe();}
   if(!pul_set_) ErrorMessage("computeRwd...-No pul set").throwMe();
   if((rwd_set_) ||(tro_set_)||(chirp_start_frequency_set_)) 
     {
       ErrorMessage("At least one of rwd, tro, csf is already set").throwMe();
     }
   if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
     {
       ErrorMessage("Sar radar mode but try to computeALT_RWD_Tro_Csf").throwMe();
     }
   
   if ( s_mod_e_ == 1  ||s_mod_e_==9)
     {//ALTH mode
       //GH 090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
       
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();

	//------------------
       //ALTH normal operation
       //-----------------
       //if (s_bem_e_ != 4) ErrorMessage("Attempt to run high alt using other than beam 3").throwMe(); 
       if (s_bem_e_ != 4) cout << "Ieb.cpp:Warning:Attempt to run altimeter using other than beam 3."  << endl;

       //set RWD
       if (s_csr_e_ !=1) {
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       //assuming number of pulses 21
       //tro -6
       //receiver window 15
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri+3);
       setTroInUnitsOfPri(tro);

       //check rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTH): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       
       //------------------------------
       //check data size
       //should be smaller than 32k
       //----------------------------------
       int receiver_window_in_pri_units= int(getPul())+getTroInPriUnits();
       unsigned int data_points = f_pri_e_*10 ;
       data_points *=(unsigned int) receiver_window_in_pri_units;
       if (data_points >=  Nmax_buffer)
	 {
	   ErrorMessage("computeALT_RWD_TRO_CSF: number of collected data points exceeds 32K bytes ").throwMe();
	 }
        setChirpStartFrequency(chirp_start_freq);
     }
   
   else if (s_mod_e_ ==0 || s_mod_e_ == 8)
     {//ALTL (SCAT) mode normal operations
       //GH 090104: if (s_csr_e_ >0  && s_csr_e_ <8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: use set method for calibration mode").throwMe();
	
       if ( s_csr_e_ == 6 || s_csr_e_ >8) ErrorMessage("Ieb.cpp::computeALT_PUL_RWD_CSF: incorrect setting of CSR").throwMe();
	
       if (s_bem_e_ != 4) ErrorMessage("Attempt to run low alt using other than beam 3").throwMe();
       
       //set rwd
       if (s_csr_e_ != 1){
       Uvar rwd =2.0 * range(2,4)/speed_light;//boresight range
       unsigned int  roundtrip_time_in_pri 
	 = (unsigned int) (rwd/getPri()).getInUnits("");
       //------------------------------------
       //To crease valid time, we will open 
       // the receiver window -3 pri earlier and 
       // close it 6 pri later
       //------------------------------------
       setRwdInUnitsOfPri(roundtrip_time_in_pri-3);
       setTroInUnitsOfPri(tro);

       // handle multiply burst_in_flight here....
      // cout << "GHbedub3: roundtrip(inpri): " <<roundtrip_time_in_pri<< " RWD(sec): " << rwd << " :RWD(pri)):"<< f_rwd_e_ << " BPD: " << getBpd()<<f_bpd_ <<  "nbursts: " << bursts_in_flight << endl;
       

       //check whether rwd >= pul+1
       int Rwd_Pul = int(getRwdInPriUnits())-int(getPul());
       if (Rwd_Pul < 1 )
	 {
	   cout << "ERROR....." << endl ;
	   cout << "roundtrip:" << roundtrip_time_in_pri << endl;
	   cout << "rwd:" << rwd << endl;
	   cout << "pul:" << getPul() << endl;
	   cout << "tro:" << getTro() << endl;
	   cout << "rwdpul:" << Rwd_Pul << endl;
	   ErrorMessage("Ieb.cpp::computeALT_Rwd_Tro_Csf(ALTL4): RWD is need to be at least 1 PRI larger than Number of Pulses transmitted ").throwMe();
	 }       
       } else { //therefore csr = 1
           setRwdInUnitsOfPri(0);
           setTroInUnitsOfPri(0);
       }

       //set chirp start freq
       setChirpStartFrequency(chirp_start_freq);
     }
   else
     {
       ErrorMessage("No valid radar mode for computeALT_PUL_RWD_CSF").throwMe();
ErrorMessage("not done yet").throwMe();
      }
}
//------------------------
//{get,set} burst period
//-----------------------
Uvar Ieb::getBpd() const
  {
    if (!bpd_set_)
      {
	ErrorMessage("Ieb::getBdp:no burst period is set");
      }
    return(f_bpd_);
  }

unsigned int Ieb::getBpdEncodedValue() const
  {
    if (!bpd_set_)
      {
	ErrorMessage("Ieb::getBdpEncodedValue():no burst period is set");
      }
    return(f_bpd_e_);
  }


void Ieb::setBpd(const Uvar&  bpd)
  {  
  if(bpd_set_) ErrorMessage("Ieb.cpp::setBpd: bpd is set").throwMe();
  if (bpd > Uvar(4095*mstos,"s")|| bpd < Uvar(10*mstos,"s"))
    {
      cout << "BPD = " <<  bpd << endl;;
      ErrorMessage("Ieb::set Bdp: burst period is out of range").throwMe();
    }
     
  f_bpd_ = bpd;
  f_bpd_e_= (unsigned short)ceil((f_bpd_.getInUnits("s") * 1000.0));
  f_bpd_ = double(f_bpd_e_) * Uvar(1.0*mstos,"s");
  bpd_set_ = true;
  }


void Ieb::increaseBpd(const Uvar& increase_bpd)
  {  
    // will increase bpd (flag for original setting should be true)
  if(bpd_set_) 
    {
    if (increase_bpd > Uvar(4095*mstos,"s")|| increase_bpd < Uvar(10*mstos,"s"))
      {
	cout << "BPD = " <<  increase_bpd.getInUnits("ms") << endl;
	ErrorMessage("Ieb::increase Bpd: burst period is out of range").throwMe();
      }

    if (f_bpd_ > increase_bpd)
      {
	cout << "Current Bpd = " << f_bpd_ << " Increase to = "<< increase_bpd << endl;
	ErrorMessage("Ieb::increaseBpd: Trying to increase BPD, but commanded a smaller value").throwMe();
      }

     
    f_bpd_ = increase_bpd;
    f_bpd_e_= (unsigned short)ceil((f_bpd_.getInUnits("s") * 1000.0));
    f_bpd_ = double(f_bpd_e_) * Uvar(1.0*mstos,"s");
    bpd_set_ = true;
    } else
      ErrorMessage("Ieb.cpp::increaseBpd: can increase if original value is not set").throwMe();
      
  }

//----------------------------------
//{get,set} receive window delay
//-----------------------------------
Uvar Ieb::getRwd()
  {
    if(!rwd_set_)
      {
	ErrorMessage("Ieb::getRwd:rwd is not set");
      }
  return(double(f_rwd_e_)*getPri());
  }

unsigned int Ieb::getRwdInPriUnits() const
  {
    if(!rwd_set_)
      {
	ErrorMessage("Ieb::getRwdEncodedValue():rwd is not set");
      }
  return(f_rwd_e_);
  }

//------------------------------
//set RWD in units of PRI
//------------------------------
void Ieb::setRwdInUnitsOfPri(const unsigned int& rwd)
  {
    if (!pri_set_)
      {
	ErrorMessage("Ieb::setRwd: pri is not set to convert rwd in units of pri");
      }
    if(rwd_set_) 
      {
	cout << "Ieb.cpp:setRwd: mod = " << s_mod_e_ <<":"<<s_csr_e_<< endl;
	ErrorMessage("Ieb.cpp::setRwdInUnitsOfPri:rwd is set").throwMe();
      }
    
     if(s_csr_e_ == 1) {
    	f_rwd_e_ =  0 ;
    	f_rwd_ = double(f_rwd_e_) * getPri();
    	rwd_set_ = true;
    } else {
        f_rwd_e_ =  rwd ;
        f_rwd_ = double(f_rwd_e_) * getPri();
        rwd_set_ = true;
    } 


  }

void Ieb::setRwd(const Uvar& rwd)
  {
    if(rwd_set_)
      {
	ErrorMessage("Ieb::setRwd: rwd is already set").throwMe();
      }
    if (!pri_set_)
      {
	ErrorMessage("Ieb::setRwd: pri is not set to convert rwd in units of pri");
      }
    if (rwd < Uvar(0,"s"))
      {
	ErrorMessage e("Ieb::setRwd: receiver delay window should be larger than 0 sec");
	e.throwMe();
      }

    // gh change 12/04/2003 this should round
    unsigned int rwd_tmp = (unsigned int) (round_double((rwd/f_pri_).getInUnits(""))); // this will truncate
    //unsigned int rwd_tmp = (unsigned int) (rwd/f_pri_).getInUnits(""); // this will truncate
    if (rwd_tmp> 1023) 
      {
	cout << "rwd(in pri) = rwd_tmp"<<":"<<" Bpd = "<< f_bpd_e_ << endl;  
	ErrorMessage("Ieb::setRwd: rwd time is too long").throwMe();
      }
    f_rwd_e_ = rwd_tmp;
    f_rwd_ = f_pri_ * double(f_rwd_e_);//reset rwd in units of pri
    rwd_set_= true;
  }



//--------------------
//{get,set,compute} PRI and dutycycle
//---------------------
Uvar Ieb::getPri() const
  {
  if (!pri_set_)
    {
      ErrorMessage("Ieb::getPri: no pri is set");
    }
  return(f_pri_);
  }

unsigned int Ieb::getPriEncodedValue() const
  {
    if(!pri_set_)
      ErrorMessage("Ieb::getPriEncodedValue: no pri is set");
    return(f_pri_e_);
  }

void Ieb::setPri(const Uvar& pri)
  {
  if(pri_set_)
    {
      ErrorMessage("Ieb::setPri: pri is already set");
    }
  if(!adc_set_)
    {
      ErrorMessage("Ieb:setPri:No adc clock speed is set to compute an encoded value for pri").throwMe();
    }
  if (pri <=Uvar(0,"s"))
    {
    ErrorMessage("Ieb::setPri: Pri should be larger than 0 sec").throwMe();
    }

  if (s_adc_ == sarh_adc || s_adc_==sarl_adc || s_adc_ ==altl_adc)
    { //12/04/2003 GH change
      unsigned int pri_number =(unsigned int)(round_double((s_adc_*pri).getInUnits("")));// this should rounds
    //unsigned int pri_number =(unsigned int)((s_adc_*pri).getInUnits("")); // this will truncate
    pri_number = pri_number/2;
    pri_number = pri_number *2;//need to be even
    f_pri_ = double(pri_number) /s_adc_;
    prf_ = 1.0/f_pri_;
    f_pri_e_ = pri_number;
    }
  else if (s_adc_ == alth_adc)
    {//12/04/2003 GH change
      unsigned int pri_number =(unsigned int)(round_double((s_adc_*pri/10.0).getInUnits("")));// this should rounds
      //unsigned int pri_number = (unsigned int)((s_adc_*pri/10.0).getInUnits(""));
    f_pri_ = double(pri_number)*10.0/s_adc_;
    prf_ = 1.0/f_pri_;
    f_pri_e_ = pri_number;
    }
  else
    {
    ErrorMessage("No valid adc clock speed: see BB4-68").throwMe();
    }
  if (f_pri_e_ < 2 || f_pri_e_ > 1023)
    {
	cout << " PRI value out of range (within 2-1023): "<< f_pri_e_ << endl;
      ErrorMessage("Ieb::setPri:pri value in units of adc clock speed is out of range (see BB 13-16").throwMe();
    }
  pri_set_ = true;
  }

void Ieb::computePri(const Uvar& altitude,const Dvec& prf_alt_fit)
  {    
    if (prf_alt_fit.size()==0) ErrorMessage("no fitting coefficients exist").throwMe();
    double alt = altitude.getInUnits("km")/1000.0;
    double prf_in_KHz = 0;
    if(s_mod_e_ ==2 || s_mod_e_ ==3 || s_mod_e_ == 10 || s_mod_e_ ==11)
      {
	//sar mode
	if (s_bem_e_ !=31) ErrorMessage("Attempt to run SAR not using all 5 beams, may need to update computePul_Rwd_Tro_Csf").throwMe();
	for (unsigned int i = 0; i < prf_alt_fit.size();++i)
	  {
	    prf_in_KHz += prf_alt_fit(i) * pow(alt,i);
	  }
	Uvar prf= Uvar(prf_in_KHz * 1000.0,"Hz");
	setPri(1.0/prf);
      }
    else
      {
	ErrorMessage("Invalid radar mode trying to run computePri using prf vs alt relationship").throwMe();
      }
  }

//---------------------------
//{get,set} duty cycle
//--------------------------
void Ieb::setDutycycle(const double& dutycycle) 
  {
  dutycycle_= dutycycle;
  if (dutycycle > 75.0/100.0)
    {
      ErrorMessage("Ieb::setDutycycle: duty cycle is out of range > 0.75").throwMe();
    }
  dutycycle_set_ = true;
  }

double Ieb::getDutycycle() const
  {
    if(!dutycycle_set_) ErrorMessage("No dutycycle is set").throwMe();
    return(dutycycle_);
  }

void Ieb::computeDutycycle(const Uvar& altitude, const Dvec& dutycycle_alt_fit)
  {
    if (dutycycle_alt_fit.size()==0) ErrorMessage("no fitting coefficients exist").throwMe();

    double alt_1000km = altitude.getInUnits("km")/1000.0;
        
    double dutycycle = 0;
    for (unsigned int i = 0; i < dutycycle_alt_fit.size();++i)
      {
	dutycycle+= dutycycle_alt_fit(i) * pow(alt_1000km,i);
      }
    setDutycycle(dutycycle);
  }
//---------------------------------------
//{get,set, compute} csf: chirp start frequency
//----------------------------------------
Uvar Ieb::getChirpStartFrequency() const
  {
    if(!chirp_start_frequency_set_)
      {
	ErrorMessage("Ieb::getChirpStartFrequency: no chirp start frequency is set").throwMe();
      }
    return (f_csf_);
  }

unsigned int Ieb::getChirpStartFrequencyEncodedValue() const
  {
    if(!chirp_start_frequency_set_)
      {
	ErrorMessage("Ieb::getChirpStartFrequencyEncodedValue: no chirp start frequency is set").throwMe();
      }
    return (f_csf_e_);
  }

//----------------------
//Force to use a specific chirp start frequency
//---------------------
void Ieb::setChirpStartFrequency(const Uvar& csf)
  {
  if(chirp_start_frequency_set_)
    {
    ErrorMessage("Ieb::setChirpStartFrequency: chirp start frequency is already set").throwMe();
    }
  if (csf <Uvar(457.764,"Hz"))
    {
    ErrorMessage e("Ieb::setChirpStartFrequency: Chirp start frequency should be larger than 457.764 Hz");
    e.throwMe();
    }
  if (csf > Uvar(29.99954 * 1e6,"Hz"))
    {
    ErrorMessage("Ieb::setChirpStartFrequency: csf value is too large").throwMe();
    }
  if (!mod_set_)
    {
    ErrorMessage("Ieb::setChirpStartFrequency: radar mode is not set").throwMe();
    }
  
  f_csf_e_ = (unsigned short)round_double( (csf/(Uvar(30*1e6,"Hz")/pow(2,16))).getInUnits(""));
  f_csf_ = Uvar(30*1e6,"Hz")/pow(2,16)*double(f_csf_e_);

  if (getAdc() == sarh_adc)
    {
      if (f_csf_ < Uvar(8.7539e6,"Hz") || f_csf_ > Uvar(9.3916e6,"Hz")) {
	cout << "Ieb.cpp:Warning:setChirpstartfrequency(sarh): starting chirp is out of HLD range (8.7539 to 9.3916): CSF is:" <<f_csf_  << endl;
	//	ErrorMessage("Ieb.cpp::setChirpstartfrequency(sarh): starting chirp is out of range").throwMe();
      }
    }
  else if (getAdc() == sarl_adc)
    {
      if (f_csf_ < Uvar(9.0631,"Hz") || f_csf_ > Uvar(10.0119e6,"Hz")) {
	cout << "Ieb.cpp:Warning:setChirpstartfrequency(sarl): starting chirp is out of HLD range (9.0631 to 10.0119): CSF is:" <<f_csf_ << endl;
	//	ErrorMessage("Ieb.cpp::setChirpstartfrequency(sarl): starting chirp is out of range").throwMe();
      }
    }
  else if (getAdc() == alth_adc)
    {
      if (f_csf_ < Uvar(4.8471e6,"Hz") || f_csf_ > Uvar(5.9029e6,"Hz")) {
	cout << "Ieb.cpp:Warning:setChirpstartfrequency(alth): starting chirp is out of HLD range (4.8471 to 5.9029): CSF is:" <<f_csf_  << endl;
	//	ErrorMessage("Ieb.cpp::setChirpstartfrequency(alth): starting chirp is out of range").throwMe();
      }
    }
  else if (getAdc() == altl_adc)
    {
      if (f_csf_ < Uvar(9.3375e6,"Hz") || f_csf_ > Uvar(10.4312e6,"Hz")) {
	cout << "Ieb.cpp:Warning:setChirpstartfrequency(altl): starting chirp is out of HLD range (9.3375 to 10.4312): CSF is:" <<f_csf_  << endl;
	//	ErrorMessage("Ieb.cpp::setChirpstartfrequency(altl): starting chirp is out of range").throwMe();
      }
    }
  else
    {
      ErrorMessage("Ieb.cpp::setChirpstartfrequency: no adc set").throwMe();
    }
  chirp_start_frequency_set_ = true;
  }

//-----------------------------
//Use, most likely beam3, doppler to compute chirt start frequency
// For reference, see memo 334YZ-2003-002
//-----------------------------
void Ieb::computeChirpStartFrequency(const Uvar& doppler)
  {
  if(chirp_start_frequency_set_)
    {
    ErrorMessage("Ieb::setChirpStartFrequency: chirp start frequency is already set").throwMe();
    }
  if (!mod_set_)
    {
    ErrorMessage("Ieb::setChirpStartFrequency: radar mode is not set").throwMe();
    }
  //BR_: chirp bandwidth
  //s_adc_: adc sample rate
  //ref frequency: oscillation frequency 10 MHz
  Uvar csf =slo_frequency - BR_/2.0 - s_adc_/4.0-doppler;
  setChirpStartFrequency(csf);
  }

//------------------------------
//compute chirp start frequency for various doppler values from 5 beams'
//3dB or 5dB points
// For reference, see memo 334YZ-2003-002
//------------------------------

void Ieb::computeChirpStartFrequency(const Uvar& min_doppler, const Uvar& max_doppler)
  {
    if(chirp_start_frequency_set_)
      {
	ErrorMessage("Ieb::setChirpStartFrequency: chirp start frequency is already set").throwMe();
      }
    if (!mod_set_)
      {
	ErrorMessage("Ieb::setChirpStartFrequency: radar mode is not set").throwMe();
      }
   if(min_doppler > max_doppler) ErrorMessage("Ieb.cpp::computeChirpStartFrequency: min doppler is larger than max doppler").throwMe();
    //BR_: chirp bandwidth   
    Uvar csf = slo_frequency - BR_/2.0 - s_adc_/4.0-(min_doppler+max_doppler)/2.0;    
    setChirpStartFrequency(csf);
  }

//--------------------------- 
//{get,set}slow field tfi
//---------------------------
Uvar Ieb::getSlowTfi() const
  {
    if (!s_tfi_set_)
      {
	ErrorMessage("getSlowTfi:tfi is not set").throwMe();
      }
  return(s_tfi_);
  }
void Ieb::setSlowTfi(const Uvar& s_tfi)
  {
    if (s_tfi_set_)
      {
	ErrorMessage("getSlowTfi:tfi is already set").throwMe();
      }
    if (s_tfi < Uvar(0,"s") || s_tfi > Uvar(18.20417*3600.0,"s"))
      {
	ErrorMessage e("Ieb::setSlowTfi:Slow field Tfi should be between 0 and 18.20417 hour");
	e.throwMe();
      }
    
    s_tfi_ =s_tfi;
    s_tfi_e_ = (unsigned short) (s_tfi_.getInUnits("s"));
    s_tfi_ = Uvar( double(s_tfi_e_),"s");
    s_tfi_set_ = true;
  }

//------------------------
//slow field typ
//----------------------------
unsigned int Ieb::getSlowTyp() const
  {
  return(s_typ_e_);
  }


//--------------------------
//{get,set} slow field dtn
//-------------------------
unsigned int Ieb::getDtn() const
  {
  return (s_dtn_e_);
  } 
void Ieb::setDtn(const unsigned int& dtn)
  {
  unsigned int dtn_tmp = dtn;
  if (dtn_tmp > 256) 
    {
      dtn_tmp = dtn_tmp %256;
    }
  s_dtn_e_ = dtn_tmp;
  dtn_set_=true;
  }

//-------------------------
//{get,set} slow field sin
//------------------------
unsigned int Ieb::getSin() const
  { 
  return(s_sin_e_);
  }

void Ieb::setSin(const unsigned int& sin)
  {
   unsigned int sin_tmp = sin;
   if (sin >= 256) 
     {
       sin_tmp = sin_tmp % 256;
     }
   s_sin_e_ = sin_tmp;
  }


//---------------------------------
//{get,set} slow field radar  mod
//--------------------------------
string Ieb::getModName() const
  {
  if(!mod_set_)
    {
    ErrorMessage("Ieb::getRadarModeName: No radar mode is set");
    }
  return(radar_mode_);
  }

unsigned int Ieb::getMod() const
  { 
    if (!mod_set_)
      {
	ErrorMessage("Ieb::getMod: radar mode is not set").throwMe();
      }
    return(s_mod_e_);
  }
//------------------------------------------
//compute radar mode for normal operation
//------------------------------------------
void Ieb::computeMod(const Uvar& adc, 
		    const unsigned int& csr, 
		    const unsigned int& baq, 
		    const unsigned int& bem,
		    const unsigned int& percentOfBW)
  {
  if (mod_set_)
    {
    ErrorMessage("Ieb::computeRadarMode: radar mode is already set").throwMe();
    }
// change by GH 09/10/04 for build 5
//  if (!(csr==0 ||csr==8||csr==6))
//    {
//      ErrorMessage("computeMod: (csr=0,8, or 6) only").throwMe();
//    }


  if (!(csr==0 ||csr==8||csr==6))
    {
      cout << "Ieb:computeMod:  Warning CSR = " << csr << " Normal ops is 0,8, or 6.  User should verify this is OK?" << endl;
    }



  s_csr_e_ = csr;
  csr_set_ = true;
  s_baq_e_ = baq;
  baq_set_ = true;
  s_bem_e_=bem;
  bem_set_ = true;
  
  if (adc == sarh_adc)
    {
      s_adc_=adc;
      s_adc_e_ = 2;
      adc_set_ = true; 
      BR_ = sarh_chirp_bandwidth;
      s_rcv_ = sarh_rcv_bandwidth;
      s_rcv_e_ = 2;
      rcv_set_ = true; 
      radar_mode_="SARH";
      s_mod_e_ = 3;
      if (s_csr_e_ == 8)
	{
	  radar_mode_="SARH_AGC";
	  s_mod_e_ = 11;
	}
      mod_set_ = true;  
    }
  else if (adc == sarl_adc)
    {
      s_adc_ = adc; 
      s_adc_e_ = 1;
      adc_set_ = true;
      BR_= sarl_chirp_bandwidth;
      s_rcv_ = sarl_rcv_bandwidth;
      s_rcv_e_ = 1;
      rcv_set_ = true;  
      radar_mode_="SARL";
      s_mod_e_ = 2;
      if(s_csr_e_ ==8)
	{
	  radar_mode_ = "SARL_AGC";
	  s_mod_e_ = 10;
	}
      mod_set_ = true;
    }
  else if (adc == alth_adc)
    {   
      s_adc_ = adc;
      s_adc_e_ = 3;
      adc_set_ = true;
      BR_ = alth_chirp_bandwidth;
      s_rcv_ = alth_rcv_bandwidth;
      s_rcv_e_ = 3;
      rcv_set_ = true;
      radar_mode_ = "ALTH";
      s_mod_e_ = 1;   
      if (s_csr_e_ ==8)
	{
	  radar_mode_ = "ALTH_AGC";
	  s_mod_e_ = 9;
	}
      mod_set_ = true;
    }
  else if (adc == altl_adc)
    {   
      s_adc_ = adc;
      s_adc_e_ = 0;
      adc_set_ = true;
      BR_ = altl_chirp_bandwidth;
      s_rcv_ = altl_rcv_bandwidth;
      s_rcv_e_ = 0;
      rcv_set_ = true;
      radar_mode_ ="ALTL";
      s_mod_e_ = 0;
      if(s_csr_e_ == 8)
	{
	  radar_mode_="ALTL_AGC";
	  s_mod_e_ = 8;
	}
      mod_set_ = true;
    }
  else
    {
      ErrorMessage("computeMod:no valid adc clock speed").throwMe();
    }

  // Adjust if BW is not 100.
  if (percentOfBW != 100) BR_ = BR_ * percentOfBW/100.0;

  //special case
  if (s_csr_e_ == 6)
    {
      radar_mode_ ="Radiometer only";
      s_mod_e_ = 4;
      mod_set_ = true;

      //for calibration mode: rwd = 0 tro = 0
      f_rwd_e_ = 0;
      rwd_set_ = true;
      s_tro_e_ = 0;
      tro_set_ = true;
      s_baq_e_ = 5;//8bit straight
      baq_set_ = true;
      setAttenuation(62,62,62);
      //chirp start frequency
      if (s_adc_==sarh_adc) setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz"));
      if (s_adc_==sarl_adc) setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz"));
      if (s_adc_==alth_adc) setChirpStartFrequency(Uvar(5375.0*1000.0,"Hz"));
      if (s_adc_==altl_adc) setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz"));
    }

  }

//-----------------------------------------
//set Radar mode for calibration
//------------------------------------------
void Ieb::setWTKCalibration(const Uvar& adc,
			    const unsigned int& csr,
			    const unsigned int& bem)
  {
  if (csr == 0 || csr == 8)
    {
      ErrorMessage("Ieb.cpp::setCalibration:This is setCalibration method only").throwMe();
    }
  if(csr >8)
    {
      ErrorMessage("Ieb.cpp::setCalibration:try to access calibration mode reserved by CTU").throwMe();
    }

  //for calibration mode: rwd = 0 tro = 0
  f_rwd_e_ = 0;
  rwd_set_ = true;
  s_tro_e_ = 0;
  tro_set_ = true;
  
  //rip : 35 msec
  setRip(Uvar(35*mstos,"s"));

  s_csr_e_ = csr;
  csr_set_ = true;
  s_bem_e_= bem;
  bem_set_ = true;
  s_baq_e_ = 5;//8bit straight
  baq_set_ = true;

  if (adc == sarh_adc)
    {
    s_adc_=adc;
    s_adc_e_ = 2;
    adc_set_ = true; 
    BR_ = sarh_chirp_bandwidth;
    s_rcv_ = sarh_rcv_bandwidth;
    s_rcv_e_ = 2;
    rcv_set_ = true;   
    radar_mode_="SARH";
    s_mod_e_ = 3;
    mod_set_ = true;   
    //attenuation setting
    if(s_csr_e_ == 4) 	setAttenuation(62,62,62);
    else setAttenuation(22,22,22);

    //bpd
    Uvar bpd = Uvar(200*mstos,"s");
    setBpd(bpd);

    //chirpstartfrequency
    setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz"));
    }
  else if (adc == sarl_adc)
    {
    s_adc_ = adc; 
    s_adc_e_ = 1;
    adc_set_ = true;
    BR_= sarl_chirp_bandwidth;
    s_rcv_ = sarl_rcv_bandwidth;
    s_rcv_e_ = 1;
    rcv_set_ = true;
    radar_mode_="SARL";
    s_mod_e_ = 2;
    mod_set_ = true;

    //attenuation setting
    if(s_csr_e_ == 4)	setAttenuation(62,62,62);
    else setAttenuation(20,20,20);

    //bpd
    Uvar bpd = Uvar(200*mstos,"s");
    setBpd(bpd);
    

    //chirpstartfrequency
    setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz"));
    }
  else if (adc == alth_adc)
    {   
    s_adc_ = adc;
    s_adc_e_ = 3;
    adc_set_ = true;
    BR_ = alth_chirp_bandwidth;
    s_rcv_ = alth_rcv_bandwidth;
    s_rcv_e_ = 3;
    rcv_set_ = true;
    radar_mode_ = "ALTH";
    s_mod_e_ = 1;   
    mod_set_ = true;
    //attenuation setting
    if(s_csr_e_ == 4) setAttenuation(62,62,62);
    else setAttenuation(23,23,23);

    //bpd
    Uvar bpd = Uvar(500*mstos,"s");
    setBpd(bpd);

    //chirpstartfrequency
    setChirpStartFrequency(Uvar(5375.0*1000.0,"Hz"));

    }
  else if (adc == altl_adc)
    {   
    s_adc_ = adc;
    s_adc_e_ = 0;
    adc_set_ = true;
    BR_ = altl_chirp_bandwidth;
    s_rcv_ = altl_rcv_bandwidth;
    s_rcv_e_ = 0;
    rcv_set_ = true;
    radar_mode_ ="ALTL";
    s_mod_e_ = 0;
    mod_set_ = true;
    //attenuation setting
    if(s_csr_e_ == 4)	setAttenuation(62,62,62);
    else setAttenuation(20,20,20);

    //bpd
    Uvar bpd = Uvar(500*mstos,"s");
    setBpd(bpd);
    

    //chirpstartfrequency
    setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz"));
    }
  else
    {
      ErrorMessage("setWKTCalibration: no valid adc clock speed").throwMe();
    }
  
  //for different calibration modes
  if (s_csr_e_ ==1) radar_mode_ +="Antenna Cal";
  else if (s_csr_e_ == 2) radar_mode_ +="Noise Diode Cal";
  else if (s_csr_e_ == 3) radar_mode_ +="Resistive Load Cal";
  else if (s_csr_e_ == 4) radar_mode_ +="Rerouted Chipr Cal";
  else if (s_csr_e_ == 5) radar_mode_ +="Leakage Signal Cal";
  else if (s_csr_e_ == 6)
    {
      radar_mode_ = "Radiometer only";
      s_mod_e_ = 4;
      mod_set_ = true;
    }
  else if (s_csr_e_ == 7) radar_mode_ ="Transimitt Only Cal";
  else ErrorMessage("setWKTCalibration: no valid calibration mode").throwMe();
  
 
  }

//-----------------------------------------
//GH:set Radar mode for calibration
//------------------------------------------
void Ieb::setGHCalibration(const Uvar& adc,
			    const unsigned int& csr,
			    const unsigned int& percentOfBW)
  {
  if (csr == 0 || csr == 8)
    {
      ErrorMessage("Ieb.cpp::setGHCalibration:This is setCalibration method only").throwMe();
    }
  if(csr >8)
    {
      ErrorMessage("Ieb.cpp::setGHCalibration:try to access calibration mode reserved by CTU").throwMe();
    }

  //for calibration mode: rwd = 0 tro = 0 during resistive load and reouted chirp
  //
  //if (csr == 4 || csr == 5) f_rwd_e_ = 0;
  //else f_rwd_e_ = 10;
  //rwd_set_ = true;
  // s_tro_e_ = 0;
  //tro_set_ = true;
  
  //rip : 35 msec
  //setRip(Uvar(35*mstos,"s"));

  s_csr_e_ = csr;
  csr_set_ = true;
  //s_bem_e_= 31;
  //bem_set_ = true;
  s_baq_e_ = 5;//8bit straight
  baq_set_ = true;

  if (adc == sarh_adc)
    {
    s_adc_=adc;
    s_adc_e_ = 2;
    adc_set_ = true; 
    BR_ = sarh_chirp_bandwidth;
    s_rcv_ = sarh_rcv_bandwidth;
    s_rcv_e_ = 2;
    rcv_set_ = true;   
    radar_mode_="SARH";
    s_mod_e_ = 3;
    mod_set_ = true;   
    s_bem_e_= 31;
    //attenuation setting
    if(s_csr_e_ == 4) 	setAttenuation(62,62,62);
    else setAttenuation(22,22,22);

    //bpd
    //Uvar bpd = Uvar(200*mstos,"s");
    //setBpd(bpd);

    //chirpstartfrequency
    //setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz"));
    }
  else if (adc == sarl_adc)
    {
    s_adc_ = adc; 
    s_adc_e_ = 1;
    adc_set_ = true;
    BR_= sarl_chirp_bandwidth;
    s_rcv_ = sarl_rcv_bandwidth;
    s_rcv_e_ = 1;
    rcv_set_ = true;
    radar_mode_="SARL";
    s_mod_e_ = 2;
    mod_set_ = true;   
    s_bem_e_= 31;

    //attenuation setting
    if(s_csr_e_ == 4)	setAttenuation(62,62,62);
    else setAttenuation(20,20,20);

    //bpd
    //Uvar bpd = Uvar(200*mstos,"s");
    //setBpd(bpd);
    

    //chirpstartfrequency
    //setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz"));
    }
  else if (adc == alth_adc)
    {   
    s_adc_ = adc;
    s_adc_e_ = 3;
    adc_set_ = true;
    BR_ = alth_chirp_bandwidth;
    s_rcv_ = alth_rcv_bandwidth;
    s_rcv_e_ = 3;
    rcv_set_ = true;
    radar_mode_ = "ALTH";
    s_mod_e_ = 1;   
    mod_set_ = true;   
    s_bem_e_= 4;
    //attenuation setting
    if(s_csr_e_ == 4) setAttenuation(62,62,62);
    else setAttenuation(62,30,62);

    //bpd
    //Uvar bpd = Uvar(500*mstos,"s");
    //setBpd(bpd);

    //chirpstartfrequency
    //setChirpStartFrequency(Uvar(5375.0*1000.0,"Hz"));

    }
  else if (adc == altl_adc)
    {   
    s_adc_ = adc;
    s_adc_e_ = 0;
    adc_set_ = true;
    BR_ = altl_chirp_bandwidth;
    s_rcv_ = altl_rcv_bandwidth;
    s_rcv_e_ = 0;
    rcv_set_ = true;
    radar_mode_ ="ALTL";
    s_mod_e_ = 0;
    mod_set_ = true;  
    s_bem_e_= 4;
    //attenuation setting
    if(s_csr_e_ == 4)	setAttenuation(62,62,62);
    else setAttenuation(62,9,62);

    //bpd
    //Uvar bpd = Uvar(500*mstos,"s");
    //setBpd(bpd);
    

    //chirpstartfrequency
    //setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz"));
    }
  else
    {
      ErrorMessage("setGHCalibration: no valid adc clock speed").throwMe();
    }


  // Adjust if BW is not 100.
  if (percentOfBW != 100) BR_ = BR_ * percentOfBW/100.0;

  //chirp start frequency
  // start + (BR difference)/2
  //   if BW = 100, then BR_diff = 0, start equal same as blue book
  //   if BW = 0, then BR_diff = BR/2, start equal center of band
  if (s_adc_==sarh_adc) setChirpStartFrequency(Uvar(9075.0*1000.0,"Hz") + (sarh_chirp_bandwidth-BR_)/2.0 ) ;
  if (s_adc_==sarl_adc) setChirpStartFrequency(Uvar(9538.0*1000.0,"Hz") + (sarl_chirp_bandwidth-BR_)/2.0 );
  if (s_adc_==alth_adc) setChirpStartFrequency(Uvar(5375.0*1000.0,"Hz") + (alth_chirp_bandwidth-BR_)/2.0 );
  if (s_adc_==altl_adc) setChirpStartFrequency(Uvar(9884.0*1000.0,"Hz") + (altl_chirp_bandwidth-BR_)/2.0 );




  
  //for calibration mode: rwd = 0 tro = 0 during resistive load and reouted chirp
  //setNumberOfPulses(8);
  if (csr == 4 || csr == 5) f_rwd_e_ = 0;
  else f_rwd_e_ = 8 + 2;
  rwd_set_ = true;
  s_tro_e_ = 0;
  tro_set_ = true;
  
  //for different calibration modes
  if (s_csr_e_ ==1) radar_mode_ +="Antenna Cal";
  else if (s_csr_e_ == 2) radar_mode_ +="Noise Diode Cal";
  else if (s_csr_e_ == 3) radar_mode_ +="Resistive Load Cal";
  else if (s_csr_e_ == 4) {radar_mode_ +="Rerouted Chirp Cal";  s_bem_e_= 0;}
  else if (s_csr_e_ == 5) radar_mode_ +="Leakage Signal Cal";
  else if (s_csr_e_ == 6)
    {
      radar_mode_ = "Radiometer only";
      s_mod_e_ = 4;
      mod_set_ = true;
    }
  else if (s_csr_e_ == 7) radar_mode_ ="Transmit Only Cal";
  else ErrorMessage("setGHCalibration: no valid calibration mode").throwMe();
  
  bem_set_ = true;

  }


//-----------------------------------------
//GH set to a default radiometer ILX.  Used when Ideal want active mode, but power not available
//------------------------------------------
void Ieb::defaultRadiometer(const Time& iebstarttime, const Time& ilxexecutiontime)

  {
    // This will set up at the time stamps
    resyncTimeTfi ( iebstarttime, ilxexecutiontime);

    unsigned int radiometer_only_baq = 5;
    unsigned int radiometer_only_csr = 6;
    unsigned int radiometer_only_bem = 4;

    // no matter what override some set parameters to become radiometer
    //clear();  this may be a better way......
    mod_set_ = false; 
    chirp_start_frequency_set_ = false;
    csd_set_ = false;
    csq_set_ = false;
    chirp_frequency_step_set_ = false;
    pul_set_ = false;
    bpd_set_ = false;
    rip_set_ = false;
    rad_set_ = false;

    computeMod(altl_adc,radiometer_only_csr,radiometer_only_baq,radiometer_only_bem, 100);
    setPri( Uvar(1.0/1000, "s") );
    setDutycycle(0.38);
    computeWaveform();
    setNumberOfPulses(0);     
    setBpd(Uvar(1100*mstos,"s"));
    setRip(Uvar(35*mstos,"s"));
    //setRL_ND_RIP(Uvar(35,"ms"), Uvar(5,"ms") );
    setRL_ND(Uvar(35*mstos,"s"), Uvar(5*mstos,"s") );
    computeRad();
   
  }


//-----------------------------------------
//GH set to a default Interleave ILX.  Used when Ideal want active mode, but power not available
//------------------------------------------
void Ieb::interleave( const int input_value )

  {
    unsigned int temp_mod;
    unsigned int old_rwd;
    unsigned int old_pul;
    Uvar old_rwd_time(0,"s");
    Uvar old_pri (0,"s");
    // ALTH (case 1,9)expected interleave is to adjust altimeter to insure receive first 9or last burst)
    // altimeter high (1 and 9) should be set to 21 PUL, TRO -6, RWD = RWD +3 (to center)
    // special case will set PUL = 9, TRO = 6, RWD = RWD -6 to get 3 pass where it should be...

    // SAR (case 2,3,10,11) increment PRF by 10% and then recalc waveform. Maintain ~70%
    
    temp_mod = s_mod_e_;
    if (temp_mod == 3 ||temp_mod == 10 || temp_mod == 11) temp_mod = 2;  // all SARs are equal! 
    if (temp_mod == 9 ) temp_mod = 1;  // all ATLHs are equal! 
    if (temp_mod == 8 ) temp_mod = 0;  // all ATLLs are equal! 

    switch ( temp_mod ){
    case (0):
	cout << "interleavestart0: Receive ONLY.  perform special interleave: CSR:RWD:PUL:TRO " <<getCsr()<<":"<<getRwdInPriUnits()<<":"<<getPul()<<":"<<getTroInPriUnits()<< endl;
	//rwd_set_ = false; setRwdInUnitsOfPri( 0 );
	s_csr_e_ = 1; csr_set_ = true; // set to antenna calibration
	cout << "interleavestart0: Receive ONLY ........completed interleave: CSR:RWD:PUL:TRO " <<getCsr()<<":"<<getRwdInPriUnits()<<":"<<getPul()<<":"<<getTroInPriUnits()<< endl;
      break;
    case (1): // copy to case of 9 as well
      if (getTroInPriUnits() == -6 && getPul() >= 14) {
	cout << "interleavestart1: TRO is -6, PUL >= 14, perform special interleave:PUL:RWD:TRO " <<getPul()<<":"<<getRwdInPriUnits()<<":"<<getTroInPriUnits()<< endl;
	tro_set_ = false; setTroInUnitsOfPri(6);
	old_pul = getPul(); pul_set_ = false; setNumberOfPulses(old_pul-12);
	old_rwd = getRwdInPriUnits();
	rwd_set_ = false; setRwdInUnitsOfPri(old_rwd - 6 );
	//s_pul_e_ = 9;
	//s_rwd_e_ = s_rwd_e_ - 6;
	cout << "interleavestart1: TRO is -6, perform special interleave:PUL:RWD:TRO " <<getPul()<<":"<<getRwdInPriUnits()<<":"<<getTroInPriUnits()<< endl;
      } else { 
	cout <<"Warning: interleave: interleave for altimeterhigh called but TRO original not -6 & PUL not >= 14.  NO special interleave performed" << endl;
      }
      break;
    case (2):
      cout << "interleavestart(SARs)"<<input_value<<":interleave:PRF Hopping:TFI: " <<getSlowTfi() << ":PUL:"<<getPul()<<":RWD:"<<getRwdInPriUnits()<< ":RWD:" <<getRwd()<<":PRI:"<<getPriEncodedValue()<<":PRI:"<<getPri();
      old_rwd = getRwdInPriUnits();
      old_rwd_time = getRwd();
      old_pri = getPri();
      pri_set_ = false; setPri(old_pri * (1.0 - input_value/100.0) ); // typically caused 10% lower PRI (higher PRF) (duty cycle will go up)
      rwd_set_ = false; setRwd(old_rwd_time);
      csd_set_ = false; csq_set_ = false; chirp_frequency_step_set_ = false; computeWaveform();
      if (input_value < 0) { // increasing PRI, may also violate 7% duty cycle and create PUL to be larger than RWD.
	// change 7%
	Uvar bpd_based_on_burst_dutycycle = double (f_pul_e_)* 
	  double(s_csq_e_ +1) * s_csd_ /(6.95/100.0);
	if (bpd_based_on_burst_dutycycle > getBpd() ) {
	  cout << " Warning: Ieb.cpp 7% dutycycle violation.  Increasing Bpd: was:" << getBpd() << " :Now = " ;
	  bpd_set_ = false; setBpd(bpd_based_on_burst_dutycycle);
	  cout << getBpd() << endl ;
	}
	// change RWD
	if ( getPul() >= getRwdInPriUnits() ) {
	  cout << "Warning Ieb.cpp PUL>RWD violation.  Changing PUL to = RWD -1.  Pul was:" << getPul() << " :Now =" ;
	  pul_set_ = false ; setNumberOfPulses ( getRwdInPriUnits() - 1 );
	  cout << getPul () << endl;
	}

      }

      
      //pri_set_ = false; setTroInUnitsOfPri(6);
      //pul_set_ = false; setNumberOfPulses(9);
      //rwd_set_ = false; setRwdInUnitsOfPri(getRwdInPriUnits() -6 );
      //s_pul_e_ = 9;
      //s_rwd_e_ = s_rwd_e_ - 6;
      cout << " NEW> RWD:"<<getRwdInPriUnits()<< ":RWD:" <<getRwd()<<":PRI:"<<getPriEncodedValue()<<":PRI:"<<getPri()<<endl;
      //} else { 
      //	cout <<"Warning: interleave11: interleave for SARH called but TRO original not 6.  NO special interleave performed" << endl;
      // }
      break;
    case (3): // actually all SARs are equal, so should be handle by case2
      cout << "ERROR:(should never get here3)  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (4):
      cout << "Error:  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (5):
      cout << "Error:  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (6):
      cout << "Error:  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (7):
      cout << "Error:  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (8):// actually all ALTLs are equal, so should be handle by case0
      cout << "ERROR:(should never get here8)  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (9):// actually all ALTHs are equal, so should be handle by case1
      cout << "ERROR:(should never get here9)  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (10): // actually all SARs are equal, so should be handle by case2
      cout << "ERROR:(should never get here10)  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    case (11): // actually all SARs are equal, so should be handle by case2
      cout << "ERROR:(should never get here11)  No interleave function defined for this mode: " << s_mod_e_ << endl;
      break;
    default:
      ErrorMessage("ERROR: Ieb::interleave: decoded radar mode does not have an appropriate value").throwMe();
    } // end of s_mod_e switch

   
  }

//-----------------------------------------
//{get,set}slow field calibration source
//---------------------------------------
void Ieb::disableAutogain()
  {
    if(!csr_set_)
      {
	ErrorMessage("Ieb::disableAutogain: csr is not set.").throwMe();
      }else {
	if (s_csr_e_ >= 8){
	  s_csr_e_ -= 8; // subtract 8 from an autogain mode, will give same mode, but autogain disabled
	} else {
	  ErrorMessage("Ieb::disableAutogain: Cannot change CSR.  Autogain was not set.").throwMe();
	}
	if (s_mod_e_ >= 8){
	  s_mod_e_ -= 8; // subtract 8 from an autogain mode, will give same mode, but autogain disabled
	} else {
	  ErrorMessage("Ieb::disableAutogain: Cannot change MOD.  Autogain was not set.").throwMe();
	}
      }
  }

//-----------------------------------------
//{get,set}slow field calibration source
//---------------------------------------
unsigned int Ieb::getCsr() const
  {
    if(!csr_set_)
      {
	ErrorMessage("Ieb::getCsr: no csr is set").throwMe();
      }
    return(s_csr_e_);
  }

//---------------------------
//{get} Adc rate and Receiver bandwidth
// only get method is provided because they are determined by radar mode
//---------------------------
Uvar Ieb::getAdc() const
{
  if(!mod_set_)
    {
      ErrorMessage("Ieb::getAdc: radar mode is not set").throwMe();
    }
  return(s_adc_);
}

unsigned int Ieb::getAdcEncodedValue()  const
{
  if(!mod_set_)
    {
      ErrorMessage("Ieb::getEncodedAdcValue: radar mode is not set").throwMe();
    }
  return(s_adc_e_);
}

//--------------------------------
//{get} Rcv value(eng unit or encoded value)
//------------------------------
Uvar Ieb::getRcv()  const
  {
  if(!mod_set_)
    {
      ErrorMessage("Ieb::getRcv: radar mode is not set").throwMe();
    }
  return(s_rcv_);  
  }

unsigned int Ieb::getRcvEncodedValue()  const
  {
  if(!mod_set_)
    {
      ErrorMessage("Ieb::getRcvEncodedValue: radar mode is not set").throwMe();
    }
  return(s_rcv_e_);  
  }
//----------------------------
//{get,set} Tro
//-----------------------------
Uvar Ieb::getTro()  const
{
  if(!tro_set_)
    {
      ErrorMessage("Ieb::getTro: tro is not set").throwMe();
    }
  return(s_tro_);
}

int  Ieb::getTroInPriUnits()
{
  if(!tro_set_)
    {
      ErrorMessage("Ieb::getTro: tro is not set").throwMe();
    }
  
  int tro_in_pri = (int) s_tro_e_;
  if (tro_in_pri >= 8 && tro_in_pri <16)
    {
      tro_in_pri = tro_in_pri - 16;
    }
  if (tro_in_pri >= 16)
    {
      ErrorMessage("Ieb.cpp: encoded tro value is larger than 15 -> out of range").throwMe();
    }
  return(tro_in_pri);
}


void Ieb::setTroTime(const Uvar& tro_time)
  {
  if (tro_set_)
    {
      ErrorMessage("Ieb::setTro: tro is already set").throwMe();
    }
  if (!pri_set_)
    {
      ErrorMessage("Ieb:setTro: pri is not set to convert tro in units of pri").throwMe();
    }
  s_tro_ = tro_time;
  int s_tro_in_pri =  int ((s_tro_ /getPri()).getInUnits(""));
 
  if (s_tro_in_pri> 7 || s_tro_in_pri < -8) 
    {
      cout<<"stro time "<<s_tro_<<endl;
      cout<<"pri "<<getPri()<<endl;
      cout<<"s_stro_in_pri_ "<<s_tro_in_pri<<endl;
      cout << "Ieb: Warning --> TRO value is out of range:see BB12-27. " ;
      s_tro_in_pri = 7;
      cout << " --> TRO force to " << s_tro_in_pri  << endl; 
	//ErrorMessage("TRO value is out of range:see BB12-27").throwMe();
    }  
  
  if (s_tro_in_pri >= 0 && s_tro_in_pri < 8 ) 
    {
      s_tro_e_ = (unsigned short) s_tro_in_pri;
    }
  else if (s_tro_in_pri < 0 && s_tro_in_pri > -9)
    {
      s_tro_e_ = (unsigned short)(s_tro_in_pri+16);
    }  
  s_tro_ =  f_pri_ * double(s_tro_e_);
  if (s_tro_in_pri < 0 && s_tro_in_pri > -9)
    {
      s_tro_ = f_pri_ * (double(s_tro_e_) - 16.0);
    }  
  tro_set_=true;  
  }



void Ieb::setTroInUnitsOfPri(const int& tro_number)
  {
  if (tro_set_)
    {
      ErrorMessage("Ieb::setTro: tro is already set").throwMe();
    }
  if (!pri_set_)
    {
      ErrorMessage("Ieb:setTro: pri is not set to convert tro in units of pri").throwMe();
    }
 
  if (tro_number > 7 || tro_number < -8) 
    {
      cout << "Ieb(setTRO): Warning --> TRO value is out of range:see BB12-27. Was: " << tro_number ;
      s_tro_e_ = 7 ; //tro_number = 7;
      cout << " --> TRO force to " << s_tro_e_  << endl;  
        //ErrorMessage("TRO value is out of range:see BB12-27").throwMe();
    }  
  
  if (tro_number >= 0 && tro_number < 8 ) 
    {
      s_tro_e_ = (unsigned short) (tro_number);
    }
  else if (tro_number < 0 &&tro_number > -9)
    {
      s_tro_e_ = (unsigned short)(tro_number+16);
    }  
  s_tro_ =  getPri() * double(tro_number);
  tro_set_=true;  
  }



void Ieb::computeTro(const Uvar& min_range, const Uvar& max_range)
{
  if (tro_set_)
    {
      ErrorMessage("Ieb::computeTro: tro is already set").throwMe();
    }
  if (!pri_set_)
    {
      ErrorMessage("Ieb:setTro: pri is not set to convert tro in units of pri").throwMe();
    }
  if (min_range > max_range) ErrorMessage("Ieb.cpp::setTro: min range is larger than max range").throwMe();
  Uvar tro = (max_range - min_range) * 2.0/speed_light;
  setTroTime(tro);
}


//------------------
//adjust pul based on ctrx number
// this happens only when long term dutycycle limit is violated
// and the result is to have less number of active pulses
//------------------

void Ieb::resetTroBasedOnCTRX(const unsigned int& ctrx)
  {
    if (!tro_set_)
      {
	ErrorMessage("Ieb::resetTroBasedOnCTRX: tro should have been set").throwMe();
      }   
    if (!pri_set_)
      {
	ErrorMessage("Ieb:setTro: pri is not set to convert tro in units of pri").throwMe();
      }
    if(!pul_set_)
      {
	ErrorMessage("Ieb.cpp:pul is not set").throwMe();
      }
 
    int pul = int(getPul());
    int tro_in_pri = getTroInPriUnits();
    int original_window= pul + tro_in_pri;
    int new_window = (int) ctrx;
    int diff = new_window - original_window;
    if(diff==0) return;
    else
      {
	cout<<"special case: Tro will be reset based on CTRX"<<endl;
	tro_set_= false;//will recompute tro based on ctrx
	tro_in_pri += diff;
	setTroInUnitsOfPri(tro_in_pri);
      }
  }

//----------------
//{get,set} BAQ
//----------------
unsigned int Ieb::getBaq() const
  {
    if(!baq_set_)
      {
	ErrorMessage("Ieb::getBaqMode:no baq is set").throwMe();
      }
  return(s_baq_e_);
  }



//-------------------------------------
//{get} beam mask and number of beams used
//-----------------------------------

unsigned int Ieb::getBem() const
{
  if(!bem_set_) 
    {
      ErrorMessage("Ieb::getBem: beam mask is not set").throwMe();
    }
  return(s_bem_e_);
}

unsigned int Ieb::getNumberOfBeams()
{
  if(!bem_set_) 
    {
      ErrorMessage("Ieb::getBem: beam mask is not set").throwMe();
    }
  unsigned short beam_mask;
  bitset(beam_mask,0,4,s_bem_e_);
  unsigned int N = 0;
  for (unsigned int i = 0; i < 5;++i)
    {
      if (bitget(beam_mask,i,i)==1) N++;
    }
  // in Rerouted chirp mode (CSR = 4) bems should = 0
  if (N == 0 && s_csr_e_ != 4) ErrorMessage("Ieb.cpp:getNumberOfBeams: no beam is selected").throwMe();
  if (N > 5) ErrorMessage("Ieb.cpp:getNumberOfBeams: more than 5 beams are selected").throwMe();
  return(N);
}


//-----------------------------------------
//compute all three attenuations and set
// for a given input of all 5 beams' echoes
//-----------------------------------------
//void Ieb::computeAttenuationSettings(const Uvar& Pn, const Uvec& Echo_power,const Uvar& squared_deviation_for_known_source)
//{
//if(!bem_set_) ErrorMessage("Ieb.cpp::computeAttnuationSetting:computeAttenuationSettings: no beam setting").throwMe();

  //determine total att for at1, at3 and at4
  //unsigned int att1,att3,att4;
  //unsigned short bem_mask;
  //bitset(bem_mask,0,4,s_bem_e_);
  //
  //Array1D<unsigned int> beam_on("beam on mask",5);
  //for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
  //  {
//beam_on(i_beam)=bitget(bem_mask,i_beam,i_beam);
//  }
//
//if (beam_on(0) == 0 && beam_on(1) == 0)
//  { 
//  att1 = 62;
//  }
//else
//  {
//  att1 =  getSingleBeamAttenuationSetting( Pn,Echo_power(0),
//		  squared_deviation_for_known_source);
//  }
//
//if (beam_on(2) == 0)
//  {
//  att3 = 62;
//  }
//else
//  {
//  att3 =   getSingleBeamAttenuationSetting( Pn,Echo_power(2),
//		   squared_deviation_for_known_source);
//  } 
//
//
//if (beam_on(3) == 0 && beam_on(4)==0) 
//  {
//  att4 =62;
//  }
//else
//  {
//  att4 = getSingleBeamAttenuationSetting(Pn,Echo_power(4),
//		 squared_deviation_for_known_source);
//  }  
//setAttenuation(att1,att3,att4);
//}



//-----------------------
//set attenuation by adjusting
// thermal noise equal to
// power corresponding to the noise bit
//------------------------
void Ieb::computeAttenuationSettings(const double&  noise_bit,
				     const Uvar& Pn,
				    const Uvar& power_to_variance)
  {
  if(!bem_set_) ErrorMessage("Ieb.cpp::computeAttnuationSetting:computeAttenuationSettings: no beam setting").throwMe();

  //determine total att for at1, at3 and at4
  unsigned int att1,att3,att4;
  unsigned short bem_mask;
  bitset(bem_mask,0,4,s_bem_e_);

  //if(noise_bit>7.0) ErrorMessage("noise bit set is too hight : "+toStr(noise_bit)).throwMe();
  // if(noise_bit<1.0) ErrorMessage("noise bit set is too low").throwMe();

  double bit_to_number = pow(2,noise_bit);
  double max_variance = bit_to_number * bit_to_number;
  double gain= max_variance/(power_to_variance*Pn).getInUnits("");
 
  //how to get db for input gain
  //int gain_dB =cag_table_.findGaindB(gain);
  if(gain>=1) ErrorMessage("attenuation gain is larger than 1  "+toStr(gain)).throwMe();
  double gain_dB_d = -10.0*log(gain)/log(10.0);
  if(gain_dB_d <0) ErrorMessage("gain db value is negative "+toStr(gain_dB_d)).throwMe();
  unsigned int gain_dB =(unsigned int) round_double(gain_dB_d);   
  unsigned int gain_set_dB = (unsigned int) gain_dB;

  Array1D<unsigned int> beam_on("beam on mask",5);
  for (unsigned int i_beam = 0; i_beam < 5;++i_beam)
    {
      beam_on(i_beam)=bitget(bem_mask,i_beam,i_beam);
    }
 
  // if any of beam 1,2,4,5 is used, set both att1 and att4 to the same gain
  if (beam_on(0) == 0 && beam_on(1) == 0 && beam_on(3)==0 && beam_on(4)==0)
    { 
      att1 = 62;
      att4 = 62;
    }
  else
    {
      att1 = gain_set_dB;
      att4 = gain_set_dB;
    }
  
  //beam 3
  if (beam_on(2) == 0)
    {
      att3 = 62;
    }
  else
    {
      att3 =   gain_set_dB;
    } 
  
  setAttenuation(att1,att3,att4);

  }
//----------------------------------
//Compute single beam attenuation setting
//------------------------------------
//unsigned int  Ieb::getSingleBeamAttenuationSetting(const Uvar& Pn, 	
//				     const Uvar& Echo_power,
//		     const Uvar& squared_deviation_for_known_source)
//{
//unsigned int att_set=62;
//Dvec echo_squared_deviation("squared deviation buffer",75);
//Dvec noise_squared_deviation("noise term ",75);
//    
////squared deviation from noise
//for (unsigned int i_att = 0; i_att < 75;++i_att)
//  {
//  noise_squared_deviation(i_att)
//    =radar_buffer_squared_deviation(Pn,
//			      i_att,
//			      squared_deviation_for_known_source);
//    }
// 
//for (unsigned int i_att = 0; i_att < 75;++i_att)
//  {
//    echo_squared_deviation(i_att)
//=radar_buffer_squared_deviation(Echo_power,
//				i_att,
//				squared_deviation_for_known_source);
//  }
// 
//for (unsigned int i_att = 0; i_att < 75;++i_att)
//  {
//  //saturation limit: 3sigma of input signal will be captured withtout
//  //saturation   
//  double radar_buffer_squared_deviation = 
//    echo_squared_deviation(i_att)*9.0
//    + noise_squared_deviation(i_att)+ 0.29*0.29;
//  if ( radar_buffer_squared_deviation < (255.0/2.0)*(255.0/2.0))
//    {
//att_set = i_att;
//break;
//    }
//  }
//return(att_set);
//}

//-------------------------------------------------
//Set attenuation setting using attenuation mapping
//----------------------------------------------------
 void Ieb::setAttenuation(const unsigned int& at1, const unsigned int& at3, const unsigned int& at4)
  {
    if (at1 > 74 || at3 > 74 || at4 > 74) ErrorMessage("One of total attenuation value is larger than 74dB").throwMe();
    computeAttenuationMap(s_at1N1_dB_, s_at1N2_dB_, s_at1N3_dB_,at1);
    computeAttenuationMap(s_at3N1_dB_, s_at3N2_dB_, s_at3N3_dB_,at3);
    computeAttenuationMap(s_at4N1_dB_, s_at4N2_dB_, s_at4N3_dB_,at4);
    at1N1_set_ = true;
    at1N2_set_ = true;
    at1N3_set_ = true;

    at3N1_set_ = true;
    at3N2_set_ = true;
    at3N3_set_ = true;

    at4N1_set_ = true;
    at4N2_set_ = true;
    at4N3_set_ = true;
  }


//----------------------------------------------------------------------
// Return expected radar buffer's squared deviation
//   for a given input power and attenuation setting in dB scale.
// This calculation needs "measured" squared deviation with a known noise 
//   source such as system noise
//
//--------------------------------------------------------------------
//double Ieb::radar_buffer_squared_deviation(const Uvar& Power_input, 
//       const unsigned int& attenuation_in_dB,
//       const Uvar& squared_deviation_for_known_input_noise_power)
//{
//  double a = (Power_input *squared_deviation_for_known_input_noise_power).getInUnits("");
//  double gain =  cag_table_.gainLinearScale((int)at1_db);
//  return( 0.29*0.29+ a*gain);
//}

//------------------------------
//{get} receiver attenuation setting
//------------------------------
unsigned int Ieb::getAt1()
{
  if(!at1N1_set_ || !at1N2_set_ || !at1N3_set_)
    {
      ErrorMessage("Ieb::getAt1: no attenuation set").throwMe();
    }
  return(s_at1N1_dB_ + s_at1N2_dB_ + s_at1N3_dB_);
}

unsigned int Ieb::getAt1N1dB()  const
{
  if(!at1N1_set_)
    {
      ErrorMessage("Ieb::getAt1N1dB(): no attenuation is set").throwMe();
    }
  return(s_at1N1_dB_);
}

unsigned int Ieb::getAt1N2dB()  const
{
  if(!at1N2_set_)
    {
      ErrorMessage("Ieb::getAt1N2dB(): no attenuation is set").throwMe();
    }
  return(s_at1N2_dB_);
}

unsigned int Ieb::getAt1N3dB()  const
{
  if(!at1N3_set_)
    {
      ErrorMessage("Ieb::getAt1N3dB(): no attenuation is set").throwMe();
    }
  return(s_at1N3_dB_);
}


unsigned int Ieb::getAt3()
{
  if(!at3N1_set_ || !at3N2_set_ || !at3N3_set_)
    {
      ErrorMessage("Ieb::getAt3: no attenuation set").throwMe();
    }
  return(s_at3N1_dB_ + s_at3N2_dB_ + s_at3N3_dB_);
}


unsigned int Ieb::getAt3N1dB()  const
{
  if(!at3N1_set_)
    {
      ErrorMessage("Ieb::getAt3N1dB(): no attenuation is set").throwMe();
    }
  return(s_at3N1_dB_);
}


unsigned int Ieb::getAt3N2dB() const
{
  if(!at3N2_set_)
    {
      ErrorMessage("Ieb::getAt3N2dB(): no attenuation is set").throwMe();
    }
  return(s_at3N2_dB_);
}

unsigned int Ieb::getAt3N3dB() const
{
  if(!at3N3_set_)
    {
      ErrorMessage("Ieb::getAt3N3dB(): no attenuation is set").throwMe();
    }
  return(s_at3N3_dB_);
}

unsigned int Ieb::getAt4()
{
  if(!at4N1_set_ || !at4N2_set_ || !at4N3_set_)
    {
      ErrorMessage("Ieb::getAt4: no attenuation set").throwMe();
    }
  return(s_at4N1_dB_ + s_at4N2_dB_ + s_at4N3_dB_);
}

unsigned int Ieb::getAt4N1dB() const
{
  if(!at4N1_set_)
    {
      ErrorMessage("Ieb::getAt4N1dB(): no attenuation is set").throwMe();
    }
  return(s_at4N1_dB_);
}



unsigned int Ieb::getAt4N2dB() const
{
  if(!at4N2_set_)
    {
      ErrorMessage("Ieb::getAt4N2dB(): no attenuation is set").throwMe();
    }
  return(s_at4N2_dB_);
}

unsigned int Ieb::getAt4N3dB() const
{
  if(!at4N3_set_)
    {
      ErrorMessage("Ieb::getAt4N3dB(): no attenuation is set").throwMe();
    }
  return(s_at4N3_dB_);
}


void Ieb::getAttenuationBitMask(unsigned int& at1_bit_mask,
				unsigned int& at3_bit_mask,
				unsigned int& at4_bit_mask)
  {
    if(!slowfastfield_loaded_)
      ErrorMessage("Ieb.cpp:getAttenuationBitMask:no slow field is loaded").throwMe();
    at1_bit_mask=0;
    at3_bit_mask=0;
    at4_bit_mask=0;
    
    //get at1
    unsigned short tmp=0;
    tmp = bitget(slowfield_[4],0,11);
    at1_bit_mask= tmp;

    //get at3
    tmp=0;
    tmp=bitget(slowfield_[5],0,11);
    at3_bit_mask=tmp;

    //get at4
    tmp=0;
    tmp = bitget(slowfield_[6],4,15);
    at4_bit_mask=tmp;
      
    //see the following decoding method
    //s_at1N3_e_ = bitget(slowfield_[4],0,1);
    //s_at1N3_dB_ = s_at1N3_e_*4;
    //s_at1N2_e_= bitget(slowfield_[4],2,6);
    //s_at1N2_dB_ = s_at1N2_e_;
    //s_at1N1_e_ = bitget(slowfield_[4],7,11);
    //s_at1N1_dB_ = s_at1N1_e_;
  
 
    
    //s_at3N3_e_ = bitget(slowfield_[5],0,1);
    //s_at3N3_dB_ = s_at3N3_e_ *4 ;
    //s_at3N2_e_ = bitget(slowfield_[5],2,6);
    //s_at3N2_dB_ = s_at3N2_e_;
    //s_at3N1_e_ = bitget(slowfield_[5],7,11);
    //s_at3N1_dB_ = s_at3N1_e_;
    
    //s_at4N3_e_ = bitget(slowfield_[6],4,5);
    //s_at4N3_dB_ = s_at4N3_e_ *4;
    //s_at4N2_e_ =bitget(slowfield_[6],6,10);
    //s_at4N2_dB_ = s_at4N2_e_;
    //s_at4N1_e_ = bitget(slowfield_[6],11,15);
    //s_at4N1_dB_ = s_at4N1_e_;

  }
//-----------------------------
//get AttenuationdBofBeamNumber()
//-----------------------------
unsigned int Ieb::getAttenuationdB(const unsigned int& beam_number)
  {
    if(beam_number < 1 || beam_number > 5) ErrorMessage("beam number is out of range: should be between 1 and 5").throwMe();
    unsigned int att_db;
    if(beam_number ==1 || beam_number==2) att_db=getAt1();
    else if(beam_number==3) att_db=getAt3();
    else if (beam_number==4 || beam_number==5) att_db=getAt4();
    else ErrorMessage("this should not happen").throwMe();
    return(att_db);
  }
double Ieb::getAttenuation(const unsigned int& beam_number) 
{
  double att_db;
  att_db = double( getAttenuationdB(beam_number));
  return(  pow(10, -att_db/10.0));
}
//---------------------------
//{get,set} Radiometer integration time
//---------------------------    
Uvar Ieb::getRip() const
  {
  if(!rip_set_)
    {
      ErrorMessage("Ieb::getRip: no rip is set").throwMe();
    }
  return(s_rip_);
  }

unsigned int Ieb::getRipEncodedValue() const
  {
  if(!rip_set_)
    {
      ErrorMessage("Ieb::getRipEncodedValue(): no rip is set").throwMe();
    }
  return(s_rip_e_);
  }

void Ieb::setRip(const Uvar& rip)
{
  if(rip_set_)
    {
      ErrorMessage("Ieb::setRip: rip is already set");
    }
  s_rip_= rip;
  if (s_rip_ < Uvar(10*mstos,"s") || s_rip_ > Uvar(75*mstos,"s"))
    {
      ErrorMessage("Ieb::setRip: rip value is out of range").throwMe();
    }
  s_rip_e_ = (unsigned short)round_double((s_rip_/Uvar(5*mstos,"s")).getInUnits(""));
  s_rip_ = Uvar(5.0*double(s_rip_e_)*mstos,"s") ;
  rip_set_ = true;
}

void Ieb::computeRip()
  {
    cout<<"----------------------------------------------------------"<<endl;
    cout<<"default value for active mode and 35 msec for calibration"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    setRip(Uvar(35*mstos,"s"));   
  }

//------------------------------
//{get,set} radiometer window count
//------------------------------
unsigned int Ieb::getRad() const
  {
  if(!rad_set_)
    {
      ErrorMessage("Ieb::getRad: no radiometer window counter is set").throwMe();
    }
  return(s_rad_e_);
  } 


void Ieb::setRad(const unsigned int& rad)
{
  if(rad_set_)
    {
      ErrorMessage("Ieb::setRad: rad is already set").throwMe();
    }
  if(rad==0) ErrorMessage("Ieb::setRad: radiometer window count is zero").throwMe();
  if (rad > 255)
    {
      ErrorMessage("Ieb::setRad: rad value is out of range").throwMe();
    }
  s_rad_e_ = rad;
  rad_set_ = true;
}

void Ieb::computeRad()
  {//assume bpd is set
    if(!pri_set_) ErrorMessage("Ieb.cpp::computeRad:no pri set").throwMe();
    if(!tro_set_) ErrorMessage("Ieb.cpp::computeRad:no tro set").throwMe();
    if(!bpd_set_) ErrorMessage("Ieb.cpp::computeRad:no bpd is set ").throwMe();
    Uvar TA = Uvar(1.0*mstos,"s")+getPri();
    Uvar TRDSUND = Uvar(2.0*mstos,"s");
    Uvar Buffer =Uvar(2*mstos,"s");
    Uvar time_left = getBpd() - TA - TRDSUND - Buffer;


    //hamilton 12/02/2003
    //not dependant on csr, use smae equation for calibration modes????
    // if (s_csr_e_ == 0 || s_csr_e_ == 8)
    // {//either normal operation or auto gain control    
	time_left -= getRwd()+getPul()*getPri()+getTroInPriUnits()*getPri()+getRL()+getND();

	//	if (time_left < getRip() ) {
	//cout << "Warning:  need to increase BPD "<< getBpd << " for RAD.  Pad(withoutRAD)= " << time_left << " RIP+ " <<getRip() << endl;
	  // setBpd = getBpd + (getRip()-time_left);
	  //time_left =  getBpd() - TA - TRDSUND - Buffer - (getRwd()+getPul()*getPri()+getTroInPriUnits()*getPri()+getRL()+getND());
	  
	//	}


	if (time_left < Uvar(0,"s"))
	  {
	    cout<<"ComputeRad: bpd rwd pri+tor rl nd time_left"<<getBpd()<<" "<<getRwd()<<" "<<
	      getPul()*getPri()+getTroInPriUnits()*getPri()<<" "<<getRL()<<" "
		<<getND()<<time_left<<endl;

	    ErrorMessage("Ieb.cpp-computeRad:no time for Rad").throwMe();
	  }
	unsigned int rad = (unsigned int) (time_left/getRip()).getInUnits("");
	if(rad==0) ErrorMessage("Ieb.cpp::computeRad(): bpd is too short").throwMe();
	setRad(rad);
	
	Uvar pad =getBpd()
	  -( TA + getRwd() + double(getPul())*getPri() 
	     + double(getTroInPriUnits())*getPri()
	     +TRDSUND+getND() + getRL() + getRip()*getRad());
	if (pad < Uvar(0,"s")) ErrorMessage("Ieb.cpp::computeRad:no padding time ").throwMe();
	/*	 }
	   else
	     {//calibration mode
	//fill up the burst period
	time_left -= getRL()+getND();
	unsigned int rad = (unsigned int) (time_left/getRip()).getInUnits("");
	if(rad==0) ErrorMessage("Ieb.cpp::computeRad(): bpd is too short").throwMe();
	setRad(rad);
	Uvar pad = getBpd()
	  -( TA+ getND() + getRL() + getRip()*getRad());
	if (pad < Uvar(0,"s")) ErrorMessage("Ieb.cpp::computeRad:no padding time ").throwMe();
	}*/
  }


void Ieb::computeRadandBpd(const unsigned int& bursts_in_flight)
  {
    Uvar bpd = Uvar(1,"s");//default
    Uvar orig_bpd = bpd ;
    //need data rate
    if(!data_rate_set_) ErrorMessage("Ieb.cpp::computeRadandBpd: data rate is not set").throwMe();
    
    //active mode
    if (s_mod_e_ < 4 || (s_mod_e_ > 7 && s_mod_e_< 12))
      {
	if(!pul_set_) 
	  {
	    ErrorMessage("Ieb.cpp::computeBpd: Number of pulses are not set").throwMe();
	  }  
	if(!pri_set_) 
	  {
	    ErrorMessage("Ieb.cpp::computeBpd: No Pri is set").throwMe();
	  }
	if(!dutycycle_set_) 
	  {
	    ErrorMessage("Ieb.cpp::computeBpd: No dutycycle is set").throwMe();
	  }
	if(!csq_set_)
	  {
	    ErrorMessage("Ieb.cpp::computeBpd: no csq set ").throwMe();
	  }
	Uvar bpd_based_on_burst_dutycycle = double (f_pul_e_)* 
	  double(s_csq_e_ +1) * s_csd_ /(6.95/100.0);
	unsigned int Nword = getNumberOfWords();
	Uvar bpd_based_on_data_rate = Uvar(double(Nword * 16)/data_rate_,"s");
	bpd = bpd_based_on_burst_dutycycle;
	if (bpd < bpd_based_on_data_rate  ) {
	  bpd = bpd_based_on_data_rate;}

       // handle multiply burst_in_flight here....
	orig_bpd = bpd;
//cout <<"ghtemp:Nword"<< Nword<<" :bpd(datarate)"<<bpd_based_on_data_rate<<":bpd(dytycycle)"<<bpd_based_on_burst_dutycycle<<":bpd"<<bpd<<endl;
	if (bursts_in_flight > 1 && s_csr_e_ != 1) RwdAdjustForNBursts (bursts_in_flight, bpd);

      }
    else if (s_mod_e_ == 4)
      {//radiometer mode
	unsigned int Nword = getNumberOfWords();
	bpd = Uvar(double(Nword * 16)/data_rate_,"s");
      }
    else
      {
	ErrorMessage("Ieb.cpp::computeBpd(data rate):invalid radar mode to compute burst period").throwMe();
      }

   
    if(!pri_set_) ErrorMessage("Ieb.cpp::computeRad:no pri set").throwMe();
    if(!tro_set_) ErrorMessage("Ieb.cpp::computeRad:no tro set").throwMe();
    
    Uvar TA = Uvar(1.0*mstos,"s")+getPri();
    Uvar TRDSUND = Uvar(2.0*mstos,"s");
    Uvar Buffer =Uvar(2*mstos,"s");
    Uvar time_left = bpd - TA - TRDSUND - Buffer;
    unsigned int rad=1;//default

    if (s_csr_e_ == 0 || s_csr_e_ == 8 || s_csr_e_ == 1)  // should this really just be for all cases???? GH 09012004
      {//either normal operation or auto gain control    
	time_left -= getRwd()+getPul()*getPri()+getTroInPriUnits()*getPri()+getRL()+getND();
	//if (time_left < Uvar(0,"s")) {
	if (time_left < getRip() ) {
	  //need to increase bpd
	  cout << "Bpd not large enough for a RAD sample.  Bpd:" << bpd << ":timeleft:" << time_left;
	  bpd += time_left*-1+getRip();// increase by time_left plus 1 RIP
	  cout << " increase BPD (scape+RIP): " << bpd << endl;
	  time_left = bpd - TA - TRDSUND - Buffer - (getRwd()+getPul()*getPri()+getTroInPriUnits()*getPri()+getRL()+getND());
	  cout << "New time_left (for RAD) = " << time_left << endl;
	}

	//	if (time_left < Uvar(0,"s"))
	  //	  {
	    //	    cout<<" bpd rwd pri+tor rl nd time_left "<<bpd<<" "<<getRwd()<<" "<<
	      //	      getPul()*getPri()+getTroInPriUnits()*getPri()<<" "<<getRL()<<" "
	      //		<<getND()<<time_left<<endl;
	    //	    ErrorMessage("Ieb.cpp-computeRad&Bpd:no time for Rad").throwMe();
	    //	  }

	rad = (unsigned int) (time_left/getRip()).getInUnits("");
	if (rad ==0)
	  {//need to increase bpd
	    cout << "Warning(Ieb.cpp): Bpd not large enough for a RAD sample.  Bpd:" << bpd ;
	    bpd += getRip() + time_left;
	    rad++;
	    cout << " increase BPD(for RIP): " << bpd << endl;
	  }
	
	Uvar pad = bpd
	  -( TA + getRwd() + double(getPul())*getPri() 
	     + double(getTroInPriUnits())*getPri()
	     +TRDSUND+getND() + getRL() + getRip()*double(rad));
	if (pad < Uvar(0,"s")) ErrorMessage("Ieb.cpp::computeRad&Bpd:no padding time ").throwMe();
      }
    else
      {//calibration mode
	//fill up the burst period
	time_left -= getRL()+getND();
	rad = (unsigned int) (time_left/getRip()).getInUnits("");
	if(rad==0)
	  {
	    bpd += getRip();
	    rad++;
	  }
	Uvar pad = bpd
	  -( TA+ getND() + getRL() + getRip()*double(rad));
	if (pad < Uvar(0,"s")) ErrorMessage("Ieb.cpp::computeRad&Bpd:no padding time ").throwMe();
      }

    if (orig_bpd != bpd && bursts_in_flight > 1) {
      cout << "Ghdebug: bpd was re-adjusted, recalculate RWD for multiple bursts in flight" << endl;
      ErrorMessage("Ieb.cpp-cannot double adjust.  Need to find valid Data rate and pulse in air calc...").throwMe();
	//RwdAdjustForNBursts (bursts_in_flight, bpd);
    }

	
    setBpd(bpd);
    setRad(rad);
  }

//------------------------------------------
//{get,set, compute} chirp step duration, no set method
//-----------------------------------------
void Ieb::computeWaveform()
  {
    unsigned int csd_value = 5;
  while (!csq_set_){
  //--------------------
  //csd set
  //-----------------
  if(csd_set_)
    {
      ErrorMessage("Ieb::computeCsd: chirp step duration is already set").throwMe();
    }
  s_csd_e_ = csd_value;
  s_csd_ = double(s_csd_e_) * 4.0 / Uvar(30*1e6,"Hz");
  //csd_set_ = true;

  //-------------
  //csq set
  //-------------
  if(csq_set_)
    {
      ErrorMessage("Ieb::computeCsq: csq is already set").throwMe();
    }
  if (!pri_set_|| !dutycycle_set_)
    {
      ErrorMessage("Ieb::computeCsq:Either pri or dutycyle is not set").throwMe();
    }
  
  tau_p_ = getPri() * dutycycle_;//pri*dutycycle
  int csq_tmp = int((tau_p_ /s_csd_).getInUnits(""));
  csq_tmp -=1;

  if (csq_tmp < 2)
    {
      cout << "Ieb.cpp:CSQ  must be greater than 2.  CSQ is: " << csq_tmp << endl; 
      csq_tmp = 2;
      cout << "Ieb.cpp:CSQ Forcing CSQ to minimum of 2" << endl;
      //ErrorMessage e("Ieb::setCsq: Chirp step quantity is out of range: 2..750");
      //e.throwMe();
	
          s_csq_e_ = (unsigned int) csq_tmp;
          tau_p_ = s_csd_ * double(s_csq_e_+1);
          csq_set_ = true;
          csd_set_ = true;
        }

  if (csq_tmp > 750)
    {
      cout << "Ieb.cpp:CSQ  must be less than 750.  CSQ is: " << csq_tmp << endl; 
      csd_value++;
      cout << "Ieb.cpp:computeWaveform: Increase CSD to try and get longer duration..CSD to " << csd_value << endl;
      //ErrorMessage e("Ieb::setCsq: Chirp step quantity is out of range: 2..750");
      //e.throwMe();
	} else { // assign csq parameter and set to true and get out of while loop
	  s_csq_e_ = (unsigned int) csq_tmp;
	  tau_p_ = s_csd_ * double(s_csq_e_+1);
	  csq_set_ = true;
	  csd_set_ = true;
	}

  } // End while !csq_set

  //----
  //Cfs
  //----
  if(chirp_frequency_step_set_)
    {
      ErrorMessage("Ieb::computeCfs(): cfs is already set").throwMe();
    }
  if(!mod_set_)
    {
      ErrorMessage("Ieb::computeCfs: no radar mode is set to determine chirp bandwidth").throwMe();
    }
    
  s_cfs_ = BR_ / double(s_csq_e_); 
  //Uvar cfs_min= Uvar(30.0*1e6,"Hz")/pow(2,24);
  //Uvar cfs_max = cfs_min * pow(2,16); 
  Uvar cfs_min= Uvar(0,"Hz");
  Uvar cfs_max = Uvar(30.0*1e6,"Hz")/pow(2,24) * pow(2,16); 
  if (s_cfs_ < cfs_min)
    {
      ErrorMessage("Chirp frequency step <  min value of (1.788*0) Hz").throwMe();
    }
  if (s_cfs_ > cfs_max)
    {
      ErrorMessage("Chirp frequency step  > 117,187 KHz").throwMe();
    }
  s_cfs_e_ = (unsigned short) (s_cfs_/(Uvar(30*1e6,"Hz")/pow(2,24))).getInUnits(""); 
  s_cfs_ = Uvar(30*1e6,"Hz")/pow(2,24)*double(s_cfs_e_);
  BR_ = s_cfs_ * double(s_csq_e_);
  chirp_frequency_step_set_ = true;
  }

void Ieb::constructWaveform(const Uvar& csd, const unsigned int& csq, const Uvar& cfs)
  {
  //--------------------
  //csd set
  //-----------------
  if(csd_set_)
    {
      ErrorMessage("Ieb::computeCsd: chirp step duration is already set").throwMe();
    }
  double x = (csd/(4.0/Uvar(30*1e6,"Hz"))).getInUnits("");
  s_csd_e_ = (unsigned int) round_double(x);
  if(s_csd_e_ != 5) ErrorMessage("Ieb.constructWaveform(): encoded csd is not 5").throwMe();
  s_csd_ = double(s_csd_e_) * 4.0 / Uvar(30*1e6,"Hz");
  csd_set_ = true;

  //-------------
  //csq set
  //-------------
  if(csq_set_)
    {
      ErrorMessage("Ieb::computeCsq: csq is already set").throwMe();
    }
  if (csq < 2 || csq > 750)
    {
      ErrorMessage e("Ieb::setCsq: Chirp step quantity is out of range: 2..750");
      e.throwMe();
    }
  s_csq_e_ = csq;
  tau_p_ = double(csq+1)*csd;
  csq_set_ = true;

  //----
  //cfs
  //----
  if(chirp_frequency_step_set_)
    {
      ErrorMessage("Ieb::computeCfs(): cfs is already set").throwMe();
    }
  if(!mod_set_)
    {
      ErrorMessage("Ieb::computeCfs: no radar mode is set to determine chirp bandwidth").throwMe();
    }

  s_cfs_ =cfs; 
  //Uvar cfs_min= Uvar(30.0*1e6,"Hz")/pow(2,24);
  //Uvar cfs_max = cfs_min * pow(2,16); 
  Uvar cfs_min= Uvar(0,"Hz");
  Uvar cfs_max = Uvar(30.0*1e6,"Hz")/pow(2,24) * pow(2,16); 
  if (s_cfs_ < cfs_min)
    {
      ErrorMessage("Chirp frequency step <  min value of (1.788*0) Hz").throwMe();
    }
  if (s_cfs_ > cfs_max)
    {
      ErrorMessage("Chirp frequency step  > 117,187 KHz").throwMe();
    }
  s_cfs_e_ = (unsigned short) round_double((s_cfs_/(Uvar(30*1e6,"Hz")/pow(2,24))).getInUnits("")); 
  s_cfs_ = Uvar(30*1e6,"Hz")/pow(2,24)*double(s_cfs_e_);
  BR_ = s_cfs_ * double(s_csq_e_);
  chirp_frequency_step_set_ = true;
  }


//-------------------------------------------------
//GH:void changeWaveform used to merge slowfields instructions
//-------------------------------------------------
void Ieb::changeWaveform( const unsigned int& newcsq, const Uvar& newcfs )
{
  // not set flag checking.  Actually the flag should be set.  We're overwritting regardless
 //-------------
  //csq set
  //-------------
 
  if (newcsq < 2 || newcsq > 750)
    {
      ErrorMessage e("Ieb::setCsq: Chirp step quantity is out of range: 2..750." );
      e.throwMe();
    }
  // accepting new csq
  s_csq_e_ = newcsq;
  tau_p_ = double(newcsq+1)*s_csd_ ;
  csq_set_ = true;

  //----
  //cfs set
  //----

  s_cfs_ = newcfs;
  //Uvar cfs_min= Uvar(30.0*1e6,"Hz")/pow(2,24);
  //Uvar cfs_max = cfs_min * pow(2,16); 
  Uvar cfs_min= Uvar(0,"Hz");
  Uvar cfs_max = Uvar(30.0*1e6,"Hz")/pow(2,24) * pow(2,16); 
  if (s_cfs_ < cfs_min)
    {
      ErrorMessage("Chirp frequency step <  min value of (1.788*0) Hz").throwMe();
    }
  if (s_cfs_ > cfs_max)
    {
      ErrorMessage("Chirp frequency step  > 117,187 KHz").throwMe();
    }
  s_cfs_e_ = (unsigned short) round_double((s_cfs_/(Uvar(30*1e6,"Hz")/pow(2,24))).getInUnits(""));
  s_cfs_ = Uvar(30*1e6,"Hz")/pow(2,24)*double(s_cfs_e_);
  BR_ = s_cfs_ * double(s_csq_e_);
  chirp_frequency_step_set_ = true;
}

//--------------------
//{get} Chirp step quantity
//------------------
Uvar Ieb::getCsd() const
  {
    if(!csd_set_)
      {
	ErrorMessage("Ieb::getCsd: chirp step duration is not set").throwMe();
      }
  return(s_csd_);
  }

unsigned int Ieb::getCsdEncodedValue() const
{
  if(!csd_set_)
    {
      ErrorMessage("Ieb::getCsdEncodedValue: chirp step duration is not set").throwMe();
    }
  return(s_csd_e_);
}


//---------------------------------
//{get} chirp step quantity
//--------------------------------
unsigned int Ieb::getCsq() const
  {
    if (!csq_set_)
      {
	ErrorMessage("Ieb::getCsq: csq is not set");
      }
    return(s_csq_e_);
  }


//----------------------------------------
//{get} chirp frequency step size
//----------------------------------------
Uvar Ieb::getCfs() const
  {
    if(!chirp_frequency_step_set_)
      {
	ErrorMessage("Ieb::getCfs:cfs is not set").throwMe();
      }
    return(s_cfs_);
  }

unsigned int Ieb::getCfsEncodedValue() const
  {
    if(!chirp_frequency_step_set_)
      {
	ErrorMessage("Ieb::getCfsEncodedValue():cfs is not set").throwMe();
      }
    return(s_cfs_e_);
  }

//------------------------------------------------
//GH:resyncTimeTfi: method to reset TFIs after time has shifted
//-------------------------------------------------

void Ieb::resyncTimeTfi (const Time& iebstarttime, const Time& ilxexecutiontime)

 {
   //cout << "HERE(resync)" << iebstarttime <<":"<< t_ <<":"<< s_tfi_ ;
   // change original time(t_) and tfi to new input value

   s_tfi_ = ilxexecutiontime  - iebstarttime  ;
   t_ = ilxexecutiontime;

   //cout << "(resync TFI: " <<  s_tfi_ << endl;
	
   // set all TFIs to same value
   s_tfi_e_ = (unsigned short)(s_tfi_.getInUnits("s"));
   s_tfi_ = Uvar(double(s_tfi_e_),"s");

   f_tfi_e_ = (unsigned short)(s_tfi_.getInUnits("s"));
   f_tfi_ = Uvar(double(f_tfi_e_),"s");
   p_tfi_e_ = (unsigned short)(s_tfi_.getInUnits("s"));
   p_tfi_ = Uvar(double(p_tfi_e_),"s");
   t_tfi_e_ = (unsigned short)(s_tfi_.getInUnits("s"));
   t_tfi_ = Uvar(double(t_tfi_e_),"s");

    s_tfi_set_ = true;
    f_tfi_set_ = true;
    p_tfi_set_ = true;
    t_tfi_set_ = true;
    time_set_ = true;  
 }

//----------------------------------
//get pri*dutycycle
//---------------------------------
Uvar Ieb::getTaup() const
  {
  if(!csq_set_)
    {
      ErrorMessage("Ieb::getTaup: Chirp step quantity is not set").throwMe();
    }
  return(tau_p_);
  }

//--------------------------
//get chirp bandwidth
//---------------------------
Uvar Ieb::getBr() const
  {
  if(!chirp_frequency_step_set_)
    {
      ErrorMessage("Ieb::getBr:chirp frequency step is not set").throwMe();
    }
  return(BR_);
  }

//-------------------------------------------
//get active mode time offset from the time of burst
//-------------------------------------------
Uvar Ieb::getOffsetTimefromBurst() const
  {
    if(!pri_set_){
      ErrorMessage("Ieb::getOffsetTimefromBurst(): no pri set ").throwMe();
    }
    if(!rwd_set_){
      ErrorMessage("Ieb::getOffsetTimefromBurst(): no rwd set ").throwMe();
    }
    return( Uvar(1*mstos,"s")+ f_pri_+f_rwd_);
  }
//--------------------------------
//{get,set} power tfi
//--------------------------------
Uvar Ieb::getPowerTfi() const
  {
  if(!p_tfi_set_) 
    ErrorMessage("Ieb.cpp:getPowerTfi: no tfi is set").throwMe();
  return(p_tfi_);
  }
void Ieb::setPowerTfi(const Uvar& tfi)
{
  if(p_tfi_set_) 
    ErrorMessage("Ieb.cpp:setPowerTfi: tfi is set").throwMe();
  int i = (unsigned int) tfi.getInUnits("s");
  if (i < 0 || i > 65536) ErrorMessage("Ieb.cpp:setPowerTfi: tfi is out of range").throwMe();
  p_tfi_e_ = (unsigned int) i;
  p_tfi_ = Uvar(p_tfi_e_,"s");
  p_tfi_set_ = true;
}

//------------------------------
//{get} power typ
//----------------------------- 
unsigned int Ieb::getPowerTyp() const
{
  return(p_typ_e_);
}

//--------------------------
//{get,set} power pmd
//---------------------------   
unsigned int Ieb::getPmd() const
{
  return(pmd_e_);
}
void Ieb::setPmd(const unsigned int& pmd)
{
  if((pmd > 5 && pmd < 17)|| (pmd>21)) ErrorMessage("Ieb.cpp:setPmd: pmd is out of range valid range: 1 2 3 4 5 17 18 19 20 21").throwMe();
  pmd_e_ = pmd;
}


//------------------------------------
//{get,set} TncTfi
//------------------------------------
Uvar Ieb::getTncTfi() const
  {
    if(!t_tfi_set_)
      ErrorMessage("Ieb.cpp::getTncTfi: no tfi value is set").throwMe();
    return(t_tfi_);
  }

void Ieb::setTncTfi(const Uvar& tfi)
  {
    if(t_tfi_set_)
      ErrorMessage("Ieb.cpp::setTncTfi: tfi is set").throwMe();
    int  a = int(tfi.getInUnits("s"));
    if (a < 0) ErrorMessage("Ieb.cpp:tnc_tfi value is negative").throwMe();
    t_tfi_e_ = (unsigned int) a;
    t_tfi_ = Uvar(t_tfi_e_,"s");
    t_tfi_set_ = true;
  }

//------------------------------
//{get} tncparamters (TYP, TNC, TCA, TCB, TCC)
//----------------------------- 
unsigned int Ieb::getTncTyp() const
{
  return(t_typ_e_);
}

unsigned int Ieb::getTnc() const
{
  return(tnc_e_);
}

unsigned short Ieb::getTca() const
{
  return(tca_);
}

unsigned short Ieb::getTcb() const
{
  return(tcb_);
}

unsigned short Ieb::getTcc() const
{
  return(tcc_);
}

//---------------------------
//{set} Resistive Load and Noise Diode Integration Radiometer Calibration
// with AutoRad On
//-----------------------------
void Ieb::setRL_ND_AutoRad(const Uvar& rl, const Uvar& max_rl,
			   const Uvar& nd, const Uvar& max_nd)
  {//AutoRad is "On"
    tnc_e_ = 64;//40 in hex (see blue book 11-15)
    if (rl < Uvar(5.0/1000.0,"s") || rl > Uvar(75.0/1000.0,"s"))
      {
	cout<<"input rip "<<rl<<endl;
	ErrorMessage("Ieb.cpp::input resistive load integeration period is out of range").throwMe();
      }
    if (max_rl < Uvar(5.0/1000.0,"s") || max_rl > Uvar(75.0/1000.0,"s"))
      {
	cout<<"input rip "<<max_rl<<endl;
	ErrorMessage("Ieb.cpp::max input resistive load integeration period is out of range").throwMe();
      }

    if (nd < Uvar(2.0/1000.0,"s") || nd > Uvar(16.0/1000.0,"s"))
      {
	cout<<"input nd "<<nd<<endl;
	ErrorMessage("Ieb.cpp::input noise diode integration period is out of range").throwMe();
      }

    if (max_nd < Uvar(2.0/1000.0,"s") || max_nd > Uvar(16.0/1000.0,"s"))
      {
	cout<<"input nd "<<max_nd<<endl;
	ErrorMessage("Ieb.cpp::max input noise diode integration period is out of range").throwMe();
      }

    rl_e_ = (unsigned int) round_double((rl*1000.0/5.0).getInUnits("s"));
    rl_  = Uvar(rl_e_ * 5.0 * mstos,"s");
    nd_e_ = (unsigned int)  round_double((nd*1000.0).getInUnits("s"));
    nd_e_ = nd_e_ -1;
    nd_ = Uvar((nd_e_ + 1)*mstos,"s");
    bitset(tca_,0,3,nd_e_);
    bitset(tca_,4,7,rl_e_);
    bitset(tca_,8,15,0);

    bitset(tcb_,0,15,101);//65 in hex (see blue book 11-15)-AutoRad On

    //store max value in parameter C, and set rl and nd to be max value
    rl_e_ = (unsigned int) round_double((max_rl*1000.0/5.0).getInUnits("s"));
    rl_  = Uvar(rl_e_ * 5.0 * mstos,"s");
    nd_e_ = (unsigned int)  round_double((max_nd*1000.0).getInUnits("s"));
    nd_e_ = nd_e_ -1;
    nd_ = Uvar((nd_e_ + 1)*mstos,"s");
    bitset(tcc_,0,3,nd_e_);
    bitset(tcc_,4,7,rl_e_);
    bitset(tcc_,8,15,0);

    tnc_set_ = true;
    tca_set_ = true;
    tcb_set_ = true;
    tcc_set_ = true;
  }
//---------------------------
//{set} Resistive Load and Noise Diode Integration Radiometer Calibration
// with AutoRad Off
//-----------------------------
void Ieb::setRL_ND(const Uvar& rl, const Uvar& nd)
  {//AutoRad is "Off"
    tnc_e_ = 64;//40 in hex (see blue book 11-15)
    if (rl < Uvar(5.0/1000.0,"s") || rl > Uvar(75.0/1000.0,"s"))
      {
	cout<<"Warning:::input rip "<<rl<<endl;
	//ErrorMessage("Ieb.cpp::input resistive load RIP value is out of range").throwMe();
      }
    if (nd < Uvar(2.0/1000.0,"s") || nd > Uvar(16.0/1000.0,"s"))
      {
	cout<<"Warning::input nd "<< nd<<endl;
	//ErrorMessage("Ieb.cpp::input noise diode RIP value is out of range").throwMe();
      }
    rl_e_ = (unsigned int) round_double((rl.getInUnits("s")*1000.0/5.0));
    rl_  = Uvar(rl_e_ * 5.0 * mstos,"s");
    nd_e_ = (unsigned int) round_double((nd.getInUnits("s")*1000.0-1.0));
    nd_ = Uvar((nd_e_ + 1)*mstos,"s");
    bitset(tca_,0,3,nd_e_);
    bitset(tca_,4,7,rl_e_);
    bitset(tca_,8,15,0);

    bitset(tcb_,0,15,169);//a9_h in hex (see blue book 11-15)-- AutoRad off

    //auto rad off, max value is same as input value
    bitset(tcc_,0,3,nd_e_);
    bitset(tcc_,4,7,rl_e_);
    bitset(tcc_,8,15,0);

    tnc_set_ = true;
    tca_set_ = true;
    tcb_set_ = true;
    tcc_set_ = true;
  }

//------------------------------
//{get} RL , ND
//------------------------------
Uvar Ieb::getRL() const
  {
    if(!(tnc_set_ && tnc_e_ == 64))
      ErrorMessage("Ieb.cpp::getRL: no tnc setting for RL ").throwMe();
    return(rl_);
  }
Uvar Ieb::getND() const
  {
    if(!(tnc_set_ && tnc_e_ == 64))
      ErrorMessage("Ieb.cpp::getND: no tnc setting for ND ").throwMe();
    return(nd_);
  }
//----------------------------
//generate slow field parameters
//----------------------------
void Ieb::encodeSlowFastfield() 
  {
  if(!f_tfi_set_){ErrorMessage("Ieb::generateSlowFastfield:f_tfi  is not set").throwMe();}
  //if(!fin_set_){ErrorMessage("Ieb::generateSlowFastfield:fin  is not set").throwMe();}
  if(!pul_set_){ErrorMessage("Ieb::generateSlowFastfield:pul  is not set").throwMe();}
  if(!bpd_set_){ErrorMessage("Ieb::generateSlowFastfield:bpd is not set").throwMe();}
  if(!rwd_set_){ErrorMessage("Ieb::generateSlowFastfield:rwd  is not set").throwMe();}
  if(!pri_set_){ErrorMessage("Ieb::generateSlowFastfield: pri is not set").throwMe();}
  if(!chirp_start_frequency_set_){ErrorMessage("Ieb::generateSlowFastfield:chirp start frequency is not set").throwMe();}
  if(!s_tfi_set_){ErrorMessage("Ieb::generateSlowFastfield:s_tfi is not set").throwMe();}
  if(!dtn_set_){ErrorMessage("Ieb::generateSlowFastfield:dtn is not set").throwMe();}
  //if(!sin_set_){ErrorMessage("Ieb::generateSlowFastfield:sin is not set").throwMe();}
  if(!mod_set_){ErrorMessage("Ieb::generateSlowFastfield:mod is not set").throwMe();}
  if(!csr_set_){ErrorMessage("Ieb::generateSlowFastfield:csr is not set").throwMe();}
  if(!adc_set_){ErrorMessage("Ieb::generateSlowFastfield:adc is not set").throwMe();}
  if(!rcv_set_){ErrorMessage("Ieb::generateSlowFastfield:rcv is not set").throwMe();}
  if(!tro_set_){ErrorMessage("Ieb::generateSlowFastfield:tro is not set").throwMe();}
  if(!baq_set_){ErrorMessage("Ieb::generateSlowFastfield:baq is not set").throwMe();}
  if(!bem_set_){ErrorMessage("Ieb::generateSlowFastfield:bem is not set").throwMe();}
  if(!at1N1_set_){ErrorMessage("Ieb::generateSlowFastfield:at1N1 is not set").throwMe();}
  if(!at1N2_set_){ErrorMessage("Ieb::generateSlowFastfield:at1N2 is not set").throwMe();}
  if(!at1N3_set_){ErrorMessage("Ieb::generateSlowFastfield:at1N3 is not set").throwMe();}
  if(!at3N1_set_){ErrorMessage("Ieb::generateSlowFastfield:at3N1 is not set").throwMe();}
  if(!at3N2_set_){ErrorMessage("Ieb::generateSlowFastfield:at3N2 is not set").throwMe();}
  if(!at3N3_set_){ErrorMessage("Ieb::generateSlowFastfield:at3N3 is not set").throwMe();}
  if(!at4N1_set_){ErrorMessage("Ieb::generateSlowFastfield:at4N1 is not set").throwMe();}
  if(!at4N2_set_){ErrorMessage("Ieb::generateSlowFastfield:at4N2 is not set").throwMe();}
  if(!at4N3_set_){ErrorMessage("Ieb::generateSlowFastfield:at4N3 is not set").throwMe();}
  if(!rip_set_){ErrorMessage("Ieb::generateSlowFastfield:rip is not set").throwMe();}
  if(!rad_set_){ErrorMessage("Ieb::generateSlowFastfield:rad is not set").throwMe();}
  if(!csd_set_){ErrorMessage("Ieb::generateSlowFastfield:csd is not set").throwMe();}
  if(!csq_set_){ErrorMessage("Ieb::generateSlowFastfield:csq is not set").throwMe();}
  if(!chirp_frequency_step_set_){ErrorMessage("Ieb::generateSlowFastfield:chirp frequency step is not set").throwMe();}
  
  //-----------------------------------
  //Check constraints involving more than two parameters
  //----------------------------------
  // RW max check
  // -----------
  int s_tro_pri = (int) s_tro_e_;
  if (s_tro_pri > 7) {s_tro_pri = s_tro_pri - 16;}
  //s_tro_ needs pri value 
  int rw = int(f_pul_e_) + s_tro_pri;
  if (rw >= 256) 
      ErrorMessage("Ieb::encodeSlowFastfield: rw= pul + tro is out of range > 255").throwMe();
    
  // -----------
  // GH_ILXCHECKER: BAQ to PRI check
  // -----------
  // pri check:  pri min = 192 for BAQ 0,1,4,5,7.  pri_min = 96 for BAQ 2,3,6
  // found that this is also clock dependant.  if 10 Mhz 192 is ok.  tested doen to 140
   unsigned int pri_min_check = 192;

  // for BAQ 5 (str8) this should still be 192 min.  However it is currently being allowed.
  if (s_baq_e_ == 5) pri_min_check = 20;

  if (f_pri_e_ < 192 && s_baq_e_ == 5) {
    cout << "Warning(ILXCHECK):  This will cause FLTSW errors.  However we are currently allowing this mode." << endl;
    cout << " \tThis should only be allowed for Distant high PRI scat.  Is it?" << endl;
  }


  if (s_baq_e_ == 2 || s_baq_e_ == 6) pri_min_check = 96; // GH edit 03/24/2006. let PRI go lower for baq = 3
  if ( s_baq_e_ == 3 ) pri_min_check = 20; 
//20110708 gh adjust if 10MHz clock
  if (s_adc_e_ == 3){
	 pri_min_check = 96;  // tested down to 80, still works, but we will set to 96 for now (for no good reason)
	cout << "Warning:  Accepting a low pri value when clock is 10 MHz: adc=" << s_adc_e_ << "  PRImin check =" << pri_min_check << endl;
   }
// 20110708


  if( f_pri_e_ < pri_min_check) {
    cout << "ERROR: PRI must be greater or equal to "<<pri_min_check<<"! Failed....Halting!  was:" << f_pri_e_<< ".  For BAQ = " << s_baq_e_ <<  endl;
    ErrorMessage("Ieb::f_pri_e_ : Less than minimum (192 or 96). out of range").throwMe();
    }


  //------------------------------------
  //data buffer < Nmax_buffer
  //-------------------------------------
  unsigned int Ndata = getRadarDataPoints();
  if (Ndata > Nmax_buffer) 
    {
      cout<<"mode "<< getModName()<<endl;
      cout<<"pulse "<<getPul()<<endl;
      cout<<"adc "<<getAdc()<<endl;
      cout<<"pri "<<getPri()<<endl;
      cout<<"encoded pri "<<getPriEncodedValue()<<endl;
      cout<<"baq "<<getBaq()<<endl;
      ErrorMessage("Ieb.cpp::encoding: data buffer is larger than Nmax_buffer= 32Kbytes").throwMe();
    }
  unsigned int no_of_beam_used = 0;
  for (unsigned int i = 0; i < 5;++i)
    {
      unsigned int mark=bitget(s_bem_e_,i,i);
      if (mark == 1) {no_of_beam_used ++;}
    }
  //if nnumber of beam used is not zero
  if(no_of_beam_used !=0 && f_bii_e_ % no_of_beam_used !=0)
    {
      cout<<"bad command:number of bursts is not an integer multiple of the number of beams used"<<endl;
      cout<<"beam mask "<<s_bem_e_<<endl;
      cout<<"number of burst"<<f_bii_e_<<endl;
      cout<<"number of beam used "<<no_of_beam_used<<endl;
      cout<<endl;
    }
  //it is possible to have no_of_beam_used=0 because of calibration mode
  if (s_adc_ <= 2.0 * s_rcv_)
    {
      ErrorMessage("Ieb::generateSlowFastfield: adc clock speed is not roughly twice the receiver bandwidth").throwMe();
    }
  if (BR_ > s_rcv_)
    {
	cout << "Error(nonhalt)Ieb::generateSlowFastfield: chirp bandwidth is larger than receiver bandwidth" << endl;

	//ErrorMessage("Ieb::generateSlowFastfield: chirp bandwidth is larger than receiver bandwidth").throwMe();
    }
  
  double HPA_duty_cycle = (s_csd_ * double(s_csq_e_ + 1)/f_pri_).getInUnits("");
  if (HPA_duty_cycle >= 75.0/100.0)
    {
      ErrorMessage("Ieb::generateSlowFastfield: HPA duty cycle is larger than 75%").throwMe();
    }
  //double ESS_duty_cycle = (getPulseTransmitTime()*dutycycle_/getBpd()).getInUnits("");
  double ESS_duty_cycle = (double (f_pul_e_)* double(s_csq_e_ +1) * 
			  s_csd_ /f_bpd_).getInUnits("");
  if (ESS_duty_cycle > 7.0/100.0)
    {
      cout<<"----------------------------------------------------"<<endl;
      cout<<"ERROR: Ieb::generateSlowFastfield: ESS duty cycle is larger than 7%"<<endl;
      cout<<"-------------------------------------------------------"<<endl;
      
      cout<<ESS_duty_cycle<<endl;
    }

 

  //-----
  //Slow field
  //----- 
  if (s_tfi_e_ >= (unsigned int) pow(2,16)) {ErrorMessage("Ieb::s_tfi_e_: out of range").throwMe();} 
  if (s_typ_e_ !=3) {ErrorMessage("Ieb::s_type_e_: out of range").throwMe();}
  if (s_dtn_e_ >= (unsigned int) pow(2,8)){ErrorMessage("Ieb::s_dtn_e_:out of range").throwMe();}

  bitset(slowfield_[0],0,15,s_tfi_e_);
  bitset(slowfield_[1],0,7,s_dtn_e_);  
  bitset(slowfield_[1],8,9,s_typ_e_);
  
  if(s_sin_e_ >= (unsigned int) pow(2,8)) {ErrorMessage("Ieb::s_sin_e_:out of range").throwMe();}
  if(s_mod_e_ >= (unsigned int) pow(2,4)) {ErrorMessage("Ieb::s_mod_e_:out of range").throwMe();}
  if(s_csr_e_ >= (unsigned int) pow(2,4)) {ErrorMessage("Ieb::s_csr_e_:out of range").throwMe();}
  bitset(slowfield_[2],8,15,s_sin_e_);
  bitset(slowfield_[2],4,7,s_mod_e_);  
  bitset(slowfield_[2],0,3,s_csr_e_);
 
 

  //------------------
  //bitset for adc,rcv,tro,baq,bem
  //-----------------
  if(s_bem_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_beam_e_:out of range").throwMe();}
  if(s_baq_e_ >= (unsigned int) pow(2,3)) {ErrorMessage("Ieb::s_baq_e_:out of range").throwMe();}
  if(s_tro_e_ >= (unsigned int) pow(2,4)) {ErrorMessage("Ieb::s_tro_e_:out of range").throwMe();}
  if(s_rcv_e_ >= (unsigned int) pow(2,2)) {ErrorMessage("Ieb::s_rcv_e_:out of range").throwMe();}
  if(s_adc_e_ >= (unsigned int) pow(2,2)) {ErrorMessage("Ieb::s_adc_e_:out of range").throwMe();}
  bitset(slowfield_[3],0,4,s_bem_e_);
  bitset(slowfield_[3],5,7,s_baq_e_);
  bitset(slowfield_[3],8,11,s_tro_e_);
  bitset(slowfield_[3],12,13,s_rcv_e_);
  bitset(slowfield_[3],14,15,s_adc_e_);
  
  //---------------------
  //attenuation setting: refer attenuation map
  //-------------------- 
  s_at1N1_e_ = s_at1N1_dB_;
  s_at1N2_e_ = s_at1N2_dB_;
  s_at1N3_e_ = s_at1N3_dB_/4;
  if ( s_at1N3_dB_% 4 != 0) ErrorMessage("attenuation 3 of att1 is not a mulitple of 4").throwMe();
  
  s_at3N1_e_ = s_at3N1_dB_;
  s_at3N2_e_ = s_at3N2_dB_;
  s_at3N3_e_ = s_at3N3_dB_/4;
  if ( s_at3N3_dB_% 4 != 0) ErrorMessage("attenuation 3 of att3 is not a mulitple of 4").throwMe();

  s_at4N1_e_ = s_at4N1_dB_;
  s_at4N2_e_ = s_at4N2_dB_;
  s_at4N3_e_ = s_at4N3_dB_/4;
  if ( s_at4N3_dB_% 4 != 0) ErrorMessage("attenuation 3 of att4 is not a mulitple of 4").throwMe();
 


  if(s_at1N3_e_ >= (unsigned int) pow(2,2)) {ErrorMessage("Ieb::s_at1N3_e:out of range").throwMe();}
  if(s_at1N2_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at1N2_e:out of range").throwMe();}
  if(s_at1N1_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at1N1_e:out of range").throwMe();}
  
  if(s_at3N3_e_ >= (unsigned int) pow(2,2)) {ErrorMessage("Ieb::s_at3N3_e:out of range").throwMe();}
  if(s_at3N2_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at3N2_e:out of range").throwMe();}
  if(s_at3N1_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at3N1_e:out of range").throwMe();}
  
  if(s_at4N3_e_ >= (unsigned int) pow(2,2)) {ErrorMessage("Ieb::s_at4N3_e:out of range").throwMe();}
  if(s_at4N2_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at4N2_e:out of range").throwMe();}
  if(s_at4N1_e_ >= (unsigned int) pow(2,5)) {ErrorMessage("Ieb::s_at4N1_e:out of range").throwMe();}
  if(s_rip_e_ >=(unsigned int) pow(2,4)){ErrorMessage("Ieb::s_rip_e_:out of range").throwMe();}
  
  bitset(slowfield_[4],0,1,s_at1N3_e_);
  bitset(slowfield_[4],2,6,s_at1N2_e_);
  bitset(slowfield_[4],7,11,s_at1N1_e_);

  bitset(slowfield_[5],0,1,s_at3N3_e_);
  bitset(slowfield_[5],2,6,s_at3N2_e_);
  bitset(slowfield_[5],7,11,s_at3N1_e_);

  bitset(slowfield_[6],0,3,s_rip_e_);
  bitset(slowfield_[6],4,5,s_at4N3_e_);
  bitset(slowfield_[6],6,10,s_at4N2_e_);
  bitset(slowfield_[6],11,15,s_at4N1_e_);

  //-------------------
  //radiometer integration period from rip: see BB12-36
  //--------------------  
  if(s_rad_e_ >= (unsigned int) pow(2,8)) {ErrorMessage("Ieb::s_rad_e_:out of range").throwMe();}
  if(s_csr_e_ >=(unsigned int) pow(2,8)) {ErrorMessage("Ieb::s_csr_e_:out of range").throwMe();}
  if(s_csq_e_ >=(unsigned int) pow(2,12)) {ErrorMessage("Ieb::s_csq_e_:out of range").throwMe();}

  bitset(slowfield_[7],8,15,s_rad_e_);
  bitset(slowfield_[7],0,7,s_csd_e_);
  bitset(slowfield_[8],0,11,s_csq_e_);

  //----------------
  //chirp start frequency
  //---------------
  if(s_cfs_e_ >= (unsigned int) pow(2,16)) {ErrorMessage("Ieb::s_cfs_e_:out of range").throwMe();}
  bitset(slowfield_[9],0,15,s_cfs_e_);


  if(f_tfi_e_>= (unsigned int) pow(2,16)) {ErrorMessage("Ieb::f_tfi_e_ :out of range").throwMe();}
  if(f_fin_e_ >= (unsigned int) pow(2,8)) {ErrorMessage("Ieb::f_fin_e_ :out of range").throwMe();}
  if(f_typ_e_!= 2) {ErrorMessage("Ieb::f_typ_e_ :out of range").throwMe();}
  if(f_bii_e_ >= (unsigned int) pow(2,8)) {ErrorMessage("Ieb::f_bii_e_ :out of range").throwMe();}
  if(f_pul_e_ >= (unsigned int) pow(2,8)) {ErrorMessage("Ieb::f_pul_e_ :out of range").throwMe();}
  if(f_bpd_e_ >= (unsigned int) pow(2,12)) {ErrorMessage("Ieb::f_bdp_e_ :out of range").throwMe();}
  if(f_rwd_e_ >= (unsigned int) pow(2,10)) {cout << "Ieb: f_rwd_e_ = "<< f_rwd_e_ << endl; ErrorMessage("Ieb::f_rwd_e_ :out of range").throwMe();}
  if(f_pri_e_ >= (unsigned int) pow(2,10)) {ErrorMessage("Ieb::f_pri_e_ :out of range").throwMe();}
  if(f_csf_e_ >= (unsigned int) pow(2,16)) {ErrorMessage("Ieb::f_csf_e_ :out of range").throwMe();}

  //fast field
  bitset(fastfield_[0],0,15,f_tfi_e_);
  bitset(fastfield_[1],0,7,f_fin_e_);
  bitset(fastfield_[1],8,9,f_typ_e_);
  bitset(fastfield_[2],0,7,f_pul_e_);
  bitset(fastfield_[2],8,15,f_bii_e_);
  bitset(fastfield_[3],0,11,f_bpd_e_);
  bitset(fastfield_[4],0,9,f_rwd_e_);
  bitset(fastfield_[5],0,9,f_pri_e_);
  bitset(fastfield_[6],0,15,f_csf_e_); 
  encode_slowfast_complete_ = true;

  }

//---------------------------------
//encode power : gh add to split power tnc
//--------------------------------
void Ieb::encodePowerfield()
{
  //-------------------------------------
  //power field
  //-------------------------------------
   
  if(!p_tfi_set_) ErrorMessage("Ieb.cpp:encodepower:no p_tfi set").throwMe(); 
  bitset(powerfield_[0],0,15,p_tfi_e_);
  bitset(powerfield_[1],0,7,pmd_e_);
  bitset(powerfield_[1],8,9,p_typ_e_);

  encode_power_complete_=true;
  // so it is backward compatible with powertnc together
  if (encode_tnc_complete_ == true ) encode_powertnc_complete_ = true;
}

//---------------------------------
//encode  tnc // GH added to split power tnc
//--------------------------------
void Ieb::encodeTncfield()
{
  //-------------------------------------
  //tnc field
  //-------------------------------------
 
  if(!t_tfi_set_) ErrorMessage("Ieb.cpp:encodetnc: no tnc_tfi set").throwMe();
  bitset(tncfield_[0],0,15,t_tfi_e_);
  bitset(tncfield_[1],0,7,tnc_e_);
  bitset(tncfield_[1],8,9,t_typ_e_);
  tncfield_[2] = tca_;
  tncfield_[3] = tcb_;
  tncfield_[4] = tcc_; 

  encode_tnc_complete_=true;
  // so it is backward compatible with powertnc together
  if (encode_power_complete_ == true ) encode_powertnc_complete_ = true;
}


//---------------------------------
//encode power tnc 
//--------------------------------
void Ieb::encodePowerTncfield()
{
  //-------------------------------------
  //power field
  //-------------------------------------
  encodePowerfield();
  encodeTncfield();

  /*
  if(!p_tfi_set_) ErrorMessage("Ieb.cpp:encodepowertnc:no p_tfi set").throwMe(); 
  bitset(powerfield_[0],0,15,p_tfi_e_);
  bitset(powerfield_[1],0,7,pmd_e_);
  bitset(powerfield_[1],8,9,p_typ_e_);
    
  encode_power_complete_=true;
 
  if(!t_tfi_set_) ErrorMessage("Ieb.cpp:encodepowertnc: no tnc_tfi set").throwMe();
  bitset(tncfield_[0],0,15,t_tfi_e_);
  bitset(tncfield_[1],0,7,tnc_e_);
  bitset(tncfield_[1],8,9,t_typ_e_);
  tncfield_[2] = tca_;
  tncfield_[3] = tcb_;
  tncfield_[4] = tcc_;

  encode_tnc_complete_=true;
  */
  if ( encode_power_complete_ && encode_tnc_complete_ ) encode_powertnc_complete_=true;
}


//---------------------------------
//decode(): it will check whether slowfield and fastfield data are loaded
// then start deocde
//---------------------------------
void Ieb::decodeSlowFastfield() 
  {
    if (decode_slowfast_complete_) return; 
    if(!(encode_slowfast_complete_ || slowfastfield_loaded_))
      {
	ErrorMessage("Ieb::decodeSlowFastfield: incomplete encoding or no slow/fast field is loaded").throwMe();
      }
  s_tfi_e_ = bitget(slowfield_[0],0,15);
  s_tfi_ = Uvar(double(s_tfi_e_),"s");

  s_dtn_e_ = bitget(slowfield_[1],0,7);
  unsigned short slow_instruction_type = bitget(slowfield_[1],8,9);

  // cout << "GH:s_typ should equal " << slow_instruction_type << ", but it is: " << s_typ_e_ ;
  //cout <<"stfi"<< s_tfi_e_<<"dtn:"<< s_dtn_e_<< endl;

  if (slow_instruction_type !=s_typ_e_) 
    {
      cout << "s_typ should equal " << slow_instruction_type << ", but it is: " << s_typ_e_ << endl;
      ErrorMessage("Ieb.cpp: SLOW fiel TYP error").throwMe();
    }
  s_csr_e_ = bitget(slowfield_[2],0,3);
  s_mod_e_ = bitget(slowfield_[2],4,7);
  switch(s_mod_e_)
    {
    case 0:
      radar_mode_ ="ALTL";
      break;
    case 1:
      radar_mode_ ="ALTH";
      break;
    case 2:
      radar_mode_="SARL";
      break;
    case 3:
      radar_mode_="SARH";
      break;
    case 4:
      radar_mode_="Radiometer only";
      break;
    case 5:
      radar_mode_="IGO calibration";
      break;
    case 6:
      radar_mode_="Earth Viewing Calibration";
      break;
    case 7:
      radar_mode_ = "Bi-static operation";
      break;
    case 8:
      radar_mode_="ALTL_AGC";
      break;
    case 9:
      radar_mode_="ALTH_AGC";
      break;
    case 10:
      radar_mode_="SARL_AGC";
      break;
    case 11:
      radar_mode_="SARH_AGC";
      break;
    default:
      ErrorMessage("Ieb::decode: decoded radar mode does not have an appropriate value").throwMe();
    }

  //if (s_mod_e_ < 4 && s_csr_e_!=0)
  //{
  //  cout<<"Ieb::decode: calibration source value "<<s_csr_e_<<endl;
  //  cout<<"Ieb::decode: radar mode "<<radar_mode_<<endl;
  //}
  //if ((s_mod_e_ >7 && s_mod_e_<12) && s_csr_e_!=8)
  //{
  //  cout<<"Ieb::decode: calibration source value "<<s_csr_e_<<endl;
  //  cout<<"Ieb::decode: radar mode "<<radar_mode_<<endl;
  //}

  s_sin_e_ = bitget(slowfield_[2],8,15);// slow field instruction number

  s_bem_e_ = bitget(slowfield_[3],0,4);

  s_baq_e_ = bitget(slowfield_[3],5,7);

  s_tro_e_ = bitget(slowfield_[3],8,11) ;
  int s_tro_pri = (int) s_tro_e_;
  if (s_tro_pri > 7) {s_tro_pri = s_tro_pri - 16;} 

  double rc_bw_map[4]={0.117, 0.468,0.935,4.675};
  s_rcv_= Uvar(rc_bw_map[bitget(slowfield_[3],12,13)]*1e6,"Hz");
  s_rcv_e_= bitget(slowfield_[3],12,13);

  double adc_map[4]={0.25, 1.0, 2.0, 10.0};
  s_adc_ = Uvar(adc_map[bitget(slowfield_[3],14,15)]*1e6,"Hz");
  s_adc_e_=bitget(slowfield_[3],14,15);

  s_at1N3_e_ = bitget(slowfield_[4],0,1);
  s_at1N3_dB_ = s_at1N3_e_*4;
  s_at1N2_e_= bitget(slowfield_[4],2,6);
  s_at1N2_dB_ = s_at1N2_e_;
  s_at1N1_e_ = bitget(slowfield_[4],7,11);
  s_at1N1_dB_ = s_at1N1_e_;
  
  if ( (int(s_at1N1_dB_) - int(s_at1N2_dB_)) < 0 ||
       (int(s_at1N1_dB_) - int(s_at1N2_dB_))>1) 
    cout << "ERROR:Ieb.cpp: attenuator1 decoding error.  Not following Map!" << endl;
    //ErrorMessage("attenuation1 decoding error").throwMe();

  s_at3N3_e_ = bitget(slowfield_[5],0,1);
  s_at3N3_dB_ = s_at3N3_e_ *4 ;
  s_at3N2_e_ = bitget(slowfield_[5],2,6);
  s_at3N2_dB_ = s_at3N2_e_;
  s_at3N1_e_ = bitget(slowfield_[5],7,11);
  s_at3N1_dB_ = s_at3N1_e_;
  s_rip_ = Uvar(bitget(slowfield_[6],0,3)*5*mstos,"s");
  s_rip_e_=bitget(slowfield_[6],0,3);

  if ( (int(s_at3N1_dB_) - int(s_at3N2_dB_)) < 0 ||
       (int(s_at3N1_dB_) - int(s_at3N2_dB_))>1)
    cout << "ERROR:Ieb.cpp: attenuator3 decoding error.  Not following Map!" << endl;
    //ErrorMessage("attenuation3 decoding error").throwMe();


  s_at4N3_e_ = bitget(slowfield_[6],4,5);
  s_at4N3_dB_ = s_at4N3_e_ *4;
  s_at4N2_e_ =bitget(slowfield_[6],6,10);
  s_at4N2_dB_ = s_at4N2_e_;
  s_at4N1_e_ = bitget(slowfield_[6],11,15);
  s_at4N1_dB_ = s_at4N1_e_;

  if ( (int(s_at4N1_dB_) - int(s_at4N2_dB_)) < 0 ||
       (int(s_at4N1_dB_) - int(s_at4N2_dB_))>1)
    cout << "ERROR:Ieb.cpp: attenuator4 decoding error.  Not following Map!" << endl;
  //ErrorMessage("attenuation4 decoding error").throwMe();

  s_csd_ = Uvar(bitget(slowfield_[7],0,7)*133.3333333*nstos,"s");
  s_csd_e_ = bitget(slowfield_[7],0,7);

  s_rad_e_ = bitget(slowfield_[7],8,15);
  
  s_csq_e_ = bitget(slowfield_[8],0,11);
  

  s_cfs_ = Uvar(bitget(slowfield_[9],0,15)*1.7881393,"Hz");
  s_cfs_e_ = bitget(slowfield_[9],0,15);

  tau_p_= Uvar( (s_csq_e_ + 1) * bitget(slowfield_[7],0,7)*133.333333*nstos ,"s");
  BR_ = s_cfs_ * double (s_csq_e_);
  //
  //fast field
  //
  f_tfi_ = Uvar(bitget(fastfield_[0],0,15),"s");
  f_tfi_e_ = bitget(fastfield_[0],0,15);

  f_fin_e_ = bitget(fastfield_[1],0,7);
  unsigned short fast_instruction_type =bitget(fastfield_[1],8,9); 
  if (fast_instruction_type != f_typ_e_) 
    {
      ErrorMessage("Fast field TYP error").throwMe();
    }
  f_pul_e_ = bitget(fastfield_[2],0,7);
 
  int rw = int(f_pul_e_) + s_tro_pri;
  if (rw >= 256) 
      ErrorMessage("Ieb::encodeSlowFastfield: rw= pul + tro is out of range > 255").throwMe();

  f_bii_e_ = bitget(fastfield_[2],8,15);  
  f_bpd_ = Uvar(bitget(fastfield_[3],0,11)*mstos,"s");
  f_bpd_e_ =bitget(fastfield_[3],0,11);

  
  //--------------------------------------------------------
  //total number of data per pri
  //According to blue book 13-16,
  //PRI is tightly associated with adc clock speed
  // fadc : 0.25 MHz, 1.0 MHz, 2.0 MHz, 10 MHz
  // pri_number: 10 bit information 2 - 1023
  //---------------------------------------------


  double pri_map[4] = {4.0, 1.0, 0.5, 0.1};
  f_pri_e_ = bitget(fastfield_[5],0,9);
  if (s_adc_e_ <3)
    {
    if (f_pri_e_ %2 != 0) 
      {
      ErrorMessage("Ieb::decodePRI Number is not an even number");
      }
    f_pri_ = Uvar(pri_map[s_adc_e_] * double(f_pri_e_)*1e-6, "s");
    s_tro_ = f_pri_ * double(s_tro_pri);
    }
  else if (s_adc_e_ == 3)
    {
    f_pri_ = Uvar(pri_map[s_adc_e_] * double(f_pri_e_)*10.0*1e-6, "s");
    s_tro_ = f_pri_ * double(s_tro_pri);
    }
  else
    {
    ErrorMessage("fadc number is not correct").throwMe();
    }

  prf_ = 1/f_pri_;
  dutycycle_ = (tau_p_/f_pri_).getInUnits("");
  if (dutycycle_ > 0.75) 
    {
      cout<<"Ieb.cpp::decodeSlowFastfield: duty cycle is higher than 75% "<<endl;
      cout<<"pri "<< f_pri_<<endl;
      cout<<"tau_p "<<tau_p_<<endl;
    }

  double ESS_dutycycle = (double (f_pul_e_)* double(s_csq_e_ +1) * 
			  s_csd_ /f_bpd_).getInUnits("");
  if (ESS_dutycycle > 0.07)
    {
      cout<<"Ieb.cpp::decodeSlowFastfield: long term duty cycle is higher than 7%"<<endl;
      cout<<"pri "<<f_pri_<<endl;
      cout<<"pul "<<f_pul_e_<<endl;
      cout<<"bdp "<<f_bpd_<<endl;
      cout<<"ESS duty cycle"<<ESS_dutycycle<<endl;
    }

  f_rwd_ = bitget(fastfield_[4],0,9) * f_pri_;
  f_rwd_e_ =bitget(fastfield_[4],0,9);

  f_csf_= Uvar(bitget(fastfield_[6],0,15)*457.764,"Hz");
  f_csf_e_ = bitget(fastfield_[6],0,15);



  dutycycle_set_ = true;
  f_tfi_set_ = true;
  //fin_set_ = true;
  pul_set_ = true;
  bpd_set_ = true;
  rwd_set_ = true;
  pri_set_ = true;
  chirp_start_frequency_set_ = true;  
  s_tfi_set_ = true;
  dtn_set_ = true;
  //sin_set_ = true;
  mod_set_ = true;
  csr_set_ = true;
  adc_set_ = true;
  rcv_set_ = true;
  tro_set_ = true;
  baq_set_ = true;
  bem_set_ = true;
  at1N1_set_ = true;
  at1N2_set_ = true;
  at1N3_set_ = true;
  at3N1_set_ = true;
  at3N2_set_ = true;
  at3N3_set_ = true;
  at4N1_set_ = true;
  at4N2_set_ = true;
  at4N3_set_ = true;
  rip_set_ = true;
  rad_set_ = true;
  csd_set_ = true;
  csq_set_ = true;
  chirp_frequency_step_set_ = true; 
  decode_slowfast_complete_ = true;
  }

//------------------------------
//decode power // GH; added to slip power tnc
//-----------------------------
void Ieb::decodePowerfield()
{
  if(!(encode_power_complete_ || (power_loaded_)))
      {
	ErrorMessage("Ieb::decodePowerfield: incomplete encoding or no power field is loaded").throwMe();
      }
  //power mode
  //-------------------------------------
  //power field
  //-------------------------------------
  p_tfi_e_= bitget(powerfield_[0],0,15);
  p_tfi_ = Uvar(p_tfi_e_,"s");      
  pmd_e_ = bitget(powerfield_[1],0,7);
  unsigned int p_typ = bitget(powerfield_[1],8,9);
  if (p_typ != p_typ_e_)
    {
      ErrorMessage("power instruction type mismatch").throwMe();
    }     
  p_tfi_set_ = true;
    
  decode_power_complete_ = true;
  // so it is backward compatible with powertnc together
  if (decode_tnc_complete_ == true ) decode_powertnc_complete_ = true;
}
//------------------------------
//decode power tnc // GH: added to split power/ tnc
//-----------------------------
void Ieb::decodeTncfield()
{
  if(!(encode_tnc_complete_ || (tnc_loaded_)))
      {
	ErrorMessage("Ieb::decodeTncfield: incomplete encoding or no tnc field is loaded").throwMe();
      }
  //tnc mode
  //-------------------------------------
  //tnc field
  //-------------------------------------
    
  t_tfi_e_=  bitget(tncfield_[0],0,15);
  t_tfi_set_ = true;
  t_tfi_ = Uvar(t_tfi_e_,"s");
  tnc_e_ =bitget(tncfield_[1],0,7);
  unsigned int t_typ_e=bitget(tncfield_[1],8,9);
  if (t_typ_e !=t_typ_e_ ) {
    ErrorMessage("tnc type instruction does not match").throwMe();
  }

  tca_= tncfield_[2] ;
  tcb_=tncfield_[3] ;
  tcc_=tncfield_[4] ;
  if (tnc_e_ == 64)
    {//resitive load and noise diode intergation time is available
      nd_e_ = bitget(tca_,0,3);
      nd_ = Uvar(double(nd_e_+1)*mstos,"s");
      rl_e_ = bitget(tca_,4,7);
      rl_ = Uvar( double(rl_e_) * 5.0 * mstos,"s");
    } 
  tnc_set_ = true;
  tca_set_ = true;
  tcb_set_ = true;
  tcc_set_ = true;

  decode_tnc_complete_ = true;
  // so it is backward compatible with powertnc together
  if (decode_power_complete_ == true ) decode_powertnc_complete_ = true;
}


//------------------------------
//decode power tnc
//-----------------------------
void Ieb::decodePowerTncfield()
{

  decodePowerfield();
  decodeTncfield();
  /*
  if(!(encode_powertnc_complete_ || (power_loaded_&&tnc_loaded_)))
      {
	ErrorMessage("Ieb::decodePowerTncfield: incomplete encodiding or no power/tnc field is loaded").throwMe();
      }
  //power mode
  //-------------------------------------
  //power field
  //-------------------------------------
  p_tfi_e_= bitget(powerfield_[0],0,15);
  p_tfi_ = Uvar(p_tfi_e_,"s");      
  pmd_e_ = bitget(powerfield_[1],0,7);
  unsigned int p_typ = bitget(powerfield_[1],8,9);
  if (p_typ != p_typ_e_) ErrorMessage("ieb.cpp:decodepowertnd:power typ set is wrong").throwMe();      
  p_tfi_set_ = true;
  decode_power_complete_ = true; // GH
    
  t_tfi_e_=  bitget(tncfield_[0],0,15);
  t_tfi_ = Uvar(t_tfi_e_,"s");
  tnc_e_ =bitget(tncfield_[1],0,7);
  unsigned int t_typ_e=bitget(tncfield_[1],8,9);
  if (t_typ_e !=t_typ_e_ ) ErrorMessage("Ieb.cpp::decodepowertnc: tnc type is not valid").throwMe();
  tca_= tncfield_[2] ;
  tcb_=tncfield_[3] ;
  tcc_=tncfield_[4] ;
  if (tnc_e_ == 64)
    {//resitive load and noise diode intergation time is available
      nd_e_ = bitget(tca_,0,3);
      nd_ = Uvar(double(nd_e_+1)*mstos,"s");
      rl_e_ = bitget(tca_,4,7);
      rl_ = Uvar( double(rl_e_) * 5.0 * mstos,"s");
    } 
  tnc_set_ = true;
  tca_set_ = true;
  tcb_set_ = true;
  tcc_set_ = true;
  decode_tnc_complete_ = true; // GH
  */
  if ( decode_power_complete_ && decode_tnc_complete_ ) decode_powertnc_complete_ = true;
}

//---------------------------------
//loadSlowfield(const std::vector<unsigned short> slow): load slow field data
// from telemetry file
//--------------------------------
void Ieb::loadSlowFastfield(const std::vector<unsigned short>& slow, 
			     const std::vector<unsigned short>& fast)
  {
  if (Nslowfield !=slow.size())
    {
    ErrorMessage("Ieb::Slow field size is set incorrectlyout of range").throwMe();
    }
  if (Nfastfield !=fast.size())
    {
    ErrorMessage("Ieb::Fast field size is set incorrectlyout of range").throwMe();
    }
  slowfield_ = slow; 
  fastfield_ = fast;
  slowfastfield_loaded_ = true;
  }

//--------------------------------
//load power
//---------------------------------
void Ieb::loadPowerfield(const sdata& powerfield)
			  
{
   if (Npowerfield !=powerfield.size())
    {
    ErrorMessage("Ieb::power field size is set incorrectlyout of range").throwMe();
    }
  powerfield_ = powerfield; 
  power_loaded_ = true;
}
//--------------------------------
//load TNC
//---------------------------------
void Ieb::loadTncfield( const sdata& tncfield ) 
  {
  if (Ntncfield !=tncfield.size())
    {
    ErrorMessage("Ieb::Tnc field size is set incorrectlyout of range").throwMe();
    }
  tncfield_ = tncfield;
  tnc_loaded_ = true;
  }
//---------------------------
//{get} slow/fast field
//--------------------------
vector<unsigned short> Ieb::getSlowfield() const
{ 
  return(slowfield_);
}

vector <unsigned short> Ieb::getFastfield() const
{
  return(fastfield_);
}


//---------------------------
//{get} power/tnc field
//--------------------------
vector<unsigned short> Ieb::getTncfield() const
{
  return(tncfield_);
}

vector<unsigned short> Ieb::getPowerfield() const
{
  return(powerfield_);
}


//--------------------------
//encodecomplete(): return(encode_set_)
//--------------------------
bool Ieb::encodeSlowFastfieldComplete()
  {
  return(encode_slowfast_complete_);
  }
bool Ieb::encodePowerTncfieldComplete()
  {
  return(encode_powertnc_complete_);
  }
//GH: added to split pwer/tnc
bool Ieb::encodePowerfieldComplete()
  {
  return(encode_power_complete_);
  }
bool Ieb::encodeTncfieldComplete()
  {
  return(encode_tnc_complete_);
  }


//----------------------------------
//decodecomplete(): return(decode_set_)
//----------------------------------
bool Ieb::decodeSlowFastfieldComplete()
  {
  return(decode_slowfast_complete_);
  }
bool Ieb::decodePowerTncfieldComplete()
  {
  return(decode_powertnc_complete_);
  }
//GH: added to split pwer/tnc
bool Ieb::decodePowerfieldComplete()
  {
  return(decode_power_complete_);
  }
bool Ieb::decodeTncfieldComplete()
  {
  return(decode_tnc_complete_);
  }


//------------------------------
//compute radar data points per pri
//------------------------------
unsigned int Ieb::getPointsInsidePRI()
{
  if(!pri_set_) 
    ErrorMessage("Ieb.cpp:getPointsInsidePri(): no pri set").throwMe();
  unsigned int Npoints = 0;
  if ((s_mod_e_ > 4 && s_mod_e_ < 8) || (s_mod_e_>11)) 
    ErrorMessage("Ieb.cpp:getPointsInsidePRI: invalid radar mode").throwMe();
  
  if (s_mod_e_ == 4) return(Npoints);//no active data
  else if (s_mod_e_ < 4 || (s_mod_e_ >7 && s_mod_e_ < 12))
    {
      if (s_adc_e_ <3)
	{
	  if (f_pri_e_%2 != 0) 
	    {
	      ErrorMessage("Ieb::decodePRI Number is not an even number");
	    }
	  Npoints= f_pri_e_ ;
	}
      else if(s_adc_e_ == 3)
	{
	  Npoints = f_pri_e_ * 10 ;
	}
      else
	{
	  ErrorMessage("Ieb.cpp::getPointsInsidePRI: no valid adc clock").throwMe();
	}
    }
  else
    {
      ErrorMessage("ieb.cpp:getRadarDataPoints(): no valid radar mode").throwMe();
    }
  return(Npoints);
}

//------------------------------
//compute radar data points per bpd
//------------------------------
unsigned int Ieb::getRadarDataPoints()
  { 
    if (!pul_set_)
      ErrorMessage("Ieb.cpp:getRadarDataPoints():no pul set").throwMe();
    if (!tro_set_)
      ErrorMessage("Ieb.cpp:getRadarDataPoints():no TRO set").throwMe();
    if ((s_mod_e_ > 4 && s_mod_e_ < 8) || (s_mod_e_>11)) 
      ErrorMessage("Ieb.cpp:getRadarDataPoints: invalid radar mode").throwMe();
    if (s_mod_e_ == 4) return(0);//no active data for radiometer

    int rw = (int)getPul();
    rw += getTroInPriUnits();
    if (rw < 0) 
      ErrorMessage("Ieb.cpp:getRadarDataPoints: negative window").throwMe();   
    return(getPointsInsidePRI() * (unsigned int) rw);
  }
     
     
//---------------------------
//compute number of words to store radar buffer
//---------------------------
unsigned int Ieb::getNumberOfWords()
  {
  if(!baq_set_) 
    ErrorMessage("Ieb.cpp::getNumberOfWords: no baq set").throwMe();
  if ((s_mod_e_ > 4 && s_mod_e_ < 8) || (s_mod_e_>11)) 
    ErrorMessage("Ieb.cpp:getNumberofWords: invalid radar mode").throwMe();

 
  if (s_mod_e_ == 4)
    {//radiometer mode
      return(Nheader+Nfooter);
    }

  unsigned int Ndat_bpd = getRadarDataPoints();
  unsigned int Nword = 0;
  if(s_baq_e_ == 0)
    {//8-2 bit, 8 data into one word
    unsigned int word = Ndat_bpd/8;
    int remain = int(Ndat_bpd)  - int(word * 8);
    if (remain == 0) Nword= word; 
    else if (remain > 0 && remain < 8) Nword= word +1;
    else ErrorMessage("ieb.encodeSlowFastfield: word counting error").throwMe();
    }
  else if (s_baq_e_ ==1 || s_baq_e_ == 2)
    {//8-1 16 data into one word
    unsigned int word = Ndat_bpd/16;
    int remain = int(Ndat_bpd) - int(word * 16);
    if (remain == 0) Nword = word;
    else if (remain > 0 && remain < 15) Nword = word +1;
    else ErrorMessage("ieb.encodeSlowFastfield: word counting error").throwMe();
    }
  else if (s_baq_e_ == 3)
    {//new scatt_compression: 32 bits (2 words) for offset and 
     //each data point is stored into each word
     //each data takes up 1 word+2words for dc offset
      Nword = 2 + getPointsInsidePRI();
    }
  else if (s_baq_e_ == 4)
    {//8-4 MSB, 4 data into one word
    unsigned int word = Ndat_bpd/4;
    int remain = int(Ndat_bpd) - int(word*4);
    if (remain == 0) Nword = word;
    else if (remain > 0 && remain < 3) Nword= word +1;
    else ErrorMessage("ieb.encodeSlowFastfield: word counting error").throwMe();
    }
  else if (s_baq_e_ == 5)
    {//8 bit straight
    unsigned int word= Ndat_bpd/2;
    int remain = int(Ndat_bpd) - int(word*2);
    if (remain !=0) ErrorMessage("ieb.encodeSlowFast:word counting error").throwMe(); 
    Nword = word;
    }
  else if (s_baq_e_ == 6 || s_baq_e_==7)
    {//either 8-4 BAQ
    unsigned int word = Ndat_bpd/4;
    int remain = int(Ndat_bpd) - int(word*4);
    if (remain == 0) Nword = word;
    else if (remain > 0 && remain < 3) Nword = word +1;
    else ErrorMessage("ieb.encodeSlowFastfield: word counting error").throwMe();
    }
  else
    {
      cout<<"used baq"<<s_baq_e_<<endl;
      ErrorMessage("ieb.encodeSlowFast: invalid baq mode to compute number or words").throwMe();
    }
  Nword += Nheader + Nfooter;
  return(Nword);
   
  }

//----------------------------------
//{get} datarate
//----------------------------------
double Ieb::getDataRate()
  {
    if(!bpd_set_) ErrorMessage("Ieb.cpp::getDataRate:no burst period is set").throwMe();
    return(double(getNumberOfWords())*16.0/getBpd().getInUnits("s"));
  }

void Ieb::setMaxDataRate(const double& data_rate)
  {
    if(data_rate_set_)
      {
	ErrorMessage("Ieb.cpp::data rate is set").throwMe();
      }
    data_rate_ = data_rate;
    data_rate_set_ = true;
  }

//----------------------------------------
//adjust RWD if there are multiple bursts in flight
// 
//-----------------------------------------
void Ieb::RwdAdjustForNBursts (const unsigned int& bursts_in_flight, Uvar& bpd)
 {
   Uvar tempRwd = getRwd(); //Uvar(0,"s") ;   
   Uvar tempBpd = bpd;
   
   tempBpd = Uvar(ceil(bpd.getInUnits("ms"))/1000,"s");

   // handle multiply burst_in_flight here....
   cout << "RWDADJUST(BurstInFlight): roundtrip(inpri): " << " RWD(sec): " << getRwd() << " :RWD(pri)):"<< getRwdInPriUnits()<<":"<<f_rwd_e_ << " BPD: "  << bpd.getInUnits("s") << " BPD(+round): "  << tempBpd.getInUnits("s") <<  "nbursts: " << bursts_in_flight << endl;

   tempRwd -=  ( tempBpd.getInUnits("s") * (bursts_in_flight -1 )) ;
   cout << "tempRwd: " <<tempRwd << endl;
   rwd_set_ = false; setRwd (tempRwd);

   // handle multiply burst_in_flight here....
   cout << "RWDADJUST(newNumbers): roundtrip(inpri): " << " RWD(sec): " << getRwd() << " :RWD(pri)):"<< getRwdInPriUnits()<<":"<<f_rwd_e_ << " BPD(+round): "  << tempBpd.getInUnits("s") <<  "nbursts: " << bursts_in_flight << endl;

	// add check to RWD calculation
	if ( getRwdInPriUnits() < (getPul()+1) ) // min rwd = pul+1
	{
		cout << "ERROR: Ieb.cpp: While adjusting RWD for BIF.  RWD is now inside TX window.  RWD must be at least PUL+1.  RWD:PUL=" << getRwdInPriUnits()<<":"<< getPul()<< endl;
        ErrorMessage("Ieb.cpp::While adjusting RWD for BIF.  RWD is now inside TX window.  RWD must be at least PUL+1.").throwMe();

	}
		    
}

//----------------------
//clear method: set all the flags to false
//----------------------
void Ieb::clear()
  {
   time_set_=false; 
   dutycycle_set_=false;
   
   f_tfi_set_=false;
   //fin_set_=false;
   pul_set_=false;
   bpd_set_=false;
   rwd_set_=false;
   pri_set_=false;
   chirp_start_frequency_set_=false;
   
   s_tfi_set_=false;
   dtn_set_=false;
   //sin_set_=false;
   mod_set_=false;
   csr_set_=false;
   adc_set_=false;
   rcv_set_=false;
   tro_set_=false;
   baq_set_=false;
   bem_set_=false;
   at1N1_set_=false;
   at1N2_set_=false;
   at1N3_set_=false;
   at3N1_set_=false;
   at3N2_set_=false;
   at3N3_set_=false;
   at4N1_set_=false;
   at4N2_set_=false;
   at4N3_set_=false;
   rip_set_=false;
   rad_set_=false;
   csd_set_=false;
   csq_set_=false;
   chirp_frequency_step_set_=false;
   
   encode_slowfast_complete_=false;
   decode_slowfast_complete_=false;
   encode_powertnc_complete_=false;
   decode_powertnc_complete_=false;
   slowfastfield_loaded_=false;
   power_loaded_=false;
   tnc_loaded_=false;
   p_tfi_set_ = false;
   t_tfi_set_ = false;
   tnc_set_=false;
   tca_set_=false;
   tcb_set_=false;
   tcc_set_=false;
   data_rate_set_ = false;
  }


//----------------------------
//self test
//---------------------------
bool Ieb::selfTest()
  {
    try
      {
	Time t("2005-301T04:30:37.840");
	Ieb ieb(t);
	ieb.setFastTfi(Uvar(0,"s"));
	ieb.setSlowTfi(Uvar(0,"s"));
	ieb.setFin(0);
	ieb.setSin(0);
	ieb.setDtn(0);
	ieb.computeMod(sarh_adc,0,0,31, 100);
	ieb.setBii(5);
	
	
	Uvar prf = Uvar(5000,"Hz");
	Uvar pri = 1/prf;
	double dutycycle = 0.7;
	ieb.setPri(pri);
	ieb.setDutycycle(dutycycle);
	ieb.computeWaveform();
	unsigned int csq = int(((pri * dutycycle)/(5.0*4.0/Uvar(30e6,"Hz"))).getInUnits(""))-1;
	if (csq != ieb.getCsq()) ErrorMessage("ieb.selfTest: csq encoding error").throwMe();
 
	
	if (ieb.getPriEncodedValue() != 400) ErrorMessage("Ieb:selfTest():Pri Encoding Error ").throwMe();   
	//Uvar transmit_time = Uvar(2000,"km")/speed_light;
	//ieb.setPulseTransmitTime(transmit_time-ieb.getPri());
	//unsigned int Np = (unsigned int) (transmit_time/pri).getInUnits("");
	//if (ieb.getPul() != 33) ErrorMessage("Ieb:selfTest: encoding Pul error").throwMe();
	//Uvar tro_time = pri * (-3.0);
	//ieb.setTroTime(tro_time);
	//if( ieb.getTroInPriUnits() != -3) ErrorMessage("ieb.selfTest(): tro encoding error").throwMe();
	
	//ieb.setRwd(pri * double(Np));
	//ieb.computeChirpStartFrequency(Uvar(0,"Hz"));
	//Uvar bpd = Uvar(60*mstos,"s");
	//ieb.setBpd(bpd);
	//if (ieb.getBpd() != Uvar(80*mstos,"s")) ErrorMessage("ieb.selfTest:bpd min_bpd error").throwMe();
	//ieb.setRip(Uvar(25*mstos,"s"));
	//ieb.computeRad();
	//ieb.setAttenuation(15,20,15);
	//ieb.encodeSlowFastfield();
	//ieb.decodeSlowFastfield();
      }
    catch(const ErrorMessage& e)
      {
	std::cout << e.msg << std::endl;
      }
    catch(...)
      {
	std::cout << "Exception caught" << std::endl;
      }    
    return(true);
  }




//-----------------------------
//construct attenuation map for a given total value of attenuation
//----------------------------
void computeAttenuationMap(unsigned int& at1, unsigned int& at2, unsigned int& at3,const unsigned int& total_att)
  {
    int att = int(total_att);
    if (att == 0)
      {
	at1 = 0;
	at2 = 0;
	at3 = 0;
      }
    else if (att > 0 && att <=62)
      {
	at3 = 0;
	int a = int((att+1)/2);
	at1 = (unsigned int) a;
	a = att - a;
	if (a < 0) ErrorMessage("at2 is out of range").throwMe();;
	at2 = (unsigned int) (a);
      }
    else if (att > 62 && att <= 74)
      {
	int a = att - 63;
	if (a < 0) ErrorMessage("at3 is out of range");
	int b = int(a/4);
	at3 = (unsigned int)((b+1) * 4);
	
	a = att - int(at3);//recycle a
	if (a < 0) ErrorMessage("at1 + at2 is out of range").throwMe();
	a = int((a+1)/2);
	at1 = (unsigned int) a;
	a = att - int(at3) - int(at1);
	if (a < 0) ErrorMessage("at2 is out of range").throwMe();
	at2 = (unsigned int) a;
      }
    else
      {
	ErrorMessage("total attenuation value is out of range").throwMe();
      }
  }

  //----------------------------------------
  //convert beam pattern to binary number
  // 11111 to 31
  //-----------------------------------------
 void setBeamMask(unsigned int& bem, const unsigned int& beam_mask)
  {
    //check size
    bem = 0;
    int beam_set = (int) beam_mask;
   
    if(beam_set > 11111) ErrorMessage("Error in reading beam mask").throwMe();

    if(beam_set >= 10000)
      {
	bem += 16;
	if (int(beam_set/10000) >= 2)
	  ErrorMessage("non valid beam option").throwMe();
	beam_set -= 10000;

      }
    if (beam_set >= 1000)
      {
	bem += 8;
	if (int(beam_set/1000) >= 2) 
	  ErrorMessage("non valid beam option").throwMe();
	beam_set -= 1000;
      }
    if (beam_set >=100)
      {
	bem += 4;
	if (int(beam_set/100) >= 2) 
	  ErrorMessage("non valid beam option").throwMe();
	beam_set -= 100;
      }
    if(beam_set >= 10)
      {
	bem += 2;
	if (int(beam_set/10) >= 2) 
	  ErrorMessage("non valid beam option").throwMe();
	beam_set -=10;
      }
    if (beam_set == 1) bem +=1;
    if (beam_set >= 2)  ErrorMessage("non valid beam option").throwMe();
  }




