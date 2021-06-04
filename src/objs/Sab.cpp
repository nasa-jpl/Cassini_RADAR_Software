#include "Sab.h"
#include "Array.h"
#include "Error.h"
#include "Io.h"
#include "Units.h"
#include "Plot.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "Constants.h"
#include "Ieb.h"
//-------------------
// Class Sab Methods
//-------------------

using std::cout;
using std::endl;
using std::min;
using std::cerr;

//constructor : initiallize baq_threshold(24), eng_tmp(12), 
//              and sub_comm_data(6,15)

Sab::Sab(FileMgr::orderE mt)
  : decoded_data("decoded_data",1),
    baq_threshold("baq_threshold"),
    eng_tmp("eng_tmp"),
    sub_comm_data("sub_comm_data"),
    tail_length_(0),Nrdat_(0),machine_type_(mt),decoded_(false),
    baq_decoded_(false),  is_sab_valid_(true),hip_hold_(0),cip_hold_(0)
{
  
  baq_threshold.resize(Nb_);
  baq_threshold = 0;
  eng_tmp.resize(Neng_tmp_); 
  eng_tmp = 0;

 
  // 6 x 16 (colum x row)
  sub_comm_data.resize(Nsub_col_,Nsub_row_);
  sub_comm_data = 0;
  decoded_data.resize(Ndat_max_);
  Ndecoded=0;
  for (unsigned int i = 0 ; i < Ndat_max_; i++) decoded_data(i) =0.0;

   //initiallize sub_comm_data variables to 0 except lpltlm_db (=1.0)
   fwdtmp= be1tmp= be2tmp= be3tmp=be4tmp =be5tmp= 0.0;
   diptmp=rlotmp=tadcal1=nsdtmp=lnatmp=evdtmp=0.0;
   mratmp=mruttm=dcgttm=cucttm=twttmp=epctmp=0.0;
   tw1ttm=ep1ttm=p_stmp=p_sttm=fguttm=tadcal4=0.0;
   esstmp=wgb1t1=wgb3t1=wgb3t2=wgb3t3=wgb5t1=0.0;
   pcutmp=adctmp=tadcal2=ecltmp=cputmp=memtmp=0.0;
   sadctmp=0.0;

   tadcal3=0.0;
   frwdpw=0.0;
   dcgmon=0;
   lpltlm_db=1.0;

   nsdcur=0.0;
   hpapsm=0.0;
   catcur=0.0;
   p_smon=0.0;
   svlsta=0.0;
   usotmp=0.0;

   cpbnkv = essvlt=tadcal5=pcu5v_72=0.0;
   pcu5i_73= 0.0;
   pcu5v_74=0.0;
   pcuti_75= 0.0;
   pcu15v_76=0.0;
   pcu15i_77=0.0;
   pcu15v_78=0.0;
   pcu15i_79=0.0;
   pcu12v_80=0.0;

   pcu12i_81=0.0;
   pcucur=0.0;
   pllmon=0.0;
   ctu5i=0.0;
   tadcal6=0;
   pcu9v_86=0.0;

   pcu9i_87=0.0;
   pcu9v_88=0.0;
   pcu9i_89=0.0;
 

   tadcl7=0.0;

   shpttm=0.0;
   
  }



// complete
bool Sab::complete() const
  {
  return(header_.size() == Sab::Nheader_bytes_ &&
         footer_.size() == Sab::Nfooter_bytes_ &&
         rdat_.size() == Nrdat_ &&
         tail_.size() == tail_length_);
  }

// partial
bool Sab::partial() const
  {
  return(!complete() && !empty());
  }

// empty
bool Sab::empty() const
  {
  return(header_.size() == 0);
  }

// decoded
bool Sab::decoded() const
  {
  return(decoded_);
  }
// baq_decoded
bool Sab::baq_decoded() const
  {
    return (baq_decoded_);
  }

//subtract 1 sec from sclk
void Sab::subtract1secfromSCLK()
  {
    if(!decoded_) ErrorMessage("Try to change sclk even before decoding, should not happen").throwMe();
    sclk -= 1;
  }

//decode sclk
void Sab::decodeSclk() throw(ErrorMessage)
  {
  if (partial()) throw IncompleteSab();
  if (empty()) return;

  //---------------------------------------
  // Setup arrays with uniform data fields.
  //---------------------------------------


  // SAB: header + science data + footer + tail
  // header_ = sync + sclk + scpr + brst+ slow + fast + header1
  // rdat_ : 0 - 16 k words
  // footer_= footer1 + footer2
  // tail_= tail

  //------------------------------
  // Bit information  from 0 to 80 
  //------------------------------

  sclk = getSabEntry(2,machine_type_,header_);
  }

//decode SAB counter
void Sab::decodeSABcounter() throw(ErrorMessage)
 {
   if (partial()) throw IncompleteSab();
   if (empty()) return;
   byte_int(sab_counter,machine_type_,FileMgr::unsigned_int,header_,66,2);
 }


//-------------------------------------------------------------------------
// decode
//
// Extract the useful telemetry channels from the header, footer, tail,
// and data.  These data channels are stored in other internal variables
// which are read and written by the I/O methods.  The raw header, footer,
// etc. data are only accepted through the addData method.
// If the current SAB is partially complete, an exception is thrown.
// If the current SAB is empty, nothing happens.
//-------------------------------------------------------------------------

//-----------------------------
//decode sab header,footer:
// Sab is expected to set the is_sab_valid_ false whenever
// there is an error/bit corruption instead of throwing an errormessage
// this is done to make sure the preprocessor runs withough crashing.
// whenever sab becomes invalid, its record will be dropped withoug
// writing into  LBDR, SBDR, IEBLIST, RASLISt
//----------------------------
void Sab::decode() 
{
  if (partial()) {
    is_sab_valid_ = false;
    error_code_ = -800;
    return;
  }
  if (empty()){
    is_sab_valid_ = false;
    error_code_= -700;
    return;
  }
  //---------------------------------------
  // Setup arrays with uniform data fields.
  //---------------------------------------


  // SAB: header + science data + footer + tail
  // header_ = sync + sclk + scpr + brst+ slow + fast + header1
  // rdat_ : 0 - 16 k words
  // footer_= footer1 + footer2
  // tail_= tail

  //------------------------------
  // Bit information  from 0 to 80 
  //------------------------------

  sync = getSabEntry(1,machine_type_,header_);
  sclk = getSabEntry(2,machine_type_,header_);
  unsigned int scpr_brst = getSabEntry(3,machine_type_,header_);
  scpr = bitget(scpr_brst,12,15);
  int_brst=bitget(scpr_brst,0,11);
  brst = Uvar(bitget(scpr_brst,0,11)*0.25*mstos, "s");
 
  //-----------------------------------------------------
  // Bit information from 80 to 352 (slow and fast field)
  //-----------------------------------------------------

  std::vector<unsigned short> slow;
  byte_int(slow,machine_type_,FileMgr::unsigned_int,header_,10,20);
  std::vector<unsigned short> fast;
  byte_int(fast,machine_type_,FileMgr::unsigned_int,header_,30,14);

  //--------------------------------
  // bit information from 352 to 640
  //--------------------------------
  std::vector<unsigned short> tnci;
  byte_int(tnci,machine_type_,FileMgr::unsigned_int,header_,44,10);
  std::vector<unsigned short> power;
  byte_int(power,machine_type_,FileMgr::unsigned_int,header_,54,4);

 
 
  //----------------------------------
  //Add tcn into ieb, check tnc instruction type
  //----------------------------------
  header_tfi=Uvar(bitget(tnci[0],0,15),"s");
  header_tnc = bitget(tnci[1],0,7);
  header_typ = bitget(tnci[1],8,9);
  
  header_tca=  tnci[2];
  header_tcb = tnci[3];
  header_tcc = tnci[4];
  byte_int(pwri,machine_type_,FileMgr::unsigned_int,header_,54,4);
  byte_int(vicc,machine_type_,FileMgr::unsigned_int,header_,58,2);
  byte_int(vimc,machine_type_,FileMgr::unsigned_int,header_,60,2);
  byte_int(tail_len,machine_type_,FileMgr::unsigned_int,header_,62,2);
  byte_int(tail_id,machine_type_,FileMgr::unsigned_int,header_,64,2);
  byte_int(sab_counter,machine_type_,FileMgr::unsigned_int,header_,66,2);
  byte_int(sab_len,machine_type_,FileMgr::unsigned_int,header_,68,2);
  byte_int(fswm,machine_type_,FileMgr::unsigned_int,header_,70,1);
  byte_int(fswc,machine_type_,FileMgr::unsigned_int,header_,71,1);
  byte_int(ctbc,machine_type_,FileMgr::unsigned_int,header_,72,1);
  byte_int(ctrx,machine_type_,FileMgr::unsigned_int,header_,73,1);
 
  byte_int(ctps_ctbe,machine_type_,FileMgr::unsigned_int,header_,74,2);
  ctps = bitget(ctps_ctbe,5,15);
  ctbe = bitget(ctps_ctbe,0,4);
  byte_int(header_end,machine_type_,FileMgr::unsigned_int,header_,76,4);

  //
  //bit information from footer1 and 2    
  //
  std::vector<unsigned short> footer1;
  byte_int(footer1,machine_type_,FileMgr::unsigned_int,footer_,0,30);
  std::vector<unsigned short> footer2;
  byte_int(footer2,machine_type_,FileMgr::unsigned_int,footer_,30,14);
 

  //-----------------------------------------
  // Extract data fields from uniform arrays using IEB
  //If an error occurs, it will make sab invalid
  //----------------------------------------- 
  ieb_.loadTncfield(tnci);
  ieb_.loadPowerfield(power);
  ieb_.loadSlowFastfield(slow,fast);
  try{
    ieb_.decodeSlowFastfield();
    slow_tfi =ieb_.getSlowTfi();
    dtn = ieb_.getDtn();
    slow_typ = ieb_.getSlowTyp();
 
    csr = ieb_.getCsr();
    r_mode = ieb_.getMod();
    slow_instruction_number= ieb_.getSin();
    bem = ieb_.getBem();
    baq_mode = ieb_.getBaq();
    tro = ieb_.getTro();
    tro_in_units_of_pri =  ieb_.getTroInPriUnits();
    
    rc_bw = ieb_.getRcv();
    adc = ieb_.getAdc();

    at1_db = ieb_.getAt1N1dB() + ieb_.getAt1N2dB() + ieb_.getAt1N3dB();
    at1 = cag_table_.gainLinearScale((int)at1_db);
  
    at3_db = ieb_.getAt3N1dB() + ieb_.getAt3N2dB() + ieb_.getAt3N3dB();
    at3 = cag_table_.gainLinearScale((int)at3_db);

    at4_db = ieb_.getAt4N1dB() + ieb_.getAt4N2dB() + ieb_.getAt4N3dB();
    at4 =cag_table_.gainLinearScale((int)at4_db);
    

    //transfer bit mask
    ieb_.getAttenuationBitMask(at1_db_mask, at3_db_mask,at4_db_mask);
    //cout<<"at1_db_mask "<<at1_db_mask<<endl;

    rip = ieb_.getRip();
    csd=ieb_.getCsd();
    rad=ieb_.getRad();
    csq = ieb_.getCsq();
    chirp_length = double(csq+1) * ieb_.getCsd();
    slow_cfs = ieb_.getCfs();
    
    fast_tfi = ieb_.getFastTfi();
    fin = ieb_.getFin();
    fast_typ = ieb_.getFastTyp();
    
    pul = ieb_.getPul();
    bii = ieb_.getBii();  
    bpd = ieb_.getBpd();
    pri = ieb_.getPri();
    rwd = ieb_.getRwd();
    fast_csf= ieb_.getChirpStartFrequency();
  }
  catch(ErrorMessage& e){
    cout<<"Ieb decoding failed, therefore decode individual field separately "<<endl;
    cout<<"this portion has not been tested because there has not been any invalid ieb from cruise data.  may fail and need update "<<endl;
    //----------------------------------------
    //Last resort: If ieb slow fast field fail
    //-----------------------------------------
    //--------------------------------------------------------------
    // This routine was made before IEB class is available so that slow and fast
    // field is decoded inside Sab.  As of now (July 22, 02), IEB class will
    // handle all the encoding and decoding processes.
    
    //-----------
    //slow field
    //------------
    slow_tfi = Uvar(bitget(slow[0],0,15),"s");
    dtn = bitget(slow[1],0,7);
    slow_typ = bitget(slow[1],8,9); 
    csr = bitget(slow[2],0,3);
    r_mode = bitget(slow[2],4,7);
    slow_instruction_number = bitget(slow[2],8,15);// slow field instruction number
    bem = bitget(slow[3],0,4);
    baq_mode = bitget(slow[3],5,7);
    tro_in_units_of_pri = bitget(slow[3],8,11) - bitget(slow[3],11,11)*16;
   
    double rc_bw_map[4]={0.117, 0.468,0.935,4.675};
    rc_bw= Uvar(rc_bw_map[bitget(slow[3],12,13)]*MHztoHz,"Hz");
    double adc_map[4]={0.25, 1.0, 2.0, 10.0};
    adc = Uvar(adc_map[bitget(slow[3],14,15)]*MHztoHz,"Hz");
    at1_db = bitget(slow[4],0,1)*4+bitget(slow[4],2,6)+bitget(slow[4],7,11);
    at1_db_mask=bitget(slow[4],0,11);
    at1 = cag_table_.gainLinearScale((int)at1_db);
   

    at3_db = bitget(slow[5],0,1)*4+bitget(slow[5],2,6)+bitget(slow[5],7,11);
    at3_db_mask=bitget(slow[5],0,11);
    at3 = cag_table_.gainLinearScale((int)at3_db);
    

    rip = Uvar(bitget(slow[6],0,3)*5*mstos,"s");
    at4_db = bitget(slow[6],4,5)*4+bitget(slow[6],6,10)+bitget(slow[6],11,15);
    at4_db_mask=bitget(slow[6],4,15);
    at4 =cag_table_.gainLinearScale((int)at4_db);
    
    csd = Uvar(bitget(slow[7],0,7)*133.3333333*nstos,"s");
    rad = bitget(slow[7],8,15);
    csq = bitget(slow[8],0,11);
    chirp_length
      = Uvar( (csq + 1) * bitget(slow[7],0,7)*133.333333*nstos ,"s");
    slow_cfs = Uvar(bitget(slow[9],0,15)*1.7881393,"Hz");
  
    //--------------------
    //fast field
    //------------------
    fast_tfi= Uvar(bitget(fast[0],0,15),"s");
    fin = bitget(fast[1],0,7);
    fast_typ=bitget(fast[1],8,9); 
    pul = bitget(fast[2],0,7); 
    bii = bitget(fast[2],8,15);  
    bpd = Uvar(bitget(fast[3],0,11)*mstos,"s");
    //see baq_decode() to understand the relation ship between pri 
    //information and actual pri in time
    //see blue page 4-68 for details
    // (pri + 2 ) * (4, 1, 0.5) micro-second 
    // or (pri(*10) + 6) * 0.1 micro-second
    //However, the FSW "will substract" 2 from pri setting when fadc != 10 MHz
    // and multiply 10 and substract 6 when fadc = 10 MHz
    // All this means that pri setting in fast field is actual pri
    
    double pri_map[4] = {4.0, 1.0, 0.5, 0.1};
    unsigned int i_fadc = bitget(slow[3],14,15);
    unsigned int i_pri = bitget(fast[5],0,9);
    if (i_fadc >= 0 && i_fadc <3){
      if (i_pri %2 != 0) i_pri++;//add one if it is not even
      pri = Uvar(pri_map[i_fadc] * i_pri*ustos, "s");
    }
    else if (i_fadc == 3){
      i_pri = i_pri * 10 ;
      pri = Uvar(pri_map[i_fadc] * i_pri*ustos, "s");
    }
    else{}
    tro = pri * tro_in_units_of_pri;
    rwd = bitget(fast[4],0,9) * pri;
    fast_csf= Uvar(bitget(fast[6],0,15)*457.764,"Hz");

    //-------------------------------------------------------------
    //when decoding slow and fast field fails, set valid flag false 
    // and return
    //-------------------------------------------------------------
    cerr<<"Error in decoding slow/fast field: "<<e.msg<<endl;
    cout<<"will make this sab invalid "<<endl;
    is_sab_valid_ = false;
    error_code_ = 5;
  }

  try{
    ieb_.decodePowerTncfield();
  }
  catch(ErrorMessage& e){
    //set valid flag false
    // and return
    cerr<<"Error in decoding power and tcn field"<<endl;
    cerr<<" ErrorMessage "<<e.msg<<endl;
    is_sab_valid_ = false;
    error_code_ = 7;
  }


 
 

  //-------------------------------------------------
  //footer1
  //
  //decoding footer1 depending on radar mode
  //r_mode <4 or r_mode>7 : ALT or Image mode --> BAQ
  //r_mode= 4,5,6,or 7: Inter-Galatic Object Calibration --> ENG Temp
  // first 192 bits (12 (=Neng_tmp_) x 16 bits) contain either
  // 12 16-bit eng tem values
  // or 24 8-bit baq_threshold values  
  // temperature values are stored as 12 bit information
  //------------------------------------------------------

  if (r_mode < 4 ||( r_mode > 7 && r_mode < 12) ) {//imaging mode 
    for (unsigned int  loop = 0; loop < Neng_tmp_; loop++){
      baq_threshold(2 * loop+1) = bitget(footer1[loop],0,7);  
      baq_threshold(2 * loop) = bitget( footer1[loop],8,15);
    }
  }
  else if (r_mode == 4) {//Radiometer only
    if (Nrdat_ != 0) {
      //if any error, set valid flag false and return
      cout<<"Error: non-zero radar data in passive mode"<<endl;
      is_sab_valid_ = false;
      error_code_ = -901;
    }
    for (unsigned int loop = 0; loop < Neng_tmp_; loop++)
      eng_tmp(loop) = bitget(footer1[loop],0,15);
    xfer1_eng_temp_data(); //transfer eng_tmp_data to tmp variables
    
    //screen output
    //cout<<"fwdtmp "<<fwdtmp<<endl;
    //cout<<"be3tmp "<<be3tmp<<endl;
    //cout<<"diptmp "<< diptmp <<endl;
    //cout<<"rlotmp "<< rlotmp <<endl;
    //cout<<"nsdtmp "<<   nsdtmp<<endl;
    //cout<<"lnatmp"<<   lnatmp <<endl;
    //cout<<"mratmp "<< mratmp<<endl;
    //cout<<"mruttm"<<   mruttm<<endl;
    //cout<<"wgb3t1 t2 t3  "<<wgb3t1<< wgb3t2<<   wgb3t3<<endl;
    //cout<<"nsdcur "<<nsdcur<<endl;
  }
  else{//this space should be set to 0
    for (unsigned int loop = 0; loop < Neng_tmp_; loop++){
      eng_tmp(loop) = bitget(footer1[loop],0,15);
      if (eng_tmp(loop) != 0) {
	//whenever this occurs, set valid flag false and 
	// return
       	cout<<"Sab::decode: BAQ block is not to set all zeros:see BB 7-5"<<endl;
	is_sab_valid_ = false;
	error_code_ = 241;
      }
    }
  }



  //-------------------------------------------------------
  // check zero filled baq mode
  //baq exponent should not be zero 
  //exception: when pri number is less than 192 in scatt mode
  //--------------------------------------------------------
  if ( (baq_mode==0 || baq_mode==1
	|| baq_mode==6 || baq_mode==7)&&
       (r_mode < 4 ||( r_mode > 7 && r_mode < 12))){//imaging mode 
    //imaging mode which requires baq threshold 
    if(baq_threshold.sum()==0){
      cout<<"baq_mode "<< baq_mode<<endl;
      cout<<"csr "<< csr<<endl;
      cout<<"Radar mode "<<r_mode<<endl;
      cout<<"all BAQ exponents are zero"<<endl;
      cout<<"pri * adc "<< pri*adc<<endl;
      is_sab_valid_ = false;
      error_code_= 242;
    }
  }//image mode

  unsigned int rads_r = getSabEntry(2,machine_type_,footer_);  
  cnt_rl = bitget(rads_r,0,11);
  cnt_radio = bitget(rads_r,12,31);

  unsigned short radd = getSabEntry(6,machine_type_,footer_);   
  cnt_nd = bitget(radd,0,11);
  eout = bitget(radd,12,14);

  subr = bitget(footer2[0],0,3);
  space_craft_time = bitget(footer2[0],4,15);
  
  if (subr == 15)
    {  // hip and cip invalid for subr == 15, use previous values instead.
    hip = hip_hold_;
    cip = cip_hold_;  
    iebtth =  footer2[1];      
    iebttl =  footer2[2]; 
    bgcalls = footer2[3];      
    delvmn =  footer2[4];     
    delvda =  footer2[5];      
    delvyr =  footer2[6];
    }
  else
    {
    hip = Uvar((bitget(footer2[2],12,15) + 1)*mstos,"s");
    cip = Uvar(bitget(footer2[1],12,15)*5*mstos,"s");
    hip_hold_ = hip;
    cip_hold_ = cip;     
    for ( int column = 0; column < Nsub_col_; column++)
      {
      sub_comm_data(column,subr) = bitget(footer2[column+1],0,11);
      }   
    xfer2_sub_comm_data();
    }      

  //-----------------
  //check validity of sab
  //-------------------
  if(subr!=15 && hip==Uvar(0,"s")){
    is_sab_valid_ = false;
    cout<<"noise diode integration time is zero"<<endl;
    cout<<"corrupted SAB data "<<endl;
    cout<<endl;
    error_code_ = 300;
  }
  else if (subr!=15 && cip==Uvar(0,"s")){
    is_sab_valid_=false;
    cout<<"resistive load integration time is zero"<<endl;
    cout<<"corrupted SAB data "<<endl;
    cout<<endl;
    error_code_ = 290;
    }
  else if(subr==15 && hip_hold_==Uvar(0,"s"))
    {
      is_sab_valid_= false;
      error_code_ = 301;
      cout<<"noise diode hold is "<< hip_hold_<<endl;
      cout<<"corrupted Sab data"<<endl;
    }
  else if(subr==15 && cip_hold_==Uvar(0,"s"))
    {
      is_sab_valid_=false;
      error_code_ = 291;
      cout<<"resitive load is "<<cip_hold_<<endl;
      cout<<"corrupted Sab data"<<endl;
    }
  decoded_ = true;
}



bool Sab::Valid() const   {
  return(is_sab_valid_);  
  //----------------------------------------------------------------------
  // is Sab valid after being decoded?
  // It is possible that ground station fails to receive 
  // whole sab packet and fills sab field with zeros.
  // In this case, although sab is complete (or all necessary data field),
  // some sab data fields (such as RL integration time) will have zero, which
  // is not valid.
  // 
  //This routine will check whether two parameters (RL and ND integration time)
  // are zero or not.  If that happens, it will return false
  // In addition, it will send out warning messages if ND, RL, and Radiometer
  // counts are zero.
  // list of things that are checked
  //  1: tnc_type
  //  2: pwr_type
  //  3: slow_typ
  //  4: fast_typ
  //  5: hip, hip_hold
  //  6: cip, cip_hold
  // -----------------------------------------------------------------------
}




// conversion of eng. tel. data when baq data are not reported
void Sab::xfer1_eng_temp_data() 
{
  fwdtmp = Uvar(tadc_conversion(0,eng_tmp(0)),"K");
  be3tmp = Uvar(tadc_conversion(0,eng_tmp(1)),"K");
  diptmp =  Uvar(tadc_conversion(0,eng_tmp(2)),"K" );
  rlotmp =  Uvar(tadc_conversion(0,eng_tmp(3)),"K" );
  nsdtmp =  Uvar(tadc_conversion(0,eng_tmp(4)),"K" );
  lnatmp =  Uvar(tadc_conversion(0,eng_tmp(5)),"K" );
  mratmp= Uvar( tadc_conversion(0,eng_tmp(6)),"K");
  mruttm= Uvar( tadc_conversion(0,eng_tmp(7)),"K");
  wgb3t1= Uvar(tadc_conversion(1,eng_tmp(8)),"K");
  wgb3t2= Uvar(tadc_conversion(1,eng_tmp(9)),"K");
  wgb3t3= Uvar(tadc_conversion(1,eng_tmp(10)),"K");
  nsdcur= Uvar(tadc_conversion(2, eng_tmp(11))*mAtoA,"A");
}


//
// conversion of sub_comm_data set
//
void Sab::xfer2_sub_comm_data() 
{
  if (subr == 15) {
    cout<<"no transfer because subr==15"<<endl;
    return;
  }
  if (subr == 0){
    fwdtmp = Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K");
    be1tmp = Uvar(tadc_conversion(0,sub_comm_data(1,subr)),"K");
    be2tmp = Uvar(tadc_conversion(0,sub_comm_data(2,subr)),"K");
    be3tmp = Uvar(tadc_conversion(0,sub_comm_data(3,subr)),"K");
    be4tmp = Uvar(tadc_conversion(0,sub_comm_data(4,subr)),"K");
    be5tmp = Uvar(tadc_conversion(0,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 1){  
    diptmp =  Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K" );
    rlotmp =  Uvar(tadc_conversion(0,sub_comm_data(1,subr)),"K" );
    tadcal1=  Uvar(tadc_conversion(0,sub_comm_data(2,subr)),"K" );
    nsdtmp =  Uvar(tadc_conversion(0,sub_comm_data(3,subr)),"K" );
    lnatmp =  Uvar(tadc_conversion(0,sub_comm_data(4,subr)),"K" );
    evdtmp =  Uvar(tadc_conversion(0,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 2){
    mratmp= Uvar( tadc_conversion(0,sub_comm_data(0,subr)),"K");
    mruttm= Uvar( tadc_conversion(0,sub_comm_data(1,subr)),"K");
    dcgttm= Uvar( tadc_conversion(0,sub_comm_data(2,subr)),"K");
    cucttm= Uvar( tadc_conversion(0,sub_comm_data(3,subr)),"K");
    twttmp= Uvar( tadc_conversion(0,sub_comm_data(4,subr)),"K");
    epctmp= Uvar( tadc_conversion(0,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 3){
    tw1ttm =Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K");
    ep1ttm= Uvar(tadc_conversion(0,sub_comm_data(1,subr)),"K");
    p_stmp= Uvar(tadc_conversion(0,sub_comm_data(2,subr)),"K");
    p_sttm= Uvar(tadc_conversion(0,sub_comm_data(3,subr)),"K");
    fguttm= Uvar(tadc_conversion(0,sub_comm_data(4,subr)),"K");
    tadcal4= Uvar(tadc_conversion(0,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 4){
    esstmp= Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K");
    wgb1t1= Uvar(tadc_conversion(1,sub_comm_data(1,subr)),"K");
    wgb3t1= Uvar(tadc_conversion(1,sub_comm_data(2,subr)),"K");
    wgb3t2= Uvar(tadc_conversion(1,sub_comm_data(3,subr)),"K");
    wgb3t3= Uvar(tadc_conversion(1,sub_comm_data(4,subr)),"K");
    wgb5t1= Uvar(tadc_conversion(1,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 5){
    pcutmp= Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K");
    adctmp= Uvar(tadc_conversion(0,sub_comm_data(1,subr)),"K");
    tadcal2= Uvar(tadc_conversion(1,sub_comm_data(2,subr)),"K");
    ecltmp= Uvar(tadc_conversion(0,sub_comm_data(3,subr)),"K");
    cputmp= Uvar(tadc_conversion(0,sub_comm_data(4,subr)),"K");
    memtmp= Uvar(tadc_conversion(0,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 6){
    sadctmp = Uvar(tadc_conversion(0,sub_comm_data(1,subr)),"K");
  }
  else if (subr == 7){
    tadcal3= Uvar(tadc_conversion(1,sub_comm_data(2,subr)),"K");
    double frwdpw_in_Watts=tadc_conversion(2,sub_comm_data(3,subr));
    frwdpw = Uvar(WtoMW*frwdpw_in_Watts,"MW");
    dcgmon= sub_comm_data(4,subr);
    //digital chirp generator monitor (see BB E-8 and 9)
    //range: 100 - 4000
    lpltlm_db= tadc_conversion(3,sub_comm_data(5,subr));
  }
  else if (subr == 8){   
    nsdcur= Uvar(tadc_conversion(2, sub_comm_data(0,subr))*mAtoA,"A");
    hpapsm= Uvar(tadc_conversion(4,sub_comm_data(1,subr)),"A");
    catcur= Uvar(tadc_conversion(5,sub_comm_data(2,subr))*mAtoA,"A");
    p_smon =Uvar(tadc_conversion(6, sub_comm_data(3,subr)),"A");
    svlsta =Uvar(sub_comm_data(4,subr)*VtoMV,"MV");
    usotmp= Uvar(tadc_conversion(7,sub_comm_data(5,subr)),"K");
  }
  else if (subr == 9){
    cpbnkv= Uvar(tadc_conversion(8,sub_comm_data(0,subr))*VtoMV,"MV");
    essvlt= Uvar(tadc_conversion(9,sub_comm_data(1,subr))*VtoMV,"MV");
    tadcal5= Uvar(sub_comm_data(2,subr)*VtoMV,"MV");
    pcu5v_72 =Uvar(tadc_conversion(10,sub_comm_data(3,subr))*VtoMV,"MV");
    pcu5i_73= Uvar(tadc_conversion(11,sub_comm_data(4,subr))*mAtoA,"A");
    pcu5v_74= Uvar(tadc_conversion(12,sub_comm_data(5,subr))*VtoMV,"MV");
  }
  else if (subr == 10){
    pcuti_75= Uvar(tadc_conversion(13, sub_comm_data(0,subr))*mAtoA,"A");
    pcu15v_76= Uvar(tadc_conversion(14, sub_comm_data(1,subr))*VtoMV,"MV");
    pcu15i_77= Uvar(tadc_conversion(15,sub_comm_data(2,subr))*mAtoA,"A");
    pcu15v_78 = Uvar(tadc_conversion(16,sub_comm_data(3,subr))*VtoMV,"MV");
    pcu15i_79= Uvar(tadc_conversion(17,sub_comm_data(4,subr))*mAtoA,"A");
    pcu12v_80= Uvar(tadc_conversion(18,sub_comm_data(5,subr))*VtoMV,"MV");
  }
  else if (subr == 11){
    pcu12i_81= Uvar(tadc_conversion(19,sub_comm_data(0,subr))*mAtoA,"A");
    pcucur= Uvar(tadc_conversion(20,sub_comm_data(1,subr))*mAtoA,"A");
    pllmon= Uvar(tadc_conversion(21,sub_comm_data(2,subr))*MHztoHz,"Hz");
    ctu5i = Uvar(tadc_conversion(22,sub_comm_data(3,subr))*mAtoA,"A");
    tadcal6= sub_comm_data(4,subr);//just number
    pcu9v_86= Uvar(tadc_conversion(23,sub_comm_data(5,subr))*VtoMV,"MV");
  }
  else if (subr == 12){
    pcu9i_87 = Uvar(tadc_conversion(24,sub_comm_data(0,subr))*mAtoA,"A");
    pcu9v_88 = Uvar(tadc_conversion(25,sub_comm_data(1,subr))*VtoMV,"MV");
    pcu9i_89 = Uvar(tadc_conversion(26,sub_comm_data(2,subr))*mAtoA,"A");
  }
  else if (subr == 13){
    tadcl7 =Uvar(sub_comm_data(5,subr)*VtoMV,"MV");
  }
  else{ 
    shpttm= Uvar(tadc_conversion(0,sub_comm_data(0,subr)),"K");
  }
}

//
// based on TADC converion table
//
double Sab::tadc_conversion(const unsigned short i,const unsigned int adc_sub_comm)
{ 
  double  x = adc_sub_comm;
  if (i == 0)
  {
    double tadc_map[]={0.029363, 23.4};
    return (  x *tadc_map[0] - tadc_map[1] + 273.15);
  }
  else if (i == 1)
  {
    double tadc_map[]={0.083382, 214.4};
    return ( x *tadc_map[0] - tadc_map[1]+273.15);
  }
  
  else if (i == 2) {return(x * x *0.000005);}
  else if( i == 3){return (pow(10,(x *0.001176+2.06)/10.0));}
  else if (i == 4){return (x * 0.00401 - 8.42);}
  else if (i == 5) {return (x *0.020173 + 8.1);}
  else if (i == 6){ return (x *0.000614+0.61);}
  else if (i == 7) {return (x *0.010452+71.1);}
  else if (i == 8){return (x *0.02867-10.65);}
  else if (i == 9) {return (x *0.0122);}
  else if (i ==10){ return (x * 0.001606-0.001);}
  else if (i == 11) { return (x *1.14623 -27.0);}
  else if (i ==12) {return (x *0.0015823-0.003);}
  else if (i ==13) {return (x *0.231158-8.5);}
  else if (i ==14) {return (x*0.00418+0.014);}
  else if (i ==15){return (x*0.060405-17.8);}
  else if (i ==16){return (x*0.004201+0.0036);}
  else if (i ==17){return (x*0.059303-11.8);}
  else if (i ==18){ return (x*0.003286+0.353);}
  else if (i ==19){return (x*0.014369-1.88);}
  else if (i ==20) {return (x*0.00101289-0.098);}
  else if (i ==21){return (x*0.000586+8.8243);}
  else if (i ==22){return (x*0.1861+10.7);}
  else if (i ==23){return (x*0.0025566+0.011);}
  else if (i ==24){return (x*0.102923-5.6);}
  else if (i==25){return (x*0.002469+0.209);}
  else if (i ==26){return (x*0.06414-3.8);}
  else{
    cout<<"out of range input value in tadc_conversion: will return 0"<<endl;
    return(0);
  }
}

//-----------------------------------------
//transfer baq threshold and raw data
// This routine is requested by l0_ras_extract.
// SH wants to use undecoded data set for his RAS tools
//-----------------------------------------
void Sab::extract_baq_raw_data(vector<unsigned short>&  baq, 
			       cdata& rdata) const
  {//will transport pure raw data without any conversion
   //raw data first
    rdata.clear();
    if(rdat_.size()==0) return;
    else
      {
	for(vector<unsigned char>::const_iterator p=rdat_.begin();
	    p != rdat_.end();p++) rdata.push_back(*p);
      }
    //baq
    baq.resize(Neng_tmp_);//twelve 16bit words
    std::vector<unsigned short> footer1;
    footer1.clear();
   
    byte_int(footer1,machine_type_,FileMgr::unsigned_int,footer_,0,30);
    for(unsigned int i=0;i<Neng_tmp_;i++) baq[i] = footer1[i];//first 12 16bit words
  }

//--------------
//transfer raw data
//-----------------
void Sab::get_raw_data(vector<unsigned char>& data){
  if(!complete() || partial())
    ErrorMessage("Sab is not ready to transfer data").throwMe();
  
  //clear the record
  data.clear();
  for(unsigned int i=0;i<header_.size();++i)
    data.push_back(header_[i]);

  
  if(rdat_.size()!=0){
    for(unsigned int i=0;i<rdat_.size();++i)
      data.push_back(rdat_[i]);
  }
  
  for(unsigned int i=0;i<footer_.size();++i){
    data.push_back(footer_[i]);
  }
  if(tail_.size()!=0)
    ErrorMessage("Sab: non-zero tail length").throwMe();
  return;
  
}

//________________________________________________________________
//
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
//--------------------------------------------------------------------------------

void Sab::baq_decode() 
{
  if(!decoded_) ErrorMessage("SAB header is not decoded");
  Ndecoded = 0;//reset Ndecoded
  rms_radar_data=0.0;//default value

  //declare variable
  int Ns_Nb, Ns, j;
  unsigned short iblock;
  double x;
  string input;
  unsigned int number_data_pri;
  int rw_tmp;
  unsigned int rw;

  //when radar is not in imaging mode
  if( (r_mode>3  && r_mode<8) || r_mode>11) {
    if (Nrdat_ != 0) {//passive data
      cout<<"Error: non-zero radar data in passive mode"<<endl;
      is_sab_valid_ = false;
      error_code_ = -901;
    }
    return;//passive mode, don't have to go down further
  }
  else {//active mode
    if(Nrdat_ == 0){
      cout<<"Warning::Nrdat_ is zero while radar is in active mode"<<endl;
      return;//sab is still valid, but no data 
    }
  }
  
  if(ieb_.decodeSlowFastfieldComplete()){
    rw_tmp = (int) (ieb_.getPul());
    rw_tmp += ieb_.getTroInPriUnits();
    rw = (unsigned int) rw_tmp;
  }
  else{
     rw_tmp = (int) pul;
     rw_tmp += round_double( (tro/pri).getInUnits(""));
     rw = (unsigned int) rw_tmp;
  }
  
  if (rw != ctrx){
    cout<<"sab number "<< sab_counter<<endl;
    cout<<"baq mode "<< baq_mode<<endl;
    cout<<"buffer size "<<rdat_.size()<<endl;
    cout<<"sab length "<<sab_len<<endl;
    cout<<"pri and adc "<< pri<<" "<<adc<<endl;
    //-----------------------------
    //In case pul+tro != ctrx: it happens in C27
    //-----------------------------
    if(ieb_.decodeSlowFastfieldComplete()){
      unsigned int expected_ctrx_from_ieb = (unsigned int) (int(ieb_.getPul())+ieb_.getTroInPriUnits());
      if(ctrx != expected_ctrx_from_ieb) ieb_.resetTroBasedOnCTRX(ctrx);
      Uvar old_tro= tro;      
      tro = ieb_.getTro();
      tro_in_units_of_pri=ieb_.getTroInPriUnits();
      cout<<"tro is adjusted from "<< old_tro<<endl;
      cout<<"new tro "<<tro<<endl;
    }
    else{
      Uvar old_tro= tro;
      int diff = int(ctrx) - int(rw);
      tro_in_units_of_pri += diff;
      tro = pri * tro_in_units_of_pri;
      cout<<"tro is adjusted from "<< old_tro<<endl;
      cout<<"new tro "<<tro<<endl;
    }
    //------------------------------------------------------------------
    //According to the blue book(12-26)
    //pul + tro = rw
    //however, I found a few cases in c27 where pul +tro is not
    // equal to ctrx
    // and ctrx produce a correct data size.  decide to use ctrx
    //pul : # of pulses in a transmitt burst in units of PRI
    //tro : offset between the transmitt burst and the receive window
    // in units
    //      of PRI
    //rw  : length of the receive window in units of PRI
    // when the baq mode is running, the system needs to know how many
    // pulses it has inside the receive window in order to process
    // the data (see C2)  "In any case, the BAQ will divide each PRI
    // into Nb blocks, each of which having Ns samples
    // From this, the total number of pulse rw = pul + tro
    // Again, see Page 12-27, during SAR, tro will tend to be negative
    //  (rw < pul)
    //----------------------------------------------------------
    cout <<"pul+tro_in_units_of_pri is not equal to ctrx ";
    cout<<"pul tro_in_units_of_pri ctrx values "<<pul<<" "<<tro_in_units_of_pri<<" "<<ctrx<<endl;
    cout<<"ctrx will be used to estimate the size of SAB body "<<endl;
    rw = ctrx;
  }
  if (rw == 0 )
    {
      cout<<"rw is zero while radar buffer is not"<<endl;
      is_sab_valid_ = false;
      error_code_ = 51;
      return;
    }    
  
  
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
  if(ieb_.decodeSlowFastfieldComplete())
    number_data_pri  = ieb_.getPointsInsidePRI();
  else 
    number_data_pri  = (unsigned int) round_double((pri*adc).getInUnits(""));

  Ndecoded = rw * number_data_pri;
  
  if ( (Ndecoded > Ndat_max_) &&
       (rdat_.size()!=Ndat_max_)){
    cout<<"--------- Warning --------------"<<endl;
    cout<<"adc "<<adc<<endl;
    cout<<"pri "<< pri<<endl;
    cout<<"prf "<<1/pri<<endl;
    cout<<"pul "<< pul<<endl;
    cout<<"baq "<<baq_mode<<endl;
    cout<<"mode "<<r_mode<<endl;
    cout<<"rw "<< rw<<endl;
    cout<<"number data pri "<< number_data_pri<<endl;
    cout<<" "<<Ndecoded<<endl;
    cout<<"tro "<< tro<<endl;
    cout<<"Ndat max "<<Ndat_max_<<endl;
    //	throw ErrorMessage("Ndecoded > Ndat_max_");
    cout<<"invalid sab "<<endl;
    cout<<"we can not handle this one "<<endl;
    cout<<"original data size: char array size "<< rdat_.size()<<endl;
    cout<<"sab length "<< sab_len<<endl;
    cout<<"will be returned to main program "<<endl;
    cout<<endl;
    Ndecoded=0;
    is_sab_valid_ = false;
    error_code_= -902;
    return;
  }
  else if( (Ndecoded>Ndat_max_) &&
	   (rdat_.size()==Ndat_max_)){
    cout<<"--------- Warning --------------"<<endl;
    cout<<"possible miscommanding "<<endl;
    cout<<"expected data size is larger than 16 kword"<<endl;
    cout<<"sab length size is full 16K words "<<endl;
    cout<<"adc "<<adc<<endl;
    cout<<"pri "<< pri<<endl;
    cout<<"prf "<<1/pri<<endl;
    cout<<"pul "<< pul<<endl;
    cout<<"baq "<<baq_mode<<endl;
    cout<<"mode "<<r_mode<<endl;
    cout<<"rw "<< rw<<endl;
    cout<<"number data pri "<< number_data_pri<<endl;
    cout<<" "<<Ndecoded<<endl;
    cout<<"tro "<< tro<<endl;
    cout<<"Ndat max "<<Ndat_max_<<endl;
    cout<<"we can not handle this one "<<endl;
    cout<<"original data size: char array size "<< rdat_.size()<<endl;
    cout<<"sab length "<< sab_len<<endl;
    cout<<"reset Ndecoded=16K words size "<<endl;
    is_sab_valid_ = false;
    error_code_= -902;
    Ndecoded=Ndat_max_;
  }
  else{}

  //----------------------------
  //Before decoding, take care of byte swapping
  //----------------------------   
  Array1D<unsigned char> radar_buffer("radar buffer",Nrdat_);
  const_CDI p= rdat_.begin();
  unsigned int i_buffer_count =0;
  for (p = rdat_.begin(); p != rdat_.end();p++)
    {
      radar_buffer(i_buffer_count) = *p;
      i_buffer_count++;
    }
  //----------------------
  //validity check
  //----------------------
  if (i_buffer_count != Nrdat_) {
    cout<<"Sab.cpp::Data number mismatch"<<endl;
    is_sab_valid_ = false;
    error_code_= -903;
    return;
  }
  if (Nrdat_ % 2 != 0) {
    cout<<"Sab buffere does not have even number of bytes"<<endl;
    is_sab_valid_ = false;
    error_code_ = -904;
    return;
  }
  
  
  //-------------------
  //Byte swapping
  //------------------
  for (unsigned int i = 0; i < Nrdat_-1;i = i+2){
    unsigned char a = radar_buffer(i);
    unsigned char b = radar_buffer(i+1);
    radar_buffer(i) = b;
    radar_buffer(i+1) = a;//byte swapping
  }
  
  
  //-------------------------------------------------------
  //cout <<"BAQ mode started"<<"baq_mode "<<baq_mode<<endl;
  //-------------------------------------------------------
  try{
    if(baq_mode == 0) {//8-to- 2 bit BAQ (SAR mode)
      unsigned int i_data_count = 0;    
      Ns_Nb = number_data_pri;
      Ns = int(Ns_Nb / Nb_);   
      double  d[4]={0.49,1.64,-0.49,-1.64};
      //00, 01, 10, 11 = 0.49 1.64 -0.49  -1.64    
      for (unsigned int i_buffer = 0; i_buffer < Nrdat_;++i_buffer)
      for (int i = 0; i < 8; i = i + 2){
	j = chrbitget(radar_buffer(i_buffer),i,i+1);
	iblock = i_data_count % Ns_Nb;
	iblock = int(iblock/Ns);
	if (iblock >= Nb_) {iblock = Nb_ -1 ;}
	x =d[j]*(float)(baq_threshold(iblock))/2.0;	
	decoded_data(i_data_count)= x;		  
	i_data_count++;
      }
      
      if (i_data_count != 4 * Nrdat_){
	throw ErrorMessage("baq=0:decoded data != 4 *Nrdat_");
      }
      if (i_data_count < Ndecoded ){
	throw ErrorMessage(" number of decodede< collected data");
      }  
      if (i_data_count >  Ndecoded ){
	//cout<<baq_mode<<endl;
	//cout<<i_data_count<<endl;
	//cout<<Ndecoded<<endl;
	//cout<<Nrdat_<<endl;
	unsigned int i_word = int(Ndecoded/8);
	unsigned int remain = Ndecoded%8;
	//number of words capable of storing Ndecoded
	if(remain!=0) 	i_word++;
	if(i_word != Nrdat_/2) throw ErrorMessage("SAB data counting error");
      }          
    }     	
    else if (baq_mode == 1) {//8-to-1 bit BAQ: sign bit and threshold only 
      unsigned int i_data_count=0;
      Ns_Nb = number_data_pri;
      Ns = int(Ns_Nb/Nb_);	  
      double  d[2]= {1.0, -1.0};//0 1 = 1.0 -1.0
      for (unsigned int i_buffer = 0; i_buffer < Nrdat_;++i_buffer)
      for (int i = 0; i < 8; i++){
	j = chrbitget(radar_buffer(i_buffer),i,i);
	iblock = i_data_count % Ns_Nb;
	iblock = int(iblock/Ns);
	if (iblock >= Nb_) {iblock = Nb_ -1 ;}
	x = d[j]*(float)(baq_threshold(iblock))/2.0;  	
	decoded_data(i_data_count) = x;
	i_data_count ++;
      }		
      
      if (i_data_count != 8 * Nrdat_)
	throw ErrorMessage("baq=1: decoded data != 8 *Nrdat_");
	
      if (i_data_count < Ndecoded )
	throw ErrorMessage(" number of decodede< collected data");
      
      if (i_data_count != Ndecoded){
	//cout<<baq_mode<<endl;
	//cout<<i_data_count<<endl;
	//cout<<Ndecoded<<endl;
	//cout<<Nrdat_<<endl;	    
	unsigned int i_word = int(Ndecoded/16);
	unsigned int remain= Ndecoded%16;
	if(remain!=0) i_word++;//number of word capable of storing Ndecoded
	if (i_word != Nrdat_/2) throw ErrorMessage("SAB data counting error");
      }  
    }  
    else if (baq_mode == 2){
      cout<<"8-to 0 (standard sci data but no data)"<< baq_mode<<endl;
      Ndecoded = 0;
      return;
    }
    else if (baq_mode == 3) {
      //ALTL DATA COMPRESSION
      //ALTL compression method produces 16 bit data
      //number of data:number_data_pri
      //data type: sum of absolute value
      //However, rdat_ is 8 bit.  So, need to add even and odd set to
      //construct new data

      //--------------------------
      //temporary holder
      //size = Nrdat_ :same as 8 bit size
      //real data size : pri_number (Nrdat_ /2)
      // Caution: first 32 bits corresponds to DC offset
      //First check data size
      //-------------------------------------------------

      if ( (Nrdat_- 4) !=  2* number_data_pri)
	ErrorMessage("Sab.cpp::ALTL compressed: number_data_pri does not match to 8-bit buffer size").throwMe();
      
      //---------------------------
      //Second, check if twos compliment was used
      //---------------------------
      scatt_dc_offset = 0.0;
      unsigned int sign = chrbitget(radar_buffer(0),7,7);
      if (sign == 0)
	{//normal representation
	scatt_dc_offset +=  chrbitget(radar_buffer(0),0,7)*pow(2,24);
	scatt_dc_offset +=  chrbitget(radar_buffer(1),0,7)*pow(2,16);
	scatt_dc_offset +=  chrbitget(radar_buffer(2),0,7)*pow(2,8);
	scatt_dc_offset +=  chrbitget(radar_buffer(3),0,7);    
	}
      else 
	{//twos compliment was used
	// (2^B - N) : B = 32, N= normal way of reading
	// change sign when the first bit is 1
	scatt_dc_offset = pow(2,32);
	scatt_dc_offset -=  chrbitget(radar_buffer(0),0,7)*pow(2,24);
	scatt_dc_offset -=  chrbitget(radar_buffer(1),0,7)*pow(2,16);
	scatt_dc_offset -=  chrbitget(radar_buffer(2),0,7)*pow(2,8);
	scatt_dc_offset -=  chrbitget(radar_buffer(3),0,7);  
	scatt_dc_offset *= -1.0;
	}
 
      //cout<<"total number of data"<<rw*number_data_pri<<endl;
      //cout<<"altl dc offset"<<scatt_dc_offset/double(rw*number_data_pri)<<endl;
	
      //for (unsigned int i_byte = 0; i_byte < 4;++i_byte)
      //{
      //  cout<<i_byte<<"th byte "<<endl;
      //  for (unsigned int i_bit =0; i_bit < 8;++i_bit)
      //    {
      //    cout<<chrbitget(radar_buffer(i_byte),i_bit,i_bit);
      //    }
      //  cout<<endl;
      //}

      Ndecoded = number_data_pri;
      for (unsigned int i = 0; i < number_data_pri;++i)
	{
	unsigned int value1 = radar_buffer(2*i+4);
	unsigned int value2 = radar_buffer(2*i+1+4);
	decoded_data(i) = float(value1) + float(value2)*256.0;
	}
      decoded_data(number_data_pri)= scatt_dc_offset;
      //cout<<"offset "<< decoded_data(number_data_pri)<<endl;
      }	
    else if (baq_mode == 4) {
      // 8-4 MSB
      unsigned int i_data_count = 0;
      double d[8]={8.0,24.0,40.0,56.0,72.0,88.0,104.0,120.0};
      // 0000 0001 0001 ....1000 1001..
      for (unsigned int i_buffer = 0; i_buffer < Nrdat_;++i_buffer)
	{
	j = chrbitget(radar_buffer(i_buffer),0,2);  
	decoded_data(i_data_count) = d[j]; 
	if (chrbitget(radar_buffer(i_buffer),3,3)==1) {decoded_data(i_data_count)= -decoded_data(i_data_count);}
	i_data_count ++;  	    
	j = chrbitget(radar_buffer(i_buffer),4,6);  
	decoded_data(i_data_count)=d[j] ; 
	if (chrbitget(radar_buffer(i_buffer),7,7) == 1)  {decoded_data(i_data_count)= -decoded_data(i_data_count);}
	i_data_count ++;
	}	        
      if (i_data_count != 2 * Nrdat_)
	throw ErrorMessage("baq=4:decoded data != 2 *Nrdat_");
      
      if (i_data_count < Ndecoded )
	throw ErrorMessage("baq_mode=4: decoded_data < actual data");
      
      if (i_data_count != Ndecoded ){
	//cout<<baq_mode<<endl;
	//cout<<i_data_count<<endl;
	//cout<<Ndecoded<<endl;
	//cout<<Nrdat_<<endl;
	unsigned int i_word = int(Ndecoded/4);
	unsigned int remain=Ndecoded%4;
	if(remain!=0) i_word++;//number of word capable of storing Ndecoded
	if (i_word != Nrdat_/2) throw ErrorMessage("SAB data counting error");
	}
    }
    else if(baq_mode ==5) {//8-to-8 straight
      unsigned int i_data_count = 0;
      double  d[2]={1.0, -1.0};	   
      for (unsigned int i_buffer = 0; i_buffer < Nrdat_;++i_buffer)
	{	     
	  j = chrbitget(radar_buffer(i_buffer),7,7);
	  x = ( (float)(chrbitget(radar_buffer(i_buffer),0,6)) + 0.5 ) *d[j];  
	  decoded_data(i_data_count) = x;
	  i_data_count ++;
	}
      if (i_data_count !=  Nrdat_)
	throw ErrorMessage("baq=5:decoded data != Nrdat_");
	
      if (i_data_count < Ndecoded )
	throw ErrorMessage("baq_mode=5: decoded_data < actual data");
      
      if (i_data_count >  Ndecoded)
	ErrorMessage("Sab data counting error:8 bit straingt: decoded data is larger than computed one based on commands").throwMe();
      
    }
    else if (baq_mode ==6 || baq_mode ==7 ){
      //low and high ALT mode 8-to-4 BAQ
      unsigned int i_data_count = 0;
      Ns_Nb= number_data_pri;
      Ns = int(Ns_Nb/Nb_);// number of data per each baq_threshold	     
      double d[8]={0.117,0.355,0.600,0.861,1.148,1.479,1.891,2.498};
      for (unsigned int i_buffer = 0; i_buffer < Nrdat_;++i_buffer)
	{ 		
	j = chrbitget(radar_buffer(i_buffer),0,2);
	iblock = i_data_count % Ns_Nb;
	iblock = int(iblock/Ns);
	if (iblock >= Nb_) {iblock = Nb_ -1 ;}
	x =  d[j]*(float)(baq_threshold(iblock))/2.0;  
	decoded_data(i_data_count) = x;
	if(chrbitget(radar_buffer(i_buffer),3,3)==1) {decoded_data(i_data_count) = -decoded_data(i_data_count);}
	i_data_count++;
	j = chrbitget(radar_buffer(i_buffer),4,6);
	iblock = i_data_count % Ns_Nb;
	iblock = int(iblock/Ns);
	if (iblock >= Nb_) {iblock = Nb_ -1 ;}
	x =  d[j]*(float)(baq_threshold(iblock))/2.0;  
	if (i_data_count > Ndat_max_) {throw ErrorMessage("i_data_count > Ndat_max_");}		
	decoded_data(i_data_count)=x;
	if(chrbitget(radar_buffer(i_buffer),7,7) == 1) {decoded_data(i_data_count) = -decoded_data(i_data_count);}
	i_data_count++;		
	}
      if (i_data_count != 2* Nrdat_)
	throw ErrorMessage("baq=6 or 7:decoded data != 2 *Nrdat_");
      
      if (i_data_count < Ndecoded )
	throw ErrorMessage("baq_mode=6 or 7: decoded_data < actual data");
      
      if (i_data_count != Ndecoded ){
	//cout<<baq_mode<<endl;
	//cout<<i_data_count<<endl;
	//cout<<Ndecoded<<endl;
	//cout<<Nrdat_<<endl;
	unsigned int i_word = int(Ndecoded/4);
	unsigned int remain=Ndecoded%4;
	if(remain!=0) i_word++;//number of word capable of storing Ndecoded
	if (i_word != Nrdat_/2) throw ErrorMessage("SAB data counting error");
      }	
    }
    else{ 
      Ndecoded = 0;
      throw ErrorMessage("no baq_mode found");
    }
  }
  catch(ErrorMessage& e){
    cerr<<"Error in decoding baq : invalid sab"<<e.msg<<endl;
    is_sab_valid_= false;
  }
  if(Ndecoded!=0) rms_radar_data = decoded_data.rms(Ndecoded);
  //cout<<"Active data set: Ndecoded "<<Ndecoded<<endl;
  baq_decoded_=true;
}

//----------------------------------------------------------------------------
// addData(record)
//
// This method accumulates data from an input record read from a file written
// with the machine type set by the constructor.
// Repeated calls add more data until the SAB is complete.
// Calling this method after the SAB is complete will cause an exception.
// The caller is responsible for passing the records in the correct order
// (header,data,footer,tail).
//----------------------------------------------------------------------------

void Sab::addData(const Sab::cdata& record) throw(ErrorMessage)
  {

  // irec indexes the next element to pull from the input record.
  unsigned irec = 0;

  if (complete() && irec < record.size()) throw OverfullSab();

  //---------------------------------
  // Accumulate data from the record.
  //---------------------------------

  // Header
  if (header_.size() < Sab::Nheader_bytes_ && irec < record.size())
    {
    // Nixxx is the minimum of what we want vs what is available.
    // We can't xfer more than what is available.
    int Niheader = min(Sab::Nheader_bytes_-header_.size(), record.size()-irec);
    copy(record.begin()+irec, record.begin()+irec + Niheader,
              back_inserter(header_));
    irec += Niheader;
    }

  if (irec > 69)
    {  // Now we have enough of the header to pull out size fields.
    //-----------------------------------------------------------------------
    // Pull out data length fields from the header.
    // Tail and data lengths are in 16 bit words, so they are multiplied by 2
    // to get the lengths in bytes.
    //-----------------------------------------------------------------------
    byte_int(tail_length_,machine_type_,FileMgr::unsigned_int,record,62,2);
    tail_length_ *= 2;
    byte_int(Nrdat_,machine_type_,FileMgr::unsigned_int,record,68,2);
    Nrdat_ *= 2;
    }

  // Echo data (not always present)
  if (header_.size() == Sab::Nheader_bytes_ &&
      rdat_.size() < Nrdat_ && irec < record.size())
    {
    // Finished the header data, not finished with the echo data,
    // and there is still more data in this record, so xfer more echo data.
    int Nidata = min(Nrdat_ - rdat_.size(), record.size() - irec);
    copy(record.begin()+irec, record.begin()+irec + Nidata,
              back_inserter(rdat_));
    irec += Nidata;
    }

  // Footer
  if (Nrdat_ == rdat_.size() && footer_.size() < Sab::Nfooter_bytes_ &&
      irec < record.size())
    {
    // Finished the echo data, not finished with the footer,
    // and there is still more data in this record, so xfer more footer data.
    int Nifooter = min(Sab::Nfooter_bytes_-footer_.size(), record.size()-irec);
    copy(record.begin()+irec, record.begin()+irec + Nifooter,
              back_inserter(footer_));
    irec += Nifooter;
    }

  // Tail (not always present)
  if (footer_.size() == Sab::Nfooter_bytes_ &&
      tail_.size() < tail_length_ && irec < record.size())
    {
    // Finished the footer data, not finished with the tail,
    // and there is still more data in this record, so xfer more tail data.
    int Nitail = min(tail_length_ - tail_.size(), record.size() - irec);
    copy(record.begin()+irec, record.begin()+irec + Nitail,
              back_inserter(tail_));
    irec += Nitail;
    }

  if (irec < record.size())
    {  // finished the SAB but there is still more data in this record!
    throw ErrorMessage("Extra data in record\n");
    }

  }


void Sab::addPartialRecord(const Sab::cdata& record) throw(ErrorMessage)
  {

  // irec indexes the next element to pull from the input record.
  unsigned irec = 0;

  if (complete() && irec < record.size()) throw OverfullSab();

  //---------------------------------
  // Accumulate data from the record.
  //---------------------------------

  // Header
  if (header_.size() < Sab::Nheader_bytes_ && irec < record.size())
    {
    // Nixxx is the minimum of what we want vs what is available.
    // We can't xfer more than what is available.
    int Niheader = min(Sab::Nheader_bytes_-header_.size(), record.size()-irec);
    copy(record.begin()+irec, record.begin()+irec + Niheader,
              back_inserter(header_));
    irec += Niheader;
    }



  if (header_.size() > 69)
    {  // Now we have enough of the header to pull out size fields.
    //-----------------------------------------------------------------------
    // Pull out data length fields from the header.
    // Tail and data lengths are in 16 bit words, so they are multiplied by 2
    // to get the lengths in bytes.
    //-----------------------------------------------------------------------
    byte_int(tail_length_,machine_type_,FileMgr::unsigned_int,header_,62,2);
    tail_length_ *= 2;
    byte_int(Nrdat_,machine_type_,FileMgr::unsigned_int,header_,68,2);
    Nrdat_ *= 2;
    }
  
  

  // Echo data (not always present)
  if (header_.size() == Sab::Nheader_bytes_ &&
      rdat_.size() < Nrdat_ && irec < record.size())
    {
    // Finished the header data, not finished with the echo data,
    // and there is still more data in this record, so xfer more echo data.
    int Nidata = min(Nrdat_ - rdat_.size(), record.size() - irec);
    copy(record.begin()+irec, record.begin()+irec + Nidata,
              back_inserter(rdat_));
    irec += Nidata;
    }

  // Footer
  if (Nrdat_ == rdat_.size() && footer_.size() < Sab::Nfooter_bytes_ &&
      irec < record.size())
    {
    // Finished the echo data, not finished with the footer,
    // and there is still more data in this record, so xfer more footer data.
    int Nifooter = min(Sab::Nfooter_bytes_-footer_.size(), record.size()-irec);
    copy(record.begin()+irec, record.begin()+irec + Nifooter,
              back_inserter(footer_));
    irec += Nifooter;
    }

  // Tail (not always present)
  if (footer_.size() == Sab::Nfooter_bytes_ &&
      tail_.size() < tail_length_ && irec < record.size())
    {
    // Finished the footer data, not finished with the tail,
    // and there is still more data in this record, so xfer more tail data.
    int Nitail = min(tail_length_ - tail_.size(), record.size() - irec);
    copy(record.begin()+irec, record.begin()+irec + Nitail,
              back_inserter(tail_));
    irec += Nitail;
    }

  if (irec < record.size())
    {  // finished the SAB but there is still more data in this record!
    throw ErrorMessage("Extra data in record\n");
    }

  //debug purpose
  //if(rdat_.size()==Nrdat_ && Nrdat_!=0){
  //cout<<"header dat footer tail "<< header_.size()<<" "<<rdat_.size()<<" "<<footer_.size()<<" "<<tail_.size()<<endl;
  //cout<<"tail length "<< tail_length_<<endl;
  //}

  }

//------------
// getSabEntry
//------------

int Sab::getSabEntry(unsigned int tablenum, FileMgr::orderE machine_type,
  const cdata& data) const throw(ErrorMessage)
  {
  int result;
  if (data.size() == Nheader_bytes_)
    {
    if (tablenum < 0 || tablenum > header_item_bytes_.size())
      {
      throw ErrorMessage("Error: tablenum out of range in getSabEntry");
      }
    unsigned int field_pos = Sab::header_item_byte_pos_[tablenum-1];
    unsigned int field_length = Sab::header_item_bytes_[tablenum-1];
    byte_int(result,machine_type,FileMgr::unsigned_int,data,
             field_pos,field_length);
    }
  else if (data.size() == Nfooter_bytes_)
    {
    if (tablenum > footer_item_bytes_.size())
      {
      throw ErrorMessage("Error: tablenum out of range in getSabEntry");
      }
    unsigned int field_pos = Sab::footer_item_byte_pos_[tablenum-1];
    unsigned int field_length = Sab::footer_item_bytes_[tablenum-1];
    byte_int(result,machine_type,FileMgr::unsigned_int,data,
             field_pos,field_length);
    }
  else
    {
    throw ErrorMessage("Error: Invalid data array passed to getSabEntry");
    }
  return(result);
  }

//------------------------------
//Clear method
//------------------------------

void Sab::clear()
  {
  header_.clear();
  footer_.clear();
  tail_.clear();
  rdat_.clear();
  Nrdat_ = 0;
  tail_length_ = 0;
  decoded_ = false;
  baq_decoded_=false; 
  is_sab_valid_ = true;//assume sab is valid when restarting
  eng_tmp = 0;
  baq_threshold = 0;
  decoded_data = 0.0;//set to 0.0
  rms_radar_data=0.0;//set to 0.0
  ieb_.clear();
  error_code_ = 0;
  }

//----------------------------
//getIeb
//----------------------------
Ieb Sab::getIeb() const
 {
 if(!decoded_) ErrorMessage("Sab.cpp::getIeb(Ieb& ieb):sab header is not decoded").throwMe();
   return(ieb_);
 }

//--------------
//get error location
//-------------
int Sab::getErrorInfo()
{
  return(error_code_);
}


/*** The show routines below need cdata stream output operators defined

//---------------
// showHeader
//---------------
void Sab::showHeader() const
  {
  std::cout << "Sab header:" << std::endl;
  std::cout << header_;
  std::cout << std::endl;
  }

//----------------
// showFooter
//----------------
void Sab::showFooter() const
  {
  std::cout << "Sab footer:" << std::endl;
  std::cout << footer_;
  std::cout << std::endl;
  }

// showTail
void Sab::showTail() const
  {
  std::cout << "Sab tail:" << std::endl;
  std::cout << tail_;
  std::cout << std::endl;
  }

// showData
void Sab::showData() const
  {
  std::cout << "Sab data:" << std::endl;
  std::cout << rdat_;
  std::cout << std::endl;
  }

**/

//--------------------------------------------------------------------
// This constructor builds a byte index table from a byte element
// size table.  The sizes are summed cumulatively (starting from 0) to
// give the corresponding zero-offset positions.
//--------------------------------------------------------------------

Sab::ByteIndexTable::ByteIndexTable(Sab::const_CDI p1, Sab::const_CDI p2)
  {
  this->push_back(0);
  copy(p1, p2-1, back_inserter(*this));
  std::partial_sum(this->begin(), this->end(), this->begin());
  }

//-------------------------------------------------------------------------
// hb_ and fb_ give the number of bytes in each element read or written as
// one unit (SAB header/footer structure).
// hb_ and fb_ are both expanded by the number of SAB table entries for each
// element with zeros filling the subsequent positions for each entry after
// the initial entry.  This ensures that cumulatively summing the table will
// give the correct byte offsets for each SAB table entry.  SAB table entries
// are counted by rows in the Cassini RADAR HLD.
// hnum_ and fnum_ are the same as hb_ and fb_ except that the zeros are
// replaced by duplicates of the byte count.
//-------------------------------------------------------------------------

const unsigned char Sab::hb_[] = {4,4,0,2,20,14,10,4,2,2,2,2,2,2,
  1,1,1,1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4};
const unsigned char Sab::fb_[] = {24,0,4,0,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2};

const unsigned char Sab::hnum_[] = {4,4,2,2,20,14,10,4,2,2,2,2,2,2,
  1,1,1,1,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
const unsigned char Sab::fnum_[] = {24,4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

//-----------------------------------------------------------------------
// The byte_pos tables index the byte position of an entry in the SAB
// header table (BB, 7-2).  The _bytes tables index the size of an entry
// in bytes.
//-----------------------------------------------------------------------

const Sab::ByteIndexTable Sab::header_item_byte_pos_ =
  Sab::ByteIndexTable((Sab::const_CDI)&Sab::hb_[0],(Sab::const_CDI)&Sab::hb_[sizeof(Sab::hb_)-1]);
const Sab::cdata Sab::header_item_bytes_ =
  Sab::cdata(&Sab::hnum_[0],&Sab::hnum_[sizeof(Sab::hnum_)-1]);

const Sab::ByteIndexTable Sab::footer_item_byte_pos_ =
  Sab::ByteIndexTable((Sab::const_CDI)Sab::fb_,(Sab::const_CDI)&Sab::fb_[sizeof(Sab::fb_)-1]);
const Sab::cdata Sab::footer_item_bytes_ =
  Sab::cdata(&Sab::fnum_[0],&Sab::fnum_[sizeof(Sab::fnum_)-1]);

//----------------------------------------//
// I/O Supporting functions and operators //
//----------------------------------------//

// stream output for cdata
std::ostream& operator<<(std::ostream& s, const Sab::cdata& c)
  {
  for (unsigned int i=0; i < c.size(); i++)
    {
    s << int(c[i]) << " ";
    }
  return s;
  }

// stream output for vector<unsigned short>
std::ostream& operator<<(std::ostream& s,
  const std::vector<unsigned short>& sint)
  {
  for (unsigned int i=0; i < sint.size(); i++)
    {
    std::cout << sint[i] << " ";
    }
  return s;
  }


//------------------------------------------------
//
// show all the data on the screen 
//---------------------------------------------------
void Sab::show_data_Onscreen()
  {
  string input; 
  cout<<"sync and sclk "<<sync << " "<<sclk<<endl;
  cout<<"scpr "<<scpr<<endl;
  cout<<"brst "<<brst <<endl;  
  cout <<"header tfi tnc typ: "<< header_tfi<< " "<<header_tnc << " "
       << header_typ<<endl;
  cout <<"header tca tcb tcc:"<<header_tca<<" "<<header_tcb <<" "
       <<header_tcc <<endl;
  cout<<"pwri vicc vimc tail_id sab_counter: "<<pwri<<" "<<vicc<<" "
      <<vimc<<" "<<tail_id<<" "<<sab_counter<<endl;
  cout <<"Nrdat "<< Nrdat_<<endl;

  cout <<"fswm ctbc ctrx ctps ctbe header_end:"<<fswm<<" "<<ctbc<<" "
       <<ctrx<<" "<<ctps<<" "<<ctbe<<" "<<header_end<<endl;

  //-----------------------------------------
  // Extract data fields from uniform arrays.
  //-----------------------------------------
  //slow field
  //for variable names, see blue book  12-2
  cout <<"slow tfi dtn slow type csr r_mode slow_instruction_number: "<<" "<<
    slow_tfi <<" "<<dtn <<" "<<slow_typ <<" "<<csr <<" "<<
    r_mode <<" "<<slow_instruction_number <<endl;
  cout << "bem baq_mode tro_in_units_of_pri rc_bw adc: "<<bem <<" "
       <<baq_mode <<" "<<tro_in_units_of_pri <<" "
       <<rc_bw<< " "<< adc <<endl;
  cout <<"rip: "<<rip<<endl;
  cout << "at1db  at3_db  at4_db csd rad: "<<" "
       <<at1_db<<" "<<at3_db<<" "<<at4_db<<" "<<csd
       <<" "<<rad<<endl;
  cout <<"csq chirp length csf:"<<csq<<" "<<chirp_length<<" "<<slow_cfs<<endl;
 
  //fast field
  cout << "fast_tfi fin fast_type pul bii rwd  fast_csf: "<<fast_tfi<<" "<<
    fin<<" "<<fast_typ<<" "<<pul<<" "<<bii<<" " <<rwd<<" "<<fast_csf<<endl;
  cout<<"bpd and pri "<< bpd<<" "<<pri<<endl;
  
  //footer1
  std::vector<unsigned short> footer1;
  byte_int(footer1,machine_type_,FileMgr::unsigned_int,footer_,0,30);
  std::vector<unsigned short> footer2;
  byte_int(footer2,machine_type_,FileMgr::unsigned_int,footer_,30,14);
  cout <<"footer1 "<< footer1<<endl;
  cout<<"footer2"<<footer2<<endl;

  //decoding footer1 depending on radar mode
  //r_mode <4 or r_mode>7 : ALT or Image mode --> BAQ
  //r_mode= 5,6,or 7: Inter-Galatic Object Calibration --> ENG Temp  
  if (r_mode <4 ||( r_mode > 7 && r_mode < 12) )
    {
      for (unsigned int loop = 0; loop < Nb_; loop++)
	{
	cout <<loop<<"th baq value: "<<baq_threshold( loop)<<endl;
	}
      unsigned int j  = 0;
      for (unsigned int i = 0; i < Ndecoded; i++)
	{ 
	if (j < 12 || j> Ndecoded - 12){
	  cout <<j<<"th value: "<<decoded_data(i)<<endl;}
	j++;
	}
    }
  else 
    {
    cout <<"eng temp"<<endl;
    cout<<"fwdtmp "<<  fwdtmp <<endl;
    cout<< "be3tmp "<<be3tmp <<endl;
    cout<< "diptmp "<<diptmp <<endl;
    cout<< "rlotmp "<<  rlotmp<<endl;
    cout<< "nsdtmp "<< nsdtmp <<endl;
    cout<< "lnatmp "<< lnatmp <<endl;
    cout<< "mratmp "<<mratmp<<endl;
    cout<< "mruttm "<<mruttm<<endl;
    cout<<  "wbg3t1 "<< wgb3t1<<endl;
    cout<<  "wgb3t2 "<<wgb3t2<<endl;
    cout<<"wgb3t3 "<< wgb3t3<<endl;
    cout<< "nsdcur "<< nsdcur<<endl;
  }
  cout <<"cntrl crnradio cntnd eout: "<<cnt_rl<<" "<<cnt_radio<<" "<<cnt_nd<<
    " "<<eout<<endl;
  if (subr == 15)
    {  // hip and cip invalid for subr == 15, use previous values instead.
      cout <<" hip cip iebtth iebttl bgcalls devlmn devlda devlyr"<< hip
	   <<endl
	   <<"cip "<<cip <<endl
	   <<"iebtth "<<iebtth<<endl 
	   <<"iebttl "<<iebttl<<endl
	   <<"bgcalls "<<bgcalls<<endl
	   <<"delvmn "<<delvmn<<endl
	   <<"delvda "<<delvda<<endl
	   <<"delvyr "<<delvyr<<endl;
    }
  else
    {
      cout<< "hip and cip hol"<<endl<<
	hip_hold_ <<endl<<
	cip_hold_ <<endl;     
    }  
  cout <<"space craft time"<<space_craft_time<<endl;
  cout <<"sub comm data out"<<subr<<endl;

 if (subr == 0)
   {
     cout <<
       "fwdtmp"<<   fwdtmp <<endl<<
       "be1tmp "<<be1tmp <<endl<<
       "be2tmp "<<be2tmp<<endl<<
       "be3tmp "<<be3tmp <<endl<<
       "be4tmp "<<be4tmp<<endl<<
       "be5tmp "<<be5tmp<<endl; 
   }
 else if (subr == 1)
 {   cout <<
      " diptmp "<<endl<<
       "rlotmp"<<rlotmp <<endl<<
       "tdacall"<<tadcal1<<endl<<
       "nsdtmp "<<nsdtmp <<endl<<
       "lnatmp "<<lnatmp<<endl<<
       "evdtmp "<<evdtmp<<endl;
 }
 else if (subr == 2)
   { cout <<
       mratmp<<endl<<
       mruttm<<endl<<
       dcgttm<<endl<<
       cucttm<<endl<<
       twttmp<<endl<<
       epctmp<<endl;
   }
 else if (subr == 3)
   { cout <<
       tw1ttm <<endl<<
       ep1ttm<<endl<<
       p_stmp<<endl<<
       p_sttm<<endl<<
       fguttm<<endl<<
       tadcal4<<endl;
   }
 else if (subr == 4)
   {cout <<
      esstmp<<endl<<
      wgb1t1<<endl<<
      wgb3t1<<endl<<
      wgb3t2<<endl<<
      wgb3t3<<endl<<
      wgb5t1<<endl;
   }
 else if (subr == 5)
   {cout<<
      pcutmp<<endl<<
      adctmp<<endl<<
      tadcal2<<endl<<
      ecltmp<<endl<<
      cputmp<<endl<<
      memtmp<<endl;
   }
 else if (subr == 6)
   {cout <<
      sadctmp <<endl;
   }
 else if (subr == 7)
   {cout<<
      tadcal3<<endl<<
      frwdpw <<endl<<
      dcgmon<<endl<<
      lpltlm_db<<endl;
   }
 else if (subr == 8)
   {cout <<
      nsdcur<< endl<<
      hpapsm<< endl<<
      catcur<< endl<<
      p_smon << endl<<
      svlsta<< endl<<
      usotmp<< endl;
   }
 else if (subr == 9)
   {cout<<
      cpbnkv<<endl<<
      essvlt<<endl<<
      tadcal5<<endl<<
      pcu5v_72 <<endl<<
      pcu5i_73<<endl<<
      pcu5v_74<<endl;
   }
 else if (subr == 10)
   {cout<<
      pcuti_75<<endl<<
      pcu15v_76<<endl<<
      pcu15i_77<<endl<<
      pcu15v_78<<endl<<
      pcu15i_79<<endl<<
      pcu12v_80<<endl;
   }
 else if (subr == 11)
   {cout<<
      pcu12i_81<<endl<<
      pcucur<<endl<<
      pllmon<<endl<<
      ctu5i <<endl<<
      tadcal6<<endl<<
      pcu9v_86<<endl;
   }
 else if (subr == 12)
   {cout<<
      pcu9i_87 <<endl<<
      pcu9v_88 <<endl<<
      pcu9i_89 <<endl;
   }
 else if (subr == 13)
   { cout<<
       tadcl7 <<endl;
   }
 else if (subr == 14)
   { cout <<
       shpttm<<endl;
   }
 else
   {cout <<"subr is 15"<<endl;
   } 
 
}


//-----------------------------------
//write Sab data into xmgr plot
//-----------------------------------
void Sab::generatePlot()
{
  if (!baq_decoded_) 
    ErrorMessage("No Sab data decoded").throwMe();

  unsigned int number_data_pri = ieb_.getPriEncodedValue();
  
  Dvec index("index",Ndecoded);
  Dvec radar_echo_decoded("echo",Ndecoded);
  for (unsigned int i = 0; i < Ndecoded;++i)
    {
      index(i) = i;
      radar_echo_decoded(i) = decoded_data(i);
    }
  
  Plot plot1;
  plot1.addXY(index,radar_echo_decoded,
	      line("solid","black"),sym("none"));
  plot1.setTitle("baq pri ctrx csr "+toStr(baq_mode)+" "+toStr(number_data_pri)+" "+toStr(ctrx)+" "+toStr(csr));
  plot1.show("x");
}

//------------------------------------------
//Calculate standard deviation and squared deviation
// of radar buffer (SAB data)
//------------------------------------------
void Sab::get_Sab_data_statistics(double& standard_deviation,double& average,Uvar& receiver_bandwidth, double& attenuation)
  {
    if (!baq_decoded_) 
      ErrorMessage("No Sab data decoded").throwMe();
    if (Nrdat_==0)
      ErrorMessage("Data buffer size is zero").throwMe();
    if (bem !=4)
      ErrorMessage("Calibration with beam 3 only!").throwMe();
    float sum = 0.0;
    for (unsigned int i = 0; i < Ndecoded;++i)
      {
      sum += decoded_data(i);
      }
    average = double(sum)/double(Ndecoded);
    float dev2= 0.0;
    for (unsigned int i = 0; i < Ndecoded;++i)
      {
      dev2 += (decoded_data(i)-average)*(decoded_data(i)-average)/float(Ndecoded-1);
      }
    standard_deviation = sqrt(double(dev2));
    receiver_bandwidth = rc_bw;
    attenuation = at3;//assume beam 3
    //debugging purposes
    //cout<<"beam3 gain setting "<<ieb_.getAt3N1dB();
    //cout<< ieb_.getAt3N2dB()<<" "<< ieb_.getAt3N3dB()<<endl;
    //cout<<"total gain "<<at3_db<<endl;
  }


void Sab::constructActiveRecord(const Time& t,
				const unsigned int& burst_number,
				const unsigned int& beam_number,
				Ieb& ieb, 
				const Array1D<float>& radar_data)
  {
    //set each field size
    header_.resize(Sab::Nheader_bytes_);
    footer_.resize( Sab::Nfooter_bytes_) ;
    rdat_.resize(Nrdat_) ;
    tail_.resize(0);

    sync =2004118378;
    Time t1=t;
    sclk = t1.sclk("Cassini");
    t1.setSclk("Cassini", sclk);
    brst=t - t1;
  
    slow_tfi =ieb.getSlowTfi();
    dtn = ieb.getDtn();
    slow_typ = ieb.getSlowTyp();
    
    csr = ieb.getCsr();
    r_mode = ieb.getMod();
    slow_instruction_number= ieb.getSin();
    bem = ieb.getBem();
    baq_mode = ieb.getBaq();
    tro = ieb.getTro();
    tro_in_units_of_pri =  ieb.getTroInPriUnits();
    
    rc_bw = ieb.getRcv();
    adc = ieb.getAdc();
    at1_db = ieb.getAt1N1dB() + ieb.getAt1N2dB() + ieb.getAt1N3dB();
    at1 = pow(10,-(float) (at1_db)/10.0);
    at3_db = ieb.getAt3N1dB() + ieb.getAt3N2dB() + ieb.getAt3N3dB();
    at3 = pow(10, -(float) (at3_db)/10.0);
    at4_db = ieb.getAt4N1dB() + ieb.getAt4N2dB() + ieb.getAt4N3dB();
    at4 = pow(10, -(float)(at4_db)/10.0);
    rip = ieb.getRip();
    csd=ieb.getCsd();
    rad=ieb.getRad();
    csq = ieb.getCsq();
    chirp_length = double(csq+1) * ieb.getCsd();
    slow_cfs = ieb.getCfs();
    
    fast_tfi = ieb.getFastTfi();
    fin = ieb.getFin();
    fast_typ = ieb.getFastTyp();
    
    pul = ieb.getPul();
    bii = ieb.getBii();  
    bpd = ieb.getBpd();
    pri = ieb.getPri();
    rwd = ieb.getRwd();
    fast_csf= ieb.getChirpStartFrequency();
   
   
     rad = ieb.getRad();
     //derived quantity
     int int_ctrx=int(ieb.getPul()) + ieb.getTroInPriUnits();
     if(int_ctrx<0) ErrorMessage("Sab.cpp:constructActiveRecord():pul+tro is less than 0").throwMe();
     ctrx =(unsigned int) int_ctrx;
     sab_len=ieb.getNumberOfWords();
    
     //burst number
     sab_counter = burst_number;
     //beam number setting
     ctbe=1;
     for(unsigned int i=1;i<beam_number;++i) ctbe *=2;

     tail_len=0;
     if(radar_data.size() != ieb.getRadarDataPoints()) ErrorMessage("radar data size does not match to ieb's computed one").throwMe();
     Ndecoded = radar_data.size();

     for(unsigned int i=0;i<radar_data.size();++i) 
       decoded_data(i)=radar_data(i);
     rms_radar_data=radar_data.rms(radar_data.size());
     decoded_ = true;
     baq_decoded_=true;
    
  }









