#include "L1A.h"
#include "Error.h"
#include <string>
#include <stdlib.h>
#include <string.h>

using std::cout;
using std::endl;
using std::cin;

//-----------------
// Methods for L1A
//-----------------

//--------------
// Constructors
//--------------

L1A::L1A(const std::string& filename, const std::string& mode,const std::string& filetype)
  : BurstData(filename,mode,filetype), decoded_data("decoded_data",1),
    radar_data_loaded_(false)
  {

  header_size_= L1A_HEADER_LENGTH;
  if (filetype =="active")
    {
      Ndat_max_=16*16*1024;//16k (16x1024x16)
      record_size_=L1A_PASSIVE_RECORD_LENGTH+sizeof(Ndecoded)+
	Ndat_max_*sizeof(float);
      parameters_.appendParameter("Ndecoded",Ndecoded);
      parameters_.appendParameter("decoded_data_full",Ndat_max_*sizeof(float),
			      Parameter::FDATA,(void*) &decoded_data);

    }
  else if (filetype =="passive"){
    Ndat_max_=1;
    record_size_=L1A_PASSIVE_RECORD_LENGTH;
  }
  else
    { 
     ErrorMessage e("no proper filetype: either active or passive"); 
     e.throwMe(); 
    }
    
  decoded_data.resize(Ndat_max_);
  for (unsigned int i = 0 ; i < Ndat_max_; i++) decoded_data(i)=0.0;


  }

L1A::~L1A()
{
  if(mode_[0]=='w') rewriteHeader();
}

//--------------
// I/O
//--------------

//------------------------------------------------------------------
// loadSab(sab)
//
// Copy relevant data from the passed sab to this L1A's data fields.
// If the passed SAB is empty or partially complete, an exception is thrown.
// If the passed SAB has not been decoded, an exception is thrown.
//------------------------------------------------------------------

void L1A::loadSab(const Sab& sab) 
  {
  if (sab.partial()){
    L1AError e(L1A::incomplete_sab);
    e.throwMe();
  }
  if (sab.empty()){
    L1AError e(L1A::empty_sab);
    e.throwMe();
  }
  if (!sab.decoded()){
    L1AError e(L1A::undecoded_sab);
    e.throwMe();
  }


 //header field
  sync = sab.sync;
  sclk= sab.sclk;

  // decode CDS pickup rate
  double value;
  if ( sab.scpr==0x03 ) value=30400.0;
  else if (sab.scpr==0x0C ) value=364800.0;
  scpr= Uvar(value,"Hz"); // bits/sec
  
  brst= sab.brst;

  header_tfi=sab.header_tfi;
  header_tnc=sab.header_tnc;
  header_typ=sab.header_typ;
  header_tca=sab.header_tca;
  header_tcb= sab.header_tcb;
  header_tcc= sab.header_tcc;

  pwri=sab.pwri;
  vicc=sab.vicc;
  vimc=sab.vimc;
  tail_len=sab.tail_len;
  tail_id=sab.tail_id;
  sab_counter=sab.sab_counter;
  sab_len=sab.sab_len;
  fswm=sab.fswm;
  fswc=sab.fswc;
  ctbc=sab.ctbc; 
  ctrx = sab.ctrx; 

  ctps= sab.ctps;
  ctbe=sab.ctbe;
  ctps_ctbe=sab.ctps_ctbe;
  header_end=sab.header_end;

  //slow field variables
  slow_tfi=sab.slow_tfi;
  dtn= sab.dtn;
  slow_typ= sab.slow_typ;
  csr=sab.csr;  
  r_mode= sab.r_mode; 
  slow_instruction_number= sab.slow_instruction_number; 
  bem=sab.bem; 
  baq_mode=sab.baq_mode ; 
  tro=sab.tro;
  
  rc_bw= sab.rc_bw;  
  adc=sab.adc;

  at1= sab.at1;  
  at3= sab.at3; 
  at4= sab.at4;
  at1_db= sab.at1_db; at3_db= sab.at3_db; at4_db= sab.at4_db;  
  rip = sab.rip; 
  csd= sab.csd;
  rad = sab.rad; 

  csq= sab.csq;
  chirp_length= sab.chirp_length;  

  slow_cfs= sab.slow_cfs;

  //fast field variables
  fast_tfi= sab.fast_tfi;
  fin= sab.fin; 
  fast_typ= sab.fast_typ; 
  pul = sab.pul; 
  bii = sab.bii;
  bpd = sab.bpd; 

  pri = sab.pri; 
  rwd = sab.rwd; 
  fast_csf = sab.fast_csf;

  //sci data field
  
  Ndecoded= sab.Ndecoded; 
  if (ft_==active && Ndecoded != 0)
  {
    decoded_data= sab.decoded_data;
    radar_data_loaded_ = true;
  }
else
  {
    radar_data_loaded_ = false;
  }

 

  cnt_rl= sab.cnt_rl; 
  cnt_radio= sab.cnt_radio; 
  cnt_nd= sab.cnt_nd;
  eout= sab.eout; 
  subr= (unsigned int) sab.subr; 
  space_craft_time= sab.space_craft_time;
  hip= sab.hip; 
  cip=sab.cip;
  
  iebtth= sab.iebtth;  
  iebttl= sab.iebttl;
  bgcalls= sab.bgcalls; 
  delvmn= sab.delvmn; 
  delvda= sab.delvda; 
  delvyr= sab.delvyr;
    
  //sub_comm_data variables 
  //subr = 0
  fwdtmp= sab.fwdtmp; 
  be1tmp=sab.be1tmp; 
  be2tmp=sab.be2tmp;
  be3tmp=sab.be3tmp;
  be4tmp=sab.be4tmp;
  be5tmp=sab.be5tmp;

 //subr = 1
  diptmp=sab.diptmp;
  rlotmp=sab.rlotmp;
  tadcal1=sab.tadcal1;
  nsdtmp=sab.nsdtmp;
  lnatmp=sab.lnatmp;
  evdtmp=sab.evdtmp;
  
  //subr = 2
  mratmp=sab.mratmp; 
  mruttm=sab.mruttm; 
  dcgttm=sab.dcgttm;
  cucttm=sab.cucttm;
  twttmp=sab.twttmp;
  epctmp=sab.epctmp;
  
  //subr = 3
  tw1ttm=sab.tw1ttm;
  ep1ttm=sab.ep1ttm;
  p_stmp=sab.p_stmp;
  p_sttm=sab.p_sttm;
  fguttm=sab.fguttm;
  tadcal4=sab.tadcal4;
  
 //subr = 4
  esstmp=sab.esstmp;
  wgb1t1=sab.wgb1t1;
  wgb3t1=sab.wgb3t1;
  wgb3t2=sab.wgb3t2;
  wgb3t3=sab.wgb3t3;
  wgb5t1=sab.wgb5t1;
  
 //subr = 5
  pcutmp=sab.pcutmp;
  adctmp=sab.adctmp;
  tadcal2=sab.tadcal2;
  ecltmp=sab.ecltmp;
  cputmp=sab.cputmp;
  memtmp=sab.memtmp;
  
 //subr = 6
  sadctmp=sab.sadctmp;

 //subr = 7
  tadcal3=sab.tadcal3; 
  frwdpw=sab.frwdpw;
  dcgmon=sab.dcgmon;
  lpltlm_db=sab.lpltlm_db;
  

 //subr = 8
  nsdcur=sab.nsdcur;
  hpapsm=sab.hpapsm; 
  catcur=sab.catcur;
  p_smon=sab.p_smon;
  svlsta=sab.svlsta;
  usotmp=sab.usotmp;
  
//subr = 9
  cpbnkv= sab.cpbnkv; 
  essvlt=sab.essvlt;
  tadcal5=sab.tadcal5;
  pcu5v_72=sab.pcu5v_72;
  pcu5i_73=sab.pcu5i_73;
  pcu5v_74=sab.pcu5v_74;

//subr = 10
  pcuti_75=sab.pcuti_75;
  pcu15v_76=sab.pcu15v_76;
  pcu15i_77=sab.pcu15i_77;
  pcu15v_78=sab.pcu15v_78;
  pcu15i_79=sab.pcu15i_79; 
  pcu12v_80=sab.pcu12v_80;

//subr = 11
  pcu12i_81=sab.pcu12i_81;
  pcucur=sab.pcucur; 
  pllmon=sab.pllmon; 
  ctu5i=sab.ctu5i;
  tadcal6=sab.tadcal6;
  pcu9v_86=sab.pcu9v_86;
 
//subr = 12
  pcu9i_87=sab.pcu9i_87;
  pcu9v_88=sab.pcu9v_88;
  pcu9i_89=sab.pcu9i_89; 

//subr = 13
  tadcl7= sab.tadcl7 ;

//subr = 14 
  shpttm=sab.shpttm;
 
  }

//------------------------------------------------------------------
// writeHeader()
//
// Write a L1A header to the current L1A file in binary.
// This method defines the L1A Header format.
//------------------------------------------------------------------

void L1A::writeHeader() 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1AError e("Can't write to input file " + filename_, L1A::write_error);
      e.throwMe();
    }

  // Should not call this routine twice
  if (header_handled_)
    {
      L1AError e("Attempt to write header twice", L1A::write_error);
      e.throwMe();
    }

  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",0,record_size_);
  header=headbuf;
  if ( header.length() != L1A_HEADER_LENGTH){
    L1AError e("Header length is incorrect",L1A::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  }

void L1A::rewriteHeader(){
  int position = file_.getPosition();
  file_.rewind();
  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",record_count_,record_size_);
  header=headbuf;
  if ( header.length() != L1A_HEADER_LENGTH){
    L1AError e("Header length is incorrect",L1A::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  file_.setPosition(position);
}

//------------------------------------------------------------------
// readHeader()
//
// Read a L1A header from the current L1A file in binary.
// This method defines the L1A Header format.
//------------------------------------------------------------------

void L1A::readHeader() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1AError e("Can't read from output file " + filename_,
	       L1A::read_error);
      e.throwMe();
    }
  if (header_handled_)
    {
      L1AError e("Attempt to read header twice",
	       L1A::read_error);
      e.throwMe();
    }
  checkEndedness();
  string header;
  header.resize(L1A_HEADER_LENGTH);
  file_.read(header);
  string token=get_next_token(header," /n/t"); // skip text
  token=get_next_token(header," /n/t"); // skip text
  
  token=get_next_token(header," /n/t"); // get num records
  record_count_=atoi(token.c_str());
  records_counted_=true;

  token=get_next_token(header," /n/t"); // skip text
  token=get_next_token(header," /n/t"); // skip text

  token=get_next_token(header," /n/t"); // check record size
  unsigned int recsize=(unsigned int) atoi(token.c_str());
  if(recsize!=record_size_){
    L1AError e("Record size incorrect in " + filename_, L1A::read_error);
    e.throwMe();
  }
  header_handled_ = true;
  }

//------------------------------------------------------------------
// writeRecord()
//
// Write the current L1A record to the end of the current L1A file in binary.
// This method defines the L1A record format.
// It is the responsibility of the user to ensure that valid data is present
// before calling this routine.  No tracking of data status is performed.
//------------------------------------------------------------------

void L1A::writeRecord(int abs_record_no) 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1AError e("Can't write to input file " + filename_, L1A::write_error);
      e.throwMe();
    }
  if (!header_handled_)
    {
      L1AError e("Can't write L1A record until header is written",
		 L1A::write_error);
      e.throwMe();
    }

  if(file_.getPosition() > MAX_FILE_LENGTH-(int)record_size_) 
    createNextFile();
  //-----------------------------------------------------------------
  // radar_data_loaded_ was set in loadSab
  // radar_data_loaded_ is true only when radar is in active mode and 
  //                           number of decoded data is not zero
  //----------------------------------------------------------------

  if (ft_==active && radar_data_loaded_==false) return;


  int start_position=file_.getPosition();
  writePassiveSABData(abs_record_no);
 
  if (ft_== active &&  radar_data_loaded_==true)
  {  
    file_.write(Ndecoded);   
    //      file_.write(decoded_data);      
    for (unsigned int i = 0; i < Ndecoded; i++){file_.write(decoded_data(i));}
    for (unsigned int i = Ndecoded; i< Ndat_max_; i++){file_.write(float(0.0));}
    
  }

 
 
  int end_position=file_.getPosition();
  int len_written=end_position-start_position;
  if((unsigned int)len_written!=record_size_){
    L1AError("Incorrect record size "+toStr(len_written)+" should be"+
	     toStr(record_size_),L1A::write_error).throwMe();
  }
  record_count_++;

  }


//------------------------------------------------------------------
// readRecord()
//
// Read a L1A record from the current L1A file.
//------------------------------------------------------------------

void L1A::readRecord() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1AError e("Can't read from output file " + filename_,
	       L1A::read_error);
//
// create next file
      e.throwMe();
    }  if (!header_handled_)
    {
      L1AError e("Can't read L1A record until header is read",
	       L1A::read_error);
      e.throwMe();
    } 

  // This check has the side effect of moving to the next  file if necessary
  if(eof())
    {
      L1AError e("Unexpected EOF in L1A file "+filename_,
	       L1A::read_error);
      e.throwMe();
    } 


  // Read in all public variables
  //convert this part to read
  //header slow fast data footer
 
  int start_position=file_.getPosition();
  
  readPassiveSABData();  
 
  if (ft_==active)
  {     
      file_.read(Ndecoded);     
      //      file_.read(decoded_data);      
      for (unsigned int i = 0; i <Ndecoded; i++){file_.read(decoded_data(i));}
      int position=file_.getPosition();
      position+=sizeof(float)*(Ndat_max_-Ndecoded);
      file_.setPosition(position);
  }
  
  
 
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    L1AError e("Incorrect record size",L1A::read_error);
    e.throwMe();
  }
  }






//-----------------------------------------------------
//
//
// this function is used to print read-in data on the screen
//
//
//-----------------------------------------------------


void L1A::show_passive_data_on_screen() 
{
 
 
 
  string input; 

  cout<<"sync"<< sync<<endl;
  cout<<"sclk"<<sclk<<endl; 
  cout<<"scpr"<<scpr<<endl;
  cout<<"brst:"<<brst <<endl;  
  cout <<"header tfi tnc typ: "<< header_tfi<< " "<<header_tnc << " "<< header_typ<<endl;
  cout <<"header tca tcb tcc:"<<header_tca<<" "<<header_tcb <<" "<<header_tcc <<endl;
  cout<<"pwri vicc vimc tail_id sab_counter: "<<pwri<<" "<<vicc<<" "<<vimc<<" "<<tail_id<<" "<<sab_counter<<endl;
  
  cout <<"fswm ctbc ctrx ctps ctbe header_end:"<<fswm<<" "<<ctbc<<" "<<ctrx<<" "<<ctps<<" "<<ctbe<<" "<<header_end<<endl;

getline(cin, input); if (input=="q"){return;}
 
  //-----------------------------------------
  // Extract data fields from uniform arrays.
  //-----------------------------------------

  //slow field
  //for variable names, see blue book  12-2
  cout <<"slow tfi dtn slow type csr r_mode sin: "<<" "<<
    slow_tfi <<" "<<dtn <<" "<<slow_typ <<" "<<csr <<" "<<
    r_mode <<" "<<slow_instruction_number<<endl;
  cout << "bem baq_mode tro rc_bw adc: "<<bem <<" "<<baq_mode <<" "<<tro <<" "<<
    rc_bw<< " "<< adc <<endl;
  cout <<"rip: "<<rip<<endl;
  cout << "at1db  at3_db  at4_db csd rad: "<<" "
    <<at1_db<<" "<<at3_db<<" "<<at4_db<<" "<<csd
       <<" "<<rad<<endl;
 
  cout <<"csq chirp length csf:"<<csq<<" "<<chirp_length<<" "<<slow_cfs<<endl;
  getline(cin, input); if (input == "q"){return;}

  
  //
  //fast field
  //
  cout << "fast tfi fin fast type pul bii bpd rwd pri:"<<fast_tfi<<" "<<
    fin<<" "<<fast_typ<<" "<<pul<<" "<<bii<<" "<<bpd<<" "<<rwd<<" "<<pri<<" "<<fast_csf
    <<endl;


    
  
  //decoding footer1 depending on radar mode
  //r_mode <4 or r_mode>7 : ALT or Image mode --> BAQ
  //r_mode= 5,6,or 7: Inter-Galatic Object Calibration --> ENG Temp  
 
    cout <<"eng temp"<<endl;
    cout<<  fwdtmp <<endl;
   cout<< be3tmp <<endl;
   cout<< diptmp <<endl;
  cout<<   rlotmp<<endl;
   cout<< nsdtmp <<endl;
  cout<<  lnatmp <<endl;
   cout<< mratmp<<endl;
   cout<< mruttm<<endl;
  cout<<  wgb3t1<<endl;
  cout<<  wgb3t2<<endl;
   cout<< wgb3t3<<endl;
   cout<< nsdcur<<endl;
 
 

  cout <<"cntrl crnradio cntnd eout: "<<cnt_rl<<" "<<cnt_radio<<" "<<cnt_nd<<
    " "<<eout<<endl;

 
  getline(cin, input); if (input=="q"){return;}
  
  if (subr == 15)
  {  // hip and cip invalid for subr == 15, use previous values instead.
    cout <<" hip cip iebtth iebttl bgcalls devlmn devlda devlyr"<< hip<<endl<<
    cip <<endl<<
    iebtth<<endl<< 
    iebttl<<endl<<
    bgcalls<<endl<<
    delvmn<<endl<<
    delvda<<endl<<
    delvyr<<endl;
  }
  else
  {
    cout<< "hip and cip hol"<<endl<<
    hip<<endl<<
    cip<<endl;     
    
  }  
  cout <<"space craft time"<<space_craft_time<<endl;
  cout <<"sub comm data out"<<subr<<endl;

 if (subr == 0)
 {cout <<
   fwdtmp <<endl<<
   be1tmp <<endl<<
   be2tmp<<endl<<
   be3tmp <<endl<<
   be4tmp<<endl<<
   be5tmp<<endl
; 
 }
 else if (subr == 1)
 {   cout <<
   diptmp<<endl<<
   rlotmp <<endl<<
   tadcal1<<endl<<
   nsdtmp <<endl<<
   lnatmp<<endl<<
   evdtmp<<endl;
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
   {     cout <<
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


//--------Parameter Access ----//
double L1A::getSpecialParam(const char* name)
{
  double retval;
  if(strcmp(name,"rms_decoded_data")==0){
    retval=decoded_data.rms(Ndecoded);
  }
  else if(strcmp(name,"norm_cnt_rl")==0){
    Uvar r=normalizeCntRL();
    retval=r.getInUnits("1/s");
  }
  else if(strcmp(name,"norm_cnt_nd")==0){
    Uvar r=normalizeCntND();
    retval=r.getInUnits("1/s");
  }
  else if(strcmp(name,"norm_cnt_radio")==0){
    Uvar r=normalizeCntRadio();
    retval=r.getInUnits("1/s");
  }
  else retval=BurstData::getSpecialParam(name);
  return(retval);
}
void L1A::enableSpecialParam(const char* name)
{
  if(strcmp(name,"rms_decoded_data")==0){
    parameters_.enable("full_record");
  }
  else if(strcmp(name,"norm_cnt_rl")==0){
    parameters_.enable("cnt_rl");
    parameters_.enable("cip");
  }
  else if(strcmp(name,"norm_cnt_nd")==0){
    parameters_.enable("cnt_nd");
    parameters_.enable("hip");
  }
  else if(strcmp(name,"norm_cnt_radio")==0){
    parameters_.enable("cnt_radio");
    parameters_.enable("rip");
    parameters_.enable("rad");
  }
  else{
    BurstData::enableSpecialParam(name);
  }
}
void L1A::disableSpecialParam(const char* name)
{
  if(strcmp(name,"rms_decoded_data")==0){
    parameters_.disable("full_record");
  }
  else if(strcmp(name,"norm_cnt_rl")==0){
    parameters_.disable("cnt_rl");
    parameters_.disable("cip");
  }
  else if(strcmp(name,"norm_cnt_nd")==0){
    parameters_.disable("cnt_nd");
    parameters_.disable("hip");
  }
  else if(strcmp(name,"norm_cnt_radio")==0){
    parameters_.disable("cnt_radio");
    parameters_.disable("rip");
    parameters_.disable("rad");
  }
  else{
    BurstData::disableSpecialParam(name);
  }
}

Uvar L1A::normalizeCntRL()
{
 Uvar ncnt_rl=(cnt_rl+radiometer_offset_)/(cip+radiometer_delta_tau_);
 return(ncnt_rl);
}

Uvar L1A::normalizeCntND()
{
  Uvar ncnt_nd=(cnt_nd+radiometer_offset_)/(hip+radiometer_delta_tau_);
  return(ncnt_nd);
}

Uvar L1A::normalizeCntRadio()
{
  Uvar ncnt_radio=((cnt_radio+radiometer_offset_*rad)/
		   (rip+radiometer_delta_tau_))/rad;
  return(ncnt_radio);
}




