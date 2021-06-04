#include "RadiometerData.h"
#include "L1B.h"

//---------------------------
// Methods for RadiometerData
//---------------------------

using std::cerr;

//--------------
// Constructors
//--------------

RadiometerData::RadiometerData(L1B& l1b, Time t1, Time t2) throw(ErrorMessage)
  : cnt_rl("cnt_rl",1),
    cnt_nd("cnt_nd",1),
    cnt_radio("cnt_radio",1),
    rad("rad",1),
    csr("csr",1),
    r_mode("r_mode",1),
    ctbe("ctbe",1),
    sclk("sclk",1),
    brst("brst",1),
    hip("hip",1),
    cip("cip",1),
    rip("rip",1),
    bpd("bpd",1),
    pul("pul",1),
    tro("tro",1),
    rwd("rwd",1),
    pri("pri",1),
    be1tmp("be1tmp",1),
    be2tmp("be2tmp",1),
    be3tmp("be3tmp",1),
    be4tmp("be4tmp",1),
    be5tmp("be5tmp",1),
    wgb1t1("wgb1t1",1),
    wgb3t1("wgb3t1",1),
    wgb3t2("wgb3t2",1),
    wgb3t3("wgb3t3",1),
    wgb5t1("wgb5t1",1),
    mratmp("mratmp",1),
    mruttm("mruttm",1),
    nsdtmp("nsdtmp",1),
    rlotmp("rlotmp",1),
    record_count(0),
    epoch(),
    time("time",1),
    beam_number("beam_number",1),
    ncnt_rl("ncnt_rl",1),
    ncnt_nd("ncnt_nd",1),
    ncnt_radio("ncnt_radio",1),
    misc_delay(3,"ms"),
    delta_tau(0,"ms"),
    offset(3550)
  {

  //---------------------------------
  // Read past L1B Header
  //---------------------------------

  l1b.readHeader();
  //---------------------------------
  // Setup data arrays.
  //---------------------------------

  unsigned int start_record=0;
  unsigned int end_record=0;
  if(t1.valid())
    {
      start_record = l1b.sclkToRecord(t1.sclk("Cassini"));
    }
  if(t2.valid())
    {
      end_record = l1b.sclkToRecord(t2.sclk("Cassini"));
    }
  else
    {
      end_record=l1b.recordCount();
    }

  record_count=end_record-start_record;
  if(record_count==0)
    {
      throw ErrorMessage("RadiometerData: No data on time range");
    }

 if(end_record<start_record)
    {
      throw ErrorMessage("RadiometerData: Time Range Misspecified t2<t1");
    }


  cnt_rl.resize(record_count);
  cnt_nd.resize(record_count);
  cnt_radio.resize(record_count);
  rad.resize(record_count);
  csr.resize(record_count);
  r_mode.resize(record_count);
  ctbe.resize(record_count);
  sclk.resize(record_count);
  brst.resize(record_count);
  hip.resize(record_count);
  cip.resize(record_count);
  rip.resize(record_count);
  bpd.resize(record_count);
  pul.resize(record_count);
  tro.resize(record_count);
  rwd.resize(record_count);
  pri.resize(record_count);

  be1tmp.resize(record_count);
  be2tmp.resize(record_count);
  be3tmp.resize(record_count);
  be4tmp.resize(record_count);
  be5tmp.resize(record_count);
  wgb1t1.resize(record_count);
  wgb3t1.resize(record_count);
  wgb3t2.resize(record_count);
  wgb3t3.resize(record_count);
  wgb5t1.resize(record_count);
  mratmp.resize(record_count);
  mruttm.resize(record_count);
  nsdtmp.resize(record_count);
  rlotmp.resize(record_count);


  time.resize(record_count);
  beam_number.resize(record_count);
  ncnt_rl.resize(record_count);
  ncnt_nd.resize(record_count);
  ncnt_radio.resize(record_count);


  //---------------------------------
  // Load L1B records one by one.
  //---------------------------------

  for (unsigned int i=0; i < start_record; ++i)
    {
      l1b.readSclkOnly();
    }

  for (unsigned int i=0; i < record_count; ++i)
    {
    l1b.readRecord();

    cnt_radio(i) = l1b.cnt_radio;
    cnt_rl(i) = l1b.cnt_rl;
    cnt_nd(i) = l1b.cnt_nd;
    rad(i) = l1b.rad;
    csr(i) = l1b.csr;
    r_mode(i) = l1b.r_mode;
    ctbe(i) = l1b.ctbe;

    // l1b.sclk is along int directly from telemetry
    // The LSB is 1.0 s.
    sclk(i) = l1b.sclk;

    brst(i) = l1b.brst;
    hip(i) = l1b.hip;
    cip(i) = l1b.cip;
    rip(i) = l1b.rip;
    bpd(i) = l1b.bpd;
    pul(i) = l1b.pul;
    tro(i) = l1b.tro;
    rwd(i) = l1b.rwd;
    pri(i) = l1b.pri;

    be1tmp(i)=l1b.be1tmp;
    be2tmp(i)=l1b.be2tmp;
    be3tmp(i)=l1b.be3tmp;
    be4tmp(i)=l1b.be4tmp;
    be5tmp(i)=l1b.be5tmp;

    wgb1t1(i)=l1b.wgb1t1;
    wgb3t1(i)=l1b.wgb3t1;
    wgb3t2(i)=l1b.wgb3t2;
    wgb3t3(i)=l1b.wgb3t3;
    wgb5t1(i)=l1b.wgb5t1;

    mratmp(i)=l1b.mratmp;
    mruttm(i)=l1b.mruttm;
    nsdtmp(i)=l1b.nsdtmp;
    rlotmp(i)=l1b.rlotmp;
    }


  processTime();

  extractBeamNumbers();

  normalize();
  }

void RadiometerData::processTime(){
  epoch = Time("Cassini",sclk(0));   
  for(unsigned int i=0;i<record_count;i++){  
    time(i) = Time("Cassini",sclk(i)) + Time(brst(i));
    Uvar receive_window=pul(i) * pri(i) +tro(i);
    Uvar time_offset=misc_delay+rwd(i)+receive_window+
      hip(i) + cip(i) + (rip(i) * rad(i))/2;
    time(i)=time(i)+Time(time_offset);
  }
}

void RadiometerData::extractBeamNumbers(){
  for(unsigned int i=0;i<record_count;i++){
    int i_ctbe=(int) ctbe(i).getValue();
    int j=0;
    while(i_ctbe){
      i_ctbe>>=1;
      j++;
    }
    beam_number(i)=j;
  }
}

void RadiometerData::normalize(){

  ncnt_rl=(cnt_rl+offset)/(cip+delta_tau);
  ncnt_nd=(cnt_nd+offset)/(hip+delta_tau);
  ncnt_radio=((cnt_radio+offset*rad)/(rip+delta_tau))/rad;

}














