//----------------------------------------------------------------------------
//
// RasList.cpp
//
// handle intermediate file format for SS's RAS tools
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------


//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>


//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <map>
#include "Error.h"
#include "Units.h"
#include "Ieb.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "RasList.h"
#include "Plot.h"
#include "SARFunctions.h"
#include "Baq.h"
#include "Frame.h"
#include "Utils.h"
using std::cout;
using std::endl;
using std::cerr;

//----------------------------------------
// Constructors: contruct with cfg only
//-------------------------------------
RasList::RasList(const string& filename, const string& mode)
  :echo_data("echo"),
   fft_echo_data("fft echo data"),
   abs_echo_data("absolute value of echo "),
   baq_decoded("baq decoded"),
   sign_bit("sign bit"),
   echo2D("echo2D"),
   echo2D_c("echo2D_c"),
   file_(filename,mode),
   filename_(filename),  
   mode_(mode),
   decoded_(false),
   echo2D_set_(false)
  
  {   
    ras_sab_size_before_buffer_ = sizeof(unsigned int) + sizeof(unsigned int);
    ras_sab_size_before_buffer_ += sizeof(double)+5*sizeof(unsigned int);
    ras_sab_size_before_buffer_ += 12*sizeof(unsigned short);
    ras_sab_size_before_buffer_ += 2*sizeof(unsigned int);

    ras_sab_size_after_sab_before_buffer_=sizeof(double)+5*sizeof(unsigned int);
    ras_sab_size_after_sab_before_buffer_+=12*sizeof(unsigned short);
    ras_sab_size_after_sab_before_buffer_+=2*sizeof(unsigned int);
   
    ras_up_to_sab_size_ = sizeof(unsigned int) + sizeof(unsigned int);
    
    baq_decoded.resize(Nblock_);//24 block
    sab_counter = 1;//initial value
    if(mode=="rb" || mode=="r")
      mapAllRecords();
  }

RasList::~RasList()
{}
//-----------------------------------------
//load sab
//-----------------------------------------
void RasList::loadSab(const  Sab& sab)
{
  //write only mode
  if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
  //import data from sab    
  t.setSclk("Cassini",sab.sclk);
  t+= sab.brst;
  
  sync=sab.sync;
  sab_counter = sab.sab_counter;
  sclk=sab.sclk;
  brst = sab.int_brst;
  int i_ctbe= sab.ctbe;
  beam_number = 0;
  while(i_ctbe)
    {
      i_ctbe>>=1;
      beam_number++;
    }
  baq_mode=sab.baq_mode;
  radar_mode=sab.r_mode;
 
  Ieb ieb;
  ieb=sab.getIeb();
  valid_data_per_pri=ieb.getPointsInsidePRI();
  //----------------------------------------------------
  //  valid_data_per_burst=ieb.getRadarDataPoints();
  // There was a case in which pul+tro was not an actual window count
  // because of bad command.  In that case, sab.ctrx is the true value
  // of radar window count
  //----------------------------------------------------------
  valid_data_per_burst = ieb.getRadarDataPoints();
  sab.extract_baq_raw_data(baq_,rdata_);
  radar_data_size_in_bytes = rdata_.size();
  writeRecord();
}

//----------------------------
//load simulated data
//---------------------------
void RasList::loadEncodedData(Ieb& ieb, 
			      const Time& transmit_time,
			      const unsigned int& burst_id,
			      const unsigned int& beam_id,
			      const Array1D<unsigned short>& baq,
			      const Array1D<unsigned char>& data)
{ 
  //write only mode
  if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
  //import data from simulation
  //(0) sync 
  sync=2004118378;
  //(1) sab counter
  sab_counter = burst_id;
  //(2) ephemeris time
  t = transmit_time;
  if(t < ieb.getTime()) ErrorMessage("RasList.cpp:: transmit time is earlier than valid ieb time").throwMe();
  //(3) sclk
  sclk = t.sclk("Cassini");
  //(4) brst
  Time t2;
  t2.setSclk("Cassini",sclk);
  Uvar brst_time = t - t2;
  if(brst_time<Uvar(0,"s")) brst_time=Uvar(0,"s");
  brst =(unsigned int) (brst_time.getInUnits("s")/(0.25*0.001));
  //(5) beam number    
  beam_number = beam_id;
  //(6) baq mode
  baq_mode=ieb.getBaq();
  //(7) radar mode
  radar_mode = ieb.getMod();
  //(8) baq threshold
  baq_.resize(baq.size());
  if(baq.size()!=12) ErrorMessage("Baq size mismatch").throwMe();
  for(unsigned int i=0;i<baq.size();++i) baq_[i]=baq(i);
 
  //(9) valid data per pri
  valid_data_per_pri=ieb.getPointsInsidePRI();
  //(10) valid data per burst
  valid_data_per_burst=ieb.getRadarDataPoints();
  //(11)
  radar_data_size_in_bytes= (ieb.getNumberOfWords()-ieb.Nheader-ieb.Nfooter)*2;
  //(12) data
  if(radar_data_size_in_bytes!=data.size()) ErrorMessage("encoded data size doe not match to ieb's expected data size").throwMe();
  rdata_.clear();
  for(unsigned int i=0;i<data.size();++i) rdata_.push_back(data(i));
   writeRecord();
}

//---------------------------
//load ieb
//---------------------------
void RasList::decodeEcho(Ieb& ieb)
 {
   if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
   if(ieb.getTime() > t) ErrorMessage("RasList.cpp:: ieb time is later than burst time").throwMe();
   ieb.decodeSlowFastfield();
   
   //decode
   //copied from sab decode
    pos=0;
    neg=0;
    if(radar_data_size_in_bytes==0) return;//no data to decode
    //--------------
    //decode baq
    //--------------
    for (unsigned int  loop = 0; loop < 12; loop++)
      {	//byte swapped reading
	baq_decoded(2 * loop+1 ) = bitget(baq_[loop],0,7);  
	baq_decoded(2 * loop) = bitget( baq_[loop],8,15);	   	
      }
    decodeCassiniData(rdata_, 
		      baq_decoded,
		      ieb, 
		      echo_data,
		      valid_data_per_burst,offset);
    decoded_ = true; 	    


    //debug purpose
    //Uvec y1("",echo_data.size());
    //y1=echo_data;// original echo
    //
    //Baq baq;
    //
    //baq.setParams( ieb.getPri(), ieb.getPul(),ieb.getBaq(), ieb.getAdc(), ieb.getTroInPriUnits());
    //
    //Charvec words("",rdata_.size());
    //for(unsigned int i=0;i<rdata_.size();++i)
    //words(i)=rdata_[i];//transfer coded word 
    //
    //Array1D<float> echo_data2("");
    //
    //store baq into integer array
    //Ivec sab_baq("",Nblock_);
    //sab_baq=baq_decoded;
    //Array1D<unsigned short> thresh_encoded("",Nblock_);
    //baq.Encode_thresh(sab_baq,thresh_encoded);
    //
    //Fvec fout("");
    //baq.Decode_Nbit(words,sab_baq,fout);
    //
    //Uvec y3("",fout.size());
    //y3=fout;
    //
    //
    //
    //
    //plot data
    //Plot a;
    //a.addY(y1,"",line("solid","black",2),sym("none"));
    //
    //a.addY(y3,"",line("solid","magenta",2),sym("none"));
    //a.setTitle("original-black, baq encoded and decoded-red, baq_encode_decode-mag");
    //a.show("x");
 }

void RasList::decodeBaqThresholds()
 {
   if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::decodeBaqThresholds works only reading mode").throwMe();
 
   //--------------
   //decode baq
   //--------------
   for (unsigned int  loop = 0; loop < 12; loop++)
     {	//byte swapped reading
       baq_decoded(2 * loop+1 ) = bitget(baq_[loop],0,7);  
       baq_decoded(2 * loop) = bitget( baq_[loop],8,15);	   	
     }
 }


//-----------------------------
//find and read Record(const int n_sab)
//------------------------------
bool RasList::findRecord(const unsigned int& n_sab)
{//read record whose sab number is n_sab
  if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
  bool foundRecord= false;
  sabnumber_vs_fileposition_ptr_= sabnumber_vs_fileposition_.find(n_sab);
  if(sabnumber_vs_fileposition_ptr_ != sabnumber_vs_fileposition_.end())
    {
      foundRecord = true;
      file_.setPosition(  sabnumber_vs_fileposition_ptr_ -> second);
      readRecord();
    }
  return(foundRecord);
}

//-----------------------------------
//find record using time
//-----------------------------------
bool RasList::findRecord(const Time& t)
{
  if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
  bool foundRecord= false;
  bool readRecord= false;
  time_vs_sabnumber_ptr_= time_vs_sabnumber_.find(t);
  if(time_vs_sabnumber_ptr_ != time_vs_sabnumber_.end())
    {
      foundRecord= true;
      readRecord= findRecord( time_vs_sabnumber_ptr_ ->second);
      if(!readRecord) ErrorMessage("Possible mapping error  ").throwMe();
    }
  return(foundRecord);
}
unsigned int& RasList::operator[] (const unsigned int& n)//return sab number of nth recor
 {
   if(n > sabnumber_list_.size()) ErrorMessage("RasList.cpp:: no record ").throwMe();
   return(sabnumber_list_[n]);
 }

//----------------------------------
//read first record
//----------------------------------
void RasList::readFirstRecord()
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    file_.setPosition(0);//set to beginning
    readRecord();
  }

//--------------------
//read next record
//--------------------
void RasList::readNextRecord()
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    if(file_.eof())
      {
	ErrorMessage("no more record").throwMe();
      }
    else
      {
	readRecord();
      }
  }
//------------------------------------
//skipOneRecord()
//------------------------------------
void RasList::skipOneRecord()
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    //current record is stored in sab_counter
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    if(file_.eof()) ErrorMessage("RasList::skipOneRecord(): no more record").throwMe();
    //skip one record
    int start = file_.getPosition();
    file_.setPosition(start+int(ras_sab_size_before_buffer_));
    file_.read(radar_data_size_in_bytes);
    start = file_.getPosition();
    file_.setPosition(start+int(radar_data_size_in_bytes*sizeof(unsigned char)));
  }


//----------------------------------
//computeRecordNumbers()
//-----------------------------------
unsigned int RasList::computeNumberOfRecords()
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    return(sabnumber_list_.size());
  }

//------------------------
//compute first and last sab number
//------------------------
 void RasList::getFirstLastSabNumber(unsigned int& first_sab, unsigned int& last_sab)
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    first_sab = sabnumber_list_.front();//first element
    last_sab= sabnumber_list_.back();
  }


//-----------------
//rewind
//---------------
void RasList::rewind()
  {
    file_.setPosition(0);
  }

//--------------------
//EOF
//--------------------
bool RasList::eof()
  {
    return(file_.eof());
  }




//---------------------------------------------------------------
// Method Name: shiftEcho
// Author: YZ & YG
// Function Description:
// Shift 1D burst echo (echo_data) by dt, so that rwd starts at
// the beginning of the first pulse.
// Notes:
//---------------------------------------------------------------
//void RasList::shiftEcho(const unsigned int& dn)
//{
//  if(!decoded_) ErrorMessage("RasList.cpp:: echo is not decoded").throwMe();
//  Array1D<float> shifted_echo_data("se",echo_data.size());
//  shifted_echo_data = 0.0;//default
//  for(unsigned int i=0;i<echo_data.size();++i) 
//    {
//     	if(i<dn) continue;//go to next i
//unsigned int ip =(unsigned int) (i - dn);
//shifted_echo_data(ip)=echo_data(i);
//    }
//  echo_data = shifted_echo_data;
//}



//-----------------------------------------------------------------------
// return decoded data
//------------------------------------------------------------------------
void RasList::getEchoData(Array1D<float>& echo1D)
{
  if(!decoded_) ErrorMessage("RasList.cpp:: echo is not decoded").throwMe();
  unsigned int N = echo_data.size();
  echo1D.resize(N);  
  echo1D = echo_data;
}


//------------------------------------------------------------------------
// Method Name: get2Decho
// Author: YZ
// Function Description:
// Change 1D data array into 2D matrix with the size of the power of 2 for
// both dimensions.
// Notes:
//------------------------------------------------------------------------
void RasList::get2Decho(const unsigned int& start_index,
			 const unsigned int& end_index,
			 const unsigned int& N_samples_per_window,
			 const unsigned int& N_pulses_inside_echo)
{
  if(!decoded_) ErrorMessage("RasList.cpp:: echo is not decoded").throwMe();
  if(end_index>echo_data.size()) ErrorMessage("end index is larger than data size").throwMe();

   Convert2DEcho(echo_data, start_index,
		     end_index, 
		     N_samples_per_window,
		     N_pulses_inside_echo,
		     valid_data_per_pri,
		     echo2D);
   Convert2DReal2Complx( echo2D,  echo2D_c,-1);//choose negative bandwidth
  
  /*
  Ndop_bins_ = (unsigned int)(N_pulses_inside_echo);
  Convert2Pow(Ndop_bins_);
    
  Nrange_bins_ = (unsigned int)(N_samples_per_window);
  Convert2Pow(Nrange_bins_);

  echo2D.resize(Ndop_bins_,Nrange_bins_);
  echo2D=0.0;//0
  
  
  
  for(unsigned int i=0; i < N_pulses_inside_echo; i++) 
    for(unsigned j=0; j < N_samples_per_window; j++){
      unsigned int n = i*valid_data_per_pri + j+start_index;
      if(n<end_index)
	echo2D(i,j) = echo_data(n); 
    }
  getComplex2Decho();   
  */

  echo2D_set_=true;
  //debug
  //cout<<"Nrange and Ndop "<< Nrange_bins_<<" "<< Ndop_bins_<<endl;
}

//-----------------
//direct fft
//-----------------
void RasList::directFFT()
 {
   if(!decoded_) ErrorMessage("data not decoded").throwMe();
   unsigned int N= echo_data.size();
   unsigned int M=   get_next_power_of_2(N);
   CFvec  echo_data_tmp("",M);
   for(unsigned int i=0;i<M;++i){
     if(i<N) echo_data_tmp(i)=complex<float>(echo_data(i),0);
     else echo_data_tmp(i)=complex<float>(0,0);
   }
   fft_echo_data.resize(M);
   fft(echo_data_tmp, fft_echo_data);
   //scaling
   for(unsigned int i=0;i<M;++i) fft_echo_data(i)/=float(N);

   Plot a;
 
   Uvec y("",M/2);
   for(unsigned int i=0;i<M/2;++i){
     y(i)=norm(fft_echo_data(i));
   }
   a.addY(y,"",line("solid","red",2),sym("none"));
   a.show("x");
 }
//----------------------------
//display echo
//----------------------------
void RasList::displayEcho()
{
  if(!decoded_) ErrorMessage("data not decoded").throwMe();
  Uvec a_list("a",Nblock_);
  for(unsigned int i=0;i<a_list.size();++i) a_list(i)= baq_decoded(i);
  Uvec b_list("b",echo_data.size());
  for(unsigned int i=0;i<b_list.size();++i) b_list(i)=echo_data(i);
 

  cout<<"SAB information "<<endl;
  cout<<"sab number "<<sab_counter<<endl;
  cout<<"beam number "<<beam_number<<endl;
  cout<<"baq mode "<<baq_mode<<endl;
  cout<<"radar mode "<<radar_mode<<endl;
  cout<<"data inside pri "<< valid_data_per_pri<<endl;
  cout<<"data inside burst "<<valid_data_per_burst<<endl;
  cout<<"radar data size in bytes "<< radar_data_size_in_bytes<<endl;
  cout<<"offset "<<offset<<endl;

  Plot a;
  a.addY(a_list,"",line("solid","red",2),sym("none"));
  a.setTitle("loaded BAQ threshold ");
  a.setSubTitle("Sab counter:"+toStr(sab_counter));
  a.show("x");

  Plot b;
  b.addY(b_list,"",line("solid","red",2),sym("none"));
  b.setTitle(" loaded echo sab number:" +toStr(sab_counter));
  b.setSubTitle("baq: "+toStr(baq_mode));
  b.show("x");    

}

//---------------------------
//compare baq exponent: downlinked one and theoretical one
//---------------------------
void RasList::compare_baq_exponent(const Uvar& pri, const unsigned int& pul,
				   const unsigned int& tmp_baq,
				   const Uvar& adc,
				   const int& tro_int)
{
  if(!decoded_) ErrorMessage("data not decoded").throwMe();

  //----
  for(unsigned int i=0;i<echo_data.size();++i){
    if(echo_data(i)>0.0) echo_data(i)-=0.5;
    else echo_data(i) += 0.5;
  }

  //store SAB footter baq values
  Uvec a_list("a",Nblock_);
  for(unsigned int i=0;i<a_list.size();++i) a_list(i)= baq_decoded(i);
  
  //set baq parameter
  Baq baq;
  baq.setParams( pri, pul, tmp_baq, adc, tro_int);

  //store data into character array
  Charvec echo_char("");
  echo_char=baq.ADC_without_normalization(echo_data);
  
  //store baq exponent into integer array
  cout<<"compute threshold using baq mode of  "<< tmp_baq<<endl;
  Ivec threshold("",Nblock_);
  baq.compuThreshold(echo_char, threshold);
  
  //store baq threshold
  Uvec b_list("a",Nblock_);
  for(unsigned int i=0;i<a_list.size();++i) b_list(i)= threshold(i);

  //print out i_block_th baq exponents: sab footer and computed
  int i_block=12;
  cout<<"block number "<<i_block<<endl;
  cout<<"from downlink, baq threshold "<< baq_decoded(i_block)<<endl;
  cout<<"baq method "<< threshold(i_block)<<endl;
  
  //debug: print on the screen
  int Npri= round_double( (adc*pri).getInUnits(""));
  int rcv_window = (unsigned int) ( (int)pul + tro_int);
  cout<<"Npri and rev "<< Npri<<" "<< rcv_window<<endl;
  int Necho = Npri*rcv_window;
  int Nblock=24;
  int Ns = int(Npri/Nblock);
  cout<<"Ns : points inside block "<< Ns<<endl;
  cout<<"total number of data "<<Necho<<endl;

  
 
  //store echo into 2D array form: number_of_points_inside_pri * number_of_rcv_window
  // in units of pri
  Array2D<float> echo2D("",Npri,rcv_window);
  for(int i=0;i<rcv_window;++i)
    for( int j=0;j<Npri;++j){
      echo2D(j,i)= int(echo_data(j + i*Npri));
    }

  //compute average of absolute magnitude of echoes
  float sum=0.0;
  int index=0;
  int Npoints_inside_block;
  int Npulse_inside_echo;
  int avg_factor;
  if(tmp_baq==0){
    Npoints_inside_block=8;
    Npulse_inside_echo=8;
    avg_factor=32;
  }
  else if(tmp_baq==6){
    Npoints_inside_block=4;
    Npulse_inside_echo=8;
    avg_factor=16;
  }
  else if(tmp_baq==7){
    Npoints_inside_block=8;
    Npulse_inside_echo=4;
    avg_factor=16;
  }
  else ErrorMessage("no baq").throwMe();
  
  //do echo summation: beginning of echoes
  for( int j=0;j<Npulse_inside_echo;++j){
    for(int i=0;i<Npoints_inside_block;++i){
      //cout<<"index "<< index<<" "<< echo2D(i+i_block*Ns,j)<<endl;
      index++;
      sum += fabs(echo2D(i+i_block*Ns,j));
    }
    for(int i=Ns-Npoints_inside_block;i<Ns;++i){
      //cout<<"index "<< index<<" "<< echo2D(i+i_block*Ns,j)<<endl;
      index++;
      sum += fabs(echo2D(i+i_block*Ns,j));
    }
  }

  //do summation: end of echoes
  for( int j=rcv_window-Npulse_inside_echo;j<rcv_window;++j){
    for(int i=0;i<Npoints_inside_block;++i){
      //cout<<"index "<< index<<" "<< echo2D(i+i_block*Ns,j)<<endl;
      index++;
      sum += fabs(echo2D(i+i_block*Ns,j));
    }
    for(int i=Ns-Npoints_inside_block;i<Ns;++i){
      //cout<<"index "<< index<<" "<< echo2D(i+i_block*Ns,j)<<endl;
      index++;
      sum += fabs(echo2D(i+i_block*Ns,j));
    }
  }

  cout<<"total sum "<< sum<<endl;
  cout<<"look up table "<< sum/avg_factor<<endl;
  
  
  
  //plot loaded baq threshold and computed baq threshold
  Plot a;
  a.addY(a_list,"",line("solid","red",2),sym("none"));
  a.addY(b_list,"",line("solid","black",2),sym("none"));
  a.setTitle("loaded BAQ threshold: red(downlink) black(computed) ");
  a.setSubTitle("Sab counter:"+toStr(sab_counter));
  a.show("x");


  cout<<"baq threshold: down "<< baq_decoded<<endl;
  cout<<"baq trheshold: baq  "<<threshold<<endl;
  //plot loaded baq and computed baq with echoes used to compute baqs
  //add original echo
  Plot b;
  Uvec x("",Npri);
  Uvec y("",Npri);
  for(int i=0;i<Npulse_inside_echo;++i){
    for(int j=0;j<Npri;++j){
      x(j)=j;
      y(j)=fabs(echo2D(j,i));
    }
    b.addXY(x,"",y,"",line("solid","black"),sym("none"));
  }

  for( int i=rcv_window-Npulse_inside_echo;i<rcv_window;++i){
    for(int j=0;j<Npri;++j){
      x(j)=j;
      y(j)=fabs(echo2D(j,i));
    }
    b.addXY(x,"",y,"",line("solid","black"),sym("none"));
  }
  //add original baq
  for( int j=0;j<Npri;++j){
    x(j)=j;
    int k= int(j/Ns);
    if(k>=Nblock) k=Nblock-1;
    //cout<<"j and k "<<j<<" "<<k<<endl;
    y(j)=baq_decoded(k);
  }
  b.addXY(x,"",y,"",line("solid","red"),sym("none"));
  //add computed baq
  for( int j=0;j<Npri;++j){
    x(j)=j;
    int k= int(j/Ns);
    if(k>=Nblock) k=Nblock-1;
    y(j)=threshold(k);
  }
  b.addXY(x,"",y,"",line("solid","blue"),sym("none"));
  b.setTitle("echo-black, blue-compued baq, red-sab header baq");
  b.show("x");
  


 
 
  Charvec words("words");//encoded value
  Array1D<float> f_out("");
  baq.setParams( pri, pul, tmp_baq, adc, tro_int);
  cout<<"echo char size "<< echo_char.size()<<endl;
  cout<<"encode echo using baq mode of "<< tmp_baq<<endl;
  baq.Encode_Nbit(echo_char, threshold, words);//encode echo using computed threshold into words
  cout<<"encoded word size "<< words.size()<<endl;
  baq.Decode_Nbit(words,threshold,f_out);//decode back into float array
  cout<<"decoded size "<<f_out.size()<<endl;
  Uvec original_echo("",echo_data.size());
  Uvec decoded_echo("",f_out.size());
  original_echo= echo_data;
  decoded_echo=f_out;

  Uvec decoded_sab_baq("",f_out.size());
  Ivec sab_baq("",Nblock_);
  sab_baq=baq_decoded;//store original baq exponent into integer array
  
  
 
  f_out=0.0;
  Baq baq2;
  unsigned int baq_used=tmp_baq;
  if(baq_mode==5) baq_used=6;
  baq2.setParams( pri, pul, baq_used ,adc, tro_int);//original decoding
  Charvec echo_char2=baq.ADC_without_normalization(echo_data);
  Charvec words2("words");//encoded value
  cout<<"echo char size "<< echo_char2.size()<<endl;
  cout<<"encode with baq mode of "<< baq_used<<endl;
  baq2.Encode_Nbit(echo_char2, sab_baq, words2);//encode using original baq exp
  cout<<"encoded word size "<< words2.size()<<endl;
  baq2.Decode_Nbit(words2,sab_baq,f_out);//decode back
  cout<<"decoded size "<<f_out.size()<<endl;
  decoded_sab_baq = f_out;
  //plot original data and decoded data
  Plot c;
  c.addY(original_echo,"",line("solid","black",2),sym("none"));
  c.addY(decoded_echo,"",line("solid","red",2),sym("none"));
  c.addY(decoded_sab_baq,"",line("solid","magenta",2),sym("none"));
  c.setTitle("orginal data(black),first_encoded-and-then_decoded (red), using original baq(mag)");
  c.show("x");
}



//----------------------------
// plot 2D sampled echo
//----------------------------
void RasList::display2Decho(const unsigned int& n_pulse)
{
  if(!echo2D_set_) ErrorMessage("RasList.cpp:: echo2D is not set").throwMe();
  unsigned int M, N;
  echo2D.size(M,N);
  if(n_pulse ==0) ErrorMessage("RasList.cpp: pulse starts from 1").throwMe();
  if(n_pulse>(M+1)) ErrorMessage("RasList.cpp:: max pulse is "+toStr(M+1)).throwMe();
  unsigned int ii = n_pulse-1;
  Uvec y("",echo2D.cols());
  y = echo2D.getRow(ii);
  Plot a;
  a.addY(y,"",line("solid","red",2),sym("none"));
  a.setTitle("pulse number "+toStr(n_pulse));
  a.show("x");

}
//display all echoes
void RasList::displayAll2Dechoes()
  {
    if(!echo2D_set_) ErrorMessage("RasList.cpp:: echo2D is not set").throwMe();
    unsigned int M, N;
    echo2D.size(M,N);
    Uvec y("",echo2D.cols());
    Plot a;

    for(unsigned int ii=0;ii<M;ii++){
      y=echo2D.getRow(ii);
      y += float(125.0*(float)ii);
      a.addY(y,"",line("solid","red",2),sym("none"));
    }
    a.setTitle("plot all echoes");
    a.show("x");
  }


//----------------------------
// plot complex pulse
//----------------------------
void RasList::displayComplexPulse(const unsigned int& n_pulse)
  {
    if(!echo2D_set_) ErrorMessage("RasList.cpp:: echo2D is not set").throwMe();
    unsigned int M, N;
    echo2D_c.size(M,N);
    if(n_pulse ==0) ErrorMessage("RasList.cpp: pulse starts from 1").throwMe();
    if(n_pulse>(M+1)) ErrorMessage("RasList.cpp:: max pulse is "+toStr(M+1)).throwMe();
    unsigned int ii = n_pulse-1;

   
    Uvec data_real("data_real",N);
    Uvec data_imag("data_imag",N);
    for(unsigned int jj = 0;jj < N; jj++){
      data_real(jj) = real(echo2D_c(ii,jj));
      data_imag(jj) = imag(echo2D_c(ii,jj));
    }
    
    Plot a;
    a.addY(data_real,"",line("solid","black",1),sym("none"));
    a.addY(data_imag,"",line("solid","red",1),sym("none"));
    a.setTitle("Complex Echo");
    a.show("x");    

  
    
    Uvec complex_real("",2*N);
    complex_real=0.0;
    Uvec original_echo("original echo",2*N);
    for(unsigned int jj = 0;jj < 2*N; jj++)
      original_echo(jj)= echo2D(n_pulse,jj);
    for(unsigned int jj=0;jj<N;++jj)
      complex_real(2*jj)= real( echo2D_c(n_pulse,jj));
    
    Plot b;
    b.addY(complex_real,"",line("none"),sym("circle","red",1));
    b.addY(original_echo,"",line("solid","red",1),sym("none"));
    b.setTitle("Complex and original Echo");
    b.show("x");    
    
  }



//---------------------
//show sign bit statistics
//---------------------
void RasList::showSignBitStatistics(const unsigned int& n_sab)
  {
    if(findRecord(n_sab))
      {
	Plot a;
	a.addY(sign_bit,"",line("none"),sym("circle","red",1));
	a.setTitle("burst number "+toStr(n_sab)+" "+"pos: "+toStr(pos)+" neg: "+toStr(neg));
	a.show("x");
      }
    else
      {
	cout<<"No Record to Plot "<<endl;
      }
  }

//------------------------------
//get data statistics
//-------------------------------
 void RasList::getEchoStatistics(double& standard_deviation,
			   double& average)
  {
    if(!decoded_) ErrorMessage("data is not decoded").throwMe();
    double sum=0.0;
    unsigned int N=echo_data.size();
    for(unsigned int i=0;i<N;++i)
      sum +=double( echo_data(i));
    
    average = sum/double(N);
    double dev2=0.0;
    for(unsigned int i=0;i<echo_data.size();++i)
      dev2 += (double(echo_data(i))-average)*(double(echo_data(i))-average)/double(N-1);
    
    standard_deviation = sqrt(dev2);
  }
			   
//-------------------------------------
//Private Function Declarations
//------------------------------------

//------------------------------------------------
//mapAllRecords(): 
// 
// these containers are set onece when rasfile is opened
//    vector<unsigned int> sabnumber_list_;
//   map<unsigned int, int> sabnumber_vs_fileposition_;
//  vector<Time> time_list;
//  map<Time, unsigned int> time_vs_sabnumber_;
//---------------------------------------------------
void RasList::mapAllRecords()
  {
    sabnumber_list_.clear();
    sabnumber_vs_fileposition_.clear();
    time_list_.clear();
    time_vs_sabnumber_.clear();
    cout<<"start Mapping "<<endl;
    //extract file pointer, sab number and time information from all the records
    if(mode_=="wb"||mode_=="w") ErrorMessage("file mode is not for reading").throwMe();
    //rewind first
    int dumm=0;
    file_.rewind();
    while(!file_.eof())
      {
	//save file pointer position of the record before read
	dumm = file_.getPosition();
	readRecord();
	//(1) sab number
	sabnumber_list_.push_back(sab_counter);
	//(2) sab number vs file position
	sabnumber_vs_fileposition_[sab_counter]=dumm;
	//(3) time list
	time_list_.push_back(t);
	//(4) time vs sab number
	time_vs_sabnumber_[t]= sab_counter;
      }
    cout<<"All the maps are complete "<<endl;
  }


//-----------------------------------
//read next record
//-----------------------------------
void RasList::readRecord()
  {
    if(file_.eof()) ErrorMessage("RasList::readRecord(): no more record").throwMe();
    if(mode_=="wb"||mode_=="w") ErrorMessage("file mode is not for reading").throwMe();
    //(0) read sync
    file_.read(sync);
    if(sync!=2004118378) ErrorMessage("Incorrect sync").throwMe();

    //(1)sab counter
    file_.read(sab_counter); 

    //(2)ephemeris time in second
    double t_et;
    file_.read(t_et);
    t.setEt(Uvar(t_et,"s"));
    
    //(3)sclk : long unsigned int
    file_.read(sclk);//unsigned int
    
    //(4)brst : long unsigned int
    file_.read(brst);//unsigned int
    
    //(5)compute beam number
    file_.read(beam_number);//unsigned int
    
    //(6)baq mode
    file_.read(baq_mode);//unsigned int
    
    //(7)radar mode
    file_.read(radar_mode);//unsigned int
    
    baq_.clear();
    //(8)get baq threshold: twelve 16bit words
    unsigned short dummy;
    for(unsigned int i=0;i<Nblock_/2;++i)
      {
	file_.read(dummy);//unsigned short
	baq_.push_back(dummy);
      } 
    if(baq_.size()!=Nblock_/2) ErrorMessage("baq size error").throwMe();
    //(9) valid data number per pri
    file_.read(valid_data_per_pri);//unsigned int
    
    //(10) valid data number per burst
    file_.read(valid_data_per_burst);//unsigned int
    
    //(11)raw data buffer size
    file_.read(radar_data_size_in_bytes);//unsigned int
    
    //(12)write raw data
    unsigned char dumm2;
    rdata_.clear();
    if(radar_data_size_in_bytes!=0)
      {//there is some case when buffer size is zero while the 
	// radar is in active mode
	for(unsigned int i=0;i<radar_data_size_in_bytes;++i)
	  {
	    file_.read(dumm2);//unsigned char
	    rdata_.push_back(dumm2);
	  }
      }   
    if(rdata_.size()!=radar_data_size_in_bytes) ErrorMessage("radar size error").throwMe();
    //if read , decode
    //decodeEcho();    
    //set previously loaded ieb invalid
    decoded_ = false;
    echo2D_set_ = false;
  }

//------------------------
//{write) one record
//------------------------
void RasList::writeRecord()
  { 
    if(mode_=="rb"||mode_=="r") ErrorMessage("file mode is not for writing").throwMe();
    
    //(0)write sync value
    file_.write(sync);//unsigned int

    //(1)sab counter
    file_.write(sab_counter);//unsigned int
		
    //(2)ephemeris time in second
    file_.write(t.et().getInUnits("s"));//double
		
		
    //(3)sclk : long unsigned int
    file_.write(sclk);//unsigned int
		
    //(4)brst : long unsigned int
    file_.write(brst);//unsigned int
				
    //(5)compute beam number
    file_.write(beam_number);//unsigned int
		
    //(6)baq mode
    file_.write(baq_mode);//unsigned int
		
    //(7)radar mode
    file_.write(radar_mode);//unsigned int
		
    //(8)get baq threshold: twelve 16bit words
   
    unsigned short dumm=0;
    if(baq_.size()!=12) ErrorMessage("baq size is not 12 unsigned short").throwMe();
    for(unsigned int i=0;i<baq_.size();++i)
      {
	dumm = baq_[i];
	file_.write(dumm);//unsigned short
      } 
		
    //(9) valid data number per pri
    file_.write(valid_data_per_pri);//unsigned int
		
    //(10) valid data number per burst
    file_.write(valid_data_per_burst);//unsigned int
		
    //(11)raw data buffer size
    file_.write(radar_data_size_in_bytes);//unsigned int
		
    //(12)write raw data
    unsigned char  dumm2=0;
    if(radar_data_size_in_bytes!=rdata_.size()) ErrorMessage("buffer size mismatch").throwMe();
    if(radar_data_size_in_bytes!=0)
      {//there is some case when buffer size is zero while the 
	// radar is in active mode
	for(unsigned int i=0;i<rdata_.size();++i)
	  {
	    dumm2 = rdata_[i];
	    file_.write(dumm2);//unsigned char
	  }
      }
  }

//--------------------------------------------------------------------------
// Method Name: getComplex2Decho
// Author: YZ
// Function Description:
// Double sample the real echo2D and convert it into complex echo2D_c by using
// the same method as getComplex2Decho. The size of echo2D_c 
// is the same as echo2D.
// Notes:
//-------------------------------------------------------------------------
void RasList::getComplex2Decho()
{
  unsigned int M, N;
  echo2D.size(M,N);
  echo2D_c.resize(M,N/2);

  CFvec pulse_fft("pulse_fft",N);
  CFvec echo_cmplx("p_cmplx",N);
  CFvec pulse_fft_half("pulse_fft_half",N/2);
  CFvec cmplx_pulse("cmplx_pulse",N/2);

  for(unsigned int ii=0; ii < M; ii++)
    {
      pulse_fft = complex<float> (0.0,0.0);
      echo_cmplx= complex<float> (0.0,0.0);
      pulse_fft_half =  complex<float> (0.0,0.0);
      cmplx_pulse= complex<float> (0.0,0.0);

      //for each pulse
      for(unsigned int jj = 0; jj < N; jj++)
	echo_cmplx(jj) = complex<float>(echo2D(ii,jj),0.);
      
      //fft
      fft(echo_cmplx, pulse_fft);
    
      // select negative spectrum
      for(unsigned int n = N/2; n < N; n++)
	pulse_fft_half(n-N/2) = pulse_fft(n);

      // inverse fft 
      ifft(pulse_fft_half, cmplx_pulse);

      //save it!
      for(unsigned int n = 0; n < N/2; n++)
	echo2D_c(ii,n) = cmplx_pulse(n);
    }    
}



//-----------------------------------------------
//End of RASLIST.CPP
//------------------------------------------------


//------------------------------
//Originally developed by Yan Zhang to study doppler centroid
// most of them are moved to SARFunctions.h and cpp
//-------------------------------



//-----------------------------
// Method Name: convert2pow
// Author: YZ
// Function Description:
// Convert n to the power of 2
// Notes:
//-----------------------------
//void RasList::convert2pow(unsigned int& n)
//{
//unsigned int m = (unsigned int)(log((double)n)/log(2.));
//if(n == (unsigned int)pow(2,m))
//  n = (unsigned int)pow(2.,m);
//else
//  n = (unsigned int)pow(2.,m+1);
//}
//------------------------------------------------------------------------
// Method Name: avgPulseAmp
// Author: YZ
// Function Description:
// Get averaged amplitude of pulses 
// Notes:
//------------------------------------------------------------------------
/*
void RasList::avgPulseAmp()
{
  if(!echo2D_set_) ErrorMessage("RasList.cpp: echo2D is zero. ").throwMe();
  unsigned int M,N;
  echo2D.size(M,N);

  pulse_amp.resize(N);
  pulse_amp=0.0;//0
 
  for(unsigned int j=0; j < N; j++) 
  for(unsigned i=0; i < M; i++)
    {
      pulse_amp(j) = pulse_amp(j) + abs(echo2D(i,j))/M; 
    }

  pulse_amp_set_ = true;
  
  //plotAvgPulseAmp(); //Yan
}
*/


/*
void RasList::getComplex2Decho_d()
{

  unsigned int M, N;
  echo2D.size(M,N);
  echo2D_c.resize(M,N);
  echo2D_c = 0.0;

  for(unsigned int ii=0; ii < M; ii++)
    {
      CFvec echo_cmplx("p_cmplx",2*N);
      for(unsigned int jj = 0; jj < N; jj++)
	  echo_cmplx(2*jj) = complex<float>(echo2D(ii,jj),0.);
      for(unsigned int jj = 0; jj < N-1; jj++)
	  echo_cmplx(2*jj+1) = (echo_cmplx(2*jj) + echo_cmplx(2*jj+2))/2;

      CFvec pulse_fft("pulse_fft",2*N);
      fft(echo_cmplx, pulse_fft);

      // select negative spectrum
      CFvec pulse_fft_half("pulse_fft_half",N);
      CFvec cmplx_pulse("cmplx_pulse",N);
      for(unsigned int n = N; n < 2*N; n++)
	pulse_fft_half(n-N) = pulse_fft(n);

      // inverse fft 
      ifft(pulse_fft_half, cmplx_pulse);

      for(unsigned int n = 0; n < N; n++)
	echo2D_c(ii,n) = cmplx_pulse(n);
    }    

  //plotComplexPulse();
}
*/

//--------------------------------------------------------------------------------
// Method Name: pulseCorrelation
// Author: YZ
// Function Description:
// Calculate pulse-to-pulse correlation
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::pulseCorrelation()
{

  unsigned int M, N;
  echo2D_c.size(M,N);
  pulse_corr.resize(M,N);
  pulse_corr = 0.0;

  for(unsigned int ii=1; ii < M; ii++)
    for(unsigned int jj = 0; jj < N; jj++)
      pulse_corr(ii-1,jj) = echo2D_c(ii,jj)*conj(echo2D_c(ii-1,jj));
    
}
*/

//--------------------------------------------------------------------------------
// Method Name: getPhaseDiff
// Author: YZ
// Function Description:
// Calculate phase difference between neighboring pulses
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::getPhaseDiff()
{

  unsigned int M, N;
  pulse_corr.size(M,N);
  phase_diff.resize(M,N);
  phase_diff = 0.0;
  
  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
 
  Uvar Pri=ieb_.getPri();

  for(unsigned int ii=0; ii < M; ii++)
    for(unsigned int jj = 0; jj < N; jj++)
      phase_diff(ii,jj) = atan2(imag(pulse_corr(ii,jj)), real(pulse_corr(ii,jj)));
    
}
*/

//--------------------------------------------------------------------------------
// Method Name: get1DPhaseDiff
// Author: YZ
// Function Description:
// Calculate phase difference from averaged pulse correlation data
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::get1DPhaseDiff()
{

  unsigned int N = avg_pulse_corr.size();
  avg_phase_diff.resize(N);
  avg_phase_diff = 0.0;
  
  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
 
  Uvar Pri=ieb_.getPri();

    for(unsigned int jj = 0; jj < N; jj++)
      avg_phase_diff(jj) = atan2(imag(avg_pulse_corr(jj)), real(avg_pulse_corr(jj)));
    
}
*/

//--------------------------------------------------------------------------------
// Method Name: phaseUnwrapp
// Author: YZ
// Function Description:
// Phase Unwrapping
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::phaseUnwrapp()
{

  unsigned int M, N;
  phase_diff.size(M,N);
  Array2D<float> temp("temp",M,N);

  temp = phase_diff;

  float pi = 4.*atan(1.);

  for(unsigned int ii=0; ii < M; ii++)
    {
      unsigned int count = 0;

      for(unsigned int jj = 1; jj < N; jj++)
	{
	  if(phase_diff(ii,jj)-phase_diff(ii,jj-1) >= pi)
	    count ++;
	  if(phase_diff(ii,jj)-phase_diff(ii,jj-1) <= -pi)
	    count --;
	  
	  temp(ii,jj) = phase_diff(ii,jj) - count*2*pi;  
	}       
    }
  
  phase_diff = temp;
}
*/

//--------------------------------------------------------------------------------
// Method Name: averagePhaseDiff
// Author: YZ
// Function Description:
// Averaging phase difference, the first and last 5 pulses are not used.
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::averagePhaseDiff()
{

  unsigned int M, N;
  phase_diff.size(M,N);
  avg_phase_diff.resize(N);
  avg_phase_diff = 0.0;

  unsigned int num_pulse = (unsigned int)(valid_data_per_burst/valid_data_per_pri);
  unsigned int dM = 5;

  //cout<<"num_pulse"<<num_pulse<<endl;

  for(unsigned int jj = 0; jj < N; jj++)
    {
      avg_phase_diff(jj) = 0;

      for(unsigned int ii = dM; ii < num_pulse - dM ; ii++)
	{
	  avg_phase_diff(jj) = avg_phase_diff(jj) + phase_diff(ii,jj)/(num_pulse - 2*dM);
	}       
    }
  
}
*/

//--------------------------------------------------------------------------------
// Method Name: avgPulseCorrelaton
// Author: YZ
// Function Description:
// Averaging pulse correlation. The first and last 5 pulses are not used.
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::avgPulseCorrelation()
{
  unsigned int M, N;
  pulse_corr.size(M,N);
  avg_pulse_corr.resize(N);
  avg_pulse_corr = 0.0;

  unsigned int num_pulse = (unsigned int)(valid_data_per_burst/valid_data_per_pri);
  unsigned int dM = 0;

  //cout<<"num_pulse"<<num_pulse<<endl;

  for(unsigned int jj = 0; jj < N; jj++)
    {
      avg_pulse_corr(jj) = 0;

      for(unsigned int ii = dM; ii < num_pulse - dM ; ii++)
	  avg_pulse_corr(jj) = avg_pulse_corr(jj) + pulse_corr(ii,jj)/(num_pulse - 2*dM);   
    }

  avgTotalCorr = 0;
  for(unsigned int jj = 0; jj < N; jj++)
    avgTotalCorr = avgTotalCorr + avg_pulse_corr(jj)/N;
  
}
*/

//-----------------------------------------------------------
// Method Name: fractDopCentroid
// Author: YZ
// Function Description:
// Calculate fractional Doppler centroid with respect to prf
// Notes:
//-----------------------------------------------------------
/*
void RasList::fractDopCentroid()
{

  unsigned int N = avg_phase_diff.size();
  fDOPrePRF.resize(N);
  fDOPrePRF = 0.0;
  float pi = 4.*atan(1.);

  for(unsigned int jj = 0;jj < N; jj++)
    fDOPrePRF(jj) = avg_phase_diff(jj)/2/pi;   
 
}
*/


//---------------------------------------------
// Method Name: IntegDopCentroid
// Author: YZ
// Function Description:
// Get integer doppler centroid by calculation
// Notes: 
//---------------------------------------------
/*
void RasList::IntegDopCentroid(const Uvar& bore_dop)
{
  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
 
  Uvar Pri=ieb_.getPri();
  
  if(bore_dop > 0)
    int_m = (int)((bore_dop*Pri).getInUnits("") + 0.501);
  else if(bore_dop < 0)
    int_m = (int)((bore_dop*Pri).getInUnits("") - 0.499);
  else 
  int_m = 0;
 
}
*/

//-----------------------------------------------------
// Method Name: TotalDopCentroid
// Author: YZ
// Function Description:
// Calculate total Doppler centroid: fractal + integer
// Notes:
//-----------------------------------------------------
/*
void RasList::TotalDopCentroid()
{

  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
 
  Uvar Pri=ieb_.getPri();

  unsigned int N = fDOPrePRF.size();
  totalDop.resize(N);
  totalDop = 0.0;

  for(unsigned int jj = 0;jj < N; jj++)
    totalDop(jj) = (fDOPrePRF(jj)+int_m)/Pri;   

  totaldop_set_ = true;
  //plotTotalDop(); //Yan
}
*/
//---------------------------------------------------------
// Method Name: avgTotalDopCentroid
// Author: YZ
// Function Description:
// Weighted averaging for the total Doppler centroid along
// the fast time bin.
// Notes:
//---------------------------------------------------------
/*
void RasList::avgTotalDopCentroid()
{

  if(!totaldop_set_) ErrorMessage("RasList.cpp: total dopper centroid is not calculated yet.").throwMe();

  Uvar temp1 = 0.0;
  Uvar temp2 = 0.0;

  unsigned int N = totalDop.size();
  for(unsigned int jj = 0;jj < N; jj++)
    {
       temp1 = temp1 + totalDop(jj)*pulse_amp(jj);
       temp2 = temp2 + pulse_amp(jj);
    }   
    avgTotalDop = temp1/temp2; 

  Uvar Pri = ieb_.getPri();
  double pi = 4.*atan(1.);

  fractAvgTotalDop = atan2(imag(avgTotalCorr),real(avgTotalCorr))/2./pi;
  avgTotalDop = (fractAvgTotalDop + int_m)/Pri;

  cout<<"  Doppler centroid by data: "<<avgTotalDop<<endl;
}
*/
//--------------------------------------------------------------------------------
// Method Name: FractDopCentroid
// Author: YZ
// Function Description:
// Main code for fractional Doppler centroid
// Notes:
//--------------------------------------------------------------------------------
/*
void RasList::FractDopCentroid()
{

  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
 
  Uvar Pri=ieb_.getPri();

  get2Decho();
  avgPulseAmp();
  getComplex2Decho();
  pulseCorrelation();
  avgPulseCorrelation();
  get1DPhaseDiff();
  //getPhaseDiff();
  //averagePhaseDiff();
  fractDopCentroid();
  //plotDopCentroid();  //Yan
  
}
*/


//----------------------------
// plot pulse correlation
//----------------------------
/*
void RasList::plotPulseCorrelation()
{

  unsigned int M, N;
  pulse_corr.size(M,N);
  Uvec data_real("data_real",N);
  Uvec data_imag("data_imag",N);
  for(unsigned int ii = 0; ii < M; ii++)
    {
      for(unsigned int jj = 0;jj < N; jj++)
	{
	  data_real(jj) = real(pulse_corr(ii,jj));
	  data_imag(jj) = imag(pulse_corr(ii,jj));
	}

      Plot a;
      a.addY(data_real,"",line("solid","black",1),sym("none"));
      a.addY(data_imag,"",line("solid","red",1),sym("none"));
      a.setTitle("Pulse Correlation");
      a.show("x");    
    }
}
*/



//----------------------------
// plot phase difference
//----------------------------
/*
void RasList::plotPhaseDiff()
{

  unsigned int M, N;
  phase_diff.size(M,N);
  Uvec data_in_range("data_range",N);
  for(unsigned int ii = 0; ii < M; ii++)
    {
      for(unsigned int jj = 0;jj < N; jj++)
	data_in_range(jj) = phase_diff(ii,jj);

      Plot a;
      a.addY(data_in_range,"",line("solid","black",1),sym("none"));
      a.setTitle("Phase Difference");
      a.show("x");    
    }
}
*/

//--------------------------------
// plot averaged phase difference
//--------------------------------
/*
void RasList::plotAvgPhaseDiff()
{

  unsigned int N = avg_phase_diff.size();
  Uvec data_in_range("data_range",N);

  for(unsigned int jj = 0;jj < N; jj++)
    data_in_range(jj) = avg_phase_diff(jj);
  
  Plot a;
  a.addY(data_in_range,"",line("solid","black",1),sym("none"));
  a.setTitle("Averaged Phase Difference");
  a.show("x");    
 
}
*/

//--------------------------------
// plot averaged pulse amplitude
//--------------------------------
/*
void RasList::plotAvgPulseAmp()
{

  unsigned int N = pulse_amp.size();
  Uvec data_in_range("data_range",N);

  for(unsigned int jj = 0;jj < N; jj++)
    data_in_range(jj) = pulse_amp(jj);
  
  Plot a;
  a.addY(data_in_range,"",line("solid","black",1),sym("none"));
  a.setTitle("Averaged Pulse Amplitude");
  a.show("x");    
 
}
*/
//----------------------------------
// plot fractional Doppler centroid
//----------------------------------
/*
void RasList::plotDopCentroid()
{

  unsigned int N = fDOPrePRF.size();
  Uvec data_in_range("data_range",N);

  for(unsigned int jj = 0;jj < N; jj++)
    data_in_range(jj) = fDOPrePRF(jj);
  
  Plot a;
  a.addY(data_in_range,"",line("solid","black",1),sym("none"));
  a.setTitle("Fractional Doppler Centroid");
  a.show("x");    
 
}
*/
//----------------------------------
// plot total Doppler centroid
//----------------------------------
/*
void RasList::plotTotalDop()
{
  
  Plot a;
  a.addY(totalDop,"KHz",line("solid","black",1),sym("none"));
  a.setTitle("Total Doppler Centroid");
  a.show("x");    
 
}
*/
//------------------------------------------
// compare with boresight Doppler frequency
//------------------------------------------
/*
void RasList::compareBoreSightDop(const Uvar& bore_dop)
{
  
  Uvec Bore_dop("Bore_dop",totalDop.size());

  Bore_dop = bore_dop;

  Plot a;
  a.addY(totalDop,"KHz",line("solid","black",1),sym("none"));
  a.addY(Bore_dop,"KHz",line("solid","red",1),sym("none"));
  a.setTitle("Total Doppler Centroid");
  a.setXlabelandUnits("Fast Time Index");
  a.setYlabelandUnits("Doppler Centroid");
  a.show("x");    
 
}
*/
