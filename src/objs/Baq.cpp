//------------------------------------------------------------------
// Baq.cpp
//
// This file contains the Baq methods.
// The Baq class provides following processes:
//   (1) load parameters from Ieb:                Baq::setIeb()
//   (2) Digitizing anologue SAR data:            Baq::ADC()
//   (3) Calculate threshold:                     Baq::setThreshold()
//   (4) Encoding digitized SAR:                  Baq::Encode_Nbit()
//   (5) Decoding SAR data:                       Baq::Decode_Nbit()
//   (6) Encoding threshold:                      Baq::Encode_thresh()
//
//------------------------------------------------------------------


#include <stdlib.h>
#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Baq.h"
#include "Io.h"
#include "Ieb.h"
#include "Sab.h"
#include "Utils.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//-------------------------
// Baq Methods
//-------------------------
//-------------------------
// Constructors
//-------------------------
Baq::Baq()
{
  N_=0;
  Np_=0; Nt_=0; Nb_=0; Ns_=0; Ns_p_=0; Ne_=0;
  dbg=0;
}

//------------------------------------------------------------------------
// Debug
//------------------------------------------------------------------------
void Baq::debug_mode()
{
  dbg=1;
}


//------------------------------------------------------------------------
// setIeb(Ieb& ieb)
// 
// This method sets parameters for Baq class.
//
// Parameters set for Baq:
// Np_:    number of points per echo (pulse)
// Nt_:    number of points per burst
// Ne_:    number of echoes per burst (same as pul, see below)
// Nb_:    number of blocks per echo (fixed to 24)
// Ns_:    number of points per block, except the last block
// Ns_p_:  number of points in the last block
// N_mod_: Baq mode
//
// Parameters from Ieb useful for Baq
// SR:    sampling rate
// pri:   pulse repetition inteval
// pul:   number of pulses per burst
//
// Baq modes:
// N_mod_ = 0: 8-to-2 bit BAQ
// N_mod_ = 1: 8-to-1 bit BAQ
// N_mod_ = 2: 8-to-0 bit BAQ
// N_mod_ = 3: 8 bit to 2 MSB's
// N_mod_ = 4: 8 bit to 4 MSB's
// N_mod_ = 5: 8 bits straight
// N_mod_ = 6: 8-to-4 bit BAQ (Low Res ALT mode)
// N_mod_ = 7: 8-to-4 bit BAQ (Hi Res ALT mode)
//------------------------------------------------------------------------
void Baq::setParams(const Uvar& pri, const unsigned int& pul,
		    const unsigned int& baq_mode, const Uvar& SR, 
		    const int& tro_int)
{
  Uvar prf = 1/pri;

  N_mod_=baq_mode;


  int N_mod_temp=N_mod_;

  if(dbg==1)
    std::cout<<"--------------- set Ieb Parameters ---------------"<<std::endl;
  //---------------------------------------------------
  // for debugging purpose, set N_mod_ to a fixed value
  //---------------------------------------------------
  //N_mod_=5; 
  if(N_mod_temp!=N_mod_) 
    std::cout << "   Warning!!! setIeb is in debugging mode."<<std::endl;

  if(dbg==1)
  {
    //std::cout << "   prf: " << prf << std::endl;
    //std::cout << "   tau_p: " << tau_p << std::endl;
    //std::cout << "   dutycycle: " << dutycycle << std::endl;
    //std::cout << "   bpd: " << bpd << std::endl;
    //std::cout << "   chirp_start_frequency: " << chirp_start_frequency << std::endl;
        
    //std::cout << "   BR: " << BR << std::endl;
    std::cout << "   pulse repetition frequency (prf): " << prf << std::endl;
    std::cout << "   number of echoes per burst (pul): " << pul << std::endl;
    std::cout << "   Baq mode (N_mod_): " << N_mod_ << std::endl;
    std::cout << "   sampling rate (SR): " << SR << std::endl;
  }
  //Np_=int((SR*pri).getInUnits(""));
  // Np_=int((SR/prf).getInUnits("")+0.5);
  Np_ = (unsigned int) round_double( (SR*pri).getInUnits(""));
//Nt_=Np_*ieb.getPul(); //a direct way to get Nt_

  // compute number of samples including effect of tro
  Nt_=Np_*pul+Np_*tro_int;

  Nb_=24;
  Ns_=int(Np_/Nb_);
  Ns_p_=Np_-(Nb_-1)*Ns_;
  //Ne_=Nt_/Np_; //originally defined in Blue Book
  int receive_length = int(pul)+tro_int;
  if(receive_length<0) ErrorMessage("pul + tro <0").throwMe();
  Ne_=(unsigned int) receive_length; //original in Blue Book is Ne_Nt_/Np_

  if(dbg==1)
  {
    std::cout << "   points per echo (Np_): " << Np_ << std::endl;	
    std::cout << "   points per burst (Nt_): " << Nt_ << std::endl;
    std::cout << "   echoes per burst (Ne_): " << Ne_ << std::endl;
    std::cout << "   blocks per burst (Nb_, fixed to 24): " << Nb_ << std::endl;
    std::cout << "   points per block, except the last block (Ns_): " << Ns_ << std::endl;
    std::cout << "   points in the last block (Ns_p_): " << Ns_p_ << std::endl;
  }
  //----------------------
  // validation 
  //----------------------
  if(dbg==1)
  {
    if(Nt_==Np_*Ne_)
      std::cout<<"   check Np_, Nt_, and Ne_: passed"<<std::endl;
    else ErrorMessage("Nt is not equal to Np * Ne !").throwMe();
    if(Np_==Ns_*(Nb_-1)+Ns_p_)
      std::cout<<"   check Nb_, Ns_, and Ns_p_: passed"<<std::endl;
    else ErrorMessage("Np is not equal to Ns * (Nb-1) + Ns_p !").throwMe();
  } 
}


//------------------------------------------------------------------------
// setIeb(Ieb& ieb)
// 
// This method sets parameters for Baq class.
//
// Parameters set for Baq:
// Np_:    number of points per echo (pulse)
// Nt_:    number of points per burst
// Ne_:    number of echoes per burst (same as pul, see below)
// Nb_:    number of blocks per echo (fixed to 24)
// Ns_:    number of points per block, except the last block
// Ns_p_:  number of points in the last block
// N_mod_: Baq mode
//
// Parameters from Ieb useful for Baq
// SR:    sampling rate
// pri:   pulse repetition inteval
// pul:   number of pulses per burst
//
// Baq modes:
// N_mod_ = 0: 8-to-2 bit BAQ
// N_mod_ = 1: 8-to-1 bit BAQ
// N_mod_ = 2: 8-to-0 bit BAQ
// N_mod_ = 3: 8 bit to 2 MSB's
// N_mod_ = 4: 8 bit to 4 MSB's
// N_mod_ = 5: 8 bits straight
// N_mod_ = 6: 8-to-4 bit BAQ (Low Res ALT mode)
// N_mod_ = 7: 8-to-4 bit BAQ (Hi Res ALT mode)
//------------------------------------------------------------------------
void Baq::setIeb(Ieb& ieb)
{

  Uvar pri = ieb.getPri();
  //Uvar tau_p = ieb.getTaup();
  //Uvar dutycycle = tau_p/pri;
  //Uvar bpd =ieb.getBpd();
  //Uvar chirp_start_frequency=ieb.getChirpStartFrequency();
  //Uvar BR=ieb.getBr();
  unsigned int pul = ieb.getPul();
  N_mod_=ieb.getBaq();
  //std::cout<<"BAQ mode: "<<N_mod_<<std::endl;
  Uvar SR=ieb.getAdc();
  int tro_int=ieb.getTroInPriUnits();
  unsigned int baq_mode=ieb.getBaq();
  setParams(pri,pul,baq_mode,SR,tro_int);
}

//------------------------------------------------------------------------
// Baq.loadEncodedata(SUvec)
//
// This method reads encoded SAR data (unsigned short) from a file into 
// memory. This method is not used in the Baq class. 
//------------------------------------------------------------------------
void Baq::LoadEncodedata(SUvec& words)
{
   FileMgr encoded_data("bit_file.dat","rb");
   int i=0;

   std::cout<<"--------------- load data in words ---------------"<<std::endl;
   
   if(N_mod_==0) MN_words_=Nt_/8;       //8-2 bit (normally used for SAR mode)
   else if(N_mod_==1) MN_words_=Nt_/16; //8-1 bit (need threshold !!!)
   else if(N_mod_==2) MN_words_=0;      //8-0 bit
   else if(N_mod_==3) MN_words_=Nt_/8;  //MSB 8-2 bit
   else if(N_mod_==4) MN_words_=Nt_/4;  //MSB 8-4 bit
   else if(N_mod_==5) MN_words_=Nt_/2;  //8-8 bit straight
   else if((N_mod_==6)||(N_mod_==7)) MN_words_=Nt_/4;  //8-4 bit (Lo or Hi-Re ALT mode)
   else ErrorMessage("Bit number is wrong.").throwMe(); 

   if(MN_words_!=0)
   {
     words.resize(MN_words_);

     while(!encoded_data.eof())
     {
       encoded_data.read(words(i));
       i++;	   
     }
   }
   else
     std::cout<<"   8-to-0 bit BAQ, no buffer created." <<std::endl;

   std::cout << "   Total number of words (16-bits): "<<MN_words_<<std::endl;
   //std::cout << "Total number of words (16-bits) verified: "<<i<<std::endl;
   //if(i==MN_words_) 
   //  std::cout<<"   check number of words: passed"<<std::endl;
   //else 
   //  ErrorMessage("Number of words is wrong!").throwMe();

}

//------------------------------------------------------------------------
// addHeaderFooter

void Baq::addHeaderFooter(const Charvec& words, const Ivec& Thresh)
{

   FileMgr baq_dat("yan_msb8t4_7.10Dec97_1947","rb");
   FileMgr baq_out("yan_msb8t4_7.10Dec97_1947_encoded","wb");
    
   unsigned short header;
   unsigned char baq_thresh;

   //-------------------------------------
   // Determine the size of encoded data
   //-------------------------------------
   if(N_mod_==0) MN_words_=Nt_/4;       //8-2 bit (normally used for SAR mode)
   else if(N_mod_==1) MN_words_=Nt_/8; //8-1 bit (need threshold !!!)
   else if(N_mod_==2) MN_words_=0;      //8-0 bit
   else if(N_mod_==3) MN_words_=Nt_/4;  //MSB 8-2 bit
   else if(N_mod_==4) MN_words_=Nt_/2;  //MSB 8-4 bit
   else if(N_mod_==5) MN_words_=Nt_;  //8-8 bit straight
   else if((N_mod_==6)||(N_mod_==7)) MN_words_=Nt_/2;  //8-4 bit (Lo or Hi-Re ALT mode)
   else ErrorMessage("Bit number is wrong.").throwMe(); 

   //---------------------
   // write header
   //---------------------
   for(int i=0;i<90;i++)
   {
     baq_dat.read(header);
     baq_out.write(header);
   }

   for(int i=0;i<MN_words_;i++)
     baq_out.write(words(i));

   //--------------------------------------------
   // write footer (threshold): 12 (16-bit) words
   //--------------------------------------------
   for(int i=0;i<Nb_;i++)
   {
     baq_thresh=(unsigned char)(Thresh(i));
     std::cout<<"  Thresh: "<<int(baq_thresh)<<std::endl; 
     baq_out.write(baq_thresh);
   }     
   //------------------------------------------------
   // filling blanks for the rest 10 (16-bit) words
   //------------------------------------------------
   for(int i=0;i<10;i++)
   {
     header=0x05ff;
     baq_out.write(header);
   }

   std::cout<<" MN_words_: "<< MN_words_<<std::endl;
  
}

float Baq::minVariance(const Ivec& th){
  // initialize up in/out arrays
  if(Nt_==0) return(0.0);
  Charvec z("zeros",Nt_);
  Charvec w("words_name");
  Fvec d("data",Nt_);

  // set up least threshold bits
  for(int c=0;c<Nt_;c++) z(c)=(char)0;
  Encode_Nbit(z,th,w);
  Decode_Nbit(w,th,d);
  float var=0;
  for(int c=0;c<Nt_;c++){
    var+=d(c)*d(c);
  }
  return(var/Nt_);
    
}
//------------------------------------------------------------------------
// Baq.Decode_Nbit(Charvec, Ivec, SUvec, Fvec)
//
// This method decodes SAR data stored in words (char array) based on the 
// threshold Th8_N, and writes the decoded SAR data into echo_N. dbg is a
// flag for debug (1-print out notes, 0-without printing out notes).
//------------------------------------------------------------------------
void Baq::Decode_Nbit(const Charvec& words, const Ivec& Th8_N, Fvec& echo_N)
{
  //FileMgr encoded_data("bit_file.dat","rb");
   unsigned short bits;
   unsigned short sign;
   int N_echo=0;

   Fvec d2("d2_name",2);
   Fvec d4("d4_name",4);
   Fvec dMSB4("dMSB4_name",4);
   Fvec d16("d16_name",16);
   Fvec dMSB16("d16_name",16);

   if(Nt_!=0) echo_N.resize(Nt_);

   if(dbg==1)
   {
     if(Nt_==int(echo_N.size())) 
       std::cout<<"   check number of data: passed"<<std::endl;
     else 
       ErrorMessage("Number of data is wrong!").throwMe();
   }

   d2(1)=-64; //1                      // 8-to-1 BAQ
   d2(0)=64;  //0

   d4(3)=-1.64; //11                   // 8-to-2 BAQ
   d4(2)=-0.49; //10
   d4(1)= 1.64; //01
   d4(0)= 0.49; //00

   dMSB4(3)=-96;  //11                 // 8-to-2 MSB
   dMSB4(2)=-32;  //10
   dMSB4(1)= 96;  //01
   dMSB4(0)= 32;  //00

   d16(15)=-2.498; //1111              // 8-to-4 BAQ
   d16(14)=-1.891; //1110
   d16(13)=-1.479; //1101 
   d16(12)=-1.148; //1100
   d16(11)=-0.861; //1011
   d16(10)=-0.600; //1010
   d16(9)= -0.355; //1001
   d16(8)= -0.117; //1000
   d16(7)=  2.498; //0111
   d16(6)=  1.891; //0110
   d16(5)=  1.479; //0101
   d16(4)=  1.148; //0100
   d16(3)=  0.861; //0011
   d16(2)=  0.600; //0010
   d16(1)=  0.355; //0001
   d16(0)=  0.117; //0000

   dMSB16(15)=-120; //1111             // 8-to-4 MSB
   dMSB16(14)=-104; //1110
   dMSB16(13)=-88; //1101 
   dMSB16(12)=-72; //1100
   dMSB16(11)=-56; //1011
   dMSB16(10)=-40; //1010
   dMSB16(9)= -24; //1001
   dMSB16(8)= -8; //1000
   dMSB16(7)=  120; //0111
   dMSB16(6)=  104; //0110
   dMSB16(5)=  88; //0101
   dMSB16(4)=  72; //0100
   dMSB16(3)=  56; //0011
   dMSB16(2)=  40; //0010
   dMSB16(1)=  24; //0001
   dMSB16(0)=  8; //0000
     
   //-------------------------------------
   // Determine the size of encoded data
   //-------------------------------------
   if(N_mod_==0) MN_words_=Nt_/4;       //8-2 bit (normally used for SAR mode)
   else if(N_mod_==1) MN_words_=Nt_/8; //8-1 bit (need threshold !!!)
   else if(N_mod_==2) MN_words_=0;      //8-0 bit
   else if(N_mod_==3) MN_words_=Nt_/4;  //MSB 8-2 bit
   else if(N_mod_==4) MN_words_=Nt_/2;  //MSB 8-4 bit
   else if(N_mod_==5) MN_words_=Nt_;  //8-8 bit straight
   else if((N_mod_==6)||(N_mod_==7)) MN_words_=Nt_/2;  //8-4 bit (Lo or Hi-Re ALT mode)
   else ErrorMessage("Bit number is wrong.").throwMe(); 

   //-------------------------------------
   // decode main routine
   //-------------------------------------
   if(N_mod_==0) 
   {
     if(dbg==1)
       std::cout<<"--------------- Decode 8-2 bit BAQ ---------------"<<std::endl;
     for(int w=0; w<MN_words_; w++)
     {	   
       int ie=N_echo/Np_;
       int ib=(N_echo-ie*Np_)/Ns_;
       if(ib>Nb_-1) ib=Nb_-1;

       if((N_echo<4)&&(dbg==1)) print_bits_endl(words(w),8);

       for(int i=0;i<4;i++)
       {
	 if(int(fmod(double(w),2))==1)
           bits=bitget((unsigned char)(words(w-1)),2*i,2*i+1);
         else
           bits=bitget((unsigned char)(words(w+1)),2*i,2*i+1);

         
         echo_N(N_echo)=d4(bits)*Th8_N(ib)/2;
	 if((N_echo<4)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }
     }
   }
   else if(N_mod_==1) 
   {
     if(dbg==1)
     {
       std::cout<<"--------------- Decode 8-1 bit BAQ ---------------"<<std::endl;
       std::cout<<"                (missing threshold)               "<<std::endl;
     }
     for(int w=0; w<MN_words_; w++)
     {	   
       int ie=N_echo/Np_;
       int ib=(N_echo-ie*Np_)/Ns_;
       if(ib>Nb_-1) ib=Nb_-1;

       if((N_echo<8)&&(dbg==1))  print_bits_endl(words(w),8);
          
       for(int i=0;i<8;i++)
       {
         if(int(fmod(double(w),2))==1)
           bits=bitget((unsigned char)(words(w-1)),i,i);
         else
	   bits=bitget((unsigned char)(words(w+1)),i,i);

         echo_N(N_echo)=d2(bits); //threshold??
	 if((N_echo<8)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }
     }
   }
   else if(N_mod_==2) 
   {
     if(dbg==1)
     {
       std::cout<<"--------------- Decode 8-0 bit BAQ ---------------"<<std::endl;
       std::cout<<"   Empty data buffer"<<std::endl;
     }
   }
   else if(N_mod_==3) 
   {
     if(dbg==1)
       std::cout<<"--------------- Decode 8-2 bit MSB ---------------"<<std::endl;
     for(int w=0; w<MN_words_; w++)
     {	  
       if((N_echo<4)&&(dbg==1)) print_bits_endl(words(w),8);
          
       for(int i=0;i<4;i++)
       {
         if(int(fmod(double(w),2))==1)
           bits=bitget((unsigned char)(words(w-1)),2*i,2*i+1);
         else
           bits=bitget((unsigned char)(words(w+1)),2*i,2*i+1);

         echo_N(N_echo)=dMSB4(bits);
	 if((N_echo<4)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }
     }
   }
   else if(N_mod_==4) 
   {
     if(dbg==1)
       std::cout<<"--------------- Decode 8-4 bit MSB ---------------"<<std::endl;
     for(int w=0; w<MN_words_; w++)
     {	  
       if((N_echo<2)&&(dbg==1)) print_bits_endl(words(w),8);
          
       for(int i=0;i<2;i++)
       {
         if(int(fmod(double(w),2))==1)
           bits=bitget((unsigned char)(words(w-1)),4*i,4*i+3);
         else
           bits=bitget((unsigned char)(words(w+1)),4*i,4*i+3);
           
         echo_N(N_echo)=dMSB16(bits);
	 if((N_echo<2)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }
     }
   }
   else if(N_mod_==5) 
   {
     if(dbg==1)
       std::cout<<"--------------- Decode 8-8 bit BAQ ---------------"<<std::endl;
     for(int w=0; w<MN_words_; w++)
     {	  
       if((N_echo<1)&&(dbg==1))  print_bits_endl(words(w),8);
          
       for(int i=0;i<1;i++)
       {
         if(int(fmod(double(w),2))==1) 
	 {        
           bits=bitget((unsigned char)(words(w-1)),8*i,8*i+6);
           sign=bitget((unsigned char)(words(w-1)),8*i+7,8*i+7);
	 }
         else
	 {        
           bits=bitget((unsigned char)(words(w+1)),8*i,8*i+6);
           sign=bitget((unsigned char)(words(w+1)),8*i+7,8*i+7);
	 }

         if(sign==0) 
            echo_N(N_echo)=bits;
         else
	    echo_N(N_echo)=-bits;

	 if((N_echo<1)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }
     }
   }
   else if((N_mod_==6)||(N_mod_==7)) 
   {
     if(dbg==1)
       std::cout<<"--------------- Decode 8-4 bit BAQ ---------------"<<std::endl;
     for(int w=0; w<MN_words_; w++)
     { 
       int ie=N_echo/Np_;
       int ib=(N_echo-ie*Np_)/Ns_;
       if(ib>Nb_-1) ib=Nb_-1;

       if((N_echo<2)&&(dbg==1)) print_bits_endl(words(w),8);

       for(int i=0;i<2;i++)
       {
         if(int(fmod(double(w),2))==1) 
           bits=bitget((unsigned char)(words(w-1)),4*i,4*i+3);
         else
           bits=bitget((unsigned char)(words(w+1)),4*i,4*i+3);

         echo_N(N_echo)=d16(bits)*Th8_N(ib)/2;
	 if((N_echo<2)&&(dbg==1)) std::cout<<"  ecoded data: "<<echo_N(N_echo)<<std::endl;
         N_echo++;
       }	   
     }
   }
   else
     ErrorMessage("Bit number is wrong.").throwMe(); 
}

//---------------------------------------------------------------------------
// Encode threshold from Ivec to Array1D<unsigned short>
//---------------------------------------------------------------------------
void Baq::Encode_thresh(const Ivec& Thresh, Array1D<unsigned short>& thresh)
{
  //int N_thresh = Thresh.size();
  //cout<<"N_thresh: "<<N_thresh<<endl;

  for(unsigned int i = 0; i < 12; i++)
    {
      
      //byte swapped
      bitset(thresh(i),8,15,Thresh(2*i));
      //cout<<"  Thresh(i): "<<Thresh(2*i)<< "thresh(i): "<<thresh(i)<<endl;
      bitset(thresh(i),0,7,Thresh(2*i+1));
      //cout<<" Thresh(i): "<<Thresh(2*i+1)<< "thresh(i): "<<thresh(i)<<endl;

    }
}

//---------------------------------------------------------------------------
// Baq.Encode_Nbit(const Charvec, const Ivec, Cdata)
//
// This method encodes SAR data stored in "echo_8" (char array) based on the 
// threshod "Thresh" and writes the encoded data into "words". dbg is a flag 
// for debug (1-print out notes, 0-without printing out notes).
// 
//---------------------------------------------------------------------------
void Baq::Encode_Nbit(const Charvec& echo_8, const Ivec& Thresh, Charvec& words)
{

  int NP_echo_8 = echo_8.size();
  if(NP_echo_8!=Nt_)
    ErrorMessage("The size of SAR data file is not correct.").throwMe();

  //FileMgr encoded_data("bit_file.dat","wb");

  int count=0;  
  int W_count=0;
  char bits;
  unsigned char sign;
  unsigned char number;
  //unsigned short p;


  //-------------------------------------
  // Determine the size of encoded data
  //-------------------------------------

  
  if(N_mod_==0){
    MN_words_=Nt_/4;       //8-2 bit (normally used for SAR mode)
    if(Nt_%4!=0) MN_words_++;
  }
  else if(N_mod_==1){
    MN_words_=Nt_/8; //8-1 bit (need threshold !!!)
    if(Nt_%8!=0) MN_words_++;
  }
  else if(N_mod_==2){
    MN_words_=0;      //8-0 bit
  }
  else if(N_mod_==3){
    MN_words_=Nt_/4;  //MSB 8-2 bit
    if(Nt_%4!=0) MN_words_++;    
  }
  else if(N_mod_==4){
    MN_words_=Nt_/2;  //MSB 8-4 bit
    if(Nt_%2!=0) MN_words_++;
  }
  else if(N_mod_==5){
    MN_words_=Nt_;  //8-8 bit straight
  }
  else if((N_mod_==6)||(N_mod_==7)){
    MN_words_=Nt_/2;  //8-4 bit (Lo or Hi-Re ALT mode)
    if(Nt_%2!=0) MN_words_++;  
  }

  else ErrorMessage("Bit number is wrong.").throwMe(); 

  // add extra byte if MN_words is odd.
  if(MN_words_%2==1) MN_words_++;

  if(MN_words_!=0) words.resize(MN_words_);

  //-------------------------------
  // 2 bit encoding
  //-------------------------------
  if(N_mod_==0)
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-2 bit BAQ ----------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {
      int ie=ii/Np_;
      int ib=(ii-ie*Np_)/Ns_;
      if(ib>Nb_-1) ib=Nb_-1;    // get number of blocks
      
      if(echo_8(ii)<-Thresh(ib)/2) 
	  bits=3; //11
      else if((echo_8(ii)>=-Thresh(ib)/2)&&(echo_8(ii)<=0))
          bits=2; //10
      else if((echo_8(ii)>0)&&(echo_8(ii)<=Thresh(ib)/2))
          bits=0; //00
      else if(echo_8(ii)>Thresh(ib)/2)
	  bits=1; //01
      else
            ErrorMessage("Encoding Error!").throwMe();

      //--------------------
      // packing 2 bit data
      //--------------------
      //if(count==8) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<2*count);
      //if((ii<8)&&(dbg==1))
      //{
      //  std::cout<<"data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   ";  
      //  print_bits_endl(p,16);
      //}     
      //if(count==7)
      //   encoded_data.write(p);
      //count++;


      if(count==4) 
      {
        count=0; 
        W_count++;
      }

      if(int(fmod(double(W_count),2))==1)
        words(W_count-1)=words(W_count-1)|(bits<<2*count);
      else
	words(W_count+1)=words(W_count+1)|(bits<<2*count);
      

      if((ii<4)&&(dbg==1))
      {
        std::cout<<"data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   ";  
        print_bits_endl(words(W_count),8);
      }     
      count++;
      
    }

    //std::cout<<"number of words: "<<W_count<<std::endl;
  } 
  //-------------------------------
  // 1-bit encoding
  //-------------------------------
  else if(N_mod_==1)
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-1 bit BAQ ---------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {
      int ie=ii/Np_;
      int ib=(ii-ie*Np_)/Ns_;
      if(ib>Nb_-1) ib=Nb_-1;    // get number of blocks
      
      if(echo_8(ii)<0) 
	  bits=1; //1
      else if(echo_8(ii)>=0)
	  bits=0; //0
      else
            ErrorMessage("Encoding Error!").throwMe();

      //--------------------
      // packing 1-bit data
      //--------------------
      //if(count==16) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<count);
      //if((ii<16)&&(dbg==1))
      //{
      //  std::cout<<"   data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   "; 
      //  print_bits_endl(p,16);
      //}     
      //if(count==15)
      //   encoded_data.write(p);
      //count++;

      if(count==8) 
      {
        count=0; 
        W_count++;
      }

      if(int(fmod(double(W_count),2))==1)
        words(W_count-1)=words(W_count-1)|(bits<<count);
      else
        words(W_count+1)=words(W_count+1)|(bits<<count);

      if((ii<8)&&(dbg==1))
      {
        std::cout<<"data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   ";  
        print_bits_endl(words(W_count),8);
      }     
      count++;

    }
  } 
  //-------------------------------
  // 0-bit BAQ encoding
  //-------------------------------
  else if(N_mod_==2)
  {
    if(dbg==1)
    {
      std::cout <<"--------------- Encode 8-0 bit BAQ ---------------"<<std::endl;
      std::cout<<"   Empty data buffer"<<std::endl;
    }
  }
  //-------------------------------
  // 2-bit MSB encoding
  //-------------------------------
  else if(N_mod_==3)
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-2 bit MSB ---------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {
      
      sign=bitget(unsigned(echo_8(ii)),7,7);
      number=bitget(echo_8(ii),6,6);

      if(sign==0) bits=number;
      if(sign==1) bits=2+(1-number); 

      //--------------------
      // packing 2 bit MSB data
      //--------------------
      //if(count==8) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<2*count);
      //if((ii<8)&&(dbg==1))
      //  {
      //    std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
      //    for(int t=128; t>0; t=t/2)
      //      if(echo_8(ii) & t) std::cout << "1";
      //      else std::cout << "0";
      //    std::cout<<"   code: "<<int(bits)<<"   ";   
      //    print_bits_endl(p,16);
      //  }     
      //if(count==7)
      //   encoded_data.write(p);
      //count++;

      if(count==4) 
      {
        count=0; 
        W_count++;
      }
      if(int(fmod(double(W_count),2))==1)
        words(W_count-1)=words(W_count-1)|(bits<<2*count);
      else
        words(W_count+1)=words(W_count+1)|(bits<<2*count);

      if((ii<4)&&(dbg==1))
        {
          std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
          for(int t=128; t>0; t=t/2)
            if(echo_8(ii) & t) std::cout << "1";
            else std::cout << "0";
          std::cout<<"   code: "<<int(bits)<<"   ";   
          print_bits_endl(words(W_count),8);
        }     
      count++;

    }
  } 
  //-------------------------------
  // 4-bit MSB encoding
  //-------------------------------
  else if(N_mod_==4)
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-4 bit MSB ---------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {
      
      sign=bitget(unsigned(echo_8(ii)),7,7);
      number=bitget(echo_8(ii),4,6);

      if(sign==0) bits=number;
      if(sign==1) bits=8+(7-number); 

      //------------------------
      // packing 4 bit MSB data
      //------------------------
      //if(count==4) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<4*count);
      //if((ii<4)&&(dbg==1))
      //{
      //  std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
      //  for(int t=128; t>0; t=t/2)
      //    if(echo_8(ii) & t) std::cout << "1";
      //  else std::cout << "0";
      //  std::cout<<"   code: "<<int(bits)<<"   ";   
      //  print_bits_endl(p,16);
      //}     
      //if(count==3)
      //   encoded_data.write(p);
      //count++;

      if(count==2) 
      {
        count=0; 
        W_count++;
      }
      if(int(fmod(double(W_count),2))==1)
        words(W_count-1)=words(W_count-1)|(bits<<4*count);
      else
        words(W_count+1)=words(W_count+1)|(bits<<4*count);

      if((ii<2)&&(dbg==1))
      {
        std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
        for(int t=128; t>0; t=t/2)
          if(echo_8(ii) & t) std::cout << "1";
        else std::cout << "0";
        std::cout<<"   code: "<<int(bits)<<"   ";   
        print_bits_endl(words(W_count),8);
      }     
      count++;

    }
  }
  else if(N_mod_==5)
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-8 bit BAQ ---------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {      
      sign=bitget(unsigned(echo_8(ii)),7,7);
      number=bitget(echo_8(ii),0,6);

      if(sign==0) bits=number;
      else bits=128-number;

      //--------------------
      // packing 8 bit data
      //--------------------
      //if(count==2) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<8*count);
      //p=p|(sign<<8*count+7);
      //if((ii<2)&&(dbg==1))
      //{
      //  std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
      //  for(int t=128; t>0; t=t/2)
      //    if(echo_8(ii) & t) std::cout << "1";
      //  else std::cout << "0"; 
      //  std::cout<<"   code: "<<int(bits)<<"   ";   
      //  print_bits_endl(p,16);
      //}     
      //if(count==1)
      //   encoded_data.write(p);
      //count++;

      if(count==1) 
      {
        count=0; 
        W_count++;
      }
      if(int(fmod(double(W_count),2))==1)
      {
        words(W_count-1)=words(W_count-1)|(bits<<8*count);
        words(W_count-1)=words(W_count-1)|(sign<<8*count+7);
      }
      else
      {
        words(W_count+1)=words(W_count+1)|(bits<<8*count);
        words(W_count+1)=words(W_count+1)|(sign<<8*count+7);
      }

      if((ii<1)&&(dbg==1))
      {
        std::cout<<"data:"<<int(echo_8(ii))<<"   two's:";
        for(int t=128; t>0; t=t/2)
          if(echo_8(ii) & t) std::cout << "1";
        else std::cout << "0"; 
        std::cout<<"   code: "<<int(bits)<<"   ";   
        print_bits_endl(words(W_count),8);
      }     
      count++;

    } //for
  }  
  else if((N_mod_==6)||(N_mod_==7))
  {
    if(dbg==1)
      std::cout <<"--------------- Encode 8-4 bit BAQ ---------------"<<std::endl;
    for(int ii=0;ii<NP_echo_8;ii++)
    {      
      int ie=ii/Np_;
      int ib=(ii-ie*Np_)/Ns_;
      if(ib>Nb_-1) ib=Nb_-1;    // get number of blocks

      Fvec z("z_name",7);
      z(0)=0.235*Thresh(ib);
      z(1)=0.475*Thresh(ib);
      z(2)=0.730*Thresh(ib);
      z(3)=1.000*Thresh(ib);
      z(4)=1.310*Thresh(ib);
      z(5)=1.680*Thresh(ib);
      z(6)=2.200*Thresh(ib);

      
      if(echo_8(ii)<-z(6)/2) 
	  bits=15; // 1111
      else if((echo_8(ii)>=-z(6)/2)&&(echo_8(ii)<-z(5)/2))
          bits=14; // 1110
      else if((echo_8(ii)>=-z(5)/2)&&(echo_8(ii)<-z(4)/2))
          bits=13; // 1101
      else if((echo_8(ii)>=-z(4)/2)&&(echo_8(ii)<-z(3)/2))
	  bits=12; // 1100
      else if((echo_8(ii)>=-z(3)/2)&&(echo_8(ii)<-z(2)/2))
	  bits=11; // 1011
      else if((echo_8(ii)>=-z(2)/2)&&(echo_8(ii)<-z(1)/2))
	  bits=10; //1010
      else if((echo_8(ii)>=-z(1)/2)&&(echo_8(ii)<-z(0)/2))
	  bits=9; //1001
      else if((echo_8(ii)>=-z(0)/2)&&(echo_8(ii)<0))
	  bits=8; //1000
      else if((echo_8(ii)>=0)&&(echo_8(ii)<z(0)/2))
	  bits=0; //0000
      else if((echo_8(ii)>=z(0)/2)&&(echo_8(ii)<z(1)/2))
	  bits=1; //0001
      else if((echo_8(ii)>=z(1)/2)&&(echo_8(ii)<z(2)/2))
	  bits=2; //0010
      else if((echo_8(ii)>=z(2)/2)&&(echo_8(ii)<z(3)/2))
	  bits=3; //0011
      else if((echo_8(ii)>=z(3)/2)&&(echo_8(ii)<z(4)/2))
	  bits=4; //0100
      else if((echo_8(ii)>=z(4)/2)&&(echo_8(ii)<z(5)/2))
	  bits=5; //0101
      else if((echo_8(ii)>=z(5)/2)&&(echo_8(ii)<z(6)/2))
	  bits=6; //0110
      else if(echo_8(ii)>=z(6)/2)
	  bits=7; //0111
      else 
          ErrorMessage("Encoding Error!").throwMe();

      //--------------------
      // packing 4 bit data
      //--------------------
      //if(count==4) count=0;
      //if(count==0) p=0;
      //p=p|(bits<<4*count);
      //if((ii<4)&&(dbg==1))
      //{
      //  std::cout<<"data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   ";
      //  print_bits_endl(p,16);
      //}     
      //if(count==3)
      //   encoded_data.write(p);
      //count++;

      if(count==2) 
      {
        count=0; 
        W_count++;
      }
      if(int(fmod(double(W_count),2))==1)
        words(W_count-1)=words(W_count-1)|(bits<<4*count);
      else
        words(W_count+1)=words(W_count+1)|(bits<<4*count);

      if((ii<2)&&(dbg==1))
      {
        std::cout<<"data:"<<int(echo_8(ii))<<"   code: "<<int(bits)<<"   ";
        print_bits_endl(words(W_count),8);
      }     
      count++;

    } //for
  } 
  else
    ErrorMessage("There is no such bit number for encoding.").throwMe();
}

//---------------------------------------------------------------------------
// Baq.Encode_Nbit(const Fvec, const Ivec, Charvec)
//
// This method is overloaded with the input SAR data being float array instead
// of char array. It encodes SAR data stored in "echo" (float array) based on 
// the threshod "Thresh" and writes the encoded data into "words". dbg is a 
// flag for debug (1-print out notes, 0-without printing out notes).
// 
//---------------------------------------------------------------------------
void Baq::Encode_Nbit(const Fvec& echo, const Ivec& Thresh, Charvec& words)
{
  
  int NP_echo_8 = echo.size();
  Charvec echo_8("echo_8_name",NP_echo_8);
  echo_8=Baq::ADC(echo);
  if(NP_echo_8!=Nt_)
    ErrorMessage("The size of SAR data file is not correct.").throwMe();
    
  Encode_Nbit(echo_8,Thresh,words);
}

//---------------------------------
// print_Bits(unsigned short, int)
//---------------------------------
void Baq::print_bits_endl(const unsigned short p, const int a)
{
  int t_=1;
  for(int i=1;i<a;i++)
    t_=t_*2;
  for(int t=t_; t>0; t=t/2)
    if(p & t) std::cout << "1";
    else std::cout << "0";
  std::cout << std::endl;
}

//---------------------------------
// print_Bits(unsigned short, int)
//---------------------------------
void Baq::print_bits(const unsigned short p, const int a)
{
  int t_=1;
  for(int i=1;i<a;i++)
    t_=t_*2;
  for(int t=t_; t>0; t=t/2)
    if(p & t) std::cout << "1";
    else std::cout << "0";

}

//------------------------------------------------------------------
// compuThreshold(const Charvec, Ivec)
//
// This method is to calculate the threshold of the simulated burst SAR
// data based on the look-up tables. 
//
// The inputs of setThreshold is
//
// echo: SAR data (8-bit char array) with overload data type of float array
// dbg:  flag for debug mode, 1-print out notes, 0-without printing
//                                                 out notes
// The output of setThreshold is 
//
// Threshold: threshold values (integer array with size fixed at 24) 
//------------------------------------------------------------------
void Baq::compuThreshold(const Charvec& echo_8, Ivec& Thresh)
{

  Thresh.resize(Nb_);
  if(dbg==1)
  {
    std::cout<<"------------------ set threshold -----------------"<<std::endl;
    std::cout << "   Baq mode (N_mod_): " << N_mod_ << std::endl;
  }
  //-------------------------------
  // SAR echo digitization and 2D mapping 
  //-------------------------------

  Imat Xecho("Xecho_name",Np_,Ne_);
  unsigned int echo_count=0;
  for(int jj=0;jj<Ne_;jj++)
    for(int ii=0;ii<Np_;ii++)
    {
      Xecho(ii,jj)=echo_8(echo_count);
      echo_count++;
    }
  //std::cout<<"   check number of points: "<<echo_count<< std::endl;

  //---------------------------------
  // calculate vk
  //---------------------------------
  Ivec vk("vk_name",Nb_);
  Ivec nuk("nuk_name",Nb_);

  if(N_mod_==0)
  {
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      nuk(k)= vk(k)/32;
      Thresh(k)=Table8_2(nuk(k));
    }
  }
  else if(N_mod_==6)
    {//scatt mode:first and last  4 samples of the first  and last 8PRI
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<4;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-4;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<4;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-4;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      nuk(k)=vk(k)/16;
      Thresh(k)=Table8_4(nuk(k));
    }
  }
 else if(N_mod_==7)
  {//altimeter: first and last 8 samples of first and last 4 PRIs
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<4;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      for(int jj=Ne_-4;jj<Ne_;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      nuk(k)=vk(k)/16;
      Thresh(k)=Table8_4(nuk(k));
    }
  }
  else if(N_mod_==1)
  {
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
          }
        }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      
      nuk(k)=vk(k)/32;
      Thresh(k)=Table8_2(nuk(k));
    }
  }
  else if((N_mod_==2)||(N_mod_==3)||(N_mod_==4)||(N_mod_==5))
  {
    for(int k=0; k<Nb_; k++)
      Thresh(k)=0;        
  }
  else
    {
      ErrorMessage("There is no such bit number to get threshold.").throwMe();
    }
}


//------------------------------------------------------------------
// compuThreshold(const Fvec, Ivec)
//
// This method is to calculate the threshold of simulated burst SAR
// data based on the look-up tables. 
//
// The inputs of setThreshold is
//
// echo: SAR data (float array) with overloaded data type of char array
// dbg:  flag for debug mode, 1-print out notes, 0-without printing
//                                                 out notes
// The output of setThreshold is 
//
// Threshold: threshold values (integer array with size fixed at 24) 
//------------------------------------------------------------------
void Baq::compuThreshold(const Fvec& echo, Ivec& Thresh)
{

  int NP_echo_8 = echo.size();
  //std::cout << NP_echo_8 << " " << Nt_ << std::endl;
  Charvec echo_8("echo_8_name",NP_echo_8);
  echo_8=Baq::ADC(echo);
  if(NP_echo_8!=Nt_)
    ErrorMessage("The size of SAR data file is not correct.").throwMe();

  Thresh.resize(Nb_);
  if(dbg==1)
  {
    std::cout<<"------------------ set threshold -----------------"<<std::endl;
    std::cout << "   Baq mode (N_mod_): " << N_mod_ << std::endl;
  }
  //-------------------------------
  // SAR echo digitization and 2D mapping 
  //-------------------------------

  Imat Xecho("Xecho_name",Np_,Ne_);
  unsigned int echo_count=0;
  for(int jj=0;jj<Ne_;jj++)
    for(int ii=0;ii<Np_;ii++)
    {
      Xecho(ii,jj)=echo_8(echo_count);
      echo_count++;
    }
  //std::cout<<"   check number of points: "<<echo_count<< std::endl;

  //---------------------------------
  // calculate vk
  //---------------------------------
  Ivec vk("vk_name",Nb_);
  Ivec nuk("nuk_name",Nb_);

  if(N_mod_==0)
  {
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      nuk(k)=vk(k)/32;
      //if(nuk(k)>670) 
      //cout<<nuk(k)<<" cc: "<<vk(k)<<endl;
      Thresh(k)=Table8_2(nuk(k));
    }
  }
  else if(N_mod_==6)
    {//ALTL: 4 samples, 8 echoes
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<4;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-4;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<4;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-4;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
     
      nuk(k)=vk(k)/16;
      Thresh(k)=Table8_4(nuk(k));
    }
  }
  else if(N_mod_==7)
    {//ALTH
      for(int k=0; k<Nb_; k++)
	{
	  echo_count=0;
	  vk(k)=0;
	  
	  for(int jj=0;jj<4;jj++)
	    {
	      for(int ll=0;ll<8;ll++)
		{
		  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
		  echo_count++;
		}
	      for(int ll=Ns_-8;ll<Ns_;ll++)
		{
		  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
		  echo_count++;
		}
	    }
	  for(int jj=Ne_-4;jj<Ne_;jj++)
	    {
	      for(int ll=0;ll<8;ll++)
		{
		  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
		  echo_count++;
		}
	      for(int ll=Ns_-8;ll<Ns_;ll++)
		{
		  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
		  echo_count++;
		}
	    }
	
	  nuk(k)=vk(k)/16;
	  Thresh(k)=Table8_4(nuk(k));
	}
    }
  else if(N_mod_==1)
  {
    for(int k=0; k<Nb_; k++)
    {
      echo_count=0;
      vk(k)=0;
    
      for(int jj=0;jj<8;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
          }
        }
      for(int jj=Ne_-8;jj<Ne_;jj++)
      {
	for(int ll=0;ll<8;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
	for(int ll=Ns_-8;ll<Ns_;ll++)
	{
	  vk(k)+=abs(Xecho(ll+k*Ns_,jj));
          echo_count++;
        }
      }
     
      nuk(k)=vk(k)/32;
     Thresh(k)=Table8_2(nuk(k));
    }
  }
  else if((N_mod_==2)||(N_mod_==3)||(N_mod_==4)||(N_mod_==5))
  {
    for(int k=0; k<Nb_; k++)
      Thresh(k)=0;        
  }
  else
    {
      ErrorMessage("There is no such bit number to get threshold.").throwMe();
    }
}

//--------------------------------------
// 8 to 1 bit threshold lookup table ??
//--------------------------------------

//--------------------------------------
// 8 to 4 bit threshold lookup table   
//--------------------------------------
unsigned int Baq::Table8_4(const unsigned int a)
{
  if((a>=0)&&(a<201))
    {
      if((a>=0)&&(a<=1)) return(1);
      if((a>=2)&&(a<=3)) return(2);
      if((a>=4)&&(a<=6)) return(3);
      if((a>=7)&&(a<=9)) return(4);
      if((a>=10)&&(a<=12)) return(5);
      if((a>=13)&&(a<=15)) return(6);
      if((a>=16)&&(a<=17)) return(7);
      if((a>=18)&&(a<=20)) return(8);
      if((a>=21)&&(a<=23)) return(9);
      if((a>=24)&&(a<=26)) return(10);
      if((a>=27)&&(a<=29)) return(11);
      if((a>=30)&&(a<=32)) return(12);
      if((a>=33)&&(a<=35)) return(13);
      if((a>=36)&&(a<=38)) return(14);
      if((a>=39)&&(a<=41)) return(15);
      if((a>=42)&&(a<=43)) return(16);
      if((a>=44)&&(a<=46)) return(17);
      if((a>=47)&&(a<=49)) return(18);
      if((a>=50)&&(a<=52)) return(19);
      if((a>=53)&&(a<=55)) return(20);
      if((a>=56)&&(a<=58)) return(21);
      if((a>=59)&&(a<=61)) return(22);
      if((a>=62)&&(a<=64)) return(23);
      if((a>=65)&&(a<=67)) return(24);
      if((a>=68)&&(a<=70)) return(25);
      if((a>=71)&&(a<=72)) return(26);
      if((a>=73)&&(a<=75)) return(27);
      if((a>=76)&&(a<=78)) return(28);
      if((a>=79)&&(a<=81)) return(29);
      if((a>=82)&&(a<=84)) return(30);
      if((a>=85)&&(a<=87)) return(31);
      if((a>=88)&&(a<=90)) return(32);
      if((a>=91)&&(a<=93)) return(33);
      if((a>=94)&&(a<=96)) return(34);
      if((a>=97)&&(a<=99)) return(35);
      if((a>=100)&&(a<=101)) return(36);
      if((a>=102)&&(a<=104)) return(37);
      if((a>=105)&&(a<=107)) return(38);
      if((a>=108)&&(a<=110)) return(39);
      if((a>=111)&&(a<=113)) return(40);
      if((a>=114)&&(a<=116)) return(41);
      if((a>=117)&&(a<=119)) return(42);
      if((a>=120)&&(a<=122)) return(43);
      if((a>=123)&&(a<=125)) return(44);
      if((a>=126)&&(a<=128)) return(45);
      if((a>=129)&&(a<=131)) return(46);
      if((a>=132)&&(a<=133)) return(47);
      if((a>=134)&&(a<=136)) return(48);
      if((a>=137)&&(a<=139)) return(49);
      if((a>=140)&&(a<=142)) return(50);
      if((a>=143)&&(a<=145)) return(51);
      if((a>=146)&&(a<=148)) return(52);
      if((a>=149)&&(a<=151)) return(53);
      if((a>=152)&&(a<=154)) return(54);
      if((a>=155)&&(a<=157)) return(55);
      if((a>=158)&&(a<=160)) return(56);
      if((a>=161)&&(a<=162)) return(57);
      if((a>=163)&&(a<=165)) return(58);
      if((a>=166)&&(a<=168)) return(59);
      if((a>=169)&&(a<=171)) return(60);
      if((a>=172)&&(a<=174)) return(61);
      if((a>=175)&&(a<=177)) return(62);
      if((a>=178)&&(a<=180)) return(63);
      if((a>=181)&&(a<=183)) return(64);
      if((a>=184)&&(a<=186)) return(65);
      if((a>=187)&&(a<=189)) return(66);
      if((a>=190)&&(a<=191)) return(67);
      if((a>=192)&&(a<=194)) return(68);
      if((a>=195)&&(a<=197)) return(69);
      if((a>=198)&&(a<=200)) return(70);
    }
  else if((a>=201)&&(a<401))
    {
      if((a>=201)&&(a<=203)) return(71);
      if((a>=204)&&(a<=206)) return(72);
      if((a>=207)&&(a<=209)) return(73);
      if((a>=210)&&(a<=212)) return(74);
      if((a>=213)&&(a<=215)) return(75);
      if((a>=216)&&(a<=218)) return(76);
      if((a>=219)&&(a<=220)) return(77);
      if((a>=221)&&(a<=223)) return(78);
      if((a>=224)&&(a<=226)) return(79);
      if((a>=227)&&(a<=229)) return(80);
      if((a>=230)&&(a<=232)) return(81);
      if((a>=233)&&(a<=235)) return(82);
      if((a>=236)&&(a<=238)) return(83);
      if((a>=239)&&(a<=241)) return(84);
      if((a>=242)&&(a<=244)) return(85);
      if((a>=245)&&(a<=247)) return(86);
      if((a>=248)&&(a<=249)) return(87);
      if((a>=250)&&(a<=252)) return(88);
      if((a>=253)&&(a<=255)) return(89);
      if((a>=256)&&(a<=258)) return(90);
      if((a>=259)&&(a<=261)) return(91);
      if((a>=262)&&(a<=264)) return(92);
      if((a>=265)&&(a<=267)) return(93);
      if((a>=268)&&(a<=270)) return(94);
      if((a>=271)&&(a<=272)) return(95);
      if((a>=273)&&(a<=275)) return(96);
      if((a>=276)&&(a<=278)) return(97);
      if((a>=279)&&(a<=281)) return(98);
      if((a>=282)&&(a<=284)) return(99);
      if((a>=285)&&(a<=287)) return(100);
      if((a>=288)&&(a<=290)) return(101);
      if((a>=291)&&(a<=292)) return(102);
      if((a>=293)&&(a<=295)) return(103);
      if((a>=296)&&(a<=298)) return(104);
      if((a>=299)&&(a<=301)) return(105);
      if((a>=302)&&(a<=304)) return(106);
      if((a>=305)&&(a<=306)) return(107);
      if((a>=307)&&(a<=309)) return(108);
      if((a>=310)&&(a<=312)) return(109);
      if((a>=313)&&(a<=315)) return(110);
      if((a>=316)&&(a<=318)) return(111);
      if((a>=319)&&(a<=320)) return(112);
      if((a>=321)&&(a<=323)) return(113);
      if((a>=324)&&(a<=326)) return(114);
      if((a>=327)&&(a<=329)) return(115);
      if((a>=330)&&(a<=331)) return(116);
      if((a>=332)&&(a<=334)) return(117);
      if((a>=335)&&(a<=337)) return(118);
      if((a>=338)&&(a<=340)) return(119);
      if((a>=341)&&(a<=342)) return(120);
      if((a>=343)&&(a<=345)) return(121);
      if((a>=346)&&(a<=348)) return(122);
      if((a>=349)&&(a<=350)) return(123);
      if((a>=351)&&(a<=353)) return(124);
      if((a>=354)&&(a<=356)) return(125);
      if((a>=357)&&(a<=358)) return(126);
      if((a>=359)&&(a<=361)) return(127);
      if((a>=362)&&(a<=364)) return(128);
      if((a>=365)&&(a<=366)) return(129);
      if((a>=367)&&(a<=369)) return(130);
      if((a>=370)&&(a<=372)) return(131);
      if((a>=373)&&(a<=374)) return(132);
      if((a>=375)&&(a<=377)) return(133);
      if((a>=378)&&(a<=379)) return(134);
      if((a>=380)&&(a<=382)) return(135);
      if((a>=383)&&(a<=384)) return(136);
      if((a>=385)&&(a<=387)) return(137);
      if((a>=388)&&(a<=390)) return(138);
      if((a>=391)&&(a<=392)) return(139);
      if((a>=393)&&(a<=395)) return(140);
      if((a>=396)&&(a<=397)) return(141);
      if((a>=398)&&(a<=400)) return(142);
    }
  else if((a>=401)&&(a<570))
    {
      if((a>=401)&&(a<=402)) return(143);
      if((a>=403)&&(a<=405)) return(144);
      if((a>=406)&&(a<=407)) return(145);
      if((a>=408)&&(a<=409)) return(146);
      if((a>=410)&&(a<=412)) return(147);
      if((a>=413)&&(a<=414)) return(148);
      if((a>=415)&&(a<=417)) return(149);
      if((a>=418)&&(a<=419)) return(150);
      if((a>=420)&&(a<=421)) return(151);
      if((a>=422)&&(a<=424)) return(152);
      if((a>=425)&&(a<=426)) return(153);
      if((a>=427)&&(a<=429)) return(154);
      if((a>=430)&&(a<=431)) return(155);
      if((a>=432)&&(a<=433)) return(156);
      if((a>=434)&&(a<=436)) return(157);
      if((a>=437)&&(a<=438)) return(158);
      if((a>=439)&&(a<=440)) return(159);
      if((a>=441)&&(a<=442)) return(160);
      if((a>=443)&&(a<=445)) return(161);
      if((a>=446)&&(a<=447)) return(162);
      if((a>=448)&&(a<=449)) return(163);
      if((a>=450)&&(a<=451)) return(164);
      if((a>=452)&&(a<=454)) return(165);
      if((a>=455)&&(a<=456)) return(166);
      if((a>=457)&&(a<=458)) return(167);
      if((a>=459)&&(a<=460)) return(168);
      if((a>=461)&&(a<=462)) return(169);
      if((a>=463)&&(a<=465)) return(170);
      if((a>=466)&&(a<=467)) return(171);
      if((a>=468)&&(a<=469)) return(172);
      if((a>=470)&&(a<=471)) return(173);
      if((a>=472)&&(a<=473)) return(174);
      if((a>=474)&&(a<=475)) return(175);
      if((a>=476)&&(a<=477)) return(176);
      if((a>=478)&&(a<=479)) return(177);
      if((a>=480)&&(a<=481)) return(178);
      if((a>=482)&&(a<=483)) return(179);
      if((a>=484)&&(a<=485)) return(180);
      if((a>=486)&&(a<=487)) return(181);
      if((a>=488)&&(a<=489)) return(182);
      if((a>=490)&&(a<=491)) return(183);
      if((a>=492)&&(a<=493)) return(184);
      if((a>=494)&&(a<=495)) return(185);
      if((a>=496)&&(a<=497)) return(186);
      if((a>=498)&&(a<=499)) return(187);
      if((a>=500)&&(a<=501)) return(188);
      if((a>=502)&&(a<=503)) return(189);
      if((a>=504)&&(a<=505)) return(190);
      if((a>=506)&&(a<=507)) return(191);
      if((a>=508)&&(a<=509)) return(192);
      if((a>=510)&&(a<=511)) return(193);
      if((a>=512)&&(a<=513)) return(194);
      if((a>=514)&&(a<=515)) return(195);
      if((a>=516)&&(a<=517)) return(196);
      if((a>=518)&&(a<=518)) return(197);
      if((a>=519)&&(a<=520)) return(198);
      if((a>=521)&&(a<=522)) return(199);
      if((a>=523)&&(a<=524)) return(200);
      if((a>=525)&&(a<=526)) return(201);
      if((a>=527)&&(a<=527)) return(202);
      if((a>=528)&&(a<=529)) return(203);
      if((a>=530)&&(a<=531)) return(204);
      if((a>=532)&&(a<=533)) return(205);
      if((a>=534)&&(a<=535)) return(206);
      if((a>=536)&&(a<=536)) return(207);
      if((a>=537)&&(a<=538)) return(208);
      if((a>=539)&&(a<=540)) return(209);
      if((a>=541)&&(a<=541)) return(210);
      if((a>=542)&&(a<=543)) return(211);
      if((a>=544)&&(a<=545)) return(212);
      if((a>=546)&&(a<=547)) return(213);
      if((a>=548)&&(a<=548)) return(214);
      if((a>=549)&&(a<=550)) return(215);
      if((a>=551)&&(a<=552)) return(216);
      if((a>=553)&&(a<=553)) return(217);
      if((a>=554)&&(a<=555)) return(218);
      if((a>=556)&&(a<=556)) return(219);
      if((a>=557)&&(a<=558)) return(220);
      if((a>=559)&&(a<=560)) return(221);
      if((a>=561)&&(a<=561)) return(222);
      if((a>=562)&&(a<=563)) return(223);
      if((a>=564)&&(a<=564)) return(224);
      if((a>=565)&&(a<=566)) return(225);
      if((a>=567)&&(a<=568)) return(226);
      if((a>=569)&&(a<=569)) return(227);
    }
  else if((a>=570)&&(a<=607))
    {
      if((a>=570)&&(a<=571)) return(228);
      if((a>=572)&&(a<=572)) return(229);
      if((a>=573)&&(a<=574)) return(230);
      if((a>=575)&&(a<=575)) return(231);
      if((a>=576)&&(a<=577)) return(232);
      if((a>=578)&&(a<=578)) return(233);
      if((a>=579)&&(a<=580)) return(234);
      if((a>=581)&&(a<=581)) return(235);
      if((a>=582)&&(a<=583)) return(236);
      if((a>=584)&&(a<=584)) return(237);
      if((a>=585)&&(a<=586)) return(238);
      if((a>=587)&&(a<=587)) return(239);
      if((a>=588)&&(a<=588)) return(240);
      if((a>=589)&&(a<=590)) return(241);
      if((a>=591)&&(a<=591)) return(242);
      if((a>=592)&&(a<=593)) return(243);
      if((a>=594)&&(a<=594)) return(244);
      if((a>=595)&&(a<=595)) return(245);
      if((a>=596)&&(a<=597)) return(246);
      if((a>=598)&&(a<=598)) return(247);
      if((a>=599)&&(a<=600)) return(248);
      if((a>=601)&&(a<=601)) return(249);
      if((a>=602)&&(a<=602)) return(250);
      if((a>=603)&&(a<=604)) return(251);
      if((a>=605)&&(a<=605)) return(252);
      if((a>=606)&&(a<=606)) return(253);
      if((a>=607)&&(a<=607)) return(254);
    }
  else
    {
      cout<<"a is not in the list: "<<a<<endl;
      ErrorMessage("Input value is out of range").throwMe();
    }
  return(0);
}
//--------------------------------------
// 8 to 2 bit threshold lookup table
//--------------------------------------
unsigned int Baq::Table8_2(const unsigned int a)
{
  if((a>=0)&&(a<200))
    {
      if(a==0) return(0);
      if(a==1) return(1);
      if((a>=2)&&(a<=4)) return(2);
      if((a>=5)&&(a<=7)) return(3);
      if((a>=8)&&(a<=10)) return(4);
      if((a>=11)&&(a<=14)) return(5);
      if((a>=15)&&(a<=17)) return(6);
      if((a>=18)&&(a<=20)) return(7);
      if((a>=21)&&(a<=23)) return(8);
      if((a>=24)&&(a<=27)) return(9);
      if((a>=28)&&(a<=30)) return(10);
      if((a>=31)&&(a<=33)) return(11);
      if((a>=34)&&(a<=36)) return(12);
      if((a>=37)&&(a<=40)) return(13);
      if((a>=41)&&(a<=43)) return(14);
      if((a>=44)&&(a<=46)) return(15);
      if((a>=47)&&(a<=49)) return(16);
      if((a>=50)&&(a<=53)) return(17);
      if((a>=54)&&(a<=56)) return(18);
      if((a>=57)&&(a<=59)) return(19);
      if((a>=60)&&(a<=62)) return(20);
      if((a>=63)&&(a<=66)) return(21);
      if((a>=67)&&(a<=69)) return(22);
      if((a>=70)&&(a<=72)) return(23);
      if((a>=73)&&(a<=75)) return(24);
      if((a>=76)&&(a<=79)) return(25);
      if((a>=80)&&(a<=82)) return(26);
      if((a>=83)&&(a<=85)) return(27);
      if((a>=86)&&(a<=88)) return(28);
      if((a>=89)&&(a<=92)) return(29);
      if((a>=93)&&(a<=95)) return(30);
      if((a>=96)&&(a<=98)) return(31);
      if((a>=99)&&(a<=101)) return(32);
      if((a>=102)&&(a<=105)) return(33);
      if((a>=106)&&(a<=108)) return(34);
      if((a>=109)&&(a<=111)) return(35);
      if((a>=112)&&(a<=114)) return(36);
      if((a>=115)&&(a<=118)) return(37);
      if((a>=119)&&(a<=121)) return(38);
      if((a>=122)&&(a<=124)) return(39);
      if((a>=125)&&(a<=127)) return(40);
      if((a>=128)&&(a<=131)) return(41);
      if((a>=132)&&(a<=134)) return(42);
      if((a>=135)&&(a<=137)) return(43);
      if((a>=138)&&(a<=140)) return(44);
      if((a>=141)&&(a<=144)) return(45);
      if((a>=145)&&(a<=147)) return(46);
      if((a>=148)&&(a<=150)) return(47);
      if((a>=151)&&(a<=153)) return(48);
      if((a>=154)&&(a<=157)) return(49);
      if((a>=158)&&(a<=160)) return(50);
      if((a>=161)&&(a<=163)) return(51);
      if((a>=164)&&(a<=166)) return(52);
      if((a>=167)&&(a<=170)) return(53);
      if((a>=171)&&(a<=173)) return(54);
      if((a>=174)&&(a<=176)) return(55);
      if((a>=177)&&(a<=180)) return(56);
      if((a>=181)&&(a<=183)) return(57);
      if((a>=184)&&(a<=186)) return(58);
      if((a>=187)&&(a<=189)) return(59);
      if((a>=190)&&(a<=193)) return(60);
      if((a>=194)&&(a<=196)) return(61);
      if((a>=197)&&(a<=199)) return(62);
    }
  else if((a>=200)&&(a<402))
    {
      if((a>=200)&&(a<=202)) return(63);
      if((a>=203)&&(a<=206)) return(64);
      if((a>=207)&&(a<=209)) return(65);
      if((a>=210)&&(a<=212)) return(66);
      if((a>=213)&&(a<=215)) return(67);
      if((a>=216)&&(a<=219)) return(68);
      if((a>=220)&&(a<=222)) return(69);
      if((a>=223)&&(a<=225)) return(70);
      if((a>=226)&&(a<=228)) return(71);
      if((a>=229)&&(a<=232)) return(72);
      if((a>=233)&&(a<=235)) return(73);
      if((a>=236)&&(a<=238)) return(74);
      if((a>=239)&&(a<=241)) return(75);
      if((a>=242)&&(a<=245)) return(76);
      if((a>=246)&&(a<=248)) return(77);
      if((a>=249)&&(a<=251)) return(78);
      if((a>=252)&&(a<=254)) return(79);
      if((a>=255)&&(a<=258)) return(80);
      if((a>=259)&&(a<=261)) return(81);
      if((a>=262)&&(a<=264)) return(82);
      if((a>=265)&&(a<=267)) return(83);
      if((a>=268)&&(a<=270)) return(84);
      if((a>=271)&&(a<=274)) return(85);
      if((a>=275)&&(a<=277)) return(86);
      if((a>=278)&&(a<=280)) return(87);
      if((a>=281)&&(a<=283)) return(88);
      if((a>=284)&&(a<=286)) return(89);
      if((a>=287)&&(a<=290)) return(90);
      if((a>=291)&&(a<=293)) return(91);
      if((a>=294)&&(a<=296)) return(92);
      if((a>=297)&&(a<=299)) return(93);
      if((a>=300)&&(a<=302)) return(94);
      if((a>=303)&&(a<=305)) return(95);
      if((a>=306)&&(a<=309)) return(96);
      if((a>=310)&&(a<=312)) return(97);
      if((a>=313)&&(a<=315)) return(98);
      if((a>=316)&&(a<=318)) return(99);
      if((a>=319)&&(a<=321)) return(100);
      if((a>=322)&&(a<=324)) return(101);
      if((a>=325)&&(a<=327)) return(102);
      if((a>=328)&&(a<=330)) return(103);
      if((a>=331)&&(a<=333)) return(104);
      if((a>=334)&&(a<=337)) return(105);
      if((a>=338)&&(a<=340)) return(106);
      if((a>=341)&&(a<=343)) return(107);
      if((a>=344)&&(a<=346)) return(108);
      if((a>=347)&&(a<=349)) return(109);
      if((a>=350)&&(a<=352)) return(110);
      if((a>=353)&&(a<=355)) return(111);
      if((a>=356)&&(a<=358)) return(112);
      if((a>=359)&&(a<=361)) return(113);
      if((a>=362)&&(a<=364)) return(114);
      if((a>=365)&&(a<=367)) return(115);
      if((a>=368)&&(a<=369)) return(116);
      if((a>=370)&&(a<=372)) return(117);
      if((a>=373)&&(a<=375)) return(118);
      if((a>=376)&&(a<=378)) return(119);
      if((a>=379)&&(a<=381)) return(120);
      if((a>=382)&&(a<=384)) return(121);
      if((a>=385)&&(a<=387)) return(122);
      if((a>=388)&&(a<=390)) return(123);
      if((a>=391)&&(a<=393)) return(124);
      if((a>=394)&&(a<=395)) return(125);
      if((a>=396)&&(a<=398)) return(126);
      if((a>=399)&&(a<=401)) return(127);
    }
  else if((a>=402)&&(a<570))
    {
      if((a>=402)&&(a<=404)) return(128);
      if((a>=405)&&(a<=406)) return(129);
      if((a>=407)&&(a<=409)) return(130);
      if((a>=410)&&(a<=412)) return(131);
      if((a>=413)&&(a<=415)) return(132);
      if((a>=416)&&(a<=417)) return(133);
      if((a>=418)&&(a<=420)) return(134);
      if((a>=421)&&(a<=423)) return(135);
      if((a>=424)&&(a<=425)) return(136);
      if((a>=426)&&(a<=428)) return(137);
      if((a>=429)&&(a<=431)) return(138);
      if((a>=432)&&(a<=433)) return(139);
      if((a>=434)&&(a<=436)) return(140);
      if((a>=437)&&(a<=438)) return(141);
      if((a>=439)&&(a<=441)) return(142);
      if((a>=442)&&(a<=443)) return(143);
      if((a>=444)&&(a<=446)) return(144);
      if((a>=447)&&(a<=449)) return(145);
      if((a>=450)&&(a<=451)) return(146);
      if((a>=452)&&(a<=454)) return(147);
      if((a>=455)&&(a<=456)) return(148);
      if((a>=457)&&(a<=458)) return(149);
      if((a>=459)&&(a<=461)) return(150);
      if((a>=462)&&(a<=463)) return(151);
      if((a>=464)&&(a<=466)) return(152);
      if((a>=467)&&(a<=468)) return(153);
      if((a>=469)&&(a<=471)) return(154);
      if((a>=472)&&(a<=473)) return(155);
      if((a>=474)&&(a<=475)) return(156);
      if((a>=476)&&(a<=478)) return(157);
      if((a>=479)&&(a<=480)) return(158);
      if((a>=481)&&(a<=482)) return(159);
      if((a>=483)&&(a<=485)) return(160);
      if((a>=486)&&(a<=487)) return(161);
      if((a>=488)&&(a<=489)) return(162);
      if((a>=490)&&(a<=491)) return(163);
      if((a>=492)&&(a<=494)) return(164);
      if((a>=495)&&(a<=496)) return(165);
      if((a>=497)&&(a<=498)) return(166);
      if((a>=499)&&(a<=500)) return(167);
      if((a>=501)&&(a<=502)) return(168);
      if((a>=503)&&(a<=505)) return(169);
      if((a>=506)&&(a<=507)) return(170);
      if((a>=508)&&(a<=509)) return(171);
      if((a>=510)&&(a<=511)) return(172);
      if((a>=512)&&(a<=513)) return(173);
      if((a>=514)&&(a<=515)) return(174);
      if((a>=516)&&(a<=517)) return(175);
      if((a>=518)&&(a<=519)) return(176);
      if((a>=520)&&(a<=521)) return(177);
      if((a>=522)&&(a<=523)) return(178);
      if((a>=524)&&(a<=525)) return(179);
      if((a>=526)&&(a<=527)) return(180);
      if((a>=528)&&(a<=529)) return(181);
      if((a>=530)&&(a<=531)) return(182);
      if((a>=532)&&(a<=533)) return(183);
      if((a>=534)&&(a<=535)) return(184);
      if((a>=536)&&(a<=537)) return(185);
      if((a>=538)&&(a<=539)) return(186);
      if((a>=540)&&(a<=541)) return(187);
      if((a>=542)&&(a<=543)) return(188);
      if((a>=544)&&(a<=545)) return(189);
      if((a>=546)&&(a<=547)) return(190);
      if((a>=548)&&(a<=549)) return(191);
      if((a>=550)&&(a<=551)) return(192);
      if((a>=552)&&(a<=552)) return(193);
      if((a>=553)&&(a<=554)) return(194);
      if((a>=555)&&(a<=556)) return(195);
      if((a>=557)&&(a<=558)) return(196);
      if((a>=559)&&(a<=560)) return(197);
      if((a>=561)&&(a<=562)) return(198);
      if((a>=563)&&(a<=563)) return(199);
      if((a>=564)&&(a<=565)) return(200);
      if((a>=566)&&(a<=567)) return(201);
      if((a>=568)&&(a<=569)) return(202);
    }
  else if((a>=570)&&(a<=645))
    {
      if((a>=570)&&(a<=570)) return(203);
      if((a>=571)&&(a<=572)) return(204);
      if((a>=573)&&(a<=574)) return(205);
      if((a>=575)&&(a<=575)) return(206);
      if((a>=576)&&(a<=577)) return(207);
      if((a>=578)&&(a<=579)) return(208);
      if((a>=580)&&(a<=580)) return(209);
      if((a>=581)&&(a<=582)) return(210);
      if((a>=583)&&(a<=584)) return(211);
      if((a>=585)&&(a<=585)) return(212);
      if((a>=586)&&(a<=587)) return(213);
      if((a>=588)&&(a<=589)) return(214);
      if((a>=590)&&(a<=590)) return(215);
      if((a>=591)&&(a<=592)) return(216);
      if((a>=593)&&(a<=593)) return(217);
      if((a>=594)&&(a<=595)) return(218);
      if((a>=596)&&(a<=596)) return(219);
      if((a>=597)&&(a<=598)) return(220);
      if((a>=599)&&(a<=600)) return(221);
      if((a>=601)&&(a<=601)) return(222);
      if((a>=602)&&(a<=603)) return(223);
      if((a>=604)&&(a<=604)) return(224);
      if((a>=605)&&(a<=606)) return(225);
      if((a>=607)&&(a<=607)) return(226);
      if((a>=608)&&(a<=609)) return(227);
      if((a>=610)&&(a<=610)) return(228);
      if((a>=611)&&(a<=612)) return(229);
      if((a>=613)&&(a<=613)) return(230);
      if((a>=614)&&(a<=614)) return(231);
      if((a>=615)&&(a<=616)) return(232);
      if((a>=617)&&(a<=617)) return(233);
      if((a>=618)&&(a<=619)) return(234);
      if((a>=620)&&(a<=620)) return(235);
      if((a>=621)&&(a<=621)) return(236);
      if((a>=622)&&(a<=623)) return(237);
      if((a>=624)&&(a<=624)) return(238);
      if((a>=625)&&(a<=626)) return(249);
      if((a>=627)&&(a<=627)) return(240);
      if((a>=628)&&(a<=628)) return(241);
      if((a>=629)&&(a<=630)) return(242);
      if((a>=631)&&(a<=631)) return(243);
      if((a>=632)&&(a<=632)) return(244);
      if((a>=633)&&(a<=634)) return(245);
      if((a>=635)&&(a<=635)) return(246);
      if((a>=636)&&(a<=636)) return(247);
      if((a>=637)&&(a<=638)) return(248);
      if((a>=639)&&(a<=639)) return(249);
      if((a>=640)&&(a<=640)) return(250);
      if((a>=641)&&(a<=641)) return(251);
      if((a>=642)&&(a<=643)) return(252);
      if((a>=644)&&(a<=644)) return(253);
      if((a>=645)&&(a<=645)) return(254);
    }
  else
    {
      return(254);
      //cout<<"a is not in the list: "<<a<<endl;
      //ErrorMessage("Input value is out of range").throwMe();
    }
  return(0);
}

//-----------------------------------------------------------------------
// ADC(Fvec)
//
// This method is to (1) normalize the SAR echo data with maximum value
// being 127, and (2) convert the data into char data type with the range
// of -127 to 127.
//-----------------------------------------------------------------------
Charvec Baq::ADC(const Fvec &c)
{
  N_ = c.size();
  unsigned int max_indx;
  Charvec b("digit_name",N_); 
  //Fvec d("temp_name",N_);
  //d=fabs(c);

  float max_c=max(fabs(c),max_indx);
  //std::cout << "max_c" << max_c << std::endl;
  //std::cout << d(max_indx) << std::endl;
  //std::cout << c(max_indx) << std::endl;

  //-------------------------------------------------
  // the following subtitued code is to find max_c
  //-------------------------------------------------
  //float max_c=-1.e20;
  //for(int i=0; i<N_; i++)
  //  {
  //    if(max_c < d(i))
  //      max_c = d(i);
  //  } 

  for(int i=0; i<N_; i++)
    {
      if(c(i)==0)
        b(i)=0;
      else if(c(i)>0)
        b(i)=int(c(i)/max_c*127+0.5);
      else 
        b(i)=int(c(i)/max_c*127-0.5);
      
      //std::cout <<"  b(i)=" <<int(b(i))<<std::endl;
    }
  //std::cout << "   Maximum value (max_c): " << max_c << std::endl;
  return(b);
}


//-----------------
//do not normalize echo data
//------------------
Charvec Baq::ADC_without_normalization(const Fvec& f)
  {
    N_ = f.size();
  
    Charvec b("digit_name",N_); 
    for(int i=0; i<N_; i++)
      b(i)=int(f(i));
    return(b);
  }














































