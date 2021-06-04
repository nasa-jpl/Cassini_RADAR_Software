//------------------------------------------------------------------
// Baq.h
//
// This file contains the Baq methods.
// The Baq class provides following processes:
//   (1) load parameters from Ieb:                Baq::setIeb()
//   (2) Digitizing anologue SAR data:            Baq::ADC()
//   (3) Calculate threshold:                     Baq::setThreshold()
//   (4) Encoding digitized SAR:                  Baq::Encode_Nbit()
//   (5) Decoding SAR data:                       Baq::Decode_Nbit()
//
//------------------------------------------------------------------



#ifndef Baq_H
#define Baq_H

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
class Baq;
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Config.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Beam.h"
#include "Ieb.h"
#include "Sab.h"

using std::string;


//-------------------------
// Class Baq declaration
//-------------------------
class Baq
  {
  public:
    Baq();

    void debug_mode();
    Charvec ADC(const Fvec& f);
    Charvec ADC_without_normalization(const Fvec& f);
    void setIeb(Ieb& ieb);
    void setParams(const Uvar& pri, const unsigned int& pul, 
		   const unsigned int& baq_mode, const Uvar& SR, 
		   const int& tro_int);
    void compuThreshold(const Charvec& c, Ivec& th);  
    void compuThreshold(const Fvec& c, Ivec& th);  
 
    void Encode_Nbit(const Charvec& c,const Ivec& th,Charvec& words);
    void Encode_Nbit(const Fvec& c, const Ivec& th,Charvec& words);

    void Encode_thresh(const Ivec& th, Array1D<unsigned short>& thresh);


    void Decode_Nbit(const Charvec& w, const Ivec& th, Fvec& c);

    float minVariance(const Ivec& th);
    //methods not being used
    void LoadEncodedata(SUvec& w);
    void addHeaderFooter(const Charvec& words, const Ivec& th);

  private:
    int N_,Np_,Nt_,Nb_,Ns_,Ns_p_,Ne_,N_mod_,MN_words_,dbg;
    unsigned int Table8_2(const unsigned int a);
    unsigned int Table8_4(const unsigned int a);
    void print_bits(const unsigned short p, const int a);
    void print_bits_endl(const unsigned short p, const int a);
  };
#endif


















































