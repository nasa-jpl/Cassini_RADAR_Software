//----------------------------------------------------------------------------
// Ckernel.cpp
//
// This file contains method definitions for the Ckernel class.
// Ckernel provides support for Ckernel manipulations.
// In particular, it supports writing new ckernel files, and determining
// the valid time range of existing ckernel files.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char ckernel_c[] =
  "@(#) $Id: Ckernel.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "Ivd.h"
#include "Array.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Frame.h"
#include "Time.h"
#include "Constants.h"
#include "SpiceUsr.h"
#include "Ckernel.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//------------------------//
// Class Ckernel  Methods //
//------------------------//

// Constructors

Ckernel::Ckernel() 
  :quats("quaternions"), 
   avvs("avvs"),
   attitude_x("attitude_x"),
   attitude_mz("attitude_mz"),
   et_x("et_x"), 
   et_mz("et_mz"),
   sclkdp_x("sclkdp_x"),
   sclkdp_mz("sclkdp_mz"),
   ivd_data_loaded_(false),
   quats_data_loaded_(false),
   INST_(-82000),
   REF_("J2000"),
   CENTERNAME_("Titan"),
   SCNAME_("CASSINI"),
   SEGID_("CASSINI PREDICT TYPE 3 SEGMENT")
  {
  quats.resize(maxrec,4);
  avvs.resize(maxrec,3);
  attitude_x.resize(maxrec,3);
  attitude_mz.resize(maxrec,3);
  et_x.resize(maxrec);
  et_mz.resize(maxrec);
  sclkdp_x.resize(maxrec);
  sclkdp_mz.resize(maxrec);
  quats = 0.0;
  avvs = 0.0;
  attitude_x = 0.0;
  attitude_mz = 0.0;
  et_x = Uvar(0.0,"s");
  et_mz = Uvar(0.0,"s");
  sclkdp_x = 0.0;
  sclkdp_mz = 0.0;
  N_dat = 0;
  AVFLAG_ = SPICETRUE;
  }


//-----------------------------------------------
// LoadIvd(const Ivd& ivd) 
//_______________________________________________

void  Ckernel::LoadIvd(const Ivd& ivd) 
  {
  if (ivd.ivd_data_loaded == false) 
    {
    ErrorMessage e("Ivd data is not loaded"); e.throwMe();
    return;
    }
  quats= ivd.quats ;
  avvs= ivd.avvs;
  attitude_x= ivd.attitude_x;
  attitude_mz= ivd.attitude_mz;
  et_x= ivd.et_x;
  et_mz =ivd.et_mz;
  sclkdp_x= ivd.sclkdp_x;
  sclkdp_mz= ivd.sclkdp_mz ;
  N_dat= ivd.N_dat ;
   
  start_x=ivd.start_x;
  end_x=ivd.end_x;
  start_mz= ivd.start_mz;
  end_mz= ivd.end_mz;
 
  if (start_x != start_mz) 
    {
    ErrorMessage e("Start time mismatch"); e.throwMe();
    return;
    }
  if (end_x != end_mz) 
    {
    ErrorMessage e("End time mismatch"); e.throwMe();
    return;
    }

  N_dat= ivd.N_dat;
  ivd_data_loaded_ = true;
  }

//----------------------------------------------
// void GenerateQuats() 
// 
//----------------------------------------------

void Ckernel::GenerateQuats() 
  {
  if (ivd_data_loaded_ == false) 
    {
    ErrorMessage e("ivd data is not loaded"); e.throwMe();
    return;
    }


  for (unsigned int i = 0; i <N_dat; i++)
    {
    Frame rotation("J2000","Titan");
    Time t_x,t_mz;
    t_x.setEt(et_x(i));
    t_mz.setEt(et_mz(i));
    Dvec quats_index("Quaternion values");
    quats_index.resize(4);
    DirectionVector dir_x("dir_x",rotation,t_x,
      attitude_x(i,0),attitude_x(i,1),attitude_x(i,2));
    DirectionVector dir_z("dir_z",rotation,t_mz,
      -attitude_mz(i,0),-attitude_mz(i,1),-attitude_mz(i,2));
    DirectionVector dir_y("dir_y",cross(dir_z,dir_x));
    quats_index=rotation.quaternion(dir_x,dir_y,dir_z);
      
    quats(i,0) = quats_index(0);
    quats(i,1) = quats_index(1);
    quats(i,2) = quats_index(2);
    quats(i,3) = quats_index(3);
    avvs(i,0) = avvs(i,1)=avvs(i,2) = 0.0;
    //cout<<i<<" "<<quats(i,0)<<" "<<quats(i,1)<<" "<<quats(i,2)<<" "<<quats(i,3)<<endl;
    quats_data_loaded_ = true;
    }

  }

//------------------------------------/--------------------
//void Ckernel::WriteData() 
//--------------------------------------------------------

void Ckernel::WriteData(string& ckernel_file) 
  {
  if(quats_data_loaded_ == false)
    {
    ErrorMessage e("quaternions are not loaded"); e.throwMe();
    return;
    }

  //--------------------------------------------------------
  //C-kernel writing routine
  //preset variables
  //nint = 1, 
  //-------------------------------------------------------
    
  //The following definitions are based on the information regarding 
  //data type3 writing method described in ck.req
  // no character in comment area ->NCOMCH = 0
  //Input file name ->    IFNAME = ck
   
  SpiceDouble begtime,endtime;
  SpiceDouble spice_start_x_time[maxrec],spice_start_mz_time[maxrec];
  SpiceDouble spice_sclkdp[maxrec];
  SpiceDouble spice_quats[maxrec][4];
  SpiceDouble spice_avvs[maxrec][3];
  
  begtime=start_x.encodedSclk(SCNAME_); //encoded Sclk
  endtime=end_x.encodedSclk(SCNAME_); //encoded Sclk
  int nrec = N_dat; //total number of record
  unsigned int NCOMCH = 0; // no comments in the c-kernel file
   
  for (unsigned int i = 0; i <N_dat; i++)
    {
      spice_sclkdp[i] = sclkdp_x(i);     
      spice_quats[i][0]=quats(i,0) ;
      spice_quats[i][1]=quats(i,1) ;
      spice_quats[i][2]=quats(i,2) ;
      spice_quats[i][3]=quats(i,3) ;
      spice_avvs[i][0]= avvs(i,0);
      spice_avvs[i][1]=avvs(i,1);
      spice_avvs[i][2]=avvs(i,2);  
    }

  AVFLAG_ = SPICETRUE;

  spice_start_x_time[0]=sclkdp_x(0);
  spice_start_mz_time[0] = sclkdp_mz(0);
   
  spice_write_ckernel(ckernel_file, NCOMCH, begtime,endtime,INST_,REF_,
                      AVFLAG_,SEGID_,nrec,spice_sclkdp,spice_quats,
                      spice_avvs,spice_start_x_time);
   
    /********************************************** 
    cout<<"ncom "<< NCOMCH<<endl;
    cout<<"beg time "<<begtime<<endl;
    cout<<"end time "<< endtime<<endl;
    cout<<"Ins "<<INST_<<endl;
    cout<<"Ref "<<REF_<<endl;
    cout<<"Av flag"<<AVFLAG_<<endl;
    cout<<"Seg id "<<SEGID_<<endl;
    cout<<"nrec "<<nrec<<endl;
    for (unsigned int i = 0; i < N_dat; i++)
    {
      cout<< "sclkdp "<<spice_sclkdp[i]<<endl;
      cout<<setw(15)<<spice_quats[i][0]<<setw(15)<<" "<<spice_quats[i][1]<<endl;
      cout<<setw(15)<<spice_quats[i][2]<<setw(15)<<" "<<spice_quats[i][3]<<endl;
      cout<<setw(15)<<spice_avvs[i][0]<<setw(15)<<spice_avvs[i][1]<<setw(15)<<spice_avvs[i][2]<<endl;
    }
    cout<<"nint "<<nint<<endl;
    cout<<"spice start x time "<<spice_start_x_time[0]<<endl;
    cout<<"spice true "<<SPICETRUE<<endl;
    **************************************************/

  }  

//-----------------------------------------------------------------------------
// validTimeRange(filename,sc_name,start_time,end_time)
//
// Static method to determine the valid time range (inclusive) for the ckernel
// identified by filename.  The ckernel is defined for the spacecraft
// identified by sc_name.
//-----------------------------------------------------------------------------

void Ckernel::validTimeRange(const string& filename, const string& sc_name,
  Time& start_time, Time& end_time)
  {
  SpiceInt handle;
  SpiceBoolean found = false;
  dafopr_c(filename.c_str(),&handle);
  dafbfs_c(handle);
  daffna_c(&found);
  if (failed_c())
    {
    GeomError e(GeomError::spice_misc);
    e.throwMe();
    }
  SpiceDouble summary[5];
  SpiceDouble dsum[2];
  SpiceInt isum[6];
  SpiceDouble start_t,end_t;
  bool first = true;
  while (found)
    {
    dafgs_c(summary);
    dafus_c(summary,2,6,dsum,isum);
    if (first)
      {
      start_t = dsum[0];
      end_t = dsum[1];
      }
    if (dsum[0] < start_t) start_t = dsum[0];
    if (dsum[1] > end_t) end_t = dsum[1];
    daffna_c(&found);
    first = false;
    if (failed_c())
      {
      GeomError e(GeomError::spice_misc);
      e.throwMe();
      }
    }
  start_time.setEncodedSclk(sc_name,start_t);
  end_time.setEncodedSclk(sc_name,end_t);
  }

