//----------------------------------------------------------------------------
//
// PerfFile.cpp
//
// This file contains method definitions for the PerfFile handling classes
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------


//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <fstream>
#include "Error.h"
#include "Units.h"
#include "Ambiguity.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "PerfFile.h"
#include "Utils.h"
using std::cout;
using std::endl;
using std::cerr;

//-------------------------
// Methods for class PerfFile
//-------------------------

//--------------
// Constructors: contruct with a filename and file type
//  to write/read record
//--------------
PerfFile::PerfFile(const string& filename, const string& mode) 
  :file_opened_(true),
   file_(filename,mode),
   filetype_(mode) 
  {
  }

//-----------------------
//destructor
//----------------------
PerfFile::~PerfFile()
  {
  if (file_opened_)
    {
    file_.close();
    }
  }


//------------------
//close
//------------------
void PerfFile::close()
  {
    file_.close();
    file_opened_ = false;
  }


//-----------------------------
//readRecord(AmbiguityFNN& amb)
//-----------------------------
void PerfFile::readRecord(AmbiguityFNN& amb)
  {
  if (!file_opened_)
    {
      ErrorMessage("No file is opened to read in").throwMe();
    }
  if (filetype_=="wb"|| filetype_=="w" )
    {
      ErrorMessage("Can't read from outputfile").throwMe();
    }

  //read time
  double read_double;
  file_.read(read_double);
  Time t;
  t.setEt(Uvar(read_double,"s"));
  amb.setTime(t);

  //read beam number
  unsigned int beam_num;
  file_.read(beam_num);
  amb.setBeamNumber(beam_num);

  //read range doppler grid
  unsigned int Ngrid = amb.getGridsize();
  for(unsigned int i = 0; i < Ngrid ; ++i)
    {
      file_.read(read_double);
      amb.range(i) = Uvar(read_double,"km");
    }
  for(unsigned int i = 0; i < Ngrid ; ++i)
    {
      file_.read(read_double);
      amb.doppler(i) = Uvar(read_double,"Hz");
    }
  
  //read number of looks, usable areas
  unsigned int nlook;
  file_.read(nlook);
  amb.setNumberofLooks(nlook);

  file_.read(read_double);
  amb.usable_area_SL = Uvar(read_double,"km km");
  file_.read(read_double);
  amb.usable_area_ML=Uvar(read_double,"km km");
  
  //read multilook amb
  for (unsigned int i = 0; i < Ngrid; ++i)
    {
      file_.read(read_double);
      for (unsigned int j = 0; j < Ngrid; ++j)
	 {
	   amb.amb_ratio_ML(i,j) = Uvar(read_double,"");
	 }
    }

  //read single look snr, amb, and multilook snr
  for (unsigned int i = 0; i < Ngrid; ++i)
  for (unsigned int j = 0; j < Ngrid; ++j)
    {
      file_.read(read_double);
      amb.thermal_snr_SL(i,j) = Uvar(read_double," ");    
      amb.thermal_snr_ML(i,j) = amb.thermal_snr_SL(i,j) * sqrt(double(nlook));

      file_.read(read_double);
      amb.amb_ratio_SL(i,j) = Uvar(read_double,"");

      file_.read(read_double);
      amb.alongtrack(i,j) = Uvar(read_double,"km");
      
      file_.read(read_double);
      amb.crosstrack(i,j) = Uvar(read_double,"km");

      file_.read(read_double);
      amb.oneway_beam_gaindB(i,j) = read_double;
     
      file_.read(read_double);
      amb.echo_power(i,j) = Uvar(read_double,"kg km km/(s s s)");
    }
  for (unsigned int i = 0; i < Ngrid+1; ++i)
  for (unsigned int j = 0; j < Ngrid+1; ++j)
    {
      //read in alongtrack and cross track grid values (pixel+1)
      file_.read(read_double);
      amb.alongtrack(i,j) = Uvar(read_double,"km");
      file_.read(read_double);
      amb.crosstrack(i,j) = Uvar(read_double,"km");
    }
  }


//------------------------
//writeRecord(const AmbiguityFNN& amb): write fast/slow field into file
//----------------------
void PerfFile::writeRecord(const AmbiguityFNN& amb) 
  {
  if (!file_opened_)
    {
    throw ErrorMessage("No file is opened to write");
    }
  if (filetype_=="rb"|| filetype_=="r" )
    {
    throw ErrorMessage("Can't write to inputfile");
    }
 
  //write time
  Time t = amb.getTime();
  file_.write(t.et().getInUnits("s"));

  //write beam number
  unsigned int beam_num = amb.getBeamNumber();
  file_.write(beam_num);

  //write range,doppler grid
  unsigned int Ngrid = amb.getGridsize();
  for(unsigned int i = 0; i < Ngrid ; ++i)
    {
      file_.write(amb.range(i).getInUnits("km"));
    }
  for(unsigned int i = 0; i < Ngrid ; ++i)
    {
      file_.write(amb.doppler(i).getInUnits("Hz"));
    }

  //write number of looks, usable areas
  unsigned int nlook = amb.getNumberofLooks();
  file_.write(nlook);
  file_.write(amb.usable_area_SL.getInUnits("km km"));
  file_.write(amb.usable_area_ML.getInUnits("km km"));
  
  //write multilook amb
  for (unsigned int i = 0; i < Ngrid; ++i)
    {
      unsigned int j = (unsigned int) Ngrid/2;
      file_.write(amb.amb_ratio_ML(i,j).getInUnits(""));
    }

  //write single look snr and amb
  for (unsigned int i = 0; i < Ngrid; ++i)
  for (unsigned int j = 0; j < Ngrid; ++j)
    {
      file_.write(amb.thermal_snr_SL(i,j).getInUnits(""));
      file_.write(amb.amb_ratio_SL(i,j).getInUnits(""));
      file_.write(amb.alongtrack(i,j).getInUnits("km"));
      file_.write(amb.crosstrack(i,j).getInUnits("km"));
      file_.write(amb.oneway_beam_gaindB(i,j).getInUnits(""));
      file_.write(amb.echo_power(i,j).getInUnits("kg km km/(s s s)"));
    }
  for (unsigned int i = 0; i < Ngrid+1; ++i)
  for (unsigned int j = 0; j < Ngrid+1; ++j)
    {
      //write alongtrack and cross track grid values (pixel+1)
      file_.write(amb.alongtrack(i,j).getInUnits("km"));
      file_.write(amb.crosstrack(i,j).getInUnits("km"));
    }
  } 


//---------------------------
//return eof of FileMgr
//--------------------------
bool PerfFile::eof()
  {
  return(file_.eof());
  }



