//----------------------------------------------------------------------------
//
// AmbiguityFile.cpp
//
// This file contains method definitions that handle ambiguity data 
//----------------------------------------------------------------------------



//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>



//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <fstream>
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "AmbiguityFile.h"
#include "Io.h"
#include "Utils.h"

using std::cout;
using std::endl;
using std::cerr;


//-------------------------
// Methods for class AmbiguityFile
//-------------------------

//--------------
// Constructors: contruct with a filename and file type
//  to write/read record
//--------------
AmbiguityFile::AmbiguityFile(const string& filename,const string& mode)
  :ambiguity_data_loaded_(false),
   file_opened_(true),
   file_(filename,mode),
   filetype_(mode)
  {
  }



//-----------------------
//destructor
//----------------------
AmbiguityFile::~AmbiguityFile()
  {
  if (file_opened_)
    {
    file_.close();
    }
  }

//--------------------
//void load(const Ambiguity& amb)
//------------------
void AmbiguityFile::load(const Ambiguity& amb)
  {
  ambiguity_data_loaded_ = true;
  amb_ = amb;
  }

//-----------------------------
//readRecord(Ambiguity& amb) 
//-----------------------------
void AmbiguityFile::readRecord(Ambiguity& amb) 
  {
  if (!file_opened_)
    {
    ErrorMessage e("No file is opened to read in");
    e.throwMe();
    }
  if (filetype_=="wb"|| filetype_=="w" )
    {
    ErrorMessage e("Can't read from outputfile");
    e.throwMe();
    } 
  
  string time_string;
  time_string.resize(21);
  unsigned int beam_num;
  double beam_max_gain;
  double x_input;
  unsigned int Ndata;
  unsigned int  mirror_geom_indicator=0;
  Uvar lambda;

  file_.read(time_string);
  Time t = Time(time_string);
  amb.setTime(t);

  file_.read(beam_num);
  amb.setBeamNumber(beam_num);

  file_.read(beam_max_gain);
  amb.setBeamMaxGain(beam_max_gain);

  file_.read(x_input);
  amb.range_bore=Uvar(x_input,"km");

  file_.read(x_input);
  amb.doppler_bore=Uvar(x_input,"Hz");

  file_.read(x_input);
  amb.altitude = Uvar(x_input,"km");

  file_.read(x_input);
  amb.boresight_incidence_angle = Uvar(x_input,"rad");

  file_.read(x_input);
  lambda = Uvar(x_input,"km");
  amb.setWavelength(lambda);

  file_.read(Ndata);
  amb.setGridsize(Ndata);

  file_.read(mirror_geom_indicator);
  if (mirror_geom_indicator == 0)
    {
      amb.setMirrorAmbiguity(false);
    }
  else if (mirror_geom_indicator ==1)
    {
      amb.setMirrorAmbiguity(true);
    }

  //read all the public varaibles into output file
  
  for (unsigned int i = 0; i < Ndata; ++i)
    {
    file_.read(x_input);
    amb.range(i)= Uvar(x_input,"km");
    file_.read(x_input);
    amb.doppler(i)= Uvar(x_input,"Hz");
    }
  
  
  //area incidenceangle along/cross data set are initiallized with base unit
  for (unsigned int i = 0; i < Ndata; ++i)
  for (unsigned int j = 0; j < Ndata;++j)
    {
    file_.read(amb.oneway_beam_gain(i,j));

    file_.read(x_input);
    amb.incidenceangle(i,j).setValue(x_input);
    
    file_.read(x_input);
    amb.area(i,j).setValue(x_input);
   
    file_.read(amb.number_of_solution(i,j));
    }

 
  for (unsigned int i = 0; i < Ndata+1; ++i)
  for (unsigned int j = 0; j < Ndata+1;++j)
    {    
      file_.read(x_input);
      amb.alongtrack(i,j).setValue(x_input);
      file_.read(x_input);
      amb.crosstrack(i,j).setValue(x_input);
    } 

  if(mirror_geom_indicator!= 1) {return;}
  // if there is no data for mirror geom, stop reading
 
  for (unsigned int i = 0; i < Ndata; ++i)
  for (unsigned int j = 0; j < Ndata;++j)
    {
    file_.read(amb.oneway_beam_gain_mirror(i,j));

    file_.read(x_input);
    amb.incidenceangle_mirror(i,j).setValue(x_input);
   
    file_.read(x_input);
    amb.area_mirror(i,j).setValue(x_input);
    
    file_.read(amb.number_of_solution_mirror(i,j));
    } 
  }

//------------------------
//writeRecord(): write public data into file
//----------------------
void AmbiguityFile::writeRecord() 
  {
  if (!file_opened_)
    {
    ErrorMessage e("No file is opened to read in");
    e.throwMe();
    }
  if (filetype_=="rb"|| filetype_=="r" )
    {
    ErrorMessage e("Can't write to inputfile");
    e.throwMe();
    }
  if (!ambiguity_data_loaded_)
    {
    ErrorMessage e("No ambiguity data set is loaded ");
    e.throwMe();
    }

  Time t = amb_.getTime();
  unsigned int beam_num = amb_.getBeamNumber();
  double beam_max_gain = amb_.getBeamMaxGain();
  unsigned int Ndata = amb_.getGridsize();
  Uvar lambda = amb_.getWavelength();
  bool mirror_geom_set=amb_.getMirrorAmbiguity();
  unsigned int mirror_geom_indicator =0;
  if (mirror_geom_set) 
    {
    mirror_geom_indicator = 1;
    }

 
  string time_string;
  time_string.resize(21);
  time_string = t.utc("ISOD");

  file_.write(time_string);
  file_.write(beam_num);
  file_.write(beam_max_gain);
  file_.write(amb_.range_bore.getInUnits("km"));
  file_.write(amb_.doppler_bore.getInUnits("Hz"));
  file_.write(amb_.altitude.getInUnits("km"));
  file_.write(amb_.boresight_incidence_angle.getInUnits("rad"));
  file_.write(lambda.getInUnits("km"));
  file_.write(Ndata);
  file_.write(mirror_geom_indicator);

  //write all the public varaibles into output file
  for (unsigned int i = 0; i < Ndata; ++i)
    {
    file_.write(amb_.range(i).getInUnits("km"));
    file_.write(amb_.doppler(i).getInUnits("Hz"));
    }

  //force to use base units
  for (unsigned int i = 0; i < Ndata; ++i)
  for (unsigned int j = 0; j < Ndata;++j)
    {
    file_.write(amb_.oneway_beam_gain(i,j));
    file_.write(get_in_base_units(amb_.incidenceangle(i,j)));
    file_.write(get_in_base_units(amb_.area(i,j)));
    file_.write(amb_.number_of_solution(i,j));
    }
  for (unsigned int i = 0; i < Ndata+1; ++i)
  for (unsigned int j = 0; j < Ndata+1;++j)
    {    
      file_.write(get_in_base_units(amb_.alongtrack(i,j)));
      file_.write(get_in_base_units(amb_.crosstrack(i,j)));
    }   

  if(mirror_geom_indicator!= 1) {return;}
  // if there is no data for mirror geom, stop writing and return
  for (unsigned int i = 0; i < Ndata; ++i)
  for (unsigned int j = 0; j < Ndata;++j)
    {
    file_.write(amb_.oneway_beam_gain_mirror(i,j));
    file_.write(get_in_base_units(amb_.incidenceangle_mirror(i,j)));
    file_.write(get_in_base_units(amb_.area_mirror(i,j)));
    file_.write(amb_.number_of_solution_mirror(i,j));
    }
  } 



//---------------------------
//return eof of FileMgr
//--------------------------
bool AmbiguityFile::eof()
  {
  return(file_.eof());
  }






