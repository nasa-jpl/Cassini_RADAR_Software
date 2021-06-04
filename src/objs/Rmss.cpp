//----------------------------------------------------------------------------
//
// Rmss.cpp
//
// This file contains method definitions for the Rmss  handling classes
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
#include <fstream>
#include "Error.h"
#include "Utils.h"
#include "Units.h"
#include "Ieb.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "Io.h"
#include "Rmss.h"



//-------------------------
// Methods for class 
//-------------------------

//--------------
// Constructors
//--------------
 

//--------------------------------------------------------------------------
// Rmss()
//---------------------------------------------------------------------------
Rmss::Rmss(const string& filename, const string& filetype) throw(ErrorMessage)
  :ieb_sequence("ieb_sequence"),
   filename_(filename),
   ft_(filetype)
 {
   if (filetype =="wb" || filetype =="rb")
     {
       FileMgr fp_(filename,filetype); 
     }
   else
     {
       throw ErrorMessage("File type should be either wb or rb");
     }
   ieb_sequence.resize(Nmax);
   ieb_counter_ = 0;
 }

//-----------------------
//Load Ieb
//-----------------------
void Rmss::loadIeb(const Ieb& ieb) throw(ErrorMessage)
  {
  if (ieb_counter_ > Nmax) throw ErrorMessage("Too many records");
  //-----------------
  //transfer data
  //----------------
  ieb_sequence(ieb_counter_) = ieb;
  ieb_counter_=ieb_counter_ + 1;
  }
//---------------------------
//readRecord
//--------------------------
void Rmss::readRecordsfromFile()throw(ErrorMessage)
  {
  if (ft_=="wb")
    throw ErrorMessage("Can't read in output file: "+filename_); 
  } 

//------------------------------
//writeRecord
//------------------------------
void Rmss::writeRecordstoFile()throw(ErrorMessage)  
  {
  if (ft_=="rb")
    throw ErrorMessage("Can't write to inputfile: "+filename_);
 
  for (unsigned int i_ieb=0; i_ieb<ieb_counter_;++i_ieb)
    {
      Time t= ieb_sequence(i_ieb).getTime();
      string time_string;
      time_string.resize(23);
      
      time_string=t.utc("ISOD");
      cout<<"time_string "<<time_string<<endl;

    }
  }  


//----------------------
//close file
//---------------------
void Rmss::close()
{
  fp_.close();
}



