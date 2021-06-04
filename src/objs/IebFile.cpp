//----------------------------------------------------------------------------
//
// IebFile.cpp
//
// This file provides Ieb read/write services
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
#include "Ieb.h"
#include "IebFile.h"

using std::cout;
using std::endl;
using std::cerr;

//-------------------------
// Methods for class Ieb
//-------------------------

//--------------
// Constructors: contruct with a filename and file type
//  to write/read record
//--------------
IebFile::IebFile(const string& filename, const string& mode) 
  :slowfastfield_loaded_(false),
   file_opened_(true),
   ieb_all_written_(false),
   ieb_all_read_(false),
   file_(filename,mode),
   filetype_(mode)
  {
  slowfield_.resize(Nslowfield);
  fastfield_.resize(Nfastfield);
  ieb_list_.clear();
  }

//-----------------------
//destructor
//----------------------
IebFile::~IebFile()
  {
  if (file_opened_)
    {
    file_.close();
    }
  }

//--------------------
//void load(const Ieb& ieb)
//------------------
void IebFile::load(const Ieb& ieb)
  {
    ieb_list_.push_back(ieb);
  }

//-----------------------------
//readAllIebRecords(): read all the record in the file
//-----------------------------
void IebFile::readAllIebRecords() 
  {
  if (!file_opened_)
    {
      ErrorMessage("No file is opened to read in").throwMe();
    }
  if (filetype_=="wb"|| filetype_=="w" )
    {
      ErrorMessage("Can't read from outputfile").throwMe();
    }
  if(ieb_all_read_)
    {
      ErrorMessage("All the records have been already read in").throwMe();
    }
  
  while(!file_.eof())
    {
      double et_time_in_sec;
      file_.read(et_time_in_sec);
      
      Ieb ieb(Uvar(et_time_in_sec,"s"));
      
      for(unsigned int i_slow=0; i_slow < Nslowfield;++i_slow){
	file_.read(slowfield_[i_slow]);}
            
      for(unsigned int i_fast=0;i_fast<Nfastfield;++i_fast){
	file_.read(fastfield_[i_fast]);}
      ieb.loadSlowFastfield(slowfield_,fastfield_);  
      ieb_list_.push_back(ieb);    
    }
 
  //check whether time information is stored from past to future order
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < -- ieb_list_.end();++ieb_pt)
    {
      if ( (*ieb_pt).getTime() > (*(ieb_pt+1)).getTime()) ErrorMessage("IebFile: time is not ordered from past to future").throwMe();
    }
  

  ieb_all_read_= true;
  cout<<"read all the records"<<endl;
  }

//------------------------
//writeRecord(): write time ordered fast/slow field into file
//----------------------
void IebFile::writeAllIebRecords() throw(ErrorMessage)
  {
  if (!file_opened_)
    {
    throw ErrorMessage("No file is opened to read in");
    }
  if (filetype_=="rb"|| filetype_=="r" )
    {
    throw ErrorMessage("Can't write to inputfile");
    }
  if (ieb_all_written_)
    {
      ErrorMessage("ieb records have been written").throwMe();
    }
  //--------------------------------
  //format: time, slowfield, fastfield
  //---------------------------------

  //check time order
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < --ieb_list_.end();++ieb_pt)
    {
      if ((*ieb_pt).getTime() > (*(ieb_pt+1)).getTime()) ErrorMessage("IebFile: time is not ordered from past to future").throwMe();
    }
  
  for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < ieb_list_.end();++ieb_pt)
    {
      
      file_.write(((*ieb_pt).getTime()).getInUnits("s"));
     
      slowfield_ = (*ieb_pt).getSlowfield();
      fastfield_ = (*ieb_pt).getFastfield();
      
      for(unsigned int i_slow=0; i_slow < Nslowfield;++i_slow){
	file_.write(slowfield_[i_slow]);}
	
      for(unsigned int i_fast=0;i_fast<Nfastfield;++i_fast){
	file_.write(fastfield_[i_fast]);}
    }
  ieb_all_written_ = true;
  } 

//--------------------------
//Return Ieb for a given time
//---------------------------
Ieb IebFile::getIeb(const Time& t)
  {
    if (!ieb_all_read_)
      {
	ErrorMessage("ieb records have not  been read in").throwMe();
      }
    
    vector<Ieb>::iterator ieb_return;
    
    //if the request time is earlier than the first record, throw an error
    ieb_return = ieb_list_.begin();
    if (t < (*ieb_return).getTime()) ErrorMessage("Requested time is earlier than first ieb record in the file").throwMe();

    //if the request time is later than the last record, return the last record
    ieb_return = --ieb_list_.end();
    if (t > (*ieb_return).getTime())
      {
	return(*ieb_return);
      }
   
    //if the request time is between the first and last, find it
    ieb_return = ieb_list_.begin();
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < --ieb_list_.end();++ieb_pt)
      {
	ieb_return = ieb_pt;
	ieb_return++;
	if ( t >= (*ieb_pt).getTime() && t < (*ieb_return).getTime()) 
	  {
	    break;
	  }  
      }
    ieb_return--;
    return(*ieb_return);
  }

//------------------------------------
//{add,delete} ieb
//------------------------------------
void IebFile::addIeb(const Ieb& ieb)
  {
    if (!ieb_all_read_)
      {
	ErrorMessage("ieb records have not  been read in").throwMe();
      } 

    //when no record
    if (ieb_list_.size() == 0)
      {
	ieb_list_.push_back(ieb);
	return;
      }
  
    //when ieb.getTime is earlier than ieb_list ime
    vector<Ieb>::iterator  ieb_add;
    ieb_add = ieb_list_.begin();
    if (ieb.getTime() < (*ieb_add).getTime())
      {
	ieb_list_.insert(ieb_add,ieb);//insert one before begin
	return;
      }
      
    //when ieb.getTime is later than ieb_list time
    ieb_add = --ieb_list_.end();
    if (ieb.getTime() > (*ieb_add).getTime())
      {
	ieb_list_.push_back(ieb);
	return;
      }

   
    //when ieb.getTime is between record [0] and [N-1]
    ieb_add = ieb_list_.begin();
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < --ieb_list_.end();++ieb_pt)
      {
	if (ieb.getTime() == (*ieb_pt).getTime()) ErrorMessage("There is an ieb record at the same time").throwMe();
	ieb_add = ieb_pt;
	ieb_add++;
	if (ieb.getTime() > (*ieb_pt).getTime() && ieb.getTime() < (*ieb_add).getTime())
	  {
	    break;
	  }
      }
    ieb_list_.insert(ieb_add,ieb);//add a record before ieb_add   
  }

void IebFile::deleteIeb(const Ieb& ieb)
  {
    if (!ieb_all_read_)
      {
	ErrorMessage("ieb records have not  been read in").throwMe();
      } 

    //check there is a record at ieb.getTime();
    unsigned int record_found = 0;
    vector<Ieb>::iterator ieb_delete;
  
    for (vector<Ieb>::iterator ieb_pt = ieb_list_.begin(); ieb_pt < ieb_list_.end();++ieb_pt)
      {
	if (ieb.getTime() == (*ieb_pt).getTime())
	  {
	    record_found = 1;
	    ieb_delete = ieb_pt;
	    break;
	   
	  }
      }
    if (record_found != 1) ErrorMessage("There is no record to delete");
    
    //if there is a record, delete
    ieb_list_.erase(ieb_delete);
  }



//----------------------------
//Return the size of ieb records, start and end time
//-----------------------------

Time IebFile::getStartTime()
  {
    if (ieb_list_.size() == 0) ErrorMessage("IebFile::getStartTime: No record avail").throwMe();
    return( (*ieb_list_.begin()).getTime());
  }

Time IebFile:: getEndTime()
  {
  if (ieb_list_.size() == 0) ErrorMessage("IebFile::getEndTime: No record avail").throwMe();
  return( (*(--ieb_list_.end())).getTime());
  }

unsigned int IebFile::getNumberofRecords()
  {
    return(ieb_list_.size());
  }
