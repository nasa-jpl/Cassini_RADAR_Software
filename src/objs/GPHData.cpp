#include <string.h>
#include <stdlib.h>
#include"GPHData.h"

using std::cout;
using std::cerr;
using std::endl;

//--------------------
// Methods for GPHData
//--------------------

//--------------
// Constructors
//--------------
GPHData::GPHData()
: filename_(""),
  records_counted_(0),
  file_read_(false),
  header_handled_(false),
  temp_("temp_",1),
  time_("time_",1)
{}

GPHData::GPHData(const string& filename) 
  : filename_(filename),
    file_(filename,"r"),
    records_counted_(0),
    file_read_(false),
    header_handled_(false),
    temp_("temp_",1),
    time_("time_",1)
    
{
  recordCount();
  temp_.resize(records_counted_);
  time_.resize(records_counted_);
  read();
}

Uvar GPHData::interpolate(Time t, const Uvar& interpolation_valid_time) throw(ErrorMessage)
{

  unsigned int i=0;
  while(i<records_counted_ && t>time_(i)) i++;

  if(i==0 || i==records_counted_){
    throw ErrorMessage("GPHData::Interpolate Time out of range");
  }

  // linear interpolation
  Uvar tmid=t.et();
  Uvar t1=time_(i-1).et();
  Uvar t2=time_(i).et();
  if((t2-t1)>interpolation_valid_time) 
    ErrorMessage("GPHData:Two valid points are separated by more than valid time interval").throwMe();

  Uvar tmp1=temp_(i-1);
  Uvar tmp2=temp_(i);
  

  Uvar c1=(t2-tmid)/(t2-t1);
  Uvar c2=(tmid-t1)/(t2-t1);
  Uvar tmp=c1*tmp1+c2*tmp2;
  return(tmp);
}

Uvar GPHData::interpolate(Time t) throw(ErrorMessage)
{

  unsigned int i=0;
  while(i<records_counted_ && t>time_(i)) i++;

  if(i==0 || i==records_counted_){
    throw ErrorMessage("GPHData::Interpolate Time out of range");
  }

  // linear interpolation
  Uvar tmid=t.et();
  Uvar t1=time_(i-1).et();
  Uvar t2=time_(i).et();
  
  Uvar tmp1=temp_(i-1);
  Uvar tmp2=temp_(i);
  

  Uvar c1=(t2-tmid)/(t2-t1);
  Uvar c2=(tmid-t1)/(t2-t1);
  Uvar tmp=c1*tmp1+c2*tmp2;
  return(tmp);
}


unsigned int  GPHData::recordCount() throw(ErrorMessage)
{
  if (!records_counted_)
    {
    int position = file_.getPosition();  // remember current position
    file_.rewind();                      // position at file beginning

    bool old_flag= header_handled_;
    header_handled_=false;               // read header
    readHeader();             
    header_handled_=old_flag;           
           
    records_counted_ = 0;
    while (1)
      {  // read and count records
       string s;
       s.resize(100);
       file_.readLine(s);
       if(file_.eof()) break;
       records_counted_++;
      }
    file_.setPosition(position);         // restore original position
    }
  return(records_counted_);
}

void GPHData::read(const string& filename) 
{
  file_.reopen(filename,"r");
  recordCount();
  temp_.resize(records_counted_);
  time_.resize(records_counted_);
  read();
}

void GPHData::read() throw(ErrorMessage)
{
  if(file_read_)
    {
      throw ErrorMessage("GPHData:: Attempt to reread file");
    }
  if(records_counted_==0)
    {
      throw ErrorMessage("GPHData:: Reading to Unallocated Arrays");
    }
  readHeader();
  unsigned int i=0;
  while(1){
    string s;
    s.resize(100);
    file_.readLine(s);
    if(file_.eof()) break;
    unsigned int idx=s.find_first_of(",");
    string utc_str(s.c_str(),idx);
    float tmp=atof(&s[idx+1]);
    
    time_(i).setUtc(utc_str);
    temp_(i)=Uvar(tmp+273.15,"K");
    i++;
  }
  file_read_=true;
   
}


void GPHData::readHeader() throw(ErrorMessage,FileMgr::IoError)
{
  if(header_handled_){
    throw ErrorMessage("GPHData: Attempt to reread header");
  }
  string s;
  s.resize(100);
  file_.readLine(s);
  while(strncmp(s.c_str(),"@@GPH_DATA",10)!=0 && !file_.eof())
  {
    file_.readLine(s);
  }
  if(file_.eof()){
    throw FileMgr::IoError(file_.name(),"Invalid Header or Null Body");
  }
  header_handled_=true;
}



