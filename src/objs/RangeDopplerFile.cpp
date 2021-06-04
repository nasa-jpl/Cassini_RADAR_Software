
//-------------------------
//Class RangeDopplerFile
//---------------------------


#include <string>
#include <math.h>
#include <map>
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Array.h"
#include "Constants.h"
#include "Config.h"
#include "Plot.h"
#include "Utils.h"
#include "RangeDopplerFile.h"
using std::cout;
using std::endl;
using std::cerr;

//----------------------------------------
// Constructors: contruct with cfg only
//-------------------------------------
RangeDopplerFile::RangeDopplerFile(const string& filename, const string& mode)
  :Ndata(0),
   range("range"),
   fract_doppler("doppler"),
   Ndata_meas(0),
   range_meas("masured range"),
   fract_doppler_meas("measured fractional doppler"),
   file_(filename,mode),
   filename_(filename),  
   mode_(mode)
  {   
    sync_=2004118378;
    Nmax_ = 10000;
    range.resize(Nmax_);
    fract_doppler.resize(Nmax_);
    range_meas.resize(Nmax_);
    fract_doppler_meas.resize(Nmax_);
    if(mode=="wb" || mode=="w")
      file_.write(sync_);
    else if(mode=="rb" || mode=="r" ||mode=="r+"){
      file_.readKnown(sync_);
      mapAllRecords();
    }
    else ErrorMessage("RangeDopplerFile.cpp::Invalid file type").throwMe();
  }

RangeDopplerFile::~RangeDopplerFile()
{}

//------------------------
//load burst data
//-----------------------
void RangeDopplerFile::loadBurstData(const Time& burst_time,
				     const Uvar& active_time_offset,
				     const unsigned int& n_sab,
				     const unsigned int beam_id_input,
				     const Uvar& prf_input,
				     const Uvar& lambda_input)
  {
    //write only mode
    if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
    //(0) time
    t = burst_time;
    //(1) time offset
    time_offset = active_time_offset;
    //(2) sab number
    sab_counter= n_sab;
    //(3) beam id
    beam_id = beam_id_input;
    //(4) prf
    prf = prf_input;
    bitset(data_type,0,0,1);//burst data bit set
    //(5) lambda
    lambda= lambda_input;
  }

//------------------------
//load burst data
//-----------------------
void RangeDopplerFile::loadBurstData(const Time& burst_time,
				     const Uvar& active_time_offset,
				     const unsigned int& n_sab,
				     const unsigned int beam_id_input,
				     const Uvar& prf_input)
  {
    //write only mode
    if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
    //(0) time
    t = burst_time;
    //(1) time offset
    time_offset = active_time_offset;
    //(2) sab number
    sab_counter= n_sab;
    //(3) beam id
    beam_id = beam_id_input;
    //(4) prf
    prf = prf_input;
    bitset(data_type,0,0,1);//burst data bit set
    //(5) lambda
    lambda= speed_light/Uvar(13.78e9,"Hz");
    cout<<"transmitted wavelength is fixed : 13.78 GHz"<<endl;
  }

//-------------------
//load data
//---------------------------
void RangeDopplerFile::loadComputedRangeDoppler(const unsigned int& Npoint,
						const Uvec& range_values,
						const Dvec& fractional_doppler,
						const Uvar& bore_range_input,
						const double& fract0D,
						const int& bore_doppler_integer_input)
{  
  //write only mode
  if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
  
  //(5) number of data
  Ndata = Npoint;
  if(Ndata > Nmax_) {
    Nmax_ = Ndata;
    range.resize(Nmax_);
    fract_doppler.resize(Nmax_);
    cout<<"container size is resized"<<endl;
  }
  for(unsigned int i=0;i<Npoint;++i){
    range(i) = range_values(i);
    fract_doppler(i) = fractional_doppler(i);
  }
  bore_range = bore_range_input;
  bore_int_doppler = bore_doppler_integer_input;
  bore_fract_doppler=fract0D;
  bitset(data_type,1,1,1);//geom data bit set
}

//---------------------------------
//load measured geometry
//--------------------------------
void RangeDopplerFile::loadMeasuredRangeDoppler(const unsigned int& Npoint,
						const Uvec& range_values,
						const Dvec& fractional_doppler,
						const Uvar& bore_range_input,
						const double& fract0D,
						const int& bore_doppler_integer_input)
{  
  //write only mode
  if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
  
  //(5) number of data
  Ndata_meas = Npoint;
  if(Ndata_meas > Nmax_) {
    Nmax_ = Ndata_meas;
    range_meas.resize(Nmax_);
    fract_doppler_meas.resize(Nmax_);
    cout<<"container size is resized"<<endl;
  }
  for(unsigned int i=0;i<Npoint;++i){
    range_meas(i) = range_values(i);
    fract_doppler_meas(i) = fractional_doppler(i);
  }
  bore_range_meas = bore_range_input;
  bore_int_doppler_meas = bore_doppler_integer_input;
  bore_fract_doppler_meas=fract0D;
  bitset(data_type,2,2,1);//measurement data bit set
  bitset(data_type,3,3,1);//valid integer doppler
}

void RangeDopplerFile::loadMeasuredRangeFractDoppler( const unsigned int& Npoint,
					const Uvec& range_values,
					const Dvec& fractional_doppler,
					const Uvar& bore_range_input,
					const double& fract0D)
  {
    //write only mode
    if(mode_=="r" || mode_=="rb") ErrorMessage("RasList::works with only writing  mode").throwMe();
    
    //(5) number of data
    Ndata_meas = Npoint;
    if(Ndata_meas > Nmax_) {
      Nmax_ = Ndata_meas;
      range_meas.resize(Nmax_);
      fract_doppler_meas.resize(Nmax_);
      cout<<"container size is resized"<<endl;
    }
    for(unsigned int i=0;i<Npoint;++i){
      range_meas(i) = range_values(i);
      fract_doppler_meas(i) = fractional_doppler(i);
    }
    bore_range_meas = bore_range_input;
    bore_int_doppler_meas = 1000000;//invalid number, not calibrated yet
    bore_fract_doppler_meas=fract0D;
    bitset(data_type,2,2,1);//measurement data bit set
    bitset(data_type,3,3,0);//invalid integer doppler
  }
//-------------------------------
//Save loaded data
//------------------------------
void RangeDopplerFile::saveLoadedData()
  {
    if(bitget(data_type,0,0)==0) ErrorMessage("RangeDopplerFile.cpp:Burst data information is not avail ").throwMe();
    if( bitget(data_type,1,1)==0 &&
	bitget(data_type,2,2)==0) ErrorMessage("RangeDopplerFile.cpp:No range and doppler data is avail").throwMe();
    writeRecord();
  }


//-------------------------
//write int doppler back to file
//-------------------------
 void RangeDopplerFile::writeIntDoppler(const unsigned int& n_sab,
					const int& int_doppler)
  {
    //read record whose sab number is n_sab
    if(mode_ != "r+") ErrorMessage("RasList::can not write back into file").throwMe();
    sabnumber_vs_fileposition_ptr_= sabnumber_vs_fileposition_.find(n_sab);
    if(sabnumber_vs_fileposition_ptr_ != sabnumber_vs_fileposition_.end())
      {
	file_.setPosition(  sabnumber_vs_fileposition_ptr_ -> second);
	readRecord();
	if(bitget(data_type,3,3)==1) ErrorMessage("RangeDopplerFile:: already has a valid integer doppler value").throwMe();
	int file_position= file_.getPosition();//end position of this record
	int new_file_position= file_position - sizeof(int);
	file_.setPosition(new_file_position);
	file_.write(int_doppler);
	new_file_position=file_.getPosition();
	if(new_file_position!=file_position) ErrorMessage("RangeDoppler:: error in file position after writing int doppler").throwMe();
      }
    else{
      ErrorMessage("RangeDopplerFile: can not find the record "+toStr(n_sab)).throwMe();
    }
  }





//-----------------------------
//find and read Record(const int n_sab)
//------------------------------
bool RangeDopplerFile::findRecord(const unsigned int& n_sab)
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
//check availability of geometry and measurement data
//------------------------------------
bool RangeDopplerFile::validGeomData()
  {
    bool return_value=false;
    if(bitget(data_type,1,1)==1) return_value = true;
    return(return_value);
    
  }

bool RangeDopplerFile::validMeasData()
  {
    bool return_value=false;
    if(bitget(data_type,2,2)==1) return_value = true;
    return(return_value);
  }

bool RangeDopplerFile::validIntDoppler()
 {
   bool return_value=false;
   if(bitget(data_type,3,3)==1) return_value = true;
   return(return_value);
 }
//-----------------------------------
//find record using time
//-----------------------------------
bool RangeDopplerFile::findRecord(const Time& t)
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

//---------------------
// return sab number for nth record
//----------------------------
unsigned int& RangeDopplerFile::operator[] (const unsigned int& n)//return sab number of nth recor
 {
   if(n > sabnumber_list_.size()) ErrorMessage("DopplerRangeFile:: no record ").throwMe();
   return(sabnumber_list_[n]);
 }




//----------------------------------
//computeRecordNumbers()
//-----------------------------------
unsigned int RangeDopplerFile::computeNumberOfRecords()
  {
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    return(sabnumber_list_.size());
  }

//------------------------
//compute first and last sab number
//------------------------
 void RangeDopplerFile::getFirstLastSabNumber(unsigned int& first_sab, unsigned int& last_sab)
{
  if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
  first_sab = sabnumber_list_.front();//first element
  last_sab= sabnumber_list_.back();
}

//------------
//get beam id of sab
//-------------
unsigned int RangeDopplerFile::getBeamId(const unsigned int& sab_input)
  { 
    sabnumber_vs_beamid_ptr_= sabnumber_vs_beamid_.find(sab_input);
    if(sabnumber_vs_beamid_ptr_ == sabnumber_vs_beamid_.end()){
      ErrorMessage("RangeDopplerFile:: no such record "+toStr(sab_input)).throwMe();
    }
    return(sabnumber_vs_beamid_ptr_->second);
  }
		
//------------------
//display range doppler
//------------------
void RangeDopplerFile::displayRangeDoppler(const unsigned int& n_sab)
  {  
    if(mode_=="w" || mode_=="wb") ErrorMessage("RasList::works with only reading mode").throwMe();
    if(!findRecord(n_sab)) ErrorMessage("RangeDopplerFile.cpp:: no record ").throwMe();
    Uvec x("",Ndata);
    Uvec y("",Ndata);
    for(unsigned int i=0;i<Ndata;++i){
      x(i) = range(i);
      y(i)=fract_doppler(i);
    }
    Plot a;
    a.addXY(x,"km",y,"KHz",line("solid","black",1),sym("none"));
    a.setTitle("doppler vs range");
    a.setXlabelandUnits("range");
    a.setYlabelandUnits("doppler");
    a.show("x");
  }

//------------------------
//close the file
//--------------------------	
void RangeDopplerFile::close()
{
  file_.close();
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
void RangeDopplerFile::mapAllRecords()
  {
    sabnumber_list_.clear();
    sabnumber_vs_fileposition_.clear();
    time_list_.clear();
    time_vs_sabnumber_.clear();
    cout<<"start Mapping "<<endl;
    //extract file pointer, sab number and time information from all the records
    if(mode_=="wb"||mode_=="w") ErrorMessage("file mode is not for reading").throwMe();
    int dumm;
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
	//(5) sab number vs beam number
	sabnumber_vs_beamid_[sab_counter]=beam_id;
      }
    cout<<"All the maps are complete "<<endl;
  }


//-----------------------------------
//read next record
//-----------------------------------
void RangeDopplerFile::readRecord()
  {
    if(file_.eof()) ErrorMessage("RangeDopplerFile::readRecord(): no more record").throwMe();
    if(mode_=="wb"||mode_=="w") ErrorMessage("RangeDopplerFile::file mode is not for reading").throwMe();
       
    //(0) time
    double t_dumm;
    file_.read(t_dumm);
    t.setEt(Uvar(t_dumm,"s"));

    //(1) time offset
    file_.read(t_dumm);
    time_offset =Uvar(t_dumm,"s");

    //(2) data type
    file_.read(data_type);
    if(bitget(data_type,1,1)==0 &&
       bitget(data_type,2,2)==0) 
      ErrorMessage("RangeDopplerFile::readRecord: no avail data").throwMe();
   
    //(3)sab counter
    file_.read(sab_counter); 

    //(4) beam number
    file_.read(beam_id);//unsigned int
    
    double x_dumm;
    //(5) prf
    file_.read(x_dumm);
    prf = Uvar(x_dumm,"Hz");

    //(5)-1 lambda
    file_.read(x_dumm);
    lambda=Uvar(x_dumm,"km");

    if(bitget(data_type,1,1)==1){
      //(5) valid data number per burst
      file_.read(Ndata);//unsigned int
      
      //(6) read ata
      for(unsigned int i=0;i<Ndata;++i){
	file_.read(x_dumm);
	range(i)=Uvar(x_dumm,"km");
	file_.read(x_dumm);
	fract_doppler(i)=x_dumm;
      }
      file_.read(x_dumm);
      bore_range = Uvar(x_dumm,"km");
      file_.read(bore_fract_doppler);
      file_.read(bore_int_doppler);
    }
    if(bitget(data_type,2,2)==1){
      //(7) read number of data points
      file_.read(Ndata_meas);
      //(8)
      for(unsigned int i=0;i<Ndata_meas;++i){
	file_.read(x_dumm);
	range_meas(i)=Uvar(x_dumm,"km");
	file_.read(x_dumm);
	fract_doppler_meas(i)=x_dumm;
      }
      file_.read(x_dumm);
      bore_range_meas=Uvar(x_dumm,"km");
      file_.read(bore_fract_doppler_meas);
      file_.read(bore_int_doppler_meas);
    }
  }

//------------------------
//{write) one record
//------------------------
void RangeDopplerFile::writeRecord()
  { 
    if(mode_=="rb"||mode_=="r") ErrorMessage("file mode is not for writing").throwMe();
    
    //(0) time
    file_.write(t.et().getInUnits("s"));//double

    //(1) offset
    file_.write(time_offset.getInUnits("s"));//double


    //(2) data type
    file_.write(data_type);
   
    //(3)sab counter
    file_.write(sab_counter); 

    //(4) beam number
    file_.write(beam_id);//unsigned int
    
    double x_dumm;
    //(5) prf
    x_dumm = prf.getInUnits("Hz");
    file_.write(x_dumm);
   
    //(5)-1 lambda
    x_dumm=lambda.getInUnits("km");
    file_.write(x_dumm);

    if(bitget(data_type,1,1)==1){
      //(5) valid data number per burst
      file_.write(Ndata);//unsigned int
      
      //(6) write data
      for(unsigned int i=0;i<Ndata;++i){
	x_dumm = range(i).getInUnits("km");
	file_.write(x_dumm);
       	file_.write(fract_doppler(i));
      }
      x_dumm = bore_range.getInUnits("km");
      file_.write(x_dumm);
      file_.write(bore_fract_doppler);
      file_.write(bore_int_doppler);
    }
    if(bitget(data_type,2,2)==1){
      //(7) write number of data points
      file_.write(Ndata_meas);
      //(8)
      for(unsigned int i=0;i<Ndata_meas;++i){
	x_dumm = range_meas(i).getInUnits("km");
	file_.write(x_dumm);
	file_.write(fract_doppler_meas(i));
      }
      x_dumm = bore_range_meas.getInUnits("km");
      file_.write(x_dumm);
      file_.write(bore_fract_doppler_meas);
      file_.write(bore_int_doppler_meas);
    }
    data_type=0;
  }


   
