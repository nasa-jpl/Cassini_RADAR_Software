static const char rcs_id_l1b_c[] =
  "@(#) $Id: L1B.cpp,v 11.6 2012/09/27 20:01:58 richw Exp $";


//--------------
//Code update/change log
// Y. Gim on March 17, 2005
// Adjust rip, cip, hip based on config keyword delta tau
// parameters
//------------

#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include "L1B.h"
#include "Error.h"
#include "config_keywords.h"
#include <string>
#include "L1I.h"
#include "PDSLabel.h"
#include "Constants.h"
#include "IebProfile.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::setfill;
using std::setprecision;
using std::setw;
using std::stringstream;

//-------------------------
// Static initializations
//-------------------------

string L1B::bodp_type_name_[BODPTYPECOUNT] =
{
  "SBDR",
  "LBDR",
  "ABDR",
  ""
};
int L1B::RADIOMETER_MODE    = 1 << 0;
int L1B::SCATTEROMETER_MODE = 1 << 1;
int L1B::ALTIMETER_MODE     = 1 << 2;
int L1B::SAR_MODE           = 1 << 3;

string L1B::pds_version_id_ = "PDS3";
string L1B::record_type_ = "FIXED_LENGTH";

string L1B::data_set_version_id_ = "V1.0";
string L1B::data_set_id_[BODPTYPECOUNT] =
{
  "CO-V/E/J/S-RADAR-3-SBDR-" + data_set_version_id_,
  "CO-V/E/J/S-RADAR-3-LBDR-" + data_set_version_id_,
  "CO-SSA-RADAR-3-ABDR-" + data_set_version_id_,
  ""
};
string L1B::data_set_name_[BODPTYPECOUNT] =
{
  "CASSINI RADAR SHORT BURST DATA RECORD",
  "CASSINI RADAR LONG BURST DATA RECORD",
  "CASSINI RADAR ALTIMETER BURST DATA RECORD",
  ""
};
string L1B::producer_id_ = "JPL";
string L1B::producer_institution_name_ = "JET PROPULSION LABORATORY";
string L1B::instrument_host_name_ = "CASSINI ORBITER";
string L1B::instrument_host_id_ = "CO";
string L1B::instrument_name_ = "CASSINI RADAR";
string L1B::instrument_id_ = "RADAR";
string L1B::mission_name_ = "CASSINI-HUYGENS";
string L1B::description_prefix_[BODPTYPECOUNT] =
{
"CASSINI RADAR SHORT BURST DATA RECORD FOR THE ",
"CASSINI RADAR LONG BURST DATA RECORD FOR THE ",
"CASSINI RADAR ALTIMETER BURST DATA RECORD FOR THE ",
""
};
string L1B::processing_history_text_[BODPTYPECOUNT] =
{
"NONE",
"NONE",
"NONE",
""
};
string L1B::interchange_format_ = "BINARY";
string L1B::table_description_[BODPTYPECOUNT] =
{
"This is the table definition for a Cassini Radar Short Burst Data Record, \
which contains engineering telemetry, spacecraft geometry, and calibrated \
science data for each burst in the pass.",
"This is the table definition for a Cassini Radar Long Burst Data Record, \
which includes a Short Burst Data Record (engineering telemetry, \
spacecraft geometry, and calibrated science data) plus the raw counts of \
the sampled echo data.",
"This is the table definition for a Cassini Radar Altimeter Burst Data Record, \
which includes a Short Burst Data Record (engineering telemetry, \
spacecraft geometry, and calibrated science data) plus the altimeter profile.",
""
};

int L1B::table_column_count_[BODPTYPECOUNT] = { 0 };

//-----------------
// Methods for L1B
//-----------------

//--------------
// Constructors
//--------------



L1B::L1B(const std::string& filename, const std::string& mode,const std::string& filetype, bool abdr_flag)
  : BurstData(filename,mode,filetype), radar_data("radar_data",1),
    abdr_flag_(abdr_flag),data_complete_(false), map_complete_(false),
    pds_radar_mode_(0),
    check_pds_string_lengths_(true), config_set_(false)
  {

  header_size_=L1B_HEADER_LENGTH;
  if (filetype =="active")
    {
      Nradar_data_max_=32*1024;
      record_size_=L1B_PASSIVE_RECORD_LENGTH+ Nradar_data_max_*sizeof(float);
      parameters_.appendParameter("radar_data",Nradar_data_max_*sizeof(float),
			      Parameter::FDATA,(void*) &radar_data);
    }
  else if (filetype =="passive"){
    Nradar_data_max_=1;
    record_size_=L1B_PASSIVE_RECORD_LENGTH;
  }
  else
    { 
     ErrorMessage e("no proper filetype: either active or passive"); 
     e.throwMe(); 
    }
    
  radar_data.resize(Nradar_data_max_);
  for (unsigned int i = 0 ; i < Nradar_data_max_; i++) radar_data(i)=0.0;

  }

L1B::~L1B()
{
//  if(mode_[0]=='w') rewriteHeader();
}


//--------------
// I/O
//--------------

void L1B::config(Config& cfg)
{
  // Read configuration file parameters that are used to set values
  // in the PDS labels.
  producer_full_name_ =
    PDSLabel::replaceSpaces(cfg.str("PRODUCER_FULL_NAME"));
  product_version_ = cfg.getInt("BODP_PRODUCT_VERSION_ID");
  software_version_ = cfg.str("BODP_SOFTWARE_VERSION_ID");

  if (cfg.keywordExists("IEB_Trigger") && cfg.keywordExists("IDAPT_delta") &&
      cfg.keywordExists("CDS_Cmd_Delay"))
  {
    // Work with a copy of the original Config object
    Config cfg_tmp(cfg.filename());
    IebProfile iebprofile(cfg_tmp);
    pds_trigger_time_ = iebprofile.getTriggerTime();
  }
  else
  {
    cerr << "[L1B::config] Warning:  trigger time calculation " <<
      "in PDS label will be" << endl;
    cerr << "skipped because at least one of the following parameters " <<
      "is missing from " << endl;
    cerr << cfg.filename() << ":" << endl;
    cerr << "    IEB_Trigger IDAPT_delta CDS_Cmd_Delay" << endl;
  }
  flyby_id_ = cfg.str("FLYBY_ID");

  if (cfg.keywordExists("CHECK_PDS_STRING_LENGTHS"))
  {
    check_pds_string_lengths_ = (bool) cfg.getInt("CHECK_PDS_STRING_LENGTHS");
  }
  if (!check_pds_string_lengths_)
  {
    // Display a warning message.  Normally, this option should be enabled.
    cerr << "Warning:  PDS string length limit checks have been disabled." <<
      endl;
    cerr << "To enable the checks, set CHECK_PDS_STRING_LENGTHS to 1 in the " <<
      "config file." << endl;
  }

  // Set other values used in the PDS label.
  data_take_id_ = BurstData::data_take_id_;

  config_set_ = true;
}

//------------------------------------------------------------------
// writeHeader()
//
// Write a PDS label to the current L1B file.
//------------------------------------------------------------------

void L1B::writeHeader()
{
  string label;
  char headbuf[100];

  if (mode_ == "r" || mode_ == "rb")
  {
    L1BError e("Can't write to input file " + filename_, L1B::write_error);
    e.throwMe();
  }

  // Should not call this routine twice
  if (header_handled_)
  {
    L1BError e("Attempt to write header twice", L1B::write_error);
    e.throwMe();
  }

  header_handled_ = true;

//  For reference, the following code segment defines the format of the
//  old binary label for L1B output products.
//
//  string header;
//  char headbuf[100];
//  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",0,record_size_);
//  header=headbuf;
//  if ( header.length() != header_size_){
//    char msg[100];
//    sprintf(msg,"Incorrect header size (%d bytes) should be %d",
//	    header.length(),header_size_);
//    L1BError e(msg,L1B::write_error);
//    e.throwMe();
//  }
//  file_.write(header);

  if (ft_ == active)
  {
    if(abdr_flag_){
      label = constructLabel(ABDR);
    }
    else{
      label = constructLabel(LBDR);
    }
  }
  else if (ft_ == passive)
  {
    label = constructLabel(SBDR);
  }
  else
  {
    sprintf(headbuf, "L1B::writeHeader:  unknown output file type %d\n", ft_);
    L1BError e(headbuf, L1B::internal_error);
    e.throwMe();
  }
  file_.write(label);
  header_handled_ = true;
}

void L1B::rewriteHeader()
{
  string label;
  int position = file_.getPosition();
  char headbuf[100];

  file_.rewind();
  if (ft_ == active)
  {
    if(abdr_flag_){
      label = constructLabel(ABDR);
    }
    else{
      label = constructLabel(LBDR);
    }
  }
  else if (ft_ == passive)
  {
    label = constructLabel(SBDR);
  }
  else
  {
    sprintf(headbuf, "L1B::rewriteHeader:  unknown output file type %d\n", ft_);
    L1BError e(headbuf, L1B::internal_error);
    e.throwMe();
  }
  file_.write(label);
  header_handled_ = true;
  file_.setPosition(position);
}

//------------------------------------------------------------------
// readHeader()
//
// Read a L1B header in PDS or old format from the current L1B file.
// This method defines the L1B Header format.
//------------------------------------------------------------------

void L1B::readHeader()
{
  if (mode_ == "w" || mode_ == "wb")
    {
      L1BError e("Can't read from output file " + filename_,
	       L1B::read_error);
      e.throwMe();
    }
  if (header_handled_)
    {
      L1BError e("Attempt to read header twice",
	       L1B::read_error);
      e.throwMe();
    }
  PDSLabel label(file_);
  if (label.label() == EMPTY_STRING)
  {
    // File does not appear to contain a PDS label.
    // Assume that we have an old-style (40-byte) header.
    cerr << "[L1B::readHeader] Warning:  PDS label not found in " <<
      file_.name() << endl;
    cerr << "Assuming 40-byte header." << endl;

    checkEndedness();
    string header;
    header.resize(L1B_HEADER_LENGTH);
    file_.read(header);
    string token=get_next_token(header," /n/t"); // skip text
    token=get_next_token(header," /n/t"); // skip text
    
    token=get_next_token(header," /n/t"); // get num records
    record_count_=atoi(token.c_str());
    records_counted_=true;
  
    token=get_next_token(header," /n/t"); // skip text
    token=get_next_token(header," /n/t"); // skip text
  
    token=get_next_token(header," /n/t"); // check record size
    unsigned int recsize=(unsigned int) atoi(token.c_str());
    if(recsize!=record_size_){
      L1BError e("Record size incorrect in " + filename_, L1B::read_error);
      e.throwMe();
    }
  }
  else
  {
    // PDS label found.  Update the same private members as is done
    // for the 40-byte headers as well as header_size_, which is set
    // in the constructor.
    record_size_ = label.recordLength();
    header_size_ = label.labelLength();
    int file_records = label.fileRecords();
    int label_records = label.labelRecords();
    record_count_ = file_records - label_records;
    records_counted_ = true;

    // Set additional private members.
    if (!label.getString("PRODUCT_ID", &product_id_))
    {
      L1BError e("PRODUCT_ID missing from " + filename_, L1B::unknown_parameter);
      e.throwMe();
    }

    // header_size_ must be set (to the correct value) before
    // checkEndedness() is called.
    checkEndedness();
    file_.setPosition(header_size_);
  }
  header_handled_ = true;
}

//------------------------------------------------------------------
// writeRecord()
//
// Write the current L1B record to the end of the current L1B file in binary.
// This method defines the L1B record format.
// It is the responsibility of the user to ensure that valid data is present
// before calling this routine.  No tracking of data status is performed.
//------------------------------------------------------------------

void L1B::writeRecord(int abs_record_no) 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1BError e("Can't write to input file " + filename_, L1B::write_error);
      e.throwMe();
    }
  if (!header_handled_)
    {
      L1BError e("Can't write L1B record until header is written",
		 L1B::write_error);
      e.throwMe();
    }
  if (!data_complete_)
    {
      L1BError e("Can't write L1B record until data is complete",
		 L1B::write_error);
      e.throwMe();
    }


  
  if(file_.getPosition() > MAX_FILE_LENGTH-(int)record_size_) 
    createNextFile();


  int start_position=file_.getPosition();

  //passive
  writePassiveSABData(abs_record_no); 
  writeGeometry(); 
  writeMeasurementGeometry();
  if(ft_==active){
    for (unsigned int i = 0; i < Nradar_data; i++){
      file_.write(radar_data(i));
    }

    if(baq_mode==3){
      //write dc offset of compressed scatt mode
      file_.write(radar_data(Nradar_data));
      for (unsigned int i = Nradar_data+1; i< Nradar_data_max_; i++){
	file_.write(float(0.0));
      }
    }
    else{
      for (unsigned int i = Nradar_data; i< Nradar_data_max_; i++){
	file_.write(float(0.0));
      }
    }
  }
  
  int end_position=file_.getPosition();
  int len_written=end_position-start_position;
  if((unsigned int)len_written!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_written,record_size_);
    L1BError(msg,L1B::write_error).throwMe();
  }
  record_count_++;

  }


//------------------------------------------------------------------
// readRecord()
//
// Read a L1B record from the current L1B file.
//------------------------------------------------------------------

void L1B::readRecord() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1BError e("Can't read from output file " + filename_,
	       L1B::read_error);
//
// create next file
      e.throwMe();
    }  if (!header_handled_)
    {
      L1BError e("Can't read L1B record until header is read",
	       L1B::read_error);
      e.throwMe();
    } 

  // This check has the side effect of moving to the next  file if necessary
  if(eof())
    {
      L1BError e("Unexpected EOF in L1B file "+filename_,
	       L1B::read_error);
      e.throwMe();
    } 


  // Read in all public variables
  //convert this part to read
  //header slow fast data footer
 
  int start_position=file_.getPosition();
  
  readPassiveSABData();  
  readGeometry(); 
  readMeasurementGeometry();
  if (ft_==active)
  {     
    for (unsigned int i = 0; i <Nradar_data; i++){file_.read(radar_data(i));}
    if(baq_mode==3) file_.read(radar_data(Nradar_data));
    int position=file_.getPosition();
    
    if(baq_mode==3)
      position+=sizeof(float)*(Nradar_data_max_-Nradar_data-1);
    else
      position+=sizeof(float)*(Nradar_data_max_-Nradar_data);
    file_.setPosition(position);
  }
  
  
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,record_size_);
    L1BError e(msg,L1B::read_error);
    e.throwMe();
  }
  data_complete_=true;
  data_read_=true;  
  }


//-------------------------------------------
// Routine for computing geometry
// for a previously extant L1B file
// does not read echo data buffer
//-------------------------------------------
// obsolete
/*****
void
L1B::locate(L1B& l1b){
  resetRecord();
  try{
  copyPassive(l1b);
  computeTime();
  computeGeometry();
  setQualityFlag();
  }
  catch(GeomError e){
    if(e.error_type==GeomError::spice_bad_ckernel)
      {
	cerr<<"L1B.cpp::locate(L1B& l1b): bad ckernel "<<endl;
	quality_flag+=1;
      }
    else 
      {
	cerr<<"L1B.cpp::locate(L1B& l1b): misc spice error: possibly no spice data avail "<<endl;
	quality_flag+=2;
      }
  }
 
 data_complete_=true;
}
***/

// Routine for computing geometry in current file
// SAB data must have already been loaded

void
L1B::locate(const Umat& azim_1way3dB_ellipse_fit, 
	    const Umat& elev_1way3dB_ellipse_fit,
	    const Umat& azim_2way3dB_ellipse_fit, 
	    const Umat& elev_2way3dB_ellipse_fit){
  science_qual_flag=0;
  data_complete_=false;
  if(!data_read_){
    L1BError e("locate(): requires SAB data to have been loaded",file_not_read);
    e.throwMe();
  }
  try{
  computeTime();
  computeGeometry();
  setQualityFlag();
  computeMeasurementGeometry(azim_1way3dB_ellipse_fit,
			     elev_1way3dB_ellipse_fit,
			     azim_2way3dB_ellipse_fit,
			     elev_2way3dB_ellipse_fit);
  }
  catch(GeomError e){
    if(e.error_type==GeomError::spice_bad_ckernel)
      {
	cerr<<"L1B.cpp::locate(): spice bad ckernel "<<endl;
	quality_flag+=1;
      }
    else
      {
	cerr<<"L1B.cpp::locate(): misc spice error: possibly no spice data avail "<<endl;
	quality_flag+=2;
      }
  }
  
 
 data_complete_=true;
}

//-------------------------------------------
// Routine for converting L1B file with modified time:
// used by PointTargetSim::run
// reads L1B record directly from file,
// and then computes new geometry
// sets burst start time to parameter t0 instead of
// computing from L1B file values.
// Does not read echo buffer
// Rereads from beginning if L1B EOF is reached
//-------------------------------------------
// Obsolete
/****
void
L1B::locate(L1B& l1b, const Time& t0){
  resetRecord();
  try{
    copyPassive(l1b);
  }
  catch(L1B::L1BError e){
    if(l1b.eof()){
      
      // Go to first record of current file
      // If data take includes multiple files, this
      // throws an error. The current
      // routine is primarily intended for use by the
      // Point Target Simulator, In which case 1 L1A file
      // should be sufficient to describe desired RADAR behavior
      l1b.gotoFirstRecord();
      
      copyPassive(l1b);
    }
    else e.throwMe();
  }
  try{
    editTime(t0);
    setQualityFlag();
  }
  catch(GeomError e){
    if(e.error_type==GeomError::spice_bad_ckernel)
      {
	cerr<<"L1B.cpp::locate(L1B& l1b, Time t): bad ckernel "<<endl;
	quality_flag+=1;
      }
    else
      {
	cerr<<"L1B.cpp::locate(L1B& l1b, Time t): misc spice error: possibly no spice data avail "<<endl;
      quality_flag+=2;
      }
  }
  data_complete_=true;
}

***/

void L1B::createRecord(const Time& t0){
  resetRecord();
  pretendDataRead();

  try{
    editTime(t0);
  }
  catch(GeomError e){
    if(e.error_type==GeomError::spice_bad_ckernel) 
      {
	cerr<<"L1B.cpp::createRecord(Time t) spice kernel is bad "<<endl;
	quality_flag+=1;
      }
    else
      {
	cerr<<"L1B.cpp::createRecord(Time t): misc spice error: possibly no spice data avail "<<endl;
	quality_flag+=2;
      }
  }
  data_complete_=true;
}



void L1B::resetRecord(){
  quality_flag=0;
  data_read_=false;
  data_complete_=false;
}

void
L1B::copy(L1B& l1b){
  resetRecord();
  int start_position=BurstData::copy(l1b);
  if(l1b.Nradar_data_max_ != Nradar_data_max_){
   L1BError e("Cannot copy L1B, Ndata_max incompatible");
   e.throwMe();
 } 
   if (l1b.ft_==active)
  {     
    for (unsigned int i = 0; i <Nradar_data; i++){
      l1b.file_.read(radar_data(i));
    }
    int position=l1b.file_.getPosition();
    position+=sizeof(float)*(Nradar_data_max_-Nradar_data);
    l1b.file_.setPosition(position);
  }
  
  
 
  int end_position=l1b.file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=l1b.record_size_){
     char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,l1b.record_size_);
    L1BError e(msg,L1B::read_error);
    e.throwMe();
  }

  data_complete_=true;
  data_read_=true;
}

// Routine used to output SAR ancillaty data to LBDR file by sar_proc
// Assumes we start at the end of the record (beginning of next record)
void
L1B::outputSARAncillaryData(const L1I& l1i){


  int start_pos = file_.getPosition();
  if((start_pos-header_size_)%record_size_!=0 ){
    ErrorMessage e("L1B::outputSARAncillaryData Not at end of Record");
    e.throwMe();
  }

  int pos=start_pos-record_size_;
  file_.setPosition(pos);
  writeParameter("science_qual_flag",l1i.science_qual_flag);
  writeParameter("sar_azimuth_res",l1i.sar_azimuth_res,"km");
  writeParameter("sar_range_res",l1i.sar_range_res,"km");
  writeParameter("sar_centroid_bidr_lat",l1i.sar_centroid_bidr_lat,"deg");
  writeParameter("sar_centroid_bidr_lon",l1i.sar_centroid_bidr_lon,"deg");

  if(l1i.output_scatterometer_info){
    writeParameter("sigma0_uncorrected",l1i.sigma0_uncorrected,"");
    writeParameter("sigma0_corrected",l1i.sigma0_corrected,"");
    writeParameter("sigma0_uncorrected_std",l1i.sigma0_uncorrected_std,"");
    writeParameter("x_factor",l1i.x_factor,"J");
    writeParameter("total_echo_energy",l1i.total_echo_energy,"J");
    writeParameter("noise_echo_energy",l1i.noise_echo_energy,"J");
    writeParameter("num_pulses_received",l1i.num_pulses_received);
  }
  file_.setPosition(start_pos);
}

void
L1B::copyTo(L1I& l1i){
  l1i.resetRecord();
  unsigned int dummy;
  copyTo(l1i,l1i.radar_data,dummy);
  l1i.burst_no=l1i.record_id%1000000;
}

// for backwards compatibility with L1IOld only
void
L1B::copyTo(BurstData& bd, fdata&r, unsigned int& N){
  int start_position=BurstData::copyTo(bd);
  N=bd.Nradar_data;
  if(ft_==active)
  {     
      for (unsigned int i = 0; i <N; i++){
	file_.read(r(i));
      }
      int position=file_.getPosition();
      position+=sizeof(float)*(Nradar_data_max_-N);
      file_.setPosition(position);
  }
   
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,record_size_);
    L1BError e(msg,L1B::read_error);
    e.throwMe();
  }
}

void
L1B::copyTo(BurstData& bd, float* r, unsigned int& N){
  int start_position=BurstData::copyTo(bd);
  N=bd.Nradar_data;
  if(ft_==active)
  {     
      for (unsigned int i = 0; i <N; i++){
	file_.read(r[i]);
      }
      int position=file_.getPosition();
      position+=sizeof(float)*(Nradar_data_max_-N);
      file_.setPosition(position);
  }
   
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,record_size_);
    L1BError e(msg,L1B::read_error);
    e.throwMe();
  }
}

void 
L1B::checkNradar_data()
{
  // throws Error Message if Nradar_data is too large for buffer or 0
  // This will fail if called for passive data
  if(Nradar_data<1 || Nradar_data>Nradar_data_max_){
    ErrorMessage e("L1B::checkNradar_data failed for Nradar_data="+toStr(Nradar_data));
    e.throwMe();
  }
}

void
L1B::copyPassive(L1B& l1b){
  resetRecord();
  int start_position=BurstData::copy(l1b);
  if (l1b.ft_==active)
  {   
      int position=l1b.file_.getPosition();
      position+=sizeof(float)*(l1b.Nradar_data_max_);
      l1b.file_.setPosition(position);
  }
  
  
 
  int end_position=l1b.file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=l1b.record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,l1b.record_size_);
    L1BError e(msg,L1B::read_error);
    e.throwMe();
  }
  data_complete_=true;
  data_read_=true;
}

//=====================================================================
// constructLabel()
//
// This method defines the actual keywords and values that appear in
// the PDS labels of the SBDR and LBDR products.  The keywords and values
// are specific to these two products.  The string returned by this method
// contains a properly formatted PDS label for the specified type of
// data and is to be prepended to the output data.
//=====================================================================
string L1B::constructLabel(const BODPTypeE bodp_type)
{
  if (!config_set_)
  {
    ErrorMessage e("L1B::constructLabel:  configuration file has not been read");
    e.throwMe();
  }

  if (check_pds_string_lengths_)
  { 
    // Populate the table that holds the limits on the lengths of PDS
    // character values.  (The actual checking is done within the
    // constructor for KeyValPair().)
    PDSLabel::setPDSCharValLimits();
  }

  //
  // Construct list of keyword-value pairs for BODP PDS label
  //
  KeyValSeq list;
  int data_records = record_count_;
  int record_length = record_size_; // in bytes
  int sclk;
  string target_name = BurstData::target_name_;
  string product_creation_time = Time::getCurrentUtcDOYTime();
  string table_name = bodp_type_name_[bodp_type] + "_TABLE";
  string description = description_prefix_[bodp_type];
  string flyby_desc = PDSLabel::toUpper(flyby_id_);
  stringstream buf(ios::out);
  Time epoch_time = BurstData::epoch_;
  Time closest_approach_time = BurstData::closest_approach_time_;

  description += PDSLabel::replaceSpaces(flyby_desc);
  description += " PASS WITH CLOSEST APPROACH TIME ";
  if (! closest_approach_time.valid())
  {
    description += TIMESTAMP_PLACEHOLDER;
  }
  else
  {
    description += closest_approach_time.utc("ISOD");
  }
  description += " TRIGGER TIME ";
  if (! pds_trigger_time_.valid())
  {
    description += TIMESTAMP_PLACEHOLDER;
  }
  else
  {
    description += pds_trigger_time_.utc("ISOD");
  }
  description += " EPOCH TIME ";
  if (! epoch_time.valid())
  {
    description += TIMESTAMP_PLACEHOLDER;
  }
  else
  {
    description += epoch_time.utc("ISOD");
  }
  description += ".";

  table_column_count_[SBDR] = 235;
  table_column_count_[LBDR] = table_column_count_[SBDR] + 1;
  table_column_count_[ABDR] = table_column_count_[SBDR] + 1;

  //
  // Create the list of keyword-value pairs in the order in which they
  // are to be written to the label.
  //

  list.push_back(KeyValPair("PDS_VERSION_ID", pds_version_id_, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("/*       FILE FORMAT AND LENGTH */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("RECORD_TYPE", record_type_, !IS_ODL_STRING));
  list.push_back(KeyValPair("RECORD_BYTES", record_length));
  //
  // Initialize the next two values to their smallest possible values.
  //
  list.push_back(KeyValPair("FILE_RECORDS", data_records+1));
  list.push_back(KeyValPair("LABEL_RECORDS", 1));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair(
    "/*       POINTERS TO START RECORDS OF OBJECTS IN FILE */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  //
  // Initialize the next value to its smallest possible value.
  //
  list.push_back(KeyValPair("^" + table_name, 2));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair("/*       FILE DESCRIPTION */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("DATA_SET_ID", data_set_id_[bodp_type],
    IS_ODL_STRING));
  list.push_back(KeyValPair("DATA_SET_NAME", data_set_name_[bodp_type],
    IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_INSTITUTION_NAME",
    producer_institution_name_, IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_ID", producer_id_, !IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCER_FULL_NAME", producer_full_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCT_ID", createProductID(bodp_type),
    !IS_ODL_STRING));
  list.push_back(KeyValPair("PRODUCT_VERSION_ID", product_version_,
    2, NO_UNITS));
  list.push_back(KeyValPair("INSTRUMENT_HOST_NAME", instrument_host_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_HOST_ID", instrument_host_id_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_NAME", instrument_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("INSTRUMENT_ID", instrument_id_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("TARGET_NAME", PDSLabel::toUpper(target_name),
    !IS_ODL_STRING));
  if (! sclk_start_.valid())
  {
    list.push_back(KeyValPair("START_TIME", TIMESTAMP_PLACEHOLDER,
      !IS_ODL_STRING));
  }
  else
  {
    list.push_back(KeyValPair("START_TIME", sclk_start_.utc("ISOD"),
      !IS_ODL_STRING));
  }
  if (! sclk_stop_.valid())
  {
    list.push_back(KeyValPair("STOP_TIME", TIMESTAMP_PLACEHOLDER,
      !IS_ODL_STRING));
  }
  else
  {
    list.push_back(KeyValPair("STOP_TIME", sclk_stop_.utc("ISOD"),
      !IS_ODL_STRING));
  }
  //
  // Since PDS defines SPACECRAFT_CLOCK_START_COUNT and
  // SPACECRAFT_CLOCK_STOP_COUNT as character-valued, format both values as
  // strings so that their lengths can be limit-checked.
  //
  if (! sclk_start_.valid())
  {
    sclk = 0; // placeholder
  }
  else
  {
    sclk = sclk_start_.sclk(cassini_str);
  }
  buf << std::setw(9) << std::setfill('0') << sclk;
  list.push_back(KeyValPair("SPACECRAFT_CLOCK_START_COUNT", buf.str(),
    !IS_ODL_STRING));
  buf.seekp(0, ios::beg);
  if (! sclk_stop_.valid())
  {
    sclk = 0; // placeholder
  }
  else
  {
    sclk = sclk_stop_.sclk(cassini_str);
  }
  buf << std::setw(9) << std::setfill('0') << sclk;
  list.push_back(KeyValPair("SPACECRAFT_CLOCK_STOP_COUNT", buf.str(),
    !IS_ODL_STRING));
  buf.seekp(0, ios::beg);
  list.push_back(KeyValPair("PRODUCT_CREATION_TIME", product_creation_time,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("MISSION_NAME", mission_name_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("SOFTWARE_VERSION_ID", software_version_,
    IS_ODL_STRING));
  list.push_back(KeyValPair("DESCRIPTION", description, IS_ODL_STRING));
  list.push_back(KeyValPair("PROCESSING_HISTORY_TEXT",
    processing_history_text_[bodp_type], IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));

  list.push_back(KeyValPair(
    "/*       DESCRIPTION OF OBJECTS CONTAINED IN FILE */",
    EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair(EMPTY_STRING, EMPTY_STRING, !IS_ODL_STRING));
  list.push_back(KeyValPair("OBJECT", table_name, !IS_ODL_STRING));
  list.push_back(KeyValPair("INTERCHANGE_FORMAT", interchange_format_,
    !IS_ODL_STRING));
  list.push_back(KeyValPair("ROWS", data_records));
  list.push_back(KeyValPair("COLUMNS", table_column_count_[bodp_type]));
  list.push_back(KeyValPair("ROW_BYTES", record_length));
  list.push_back(KeyValPair("^STRUCTURE", bodp_type_name_[bodp_type] + ".FMT",
    IS_ODL_STRING));
  list.push_back(KeyValPair("DESCRIPTION", table_description_[bodp_type],
    IS_ODL_STRING));
  list.push_back(KeyValPair("END_OBJECT", table_name, !IS_ODL_STRING));

  list.push_back(KeyValPair("END", EMPTY_STRING, !IS_ODL_STRING));

  //
  // Create the PDS label in the proper format.
  //
  PDSLabel label(list, record_length, data_records, "^" + table_name);
  return(label.label());
}

//=====================================================================
// createProductID()
//
// Determine the value of the PRODUCT_ID keyword in the BODP PDS label.
// This value is dependent upon the BODP data type, the radar mode,
// the radar observation counter (data take ID), and the product version.
//=====================================================================
string L1B::createProductID(const BODPTypeE bodp_type)
{
  string product_id = EMPTY_STRING;
  stringstream buf_stream(ios::app|ios::out);

  switch (bodp_type)
  {
    case SBDR:
    case LBDR:
    case ABDR:
      product_id += bodp_type_name_[bodp_type];
      break;
    default:
      ErrorMessage e("L1B::createProductID:  unknown BODP type " +
        toStr(bodp_type));
      e.throwMe();
      break;
  }

  // Radar mode
  buf_stream << "_" << std::setw(2) << std::setfill('0') << pds_radar_mode_;

  // Radar observation counter
  buf_stream << "_D" << std::setw(3) << std::setfill('0') << data_take_id_;

  // Version number
  buf_stream << "_V" << std::setw(2) << std::setfill('0') << product_version_;
  product_id += buf_stream.str();
  return product_id;
}

// Parameter access methods

double L1B::getSpecialParam(const char* name){
  double retval;
  
  retval=BurstData::getSpecialParam(name);
  return(retval);
}
void L1B::enableSpecialParam(const char* name)
{
  BurstData::enableSpecialParam(name);
}
void L1B::disableSpecialParam(const char* name)
{
  BurstData::disableSpecialParam(name);
}

//=====================================================================
// productID()
//
// Return the value of the private member product_id_.  Its value should
// have been set to the value of the PRODUCT_ID keyword in the BODP PDS
// label by a call to L1B::readHeader().
//=====================================================================
string L1B::productID()
{
  return product_id_;
}
 
//------------------------------------------------------------------
// loadSab(sab)
//
// Copy relevant data from the passed sab to this L1A's data fields.
// If the passed SAB is empty or partially complete, an exception is thrown.
// If the passed SAB has not been decoded, an exception is thrown.
//------------------------------------------------------------------

void L1B::loadSab(const Sab& sab) 
  {
  if (sab.partial()){
    L1BError e(L1B::incomplete_sab);
    e.throwMe();
  }
  if (sab.empty()){
    L1BError e(L1B::empty_sab);
    e.throwMe();
  }
  if (!sab.decoded()){
    L1BError e(L1B::undecoded_sab);
    e.throwMe();
  }

  quality_flag=0;//assume good data set
 
  //validity of sab
  if(!sab.Valid()) bitset(quality_flag,5,1);
  
    
 //header field
  sync = sab.sync;
  sclk= sab.sclk;

  // decode CDS pickup rate
  double value;
  if ( sab.scpr==0x03 ) value=30400.0;
  else if (sab.scpr==0x0C ) value=364800.0;
  scpr= Uvar(value,"Hz"); // bits/sec
  
  brst= sab.brst;

  header_tfi=sab.header_tfi;
  header_tnc=sab.header_tnc;
  header_typ=sab.header_typ;
  header_tca=sab.header_tca;
  header_tcb= sab.header_tcb;
  header_tcc= sab.header_tcc;

  pwri=sab.pwri;
  vicc=sab.vicc;
  vimc=sab.vimc;
  tail_len=sab.tail_len;
  tail_id=sab.tail_id;
  sab_counter=sab.sab_counter;
  sab_len=sab.sab_len;
  fswm=sab.fswm;
  fswc=sab.fswc;
  ctbc=sab.ctbc; 
  ctrx = sab.ctrx; 

  ctps= sab.ctps;
  ctbe=sab.ctbe;
  ctps_ctbe=sab.ctps_ctbe;
  header_end=sab.header_end;

  //slow field variables
  slow_tfi=sab.slow_tfi;
  dtn= sab.dtn;
  if(dtn<32) dtn= data_take_id_;
  slow_typ= sab.slow_typ;
  csr=sab.csr;  
  r_mode= sab.r_mode; 
  slow_instruction_number= sab.slow_instruction_number; 
  bem=sab.bem; 
  baq_mode=sab.baq_mode ; 
  tro=sab.tro;
  
  rc_bw= sab.rc_bw;  
  adc=sab.adc;

  //norminal db value
  at1_db= sab.at1_db;
  at3_db= sab.at3_db;
  at4_db= sab.at4_db; 

  //calibrated db value
  
  at1= sab.at1;
  at3= sab.at3;
  at4= sab.at4;  
  //cout<<"at1 at3 at4"<< at1<<" "<<at3<<" "<<at4<<endl;
  //cout<<"at1 at3 at4"<< at1_db<<" "<<at3_db<<" "<<at4_db<<endl;
  //bit map
  at1_each = sab.at1_db_mask;
  at3_each = sab.at3_db_mask;
  at4_each = sab.at4_db_mask;
  
  rip = sab.rip + radiometer_delta_tau_;
  csd= sab.csd;
  rad = sab.rad; 

  csq= sab.csq;
  chirp_length= sab.chirp_length;  

  slow_cfs= sab.slow_cfs;

  //fast field variables
  fast_tfi= sab.fast_tfi;
  fin= sab.fin; 
  fast_typ= sab.fast_typ; 
  pul = sab.pul; 
  bii = sab.bii;
  bpd = sab.bpd; 

  pri = sab.pri; 

  if(r_mode==0 || r_mode==8){
    rwd = sab.rwd + zero_range_time_delay_scatt_; 
  }
  else if(r_mode==1 || r_mode==9){
    rwd= sab.rwd + zero_range_time_delay_alth_;
  }
  else if(r_mode==2 || r_mode==10){
    rwd= sab.rwd + zero_range_time_delay_sarl_;
  }
  else if(r_mode==3 || r_mode==11){
    rwd=sab.rwd + zero_range_time_delay_sarh_;
  }
  else{
    rwd=sab.rwd;
  }

  fast_csf = sab.fast_csf;

  //sci data field
  
  Nradar_data= sab.Ndecoded; 
  rms_radar_data = sab.rms_radar_data;
  if (ft_==active && Nradar_data != 0)
  {
    radar_data= sab.decoded_data;
  }

 

  cnt_rl= sab.cnt_rl; 
  cnt_radio= sab.cnt_radio; 
  cnt_nd= sab.cnt_nd;
  eout= sab.eout; 
  subr= (unsigned int) sab.subr; 
  space_craft_time= sab.space_craft_time;
  hip= sab.hip + noise_diode_delta_tau_; 
  cip=sab.cip + resistive_load_delta_tau_;
  
  iebtth= sab.iebtth;  
  iebttl= sab.iebttl;
  bgcalls= sab.bgcalls; 
  delvmn= sab.delvmn; 
  delvda= sab.delvda; 
  delvyr= sab.delvyr;
    
  //sub_comm_data variables 
  //subr = 0
  fwdtmp= sab.fwdtmp; 
  be1tmp=sab.be1tmp; 
  be2tmp=sab.be2tmp;
  be3tmp=sab.be3tmp;
  be4tmp=sab.be4tmp;
  be5tmp=sab.be5tmp;

 //subr = 1
  diptmp=sab.diptmp;
  rlotmp=sab.rlotmp;
  tadcal1=sab.tadcal1;
  nsdtmp=sab.nsdtmp;
  lnatmp=sab.lnatmp;
  evdtmp=sab.evdtmp;
  
  //subr = 2
  mratmp=sab.mratmp; 
  mruttm=sab.mruttm; 
  dcgttm=sab.dcgttm;
  cucttm=sab.cucttm;
  twttmp=sab.twttmp;
  epctmp=sab.epctmp;
  
  //subr = 3
  tw1ttm=sab.tw1ttm;
  ep1ttm=sab.ep1ttm;
  p_stmp=sab.p_stmp;
  p_sttm=sab.p_sttm;
  fguttm=sab.fguttm;
  tadcal4=sab.tadcal4;
  
 //subr = 4
  esstmp=sab.esstmp;
  wgb1t1=sab.wgb1t1;
  wgb3t1=sab.wgb3t1;
  wgb3t2=sab.wgb3t2;
  wgb3t3=sab.wgb3t3;
  wgb5t1=sab.wgb5t1;
  
 //subr = 5
  pcutmp=sab.pcutmp;
  adctmp=sab.adctmp;
  tadcal2=sab.tadcal2;
  ecltmp=sab.ecltmp;
  cputmp=sab.cputmp;
  memtmp=sab.memtmp;
  
 //subr = 6
  sadctmp=sab.sadctmp;

 //subr = 7
  tadcal3=sab.tadcal3; 
  frwdpw=sab.frwdpw;
  dcgmon=sab.dcgmon;
  lpltlm_db=sab.lpltlm_db;
  

 //subr = 8
  nsdcur=sab.nsdcur;
  hpapsm=sab.hpapsm; 
  catcur=sab.catcur;
  p_smon=sab.p_smon;
  svlsta=sab.svlsta;
  usotmp=sab.usotmp;
  
//subr = 9
  cpbnkv= sab.cpbnkv; 
  essvlt=sab.essvlt;
  tadcal5=sab.tadcal5;
  pcu5v_72=sab.pcu5v_72;
  pcu5i_73=sab.pcu5i_73;
  pcu5v_74=sab.pcu5v_74;

//subr = 10
  pcuti_75=sab.pcuti_75;
  pcu15v_76=sab.pcu15v_76;
  pcu15i_77=sab.pcu15i_77;
  pcu15v_78=sab.pcu15v_78;
  pcu15i_79=sab.pcu15i_79; 
  pcu12v_80=sab.pcu12v_80;

//subr = 11
  pcu12i_81=sab.pcu12i_81;
  pcucur=sab.pcucur; 
  pllmon=sab.pllmon; 
  ctu5i=sab.ctu5i;
  tadcal6=sab.tadcal6;
  pcu9v_86=sab.pcu9v_86;
 
//subr = 12
  pcu9i_87=sab.pcu9i_87;
  pcu9v_88=sab.pcu9v_88;
  pcu9i_89=sab.pcu9i_89; 

//subr = 13
  tadcl7= sab.tadcl7 ;

//subr = 14 
  shpttm=sab.shpttm;

 
  //num_bursts_in_flight=1
  //to change this one, the user must call
  // setNum_bursts_in_flgith separately
  num_bursts_in_flight=1;
 
  data_read_=true;
 
  }


//----------------
//number of bursts
//----------------
void L1B::setNum_bursts_in_flight(const unsigned int& num_bursts)
  {
    if(!data_read_) ErrorMessage("L1B.cpp:setNum_bursts_in_flight: no data has been loaded").throwMe();
    num_bursts_in_flight= num_bursts;
  }

//======================================================================
// recordStartTime()
//
// Update the start time to be the timestamp of the current SAB.
//======================================================================
void L1B::recordStartTime(Time& t)
{
  sclk_start_ = Time(t);
}

//======================================================================
// recordStopTime()
//
// Update the stop time to be the timestamp of the current SAB.
//======================================================================
void L1B::recordStopTime(Time& t)
{
  sclk_stop_ = Time(t);
}

//======================================================================
// recordRadarMode()
//
// Update the cumulative radar mode with the mode of the current SAB.
// Use csr and adc to determine the mode instead of r_mode; r_mode will
// be set correctly for most of the observations but may produce the
// wrong result if RMSS was not used to generate the IEB.
//======================================================================
void L1B::recordRadarMode(unsigned int csr, Uvar adc)
{
  // Disregard the auto-gain bit
  csr = AUTO_GAIN_MASK & csr;
  switch(csr)
  {
    case 0:  // Active modes; ignore invalid values for adc
      if      (adc == altl_adc) pds_radar_mode_ |= SCATTEROMETER_MODE;
      else if (adc == alth_adc) pds_radar_mode_ |= ALTIMETER_MODE;
      else if (adc == sarl_adc) pds_radar_mode_ |= SAR_MODE;
      else if (adc == sarh_adc) pds_radar_mode_ |= SAR_MODE;
      break;

    case 6:  // Radiometer
      pds_radar_mode_ |= RADIOMETER_MODE;
      break;

    default:
      // Must be in one of the calibration modes; ignore this SAB
      break;
  }
}
//-----------------------------------------------------
// readPDSValuesFromLabel
// Used to set pds values which are normally determined from telemetry
// data during preprocessing using the strings in the label of the input
// file. Assumes config has been run and does not try to regenerate pds
// values obtains from the config file.
// The only values computed are sclk_start_, sclk_stop_, pds_radar_mode_
// and data_take_id_
//----------------------------------------------------
void L1B::readPdsValuesFromLabel()
{
  int currentpos=file_.getPosition();
  file_.setPosition(0);
  string start_time_string;
  string stop_time_string;
  PDSLabel label(file_);
  if(label.label()==EMPTY_STRING){
      L1BError e("No PDS label in " + filename_, L1B::read_error);
      e.throwMe();
  }
  if (!label.getString("START_TIME", &start_time_string))
    {
      L1BError e("START_TIME missing from " + filename_, L1B::unknown_parameter);
      e.throwMe();
    }

 if (!label.getString("STOP_TIME", &stop_time_string))
    {
      L1BError e("STOP_TIME missing from " + filename_, L1B::unknown_parameter);
      e.throwMe();
    }

 sclk_start_=Time(start_time_string);
 sclk_stop_=Time(stop_time_string);

 // PRODUCT_ID, record_count etc  already read from header during readHeader()
 string prmstr=product_id_.substr(5,2);
 string dtstr=product_id_.substr(9,3);
 pds_radar_mode_=atoi(prmstr.c_str());
 data_take_id_=atoi(dtstr.c_str());
 file_.setPosition(currentpos);
}
//-----------------------------------------------------
// copyPDSValues
// Used to set pds values which are normally determined from telemetry
// data during preprocessing using another L1B object
// Assumes config has been run and does not try to regenerate pds
// values obtains from the config file.
// The only values copied are sclk_start_, sclk_stop_, pds_radar_mode_
// and data_take_id_
//----------------------------------------------------
void L1B::copyPdsValues(const L1B& dat)
{
  sclk_start_=dat.sclk_start_;
  sclk_stop_=dat.sclk_stop_;
  data_take_id_=dat.data_take_id_;
  pds_radar_mode_=dat.pds_radar_mode_;
}


//-----------------------------------------------------
// copyPdsEndValues
// Used to set pds end time values which are normally determined from telemetry
// data during preprocessing using another L1B object
// Assumes config has been run and does not try to regenerate pds
// values obtains from the config file.
// The only value copied is sclk_stop_
//----------------------------------------------------
void L1B::copyPdsEndValues(const L1B& dat)
{
  sclk_stop_=dat.sclk_stop_;
}

//---------------------------------------
// combineRecordCounts()
// computes the combined record count when merging two L1B files
//----------------------------------------
void L1B::combineRecordCounts(const L1B& dat1, const L1B& dat2)
{
  record_count_=dat1.record_count_+dat2.record_count_;
}

//----------------------------
//map the record accroding to sab counter and time
//----------------------------
void L1B::mapRecords()
{
  if (mode_ == "w" || mode_ == "wb")
    {
      L1BError e("Can't read from output file " + filename_,
		 L1B::read_error);
    }
  if(!header_handled_)
    readHeader();//read header

  //------------------
  //useful internal parameter
  // record_size_
  //map<unsigned int, int> sab_counter_map_;
  //map<Uvar, int> CA_relative_time_map_;
  //-------------------
  gotoFirstRecord();
  int start_position=file_.getPosition();//after reading header
  unsigned  int Nrecord = recordCount();//count record

  //recycled parameter
  Parameter* sab=parameters_.getParam("sab_counter");
  Parameter* sclk_tmp=parameters_.getParam("sclk");
  Parameter* brst_tmp=parameters_.getParam("brst");

  int sab_location = sab->offset;
  int sclk_location = sclk_tmp->offset;
  int brst_location = brst_tmp->offset;

  int file_position, tmp_position;
  unsigned int sab_counter;
  unsigned int tmp_sclk;
  float brst_in_s;
  Time t_tmp;

  //clear internal map
  record_counter_map_.clear();
  sab_counter_map_.clear();
  time_map_.clear();

  for(unsigned int i=0; i<Nrecord;++i){
    file_position=  start_position + record_size_*i;//beginning of record
    
    //read sab counter
    tmp_position = file_position + sab_location;
    file_.setPosition(tmp_position);
    file_.read(sab_counter);

    //read sclk
    tmp_position=file_position+sclk_location;
    file_.setPosition(tmp_position);
    file_.read(tmp_sclk);//unsigned int

    //read brst
    tmp_position= file_position+brst_location;
    file_.setPosition(tmp_position);
    file_.read(brst_in_s);

    //construct time
    t_tmp.setSclk("Cassini",tmp_sclk);
    t_tmp += Uvar((double)brst_in_s,"s");
    
    record_counter_map_[i]= sab_counter;// record number vs sab counter
    sab_counter_map_[sab_counter]= file_position;// sab counter vs position
    time_map_[t_tmp]=sab_counter;// time vs sab counter
    //cout<<"sclk brst tmp time "<<tmp_sclk<<" " << brst_in_s<<" "<<t_tmp.utc("ISOD")<<endl;   
  }

  if( (sab_counter_map_.size() != Nrecord) ||
      (time_map_.size()!=Nrecord) ||
      (record_counter_map_.size()!=Nrecord))
    {
      cout<<"Mapping error"<<endl;
      cout<<"size of record counter, sab counter and time map container "<<record_counter_map_.size()<<" "<< sab_counter_map_.size()  <<" "<< time_map_.size()<<endl;
      map_complete_=false;
    }
  else
    {
      map_complete_=true;
    }
  gotoFirstRecord();
}

//----------------------
//foundRecord: return true and move file pointer to 
// the beginning of record
//-----------------------
bool L1B::foundSABRecord(const unsigned int& sab_counter)
{
  if(!map_complete_) ErrorMessage("L1B.cpp::foundSABRecord: no mapping is done").throwMe();
  bool findRecord=false;
  sab_counter_map_ptr_= sab_counter_map_.find(sab_counter);
  if(sab_counter_map_ptr_ != sab_counter_map_.end()){
    findRecord= true;//yes there is a record
    file_.setPosition(sab_counter_map_ptr_->second);
  }
  return(findRecord);
}

void L1B::findSABCounter_at_T(const Time& t, unsigned int& sab_counter,
			bool& foundRecord)
{
  if(!map_complete_) ErrorMessage("L1B.cpp::foundSABRecord: no mapping is done").throwMe();
  foundRecord=false;
  sab_counter=0;
 
  for(time_map_ptr_ = time_map_.begin(); time_map_ptr_ != time_map_.end();++time_map_ptr_){
    if( time_map_ptr_->first > (t-Uvar(1*mstos,"s")) &&
	time_map_ptr_->first < (t+Uvar(1*mstos,"s"))){
      foundRecord=true;
      break;
    }
  }

  //map.find does not work because we save brst as float although sab 
  // use double value

  if(foundRecord){
    sab_counter= time_map_ptr_->second;
  }
 
}

//-----------------------------
//access sab counter by record number
//-----------------------------
 unsigned int& L1B::operator[] (const unsigned int& n)
 {
   if(n>=record_counter_map_.size()) ErrorMessage("L1B.cpp: no record: request number is larger than total number of record").throwMe();
   return(record_counter_map_[n]);
 }

//---------------------
//copy active data 
//---------------------
void L1B::getEcho(const unsigned int& sab_counter,
		  const unsigned int& Nradar_data,
		  fdata& echo,
		  bool& valid)
  {
    if(!foundSABRecord(sab_counter)) {
      valid=false;
      echo=0.0;
      return;
    }
    //once the record, file position should be at the beginning of the record
    echo.resize(Nradar_data);
    int fp=file_.getPosition();
    file_.setPosition(fp+L1B_PASSIVE_RECORD_LENGTH);
    for(unsigned int i=0;i<Nradar_data;++i)
      file_.read(echo(i));
    file_.setPosition(fp);//return to the beginning of the current record
    valid=true;
  }

 void L1B::getEcho_fast(const unsigned int& sab_counter, 
		    const unsigned int& Nradar_data, 
		    float* echo,
		    bool& valid)
{
  if(!foundSABRecord(sab_counter)) {
    valid=false;
    return;
  }
  //once the record, file position should be at the beginning of the record
  int fp=file_.getPosition();
  file_.setPosition(fp+L1B_PASSIVE_RECORD_LENGTH);
  float x;
  for(unsigned int i=0;i<Nradar_data;++i){
    file_.read(x);
    echo[i]=x;
  }
  file_.setPosition(fp);//return to the beginning of the current record
  valid=true;
}
