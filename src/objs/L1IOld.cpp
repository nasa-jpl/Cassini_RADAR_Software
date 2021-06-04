#include "L1IOld.h"
#include "Error.h"
#include "SARFunctions.h"
#include "config_keywords.h"
#include "Constants.h"
#include <string>

using std::cout;
using std::endl;
using std::cerr;


//-----------------
// Methods for L1IOld
//-----------------

//--------------
// Constructors
//--------------



L1IOld::L1IOld(const std::string& filename, const std::string& mode)
  : 
  BurstData(filename,mode,"active"), 
  raw_sar_data("raw_sar_data",1,1),
  Xfactor("Xfactor",1,1),
  sigma0("sigma0",1,1),
  range("range",1,1),
  doppler("doppler",1,1),
  matched_filter_time("matched_filter_time",1),
  doppler_filter_freq("doppler_filter_freq",1),
  estimate_range_for_windowing_(false),
  estimate_doppler_centroid_(false),
  use_match_filter_weights_(false),
  use_doppler_filter_weights_(false),
  data_complete_(false),
  radar_data("l1i_radar_data",1), 
  use_hensley_range_compress(true),
  use_upper_band(false),
  perform_single_burst_check(false),
  check_this_burst(false),
  dechirped("dechirped",1,1)
{
  
  header_size_=L1IOld_HEADER_LENGTH;
  record_size_=L1IOld_RECORD_LENGTH;
  raw_sar_data.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  Xfactor.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  sigma0.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  range.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  doppler.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  matched_filter_time.resize(L1IOld_RANGE_WIDTH);
  doppler_filter_freq.resize(L1IOld_DOPPLER_WIDTH);
  radar_data.resize(32*1024);
  dechirped.resize(L1IOld_RANGE_WIDTH,L1IOld_DOPPLER_WIDTH);
  for (unsigned int i = 0 ; i < L1IOld_DOPPLER_WIDTH; i++){
    doppler_filter_freq(i)=0;
    for (unsigned int j = 0 ; j < L1IOld_RANGE_WIDTH; j++){
      raw_sar_data(j,i)=0;
      Xfactor(j,i)=0;
      sigma0(j,i)=0;
      range(j,i)=0;
      doppler(j,i)=0;
      }
  }
  for (unsigned int i = 0 ; i < L1IOld_RANGE_WIDTH; i++){
    matched_filter_time(i)=0;
  }
  // append parameters 
  parameters_.appendParameter("rms_radar_data",sizeof(float),
			      Parameter::FLOAT,(void*)(&rms_radar_data));
  parameters_.appendParameter("num_range_bins",num_range_bins);
  parameters_.appendParameter("num_pulses_received",num_pulses_received);

  int nelem=L1IOld_DOPPLER_WIDTH*L1IOld_RANGE_WIDTH;
  parameters_.appendParameter("raw_sar_data",nelem*sizeof(double),
			      Parameter::FMAT,(void*) &raw_sar_data);
  parameters_.appendParameter("Xfactor",nelem*sizeof(double),
			      Parameter::UMAT,(void*) &Xfactor);
  parameters_.appendParameter("sigma0",nelem*sizeof(double),
			      Parameter::UMAT,(void*) &sigma0);
  parameters_.appendParameter("range",nelem*sizeof(double),
			      Parameter::UMAT,(void*) &range);
  parameters_.appendParameter("doppler",nelem*sizeof(double),
			      Parameter::UMAT,(void*) &doppler);
  parameters_.appendParameter("window_start_delay",window_start_delay);
  parameters_.appendParameter("window_start_idx",window_start_idx);
  parameters_.appendParameter("matched_filter_time",
			      L1IOld_RANGE_WIDTH*sizeof(double),
			      Parameter::FDATA,(void*) &matched_filter_time);
  parameters_.appendParameter("doppler_filter_freq",
			      L1IOld_DOPPLER_WIDTH*sizeof(double),
			      Parameter::FDATA,(void*) &doppler_filter_freq); 
  parameters_.appendParameter("fraction_energy_in_window",
			    fraction_energy_in_window);


  // Allocate intermediate arrays
  matched_filter_hensley = new complex<float>[32*1024];
  fft_matched_filter_hensley = new complex<float>[32*1024];
  tmp_complex_array = new complex<float>[32*1024];
  tmp_complex_array2 = new complex<float>[32*1024];
}
			     
  

L1IOld::~L1IOld()
{
  // Delete intermediate arrays
  delete[] matched_filter_hensley;
  delete[] fft_matched_filter_hensley;
  delete[] tmp_complex_array;
  delete[] tmp_complex_array2;
}
 
void L1IOld::config(Config& cfg)
{
  cassini_sim.config(cfg);
  updown_shift=-cfg["stalo_frequency"];
  filter_start_freq_offset=cfg["FILTER_START_FREQ_OFFSET"];
  transmit_power=cfg["Pt"];
  range_ref_window_param=cfg.getDouble("RANGE_REF_WINDOW_PARAM");
  DC_null_bins=cfg.getInt("DC_NULL_BINS");
  if(DC_null_bins < 0 || (DC_null_bins > 0 && DC_null_bins%2==0)) {
    ErrorMessage e("L1IOld::config Bad DC_null_Bins should be (0,1,3,5,...)");
    e.throwMe();
  }
  use_hensley_range_compress=true;
  use_upper_band=(bool)cfg.getInt("USE_UPPER_BAND");
  use_linear_chirp=(bool)cfg.getInt("SARPROC_USE_LINEAR_CHIRP");
  perform_single_burst_check=(bool)cfg.getInt("PERFORM_SINGLE_BURST_CHECK"); 
  if(perform_single_burst_check){
    check_burst_idx=cfg.getInt("CHECK_BURST_NO");
    if(check_burst_idx < 1){
      ErrorMessage e("L1IOld:config bad Check Burst Number should be >= 1");
      e.throwMe();
    }
    check_burst_filename=cfg.str("CHECK_BURST_FILENAME");
    cbf.open(check_burst_filename.c_str());
  }
}

//--------------
// I/O
//--------------


//------------------------------------------------------------------
// writeHeader()
//
// Write a L1IOld header to the current L1IOld file in binary.
// This method defines the L1IOld Header format.
// Takes a num_records parameter
// various overloads follow
//------------------------------------------------------------------

void L1IOld::writeHeader(int num_records) 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1IOldError e("Can't write to input file " + filename_, L1IOld::write_error);
      e.throwMe();
    }

  // Should not call this routine twice
  if (header_handled_)
    {
      L1IOldError e("Attempt to write header twice", L1IOld::write_error);
      e.throwMe();
    }
  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",num_records,record_size_);
  header=headbuf;
  if ( header.length() != header_size_){
    char msg[100];
    sprintf(msg,"Incorrect header size (%d bytes) should be %d",
	    header.length(),header_size_);
    L1IOldError e(msg,L1IOld::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  }

void L1IOld::writeHeader(){
  writeHeader(record_count_);
}

void L1IOld::writeHeader(BurstData& l1b){
  writeHeader(l1b.recordCount());
}

void L1IOld::rewriteHeader(){
int position = file_.getPosition();
  file_.rewind();
  string header;
  char headbuf[100];
  sprintf(headbuf,"Num Records %7.7d Record Size %7.7d\n",record_count_,record_size_);
  header=headbuf;
  if ( header.length() != header_size_){
    char msg[100];
    sprintf(msg,"Incorrect header size (%d bytes) should be %d",
	    header.length(),header_size_);
    L1IOldError e("Header length is incorrect",L1IOld::write_error);
    e.throwMe();
  }
  file_.write(header);
  
  header_handled_ = true;
  file_.setPosition(position);
}

//------------------------------------------------------------------
// readHeader()
//
// Read a L1IOld header from the current L1IOld file in binary.
// This method defines the L1IOld Header format.
//------------------------------------------------------------------

void L1IOld::readHeader() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1IOldError e("Can't read from output file " + filename_,
	       L1IOld::read_error);
      e.throwMe();
    }
  if (header_handled_)
    {
      L1IOldError e("Attempt to read header twice",
	       L1IOld::read_error);
      e.throwMe();
    }
  checkEndedness();
  string header;
  header.resize(L1IOld_HEADER_LENGTH);
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
    L1IOldError e("Record size incorrect in " + filename_, L1IOld::read_error);
    e.throwMe();
  }
  header_handled_ = true;
  }

//------------------------------------------------------------------
// writeRecord()
//
// Write the current L1IOld record to the end of the current L1IOld file in binary.
// This method defines the L1IOld record format.
// It is the responsibility of the user to ensure that valid data is present
// before calling this routine.  No tracking of data status is performed.
//------------------------------------------------------------------

void L1IOld::writeRecord() 
  {
  if (mode_ == "r" || mode_ == "rb")
    {
      L1IOldError e("Can't write to input file " + filename_, L1IOld::write_error);
      e.throwMe();
    }
  if (!header_handled_)
    {
      L1IOldError e("Can't write L1IOld record until header is written",
		 L1IOld::write_error);
      e.throwMe();
    }
  if (!data_complete_)
    {
      L1IOldError e("Can't write L1IOld record until data is complete",
		 L1IOld::write_error);
      e.throwMe();
    }

  if(file_.getPosition() > MAX_FILE_LENGTH-(int)record_size_) 
    createNextFile();


  int start_position=file_.getPosition();
  writePassiveSABData();
  writeGeometry(); 
  writeSARData();
    
  int end_position=file_.getPosition();
  int len_written=end_position-start_position;
  if((unsigned int)len_written!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_written,record_size_);
    L1IOldError(msg,L1IOld::write_error).throwMe();
  }
  record_count_++;
  }


//------------------------------------------------------------------
// readRecord()
//
// Read a L1IOld record from the current L1IOld file.
//------------------------------------------------------------------

void L1IOld::readRecord() 
  {
  if (mode_ == "w" || mode_ == "wb")
    {
      L1IOldError e("Can't read from output file " + filename_,
	       L1IOld::read_error);
//
// create next file
      e.throwMe();
    }  if (!header_handled_)
    {
      L1IOldError e("Can't read L1IOld record until header is read",
	       L1IOld::read_error);
      e.throwMe();
    } 

  // This check has the side effect of moving to the next  file if necessary
  if(eof())
    {
      L1IOldError e("Unexpected EOF in L1IOld file "+filename_,
	       L1IOld::read_error);
      e.throwMe();
    } 


 
  int start_position=file_.getPosition();
  
  readPassiveSABData();  
  readGeometry(); 
  readSARData();
  
  
  int end_position=file_.getPosition();
  int len_read=end_position-start_position;
  if((unsigned int)len_read!=record_size_){
    char msg[100];
    sprintf(msg,"Incorrect record size (%d bytes) should be %d",
	    len_read,record_size_);
    L1IOldError e(msg,L1IOld::read_error);
    e.throwMe();
  }
  data_complete_=true;
  data_read_=true;  
  }


//-------------------------------------------
// Routine for computing geometry
// for a previously extant L1B file
//-------------------------------------------

void
L1IOld::convertL1B(L1B& l1b){
  
  
  // clear quality flag and reset data_read and data_complete flags
  resetRecord();

  // Copy L1B data from L1B file
  l1b.copyTo(*this,radar_data,Nradar_data);

  burst_no=(int)(record_id%1000000);
  cout << "Processing Burst " << burst_no << ", Beam " << beam_number 
       << endl;

  cbf.set("Nradar_data",Nradar_data,burst_no);

  /***** for debugging zeroed radar buffer values
  cout << "zero_valued_radar_data Indices:" << endl;
  bool last_one_zeroed=false;
  for(unsigned int c=0;c<Nradar_data;c++){
    if(last_one_zeroed && radar_data(c)!=0.0){
      last_one_zeroed=false;
      cout << "-" << c << endl;
    }
    else if (!last_one_zeroed && radar_data(c)==0.0){
      last_one_zeroed=true;
      cout << c;
    }
  }
  if(last_one_zeroed) cout << "-" << Nradar_data -1 << endl;
  *****/

  // If geometry info is bad no processing is performed
  if(! goodGeometry() || beam_number<1 || beam_number>5){
    quality_flag+=32; // sets SAR data not processed flag
    data_complete_=true;
    return; //ignores data with bad geometry (bad spice ckernels)
  }

  // set to produce debug output if this burst is to be checked
  int burst_no=record_id%1000000;
  if(perform_single_burst_check && burst_no==check_burst_idx) 
    check_this_burst=true;
  else check_this_burst=false;
  if(Nradar_data==0){
    L1IOldError e("convertL1B: No active mode data",bad_l1b);
    e.throwMe();
  }

  rms_radar_data=radar_data.rms(Nradar_data);

  
  // Compute Anti-aliasing filter start freq;
  filter_start_freq=-adc/2+filter_start_freq_offset;

  computeDoubleParams();

  cbf.set("rms_radar_data",rms_radar_data,burst_no);

  //---------- BAQ Decode step needs to go here --------------//

  //--- Update cassini_sim -----//
  // use time halfway between transmit window center and receive window
  // center to get nominal_doppler_centroid and nominal_delay
  Uvar receive_window=l1b.ctrx*l1b.pri;
  Uvar transmit_window=l1b.pul*l1b.pri;
  Uvar transmit_center=transmit_window/2;
  Uvar receive_center=rwd+receive_window/2;
  Uvar time_offset=(transmit_center+receive_center)/2;
  Time t0=t+time_offset;
  cassini_sim.setBeamTime(beam_number,t0);
  nominal_doppler_centroid=cassini_sim.boresightDoppler();
  nominal_delay=cassini_sim.boresightRange()*2/speed_light;
  cbf.set("nom_del",nominal_delay.getInUnits("s"),burst_no);
  cbf.set("nom_dop",nominal_doppler_centroid.getInUnits("Hz"),burst_no);
  computeWindows();
  if(num_samples_per_window>L1IOld_RANGE_WIDTH ||
     num_pulses_received > L1IOld_DOPPLER_WIDTH){
    L1IOldError e("L1IOld bad sizes DopWidth=" + toStr(num_pulses_received) +
	       + "Range Width="  
	       + toStr(num_samples_per_window),
	       internal_error);
    e.throwMe();
  }

  if(num_pulses_received==0){
    quality_flag+=32; // sets SAR data not processed flag
    data_complete_=true;
    return; //ignores data with no pulses received
  }

  if(use_hensley_range_compress){
    rangeCompressHensley();
  }
  else{

    //--------------------------------------------------------------//
    // Old Range compression code                                   //
    //--------------------------------------------------------------//
    computeMatchedFilter();


    for(unsigned int i=0;i<num_pulses_received;i++){
      // unsigned int step = (unsigned int)((pri*adc).getValue()+0.5);
      unsigned int start_idx=window_start_idx+i*window_step;
      CFvec dc("dc",1);
      range_compress(radar_data,start_idx,num_samples_per_window,
		     matched_filter_time,dc);
      num_range_bins=num_samples_per_window;
      for(unsigned int j=0;j<num_samples_per_window;j++){
	dechirped(j,i)=dc(j); 
      }
    }
  }

   
  //-------------------------------------------------------------//
  // Doppler compression code
  //-------------------------------------------------------------//

  for(unsigned int j=0;j<num_samples_per_window;j++){
    CFvec s=dechirped.getRow(j);
    Uvar doppler_bw=1/pri; 
    // freq shift dechirped signal to correct frequency range
    for(unsigned int i=0;i<num_pulses_received;i++){
      complex<float> shift=exp(complex<float>(0,-1)*(float)(2.0*pi*
      	   ((nominal_doppler_centroid-doppler_bw/2)*pri).getValue()
      	    *i));

      // WE NEED TO FIGURE OUT WHY THIS WORKS !!!!!!
      if(use_hensley_range_compress) s(i)=conj(s(i));

      s(i)=s(i)*shift;
		     
    }
    CFvec S("S",num_pulses_received);
    dft(s,S,num_pulses_received);
    for(unsigned int i=0;i<num_pulses_received;i++){
       raw_sar_data(j,i)=norm(S(i));  
    }
  }
  //----------------------------------------------------------------- 
  // compute doppler and range for each pixel
  // (Need to add range calibration here)
  //-----------------------------------------------------------------
  computeDopplerAndRange();

  //-----------------------------------------------------------------
  // calibrate (compute Xfactor and sigma0)
  //-----------------------------------------------------------------
  calibrate();

  /**** for debugging print out dop and delay estimates of max energy pixel
  unsigned int imax=0;
  unsigned int jmax=0;
  float smax=0.0;  
  for(unsigned int j=0;j<num_range_bins;j++){
    for(unsigned int i=0;i<num_pulses_received;i++){
      if(sigma0(j,i).getInUnits("")> smax){
	smax=sigma0(j,i).getInUnits("");
	jmax=j;
	imax=i;
      }
    }
  }
  cout << "Doppler of peak is:" << doppler(jmax,imax) << endl;
  cout << "Range of peak is:" << range(jmax,imax) << endl;
  ****/
  /***** FOR DEBUGGING WE CHECK FFT
  CFvec x("x",8);
  CFvec X("X",8);
  for(unsigned int i=0;i<8;i++){
    if(i<4) x(i)=1;
    else x(i)=0;
  }
  fft(x,X);
  ifft(X,x);
  cout << "x:" << x <<endl;
  cout << "X:" << X <<endl;
  ****/


  //------------------------------------------------------------//
  // Set Data complete
  //------------------------------------------------------------//
  data_complete_=true;
  if(check_this_burst) cbf.flush();
}



void L1IOld::resetRecord(){
  quality_flag=0;
  data_read_=false;
  data_complete_=false;
}




// Parameter access methods

double L1IOld::getSpecialParam(const char* name){
  double retval;
  
  retval=BurstData::getSpecialParam(name);
  return(retval);
}
void L1IOld::enableSpecialParam(const char* name)
{
  BurstData::enableSpecialParam(name);
}
void L1IOld::disableSpecialParam(const char* name)
{
  BurstData::disableSpecialParam(name);
}

 

void L1IOld::readSARData()
{
  // incomplete
  file_.read(rms_radar_data);
  file_.read(num_range_bins);
  file_.read(num_pulses_received);

  unsigned int pos;
  for(unsigned int i=0;i<num_range_bins;i++){
    for(unsigned int j=0;j<num_pulses_received;j++){
      file_.read(raw_sar_data(i,j));
    }
    pos=file_.getPosition();
    file_.setPosition(pos+(L1IOld_DOPPLER_WIDTH-num_pulses_received)*sizeof(float));
  }
   pos=file_.getPosition();
   unsigned int offset=L1IOld_DOPPLER_WIDTH*
     (L1IOld_RANGE_WIDTH-num_range_bins)*sizeof(float);
   file_.setPosition(pos+offset);

   file_.read(window_start_idx);
   for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
     file_.read(matched_filter_time(i));
   }

  for(unsigned int i=0;i<num_range_bins;i++){
    for(unsigned int j=0;j<num_pulses_received;j++){
      file_.read(dechirped(i,j));
    }
    pos=file_.getPosition();
    file_.setPosition(pos+(L1IOld_DOPPLER_WIDTH-num_pulses_received)*2*sizeof(float));
  }
   pos=file_.getPosition();
   offset=L1IOld_DOPPLER_WIDTH*
     (L1IOld_RANGE_WIDTH-num_range_bins)*2*sizeof(float);
   file_.setPosition(pos+offset);
}
  

void L1IOld::writeSARData()
{
  // incomplete
  file_.write(rms_radar_data);
  file_.write(num_range_bins);
  file_.write(num_pulses_received);
  for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
    for(unsigned int j=0;j<L1IOld_DOPPLER_WIDTH;j++){
      file_.write(raw_sar_data(i,j));
    }
  }
  file_.write(window_start_idx);
  for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
    file_.write(matched_filter_time(i));
  }
  for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
    for(unsigned int j=0;j<L1IOld_DOPPLER_WIDTH;j++){
      file_.write(dechirped(i,j));
    }
  }

  for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
    for(unsigned int j=0;j<L1IOld_DOPPLER_WIDTH;j++){
      range(i,j).writeFloat(file_,0);
    }
  }

  for(unsigned int i=0;i<L1IOld_RANGE_WIDTH;i++){
    for(unsigned int j=0;j<L1IOld_DOPPLER_WIDTH;j++){
      doppler(i,j).writeFloat(file_,0);
    }
  }
}

void L1IOld::computeWindows()
{
  // for now arbitrarily set nominal spread to 1 PRI 
  // nominal_pulse_spread=pri-chirp_length;
  nominal_pulse_spread=1.5*pri;

  // Commented out to allow windows to overlap
  // if(chirp_length+nominal_pulse_spread > pri){
  //  L1IOldError e("Windows Overlap.",internal_error);
  //  e.throwMe();
  //}

  window_start_delay=nominal_delay-nominal_pulse_spread/2-rwd;

  // Omit missed and partial pulses at beginning
  num_pulses_received=pul;
  first_pulse_idx=0;
  while(window_start_delay<0){
    num_pulses_received--;
    window_start_delay+=pri;
    first_pulse_idx++;
  }

  window_start_idx=(unsigned long int)((window_start_delay*adc).getValue()+0.5);

  
  // num_samples_per_window=(unsigned long int)((pri*adc).getValue()+0.5);
  num_samples_per_window=(unsigned long int)
    (((chirp_length+nominal_pulse_spread)*adc).getValue()+0.5);
  window_step=(unsigned long int)((pri*adc).getValue()+0.5);


  // Omit missed or partial pulses at the end
  unsigned long int end_sample=window_start_idx+
    (num_pulses_received-1)*window_step
    +num_samples_per_window-1;
  cout << end_sample << " " << window_start_idx << " " << num_pulses_received
       <<" "<< num_samples_per_window << endl;

  if(window_start_idx + num_samples_per_window -1 > Nradar_data)
    {
      num_pulses_received=0;
    }
  else{
    while(end_sample>Nradar_data){
      num_pulses_received--;
      end_sample=window_start_idx+
	(num_pulses_received-1)*window_step
	+num_samples_per_window-1;
      cout << end_sample << " " << Nradar_data << endl;
    }
  }


  Nfft = (int)get_next_power_of_2((unsigned int)num_samples_per_window);
  if(check_this_burst){
    cbf.set("pulse_spread",nominal_pulse_spread.getInUnits("s"));
    cbf.set("start_delay",window_start_delay.getInUnits("s"));
    cbf.set("start_idx",window_start_idx);
    cbf.set("window_len",num_samples_per_window);
    cbf.set("Nfft",Nfft);
    cbf.set("window_step",window_step);
    cbf.set("num_pulses",num_pulses_received);
    cbf.set("num_data",Nradar_data);
    cbf.set("end_idx",end_sample);
    cbf.set("pri",pri.getInUnits("s"));
  }
 
} 
 
void L1IOld::computeMatchedFilter()
{
  // This needs to be modified to handle the case in which 
  // the analog anti-aliasing filter clips part of the data
  // Clipping will vary for different ranges if doppler varies
  // significantly with range

  // compute freq_shift = updown_shift
  if(!use_linear_chirp){
    ErrorMessage
    e("L1IOld::ComputeMatchedFilter: Old rangecompress needs linear chirp");
    e.throwMe();
  }
  Uvar freq_shift=updown_shift;
  Uvar chirp_rate=slow_cfs/csd;
  Uvar f_start=freq_shift+fast_csf+nominal_doppler_centroid;
  cbf.set("updown_shift",updown_shift.getInUnits("Hz"),burst_no);
  cbf.set("chirp_rate",chirp_rate.getInUnits("Hz/s"),burst_no);
  cbf.set("fast_csf",fast_csf.getInUnits("Hz"),burst_no);
  cbf.set("f_start",f_start.getInUnits("Hz"),burst_no);
  cbf.set("rc_bw",rc_bw.getInUnits("Hz"),burst_no);
  for(unsigned int i=0;i<num_samples_per_window;i++){
    Uvar tp=i/adc;
    if(tp<chirp_length){
      cout << "MFTAG " << i << " " ;
      matched_filter_time(i)=complex_anti_aliased_linear_fm_chirp(tp,f_start,
							chirp_rate,1.0,0.0,
							filter_start_freq,
							filter_start_freq+
								  rc_bw);
      
    }
    else{ 
      matched_filter_time(i)=0.0;
     }

  }
  num_range_bins=num_samples_per_window;
}
void L1IOld::rangeCompressHensley()
{

  if(check_this_burst) 
    cbf.comment("Using Hensley's range compression ...");

  if(check_this_burst){
      cbf.comment("Number of Samples in each window");
      cbf.set("num_samples_per_window",num_samples_per_window);
      cbf.comment("Number of Pulses Received");
      cbf.set("num_pulses_received",num_pulses_received);
  }

  // compute reference function (matched filter)
  computeMatchedFilterHensley();

  // for each pulse perform range compression
  for(unsigned int i=0;i<num_pulses_received;i++){

    int pulse_offset=(int)(window_start_idx+i*window_step);
    // copy and zero-pad signal 
    for(int j=0;j<(int)num_samples_per_window;j++){
      tmp_complex_array[j]=complex<float>(radar_data(pulse_offset+j),0.0);
    }
    for(int j=(int)num_samples_per_window;j<Nfft;j++){
      tmp_complex_array[j]=complex<float>(0,0);
    }
   
    if(check_this_burst){
      cbf.comment("Offset into echo array");
      cbf.set("pulse_offset",pulse_offset,i+1);
      cbf.comment("");
      cbf.comment("Zero-pad echo data in window");
      cbf.set("hrcdata1",tmp_complex_array,Nfft,i+1,-1);
    }

    // perform FFT on signal
    fft(tmp_complex_array,tmp_complex_array2,Nfft);

    if(check_this_burst){
      cbf.comment("");
      cbf.comment("Perform FFT on echo data");
      cbf.set("hrcdata2",tmp_complex_array2,Nfft,i+1,-1);
    }

    if(!use_upper_band){
      //----------------------------------------------------------
      // Use Lower sideband
      //----------------------------------------------------------

      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Using lower sideband ...");
      }

      // multiply frequency domain signal and conjugate of reference function
      for(int j=0; j<Nfft/2; j++){ 
	tmp_complex_array2[j]*=conj(fft_matched_filter_hensley[j]);
      }

      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Multiply in frequency domain by conj(ref_func)");
	cbf.set("hrcdata3",tmp_complex_array2,Nfft,i+1,-1);
      }      

      // baseband (rotate spectra)
      for(int j=0;j<Nfft/4;j++) 
	tmp_complex_array2[j+Nfft/2]=tmp_complex_array2[j];
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Rotate Spectra");
	cbf.set("hrcdata4",tmp_complex_array2,Nfft,i+1,-1);
      }      

      // inverse FFT 
      ifft(&tmp_complex_array2[Nfft/4],tmp_complex_array,Nfft/2);
     
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Inverse FFT");
	cbf.set("hrcdata5",tmp_complex_array,Nfft/2,i+1,-1);
      }      
      
      // copy to storage array
      for(int j=0;j<Nfft/2;j++){
	dechirped(j,i)=tmp_complex_array[j]; 
      }
    }
    else{
      //------------------------------------------------------------------
      // Use Upper sideband
      //------------------------------------------------------------------
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Using lower sideband ...");
      }

      // multiply frequency domain signal and conjugate of reference function
      for(int j=Nfft/2; j<Nfft; j++){ 
	tmp_complex_array2[j]*=conj(fft_matched_filter_hensley[j]);
      }
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Multiply in frequency domain by conj(ref_func)");
	cbf.set("hrcdata3",tmp_complex_array2,Nfft,i+1,-1);
      }      

      // baseband (rotate spectra)
      for(int j=Nfft/4;j<Nfft/2;j++) 
	tmp_complex_array2[i]=tmp_complex_array2[i+Nfft/2];


      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Rotate Spectra");
	cbf.set("hrcdata4",tmp_complex_array2,Nfft,i+1,-1);
      }      
      
      // inverse FFT 
      ifft(&tmp_complex_array2[Nfft/4],tmp_complex_array,Nfft/2);
      
      if(check_this_burst){
	cbf.comment("");
	cbf.comment("Inverse FFT");
	cbf.set("hrcdata5",tmp_complex_array,Nfft/2,i+1,-1);
      }      
      
      // copy to storage array
      for(int j=0;j<Nfft/2;j++){
	dechirped(j,i)=tmp_complex_array[j]; 
      }
    }
  }
  // reduce num_samples_per_window
  // number of real samples ==> number of complex samples
  num_range_bins=num_samples_per_window/2; 
}

void L1IOld::computeMatchedFilterHensley()
{

  //    initialize matched filter and zero pad to end
  for(int i=0;i<Nfft;i++){
    matched_filter_hensley[i]=complex<float>(0.0,0.0);
  }     
      
  //    compute the chirp bandwidth and center frequency
  double chirp_bandwidth = (csq)*slow_cfs_in_Hz;  
  double r_f0 = updown_shift_in_Hz+fast_csf_in_Hz+chirp_bandwidth/2;
  r_f0=r_f0+nominal_doppler_centroid.getInUnits("Hz");

  // chirp phase at middle of pulse
  double r_phase0 = 0;

  //------------------------------------------------------------------- 
  //     generate the time domain signal 
  //------------------------------------------------------------------- 
  int i_npts=0;   
  if(!use_linear_chirp){
    dig_chirp_hensley(chirp_length_in_s,chirp_bandwidth,csd_in_s,r_f0,
       r_phase0,adc_in_Hz,i_npts,matched_filter_hensley);
  }
  else{
    dig_chirp_hensley(chirp_length_in_s,chirp_bandwidth,
    chirp_length_in_s,r_f0,
    r_phase0,adc_in_Hz,i_npts,matched_filter_hensley);
  }
  


  //-----------------------------------------------------------------     
  //     compute scale of reference function 
  //-----------------------------------------------------------------     
      
  double fftgain=0;
  tmp_complex_array[0]=complex<float>(1,0);
  for(int i=1;i<Nfft;i++) tmp_complex_array[i]=complex<float>(0,0);

  fft(tmp_complex_array,tmp_complex_array2,Nfft);
  fftgain=real(tmp_complex_array2[0]);
  fftgain=fftgain*fftgain;

  ifft(tmp_complex_array2,tmp_complex_array,Nfft/2);
  fftgain=fftgain*real(tmp_complex_array[0]);

  fftgain = fftgain*i_npts;

  //-------------------------------------------------------------
  //    Apply window to range reference function
  //-------------------------------------------------------------

  double a=range_ref_window_param;
  double b=1-a;
     

  for (int i = 0;  i < i_npts; i++){
    double r_win = a - b * cos((2.* pi*i)/(i_npts-1));
    matched_filter_hensley[i] *= r_win/fftgain;
  }
         
         
  //------------------------------------------------------
  // calculate fft of range reference function
  //-------------------------------------------------------       
  fft(matched_filter_hensley,fft_matched_filter_hensley,Nfft); 

  
  //----------------------------------------------------    
  //    zero out DC from spectrum
  //----------------------------------------------------
  
  if(DC_null_bins != 0){       
    for(int i = -DC_null_bins/2; i<= DC_null_bins/2; i++){
      double wgt = 0.5 - 0.5 * cos((float)i/(DC_null_bins/2)*pi);
      int ind = (i+Nfft)%Nfft;
      fft_matched_filter_hensley[ind] *= wgt;
    }
  
    //-------------------------------------------------------
    //    Zero out mid bandwidth DC term in offset video data
    //-------------------------------------------------------

    for(int i = -DC_null_bins/2; i<= DC_null_bins/2; i++){
      double wgt = 0.5 - 0.5 * cos((float)i/(DC_null_bins/2)*pi);
      int ind = i+Nfft/2;
      fft_matched_filter_hensley[ind] *= wgt;
    }
  } // end if DC_null_bins nonzero

  // copy to output buffer
  for(unsigned int i=0;i<num_samples_per_window;i++) 
    matched_filter_time(i)=matched_filter_hensley[i];
}

void
L1IOld::computeDoubleParams()
{
 rc_bw_in_Hz=rc_bw.getInUnits("Hz");
 adc_in_Hz=adc.getInUnits("Hz");
 chirp_length_in_s=chirp_length.getInUnits("s");
 rwd_in_s=rwd.getInUnits("s");
 fast_csf_in_Hz=fast_csf.getInUnits("Hz");
 updown_shift_in_Hz=updown_shift.getInUnits("Hz");
 pri_in_s=pri.getInUnits("s");
 slow_cfs_in_Hz=slow_cfs.getInUnits("Hz");
 csd_in_s=csd.getInUnits("s");
}

void L1IOld::computeDopplerAndRange()
{
  // Right now this is simple and stupid it just assumes that the doppler
  // shift stay in the bandwidth for all ranges
  // compute Ranges
  for(unsigned int i=0;i<num_range_bins;i++){
    int mult=num_samples_per_window/num_range_bins;
    Uvar td=mult*i/adc+window_start_idx/adc+rwd-pri*first_pulse_idx;
    Uvar r=td*speed_light/2;
    for(unsigned int j=0;j<num_pulses_received;j++){
      range(i,j)=r;
    }
  }

  // compute Doppler
  Uvar doppler_bw=1/pri;
  Uvar doppler_step=doppler_bw/num_pulses_received;
  for(unsigned int i=0;i<num_pulses_received;i++){
    for(unsigned int j=0;j<num_range_bins;j++){
      doppler(j,i)=i*doppler_step-doppler_bw/2+nominal_doppler_centroid;
    }
  }  
}
  
void L1IOld::calibrate()
{
  // Right now it stupidly set Xfactor to 1
  for(unsigned int i=0;i<num_range_bins;i++){
    for(unsigned int j=0;j<num_pulses_received;j++){
      Xfactor(i,j)=Uvar(1,"");
      sigma0(i,j)=Xfactor(i,j)*raw_sar_data(i,j);
    }
  }
}





