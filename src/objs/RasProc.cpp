//-----------------------------------------------------------------------------
// RasProc.cpp: simple range and doppler compression for doppler cetroid
// and ras tools
//-----------------------------------------------------------------------------

#include <string>
#include <fstream>
//----------------------------
//Forward declaration
//----------------------------
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "RasProc.h"
#include "SARFunctions.h"
#include "Plot.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//-------------------------
// RasProc Methods
//-------------------------
//-------------------------
// Constructors
//-------------------------
RasProc::RasProc(const string& filename, const string& mode)
  :RasList(filename, mode),
   matched_filter("matched filter"),
   fft_matched_filter("fft_matched filter"),
   echo2D_rc("echo 2D rc"),
   echo2D_rdc("echo 2D range dopp comprssed"),
   c_corr2D("2D correlation"),
   c_corr1D("1D correlation"),
   dop_frac1D("dop_cent1D"),
   range_axis("range axis"),
   doppler_axis("doppler axis"),
   ieb_loaded_(false),
   matched_filter_set_(false),
   echo2D_set_(false)
{}

RasProc::~RasProc()
  {
    
  }
//--------------------
//load Ieb
//---------------------
void RasProc::loadIeb(Ieb& ieb)
{
  ieb_ = ieb;
  ieb_loaded_ = true;
  
  decodeEcho(ieb_);
  matched_filter_set_ = false;
  echo2D_set_ = false;
}

//------------------------
//simple matched filter and range/azimuth image
//----------------------
void RasProc::computeMatchedFilter(const Uvar& doppler_centroid)
 {
  doppler_centroid_ = doppler_centroid;
  if(!ieb_loaded_) ErrorMessage("RasList.cpp:No ieb loaded").throwMe();
  doppler_centroid_ = doppler_centroid;
  unsigned int M,N;
  echo2D.size(M,N);
  matched_filter.resize(N);
  fft_matched_filter.resize(N);
  
  matched_filter=complex<float> (0.0,0.0);
  fft_matched_filter=complex<float>(0.0,0.0);
  unsigned int i_npts=0;
  double adc=ieb_.getAdc().getInUnits("Hz");
  double taup=ieb_.getTaup().getInUnits("s");
  double chirp_rate = (ieb_.getCfs()/ieb_.getCsd()).getInUnits("1/ s s");
  double cfs=ieb_.getCfs().getInUnits("Hz");
  double start_frequency= -slo_frequency.getInUnits("Hz")
  +ieb_.getChirpStartFrequency().getInUnits("Hz")
  +doppler_centroid.getInUnits("Hz");
  double f1= -adc/2.0 + 0.01*adc/2.0;
  double f2= - 0.01*adc/2.0;
  compute_matched_filter(adc, 
  		 taup, 
  		 chirp_rate,
  		 cfs,
  		 f1,
  		 f2,
  		 start_frequency,
  		 N,
  		 i_npts, 
  		 matched_filter);

  //cout<<"adc taup chirprate f1 f2 start frequency"<< adc<< " "<<taup
  //  <<" "<<chirp_rate<<" "<<cfs<<" "<<f1<<" "<<f2<<" "
  //<<start_frequency<<" "<<N<<endl;
  //cout<<"inpts "<<i_npts<<endl;
  //this part from B.S. code
  //need to understand this part more

  double fftgain=0.0;
  complex<float>* tmp_complex_array;
  complex<float>* tmp_complex_array2;
  tmp_complex_array = new complex<float> [N];
  tmp_complex_array2 = new complex<float> [N];
  tmp_complex_array[0]= complex<float> (1.0,0.0);
  for(unsigned int i=1;i<N;++i) tmp_complex_array[i]=complex<float> (0.0,0.0);
  fft(tmp_complex_array, tmp_complex_array2,N);
  fftgain = real(tmp_complex_array2[0]);
  fftgain = fftgain * fftgain;

  ifft(tmp_complex_array2, tmp_complex_array,N/2);
  fftgain= fftgain * real(tmp_complex_array[0]);
  fftgain = fftgain * double(i_npts);
  for(unsigned int i=0;i<N;++i) 
    matched_filter(i)=matched_filter(i)/(float)fftgain;
  
  cout<<"matched filter gain "<< fftgain<<endl;
  fft(matched_filter,fft_matched_filter);
  matched_filter_set_ = true;

  delete[] tmp_complex_array;
  delete[] tmp_complex_array2;
  //debug
  //displayMatchedFilter();
  }

//-------------------------------
//conver to 2D array
// compute starting point and ending point
// compute number points inside window
// compute number of pulses
//---------------------------------
void RasProc::convert2Decho(const Uvar& range_start,
			    const Uvar& range_end)
{
  if(range_end < range_start) ErrorMessage("RasProc: range end is smaller than range start").throwMe();

  //compute start index
  range_start_ = range_start;
  if(!ieb_loaded_) ErrorMessage("RasProc::ieb is not loaded").throwMe();
  Uvar dt = 2.0*range_start_/speed_light - ieb_.getRwd();
  int start_index = (int) round_double((dt*ieb_.getAdc()).getInUnits(""));
  
  //if dn is negative, discard partial pulses
  while(start_index < 0){
    cout<<"discarding one pulse of data "<<endl;
    start_index += (int) ieb_.getPointsInsidePRI();
  }

  //compute end index
  range_end_ = range_end;
  dt = double(ieb_.getPul()-1) * ieb_.getPri() + ieb_.getTaup() ;
  dt += 2.0 * range_end_/speed_light;
  dt -= ieb_.getRwd();
  int end_index = (int) round_double( (dt*ieb_.getAdc()).getInUnits(""));
  while(end_index > (int) ieb_.getRadarDataPoints()){
    cout<<"discarding one pulse of data from echo tail "<<endl;
    end_index -= (int) ieb_.getPointsInsidePRI();
  }
  if(end_index < start_index) ErrorMessage("no useful data inside echo").throwMe();

  N_pulses_inside_echo_ = ((unsigned int)(end_index-start_index)) /ieb_.getPointsInsidePRI();



  //-----------------------------------------------------------
  //determine how many points we want to hold inside each record
  // compute pulse gate
  //------------------------------------------------------------
  Uvar hold_time = 2.0*(range_end - range_start)/speed_light+ieb_.getTaup();
  if(hold_time < ieb_.getPri()) hold_time = ieb_.getPri();
  N_samples_per_window_ = (unsigned int) round_double((hold_time*ieb_.getAdc()).getInUnits(""));
  if(N_samples_per_window_%2!=0) N_samples_per_window_++;

 
  get2Decho((unsigned int) start_index,
	    (unsigned int) end_index,
	    N_samples_per_window_,
	    N_pulses_inside_echo_);//transform into 2D array
  
  echo2D_set_ = true;

  //debug
  
  cout<<"start_index end_index useful_pulse comm_pulse "<< start_index<<" "<<end_index<<" "<< N_pulses_inside_echo_ <<" "<<ieb_.getPul()<<endl;
  cout<<"data inside window "<< N_samples_per_window_<<endl;
  
}


//-------------------------------
//conver to 2D array
// compute starting point and ending point
// set number of data points inside PRI
// compute number of pulses
//---------------------------------
void RasProc::convert2Decho1PRI(const Uvar& range_start,
			    const Uvar& range_end)
  {
    if(range_end < range_start) ErrorMessage("RasProc: range end is smaller than range start").throwMe();

    //compute start index
    range_start_ = range_start;
    if(!ieb_loaded_) ErrorMessage("RasProc::ieb is not loaded").throwMe();
    Uvar dt = 2.0*range_start_/speed_light - ieb_.getRwd();
    int start_index = (int) round_double((dt*ieb_.getAdc()).getInUnits(""));
    
    //if dn is negative, discard partial pulses
    while(start_index < 0){
      cout<<"discarding one pulse of data "<<endl;
      start_index += (int) ieb_.getPointsInsidePRI();
    }
    
    //compute end index
    range_end_ = range_end;
    dt = double(ieb_.getPul()-1) * ieb_.getPri() + ieb_.getTaup() ;
    dt += 2.0 * range_end_/speed_light;
    dt -= ieb_.getRwd();
    int end_index = (int) round_double( (dt*ieb_.getAdc()).getInUnits(""));
    while(end_index > (int) ieb_.getRadarDataPoints()){
      cout<<"discarding one pulse of data from echo tail "<<endl;
      end_index -= (int) ieb_.getPointsInsidePRI();
    }
    if(end_index < start_index) ErrorMessage("no useful data inside echo").throwMe();
    
    N_pulses_inside_echo_ = ((unsigned int)(end_index-start_index)) /ieb_.getPointsInsidePRI();
    
   
    

    //---------------
    //fix number of points inside PRI as ieb.getPri()*ieb.getAdc()
    //---------------
    N_samples_per_window_ = ieb_.getPointsInsidePRI();
    if(N_samples_per_window_%2!=0) ErrorMessage("inside prf, not even number of data points ").throwMe();

    
    
    get2Decho((unsigned int) start_index,
	      (unsigned int) end_index,
	      N_samples_per_window_,
	      N_pulses_inside_echo_);//transform into 2D array
    echo2D_set_ = true;

    //debug
    
    //cout<<"start_index end_index useful_pulse comm_pulse "<< start_index<<" "<<end_index<<" "<< N_pulses_inside_echo_ <<" "<<ieb_.getPul()<<endl;
    //cout<<"data inside window "<< N_samples_per_window_<<endl;
    
  }

//---------------
//get valid number of pulses
//---------------
unsigned int RasProc::getValidNumberOfPulses()
  {
    if(!echo2D_set_) ErrorMessage("echo was not decomposed into 2d").throwMe();
    return(N_pulses_inside_echo_);
  }
//----------------------------
//range and azimuth
//----------------------------
void RasProc::computeRangeAzimuthCompression()
  {
 
  if(!ieb_loaded_) ErrorMessage("RasProc::ieb is not loaded").throwMe();
  if(!matched_filter_set_) ErrorMessage("RasProc.cpp: no matched filter available ").throwMe();
  if(!echo2D_set_) ErrorMessage("RasProc.cpp: no 2D array format").throwMe();
  Uvar prf = 1/ieb_.getPri();
 
  simpleRangeDopplerCompression(echo2D_c,
				fft_matched_filter,
				prf,
				doppler_centroid_,
				ieb_.getAdc(),
				range_start_,
				range_axis,
				doppler_axis,
				echo2D_rc,
				echo2D_rdc);

  
  //unsigned int M, N;
  //echo2D_c.size(M,N);
  //echo2D_rc.resize(M,N);
  //echo2D_rdc.resize(M,N);
  //
  //--------------------
  //range axis
  //---------------------			     
  //range_axis.resize(N);
  //for(unsigned int i=0; i<N;++i)
  //{
  //range_axis(i)=2.0*double(i)/ieb_.getAdc()*speed_light/2;
  //range_axis(i)+=range_start_;
  //}

  //----------------------------------------------
  //declare temporary variables for recycling during range compression
  //-------------------------------------------
  //CFvec x_data(" ",N);
  //CFvec fft_x_data("",N);
  //CFvec ifft_x_data("",N);
  //
  //first range compress
  //for(unsigned int ii=0; ii < M; ii++)
  //{
  //  //for each pulse
  //  x_data= echo2D_c.getRow(ii);//complex data already
  // 
  //  //fft 
  //  fft(x_data,fft_x_data);
  //
  //  // multiply conjugate of matched filter
  //  for(unsigned int jj=0;jj<N;++jj)
  //fft_x_data(jj)*= conj(fft_matched_filter(jj+N));
  //
  // 
  //  ifft(fft_x_data,ifft_x_data);
  //
  //  for(unsigned int jj=0;jj<N;++jj)
  //echo2D_rc(ii,jj)=ifft_x_data(jj);
  //}
  //
  // displayRangeCompressedData();

  //-------------------------
  // azimuth compress
  //--------------------------
  // CFvec azimuth("doppler process ",M);
  //CFvec azimuth_fft("ffted azimuth",M);
  //Uvar prf = 1/ieb_.getPri();
  //Uvar pri = ieb_.getPri();
  //
  //doppler_axis.resize(M);
  //for(unsigned int i=0;i<M;++i) doppler_axis(i)=doppler_centroid_ - prf/2.0 +prf*double(i)/double(M-1);

  //B.S. code: azimuth compression
  //complex<float> shift;
  //float pre_factor;
  //Uvar t0= ieb_.getPri()/2.0;
  //Uvar fdop_rate=Uvar(0,"Hz/s");
  //Uvar freq_rate;
  //Uvar t_midle;
  //Uvar freq_offset;
  //Uvar t_middle;
  //
  //for(unsigned int jj=0;jj<N;++jj)
  //{//for each range line
  //  azimuth=complex<float>(0,0);//reset
  //  for (unsigned int ii=0;ii<M;++ii)
  //{//for each azimuth line
  //  freq_offset= 2.0*pi*(doppler_centroid_ - prf/2.0);
  //  freq_rate = pi*fdop_rate;
  //  t_middle = t0 + ieb_.getPri()*double(ii);
  //  pre_factor= float( ( (freq_rate*t_middle + freq_offset)*t_middle).getInUnits(""));
  //  shift=exp(complex<float>(0,-1)*(float)(pre_factor));
  //  azimuth(ii) = echo2D_rc(ii,jj)*shift;
  //}
  //  //--------------
  //  //fft azimuth 
  //  //-------------
  //  fft(azimuth,azimuth_fft);
  //  for(unsigned int ii=0;ii<azimuth_fft.size();++ii)
  //echo2D_rdc(ii,jj)=azimuth_fft(ii);//save it 
  //}
  // displayAzimuthData();
  
  }


//------------------------------------------
//compute fractional doppler centroid
//------------------------------------------
double RasProc::computeFractionalDoppler()
  {
    if(!echo2D_set_) ErrorMessage("RasProc::computeFractionalDoppler: data is not in the form of 2D array").throwMe();
    unsigned int M,N;
    echo2D_c.size(M,N);
    c_corr2D.resize(M,N);
    c_corr2D=0.0;

    c_corr1D.resize(N);
    c_corr1D=0.0;

    dop_frac1D.resize(N);
    dop_frac1D=0.0;

    dop_frac0D=0.0;
    c_corr0D=0.0;

    if(valid_data_per_burst%valid_data_per_pri!=0 ) return(dop_frac0D);//bad data set
    ComputFractDopCentroid(echo2D_c, dop_frac1D, dop_frac0D);
    /*
    //compute pulse correlation
    for(unsigned int i=1; i < M; i++)
    for(unsigned int j = 0; j < N; j++)
    c_corr2D(i-1,j) = echo2D_c(i,j)*conj(echo2D_c(i-1,j));
    //
    //average pulse correlation
    for(unsigned int i = 0; i < M; i++)
    for(unsigned int j = 0; j < N; j++){
    c_corr1D(j) = c_corr1D(j) + c_corr2D(i,j)/float(M);
    c_corr0D = c_corr0D + c_corr2D(i,j)/float(N)/float(M);
    }
    //fractional doppler
    double pi = 4.*atan(1.0);
    for(unsigned int j = 0; j < N; j++)
    dop_frac1D(j) = atan2(imag(c_corr1D(j)), real(c_corr1D(j)))/2./pi;
    
    dop_frac0D =  atan2(imag(c_corr0D), real(c_corr0D))/2./pi;
    //
    */

    return(dop_frac0D);
  }

//--------------------------------------
//compute range vs fractional doppler
//----------------------------------------
void RasProc::computeRangeFractionalDoppler(Uvec& range, 
					    Dvec& fract_doppler,
					    double& mean_fract_doppler)
  {
    if(!echo2D_set_) ErrorMessage("RasProc::computeFractionalDoppler: data is not in the form of 2D array").throwMe();
    if(!ieb_loaded_) ErrorMessage("RasProc::ieb is not loaded").throwMe();
    unsigned int M, N;
    echo2D_c.size(M,N); 

    //--------------------
    //range axis
    //---------------------			     
    range_axis.resize(N);
    for(unsigned int i=0; i<N;++i){
      range_axis(i)=2.0*double(i)/ieb_.getAdc()*speed_light/2;
      range_axis(i)+=range_start_;
    }
    //compute mean value
    mean_fract_doppler=computeFractionalDoppler();
  
    //transfer range values
    range.resize(N);
    fract_doppler.resize(N);
    //transfer fractional doppler values
    range = range_axis;
    fract_doppler=dop_frac1D;
  }
//----------------------------
//display matched filter
//----------------------------
void RasProc::displayMatchedFilter()
{


  Uvec y("",matched_filter.size());
  Uvec fft_y("",fft_matched_filter.size());
  for(unsigned int i=0; i<y.size();++i)
    {
      y(i) = abs(matched_filter(i));
      fft_y(i)=abs(fft_matched_filter(i));
    }
  Plot a;
  a.addY(y,"",line("solid","black",2),sym("none"));
  a.setTitle("matched filter (time domain)");
  a.show("x");

  Plot b;
  b.addY(fft_y,"",line("solid","black",2),sym("none"));
  b.setTitle("matched filter (frequency domain)");
  b.show("x");


}

//------------------------------
//Display Range Compressed Data
//-----------------------------
void RasProc::displayRangeCompressedData()
{
 
 unsigned int M,N;
 echo2D_rc.size(M,N);//M pulse N range

 CFvec y(" ",N);
 Plot a;
 Plot b;
 for(unsigned int ii=0; ii<M;++ii)
   { 
     y = echo2D_rc.getRow(ii);
     Uvec abs_y(" ",N);
     abs_y = 0;
     for(unsigned int jj=0;jj<N;++jj)
       {
	 abs_y(jj)=abs(y(jj));
       }
     a.addY(abs_y,"",line("solid","red",2),sym("none"));
     b.addXY(range_axis,"km",abs_y,"",line("solid","red",2),sym("none"));
   }
 b.setXlabelandUnits("range");
 a.show("x");  
 b.show("x");
}


//-----------------------------
//display Azimuth compressed data
//----------------------------------
void RasProc::displayAzimuthData()
  {
    unsigned int M,N;
    echo2D_c.size(M,N);  
    Plot a;
    for(unsigned int jj=0;jj<N;++jj)
      {
	CFvec y(" ",M);
	y=echo2D_rdc.getCol(jj);
	Uvec abs_y(" ",M);
	for(unsigned int ii=0;ii<M;++ii) abs_y(ii)=abs(y(ii));
	a.addXY(doppler_axis,"KHz",abs_y,"",line("solid","red",2),sym("none"));
      }
    a.setXlabelandUnits("doppler");
    a.show("x");
  }

//------------------------
//display range and doppler image
//-------------------------
void RasProc::displayRangeDopplerImage()
  {
    unsigned int M,N;
    echo2D_rdc.size(M,N);
    Umat z("",N,M);
    
    for(unsigned int jj=0;jj<N;++jj)
      for(unsigned int ii=0;ii<M;++ii){
	z(jj,ii)=abs(echo2D_rdc(ii,jj));
      }
    //cout<<"range doppler size "<< range_axis.size()<<" "<<doppler_axis.size()<<endl;
    //cout<<"image size "<< z.rows()<<" "<<z.cols()<<endl;
    Plot b;
    b.setTool("matlab");
    b.addXYZ(range_axis,"km",doppler_axis,"KHz",z,"");
    b.setXlabelandUnits("range");
    b.setYlabelandUnits("doppler");
    b.setTitle("range doppler");
    b.show2Dimage("x"); 
  }
