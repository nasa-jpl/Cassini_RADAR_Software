#include <stdlib.h>
#include"SARFunctions.h"
#include"Constants.h"
#include "TargetGeom.h"
#include "Plot.h"
#include "Baq.h"
#include "Ieb.h"
#include "SimpleArray.h"

using std::cout;
using std::endl;



void range_compress(const fdata& x, unsigned int i, 
                    unsigned int N, const CFvec& f, 
                    CFvec& dechirped)
{
  //unsigned int M=getNextPower2(N);
  unsigned int M=N; // forces dft rather than zeropadded fft

  // append rotated signal and take fft  of  signal
  CFvec xc("xc",M);
  for(unsigned int c=0;c<M;c++){
    if(c<N) xc(c)=complex<float>(x(i+c),0);
    else xc(c)=complex<float>(0,0);
  }
  CFvec X("X",M);
  fft(xc,X);


  // zeropad and take fft of matched filter
  CFvec fc("fc",M);
  for(unsigned int c=0;c<M;c++){
    if(c<N) fc(c)=f(c);
    else fc(c)=complex<float>(0,0);
  }
  CFvec F("F",M);
  fft(fc,F);
  

 // Take complex<float> conjugate of F
  for(unsigned int i=0;i<M;i++){
    F(i)=conj(F(i));
  }

  // perform complex<float> multiply
  for(unsigned int c=0;c<M;c++){
    X(c)*=F(c);
  }

  // need to zero out DC eventually?
  

  // inverse fft
  dechirped.resize(M);
  ifft(X,dechirped);

}
// Scott Hensley algorithm for computing a real digital step chirp
// as modified by Bryan Stiles
void dig_chirp_hensley(double r_chirplen, double r_bw, double r_stepsize, 
		       double r_f0, double r_phase0, double r_sampfreq, 
		       int& i_ts, complex<float>* c_chirp){

  //--------------------------------------------------------
  //     INPUT PARAMETERS:
  //---------------------------------------------------------
  //    double r_chirplen             !chirp length in seconds
  //    double r_bw                   !bandwidth in Hz
  //    double r_stepsize             !DCG step length in seconds
  //    double r_f0                   !chirp center frequency (Hz)
  //    double r_phase0               !initial phase in chirp
  //    double r_sampfreq             !sampling frequency

  //-----------------------------------------------------------------	
  //  OUTPUT PARAMETERS:
  //-----------------------------------------------------------------
  // unsigned int i_ts                  !number of samples generated
  // complex<float*> c_chirp(*)         !complex chirp sampled at ADC rate
  //
  // Array is allocated outside of function long enough to contain the
  // i_ts samples plus whatever zero-padding may be required outside of
  // this function

  //-------------------------------------------------------------------
  //    LOCAL VARIABLES:
  //-------------------------------------------------------------------
  int i_numsteps,i_stepnum;
  int i_totsamps;
  double r_samptime,r_dfdt,r_delfreq,r_tottime,r_phasec,r_phase;
  double r_t,r_phase0p,r_ss,r_freq,r_freq_low,r_freq_high;



  //     parameters describing number of steps in DCG chirp

  i_numsteps = int(r_chirplen/fabs(r_stepsize)+0.5);
  r_ss = r_stepsize/fabs(r_stepsize);

  //    ADC sample time and total length of chirp after ADC quantization

  r_samptime = 1.0/r_sampfreq;
  i_totsamps = int(r_chirplen*r_sampfreq+0.5);
  i_ts = (unsigned int) i_totsamps;

  if(i_numsteps > 1){                         // Using DCG
    r_delfreq = r_ss*r_bw/(i_numsteps - 1);   // frequency step size in DCG
    r_dfdt = 0.0;                             // no chirpo rate
  }
  else{                                        // ideal chirp
    r_delfreq = 0.0;                           // frequency step size is zero 
    r_dfdt = r_ss*r_bw/(i_totsamps*r_samptime);// chirp rate
    r_ss = 0.0;
  }
      
  r_tottime = i_totsamps*r_samptime;    // total chirp time after quantization

  // cancel initial phase variations to guarantee desired initial chirp phase

  r_t = 0.0;
  i_stepnum = int((r_t + 0.5*r_tottime)/fabs(r_stepsize));
  r_phasec = -pi*r_delfreq*fabs(r_stepsize)*
    (i_stepnum)*(i_stepnum - i_numsteps+1);  
  r_phase0p = r_phase0 - r_phasec;

  // set the upper and lower frequency regions to avoid spectrum aliasing

  if(r_sampfreq  <  2.0*r_bw) {      // If chirp is too wide for Nyquist
    ErrorMessage e("dig_chirp_hensley: Chirp bandwidth violates Nyquist");
    e.throwMe();
  }
  else{                               
    if(r_f0 > 0){
      r_freq_low = 0.01*(r_sampfreq/2.)  + MAX((r_bw/2 - r_f0),0);
      r_freq_high = 0.99*(r_sampfreq/2.) 
	- MAX((r_f0 + r_bw/2 - r_sampfreq/2),0.); 
    }
    else{
      r_freq_low = -0.99*(r_sampfreq/2.0) 
	+ MAX((fabs(r_f0) + r_bw/2 - r_sampfreq/2.0),0.0);
      r_freq_high = -0.01*(r_sampfreq/2.0) - MAX((r_bw/2 - fabs(r_f0)),0.0);
    }
  }

//-----------------------------------------------------------
//    generate the chirp
//-----------------------------------------------------------
 for(int i=0;i<i_totsamps;i++){
        
   r_t = i*r_samptime - r_tottime/2.0;   // time of sample 
                                         // filter is 0 centered


   i_stepnum = int((r_t + 0.5*r_tottime)/fabs(r_stepsize)); // chirp step num

   r_phasec = -pi*r_delfreq*fabs(r_stepsize)*(i_stepnum)*
     (i_stepnum - i_numsteps + 1);  

   r_phase = 2.0*pi*
     (r_f0 + i_stepnum*r_delfreq - r_ss*r_bw/2 + 0.5*r_dfdt*r_t)*r_t 
     + r_phasec + r_phase0p;

   r_freq = r_f0 + i_stepnum*r_delfreq - r_ss*r_bw/2 + r_dfdt*r_t;

   // if the freq exceeds the sampling frequency then zero out to 
   // avoid aliasing the spectrum

   if(r_freq >= r_freq_low && r_freq <= r_freq_high){
     float r_chirp = cos(r_phase);          // real chirped signal
     c_chirp[i] = complex<float>(r_chirp,0.0); // complex version (imag=0)
   }         
   else{
     c_chirp[i] = complex<float>(0.0,0.0);
   }         

 } // end for
} // end dig_chirp_hensley

void dft(const CFvec& s, CFvec& S, unsigned int N){
  for(unsigned int k=0;k<N;k++){
    S(k)=0;
    for(unsigned int n=0;n<N;n++){
      S(k)+=s(n)*exp(complex<float>(0,-1)*(float)(2.0*pi*k*n/N));
    }
  }
}

void idft(const CFvec& S, CFvec& s, unsigned int N){
  for(unsigned int n=0;n<N;n++){
    s(n)=0;
    for(unsigned int k=0;k<N;k++){
      s(n)+=S(k)*exp(complex<float>(0,1)*(float)(2.0*pi*k*n/N));
    }
    s(n)=s(n)/(float)N;
  }
}

// Routine for computing echo return from a point target
// Pulse windowing and bandpass filtering
// is implemented outside of this routine
float
get_point_target_echo(double tp, double t, 
double updown_freq_shift, double carrier_freq, 
double delay, double chirp_start_freq,
double chirp_rate, double mag, double phi0 ){
  
  // tp = time in pulse; start of pulse return time =0
  // t= echo return time since start of burst; transmit_start_time=0;
  // updown_freq_shift; residual frequency left over after down conversion 
  //                    from carrier frequency
  // lambda  is carrier frequency
  // fdop is doppler shift
  // delay is round trip time
  // chirp_start_freq is the starting frequency of the chirp at tp=0
  // chirp_rate
  // mag magnitude of signal
  // phi0 phase offset equals sum of chirp phase at tp=0, 
  //           down convert reference signal phase at t=0;
  //           and carrier phase at t=0;


  double phi; 
  phi= carrier_freq*(t-delay); // phase of carrier signal from start of burst
                          // to transmit time
  phi+= chirp_start_freq*tp+chirp_rate/2*tp*tp; // add in phase due to linear
                                                // FM chirp
  phi+= (updown_freq_shift-carrier_freq)*t; // add phase at time t due to updown

                                      // convert signal
  float retval= mag*cos(phi*2*pi+phi0);
  return(retval);
}

float linear_fm_chirp(const Uvar& t, const Uvar& f0, const Uvar& cr, float mag, const Uvar& phi0)
{
  Uvar w0=Uvar(2*pi,"rad")*f0;
  Uvar a=Uvar(2*pi,"rad")*cr;
  float retval=mag*cos(w0*t+(a/2)*t*t+phi0);
  return(retval);
}

complex<float> complex_linear_fm_chirp(const Uvar& t, const Uvar& f0, const Uvar& cr, float mag, const Uvar& phi0)
{
  Uvar w0=Uvar(2*pi,"rad")*f0;
  Uvar a=Uvar(2*pi,"rad")*cr;
  double phase = (w0*t+(a/2)*t*t+phi0).getValue();
  complex<float> retval=mag*exp(complex<float>(0,1)*(float)phase);
  return(retval);
}



float anti_aliased_linear_fm_chirp(const Uvar& t, const Uvar& f0, 
				   const Uvar& cr, float mag, 
				   const Uvar& phi0, 
				   const Uvar& pass_band_start,
				   const Uvar& pass_band_end)
{
  if(mag==0) return(0);
  Uvar f_inst=f0+t*cr;
  if(f_inst<pass_band_start || f_inst > pass_band_end){
    return(0.0);
  }
  else return(linear_fm_chirp(t,f0,cr,mag,phi0));
}


complex<float> complex_anti_aliased_linear_fm_chirp(
				   const Uvar& t, const Uvar& f0, 
				   const Uvar& cr, float mag, 
				   const Uvar& phi0, 
				   const Uvar& pass_band_start,
				   const Uvar& pass_band_end)
{
  if(mag==0) return(0);
  Uvar f_inst=f0+t*cr;
  if(f_inst<pass_band_start || f_inst > pass_band_end){

    cout <<"LINEAR CHIRP OUTSIDE OF BAND!!!!!! at freq" <<  
      f_inst << "at time " << t << endl;
    cout <<"LINEAR CHIRP_OUTSIDE_OF BAND!!!!! for analog filter (" 
	 << pass_band_start << " to " << pass_band_end << ")" << endl;
    cout <<"LINEAR CHIRP_OUTSIDE_OF BAND!!!!! start_freq=" << f0
	 << " chirp rate=" << cr << endl;
    return(0.0);
  }
  else return(complex_linear_fm_chirp(t,f0,cr,mag,phi0));
}




//----------------------------------------------
// Below are some fft routines
// All of these routines must operate on arrays with a length
// that is a power of two! Otherwise they fail
void fft(const CFvec& f, CFvec& F){
  int N=f.size();
  if(N!=(int)F.size()){
    ErrorMessage e("fft: Cfvec's f and F must be the same size");
    e.throwMe();
  }
  dieIfNotPowerofTwo(N);
  complex<float> *f_array= new complex<float>[N];
  complex<float> *F_array= new complex<float>[N];
  for(int c=0;c<N;c++) f_array[c]=f(c);
  fft(f_array,F_array,N);
  for(int c=0;c<N;c++) F(c)=F_array[c];
  delete[] f_array;
  delete[] F_array;
}
void ifft(const CFvec& F, CFvec& f){
  int N=f.size();
  if(N!=(int)F.size()){
    ErrorMessage e("fft: Cfvec's f and F must be the same size");
    e.throwMe();
  }
  dieIfNotPowerofTwo(N);
  complex<float> *f_array= new complex<float>[N];
  complex<float> *F_array= new complex<float>[N];
  for(int c=0;c<N;c++) F_array[c]=F(c);
  ifft(F_array,f_array,N);
  for(int c=0;c<N;c++) f(c)=f_array[c];
  delete[] f_array;
  delete[] F_array;
}
// N must be a power of two
void fft(complex<float>*f, complex<float>*F,int N){
  dieIfNotPowerofTwo(N);
  for(int i=0;i<N;i++){
    F[i]=f[i];
  }
  fft_base((float*)F,N,-1);
}

// N must be a power of two
void ifft(complex<float>*F, complex<float>*f,int N){
  dieIfNotPowerofTwo(N);
  for(int i=0;i<N;i++){
    f[i]=F[i]/(float)N;
  }
  fft_base((float*)f,N,1);
}

// NN must be apower of two; Routine DOES NOT CHECK
// Do not use this routine use the fft or ifft wrappers
void fft_base(float* DATA, int NN, int ISIGN){

  double WR,WI,WPR,WPI,WTEMP,THETA;
  int  N=2*NN;
  int  J=0;
  for(int I=0;I<N;I+=2){
    if(J > I)
      {
	float TEMPR=DATA[J];
	float TEMPI=DATA[J+1];
	DATA[J]=DATA[I];
	DATA[J+1]=DATA[I+1];
	DATA[I]=TEMPR;
	DATA[I+1]=TEMPI;
      }
    int M=N/2;
    while((M >= 2) && (J >= M))
      {
	J=J-M;
	M=M/2;
      }
    J=J+M;
  }
  int MMAX=2;


   while (N > MMAX) {
     int ISTEP=2*MMAX;
     THETA=6.283185307179590/(ISIGN*MMAX);
     WPR=sin(0.5*THETA);
     WPR=-2*WPR*WPR;
     WPI=sin(THETA);
     WR=1.0;
     WI=0.0;
     for(int M=0;M<MMAX;M+=2){
       for(int I=M;I<N;I+=ISTEP){
	 J=I+MMAX;
	 float TEMPR=(float)WR*DATA[J]-(float)WI*DATA[J+1];
	 float TEMPI=float(WR)*DATA[J+1]+float(WI)*DATA[J];
	 DATA[J]=DATA[I]-TEMPR;
	 DATA[J+1]=DATA[I+1]-TEMPI;
	 DATA[I]=DATA[I]+TEMPR;
	 DATA[I+1]=DATA[I+1]+TEMPI;
       }
       WTEMP=WR;
       WR=WR*WPR-WI*WPI+WR;
       WI=WI*WPR+WTEMP*WPI+WI;
     }
     MMAX=ISTEP;
   }
}

void dieIfNotPowerofTwo(int N){
  // first make sure it is at least 2
  if(N>=2){
    // then see if it is a power of two by adding all the bits 
    // sum should be one.
    int sum=0;
    for(int c=0;c<31;c++){
      sum+= (0x01 & N);
      N=N >> 1;
    }
    if(sum==1) return; // It is a positive power of 2. 
  }
  ErrorMessage e("FFT needs N= a positive power of two >= 2");
  e.throwMe();
}

//--------------------------------------
//chirp generation
//--------------------------------------
float linear_chirp(const double& slow_time, 
		   const double& t, 
		   const double& fast_time,
		   const double& chirp_start_frequency_in_Hz, 
		   const double& chirp_rate,
		   const double& cfs_in_Hz,
		   const double& phase)
  {
  double fc = 13.77e9;//13.77 GHz Ku band
  double slo = 10e6;//10MHz stable oscillator
  double phi1 = fc*(t-slow_time);//time delay
  phi1 +=  chirp_start_frequency_in_Hz*fast_time 
    +chirp_rate/2.0*fast_time*fast_time;//linear chirp
  //----------------------------------------------------------
  //correction factor to take into account
  //the difference between digital chirp and linear chirp  
  //-----------------------------------------------------------
  phi1 -= cfs_in_Hz*fast_time/2.0;
  //down conversion   
  phi1 -= (fc+slo)*t;
  float value= cos(2.0*pi*(phi1+phase));
  return(value);
  }

//-------------------
//chirp generation
//------------------
complex<float> complex_linear_chirp(const double& slow_time, 
		    const double& t, 
		    const double& fast_time,
		    const double& chirp_start_frequency_in_Hz, 
		    const double& chirp_rate,
		    const double& cfs_in_Hz)
  {
  double fc = 13.77e9;//13.77 GHz Ku band
  double slo = 10e6;//10MHz stable oscillator
  double phi1 = fc*(t-slow_time);//time delay
  phi1 +=  chirp_start_frequency_in_Hz*fast_time 
    +chirp_rate/2.0*fast_time*fast_time;//linear chirp
  //----------------------------------------------------------
  //correction factor to take into account
  //the difference between digital chirp and linear chirp  
  //-----------------------------------------------------------
  phi1 -= cfs_in_Hz*fast_time/2.0;
  //down conversion   
  phi1 -= (fc+slo)*t;
  phi1 *= 2.0* pi;
  complex<float> value= exp(complex<float>(0,1.0)*(float)phi1);
  return(value);
  }

//----------------
//digital chirp
//-----------------

float digital_chirp(const double& slow_time,
		     const double& t, 
		     const double& fast_time,
		     const double& chirp_start_frequency_in_Hz, 
		     const int& nstep, 
		     const double& csd_in_s, 
		     const double& cfs_in_Hz)
  {
  double fc = 13.77e9;//13.77 GHz Ku band
  double slo = 10e6;//10MHz stable oscillator
  double phi1 = fc*(t-slow_time);//time delay
  phi1 += (chirp_start_frequency_in_Hz+ double(nstep)*cfs_in_Hz) *fast_time;
  phi1 -= double(nstep*(nstep+1))/2.0*csd_in_s * cfs_in_Hz;//step chirp
  phi1 -= (fc+slo)*t;//down conversion
  float value = cos(2.0*pi*phi1);
  return(value);
  }

//---------------
//complex digital chirp
//----------------
complex<float> complex_digital_chirp(const double slow_time,
				    const double t, 
				    const double fast_time,
				    const double chirp_start_frequency_in_Hz, 
				    const int nstep, 
				    const double csd_in_s, 
				    const double cfs_in_Hz)
  {
  double fc = 13.77e9;//13.77 GHz Ku band
  double slo = 10e6;//10MHz stable oscillator
  double phi1 = fc*(t-slow_time);//time delay
  phi1 += (chirp_start_frequency_in_Hz+ double(nstep)*cfs_in_Hz) *fast_time;
  phi1 -= double(nstep*(nstep+1))/2.0*csd_in_s * cfs_in_Hz;//step chirp
  phi1 -= (fc+slo)*t;//down conversion
  phi1 *= 2.0*pi;
  complex<float> value=exp(complex<float>(0,1.0)*(float)phi1); 
  return(value);
  }

//--------------------
//echo for re-routed chirp
//-------------------
void generate_rerouted_chirp(Ieb& ieb,
			     const double& phase,
			     fdata& echo){
  //--------------------------
  //Get ieb at burst
  //--------------------------
  double pri = ieb.getPri().getInUnits("s");
  double taup = ieb.getTaup().getInUnits("s");
  double chirp_start_frequency_in_Hz = ieb.getChirpStartFrequency().getInUnits("Hz");
  double cfs_in_Hz = ieb.getCfs().getInUnits("Hz");
  double csd_in_s = ieb.getCsd().getInUnits("s");
  //unsigned int csq = ieb.getCsq();
  double chirp_rate =cfs_in_Hz/csd_in_s;
  unsigned int Npulse = ieb.getPul();
  int Ntro = ieb.getTroInPriUnits();
  unsigned int Npoint_per_pri = ieb.getPointsInsidePRI();
  int rcv_in_pri_unit = Ntro + int(Npulse);
  if (rcv_in_pri_unit <= 0) ErrorMessage("zero receive window").throwMe();
  unsigned int N = ieb.getRadarDataPoints();//total number of data points/bpd
  if(N != echo.size()) ErrorMessage("SARFuntions.cpp:echo container size does not match to the data size").throwMe();
  echo=0.0;//reset the container
  double sample_rate = ieb.getAdc().getInUnits("Hz");
  double rcv_start = ieb.getRwd().getInUnits("s");
  double rcv_end = rcv_start + double(Npulse)*pri + double(Ntro)*pri;
  
  //---------------
  //calculate echo
  //--------------
  double delay_time;//delay time for each pulse
  double time_at_rcv;//rwd + index/sample_rate
  unsigned int  index_at_rcv;//index number 
  double time_at_transmitt;//time at transmitt
  
  double fast_time;
  //double pos_X,pos_Y,pos_Z;
  //double distance;
  //double r_dot_v;
  //double v_square;
  unsigned int pulse_id;

  //unsigned int Necho_sample_ratio = 10;
  //10 MHz sampling
  double echo_sample_rate =10*1e6;
  unsigned int Necho_sample_ratio = (unsigned int) round_double(echo_sample_rate/sample_rate);


  //tmp echo holder
  float* echo_tmp;
  echo_tmp=(float*)make_array(sizeof(float),1,Necho_sample_ratio*N);
  for(unsigned int i=0;i<Necho_sample_ratio*N;++i) echo_tmp[i]=0.0;

  for (unsigned int i=0;i<Npulse*Npoint_per_pri*Necho_sample_ratio;++i)
    {
      pulse_id = (unsigned int)(i/Npoint_per_pri/Necho_sample_ratio);
      time_at_transmitt= double(i)/echo_sample_rate;
      fast_time = time_at_transmitt - double(pulse_id)*pri;
      if(fast_time < taup)
	{
	  delay_time=0;
	  time_at_rcv = time_at_transmitt+delay_time;
	  if( time_at_rcv>rcv_start &&
	      time_at_rcv<rcv_end)
	    {
	      index_at_rcv = (unsigned int) round_double((time_at_rcv - rcv_start)*echo_sample_rate);
	      if(index_at_rcv < Necho_sample_ratio*N){
		echo_tmp[index_at_rcv]= linear_chirp(delay_time, 
						     time_at_rcv, 
						     fast_time,
						     chirp_start_frequency_in_Hz,
						     chirp_rate,
						     cfs_in_Hz,
						     phase);
	      }
	      else 
		cout<<"this point will be thrown out , outside of window "<<endl;
	      
	    }//inside receive window
	}//fast time < taup
    }//high frequency sampling
  
  //---------------------------------------------------------------------
  //from oversampled data, fill up the understample data: using average
  //-------------------------------------------------------------------
 
  //--------------------------------------
  //from oversampled data, fill up the undersampled data: no average
  //
  //-------------------------------------
  unsigned int index;
  for (unsigned int i=0;i<N;++i){
    index = i *Necho_sample_ratio;
    echo(i)=echo_tmp[index];
  }
  
  free_array((void*)echo_tmp,1,Necho_sample_ratio*N);
  
}


//-----------
//echo simulation
//-----------
void generate_echo_using_LFM(Ieb& ieb, 
			     const string& target,
			     const Frame& target_frame,
			     const Uvar& lat,
			     const Uvar& lon,
			     const Time& burst_time,
			     const double& phase,
			     fdata& echo)
			     
  {
  Uvar transmit_time_offset=Uvar(1*mstos,"s")+ieb.getPri();  
  echo=0.0; 
  double light_speed = speed_light.getInUnits("m/s");
  StateVector sc_state;
  target_frame.ephemeris(sc_state,"Cassini",burst_time+transmit_time_offset,"NONE");
  double sc_P_x, sc_P_y,sc_P_z;
  double sc_V_x, sc_V_y, sc_V_z;
  sc_P_x = sc_state.position()[PositionVector::X].getInUnits("m");
  sc_P_y = sc_state.position()[PositionVector::Y].getInUnits("m");
  sc_P_z = sc_state.position()[PositionVector::Z].getInUnits("m");
  
  sc_V_x = sc_state.velocity()[FloatVector::X].getInUnits("m/s");
  sc_V_y = sc_state.velocity()[FloatVector::Y].getInUnits("m/s");
  sc_V_z = sc_state.velocity()[FloatVector::Z].getInUnits("m/s");
  
  TargetGeom tg(burst_time+transmit_time_offset);
  tg.setState(sc_state);
  tg.setTarget(target,target_frame);
  tg.setLatLon(lat,lon);
  if(!tg.foundSurfaceIntercept()) return;//no phase added
    
  double X,Y,Z;
  X= tg.surfaceIntercept()[PositionVector::X].getInUnits("m");
  Y =tg.surfaceIntercept()[PositionVector::Y].getInUnits("m");
  Z= tg.surfaceIntercept()[PositionVector::Z].getInUnits("m");
  
  //--------------------------
  //Get ieb at burst
  //--------------------------
  double pri = ieb.getPri().getInUnits("s");
  double taup = ieb.getTaup().getInUnits("s");
  double chirp_start_frequency_in_Hz = ieb.getChirpStartFrequency().getInUnits("Hz");
  double cfs_in_Hz = ieb.getCfs().getInUnits("Hz");
  double csd_in_s = ieb.getCsd().getInUnits("s");
  //unsigned int csq = ieb.getCsq();
  double chirp_rate =cfs_in_Hz/csd_in_s;
  unsigned int Npulse = ieb.getPul();
  int Ntro = ieb.getTroInPriUnits();
  unsigned int Npoint_per_pri = ieb.getPointsInsidePRI();
  int rcv_in_pri_unit = Ntro + int(Npulse);
  if (rcv_in_pri_unit <= 0) ErrorMessage("zero receive window").throwMe();
  unsigned int N = ieb.getRadarDataPoints();//total number of data points/bpd
  if(N != echo.size()) ErrorMessage("SARFuntions.cpp:echo container size does not match to the data size").throwMe();
  echo=0.0;//reset the container
  double sample_rate = ieb.getAdc().getInUnits("Hz");
  double rcv_start = ieb.getRwd().getInUnits("s");
  double rcv_end = rcv_start + double(Npulse)*pri + double(Ntro)*pri;
  
  //---------------
  //calculate echo
  //--------------
  double delay_time;//delay time for each pulse
  double time_at_rcv;//rwd + index/sample_rate
  unsigned int  index_at_rcv;//index number 
  double time_at_transmitt;//time at transmitt
  
  double fast_time;
  double pos_X,pos_Y,pos_Z;
  double distance;
  double r_dot_v;
  double v_square;
  unsigned int pulse_id;

  v_square= sc_V_x * sc_V_x + sc_V_y*sc_V_y +  sc_V_z*sc_V_z;

  //unsigned int Necho_sample_ratio = 10;
  //10 MHz sampling
  double echo_sample_rate =10*1e6;
  unsigned int Necho_sample_ratio = (unsigned int) round_double(echo_sample_rate/sample_rate);
  
 
 
  float* echo_tmp;
  echo_tmp=(float*)make_array(sizeof(float),1,Necho_sample_ratio*N);
  for(unsigned int i=0;i<Necho_sample_ratio*N;++i) echo_tmp[i]=0.0;


  for (unsigned int i=0;i<Npulse*Npoint_per_pri*Necho_sample_ratio;++i)
    {
      pulse_id = (unsigned int)(i/Npoint_per_pri/Necho_sample_ratio);
      time_at_transmitt= double(i)/echo_sample_rate;
      fast_time = time_at_transmitt - double(pulse_id)*pri;
      if(fast_time < taup)
	{
	pos_X = sc_P_x +sc_V_x*time_at_transmitt;
	pos_Y = sc_P_y +sc_V_y*time_at_transmitt;
	pos_Z = sc_P_z +sc_V_z*time_at_transmitt;
      
	distance = sqrt(pow(pos_X-X,2)+pow(pos_Y-Y,2)+pow(pos_Z-Z,2));
	r_dot_v = (pos_X-X)*sc_V_x + (pos_Y-Y)*sc_V_y+(pos_Z-Z)*sc_V_z;
	
	delay_time =2.0*(distance*light_speed + r_dot_v);//+ for forward
	delay_time /=(light_speed*light_speed-v_square);
	time_at_rcv = time_at_transmitt+delay_time;
	if( time_at_rcv>rcv_start &&
	    time_at_rcv<rcv_end)
	  {
	  index_at_rcv = (unsigned int) round_double((time_at_rcv - rcv_start)*echo_sample_rate);
	  if(index_at_rcv < Necho_sample_ratio*N){
	    echo_tmp[index_at_rcv]= linear_chirp(delay_time, 
						 time_at_rcv, 
						 fast_time,
						 chirp_start_frequency_in_Hz,
						 chirp_rate,
						 cfs_in_Hz,
						 phase);
	  }
	  else 
	    cout<<"this point will be thrown out , outside of window "<<endl;
	  
	  }//inside receive window
	}//fast time < taup
    }//high frequency sampling
  
  //---------------------------------------------------------------------
  //from oversampled data, fill up the understample data: using average
  //-------------------------------------------------------------------
 
  //--------------------------------------
  //from oversampled data, fill up the undersampled data: no average
  //
  //-------------------------------------
  unsigned int index;
  for (unsigned int i=0;i<N;++i){
    index = i *Necho_sample_ratio;
    echo(i)=echo_tmp[index];
  }
 
  free_array((void*)echo_tmp,1,Necho_sample_ratio*N);
  
  
  //debugging
  //Plot a;
  //Uvec y("",echo_tmp.size());
  //y=echo_tmp;
  //a.addY(y,"",line("solid","red",1),sym("none"));
  //a.show("x");
  //Plot b;
  //Uvec y2("",echo.size());
  //y2=echo;
  //b.addY(y2,"",line("solid","red",1),sym("none"));
  //b.show("x");
  }

//-----------
//echo generation when bif != 1
//---------
void generate_echo_using_LFM_BIF(Ieb& ieb, 
				 const string& target,
				 const Frame& target_frame,
				 const Uvar& lat,
				 const Uvar& lon,
				 const Time& burst_time,
				 const double& phase,
				 fdata& echo,
				 const unsigned int& num_bursts_in_flight)
  {//forward solution

  //when num_bursts_in_flight is not zero, we need to compute
  // state backward
  if(num_bursts_in_flight<1) 
    ErrorMessage("sARFuntions.cpp: num bursts in flight should be larger than 1").throwMe();
  int burst_start_offset_in_units_of_bpd= num_bursts_in_flight -1;
  //--------------------------
  //Get ieb at burst
  //--------------------------
  double pri = ieb.getPri().getInUnits("s");
  double taup = ieb.getTaup().getInUnits("s");
  double chirp_start_frequency_in_Hz = ieb.getChirpStartFrequency().getInUnits("Hz");
  double cfs_in_Hz = ieb.getCfs().getInUnits("Hz");
  double csd_in_s = ieb.getCsd().getInUnits("s");
  //unsigned int csq = ieb.getCsq();
  double chirp_rate =cfs_in_Hz/csd_in_s;
  unsigned int Npulse = ieb.getPul();
  int Ntro = ieb.getTroInPriUnits();
  unsigned int Npoint_per_pri = ieb.getPointsInsidePRI();
  int rcv_in_pri_unit = Ntro + int(Npulse);
  if (rcv_in_pri_unit <= 0) ErrorMessage("zero receive window").throwMe();
  unsigned int N = ieb.getRadarDataPoints();//total number of data points/bpd
  if(N != echo.size()) ErrorMessage("SARFuntions.cpp:echo container size does not match to the data size").throwMe();
  echo=0.0;//reset the container
  double sample_rate = ieb.getAdc().getInUnits("Hz");
  double rcv_start = ieb.getRwd().getInUnits("s");
  double rcv_end = rcv_start + double(Npulse)*pri + double(Ntro)*pri;

  echo=0.0; 
  double light_speed = speed_light.getInUnits("m/s");
  


  //---------------
  //get state vector
  //-------------
  Uvar burst_start_time =  -ieb.getBpd()*(double) burst_start_offset_in_units_of_bpd;
  double d_burst_start_time =burst_start_time.getInUnits("s");


  //debug
  //cout<<"burst time "<< burst_start_time<<endl;
  Uvar transmit_time_offset=Uvar(1*mstos,"s")+ieb.getPri();
 

  Time offset_burst_time = burst_time + burst_start_time+transmit_time_offset;
  StateVector sc_state;
  target_frame.ephemeris(sc_state,"Cassini",offset_burst_time,"NONE");
  double sc_P_x, sc_P_y,sc_P_z;
  double sc_V_x, sc_V_y, sc_V_z;
  sc_P_x = sc_state.position()[PositionVector::X].getInUnits("m");
  sc_P_y = sc_state.position()[PositionVector::Y].getInUnits("m");
  sc_P_z = sc_state.position()[PositionVector::Z].getInUnits("m");
  
  sc_V_x = sc_state.velocity()[FloatVector::X].getInUnits("m/s");
  sc_V_y = sc_state.velocity()[FloatVector::Y].getInUnits("m/s");
  sc_V_z = sc_state.velocity()[FloatVector::Z].getInUnits("m/s");
  
  TargetGeom tg(offset_burst_time);
  tg.setState(sc_state);
  tg.setTarget(target,target_frame);
  tg.setLatLon(lat,lon);
  if(!tg.foundSurfaceIntercept()) return;//no phase added
    
  double X,Y,Z;
  X= tg.surfaceIntercept()[PositionVector::X].getInUnits("m");
  Y =tg.surfaceIntercept()[PositionVector::Y].getInUnits("m");
  Z= tg.surfaceIntercept()[PositionVector::Z].getInUnits("m");
  
 
  
  //---------------
  //calculate echo
  //--------------
  double delay_time;//delay time for each pulse
  double time_at_rcv;//rwd + index/sample_rate
  unsigned int  index_at_rcv;//index number 
  double time_at_transmitt;//time at transmitt
  
  double fast_time;
  double pos_X,pos_Y,pos_Z;
  double distance;
  double r_dot_v;
  double v_square;

  unsigned int pulse_id;

  v_square= sc_V_x * sc_V_x + sc_V_y*sc_V_y +  sc_V_z*sc_V_z;

  //unsigned int Necho_sample_ratio = 10;
  //10 MHz sampling
  double echo_sample_rate =10*1e6;
  unsigned int Necho_sample_ratio = (unsigned int) round_double(echo_sample_rate/sample_rate);

  float* echo_tmp;
  echo_tmp=(float*)make_array(sizeof(float),1,Necho_sample_ratio*N);
  for(unsigned int i=0;i<Necho_sample_ratio*N;++i) echo_tmp[i]=0.0;
  
 

  for (unsigned int i=0;i<Npulse*Npoint_per_pri*Necho_sample_ratio;++i)
    {
      pulse_id = (unsigned int)(i/Npoint_per_pri/Necho_sample_ratio);
      time_at_transmitt= double(i)/echo_sample_rate;
      fast_time = time_at_transmitt - double(pulse_id)*pri;
      if(fast_time < taup)
	{
	pos_X = sc_P_x +sc_V_x*time_at_transmitt;
	pos_Y = sc_P_y +sc_V_y*time_at_transmitt;
	pos_Z = sc_P_z +sc_V_z*time_at_transmitt;
	//remember: sc_p_x was computed -bpd+transmit_time_offset
      
	distance = sqrt(pow(pos_X-X,2)+pow(pos_Y-Y,2)+pow(pos_Z-Z,2));
	r_dot_v = (pos_X-X)*sc_V_x + (pos_Y-Y)*sc_V_y+(pos_Z-Z)*sc_V_z;
	
	delay_time =2.0*(distance*light_speed +r_dot_v);//+ for forward
	delay_time /=(light_speed*light_speed-v_square);
	time_at_rcv = time_at_transmitt+delay_time + d_burst_start_time;//BIF effect
	if( time_at_rcv>rcv_start &&
	    time_at_rcv<rcv_end)  {
	  index_at_rcv = (unsigned int) round_double((time_at_rcv - rcv_start)*echo_sample_rate);
	  if(index_at_rcv < Necho_sample_ratio*N){
	    echo_tmp[index_at_rcv]= linear_chirp(delay_time, 
						 time_at_rcv , 
						 fast_time,
						 chirp_start_frequency_in_Hz,
						 chirp_rate,
						 cfs_in_Hz,
						 phase);
	  }
	  else 
	    cout<<"this point will be thrown out , outside of window "<<endl;
	  
	  }//inside receive window
	}//fast time < taup
    }//high frequency sampling
  
  //---------------------------------------------------------------------
  //from oversampled data, fill up the understample data: using average
  //-------------------------------------------------------------------
 
  //--------------------------------------
  //from oversampled data, fill up the undersampled data: no average
  //
  //-------------------------------------
  unsigned int index;
  for (unsigned int i=0;i<N;++i){
    index = i *Necho_sample_ratio;
    echo(i)=echo_tmp[index];
  }
 
  
  free_array((void*)echo_tmp,1,Necho_sample_ratio*N);  
  //debugging
  //Plot a;
  //Uvec y("",echo_tmp.size());
  //y=echo_tmp;
  //a.addY(y,"",line("solid","red",1),sym("none"));
  //a.show("x");
  //Plot b;
  //Uvec y2("",echo.size());
  //y2=echo;
  //b.addY(y2,"",line("solid","red",1),sym("none"));
  //b.show("x");
  }


//--------------------------------
//compute matched filter 
//----------------------------------
void compute_matched_filter(const double& adc,
			    const double& taup,
			    const double& chirp_rate,
			    const double& cfs,
			    const double& filter_start,
			    const double& filter_end,
			    const double& start_frequency,
			    const unsigned int& N,
			    unsigned int& i_npts,
			    CFvec& matched_filter)
  {
    if(matched_filter.size()!=N)
      matched_filter.resize(N);
    i_npts=0;
    double t, f, phi;
    for (unsigned int i =0; i <N;++i){
      t = double(i)/adc;
      f= start_frequency+chirp_rate*t;
      phi=0.0;
      if(t<taup){
	if(f >= filter_start && f<=filter_end){
	  phi= start_frequency*t ;
	  phi+= chirp_rate/2.0*t*t;//linear chirp
	  //----------------------------------------------------------
	  //correction factor to take into account
	  //the difference between digital chirp and linear chirp  
	  //-----------------------------------------------------------
	  phi -= cfs*t/2.0;
	  phi *= 2.0* pi;
	  matched_filter(i)=exp(complex<float>(0,1.0)*(float)phi);
	  i_npts++;
	  //cout<<"time and phase "<<i<<" "<<phi<<" "<<matched_filter(i)<<endl;
	}
	else{
	  cout<<"Warning: compute_matched_filter: frequency is out of receiver window "<<endl;
	  matched_filter(i) = 0.0;
	}	  
      }
      else matched_filter(i) = 0.0;
    }    
  }

//-----------------------
//Wrapper function for computing matched filter
//-----------------------
void compute_matched_filter(Ieb& ieb,
			    const Uvar& doppler,
			    const Uvar& filter_start,
			    const Uvar& filter_end,
			    const unsigned int& N,
			    unsigned int& i_npts,
			    CFvec& matched_filter)
  
  {
    double adc=ieb.getAdc().getInUnits("Hz");
    double taup=ieb.getTaup().getInUnits("s");
    double chirp_rate = (ieb.getCfs()/ieb.getCsd()).getInUnits("1/ s s");
    double cfs=ieb.getCfs().getInUnits("Hz");
    double start_frequency= -slo_frequency.getInUnits("Hz")
      +ieb.getChirpStartFrequency().getInUnits("Hz")
      +doppler.getInUnits("Hz");
    double f1=filter_start.getInUnits("Hz");
    double f2=filter_end.getInUnits("Hz");

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
			    
 

  }


void compute_matched_filter_dig( double r_chirplen,
				 double r_bw,
				 double r_stepsize,
				 double r_phase0,
				 double r_sampfreq,
				 double r_f0,
				 int& i_ts,
				 const int& N,
				 CFvec& matched_filter ){

  matched_filter.resize(N);
   
  //double r_chirplen=ieb.getTaup().getInUnits("s");
  //double r_bw= ieb.getCfs().getInUnits("Hz")*ieb.getCsq();
  //double r_stepsize=ieb.getCsd().getInUnits("s");
  //double r_phase0=0;
  //double r_sampfreq= ieb.getAdc().getInUnits("Hz");
  //double r_f0= -10e6 + ieb.getChirpStartFrequency().getInUnits("Hz") + r_bw/2 + doppler.getInUnits("Hz");
  //int i_ts= i_npts;
  
  //double chirp_rate = (ieb.getCfs()/ieb.getCsd()).getInUnits("1/ s s");
  //double cfs=ieb.getCfs().getInUnits("Hz");
  //double start_frequency= -slo_frequency.getInUnits("Hz")
  //+ieb.getChirpStartFrequency().getInUnits("Hz")
  //+doppler.getInUnits("Hz");
  //double f1=filter_start.getInUnits("Hz");
  //double f2=filter_end.getInUnits("Hz");
  
  //compute_matched_filter(adc,
  //		   taup,
  //		   chirp_rate,
  //		   cfs,
  //		   f1,
  //		   f2,
  //		   start_frequency,
  //		   N,
  //		   i_npts,
  //		   matched_filter);
  
  
  //--------------------------------------------------------
  //     INPUT PARAMETERS:
  //---------------------------------------------------------
  //    double r_chirplen             !chirp length in seconds
  //    double r_bw                   !bandwidth in Hz
  //    double r_stepsize             !DCG step length in seconds
  //    double r_f0                   !chirp center frequency (Hz)
  //    double r_phase0               !initial phase in chirp
  //    double r_sampfreq             !sampling frequency
  
  //-----------------------------------------------------------------	
  //  OUTPUT PARAMETERS:
  //-----------------------------------------------------------------
  // unsigned int i_ts                  !number of samples generated
  // complex<float*> c_chirp(*)         !complex chirp sampled at ADC rate
  //
  // Array is allocated outside of function long enough to contain the
  // i_ts samples plus whatever zero-padding may be required outside of
  // this function
  
  //-------------------------------------------------------------------
  //    LOCAL VARIABLES:
  //-------------------------------------------------------------------
/*
cout << "compoute_matched_filter_dig" << endl;
cout << r_chirplen << " ";
cout << r_bw << " ";
cout << r_stepsize << " ";
cout << r_phase0 << " ";
cout << r_sampfreq << " ";
cout << r_f0 << " ";
cout << endl;
*/

  int i_numsteps,i_stepnum;
  int i_totsamps;
  double r_samptime,r_dfdt,r_delfreq,r_tottime,r_phasec,r_phase;
  double r_t,r_phase0p,r_ss,r_freq,r_freq_low,r_freq_high;
  
  
  
  //     parameters describing number of steps in DCG chirp
  
  i_numsteps = int(r_chirplen/fabs(r_stepsize)+0.5);
  r_ss = r_stepsize/fabs(r_stepsize);
  
  //    ADC sample time and total length of chirp after ADC quantization
  
  r_samptime = 1.0/r_sampfreq;
  i_totsamps = int(r_chirplen*r_sampfreq+0.5);
  i_ts = (unsigned int) i_totsamps;
  
  if(i_numsteps > 1){                         // Using DCG
    r_delfreq = r_ss*r_bw/(i_numsteps - 1);   // frequency step size in DCG
    r_dfdt = 0.0;                             // no chirpo rate
  }
  else{                                        // ideal chirp
    r_delfreq = 0.0;                           // frequency step size is zero 
    r_dfdt = r_ss*r_bw/(i_totsamps*r_samptime);// chirp rate
    r_ss = 0.0;
  }
  
  r_tottime = i_totsamps*r_samptime;    // total chirp time after quantization
  
  // cancel initial phase variations to guarantee desired initial chirp phase

  r_t = 0.0;
  i_stepnum = int((r_t + 0.5*r_tottime)/fabs(r_stepsize));
  r_phasec = -pi*r_delfreq*fabs(r_stepsize)*
    (i_stepnum)*(i_stepnum - i_numsteps+1);  
  r_phase0p = r_phase0 - r_phasec;

  // set the upper and lower frequency regions to avoid spectrum aliasing

  if(r_sampfreq  <  2.0*r_bw) {      // If chirp is too wide for Nyquist
    ErrorMessage e("dig_chirp_hensley: Chirp bandwidth violates Nyquist");
    e.throwMe();
  }
  else{                               
    if(r_f0 > 0){
      r_freq_low = 0.01*(r_sampfreq/2.)  + MAX((r_bw/2 - r_f0),0);
      r_freq_high = 0.99*(r_sampfreq/2.) 
	- MAX((r_f0 + r_bw/2 - r_sampfreq/2),0.); 
    }
    else{
      r_freq_low = -0.99*(r_sampfreq/2.0) 
	+ MAX((fabs(r_f0) + r_bw/2 - r_sampfreq/2.0),0.0);
      r_freq_high = -0.01*(r_sampfreq/2.0) - MAX((r_bw/2 - fabs(r_f0)),0.0);
    }
  }

  //-----------------------------------------------------------
  //    generate the chirp
  //-----------------------------------------------------------
  for(int i=0;i<i_totsamps;i++){
    
    r_t = i*r_samptime - r_tottime/2.0;   // time of sample 
    // filter is 0 centered
    
    
    i_stepnum = int((r_t + 0.5*r_tottime)/fabs(r_stepsize)); // chirp step num
    
    r_phasec = -pi*r_delfreq*fabs(r_stepsize)*(i_stepnum)*
      (i_stepnum - i_numsteps + 1);  
    
    r_phase = 2.0*pi*
      (r_f0 + i_stepnum*r_delfreq - r_ss*r_bw/2 + 0.5*r_dfdt*r_t)*r_t 
      + r_phasec + r_phase0p;
    
    r_freq = r_f0 + i_stepnum*r_delfreq - r_ss*r_bw/2 + r_dfdt*r_t;
    
    // if the freq exceeds the sampling frequency then zero out to 
    // avoid aliasing the spectrum
    
    if(r_freq >= r_freq_low && r_freq <= r_freq_high){
      float r_chirp = cos(r_phase);          // real chirped signal
      matched_filter(i) = complex<float>((float)cos(r_phase), (float) sin(r_phase)); // complex version (imag=0)
    }         
    else{
      matched_filter(i) = complex<float>(0.0,0.0);
    }         
    
  } // end for
 
  //   complex<float>* chirp= new complex<float>[N];
  // for(int i=0;i<N;++i)
  //   chirp[i]=complex<float>(0,0);
  // as modified by Bryan Stiles
  // dig_chirp_hensley(r_chirplen, 
  //	      r_bw,
  //	      r_stepsize, 
  //	      r_f0, 
  //	      r_phase0, 
  //	      r_sampfreq, 
  //	      i_npt, 
  //	      chirp);
  // matched_filter.resize(N);
  // for(int i=0;i<N;++i){
  //   matched_filter(i)=chirp[i];
  //   cout<<"chirp "<< real(chirp[i])<<" "<<imag(chirp[i])<<endl;
  //   cout<<"match "<< real(matched_filter(i))<<" "<<imag(matched_filter(i))<<endl;
  // }
  // delete[] chirp;
  //--------------------------------------------------------
  //     INPUT PARAMETERS:
  //---------------------------------------------------------
  //    double r_chirplen             !chirp length in seconds
  //    double r_bw                   !bandwidth in Hz
  //    double r_stepsize             !DCG step length in seconds
  //    double r_f0                   !chirp center frequency (Hz)
  //    double r_phase0               !initial phase in chirp
  //    double r_sampfreq             !sampling frequency

  //-----------------------------------------------------------------	
  //  OUTPUT PARAMETERS:
  //-----------------------------------------------------------------
  // unsigned int i_ts                  !number of samples generated
  // complex<float*> c_chirp(*)         !complex chirp sampled at ADC rate
  //
  // Array is allocated outside of function long enough to contain the
  // i_ts samples plus whatever zero-padding may be required outside of
  // this function

  //-------------------------------------------------------------------
  //    LOCAL VARIABLES:
  //-------------------------------------------------------------------

  }


double vector_dist(double x[3], double y[3]){
  double dx[3]={x[0]-y[0],x[1]-y[1],x[2]-y[2]};
  return(sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]));
}



// simple/fast dot product routine
double dot(double v1[3], double v2[3]){
  double sum=0;
  for(int c=0;c<3;c++) sum+=v1[c]*v2[c];
  return(sum);
}

//simple fast vector rotation (vector matrix multiply) routine
//multiplication direction agrees with SPICE convention ??
void rotate_vector(double vfrom[3], double rotmat[3][3], double vto[3])
{
  for(int i=0;i<3;i++){
    vto[i]=0;
    for(int j=0;j<3;j++){
      vto[i]+=vfrom[j]*rotmat[i][j];
    }
  }
}


// Fast Routine for computing area of planar triangle on a sphere.
double getFlatTriangleArea(double p1[3], double p2[3], double p3[3]){
  double p12[3],p13[3];
  for(int c=0;c<3;c++){
    p12[c]=p2[c]-p1[c];
    p13[c]=p3[c]-p1[c];
  }
  double x=p12[1]*p13[2]-p12[2]*p13[1];
  double y=p12[2]*p13[0]-p12[0]*p13[2];
  double z=p12[0]*p13[1]-p12[1]*p13[0];
  return(0.5*sqrt(x*x+y*y+z*z));
}

// Fast Routine for computing area of triangle on a sphere.
double getSphericalTriangleArea(double p1[3], double p2[3], double p3[3],
				double r, double& s12, double& s23,
				double& s31)
{
  
  double area=0.0;
  // Side lengths from the earth center angle.
  s12 = r*acos(dot(p1,p2)/(r*r));
  s23 = r*acos(dot(p2,p3)/(r*r));
  s31 = r*acos(dot(p3,p1)/(r*r));

  if (s12 != 0.0 & s23 != 0.0 & s31 != 0.0)
    {       // 3 points are distinct, so compute area
      // Cosine law
      double A1 = acos(-(s23*s23 - s12*s12 - s31*s31)/
		       s12/s31/2.0);
      area = 0.5*s12*s31*sin(A1);
    }
  return(area);
}


// Uvar/PositionVector Routine for computing area of triangle on a sphere.

Uvar getSphericalTriangleArea(const PositionVector& p1, 
			      const PositionVector& p2, 
			      const PositionVector& p3, 
			      const Uvar r,
			      Uvar& s12, Uvar& s23,
			      Uvar& s31)
{ 
  Uvar area=0.0;
  // Side lengths from the earth center angle.

  s12 = r*acos(dot(p1,p2)/(r*r));
  s23 = r*acos(dot(p2,p3)/(r*r));
  s31 = r*acos(dot(p3,p1)/(r*r));

  if (s12 != 0.0 & s23 != 0.0 & s31 != 0.0)
    {       // 3 points are distinct, so compute area
      // Cosine law
      Uvar A1 = acos(-(s23*s23 - s12*s12 - s31*s31)/
		       s12/s31/2.0);
      area = 0.5*s12*s31*sin(A1);
    }
  return(area);
}

Uvar getLatLonPixelAreaOnSphere(const Uvar& lat, const Uvar& lon,
				const Uvar& latres, const Uvar& lonres){

  Frame tbf(default_target_frame_spice_id,default_target_spice_id);
  PositionVector p[4];
  double dummy_time=0.0;
  DirectionVector normal(tbf,dummy_time,0,0,1);
  int p_idx=0;

  Uvar target_radius=(default_target_radii.magnitude())/sqrt(3.0);

  for(int i=0;i<2;i++){
    Uvar lat2=lat-latres/2+i*latres;
    for(int j=0;j<2;j++){
      Uvar lon2=lon-lonres/2+j*lonres;
      normal.setPlanetodetic(lat2,lon2);
      p[p_idx]=normal;
      p[p_idx]*=target_radius;
      p_idx++;
    }
  }


  Uvar dummy1,dummy2,dummy3;
  Uvar area=getSphericalTriangleArea(p[0],p[1],p[2],target_radius,
				     dummy1,dummy2,dummy3);
  
  area+=getSphericalTriangleArea(p[1],p[2],p[3],target_radius,
				 dummy1,dummy2,dummy3);
  return(area);
}

//-------------------------------------------------------------------------
// routine for computing look direction pairs from doppler and range
//-------------------------------------------------------------------------
int
fast_doppler_range_to_TBF_look(double doppler_in_hz, double range_in_km,
			       double target_rad,const double p[3], 
			       const double v[3],
			       double lambda, double u1[3], double u2[3]){
  //------------------------------------------------
  // Input parameters
  // doppler_in_hz = desired doppler in Hz
  // range_in_km = desired range in km
  // target_rad = radius of target body
  // p = s/c position in TBF (km)
  // v = s/c velocity in TBF (km/s)
  // lambda  = center wavelength of transmitted signal (km)
  //----------------------------------------------
  // Output parameters
  // u1,u2 pair of direction vectors in TBF pointing toward
  // desired doppler,range pair on surface of target body

  double spd=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  double pos_mag=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);  
  double pos_in_rad=pos_mag/target_rad;

  // scale pos and velocity vectors to unity
  double p1,p2,p3,v1,v2,v3;
  p1=p[0]/pos_mag;
  p2=p[1]/pos_mag;
  p3=p[2]/pos_mag;
  v1=v[0]/spd;
  v2=v[1]/spd;
  v3=v[2]/spd;

  double range_in_radius=range_in_km/target_rad;
  double A =  (1 - pos_in_rad*pos_in_rad 
		 - range_in_radius*range_in_radius)
    /(2.0 * range_in_radius * pos_in_rad);

  double B =  lambda * doppler_in_hz/(2 * spd);


  double C = (v3 * A - p3*B)/(v3*p1 - p3*v1);
  double D = -(v3*p2 - p3*v2)/(v3*p1 - p3*v1);
  double E = (A - p1 *C)/p3;
  double F = -(p2 + p1*D)/p3;
  
  double sol = (C*D + E*F)*(C*D + E*F) - (D*D + F*F + 1.0)*(C*C + E*E-1.0);

  if (sol > 0){
     // get 1 solution
     u1[1] = ( -(C*D + E*F) -sqrt( sol))/(D*D + F*F + 1.0);
     u1[0] = C + D*u1[1];
     u1[2] = E + F*u1[1];

     // get the other solution
     u2[1] = ( -(C*D + E*F) +sqrt( sol))/(D*D + F*F + 1.0);
     u2[0] = C + D*u2[1];
     u2[2] = E + F*u2[1];
     return(1);
  }
  return(0);
}

//--------------------
// Fast routine for determining a surface intercept on a triaxial body
//
// Input arguments
// scpos = observer 3-D position (km) in target body fixed coords
// look = unit vector in look direction
// radii = 3 radii of target body (km)
// 
// Output argument
// pos = location of surface intercept (km)
// r = reference for range to surface in km
// returns 0 if off limb
//---------------
int
get_surface_intercept_triaxial(double scpos[3], double look[3], double radii[3], double pos[3], double& r){
  double invrad2[3];
  invrad2[0]=1/(radii[0]*radii[0]);
  invrad2[1]=1/(radii[1]*radii[1]);
  invrad2[2]=1/(radii[2]*radii[2]);
  
  double A=look[0]*look[0]*invrad2[0]+look[1]*look[1]*invrad2[1]+look[2]*look[2]*invrad2[2];
  double B=2*(scpos[0]*look[0]*invrad2[0]+scpos[1]*look[1]*invrad2[1]+scpos[2]*look[2]*invrad2[2]);
  double C=scpos[0]*scpos[0]*invrad2[0]+scpos[1]*scpos[1]*invrad2[1]+scpos[2]*scpos[2]*invrad2[2]-1;
  double det=B*B-4*A*C;
  if(det<0) return(0);
  double denom=1/(2*A);
  double dr = (sqrt(det))*denom;
  r = -B*denom+dr;
  double r2 = -B*denom-dr;
  if(r2>0 && r2<r) r=r2;
  if(r<0) return(0);

  pos[0]=scpos[0]+look[0]*r;
  pos[1]=scpos[1]+look[1]*r;
  pos[2]=scpos[2]+look[2]*r;

  return(1);
}
//----------------------------------------------------------------------------
// Function Name: BAQEncode
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This fucntion performs calculating the BAQ threshold, and encoding the 
// burst echo, where
// (1) ieb is the IEB class.
// (2) echo is the burst echo.
// (3) thresh is the calculated BAQ threshold.
// (4) word is the encoded burst echo.
//----------------------------------------------------------------------------
void BAQEncode(Ieb&ieb, const Uvec& echo, SUvec& thresh, CharUvec& word)
{

  unsigned int N_sample = echo.size();
  Baq baq;

  Charvec echo_8("echo_8",N_sample);
  Fvec echo_n("echo",N_sample);
  
  for(unsigned int n_point = 0; n_point < N_sample; n_point++)
    {
      echo_n(n_point) = echo(n_point).getInUnits("");
    }
  
  echo_8 = baq.ADC(echo_n);
  
  baq.setIeb(ieb);
  Ivec Thresh("Thresh");
  //Array1D<unsigned short> thresh("thresh",12);
  baq.compuThreshold(echo_n, Thresh);
  
  Charvec words("words");
 
  baq.Encode_Nbit(echo_n, Thresh, words);
  baq.Encode_thresh(Thresh,thresh);
  
  word = words;
}





//----------------------------------------------------------------------------
// Function Name: ShiftEcho
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This fucntion shifts the echo array by time dt, where
// (1) echo_data is the burst echo.
// (2) dt is the time for echo shifting.
// (3) SR is the sampling rate (ADC).
//----------------------------------------------------------------------------
void ShiftEcho(Fvec& echo_data, const Uvar& dt, const Uvar& SR)
{
  unsigned int N = echo_data.size();
  int  dn = (int) round_double((dt*SR).getInUnits(""));
  Array1D<float> shifted_echo_data("se",echo_data.size());
  

  shifted_echo_data = 0.0;//default
  
  for(unsigned int i=0;i<N;++i) 
    {
      if(int(i)<dn) continue;//go to next i
      unsigned int ip =(unsigned int) (i - dn);
      if(ip < N)
	shifted_echo_data(ip)=echo_data(i);
    }
  echo_data = shifted_echo_data;
}

//----------------------------------------------------------------------------
// Function Name: Convert2Pow
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function converts an integer n into the integer of the power of 2,
// and then returns the value back to n.
//----------------------------------------------------------------------------
void Convert2Pow(unsigned int& n)
{
  unsigned int m = (unsigned int)(log((double)n)/log(2.));
  if(n == (unsigned int)pow(2,m))
    n = (unsigned int)pow(2., m);
  else
    n = (unsigned int)pow(2., m+1);
  
}

//----------------------------------------------------------------------------
// Function Name: Convert2DEcho
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function converts a 1D echo array into a 2D matrix, where
// (1) echo_data is the 1D echo array.
// (2) echo2D is the 2D echo matrix.
// (3) N_pulse is the number of pulses.
// (4) N_in_pri is the number of sampling points in PRI.
//----------------------------------------------------------------------------
void Convert2DEcho(const Fvec& echo_data, Fmat& echo2D, 
		   const unsigned int& N_pulse, const unsigned int& N_in_pri)
{
  unsigned int Ndop_bins = N_pulse;
  Convert2Pow(Ndop_bins);
  unsigned int Nrange_bins = N_in_pri;
  Convert2Pow(Nrange_bins);
  //cout<<N_pulse<<"  "<<Ndop_bins<<endl;
  //cout<<N_in_pri<<"  "<<Nrange_bins<<endl;
  echo2D.resize(Ndop_bins, Nrange_bins);
  echo2D = 0.0;

  for(unsigned int i = 0; i < N_pulse; i++)
    for(unsigned int j = 0;j < N_in_pri;j++)
      {
	unsigned int n = i*N_in_pri + j;
	echo2D(i,j) = echo_data(n);
      }
}

void Convert2DEcho(const Fvec& echo_data, 
		   const unsigned int& start_index,
		   const unsigned int& end_index, 
		   const unsigned int& N_samples_per_window,
		   const unsigned int& N_pulses_inside_echo,
		   const unsigned int& N_samples_per_pri,
		   Fmat& echo2D)
  {
    if(end_index > echo_data.size()) ErrorMessage("end index is larger than data size").throwMe();
    
    unsigned int Ndop_bins=N_pulses_inside_echo;
    Convert2Pow(Ndop_bins);
    unsigned  int Nrange_bins=N_samples_per_window;
    Convert2Pow(Nrange_bins);
    echo2D.resize(Ndop_bins,Nrange_bins);
    echo2D=0.0;
    
    for(unsigned int i=0; i < N_pulses_inside_echo; i++) 
      for(unsigned j=0; j < N_samples_per_window; j++){
	unsigned int n = i*N_samples_per_pri + j+start_index;
	if(n<end_index)
	  echo2D(i,j) = echo_data(n); 
      }
    //Plot1DReal(echo_data);
    //Plot2DReal(echo2D);
  }

//----------------------------------------------------------------------------
// Function Name: Convert2DReal2Complx
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function converts a real 2D echo matrix into a complex 2D matrix,
// where
// (1) echo2D is the 2D real echo matrix.
// (2) c_echo2D is the 2D complex echo matrix. 
// (3) band_selection: -1: negative band, +1 positive band
//----------------------------------------------------------------------------
void Convert2DReal2Complx(const Fmat& echo2D, CFmat& c_echo2D,
const int& band_selection)
{
  if(!(band_selection==1 || band_selection==-1))
    ErrorMessage("band selection should be 1 or -1 ").throwMe();

  unsigned int M, N;
  echo2D.size(M, N);
  c_echo2D.resize(M, N/2);

  for(unsigned int i=0; i < M; i++)
    {
      //for each pulse
      CFvec echo_cmplx("echo_cmplx",N);
      for(unsigned int j = 0; j < N; j++)
	echo_cmplx(j) = complex<float>(echo2D(i,j),0.);

      CFvec pulse_fft("pulse_fft",N);
      fft(echo_cmplx, pulse_fft);

      // select negative spectrum
      CFvec pulse_fft_half("pulse_fft_half",N/2);
      CFvec cmplx_pulse("cmplx_pulse",N/2);
      if(band_selection==-1){
	for(unsigned int n = N/2; n < N; n++)
	  pulse_fft_half(n-N/2) = pulse_fft(n);
      }
      else{
	for(unsigned int n = 0; n < N/2; n++)
	  pulse_fft_half(n) = pulse_fft(n);
      }
      // inverse fft 
      ifft(pulse_fft_half, cmplx_pulse);

      //save it
      for(unsigned int n = 0; n < N/2; n++)
	c_echo2D(i,n) = cmplx_pulse(n);
    }    
}


//----------------------------------------------------------------------------
// Function Name: Convert2DReal2Complx
// Author: YG
// Date:   4/21/04
// Function Description:
//      This function converts a real 2D echo matrix into a complex 2D matrix,
// convert 1D echo to 1D complex echo 
// (3) band_selection: -1: negative band, +1 positive band
//----------------------------------------------------------------------------
void ConvertRealToComplex(const Fvec& echo, CFvec& c_echo,
			  const int& band_selection)
  {
    if(!(band_selection==1 || band_selection==-1))
      ErrorMessage("band selection should be 1 or -1 ").throwMe();
    unsigned int N=echo.size();
    Convert2Pow(N);
    Fvec echo_tmp("",N);
    echo_tmp=0.0;
    for(unsigned int i=0;i<echo.size();++i) echo_tmp(i)=echo(i);
    c_echo.resize(N/2);
    c_echo=0.0;

    CFvec echo_cmplx("echo_cmplx",N);
    for(unsigned int j = 0; j < N; j++)
      echo_cmplx(j) = complex<float>(echo_tmp(j),0.);
    
    CFvec pulse_fft("pulse_fft",N);
    fft(echo_cmplx, pulse_fft);

    
    CFvec pulse_fft_half("pulse_fft_half",N/2);
    CFvec cmplx_pulse("cmplx_pulse",N/2);
    if(band_selection==-1){
      for(unsigned int n = N/2; n < N; n++)
	pulse_fft_half(n-N/2) = pulse_fft(n);
    }
    else{
      for(unsigned int n = 0; n < N/2; n++)
	pulse_fft_half(n) = pulse_fft(n);
    }
    // inverse fft 
    ifft(pulse_fft_half, cmplx_pulse);

    //save it
    for(unsigned int n = 0; n < N/2; n++)
      c_echo(n) = cmplx_pulse(n);


  }


//----------------------------------------------------------------------------
// Function Name: ComputPulseCorrelaton
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function performs the pulse-to-pulse correlation for a burst echo,
// where
// (1) c_echo2D is the 2D complex echo matrix.
// (2) c_corr2D is the 2D complex correlation matrix.
//----------------------------------------------------------------------------
//void ComputPulseCorrelation(const CFmat& c_echo2D, CFmat& c_corr2D)
//{
//unsigned int M, N;
//c_echo2D.size(M,N);
//c_corr2D.resize(M,N);
//c_corr2D = 0.0;
//
//for(unsigned int i=1; i < M; i++)
//  for(unsigned int j = 0; j < N; j++)
//    c_corr2D(i-1,j) = c_echo2D(i,j)*conj(c_echo2D(i-1,j));
//}

//----------------------------------------------------------------------------
// Function Name: AveragePulseCorrelation
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function averages the 2D correlation matrix over pulses and entire 
// values to get a 1D array and constant correlations, where
// (1) c_corr2D is the 2D complex correlaton matrix.
// (2) c_corr1D is the averaged 1D complex correlation array.
// (3) c_corr0D is the averaged complex correlation constant.
//---------------------------------------------------------------------------- 
//void AveragePulseCorrelation(const CFmat& echo2D_c, CFvec& c_corr1D, complex<float>& c_corr0D)
//{
//unsigned int M, N;
//c_corr2D.size(M, N);
//c_corr1D.resize(N);
//c_corr1D = 0.0;
//c_corr0D = 0.0;
//
//for(unsigned int i = 0; i < M; i++)
//  for(unsigned int j = 0; j < N; j++)
//    {
//c_corr1D(j) = c_corr1D(j) + c_corr2D(i,j)/(float)M;
//c_corr0D = c_corr0D + c_corr2D(i,j)/(float)N/(float)M;
//    }
//}

//----------------------------------------------------------------------------
// Function Name: ComputFractDopCentroid
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This funtion computs the fractional Doppler centroid frequency by using
// the complex correlations, where
// (1) c_corr2D is the 2D complex correlation matrix.
// (2) dop_frac1D is the fractional Doppler centroid frequency as a function of
//     fast time in PRI.
// (3) dop_frac0D is the averaged Doppler centroid frequency for a burst. It is 
//     a constant.
//----------------------------------------------------------------------------- 
void ComputFractDopCentroid(const CFmat& echo2D_c, Dvec& dop_frac1D, double& dop_frac0D)
{
  unsigned int M, N;
  echo2D_c.size(M,N);

  CFmat c_corr2D("",M,N); 
  c_corr2D=0.0;
  
  CFvec c_corr1D("c_corr1D",N);
  c_corr1D=0.0;

  if(dop_frac1D.size()!= N) dop_frac1D.resize(N);
  dop_frac1D=0.0;

  complex<float> c_corr0D;
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

}

//----------------------------------------------------------------------------
// Function Name: ComputIntDopCentroid
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This functin computes the integer Doppler centroid frequency by using
// a predicted boresight Doppler frequency, where
// (1) bore_dop is the boresight Doppler frequency.
// (2) pri is the PRI.
// (3) dop_int is the integer part of the Doppler frequency.
//----------------------------------------------------------------------------
void ComputIntDopCentroid(const Uvar& bore_dop, const Uvar& pri, int& dop_int)
{
  if(bore_dop > 0)
    dop_int = (int)((bore_dop*pri).getInUnits("") + 0.5);
  else if(bore_dop < 0)
    dop_int = (int)((bore_dop*pri).getInUnits("") - 0.5);
  else 
    dop_int = 0; 
}

//----------------------------------------------------------------------------
// Function Name: CompareDopCentroid
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function plots and compares Doppler centroid frequencies obtained 
// from data and calculation, where
// (1) exact_dop is the Doppler frequency calculated by target geometry.
// (2) dop_frac1D is the 1D fractional Doppler centroid calculated from data.
// (3) dop_int is the integer Doppler centroid.
// (4) pri is the PRI.
//----------------------------------------------------------------------------
void CompareDopCentroid(const Uvar& exact_dop, const Fvec& dop_frac1D, const int& dop_int, const Uvar& pri)
{
  unsigned int N = dop_frac1D.size();
  Uvec data_dop("data_dop",N);
  Uvec anal_dop("data_dop",N);

  for(unsigned int j = 0; j < N; j++)
    {
      data_dop(j) = dop_frac1D(j);
      anal_dop(j) = exact_dop*pri - dop_int;
    }

  Plot a;
  a.addY(data_dop,"re:PRF",line("solid","black",1),sym("none"));
  a.addY(anal_dop,"re:PRF",line("solid","red",1),sym("none"));
  a.setTitle("DC by Data (black) vs. DC by Calulation (red)");
  a.show("x");

}

//-----------------------------------------------------------------
// Function Name: PRFHopping
// Author: YZ
// Date:   10/30/03
// Function Description:
//      This function calculates absolute Doppler centroid 
// by using PRF hopping algorithm, where
//
// Inputs:
//    PRFs: Array which stores three PRFs
//    est_dop: Array which stores three estimated Doppler centroid
//    frac_dop: Array which stores three calculated fractional 
//           Doppler centroid from data
// Output:
//  hop_dop: Array which stores three absolute Doppler centroid
//-----------------------------------------------------------------
void PRFHopping(const Uvec& PRFs, const Uvec& est_dop, const Uvec& frac_dop, Uvec& hop_dop)
{
  //---------------------------
  // doppler frequency walking
  //---------------------------
  Uvar dAC = est_dop(0) - est_dop(2);
  Uvar dBC = est_dop(1) - est_dop(2);
  //-----------------------------
  // redefine fractional doppler
  // ddop(0) = Dop_B - Dop_C
  // ddop(1) = Dop_A - Dop_B
  //-----------------------------
  Uvar hat_fdA = frac_dop(0) - dAC;
  Uvar hat_fdB = frac_dop(1) - dBC;
  Uvar hat_fdC = frac_dop(2);

  //----------------------
  // doppler search range
  //----------------------
  Uvar fdmin = est_dop(2) - 3*PRFs(2);
  Uvar fdmax = est_dop(2) + 3*PRFs(2);

  unsigned int N = 200;
  Uvec e_fdC("e_fdC",N);
  for(unsigned int i = 0;i < N; i++)
    e_fdC(i) = fdmin + i*(fdmax - fdmin)/N;

  Uvec dA("dA",N);
  Uvec dB("dB",N);
  Uvec delta("delta",N);
  //-------------------
  // error calculation
  //-------------------
  for(unsigned int i = 0; i < N; i++)
    {
      dA(i)=(e_fdC(i)-hat_fdA)-((int)((e_fdC(i)-hat_fdA)/PRFs(0)).getInUnits(""))*PRFs(0);
      dB(i)=(e_fdC(i)-hat_fdC)-((int)((e_fdC(i)-hat_fdC)/PRFs(2)).getInUnits(""))*PRFs(2);
      delta(i) = sqrt(pow(dA(i),2.0) + pow(dB(i),2.0));
    }

  //-------------------
  // find minimum of d
  //-------------------
  
  //Uvar d_min = 1.e50;
  unsigned int n_indx;
  delta.min(n_indx);
  //  for(unsigned int i = 0; i < N; i++)
  //if(delta(i) <= d_min)
  //  {
  //d_min = delta(i);
  //n_indx = i;
  //  }
  
  //--------------------------------
  // doppler frequency at minimum d
  //--------------------------------
  hop_dop(2) = e_fdC(n_indx);
  hop_dop(0) = hop_dop(2) + dAC;
  hop_dop(1) = hop_dop(2) + dBC;

  /*Plot a;
  a.addY(dA,"",line("solid","black",1),sym("none"));
  a.addY(dB,"",line("solid","blue",1),sym("none"));
  a.addY(delta,"",line("solid","red",1),sym("none"));
  a.setTitle("Errors");
  a.show("x");*/

}


//----------------------------------------------------------------------------
// Function Name: Plot1DReal
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function plots a real 1D array.
//----------------------------------------------------------------------------
void Plot1DReal(const Fvec& data1D)
{
  unsigned int N = data1D.size();
  Uvec p_data1D("p_data1D",N);
  p_data1D = data1D;

  Plot a;
  a.addY(p_data1D,"",line("solid","black",1),sym("none"));
  a.setTitle("1D Real Data");
  a.show("x");
}

//----------------------------------------------------------------------------
// Function Name: Plot2DReal
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function plots a real 2D matrix.
//----------------------------------------------------------------------------
void Plot2DReal(const Fmat& data2D)
{
  unsigned int M, N;
  data2D.size(M, N);
  Uvec p_data2D("p_data2D",N);
  for(unsigned int i = 0; i < M; i++)
    {
      for(unsigned int j = 0; j < N; j++)
	p_data2D(j) = data2D(i,j);

      Plot a;
      a.addY(p_data2D,"",line("solid","black",1),sym("none"));
      a.setTitle("2D Real Data");
      a.show("x");
    }
}

//----------------------------------------------------------------------------
// Function Name: Plot1DComplx
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function plots a complex 1D array.
//----------------------------------------------------------------------------
void Plot1DComplx(const CFvec& c_data1D)
{
  unsigned int N = c_data1D.size();
  Uvec p_data1D_real("p_data2D_real",N);
  Uvec p_data1D_imag("p_data2D_imag",N);

  for(unsigned int j = 0; j < N; j++)
    {
      p_data1D_real(j) = real(c_data1D(j));
      p_data1D_imag(j) = imag(c_data1D(j));
    }
  
  Plot a;
  a.addY(p_data1D_real,"",line("solid","black",1),sym("none"));
  a.addY(p_data1D_imag,"",line("solid","red",1),sym("none"));
  a.setTitle("1D Complex Data");
  a.show("x");
    
}

//----------------------------------------------------------------------------
// Function Name: Plot2DComplx
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function plots a complex 2D matrix.
//----------------------------------------------------------------------------
void Plot2DComplx(const CFmat& c_data2D)
{
  unsigned int M, N;
  c_data2D.size(M, N);
  Uvec p_data2D_real("p_data2D_real",N);
  Uvec p_data2D_imag("p_data2D_imag",N);
  for(unsigned int i = 0; i < M; i++)
    {
      for(unsigned int j = 0; j < N; j++)
	{
	  p_data2D_real(j) = real(c_data2D(i,j));
	  p_data2D_imag(j) = imag(c_data2D(i,j));
	}

      Plot a;
      a.addY(p_data2D_real,"",line("solid","black",1),sym("none"));
      a.addY(p_data2D_imag,"",line("solid","red",1),sym("none"));
      a.setTitle("1D Complex Data");
      a.show("x");
    }

}
//----------------------------------------------------------------------------

unsigned int get_next_power_of_2(unsigned int n){
  unsigned int m=2;
  while(m<n){
    m<<=1;
  }
  return(m);
}

void simpleRangeDopplerCompression(const CFmat& echo2D_c,
			    const CFvec& fft_matched_filter,
			    const Uvar& prf,
			    const Uvar& bore_doppler,
			    const Uvar& adc,
			    const Uvar& range_start,
			    Uvec& range_axis,	
			    Uvec& doppler_axis,
			    CFmat& echo2D_rc,
			    CFmat& echo2D_rdc)
				  
  {
    unsigned int M, N;
    unsigned int Noversample=4;
    echo2D_c.size(M,N);
    echo2D_rc.resize(M,Noversample*N);
    
    //--------------------
    //range axis
    //---------------------			     
    range_axis.resize(Noversample*N);
    for(unsigned int i=0; i<Noversample*N;++i){
      range_axis(i)=2.0*double(i)/adc*speed_light/2/float(Noversample);
      range_axis(i)+=range_start;
    }

    //----------------------------------------------
    //declare temporary variables for recycling during range compression
    //-------------------------------------------
    CFvec x_data(" ",N);
    CFvec fft_x_data("",N);
    CFvec ifft_x_data("",Noversample*N);
    
    //first range compress
    for(unsigned int ii=0; ii < M; ii++){
      //for each pulse
      x_data= echo2D_c.getRow(ii);//complex data already
     
      //fft 
      fft(x_data,fft_x_data);
      
      // multiply conjugate of matched filter
      for(unsigned int jj=0;jj<N;++jj)
	fft_x_data(jj)*= conj(fft_matched_filter(jj+N));
      
      CFvec tmp_x("",Noversample*N);
      tmp_x=complex<float>(0,0);
      for(unsigned int jj=0; jj<N;++jj)
      tmp_x(jj)=fft_x_data(jj);

      ifft(tmp_x,ifft_x_data);
      
      for(unsigned int jj=0;jj<Noversample*N;++jj)
	echo2D_rc(ii,jj)=ifft_x_data(jj);
    }
   



    //------------------------
    //apply window to ffted matched filter
    //----------------------
    Array1D<double> hann_window("",M);
    for(unsigned int ii=0; ii<M;++ii){     
      hann_window(ii)= 0.5*(1.0 - cos(2.0*M_PI*float(ii)/float(M-1)));     
    }
  
  

    echo2D_rdc.resize(M*Noversample,Noversample*N);
    //-------------------------
    // azimuth compress
    //--------------------------
    CFvec azimuth("doppler process ",M);
    CFvec azimuth_fft("ffted azimuth",Noversample*M);
    Uvar pri = 1/prf;
    
    doppler_axis.resize(M*Noversample);
    for(unsigned int i=0;i<doppler_axis.size();++i) 
      doppler_axis(i)=bore_doppler- prf/2.0 +prf*double(i)/double(M-1)/float(Noversample);

    complex<float> shift;
    float pre_factor;
    Uvar t0= pri/2.0;
    Uvar fdop_rate=Uvar(0,"Hz/s");
    Uvar freq_rate;
    Uvar t_midle;
    Uvar freq_offset;
    Uvar t_middle;
    
    for(unsigned int jj=0;jj<Noversample*N;++jj){ 
      //for each range line
      azimuth=complex<float>(0,0);//reset
      for (unsigned int ii=0;ii<M;++ii){
	//for each azimuth line
	freq_offset= 2.0*pi*(bore_doppler - prf/2.0);
	freq_rate = pi*fdop_rate;
	t_middle = t0 + pri*double(ii);
	pre_factor= float( ( (freq_rate*t_middle + freq_offset)*t_middle).getInUnits(""));
	shift=exp(complex<float>(0,-1)*(float)(pre_factor));
	azimuth(ii) = echo2D_rc(ii,jj)*shift;
      }
      //--------------
      //fft azimuth  for each range line
      //-------------
      CFvec tmp_azi("",Noversample*M);
      tmp_azi=complex<float>(0,0);//zero pad for over sampling
      for(unsigned int ii=0;ii<M;++ii) tmp_azi(ii)=azimuth(ii);//*float(hann_window(ii));

      fft(tmp_azi,azimuth_fft);//fft
      for(unsigned int ii=0;ii<azimuth_fft.size();++ii)
	echo2D_rdc(ii,jj)=azimuth_fft(ii);//save it 
    }   
  }


void simpleRangeDopplerCompression2(const CFmat& echo2D_c,
				    const CFvec& fft_matched_filter,
				    const Uvar& prf,
				    const Uvar& bore_doppler,
				    const Uvar& fdop_rate,
				    const Uvar& adc,
				    const Uvar& range_start,
				    Uvec& range_axis,	
				    Uvec& doppler_axis,
				    CFmat& echo2D_rc,
				    CFmat& echo2D_rdc)
				  
  {
    unsigned int M, N;
    unsigned int Noversample=4;
    echo2D_c.size(M,N);
    echo2D_rc.resize(M,Noversample*N);
    
    //--------------------
    //range axis
    //---------------------			     
    range_axis.resize(Noversample*N);
    for(unsigned int i=0; i<Noversample*N;++i){
      range_axis(i)=2.0*double(i)/adc*speed_light/2/float(Noversample);
      range_axis(i)+=range_start;
    }

    //----------------------------------------------
    //declare temporary variables for recycling during range compression
    //-------------------------------------------
    CFvec x_data(" ",N);
    CFvec fft_x_data("",N);
    CFvec ifft_x_data("",Noversample*N);
    
    //first range compress
    for(unsigned int ii=0; ii < M; ii++){
      //for each pulse
      x_data= echo2D_c.getRow(ii);//complex data already
     
      //fft 
      fft(x_data,fft_x_data);
      
      // multiply conjugate of matched filter
      for(unsigned int jj=0;jj<N;++jj)
	fft_x_data(jj)*= conj(fft_matched_filter(jj+N));
      
      CFvec tmp_x("",Noversample*N);
      tmp_x=complex<float>(0,0);
      for(unsigned int jj=0; jj<N;++jj)
      tmp_x(jj)=fft_x_data(jj);

      ifft(tmp_x,ifft_x_data);
      
      for(unsigned int jj=0;jj<Noversample*N;++jj)
	echo2D_rc(ii,jj)=ifft_x_data(jj);
    }
   



    //------------------------
    //apply window to ffted matched filter
    //----------------------
    Array1D<double> hann_window("",M);
    for(unsigned int ii=0; ii<M;++ii){     
      hann_window(ii)= 0.5*(1.0 - cos(2.0*M_PI*float(ii)/float(M-1)));     
    }
  
  

    echo2D_rdc.resize(M*Noversample,Noversample*N);
    //-------------------------
    // azimuth compress
    //--------------------------
    CFvec azimuth("doppler process ",M);
    CFvec azimuth_fft("ffted azimuth",Noversample*M);
    Uvar pri = 1/prf;
    
    doppler_axis.resize(M*Noversample);
    for(unsigned int i=0;i<doppler_axis.size();++i) 
      doppler_axis(i)=bore_doppler- prf/2.0 +prf*double(i)/double(M-1)/float(Noversample);

    complex<float> shift;
    float pre_factor;
    Uvar t0= pri/2.0;
    Uvar freq_rate;
    Uvar t_midle;
    Uvar freq_offset;
    Uvar t_middle;
    
    for(unsigned int jj=0;jj<Noversample*N;++jj){ 
      //for each range line
      azimuth=complex<float>(0,0);//reset
      for (unsigned int ii=0;ii<M;++ii){
	//for each azimuth line
	freq_offset= 2.0*pi*(bore_doppler - prf/2.0);
	freq_rate = pi*fdop_rate;
	t_middle = t0 + pri*double(ii);
	pre_factor= float( ( (freq_rate*t_middle + freq_offset)*t_middle).getInUnits(""));
	shift=exp(complex<float>(0,-1)*(float)(pre_factor));
	azimuth(ii) = echo2D_rc(ii,jj)*shift;
      }
      //--------------
      //fft azimuth  for each range line
      //-------------
      CFvec tmp_azi("",Noversample*M);
      tmp_azi=complex<float>(0,0);//zero pad for over sampling
      for(unsigned int ii=0;ii<M;++ii) tmp_azi(ii)=azimuth(ii);//*float(hann_window(ii));

      fft(tmp_azi,azimuth_fft);//fft
      for(unsigned int ii=0;ii<azimuth_fft.size();++ii)
	echo2D_rdc(ii,jj)=azimuth_fft(ii);//save it 
    }   
  }


void simpleRangeDopplerCompression3(const CFmat& echo2D_c,
				    const CFvec& fft_matched_filter,
				    const Uvar& prf,
				    const Uvar& bore_doppler,
				    const Uvar& adc,
				    const Uvar& range_start,
				    const double& oversample_ratio,
				    Uvec& range_axis,	
				    Uvec& doppler_axis,
				    CFmat& echo2D_rc,
				    CFmat& echo2D_rdc)
  
  {
    unsigned int M, N;
    unsigned int Noversample=4;
    echo2D_c.size(M,N);
    echo2D_rc.resize(M,Noversample*N);
    
    //--------------------
    //range axis
    //---------------------			     
    range_axis.resize(Noversample*N);
    for(unsigned int i=0; i<Noversample*N;++i){
      range_axis(i)=2.0*double(i)/adc*speed_light/2/float(Noversample)*oversample_ratio;
      range_axis(i)+=range_start;
    }

    //----------------------------------------------
    //declare temporary variables for recycling during range compression
    //-------------------------------------------
    CFvec x_data(" ",N);
    CFvec fft_x_data("",N);
    CFvec ifft_x_data("",Noversample*N);
    
    //first range compress
    for(unsigned int ii=0; ii < M; ii++){
      //for each pulse
      x_data= echo2D_c.getRow(ii);//complex data already
     
      //fft 
      fft(x_data,fft_x_data);
      
      // multiply conjugate of matched filter
      for(unsigned int jj=0;jj<N;++jj)
	fft_x_data(jj)*= conj(fft_matched_filter(jj+N));
      
      CFvec tmp_x("",Noversample*N);
      tmp_x=complex<float>(0,0);
      for(unsigned int jj=0; jj<N;++jj)
      tmp_x(jj)=fft_x_data(jj);

      ifft(tmp_x,ifft_x_data);
      
      for(unsigned int jj=0;jj<Noversample*N;++jj)
	echo2D_rc(ii,jj)=ifft_x_data(jj);
    }
   



    //------------------------
    //apply window to ffted matched filter
    //----------------------
    Array1D<double> hann_window("",M);
    for(unsigned int ii=0; ii<M;++ii){     
      hann_window(ii)= 0.5*(1.0 - cos(2.0*M_PI*float(ii)/float(M-1)));     
    }
  
  

    echo2D_rdc.resize(M*Noversample,Noversample*N);
    //-------------------------
    // azimuth compress
    //--------------------------
    CFvec azimuth("doppler process ",M);
    CFvec azimuth_fft("ffted azimuth",Noversample*M);
    Uvar pri = 1/prf;
    
    doppler_axis.resize(M*Noversample);
    for(unsigned int i=0;i<doppler_axis.size();++i) 
      doppler_axis(i)=bore_doppler- prf/2.0 +prf*double(i)/double(M-1)/float(Noversample);

    complex<float> shift;
    float pre_factor;
    Uvar t0= pri/2.0;
    Uvar fdop_rate=Uvar(0,"Hz/s");
    Uvar freq_rate;
    Uvar t_midle;
    Uvar freq_offset;
    Uvar t_middle;
    
    for(unsigned int jj=0;jj<Noversample*N;++jj){ 
      //for each range line
      azimuth=complex<float>(0,0);//reset
      for (unsigned int ii=0;ii<M;++ii){
	//for each azimuth line
	freq_offset= 2.0*pi*(bore_doppler - prf/2.0);
	freq_rate = pi*fdop_rate;
	t_middle = t0 + pri*double(ii);
	pre_factor= float( ( (freq_rate*t_middle + freq_offset)*t_middle).getInUnits(""));
	shift=exp(complex<float>(0,-1)*(float)(pre_factor));
	azimuth(ii) = echo2D_rc(ii,jj)*shift;
      }
      //--------------
      //fft azimuth  for each range line
      //-------------
      CFvec tmp_azi("",Noversample*M);
      tmp_azi=complex<float>(0,0);//zero pad for over sampling
      for(unsigned int ii=0;ii<M;++ii) tmp_azi(ii)=azimuth(ii);//*float(hann_window(ii));

      fft(tmp_azi,azimuth_fft);//fft
      for(unsigned int ii=0;ii<azimuth_fft.size();++ii)
	echo2D_rdc(ii,jj)=azimuth_fft(ii);//save it 
    }   
  }



double getFastDfDt(double look[3],
		   double sc_acc[3], 
		   double lambda_in_km, 
		   double speed_in_km_s,
		   double range_in_km, 
		   double fdop){
  //dfdt=2/lamb  * [ dot(a,l)-|v|^2/range] + f^2*lambda/2range
  double adotl=sc_acc[0]*look[0]+sc_acc[1]*look[1]+sc_acc[2]*look[2];
  double dfdt=(2/lambda_in_km)*
    (adotl-speed_in_km_s*speed_in_km_s/range_in_km) +
    (lambda_in_km/2/range_in_km)*fdop*fdop;
  return(dfdt);
}

//----------------------------------------------------------------------------
// Function Name: echo_phase
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function returns the phase for a burst echo at any given time in 
// the receiving window. It is called by getBurstEcho, where
// (1) tr is the sampling time in the receiving window. Its origin is the
//     transmission time of the first pulse.
// (2) range is the distance from the target to the receiver at time tr.
// (3) cfs is the chirp frequency step from ieb.
// (4) csd is the chirp step duration from ieb.
// (5) csf is the chirp start frequency from ieb.
// (6) pri is the PRI from ieb.
//---------------------------------------------------------------------------- 
Uvar echo_phase(const Uvar& tr, const Uvar& range, const Uvar& cfs, const Uvar& csd, 
                const Uvar& csf, const Uvar& pri)
{
  Uvar te_B = 2*range/speed_light;
  double pi = 4.*atan(1.);
  Uvar fo = slo_frequency;                     // local osillating frequency
  Uvar omega_x = 2.*pi*(1378*fo);              // mixing frequency
  Uvar omega_c = 2.*pi*(1377*fo);              // carrier frequency
  Uvar D_omega = 2.*pi*cfs;                    // chirp frequency step
  Uvar omega_csf = 2.*pi*csf - D_omega/2;      // chirp starting frequency modified by CSD
  Uvar gamma = D_omega/2./csd;                 // chirp rate
  int n = int(((tr - te_B)/pri).getInUnits(""));

  Uvar phi = (omega_c - omega_x)*tr            // carrier/mixing frequency term
    - omega_c * te_B                           // time delay term (doppler)
    + omega_csf * (tr - te_B - pri * n)        // linear pulse term
    + gamma * pow(tr - te_B - pri * n, 2.);    // chirp term

  return(phi);
}

//----------------------------------------------------------------------------
// Function Name: solid_pulse
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function returns 1 if the sampling time is within tau_p, and 0 if
// the sampling time is in the gap, where
// (1) tr is the sampling time in the receiving window. Its origin is the
//     transmission time.
// (2) range is the distance from the target to the receiver at time tr.
// (3) pri is the PRI from ieb.
// (4) tau_p is the pulse width from ieb.
//----------------------------------------------------------------------------
double solid_pulse(const Uvar& tr, const Uvar& range, const Uvar& pri, const Uvar& tau_p)
{
  Uvar te_B = 2*range/speed_light;
  int n = int(((tr - te_B)/pri).getInUnits(""));

  double sp = 0;
  if((tr >= te_B + pri*n)&&(tr<te_B + pri*n + tau_p)) sp = 1;

  return(sp);
} 

//----------------------------------------------------------------------------
// Function Name: doppler4target
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function calculates the Doppler frequency at signal time tp for
// a fixed target specified by lat and lon, where
// (1) tp is the time of signal for which the Doppler frequency is calculated.
// (2) target_name is the target name, e.g., titan.
// (3) lat and lon are the latitude and longitude of a point target.
//----------------------------------------------------------------------------
Uvar doppler4target(const Uvar& tp, const string& target_name, const Uvar& lat, 
                    const Uvar& lon)
{
  Frame ftitan("IAU_TITAN",target_name);
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",tp,"NONE");
		  
  TargetGeom tg(tp);
  tg.setState(sc_state);
  tg.setTarget(target_name,ftitan);
  tg.setLatLon(lat,lon);
  DirectionVector look = tg.lookDirection();
  
  Uvar dop_freq = tg.doppler(speed_light/carrier_frequency);
  
  return(dop_freq);
}

//----------------------------------------------------------------------------
// Function Name: SCState
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function specifies the Cassini spacecraft state for a given time
// tp, where
// (1) tp is the time of interest.
// (2) target_name is the target name, e.g., titan.
// (3) Rsc is the array which stores the spacecraft position vector.
// (4) Vsc is the array which stores the spacecraft velocity vector.
//---------------------------------------------------------------------------- 
void SCState(const Uvar& tp, const string& target_name, Uvec& Rsc, Uvec& Vsc)
{
  Frame ftitan("IAU_TITAN",target_name);
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",tp,"NONE");

  Rsc(0)=sc_state.position()[PositionVector::X];
  Rsc(1)=sc_state.position()[PositionVector::Y];
  Rsc(2)=sc_state.position()[PositionVector::Z];
  
  Vsc(0)=sc_state.velocity()[FloatVector::X];
  Vsc(1)=sc_state.velocity()[FloatVector::Y];
  Vsc(2)=sc_state.velocity()[FloatVector::Z];
}

//----------------------------------------------------------------------------
// Function Name: fitRange
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function returns the range as a function of the time tp within 
// pulses by linear curve prediction, where
// (1) Rsc and Vsc are spacecraft position and velocity vectors, respectively.
// (2) Rt is the target possition vector.
// (3) tp is the time of interest.
//----------------------------------------------------------------------------
Uvar fitRange(const Uvec& Rsc, const Uvec& Vsc, const Uvec& Rt, const Uvar& tp)
{
  
  Uvar Rx = (Rsc(0) - Rt(0)) + Vsc(0)*tp;
  Uvar Ry = (Rsc(1) - Rt(1)) + Vsc(1)*tp;
  Uvar Rz = (Rsc(2) - Rt(2)) + Vsc(2)*tp;
  
  Uvar Range = sqrt(pow(Rx,2.)+pow(Ry,2.)+pow(Rz,2.));
  return(Range);
}

//----------------------------------------------------------------------------
// Function Name: SetTargets
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function sets point target at the intercept of the looking direction
// and the target surface for 5 beam, where
// (1) t is the time defines the spacecraft looking direction and state.
// (2) target_name is the target name, e.g., titan.
// (3) latb and lonb are arrays which stores the latitude and longitude, 
//     respectively, for targets on 5 beam footprints.
// (4) Xtb, Ytb, and Ztb are arrays which stores the target position for targets
//     on 5 beam footprints.
//---------------------------------------------------------------------------- 
void SetTargets(Time& t, const string& target_name, Uvec& latb, Uvec& lonb, Uvec& Xtb, Uvec& Ytb, Uvec& Ztb, const int& key)
{

  Frame ftitan("IAU_TITAN",target_name);
  StateVector sc_state("sc_state");
  ftitan.ephemeris(sc_state,"Cassini",t,"NONE");

  for (unsigned int i_beam = 0; i_beam < 5 ;++i_beam)
    {
      Frame fbeam("CASSINI_RADAR_"+toStr(i_beam+1),"Cassini");

      double deg2rad = atan(1.)/45.;
      Uvar azi = double(key)*Uvar(2.7/2*(0.5-double(rand())/RAND_MAX), "rad")*deg2rad;
      Uvar ele = double(key)*Uvar(2.7/2*(0.5-double(rand())/RAND_MAX), "rad")*deg2rad;

      DirectionVector look("bore",fbeam,t,0,0,1);
      look.setAzimuthElevation(azi,ele);
      TargetGeom tg(t);
      tg.setState(sc_state);
      tg.setTarget(target_name,ftitan);
      tg.setLookDirection(look); 


      
      latb(i_beam) = tg.lat();
      lonb(i_beam) = tg.lon();
            
      Uvar X_intc = tg.surfaceIntercept()[PositionVector::X];
      Uvar Y_intc = tg.surfaceIntercept()[PositionVector::Y];
      Uvar Z_intc = tg.surfaceIntercept()[PositionVector::Z];
      
      Xtb(i_beam) = X_intc;
      Ytb(i_beam) = Y_intc;
      Ztb(i_beam) = Z_intc;
      
    }
}

//----------------------------------------------------------------------------
// Function Name: getBurstEcho
// Author: YZ
// Date:   10/27/03
// Function Description:
//      This function calculates the echo function for a entire burst, where
// (1) ieb is the IEB class.
// (2) epoch_time is the epoch_time.
// (3) t_x is the transmission time.
// (4) target_name is the target name, e.g., titan.
// (5) Rt is the array which stores the spacecraft location.
// (6) echo is the calculated burst echo.
//---------------------------------------------------------------------------- 
void getBurstEcho(Ieb& ieb, const Time& epoch_time, const Uvar& t_x, const string& target_name, const Uvec& Rt, Uvec& rw_time, Uvec& echo)
{
  Uvar pri = ieb.getPri();
  Uvar SR = ieb.getAdc();
  Uvar tau_p = ieb.getTaup();
  Uvar csf = ieb.getChirpStartFrequency();
  Uvar cfs = ieb.getCfs();
  Uvar csq = ieb.getCsq();
  Uvar csd = ieb.getCsd();
  Uvar rwd = ieb.getRwd();

  Uvar to_TA = t_x - epoch_time;

  Uvar rw_start_time = to_TA + rwd;
  unsigned int N_sample = ieb.getRadarDataPoints();
  unsigned int N_in_pri = ieb.getPointsInsidePRI();
  unsigned int N_pulse = (unsigned int)(N_sample/N_in_pri);
  for(unsigned int n_pulses = 0; n_pulses < N_pulse; n_pulses++)
    {
      
      Uvar rw_pulse_time = rw_start_time + n_pulses*pri;
      Time tp = epoch_time + rw_pulse_time;
      
      Uvec Rsc("Rsc",3);
      Uvec Vsc("Vsc",3);
      SCState(tp, target_name, Rsc, Vsc);
      
      for(unsigned int n_in_pri = 0; n_in_pri < N_in_pri; n_in_pri++)
	{
	  unsigned int n_point = n_pulses*N_in_pri + n_in_pri; 
	  Uvar rw_sample_time = n_point/SR + rw_start_time;		  
	  Uvar rw_local_time = rw_sample_time - rw_start_time - n_pulses*pri;
	  
	  Uvar Range = fitRange(Rsc, Vsc, Rt, rw_local_time);
	  
	  Uvar tr = rw_sample_time - to_TA;  
		  
	  Uvar phi = echo_phase(tr, Range, cfs, csd, csf, pri);
	  double PI = solid_pulse(tr, Range, pri, tau_p);
	  
	  rw_time(n_point) = n_point/SR;
	  echo(n_point) = cos(phi.getInUnits("rad")) * PI;	  
	}	          
    }   
}


//---------------------------------
//compute pulse segmentation
//--------------------------------
void computeStartEndIndices(Uvar& min_range,
			    const Uvar& max_range,
			    const Uvar& pri,
			    const Uvar& taup,
			    const Uvar& adc,
			    const unsigned int& pul,
			    const Uvar& rwd,
			    const unsigned int& Ndata_pri,
			    const unsigned int& Nradar_data,
			    unsigned int& index_start,
			    unsigned int& index_end,
			    unsigned int& Npulse){
  //consistency check
  unsigned int N = (unsigned int) round_double( (pri*adc).getInUnits(""));
  if(N != Ndata_pri) ErrorMessage("Number of data inside PRI does not match to computed one based  on adc and pri").throwMe();

  //segmentation : starting point
  unsigned int pulse_skip=0;

  Uvar dt= 2.0*min_range/speed_light - rwd;
  double  r_first_index = (dt*adc).getInUnits("");
  int first_index =int( floor( r_first_index));

  double r_remain= r_first_index - float(first_index);
  if(r_remain<0) {
   
    ErrorMessage("Error in computing remainder").throwMe();
  }
  Uvar  delta_range= -r_remain/adc*speed_light/2;
  min_range= min_range + delta_range;

  while(first_index<0){
    first_index += (int) Ndata_pri;
    pulse_skip++;//miss a pulse from beginning
  }
  
  //end index
  dt = double(pul-1)*pri + taup;
  dt += 2.0* max_range/speed_light;
  dt -= rwd;
  int last_index = round_double( (dt*adc).getInUnits(""));
  while(last_index > (int)Nradar_data){
    last_index -=(int) Ndata_pri;
    pulse_skip++;//miss a pulse from end
  }

  if(first_index <0 || last_index <0) ErrorMessage("start or last index is less than 0").throwMe();
  index_start =(unsigned int) first_index;
  index_end = (unsigned int) last_index;

  if(index_end < index_start) ErrorMessage("last end is smaller than index start").throwMe();
  if(pulse_skip > pul) ErrorMessage("number of  lost pulses are larger than transmitte pulses").throwMe();
  Npulse =  pul - pulse_skip;
}
