//==============================================================//
// Copyright (C) 1997-2001, California Institute of Technology. //
// U.S. Goverment sponsorship acknowledged                      //
//==============================================================//

//==================================================================//
// Author:  Bryan Stiles   Created 9/22/97                          //
//==================================================================/

//==================================================================//
// CLASSES                                                          //
//    GenericDist, Uniform, Gaussian, GTC, AttDist                  //
//==================================================================//

//==================================================================//
// Functions: SeedFromClock                                         //
//==================================================================//

static const char rcs_id_distributions_c[] =
        "@(#) $Id: Distributions.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include"Distributions.h"
#include "Constants.h"

//=============//
// GenericDist //
//=============//

GenericDist::~GenericDist()
{
    return;
}

//=====================//
// GenericTimelessDist //
//=====================//

GenericTimelessDist::~GenericTimelessDist()
{
    return;
}

//--------------------------------//
// GenericTimelessDist::GetNumber //
//--------------------------------//

float
GenericTimelessDist::GetNumber(
    double  timex)
{
    if (timex < 0.0)    //bogus check to keep compiler quiet
        return(0.0);
    else
        return(GetNumber());
}

//=====//
// RNG //
//=====//

RNG::RNG(
    int  seed)
{
    SetSeed(seed);
    return;
}

RNG::RNG()
{
    SetRandomSeed();
    return;
}

RNG::~RNG()
{
  return;
}

//--------------//
// RNG::SetSeed //
//--------------//

void
RNG::SetSeed(
    int  seed)
{
    if (seed > 0)
        _seed = -seed;
    else if (seed == 0)
        _seed = -1;
    else
        _seed = seed;
    _Init();
}

void
RNG::SetRandomSeed()
{
    SetSeed(lrand48() * lrand48());
    return;
}

double
RNG::GetDouble()
{
    int j = (int)(97.0 * _output / RNG_M);
    _output = _tab[j];
    _seed = (RNG_IA * _seed + RNG_IC) % RNG_M;
    _tab[j] = _seed;
    return((double)(_output)/RNG_M);
}

void
RNG::_Init()
{
    if((_seed=(RNG_IC - _seed) % RNG_M) < 0)
        _seed = -_seed;
    for(int j = 0; j < 97; j++)
    {
        _seed = (RNG_IA * _seed + RNG_IC) % RNG_M;
        _tab[j] = _seed;
    }
    _seed = (RNG_IA * _seed + RNG_IC) % RNG_M;
    _output = _seed;
    return;
}

//=========//
// Uniform //
//=========//

Uniform::Uniform()
:    _radius(0.0), _mean(0.0)
{
    return;
}

Uniform::Uniform(float radius, float mean)
{
    _radius=radius;
    _mean=mean;
    return;
}

Uniform::~Uniform(){
    return;
}

//====================//
// Uniform::GetNumber //
//====================//

float
Uniform::GetNumber()
{
    float num;
    num=(float)(2*_radius*(_rng.GetDouble()-0.5)+_mean);
    return(num);
}

float Uniform::GetRadius(){
  return(_radius);
}

void Uniform::SetRadius(float r){
  _radius=r;
}

float Uniform::GetMean(){
  return(_mean);
}

void Uniform::SetMean(float m){
  _mean=m;
}

void Uniform::SetSeed(int seed){
  _rng.SetSeed(seed);
}

//============================//
// Gaussian                   //
//============================//

Gaussian::Gaussian()
:    _variance(0.0), _mean(0.0)
{
    return;
}

Gaussian::Gaussian(float variance, float mean)
{
    _variance=variance;
    _mean=mean;
    return;
}


Gaussian::~Gaussian()
{
    return;
}

float
Gaussian::GetNumber()
{
    float num;
    double v1, v2, r, fac;
        do {
             v1=2.0*_rng.GetDouble()-1.0;
             v2=2.0*_rng.GetDouble()-1.0;
             r = v1*v1+v2*v2;
        } while (r >= 1.0 || r == 0.0);
        fac=sqrt((float) -2.0*log(r)/r);
        num=(float)v2*fac;
        num=num*sqrt(_variance) + _mean;
    return(num);
}

float Gaussian::GetVariance(){
  return(_variance);
}

int Gaussian::SetVariance(float v){
  if(v < 0.0) return(0);
  _variance=v;
  return(1);
}

float Gaussian::GetMean(){
  return(_mean);
}

int Gaussian::SetMean(float m){
  _mean=m;
  return(1);
}

void Gaussian::SetSeed(int seed){
  _rng.SetSeed(seed);
}

//============================//
// Gamma                      //
//============================//

Gamma::Gamma()
:    _variance(0.0), _mean(0.0)
{
    return;
}

Gamma::Gamma(float variance, float mean)
{
    _variance=variance;
    _mean=mean;
    return;
}


Gamma::~Gamma()
{
    return;
}

float
Gamma::GetNumber()
{

  double v1, v2, r, mag, fac, lambda, x;

  // Convert mean and variance to alternate r,lambda representation.
  lambda = _mean / _variance;
  r = _mean * lambda;

  // Set coefficients of comparison function to ensure that it is always
  // greater than the gamma pdf, but as close as possible.
  double a0 = 2.0*sqrt(_variance);
  double x0 = (r-1.0)/lambda;
  double c0 = exp(log(lambda) - gammln(r) + (r-1.0)*log(r-1.0) - (r-1.0));

  do
  {  // This loop generates a gamma distributed deviate.
    do
    {  // This loop generates an x uniformly distributed under the comparison
       // function.  Integral of comparison function is:
       // a0*tan(pi*u) + x0 (u is uniform deviate).
      do
      {  // This loop generates the tangent of a random angle
        v1=2.0*_rng.GetDouble()-1.0;
        v2=2.0*_rng.GetDouble()-1.0;
        mag = v1*v1+v2*v2;
      } while (mag > 1.0 || mag == 0.0);
      double y = v2/v1; // y = tan(pi * uniform random deviate)

      // Compute x-coordinate of uniformly distributed point under the
      // comparison pdf.
      x = a0*y + x0;
    } while (x <= 0.0);

    // fac is the ratio of the gamma pdf to the comparison pdf.
    fac = exp(r*log(lambda) - gammln(r) + (r-1.0)*log(x) -lambda*x) *
          (1.0 + (x - x0)*(x - x0)/a0/a0) / c0;
  } while (_rng.GetDouble() > fac);

  return((float)x);
}

float Gamma::GetVariance(){
  return(_variance);
}

int Gamma::SetVariance(float v)
{
  if(v < 0.0) return(0);
  _variance=v;
  return(1);
}

float Gamma::GetMean(){
  return(_mean);
}

int Gamma::SetMean(float m){
  _mean=m;
  return(1);
}

void Gamma::SetSeed(int seed){
  _rng.SetSeed(seed);
}


//==================================//
// Random Velocity                  //
//==================================//

RandomVelocity::RandomVelocity()
{
    _sample_period=1.0;
    _radius=1.0;
        _mean=0.0;
    _noise= NULL;
    _position=0.0;
    _time=0.0;
    _velocity=0.0;
     return;
}

RandomVelocity::RandomVelocity(GenericTimelessDist* noise, float sample_period,
    float radius, float mean)
{
    _noise=noise;
    _sample_period=sample_period;
       _radius=radius;
    _mean=mean;
    _position=mean;
    _time=0.0;
    _velocity=noise->GetNumber();
    while(fabs(_position-_mean+_velocity*_sample_period) > _radius){
        _velocity=noise->GetNumber();
    }
    return;
}

RandomVelocity::~RandomVelocity(){
    return;
}

//================================//
// RandomVelocity::GetNumber      //
//================================//

float RandomVelocity::GetNumber(double timex){
    if (timex < 0.0){
     fprintf(stderr,"Fatal Error produced by RandomVelocity::GetNumber\n");
     fprintf(stderr,"Parameter time may not be negative.\n");
     exit(1);
    }
    if (timex < _time){
     fprintf(stderr,"Fatal Error produced by RandomVelocity::GetNumber\n");
     fprintf(stderr,"Parameter time may not decrease between \n");
     fprintf(stderr,"consecutive calls to the method. \n");
     exit(1);
    }
    while (timex >= _time + _sample_period){
      _position+=_velocity*_sample_period;
      _time+=_sample_period;
       _velocity=_noise->GetNumber();
      while(fabs(_position-_mean+_velocity*_sample_period) > _radius){
        _velocity=_noise->GetNumber();
      }

    }
    return(_position+(timex-_time)*_velocity);
}

//==================================================================//
// Class                                                            //
//  TimeCorrelatedGaussian                                          //
//==================================================================//

TimeCorrelatedGaussian::TimeCorrelatedGaussian()
  : _previousTime(0.0), _previousOutput(0.0), _correlationLength(0.0), _mean(0.0)
{
  return;
}

TimeCorrelatedGaussian::~TimeCorrelatedGaussian(){
  return;
}

int TimeCorrelatedGaussian::Initialize(){
  _previousOutput=Uncorrelated.GetNumber();
   return(1);
}

//-----------------------------------//
// TimeCorrelatedGaussian::GetNumber //
//-----------------------------------//

float
TimeCorrelatedGaussian::GetNumber(
    double  timex)
{
    /******* BIAS ONLY CASE *************/
    if (Uncorrelated.GetVariance() == 0.0)
    {
        return(_mean);
    }

    /******** Error Condition ****************/
    if (timex < _previousTime)
    {
        fprintf(stderr,
            "TimeCorrelatedGaussian requires monotonically increasing time\n");
        fprintf(stderr, "  (%.3f -> %.3f)\n", _previousTime, timex);
        exit(1);
    }

    /********** Uncorrelated case *************/
    if (_correlationLength == 0.0)
    {
        return(Uncorrelated.GetNumber() + _mean);
    }

    /******* Normal Mode *******************/
    float retval = exp( -(timex - _previousTime) / _correlationLength);
    retval = retval * _previousOutput + sqrt(1 - retval*retval)
        * Uncorrelated.GetNumber();
    _previousTime = timex;
    _previousOutput = retval;
    return(retval + _mean);
}

int TimeCorrelatedGaussian::SetVariance(float variance){
  if(! Uncorrelated.SetVariance(variance)) return(0);
  return(1);
}

int TimeCorrelatedGaussian::SetMean(float mean){
  _mean=mean;
  return(1);
}

void TimeCorrelatedGaussian::SetSeed(int seed){
  Uncorrelated.SetSeed(seed);
}

int TimeCorrelatedGaussian::SetCorrelationLength(float corrlength){
  if(corrlength < 0.0) return(0);
  _correlationLength=corrlength;
  return(1);
}

//==================================//
// AttDist                          //
//==================================//

AttDist::AttDist()
{
    return;
}


AttDist::~AttDist()
{
    return;
}

//==================================//
// SeedFromClock                    //
//==================================//
void
SeedFromClock()
{
  struct timeval now;
  gettimeofday(&now,NULL);
  srand48(now.tv_sec);
}

//======================================================================//
// gammln                                                               //
// Returns the natural logarithm of the gamma function of the argument. //
//======================================================================//

double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
