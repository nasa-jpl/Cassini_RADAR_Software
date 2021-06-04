//--------------------------------
//
//-------------------------------


#ifndef Radiometer_H
#define Radiometer_H
class Radiometer_coldsky;
class Radiometer;
#include "Array.h"
#include "Radiometer.h"

class Radiometer_coldsky
{
 public:
  Radiometer_coldsky(const Uvar& coldsky_count,
		     const Uvar& coldsky_temperature,
		     const Uvar& rl_count,
		     const Uvar& rl_temperature,
		     const Uvar& waveguide_temperature,
		     const double& loss)
    {//assume
      // cold sky temp and cold sky count :measured
      // resistive load temp and count: measured
      // lumped waveguide temperature and loss
      Ncoldsky_ = coldsky_count;
      Tcoldsky_ = coldsky_temperature;
      Nrl_ = rl_count;
      Trl_ = rl_temperature;
      waveguide_temp_=waveguide_temperature;
      loss_ = loss;
      count_ratio_ = Nrl_/Ncoldsky_;
    }
  
  void computeGainandReceiverTemp(Uvar& gain, Uvar& receiver_temp)
    {
      
      receiver_temp = Trl_ ;
      receiver_temp -= count_ratio_*(Tcoldsky_/loss_ + (1.0-1.0/loss_)*waveguide_temp_);
      receiver_temp/= (count_ratio_-1.0);
      gain = Nrl_/(receiver_temp + Trl_);
    }
 private: 
  Uvar Ncoldsky_;
  Uvar Tcoldsky_;
  Uvar Nrl_;
  Uvar Trl_;
  Uvar waveguide_temp_;
  double loss_;
  Uvar count_ratio_;
};

class Radiometer
{
 public:
  Radiometer(const Uvar& receiver_noise_temperature,
	     const Uvar& waveguide_temperature,
	     const double& loss)
    {
      receiver_ = receiver_noise_temperature;
      Tp_ = waveguide_temperature;
      loss_ = loss;
      Tsys_ = loss*receiver_ + (loss-1.0)*Tp_;
    }
  Uvar getTsys()
    {
      return(Tsys_);
    }

  void computeGainTbrightness(const Uvar& rl_count,
			      const Uvar& rlo_temp,
			      const Uvar& radio_count,
			      Uvar& gain, 
			      Uvar& Tbright)
    {
      gain = rl_count/(receiver_ + rlo_temp);
      Tbright = radio_count/gain*loss_ - Tsys_;
    }
 private:
  Uvar receiver_;
  Uvar Tp_;
  double loss_;
  Uvar Tsys_;
};
#endif


