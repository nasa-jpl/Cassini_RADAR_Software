#include"CAGTable.h"
#include<math.h>
#include "Error.h"

// Construct Cassini Radar Calibrated Attenuator Gain (CAG) Table_
CAGTable::CAGTable(){

  // The values in the table_ were calculated in pre-launch tests and
  // documented in RFES End Item Data Package Vol 10 of 13.
  // Pages 26-28 of document LOG/RFES-0485/ALS/ESM

  // The values thus obtained are incomplete for two reasons
  // (1) Multiple attenuator paths can lead to the same nominal attenuator
  //     setting. Standard commanding procedure eliminates this possibly but 
  //     it is unclear if the test procedure used the same standard.
  // (2) All possible attenuation setting were not exercised in the test.
  //     Since a different setting utilizes a different set of attenuators, it
  //     is not correct to interpolate between settings. For unexercised
  //     setting we have interpolated for lack of a better strategy. We have
  //     further assumed that the 0 dB setting correspond to exactly 0 dB.
  //     The table_ entries so estimated are flagged as EST. 

  table_[0] = 0.00; // EST
  table_[1] = 0.80; // EST
  table_[2] = 1.60;
  table_[6] = 5.42;
  table_[10]= 9.37;
  table_[14]=13.27;
  table_[18]=17.08;
  table_[22]=21.10;
  table_[26]=25.12;
  table_[30]=29.07;
  table_[34]=32.92;
  table_[38]=36.95;
  table_[42]=41.08;
  table_[46]=45.10;
  table_[50]=48.95;
  table_[54]=52.96;
  table_[58]=57.06;
  table_[62]=61.08;
  table_[66]=64.96;
  table_[70]=68.79; 
  table_[74]=72.48;


  // Interpolating other entries for lack of anything better to do
  for(int c=3;c<74;c++){
    int m=(c-2)%4;

    if(m==0) continue; // skip known entries 

    // compute weighting coeficient and indices of nearest valid entries
    float f=float(m)/4.0;
    int last_idx=c-m;
    int next_idx=last_idx+4;

    // linearly interpolate
    table_[c]= table_[last_idx]*(1-f)+table_[next_idx]*f; // EST
  }

  // convert all values from dB loss to linear scale gain
  for(int c=0;c<74;c++) table_[c]=pow(10,-0.1*table_[c]);

}

float 
CAGTable::gainLinearScale(int nominal_loss_in_db)
{
  if(nominal_loss_in_db<0 || nominal_loss_in_db> 74){
    ErrorMessage e("CAGTable:gainLinearScale nominal loss out of range");
    e.throwMe();
  }
  return table_[nominal_loss_in_db];
}

 int
CAGTable::findGaindB(const float& gain_value)
  {
    //table_: from 0 to 74,    
    //find nearest one
    if(gain_value>1) ErrorMessage("attenuation gain should be less than 1").throwMe();
    if(gain_value<0) ErrorMessage("attenuation gain should be positive ").throwMe();
    int gain_db_target=0;
    bool find_gain_db=false;
    for( int i=0; i< 73;++i){
      if(gain_value < table_[i] && gain_value>= table_[i+1]){
	find_gain_db = true;
	gain_db_target=i;
      }
    }
    if(!find_gain_db) ErrorMessage("could not find db number for a given gain v alue").throwMe();

    //if(gain_db_target%2!=0) gain_db_target++;//make even number
    //pick the closest one
    // not done yet, maybe after we got all the cal data

    return(gain_db_target);
    
  }
