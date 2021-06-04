#ifndef CAG_TABLE_H
#define CAG_TABLE_H

// Calibrated Attenuator Gain Table class
// This class is used to construct and access a constant table of
// calibrated attenuator settings

// The table is indexed by nominal attenuator setting (loss) in dB
// (0,1,2,3 ... 74).
// The output value is the linear scale calibrated attenuator gain
// (values between 0 and 1).

class CAGTable{

 public:
  // constructor copies constant (hard-coded) values into table and
  // converts them from dB loss to linear scale gain
  // the values in question were obtained from the Final Performance Test 
  // (FPT) values in the TPR/RFES/0454/ALS/CSN test procedure. This test
  // was performed 12/9/95 by I Chialastri. Values are for ambient 
  // temperature.

  CAGTable();

  // Indexes the table 
  float gainLinearScale(int nominal_loss_in_db);
  //find gain db setting for desired gain value
  int findGaindB(const float& gain_value_);

 private:
  float table_[75];
};

#endif
