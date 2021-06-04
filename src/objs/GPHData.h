//============================================================================
// GPHData.h
//
// This file contains the GPHData class declarations. This class reads
// *.gph Cassini temperature telemetry files
//============================================================================
#ifndef GPHDATA_H
#define GPHDATA_H

#include "Units.h"
#include "Io.h"
#include "Error.h"
#include "Array.h"
#include "Time.h"
#include <string>

class GPHData
{
 public:

  GPHData();
  GPHData(const string& filename); 
  Uvar interpolate(Time t, const Uvar& interpolation_valid_time) throw(ErrorMessage);
  Uvar interpolate(Time t) throw(ErrorMessage);
  unsigned int recordCount() throw(ErrorMessage);
  void read(const string& filename); 
 protected:
  void read() throw(ErrorMessage);
  void readHeader() throw(ErrorMessage,FileMgr::IoError);
  string filename_;
  FileMgr file_;
  unsigned int records_counted_;
  bool file_read_;
  bool header_handled_;
  Uvec temp_;
  Array1D<Time> time_;
};
#endif








