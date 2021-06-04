#ifndef PDSVOLUME_H
#define PDSVOLUME_H

static const char rcs_id_pdsvolume_h[] =
  "@(#) $Id: PDSVolume.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <iostream>
#include <string>

using std::ofstream;
using std::string;

#define IS_LAST_ITEM 1
#define CSV_QUOTE string("\"")
#define CSV_SEPARATOR string(",")
#define LEFT_JUSTIFY 1
#define RIGHT_JUSTIFY 2

// These four strings should be the actual BODP and BIDR data set IDs
// minus the version identifiers.
#define SBDR_BASE_DATA_SET_ID string("CO-V/E/J/S-RADAR-3-SBDR")
#define LBDR_BASE_DATA_SET_ID string("CO-V/E/J/S-RADAR-3-LBDR")
#define ABDR_BASE_DATA_SET_ID string("CO-SSA-RADAR-3-ABDR")
#define BIDR_BASE_DATA_SET_ID string("CO-SSA-RADAR-5-BIDR")

#define INVALID_LATLON -1000.0

class IndexTable {

  public:

  IndexTable();
  IndexTable(string& filename, string& mode);
  void createFileHeader();
  void readIndexData(const string &filename);
  void printRecord();
  string writeIndexField(string item, unsigned int width,
    int alignment, bool is_last_item);

  private:

  ofstream stream_out_;

  // The following fields are updated for each input file
  FileMgr file_in_;
  string file_name_;
  string path_name_;
  string data_set_id_;
  string start_time_;
  string stop_time_;
  string target_name_;
  string minimum_latitude_;  // degrees
  string maximum_latitude_;  // degrees
  string westernmost_longitude_; // degrees, positive-east (not -west)
  string easternmost_longitude_; // degrees, positive-east (not -west)
  string look_direction_;
};

#endif
