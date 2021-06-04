
static const char rcs_id_pdsvolume_c[] =
  "@(#) $Id: PDSVolume.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "DebugInfo.h"
#include "Error.h"
#include "PDSLabel.h"
#include "PDSVolume.h"
#include "Units.h"
#include "Frame.h"
#include "Utils.h"

using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::ifstream;
using std::internal;
using std::ios;
using std::left;
using std::ofstream;
using std::right;
using std::setfill;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

IndexTable::IndexTable()
: file_name_(""), path_name_(""), data_set_id_(""),
  start_time_(""), stop_time_(""), target_name_(""),
  minimum_latitude_(""), maximum_latitude_(""),
  westernmost_longitude_(""), easternmost_longitude_(""),
  look_direction_("")
{
}

IndexTable::IndexTable(string& filename, string& mode)
: stream_out_(filename.c_str(), ios::out),
  file_name_(""), path_name_(""), data_set_id_(""),
  start_time_(""), stop_time_(""), target_name_(""),
  minimum_latitude_(""), maximum_latitude_(""),
  westernmost_longitude_(""), easternmost_longitude_(""),
  look_direction_("")
{
}

//----------------------------------------------------------------------
// createFileHeader()
//
// Define the values to be written to the first record of the index table.
//----------------------------------------------------------------------
void IndexTable::createFileHeader()
{
  file_name_ = "FILE_NAME";
  path_name_ = "PATH_NAME";
  data_set_id_ = "DATA_SET_ID";
  start_time_ = "START_TIME";
  stop_time_ = "STOP_TIME";
  target_name_ = "TARGET_NAME";
  minimum_latitude_ = "MINIMUM_LATITUDE";
  maximum_latitude_ = "MAXIMUM_LATITUDE";
  westernmost_longitude_ = "WESTERNMOST_LONGITUDE";
  easternmost_longitude_ = "EASTERNMOST_LONGITUDE";
  look_direction_ = "LOOK_DIRECTION";
}

//----------------------------------------------------------------------
// readIndexData()
//
// Read the fields to be written into the index table from the current
// data product's PDS label.
//----------------------------------------------------------------------
void IndexTable::readIndexData(const string &filename)
{
  stringstream s_stream(ios::out);
  float minimum_latitude;
  float maximum_latitude;
  float westernmost_longitude;
  float easternmost_longitude;

  //
  // Read the file's PDS label.
  //
  PDSLabel file_header(filename);
  if (file_header.label() == EMPTY_STRING)
  {
    ErrorMessage e("Error:  PDS label not found in " + filename);
    e.throwMe();
  }

  //
  // Read the values that will go into the index table.
  //
  if (!file_header.getString("DATA_SET_ID", &data_set_id_))
  {
    ErrorMessage e("Keyword DATA_SET_ID not found in " + filename);
    e.throwMe();
  }

  if (!file_header.getString("START_TIME", &start_time_))
  {
    ErrorMessage e("Keyword START_TIME not found in " + filename);
    e.throwMe();
  }

  if (!file_header.getString("STOP_TIME", &stop_time_))
  {
    ErrorMessage e("Keyword STOP_TIME not found in " + filename);
    e.throwMe();
  }

  if (!file_header.getString("TARGET_NAME", &target_name_))
  {
    ErrorMessage e("Keyword TARGET_NAME not found in " + filename);
    e.throwMe();
  }

  if (!file_header.getString("LOOK_DIRECTION", &look_direction_))
  {
    ErrorMessage e("Keyword LOOK_DIRECTION not found in " + filename);
    e.throwMe();
  }

  //
  // If the file is a BIDR data product, read its latitude and longitude
  // bounds from the label.  If the file is not a BIDR data product,
  // then the latitude and longitude bounds are nonexistent; set the
  // bounds to a known invalid value instead.
  //
  if (data_set_id_.find(BIDR_BASE_DATA_SET_ID) != string::npos)
  {
    if (!file_header.getNumeric("MINIMUM_LATITUDE", &minimum_latitude))
    {
      ErrorMessage e("Keyword MINIMUM_LATITUDE not found in " + filename);
      e.throwMe();
    }

    if (!file_header.getNumeric("MAXIMUM_LATITUDE", &maximum_latitude))
    {
      ErrorMessage e("Keyword MAXIMUM_LATITUDE not found in " + filename);
      e.throwMe();
    }

    if (!file_header.getNumeric("WESTERNMOST_LONGITUDE", &westernmost_longitude))
    {
      ErrorMessage e("Keyword WESTERNMOST_LONGITUDE not found in " + filename);
      e.throwMe();
    }

    if (!file_header.getNumeric("EASTERNMOST_LONGITUDE", &easternmost_longitude))
    {
      ErrorMessage e("Keyword EASTERNMOST_LONGITUDE not found in " + filename);
      e.throwMe();
    }
  }
  else
  {
    minimum_latitude = INVALID_LATLON;
    maximum_latitude = INVALID_LATLON;
    westernmost_longitude = INVALID_LATLON;
    easternmost_longitude = INVALID_LATLON;
  }

  //
  // Round the latitude and longitude values to the nearest degree
  // and express each value as a string.
  //
  s_stream.seekp(0, ios::beg);
  s_stream << std::fixed << setprecision(0) << std::noshowpoint <<
    minimum_latitude;
  minimum_latitude_ = s_stream.str();

  s_stream.seekp(0, ios::beg);
  s_stream << std::fixed << setprecision(0) << std::noshowpoint <<
    maximum_latitude;
  maximum_latitude_ = s_stream.str();

  s_stream.seekp(0, ios::beg);
  s_stream << std::fixed << setprecision(0) << std::noshowpoint <<
    westernmost_longitude;
  westernmost_longitude_ = s_stream.str();

  s_stream.seekp(0, ios::beg);
  s_stream << std::fixed << setprecision(0) << std::noshowpoint <<
    easternmost_longitude;
  easternmost_longitude_ = s_stream.str();
}

//----------------------------------------------------------------------
// printRecord()
//----------------------------------------------------------------------
void IndexTable::printRecord()
{
  stream_out_ <<
    writeIndexField(file_name_, 21, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(path_name_, 9, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(data_set_id_, 28, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(start_time_, 21, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(stop_time_, 21, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(target_name_, 30, LEFT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(minimum_latitude_, 16, RIGHT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(maximum_latitude_, 16, RIGHT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(westernmost_longitude_, 21, RIGHT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(easternmost_longitude_, 21, RIGHT_JUSTIFY, !IS_LAST_ITEM) <<
    writeIndexField(look_direction_, 5, LEFT_JUSTIFY, IS_LAST_ITEM) <<
    TERMINATOR;
}

//----------------------------------------------------------------------
// writeIndexField()
//----------------------------------------------------------------------
string IndexTable::writeIndexField(string item, unsigned int width,
  int alignment, bool is_last_item)
{
  stringstream s_stream(ios::out);

  // Output only the first "width" characters of item.
  if (item.length() > width)
  {
    item = item.substr(0, width);
  }

  if (alignment == LEFT_JUSTIFY)
  {
    s_stream << CSV_QUOTE << setfill(' ') << setw(width) << std::left <<
      std::internal << item << CSV_QUOTE;
  }
  else if (alignment == RIGHT_JUSTIFY)
  {
    s_stream << CSV_QUOTE << setfill(' ') << setw(width) << std::right <<
      std::internal << item << CSV_QUOTE;
  }
  else
  {
    s_stream << CSV_QUOTE << setw(width) << item << CSV_QUOTE;
  }
  if (!is_last_item) s_stream << CSV_SEPARATOR;
  return s_stream.str();
}

