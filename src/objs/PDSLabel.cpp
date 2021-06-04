
static const char rcs_id_pdslabel_c[] =
  "@(#) $Id: PDSLabel.cpp,v 11.6 2012/09/27 20:01:58 richw Exp $";

#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "DebugInfo.h"
#include "Error.h"
#include "PDSLabel.h"
#include "Units.h"
#include "Frame.h"
#include "Utils.h"

using std::cerr;
using std::cout;
using std::endl;
using std::fixed;
using std::ifstream;
using std::ios;
using std::left;
using std::scientific;
using std::setfill;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::uppercase;

PDSCharValLimit PDSLabel::PDSCharValLimits;

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

KeyValPair::KeyValPair()
: keyword_(""), value_("")
{
}

// Constructor for string value; can also be used for symbol values, but
// symbol delimiters aren't supported
KeyValPair::KeyValPair(string keyword, string value, int is_ODL_string)
{
  keyword_ = keyword;

  // Check that value satisfies length constraints in PDS Data Dictionary
  int limit;
  if ((limit = PDSLabel::badCharValLength(keyword, value)))
  {
    string text = "PDS keyword " + keyword + " has value " + "\"" +
      value + "\"" + " with illegal length " + toStr(value.length()) +
      "; maximum is " + toStr(limit) + ".  Exiting.";
    ErrorMessage e(text);
    e.throwMe();
  }

  if (is_ODL_string)
  {
    value_ = STRING_QUOTE + value + STRING_QUOTE;
  }
  else
  {
    value_ = value;
  }
}

// Constructor for arrays of doubles (or a single double)
// Units are only supported if n_elements = 1
KeyValPair::KeyValPair(string keyword, double *value, int n_elements,
  int precision, string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;

  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  if (n_elements > 1)
  {
    buf_stream << OPEN_SEQUENCE + SPACE;
    for (int i = 0; i < n_elements; i++)
    {
      buf_stream << setprecision(precision) << std::fixed << value[i];
      if (i < n_elements-1) buf_stream << SEPARATOR;
      buf_stream << SPACE;
    }
    buf_stream << CLOSE_SEQUENCE;
  }
  else
  {
    buf_stream << setprecision(precision) << std::fixed << *value;
  }

  value_ = buf_stream.str();
  // Currently, units are only supported for single values.
  if (units != EMPTY_STRING && n_elements == 1) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for arrays of floats (or a single float)
// Units are only supported if n_elements = 1
KeyValPair::KeyValPair(string keyword, float *value, int n_elements,
  int precision, string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;

  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  if (n_elements > 1)
  {
    buf_stream << OPEN_SEQUENCE + SPACE;
    for (int i = 0; i < n_elements; i++)
    {
      buf_stream << setprecision(precision) << std::fixed << value[i];
      if (i < n_elements-1) buf_stream << SEPARATOR;
      buf_stream << SPACE;
    }
    buf_stream << CLOSE_SEQUENCE;
  }
  else
  {
    buf_stream << setprecision(precision) << std::fixed << *value;
  }

  value_ = buf_stream.str();
  // Currently, this method handles units only for scalar inputs
  if (units != EMPTY_STRING && n_elements == 1) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for arrays of integers (or a single integer)
// Units are only supported if n_elements = 1
KeyValPair::KeyValPair(string keyword, int *value, int n_elements,
  string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;

  if (n_elements > 1)
  {
    buf_stream << OPEN_SEQUENCE + SPACE;
    for (int i = 0; i < n_elements; i++)
    {
      buf_stream << value[i];
      if (i < n_elements-1) buf_stream << SEPARATOR;
      buf_stream << SPACE;
    }
    buf_stream << CLOSE_SEQUENCE;
  }
  else
  {
    buf_stream << *value;
  }

  value_ = buf_stream.str();
  // Currently, this method handles units only for scalar inputs
  if (units != EMPTY_STRING && n_elements == 1) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for double value without units
KeyValPair::KeyValPair(string keyword, double value, int precision)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed << value;
  value_ = buf_stream.str();
}

// Constructor for double value with support for scientific notation,
// no units
KeyValPair::KeyValPair(string keyword, double value, int precision,
  int notation)
{
  stringstream buf_stream(ios::out);

  keyword_ = keyword;
  if (notation == SCIENTIFIC)
  {
    buf_stream << scientific << setprecision(precision) << uppercase << value;
  }
  else
  {
    // This compiler interprets the setprecision() setting as the total
    // number of digits on both sides of the decimal point, not the number
    // of digits after the point.  To work around this problem, set the
    // fixed-format option before calling setprecision().  (It's not clear
    // why this works, but it does.)
    buf_stream << setprecision(precision) << std::fixed << value;
  }
  value_ = buf_stream.str();
}

// Constructor for double value with units
KeyValPair::KeyValPair(string keyword, double value, int precision,
  string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed << value;
  value_ = buf_stream.str();
  if (units != EMPTY_STRING) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for double value with units and fixed width
KeyValPair::KeyValPair(string keyword, double value, int precision,
  int width, string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed <<
    setw(width) << setfill(' ') << value;
  value_ = buf_stream.str();
  if (units != EMPTY_STRING) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for float value without units
KeyValPair::KeyValPair(string keyword, float value, int precision)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed << value;
  value_ = buf_stream.str();
}

// Constructor for float value with support for scientific notation,
// no units
KeyValPair::KeyValPair(string keyword, float value, int precision,
  int notation)
{
  stringstream buf_stream(ios::out);

  keyword_ = keyword;
  if (notation == SCIENTIFIC)
  {
    buf_stream << scientific << setprecision(precision) << uppercase << value;
  }
  else
  {
    // This compiler interprets the setprecision() setting as the total
    // number of digits on both sides of the decimal point, not the number
    // of digits after the point.  To work around this problem, set the
    // fixed-format option before calling setprecision().  (It's not clear
    // why this works, but it does.)
    buf_stream << setprecision(precision) << std::fixed << value;
  }
  value_ = buf_stream.str();
}

// Constructor for float value with units
KeyValPair::KeyValPair(string keyword, float value, int precision,
  string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed << value;
  value_ = buf_stream.str();
  if (units != EMPTY_STRING) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for float value with units and fixed width
KeyValPair::KeyValPair(string keyword, float value, int precision,
  int width, string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // This compiler interprets the setprecision() setting as the total
  // number of digits on both sides of the decimal point, not the number
  // of digits after the point.  To work around this problem, set the
  // fixed-format option before calling setprecision().  (It's not clear
  // why this works, but it does.)
  buf_stream << setprecision(precision) << std::fixed <<
    setw(width) << setfill(' ') << value;
  value_ = buf_stream.str();
  if (units != EMPTY_STRING) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Constructor for integer value without units
KeyValPair::KeyValPair(string keyword, int value)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  buf_stream << value;
  value_ = buf_stream.str();
}

// Constructor for integer value with a specified field width
// and possibly units
KeyValPair::KeyValPair(string keyword, int value, int width, string units)
{
  stringstream buf_stream(ios::app|ios::in|ios::out);

  keyword_ = keyword;
  // Since the input value is an integer, use '0' as the fill character.
  buf_stream << std::setw(width) << std::setfill('0') << value;
  value_ = buf_stream.str();
  if (units != EMPTY_STRING) {
    value_ = value_ + SPACE + OPEN_UNITS_BRACKET +
      units + CLOSE_UNITS_BRACKET;
  }
}

// Public access functions

string KeyValPair::Keyword()
{
  return keyword_;
}

string KeyValPair::Value()
{
  return value_;
}

//----------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------

PDSLabel::PDSLabel()
: label_(EMPTY_STRING), keyword_field_width_(0), file_records_(0),
  label_length_(0), label_records_(0), record_length_(0)
{
}

//----------------------------------------------------------------------
// Constructor used for reading PDS labels
//----------------------------------------------------------------------
PDSLabel::PDSLabel(const string& filename)
: label_(EMPTY_STRING), keyword_field_width_(0), file_records_(0),
  label_length_(0), label_records_(0), record_length_(0)
{
  int search_limit = 10;
  string keyword;
  string s_value;

  DebugInfo dbg("PDSLabel::PDSLabel(string)");

  if (!getStringStream(filename, "PDS_VERSION_ID", search_limit, &s_value))
  {
    // The file does not contain a (valid) PDS label.
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(string):" << endl;
      dbg.file << "No PDS label found in " << filename << endl;
    }
    return;
  }
  // This file appears to contain a PDS label.  Expand the search radius
  // and find the record length (in bytes) and the number of label records
  // in the file.  Since we don't yet know the record length, use the same
  // search method as was used for the PDS marker.
  search_limit *= 4;
  if (!getIntStream(filename, "RECORD_BYTES", search_limit, &record_length_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(string):" << endl;
      dbg.file << "Keyword RECORD_BYTES not found in " << filename << endl;
    }
    return;
  }
  if (!getIntStream(filename, "FILE_RECORDS", search_limit, &file_records_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(string):" << endl;
      dbg.file << "Keyword FILE_RECORDS not found in " << filename << endl;
    }
    return;
  }
  if (!getIntStream(filename, "LABEL_RECORDS", search_limit, &label_records_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(string):" << endl;
      dbg.file << "Keyword LABEL_RECORDS not found in " << filename << endl;
    }
    return;
  }
  // Read the entire label.
  label_length_ = record_length_ * label_records_;
  int nread;
  if ((nread = getLabelStream(filename, label_length_, &label_)) !=
       label_length_)
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(string):" << endl;
      dbg.file << "Warning:  read only " << nread << " of " <<
	label_length_ << " bytes from " << filename << endl;
      dbg.file << "PDS label is" << endl;
      dbg.file << label_ << endl;
    }
  }
  if (dbg.level)
  {
    dbg.file << "PDSLabel::PDSLabel(string):" << endl;
    dbg.file << "PDS label is" << endl;
    dbg.file << label_ << endl;
  }
}

//----------------------------------------------------------------------
// Constructor used for reading PDS labels
//----------------------------------------------------------------------
PDSLabel::PDSLabel(FileMgr& f)
: label_(EMPTY_STRING), keyword_field_width_(0), file_records_(0),
  label_length_(0), label_records_(0), record_length_(0)
{
  int search_limit = 10;
  string keyword;
  string s_value;

  DebugInfo dbg("PDSLabel::PDSLabel(FileMgr)");

  if (f.name() == EMPTY_STRING)
  {
    ErrorMessage e("PDSLabel::PDSLabel(FileMgr):  input file not open");
    e.throwMe();
  }
  if (!getStringStream(f, "PDS_VERSION_ID", search_limit, &s_value))
  {
    // The file does not contain a (valid) PDS label.
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
      dbg.file << "No PDS label found in input file" << endl;
    }
    return;
  }
  // This file appears to contain a PDS label.  Expand the search radius
  // and find the record length (in bytes) and the number of label records
  // in the file.  Since we don't yet know the record length, use the same
  // search method as was used for the PDS marker.
  search_limit *= 4;
  if (!getIntStream(f, "RECORD_BYTES", search_limit, &record_length_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
      dbg.file << "Keyword RECORD_BYTES not found in input file" << endl;
    }
    return;
  }
  if (!getIntStream(f, "FILE_RECORDS", search_limit, &file_records_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
      dbg.file << "Keyword FILE_RECORDS not found in input file" << endl;
    }
    return;
  }
  if (!getIntStream(f, "LABEL_RECORDS", search_limit, &label_records_))
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
      dbg.file << "Keyword LABEL_RECORDS not found in input file" << endl;
    }
    return;
  }
  // Read the entire label.
  label_length_ = record_length_ * label_records_;
  int nread;
  if ((nread = getLabelStream(f, label_length_, &label_)) != label_length_)
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
      dbg.file << "Warning:  read only " << nread << " of " <<
	label_length_ << " bytes from input file" << endl;
      dbg.file << "PDS label is" << endl;
      dbg.file << label_ << endl;
    }
  }
  if (dbg.level)
  {
    dbg.file << "PDSLabel::PDSLabel(FileMgr):" << endl;
    dbg.file << "PDS label is" << endl;
    dbg.file << label_ << endl;
  }
}

//----------------------------------------------------------------------
// Constructor used for writing PDS labels
//----------------------------------------------------------------------
PDSLabel::PDSLabel(KeyValSeq keyvals, int record_length, int data_records,
  string first_data_record_ptr):
  file_records_(data_records+1),
  label_records_(1),
  record_length_(record_length)
{
  int n_records;

  if (record_length <= 0)
  {
    label_ = EMPTY_STRING;
    return;
  }
  keyword_field_width_ = keywordFieldWidth(keyvals);
  n_records = writeStatements(keyvals);

  bool records_updated=false; // forces at least one record number update
                              // needed in case record number has changed
                              // e.g. subsampling performed or some such
                              // without this the data_records parameter
                              // does not get used
  while (n_records != label_records_ || !records_updated)
  {
    label_records_ = n_records;
    file_records_ = data_records + n_records;
    KeyValSeq updates;
    updates.push_back(KeyValPair("FILE_RECORDS", file_records_));
    updates.push_back(KeyValPair("LABEL_RECORDS", label_records_));
    if (first_data_record_ptr != EMPTY_STRING)
    {
      updates.push_back(KeyValPair(first_data_record_ptr, label_records_+1));
    }
    keyvals = updateValues(keyvals, updates);
    n_records = writeStatements(keyvals);
    records_updated=true;
  }
  label_length_ = label_records_ * record_length_;
}

//----------------------------------------------------------------------
// fileRecords()
//
// Return the total number of records in this file.
//----------------------------------------------------------------------
int PDSLabel::fileRecords()
{
  return file_records_;
}

//----------------------------------------------------------------------
// label()
//
// Return the (entire) PDS label.
//----------------------------------------------------------------------
string PDSLabel::label()
{
  return label_;
}

//----------------------------------------------------------------------
// labelLength()
//
// Return the length of the PDS label in bytes.
//----------------------------------------------------------------------
int PDSLabel::labelLength()
{
  return label_length_;
}

//----------------------------------------------------------------------
// labelRecords()
//
// Return the number of file records used by the PDS label.  This number
// equals the length of the PDS label in bytes divided by the length of
// one record in bytes.  (Each record is assumed to have the same length.)
//----------------------------------------------------------------------
int PDSLabel::labelRecords()
{
  return label_records_;
}

//----------------------------------------------------------------------
// recordLength()
//
// Return the length of one file record in bytes.
//----------------------------------------------------------------------
int PDSLabel::recordLength()
{
  return record_length_;
}

//----------------------------------------------------------------------
// badCharValLength()
//
// Check the length of a (character-valued) PDS keyword against the maximum
// length specified in the PDS Data Dictionary.  If the limit exists and
// is exceeded, return the limit; otherwise return 0.  Currently, this method
// only checks maximum lengths, since all the character-valued keywords used
// in the BIDR and L1B output files have limits only on their maximum lengths.
//----------------------------------------------------------------------
int PDSLabel::badCharValLength(string& keyword, string& value)
{
  PDSCharValLimit::iterator itr = PDSCharValLimits.find(keyword);
  if (itr == PDSCharValLimits.end())
  {
    // Keyword not found ==> no limit
    return 0;
  }

  unsigned int maxlen = (*itr).second;
  if (value.length() > maxlen) return maxlen;

  return 0;
}

//----------------------------------------------------------------------
// replaceSpaces()
//
// Replace all instances of the SPACE_REPLACEMENT character (or string)
// in the input string with the SPACE character (or string).
//----------------------------------------------------------------------
string PDSLabel::replaceSpaces(string s)
{
  // Replace all instances of SPACE_REPLACEMENT with SPACE in input string
  size_t i = s.find(SPACE_REPLACEMENT);
  while (i != string::npos)
  {
    s.replace(i, SPACE_REPLACEMENT.length(), SPACE);
    i = s.find(SPACE_REPLACEMENT, i+SPACE.length());
  }
  return s;
}

//----------------------------------------------------------------------
// setPDSCharValLimits()
//
// Create the reference table of maximum lengths for the character-valued
// PDS keywords in the BIDR and BODP SISes.  The lengths are expressed
// in bytes and are taken from the PDS Data Dictionary.  Keywords with no
// length limits are not included.
//----------------------------------------------------------------------
void PDSLabel::setPDSCharValLimits()
{
  // List keywords in alphabetical order to help avoid duplication.
  // Duplicate keys will supersede earlier assignments.
  //
  // Use only upper case for the keys in the table.  Although the PDS
  // standard allows lower-case letters in identifiers (keywords), our
  // keywords are all in upper-case.

  PDSCharValLimits["DATA_SET_ID"] = 40;
  PDSCharValLimits["DATA_SET_NAME"] = 60;
  PDSCharValLimits["DATA_TYPE"] = 30;
  PDSCharValLimits["INSTRUMENT_HOST_ID"] = 11;
  PDSCharValLimits["INSTRUMENT_HOST_NAME"] = 120;
  PDSCharValLimits["INSTRUMENT_ID"] = 12;
  PDSCharValLimits["INSTRUMENT_NAME"] = 60;
  PDSCharValLimits["INTERCHANGE_FORMAT"] = 6;
  PDSCharValLimits["LOOK_DIRECTION"] = 6;
  PDSCharValLimits["MAP_PROJECTION_TYPE"] = 28;
  PDSCharValLimits["MISSION_NAME"] = 60;
  PDSCharValLimits["MISSION_PHASE_NAME"] = 30;
  PDSCharValLimits["NAME"] = 61;
  PDSCharValLimits["PDS_VERSION_ID"] = 6;
  PDSCharValLimits["POSITIVE_LONGITUDE_DIRECTION"] = 4;
  PDSCharValLimits["PRODUCER_FULL_NAME"] = 60;
  PDSCharValLimits["PRODUCER_ID"] = 20;
  PDSCharValLimits["PRODUCER_INSTITUTION_NAME"] = 60;
  PDSCharValLimits["PRODUCT_ID"] = 40;
  PDSCharValLimits["PRODUCT_VERSION_ID"] = 12;
  PDSCharValLimits["RECORD_TYPE"] = 20;
  PDSCharValLimits["SAMPLE_TYPE"] = 30;
  PDSCharValLimits["SOFTWARE_VERSION_ID"] = 20;
  PDSCharValLimits["SOURCE_PRODUCT_ID"] = 40;
  PDSCharValLimits["SPACECRAFT_CLOCK_START_COUNT"] = 30;
  PDSCharValLimits["SPACECRAFT_CLOCK_STOP_COUNT"] = 30;
  PDSCharValLimits["TARGET_NAME"] = 120;
  PDSCharValLimits["UNIT"] = 40;
}

//----------------------------------------------------------------------
// toUpper()
//
// Convert all lower-case characters in a string to upper-case.
//----------------------------------------------------------------------
string PDSLabel::toUpper(string& str)
{
  int i = 0;
  int len = str.length();
  for (i = 0; i < len; i++)
  {
    str[i] = toupper(str[i]);
  }
  return(str);
}

//----------------------------------------------------------------------
// createListFromLabel()
//
// Given a string containing the text of a PDS label, construct a
// KeyValSeq structure that can be used to regenerate that label.
//----------------------------------------------------------------------
int PDSLabel::createListFromLabel(string& label, KeyValSeq& list)
{
  DebugInfo dbg("PDSLabel::createListFromLabel");

  if (label == EMPTY_STRING)
  {
    if (dbg.level)
    {
      dbg.file << "PDSLabel::createListFromLabel:  " <<
	"Input consists of an empty string" << endl;
    }
    return 0;
  }

  // Assume that the input list is empty.  If it isn't, then new
  // keyword-value pairs will be appended to the list.

  // For each line of text:
  //
  // Search for an equals sign.
  //
  // If an equals sign is found, then set the keyword to be the first
  // through the last non-blank characters to the left of the equals sign
  // and set the value to be the first through the last non-blank characters
  // to the right of the equals sign.  (Note that the code does not assume
  // that an equals sign is padded with spaces on either side.)
  //
  // If no equals sign is found, then the line is either a comment or a
  // blank line.  In either case, set the keyword to be the contents of the
  // line and leave the value set to an empty string.
  //
  // The line terminator is assumed to be a newline (\n).  (For the Cassini
  // PDS labels, the PDS label line terminator is actually a carriage return
  // followed by a newline, but the getline() function used below supports
  // only single-character line terminators.)  Trailing carriage returns
  // are stripped from the end of each line before adding the line to the list.

  // This routine returns 1 on success, 0 otherwise.

  stringstream buf_stream(label, ios::in);
  string line;
  size_t eq_index;
  size_t start, end;
  size_t index;
  string keyword;
  string value;

  if (dbg.level)
  {
    dbg.file << "PDSLabel::createListFromLabel:  " << endl;
    dbg.file << "Elements of list are " << endl;
  }
  while (getline(buf_stream, line))
  {
    eq_index = 0;
    keyword = EMPTY_STRING;
    value = EMPTY_STRING;
    int length = line.length();
    if ((eq_index = line.find(EQUALS)) != string::npos)
    {
      // Equals sign found; find the keyword
      start = line.find_first_not_of(SPACE, 0);
      end = line.find_first_of(SPACE+EQUALS, start+1) - 1;
      if (start == eq_index || start == string::npos ||
          end == eq_index || end == string::npos ||
	  start > end)
      {
	// Keyword is missing?
	if (dbg.level)
	{
	  dbg.file << "PDSLabel::createListFromLabel:  " <<
	    "Unable to parse statement" << endl;
	  dbg.file << line << endl;
	}
        return 0;
      }
      keyword = line.substr(start, end-start+1);

      // Find the value
      start = line.find_first_not_of(SPACE, eq_index+1);
      end = line.find_last_not_of(SPACE, length-1);
      if (start == eq_index || start == string::npos ||
          end == eq_index || end == string::npos ||
	  start > end)
      {
	// Value is missing?
	if (dbg.level)
	{
	  dbg.file << "PDSLabel::createListFromLabel:  " <<
	    "Unable to parse statement" << endl;
	  dbg.file << line << endl;
	}
        return 0;
      }
      value = line.substr(start, end-start+1);
      // Strip off the carriage return, if any
      if ((index = value.find(CARRIAGE_RETURN)) != string::npos)
      {
        value.erase(index, CARRIAGE_RETURN.length());
      }
    }
    else
    {
      // No equals sign found
      keyword = line;
      // Strip off the carriage return, if any
      if ((index = keyword.find(CARRIAGE_RETURN)) != string::npos)
      {
        keyword.erase(index, CARRIAGE_RETURN.length());
      }
      // Strip off any trailing spaces
      end = line.find_last_not_of(SPACE, keyword.length()-1);
      keyword = keyword.substr(0, end+1);
    }
    list.push_back(KeyValPair(keyword, value, !IS_ODL_STRING));
    if (dbg.level)
    {
      dbg.file << "Keyword = ";
      if (keyword != EMPTY_STRING)
      {
        dbg.file << keyword;
      }
      else
      {
        dbg.file << "(EMPTY STRING)";
      }
      dbg.file << " Value = ";
      if (value != EMPTY_STRING)
      {
        dbg.file << value << endl;
      }
      else
      {
        dbg.file << "(EMPTY STRING)" << endl;
      }
    }
  }
  return 1;
}

//----------------------------------------------------------------------
// getString()
//
// Fetch the value of a keyword in a PDS label.  The value is assumed
// to be a string and is returned in the argument *value.  If successful,
// this function returns the length of the string value; otherwise, it
// returns 0.
//----------------------------------------------------------------------
int PDSLabel::getString(const string& keyword, string *value)
{
  string KEYWORD = keyword;
  string buf = EMPTY_STRING;
  string line;
  string LINE;
  size_t index;
  size_t start_index, end_index;
  stringstream buf_stream(label_, ios::in);

  KEYWORD = toUpper(KEYWORD);
  *value = EMPTY_STRING;
  while (getline(buf_stream, line))
  {
    LINE = toUpper(line);
    index = 0;
    if ((index = LINE.find(KEYWORD, index)) != string::npos)
    {
      index += KEYWORD.length();
      break;
    }
  }
  if (index == string::npos)
  {
    // Keyword not found
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = LINE.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the string value.
  index += EQUALS.length();
  start_index = LINE.find(STRING_QUOTE, index);
  if (start_index != string::npos)
  {
    // Opening string delimiter found.  If the closing string delimiter
    // exists, return the characters enclosed within the delimiters.
    // If only the opening string delimiter exists, the statement is
    // malformed.
    if ((end_index = LINE.find(STRING_QUOTE, start_index+STRING_QUOTE.length()))
      != string::npos)
    {
      if (end_index > start_index+1)
      {
	*value = LINE.substr(start_index+1, end_index-start_index-1);
	return (*value).length();
      }
    }
    else
      return 0;
  }
  // No string delimiters were found to the right of the equals sign.
  // Search for the first nonspace character instead.
  start_index = LINE.find_first_not_of(SPACE+TERMINATOR, index);
  if (start_index != string::npos)
  {
    // String value found.  The string should be terminated either
    // by a space or by end-of-record; otherwise, the statement is
    // malformed.
    end_index = LINE.find_first_of(SPACE+TERMINATOR, start_index+1);
    if (end_index != string::npos)
    {
      *value = LINE.substr(start_index, end_index-start_index);
    }
    else
      return 0;
  }
  return (*value).length();
}

//----------------------------------------------------------------------
// getVector()
//
// Fetch the components of a vector from a PDS label.  The values are
// returned in the array *values.  This function returns the number of
// values successfully read.  The vector in the label is assumed to have
// exactly n components in the format
//     keyword = ( x, y, ..., z )
// where x, y, z, and the other components are all numerical values and
// the opening sequence delimiter, closing sequence delimiter, and
// element separator are as shown.  Parsing of unit strings within the
// vector is not supported.
//----------------------------------------------------------------------
int PDSLabel::getVector(const string& keyword, int n, double *values)
{
  string KEYWORD = keyword;
  string buf = EMPTY_STRING;
  string line;
  string LINE;
  int i;
  int n_found = 0;
  size_t index;
  size_t start_index, end_index;
  stringstream buf_stream(label_, ios::in);

  KEYWORD = toUpper(KEYWORD);
  while (getline(buf_stream, line))
  {
    LINE = toUpper(line);
    index = 0;
    if ((index = LINE.find(KEYWORD, index)) != string::npos)
    {
      index += KEYWORD.length();
      break;
    }
  }
  if (index == string::npos)
  {
    // Keyword not found
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = LINE.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the string value.
  index += EQUALS.length();
  start_index = LINE.find(OPEN_SEQUENCE, index);
  if (start_index != string::npos)
  {
    // Opening vector delimiter found.  Locate the closing vector delimiter.
    // If only the opening vector delimiter exists, the statement is
    // malformed.
    start_index += OPEN_SEQUENCE.length();
    if ((end_index = LINE.find(CLOSE_SEQUENCE, start_index)) != string::npos)
    {
      // Search between start_index and end_index for the elements of the
      // vector, each of which is separated by a SEPARATOR.
      for (i = 0; i < n-1; i++)
      {
	index = LINE.find(SEPARATOR, start_index);
	if (index != string::npos && index < end_index)
	{
	  values[n_found++] =
	    strtod(LINE.substr(start_index, index-start_index).c_str(),
	      (char **) NULL);
	}
	start_index = index + SEPARATOR.length();
      }
      // Last element.  The vector should not contain any more separators.
      index = LINE.find(SEPARATOR, start_index);
      if (index == string::npos)
      {
	values[n_found++] =
	  strtod(LINE.substr(start_index, end_index-start_index).c_str(),
	    (char **) NULL);
      }
    }
    else
      return 0;
  }
  return n_found;
}

//---------------------------------------------------------------------
// keywordFieldWidth()
//
// Determine the field width to use for the keywords in this PDS label
// so that all the equals signs in the label are aligned (for better
// readability).  All keywords will be left-aligned and indented by
// INDENT_STEP spaces if they are enclosed by OBJECT statements.
// OBJECT statements may be nested.  The final field width is the
// smallest value wide enough to accommodate every keyword when
// indentation is taken into account.  Keywords without values are
// assumed to be either comments or END statements and are ignored.
// No checking is done on whether OBJECT and END_OBJECT statements
// match; also, GROUP statements are not supported.
//---------------------------------------------------------------------
int PDSLabel::keywordFieldWidth(KeyValSeq keyvals)
{
  KeyValSeq::iterator i;
  string indent;
  string keyword;
  string value;
  string indent_spaces(INDENT_STEP, ' ');
  int len;
  int maxlen;

  indent = EMPTY_STRING;
  maxlen = 0;
  // Make one pass through the list to calculate the smallest field
  // width that can accommodate all the keywords, assuming that statements
  // within OBJECTs will be indented.  The level of indentation will
  // change depending upon the number of nested objects.
  for (i = keyvals.begin(); i != keyvals.end(); i++) {
    keyword = (*i).Keyword();
    value = (*i).Value();
    if (value != EMPTY_STRING)
    {
      len = indent.length() + keyword.length();
      if (len > maxlen) maxlen = len;
    }
    if (keyword == "OBJECT") {
      indent.append(indent_spaces);
    }
    if (keyword == "END_OBJECT") {
      if (indent.length() >= INDENT_STEP) {
	indent.erase(0, INDENT_STEP);
      }
    }
  }
  return maxlen;
}

//---------------------------------------------------------------------
// updateValues()
//
// Update specified values in an existing list of keyword-value pairs.
// The order of keyword-value pairs is preserved.
//
// list - list of keyword-value pairs to be updated
// updates - list of keyword-value updates
//---------------------------------------------------------------------
KeyValSeq updateValues(KeyValSeq list, KeyValSeq updates)
{
  KeyValSeq::iterator i_new, i_old;
  KeyValSeq::iterator i_tmpnew, i_tmpold;

  // If a keyword in the existing list matches one in the list of updates,
  // then replace its value with the value from the list of updates.
  //
  // No checks are done for uniqueness; if a keyword appears more than once
  // in the list of updates, then the last occurrence supersedes the others.
  // If a keyword in the list of updates does not match any keywords in the
  // existing list, then the keyword is ignored.
  for (i_new = updates.begin(); i_new != updates.end(); i_new++)
  {
    for (i_old = list.begin(); i_old != list.end(); i_old++)
    {
      if ((*i_new).Keyword() == (*i_old).Keyword())
      {
        i_tmpnew = i_new;
        i_tmpnew++;
        list.insert(i_old, i_new, i_tmpnew);
        i_tmpold = i_old;
        i_old--;
        list.erase(i_tmpold);
      }
    }
  }
  return list;
}

//---------------------------------------------------------------------
// writeStatements()
//
// Write all the statements into the label, properly formatted.
//---------------------------------------------------------------------
int PDSLabel::writeStatements(KeyValSeq keyvals)
{
  int record_offset;
  KeyValSeq::iterator i;
  stringstream statement(ios::in|ios::out);
  stringstream buf(ios::in|ios::out);
  string keyword;
  string value;
  string indent;
  string indent_spaces(INDENT_STEP, ' ');
  int len;
  int n_records;
  int buf_position;
  unsigned int chars_written = 0;

  if (record_length_ <= 0)
  {
    // This statement should not be reachable if the check remains
    // in the constructor
    ErrorMessage e("PDSLabel::writeStatements:  record length <= 0");
    e.throwMe();
  }
  indent = EMPTY_STRING;
  record_offset = 0;
  for (i = keyvals.begin(); i != keyvals.end(); i++)
  {
    statement.str(EMPTY_STRING);
    keyword = (*i).Keyword();
    value = (*i).Value();
    if (keyword == "END_OBJECT")
    {
      if (indent.length() >= INDENT_STEP) {
	indent.erase(0, INDENT_STEP);
      }
    }
    if (value != EMPTY_STRING)
    {
      statement << std::left << setw(keyword_field_width_) <<
	(indent + keyword);
      statement << SPACE << EQUALS << SPACE << value;
    }
    else
    {
      statement << keyword;
    }
    statement << TERMINATOR;
    // Determine whether statement is to be written to current record
    // or to next
    len = statement.str().length();
    if (record_offset + len > record_length_ &&
      chars_written > TERMINATOR.length())
    {
      if (record_offset > 0)
      {
	// Move (rewrite) the terminator for the last statement
	// from its current position to the end of the record so that
	// the current statement will be written to a new record.
	//
	// Note:  the call to seekp() has been revised so that a positive
	// offset is passed as the first argument.  gcc 3.4.6 apparently
	// contains a bug that causes the buffer position to be set to -1
	// when a negative offset is passed relative to either ios::cur or
	// ios::end.
	//
	// buf.seekp(-(TERMINATOR.length()), ios::end);
	buf_position = buf.tellp();
	buf.seekp(buf_position-(TERMINATOR.length()), ios::beg);
	if (buf.fail())
	{
	  ErrorMessage e("PDSLabel::writeStatements:  seekp failed at end of record");
	  e.throwMe();
	}
	buf << string(record_length_ - record_offset, ' ');
	buf << TERMINATOR;
      }
      record_offset = 0;
    }
    if (keyword == "OBJECT") {
      indent.append(indent_spaces);
    }
    buf << statement.str();
    record_offset = (record_offset + len) % record_length_;
    chars_written += len;
  }
  // See above comment about seekp().
  // buf.seekp(-(TERMINATOR.length()), ios::end);
  buf_position = buf.tellp();
  buf.seekp(buf_position-(TERMINATOR.length()), ios::beg);
  if (buf.fail())
  {
    ErrorMessage e("PDSLabel::writeStatements:  seekp failed at end of label");
    e.throwMe();
  }
  buf << string(record_length_ - record_offset, ' ');
  buf << TERMINATOR;
//  cout << "Number of records = " << (buf.str().length() / record_length_) << endl;
  label_ = buf.str();
  // At this point, buf.str().length() will be an integer multiple of
  // record_length_
  n_records = buf.str().length() / record_length_;
  return n_records;
}

//----------------------------------------------------------------------
// getStringStream()
//
// Fetch the value of a given (PDS label) keyword from a file.  The value
// is assumed to be a string and is returned in the argument *value.  This
// function opens the file via the specified filename and searches only
// the first search_limit lines of the file.  If successful, this function
// returns the length of the string value; otherwise, it returns 0.
//----------------------------------------------------------------------
int PDSLabel::getStringStream(const string& filename, const string& keyword,
  const int search_limit, string *value)
{
  string buf = EMPTY_STRING;
  char linebuf[LINEBUFLEN];
  string line = EMPTY_STRING;
  size_t index;
  size_t start_index, end_index;
  ifstream filestream(filename.c_str(), ios::in);

  *value = EMPTY_STRING;
  if (!filestream)
  {
    ErrorMessage e("PDSLabel::getStringStream:  cannot open " + filename);
    e.throwMe();
  }
  int i = 0;
  while (filestream.getline(linebuf, LINEBUFLEN) && i < search_limit)
  {
    line = linebuf;
    index = 0;
    if ((index = line.find(keyword, index)) != string::npos)
    {
      index += keyword.length();
      break;
    }
    i++;
  }
  filestream.close();
  if (i == search_limit || index == string::npos)
  {
    // Since the keyword wasn't present in the first search_limit lines
    // of the file, assume that this file does not contain the keyword.
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = line.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the string value.
  index += EQUALS.length();
  start_index = line.find(STRING_QUOTE, index);
  if (start_index != string::npos)
  {
    // Opening string delimiter found.  If the closing string delimiter
    // exists, return the characters enclosed within the delimiters.
    // If only the opening string delimiter exists, the statement is
    // malformed.
    if ((end_index = line.find(STRING_QUOTE, start_index+STRING_QUOTE.length()))
      != string::npos)
    {
      if (end_index > start_index+1)
      {
	*value = line.substr(start_index+1, end_index-start_index-1);
	return (*value).length();
      }
    }
    else
      return 0;
  }
  // No string delimiters were found to the right of the equals sign.
  // Search for the first nonspace character instead.
  start_index = line.find_first_not_of(SPACE+TERMINATOR, index);
  if (start_index != string::npos)
  {
    // String value found.  The string should be terminated either
    // by a space or by end-of-record; otherwise, the statement is
    // malformed.
    end_index = line.find_first_of(SPACE+TERMINATOR, start_index+1);
    if (end_index != string::npos)
    {
      *value = line.substr(start_index, end_index-start_index);
    }
    else
      return 0;
  }
  return (*value).length();
}

//----------------------------------------------------------------------
// getStringStream()
//
// Fetch the value of a given (PDS label) keyword from a file.  The value
// is assumed to be a string and is returned in the argument *value.  This
// function accesses the file via the FileMgr argument, which should already
// have a pointer to an open file.  This function searches only the first
// search_limit lines of the file.  If successful, this function returns the
// length of the string value; otherwise, it returns 0.
//----------------------------------------------------------------------
int PDSLabel::getStringStream(FileMgr& f, const string& keyword,
  const int search_limit, string *value)
{
  string buf = EMPTY_STRING;
  string line;
  size_t index;
  size_t start_index, end_index;

  *value = EMPTY_STRING;
  if (f.name() == EMPTY_STRING)
  {
    ErrorMessage e("PDSLabel::getStringStream(FileMgr):  input file not open");
    e.throwMe();
  }
  int i = 0;
  line.resize(LINEBUFLEN);
  int pos = f.getPosition();
  while (!f.eof() && i < search_limit)
  {
    f.read(line);
    index = 0;
    if ((index = line.find(keyword, index)) != string::npos)
    {
      index += keyword.length();
      break;
    }
    i++;
  }
  f.setPosition(pos);
  if (i == search_limit || index == string::npos)
  {
    // Since the keyword wasn't present in the first search_limit*LINEBUFLEN
    // bytes of the file, assume that this file does not contain the keyword.
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = line.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the string value.
  index += EQUALS.length();
  start_index = line.find(STRING_QUOTE, index);
  if (start_index != string::npos)
  {
    // Opening string delimiter found.  If the closing string delimiter
    // exists, return the characters enclosed within the delimiters.
    // If only the opening string delimiter exists, the statement is
    // malformed.
    if ((end_index = line.find(STRING_QUOTE, start_index+STRING_QUOTE.length()))
      != string::npos)
    {
      if (end_index > start_index+1)
      {
	*value = line.substr(start_index+1, end_index-start_index-1);
	return (*value).length();
      }
    }
    else
      return 0;
  }
  // No string delimiters were found to the right of the equals sign.
  // Search for the first nonspace character instead.
  start_index = line.find_first_not_of(SPACE+TERMINATOR, index);
  if (start_index != string::npos)
  {
    // String value found.  The string should be terminated either
    // by a space or by end-of-record; otherwise, the statement is
    // malformed.
    end_index = line.find_first_of(SPACE+TERMINATOR, start_index+1);
    if (end_index != string::npos)
    {
      *value = line.substr(start_index, end_index-start_index);
    }
    else
      return 0;
  }
  return (*value).length();
}

//----------------------------------------------------------------------
// getIntStream()
//
// Fetch the value of a given (PDS label) keyword from a file.  The value
// is assumed to be an integer.  getIntStream() opens the file via the
// specified filename and searches only the first search_limit lines of
// the file.  If successful, getIntStream() returns the value; otherwise,
// it returns 0.
//----------------------------------------------------------------------
int PDSLabel::getIntStream(const string& filename, const string& keyword,
  const int search_limit, int *value)
{
  string buf = EMPTY_STRING;
  string line = EMPTY_STRING;
  size_t index;
  ifstream filestream(filename.c_str(), ios::in);

  *value = 0;
  if (!filestream)
  {
    ErrorMessage e("PDSLabel::getIntStream:  cannot open " + filename);
    e.throwMe();
  }
  int i = 0;
  while (getline(filestream, line) && i < search_limit)
  {
    index = 0;
    if ((index = line.find(keyword, index)) != string::npos)
    {
      index += keyword.length();
      break;
    }
    i++;
  }
  filestream.close();
  if (i == search_limit || index == string::npos)
  {
    // Since the keyword wasn't present in the first search_limit lines
    // of the file, assume that this file does not contain the keyword.
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = line.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the integer value.
  index += EQUALS.length();
  stringstream buf_stream(line.substr(index), ios::app|ios::in);
  if (buf_stream >> *value)
  {
    return (*value);
  }
  else
  {
    return 0;
  }
}

//----------------------------------------------------------------------
// getIntStream()
//
// Fetch the value of a given (PDS label) keyword from a file.  The value
// is assumed to be an integer.  getIntStream() accesses the file via the
// FileMgr argument, which should already have a pointer to an open file.
// This function searches only the first search_limit lines of the file.
// If successful, getIntStream() returns the value; otherwise, it returns 0.
//----------------------------------------------------------------------
int PDSLabel::getIntStream(FileMgr& f, const string& keyword,
  const int search_limit, int *value)
{
  string buf = EMPTY_STRING;
  string line;
  size_t index;

  *value = 0;
  if (f.name() == EMPTY_STRING)
  {
    ErrorMessage e("PDSLabel::getIntStream(FileMgr):  input file not open");
    e.throwMe();
  }
  int i = 0;
  int pos = f.getPosition();
  line.resize(MAXLINELEN);
  while (!f.eof() && i < search_limit)
  {
    // Assume that line terminators will be present in the input rather
    // than attempting to read a fixed number of characters, which will
    // cause the line buffer to be misaligned with respect to the actual
    // PDS statements.  This assumption should be safe as long as the
    // calling routine has confirmed that the PDS label exists.
    f.readLine(line);
    index = 0;
    if ((index = line.find(keyword, index)) != string::npos)
    {
      index += keyword.length();
      break;
    }
    i++;
  }
  f.setPosition(pos);
  if (i == search_limit || index == string::npos)
  {
    // Since the keyword wasn't present in the first search_limit lines
    // of the file, assume that this file does not contain the keyword.
    return 0;
  }
  // The keyword is on the current line.  Find the equals sign.
  if ((index = line.find(EQUALS, index)) == string::npos)
  {
    // Malformed statement -- equals sign missing
    return 0;
  }

  // Read the integer value.
  index += EQUALS.length();
  stringstream buf_stream(line.substr(index), ios::app|ios::in);
  if (buf_stream >> *value)
  {
    return *value;
  }
  else
  {
    return 0;
  }
}

int PDSLabel::getLabelStream(const string& filename, const int length,
  string *label)
{
  string buf = EMPTY_STRING;
  ifstream filestream(filename.c_str(), ios::in);
  stringstream buf_stream(ios::out);
  char ch;

  if (!filestream)
  {
    ErrorMessage e("PDSLabel::getLabelStream:  cannot open " + filename);
    e.throwMe();
  }
  int i = 0;
  while (i < length)
  {
    filestream.get(ch);
    buf_stream << ch;
    if (filestream.eof())
    {
      break;
    }
    i++;
  }
  filestream.close();
  *label = buf_stream.str();
  return (buf_stream.str().length());
}

//----------------------------------------------------------------------
// getLabelStream()
//
// Fetch a PDS label from a file.  getLabelStream() accesses the file via
// the FileMgr argument, which should already have a pointer to an open file.
// This function reads up to "length" characters from the file and returns
// the label in the *label argument.  The function returns the number of
// characters read.
//----------------------------------------------------------------------
int PDSLabel::getLabelStream(FileMgr& f, const int length,
  string *label)
{
  string buf = EMPTY_STRING;
  stringstream buf_stream(ios::out);
  char ch;

  if (f.name() == EMPTY_STRING)
  {
    ErrorMessage e("PDSLabel::getLabelStream(FileMgr):  input file not open");
    e.throwMe();
  }
  int i = 0;
  int pos = f.getPosition();
  while (i < length)
  {
    f.read(ch);
    buf_stream << ch;
    if (f.eof())
    {
      break;
    }
    i++;
  }
  f.setPosition(pos);
  *label = buf_stream.str();
  return (buf_stream.str().length());
}

