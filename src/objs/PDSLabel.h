#ifndef PDSLABEL_H
#define PDSLABEL_H

static const char rcs_id_pdslabel_h[] =
  "@(#) $Id: PDSLabel.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <list>
#include <map>
#include <sstream>
#include <string>
#include "Config.h"

using std::ios;
using std::list;
using std::string;
using std::stringstream;

#define CARET string("^")
#define EQUALS string("=")
#define STRING_QUOTE string("\"")
#define OPEN_SEQUENCE string("(")
#define CLOSE_SEQUENCE string(")")
#define OPEN_UNITS_BRACKET string("<")
#define CLOSE_UNITS_BRACKET string(">")
#define SEPARATOR string(",")
#define CARRIAGE_RETURN string("\r")
#define TERMINATOR string("\r\n")
#define SPACE string(" ")
#define SPACE_REPLACEMENT string("_")
#define EMPTY_STRING string("")
#define NO_UNITS EMPTY_STRING
#define NO_UNITS_TEXT EMPTY_STRING
#define TIMESTAMP_PLACEHOLDER string("YYYY-DOYTHH:MM:SS.sss")

#define IS_ODL_STRING 1
#define FIXED 0
#define SCIENTIFIC 1
#define INDENT_STEP 2   // Number of spaces per indentation level
#define LINEBUFLEN 80   // Bytes
#define MAXLINELEN 2000 // Bytes

class KeyValPair {

  public:

  KeyValPair();
  KeyValPair(string keyword, string value, int is_ODL_string);
  KeyValPair(string keyword, double *value, int n_elements, int precision,
    string units);
  KeyValPair(string keyword, float *value, int n_elements, int precision,
    string units);
  KeyValPair(string keyword, int *value, int n_elements, string units);
  KeyValPair(string keyword, double value, int precision);
  KeyValPair(string keyword, double value, int precision, int notation);
  KeyValPair(string keyword, double value, int precision, string units);
  KeyValPair(string keyword, double value, int precision, int width,
    string units);
  KeyValPair(string keyword, float value, int precision);
  KeyValPair(string keyword, float value, int precision, int notation);
  KeyValPair(string keyword, float value, int precision, string units);
  KeyValPair(string keyword, float value, int precision, int width,
    string units);
  KeyValPair(string keyword, int value);
  KeyValPair(string keyword, int value, int width, string units);

  string Keyword();
  string Value();

  protected:

  string keyword_;
  string value_;
};

typedef list<KeyValPair> KeyValSeq;

typedef map<const string, int> PDSCharValLimit;

KeyValSeq updateValues(KeyValSeq old_list, KeyValSeq new_list);

class PDSLabel{

  public:

  // Constructors
  PDSLabel();
  PDSLabel(const string& filename);
  PDSLabel(FileMgr& f);
  PDSLabel(KeyValSeq keyvals, int record_length, int data_records,
    string first_data_record_ptr);

  // Access methods for data members
  int fileRecords();
  string label();
  int labelLength();
  int labelRecords();
  int recordLength();

  // Utilities
  static int badCharValLength(string& keyword, string& value);
  static string replaceSpaces(string s);
  static string toUpper(string& str);
  static void setPDSCharValLimits();
  int createListFromLabel(string& label, KeyValSeq& list);

  // Access methods for values of keywords in label
  template <class T>
  int getNumeric(const string& keyword, T *value)
  {
    string KEYWORD = keyword;
    string buf = EMPTY_STRING;
    string line;
    string LINE;
    size_t index;
    stringstream buf_stream(label_, ios::in);

    KEYWORD = toUpper(KEYWORD);
    *value = 0;
    while (getline(buf_stream, line))
    {
      index = 0;
      LINE = toUpper(line);
      if ((index = LINE.find(keyword, index)) != string::npos)
      {
	index += keyword.length();
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

    // Read the value
    index += EQUALS.length();
    stringstream buf_tmp(LINE.substr(index), ios::in);
    if (buf_tmp >> *value)
    {
      return (sizeof(*value));
    }
    else
    {
      return 0;
    }
  }
  int getString(const string& keyword, string *value);
  int getVector(const string& keyword, int n, double *values);
  static PDSCharValLimit PDSCharValLimits;

  private:

  // Methods for writing labels
  int keywordFieldWidth(KeyValSeq keyvals);
  int writeStatements(KeyValSeq keyvals);

  // Methods for reading labels
  int getStringStream(const string& filename, const string& keyword,
    const int search_limit, string *value);
  int getStringStream(FileMgr& f, const string& keyword,
    const int search_limit, string *value);
  int getIntStream(const string& filename, const string& keyword,
    const int search_limit, int *value);
  int getIntStream(FileMgr& f, const string& keyword,
    const int search_limit, int *value);
  int getLabelStream(const string& filename, const int length,
    string *label);
  int getLabelStream(FileMgr& f, const int length,
    string *label);

  string label_;
  int keyword_field_width_;
  int file_records_;
  int label_length_;
  int label_records_;
  int record_length_;
};

#endif
