//=============================================================================
// TextMgr.h
//
// This file contains the TextMgr class declarations.
// The TextMgr class provides support for integrating C++ programs and
// documents (especially TeX documents).  A TextMgr object knows how to read
// a C++ source file looking for special codes in the comments.  Document
// source (such as TeX source) is embedded in comments.
//
// Interface summary:
//
// TextMgr methods:
//
//   Construction:
//
//   TextMgr(output_filename); Construct to build indicated output file
//
//   Setup:
//
//   setSpecialString(str);    Indicate the string flagging embedded data
//   setCommentSyncStr(str);   Indicate the string identifying text in comment
//
//   I/O:
//
//   lineOut(filename,id_str); Write out one line of text from comment line
//
//   Variable handling:
//
//   addVar(varname,var); Add string rep. of var to make available for output
//   getVar(varname);     Return string rep. of var
//
//
//=============================================================================

#ifndef TEXTMGR_H
#define TEXTMGR_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include "Error.h"
//#include "Units.h"

using std::string;
using std::map;

class TextMgr
  {
  public:

  typedef std::ostringstream  OSTRINGSTREAM;
//  typedef UnitVar<double> Uvar;
  enum modeE {text, source};

  //--------------
  // Constructors
  //--------------

  TextMgr(const std::string& output_filename);

  //--------
  // Setup
  //--------

  void setSpecialStr(const string& special_str);
  void setCommentSyncStr(const string& comment_sync_str);

  //--------
  // I/O
  //--------

  void lineOut(const string& input_filename, const string& id_str)
    throw(ErrorMessage);

  //------------------
  // Variable handling
  //------------------

  void addVar(const string& varname, const string& var) throw();
//  void addVar(const string& varname, const Uvar& var);
  string getVar(const string& varname) throw(ErrorMessage);

  private:

  //-------------------
  // Internal variables
  //-------------------

  string output_filename_;
  string special_str_;
  string comment_sync_str_;
  modeE mode_;

  // Map to hold variables and associated printing string for output.
  map<string,string> variables_;

  };

#endif
