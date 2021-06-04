//============================================================================
// Config.h
//
// This file contains the declaration for class Config.
// Config provides support for configuration files that contain
// space delimited keyword, value pairs.
//============================================================================

#ifndef Config_H
#define Config_H

#include <string>
#include <map>
#include <list>
#include <set>
#include <fstream>
#include <sstream>
#include "Units.h"
#include "Error.h"

using std::string;
using std::map;
using std::set;
using std::ifstream;
using std::ofstream;

//-------------------------------
// Configuration list management
//-------------------------------

class Config
  {
  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;
  typedef map<string,string> strings;
  typedef map<string,Uvar> numbers;

  //------------------
  // Automatic testing
  //------------------

  static bool selfTest();

  //--------------
  // Constructors
  //--------------

  Config(const string& filename);
  Config(const string& caller_id, const string& filename);
  ~Config();

  //-----------
  // Operators
  //-----------

  Uvar operator[](const string& s);

  //-----------
  // Predicates
  //-----------

  bool keywordExists(const string& s) const;
  bool numberKeywordExists(const string& s) const;
  bool stringKeywordExists(const string& s) const;
  bool isDeprecated(const string& s) const;
  bool isDisabled(const string& s) const;

  //----------------
  // File I/O
  //----------------

  void readFile(const string& filename);
  void writeAutoCommentFile(const string& filename);

  //----------------
  // Access methods
  //----------------

  string str(const string& s);
  Uvar num(const string& s);
  int getInt(const string& s);
  Uvar getTime(const string& s);
  double getDouble(const string& s);
  void removeKeyword(const string& keyword);
  void addKeyword(const string& keyword, const string& value);
  void addKeyword(const string& keyword, const Uvar& value, const string& ustr);
  void suffixSet(const string& keyword, numbers& subset) const;
  void suffixSet(const string& keyword, strings& subset) const;
  void suffixSet(const string& keyword,
    list<string>& keyword_list, list<string>& suffix_list) const;
  Uvar convertWestLongitude(const string& s);
  string filename() const;

  //-------------------------
  // Other support methods
  //-------------------------

  static void suffixSort(const list<string>& suffix_list,
    list<string>& keyword_list);

  //---------------------
  // debugging support
  //---------------------

  static void numbersShow(const string& name, numbers& nn);
  static void showStringList(const string& name, const list<string>& slist);

  private:

  enum modeE {none,block_mode,keyword_mode,value_mode,keyword_cont,value_cont};

  //---------------------------
  // Private support methods
  //---------------------------

  void expandKeywordList(list<string>& keywords) const;
  void readSpecialKeywords();
  void writeLines(ofstream& output_stream, const list<string>& lines);
  void writeKeywordLines(ofstream& output_stream,
    const list<string>& keywords);
  void writeCallerKeywords(ostream& output_stream) const;
  string lineMerge(const list<string>& lines) const;

  //-----------------------
  // Private data elements
  //-----------------------

  bool auto_comments_;  // indicates if auto-comments are to be generated
  string filename_;     // Name of config file that goes with a Config object
  string caller_id_;    // Name of calling program/function for auto-comments
  string master_cfg_;   // Name of master config file

  strings string_list_;
  numbers number_list_;
  strings line_list_;
  strings line_list_copy_;
  list<string> user_header_comments_;
  map<string,list<string> > user_comments_;
  mutable set<string> caller_keywords_;
  mutable map<string,list<string> > deprecated_keywords_;
  mutable map<string,list<string> > disabled_keywords_;
  bool interactive_mode_on_; // Allow missing keywords to be entered 
                             // interactively?

  bool update_file_;         // Update the config file?
  };

#endif
