//----------------------------------------------------------------------------
// Config.cpp
//
// This file contains method definitions for the Config class.
// Config provides support for configuration files that contain
// space delimited keyword, value pairs.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char config_c[] =
  "@(#) $Id: Config.cpp,v 11.7 2016/09/20 17:10:13 richw Exp $";

#include <stdlib.h>
#include <string>
#include <fstream>
#include <map>
#include <list>
#include <set>
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Constants.h"
#include "Utils.h"
#include "DebugInfo.h"

using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::cerr;
using std::streampos;
using std::cin;

//----------------------
// Class Config Methods
//----------------------

//------------------
// Automatic testing
//------------------

bool Config::selfTest()
  {
  ofstream cfg_file("test_cfg_file.cfg");
  if (!cfg_file)
    {
     ErrorMessage e("Can't open test_cfg_file.cfg for writing");
     e.throwMe();
    }
  cfg_file
    << "% Test cfg file for use by Config::selfTest" << endl
    << "filename h:ras/ras/src/wow.dat" << endl
    << "input_name   trial_name" << endl
    << endl
    << "wavelength 3 m" << endl
    << "altitude 4000 km" << endl
    << "start_time 35 hr" << endl
    << "c1  22.3" << endl
    << "epoch_time 2005-090T20:23:29.230" << endl
    << "time_step 1 min" << endl;
  Config cfg("Config::selfTest","test_cfg_file.cfg");
  Uvar wavelength = cfg["wavelength"];
  Uvar start_time = cfg["start_time"];
  if (cfg["wavelength"] != UnitVar(3,"m")) return(false);
  if (cfg["altitude"] != UnitVar(4000,"km")) return(false);
  if (cfg["start_time"] != UnitVar(35,"hr")) return(false);
  if (cfg["c1"] != UnitVar(22.3)) return(false);
  if (cfg.str("input_name") != "trial_name") return(false);
  if (cfg.str("epoch_time") != "2005-090T20:23:29.230") return(false);
  if (cfg["time_step"] != UnitVar(1,"min")) return(false);

  return(true);
  }

//-------------
// Constructors
//-------------

//--------------------------------------------------------
// The basic construction method is from a filename.
// The file contains space delimited keyword,value pairs.
// The '%' character starts a comment line in the file.
// The '#' character starts an auto-comment line in the file.
// Blank lines are ignored.
//--------------------------------------------------------

Config::Config(const string& filename)
  : auto_comments_(false), filename_(filename),interactive_mode_on_(false),
    update_file_(true)
  {
  readFile(filename);
  // set up interactive mode if desired
  if(keywordExists("CONFIG_INTERACTIVE_MODE_ON"))
    {
    interactive_mode_on_=(bool)getInt("CONFIG_INTERACTIVE_MODE_ON");
    }
  if(keywordExists("ALLOW_CONFIG_FILE_OVERWRITE"))
    {
    update_file_=(bool)getInt("ALLOW_CONFIG_FILE_OVERWRITE");
    }
  }

//--------------------------------------------------------------------------
// If the user supplies a caller id string with the name of the
// calling program or function, then the Config object will be able to
// track which programs use which keywords.
// The basic operations reamin the same:
// The file contains space delimited keyword,value pairs.
// The '%' character starts a comment line in the file.
// The '#' character starts an auto-comment line in the file.
// Blank lines are ignored.
//--------------------------------------------------------------------------

Config::Config(const string& caller_id, const string& filename)
  : auto_comments_(false), filename_(filename), caller_id_(caller_id),
    update_file_(true)
  {
  readFile(filename);
  // set up interactive mode if desired
  if(keywordExists("CONFIG_INTERACTIVE_MODE_ON"))
    {
    interactive_mode_on_=(bool)getInt("CONFIG_INTERACTIVE_MODE_ON");
    }
  if(keywordExists("ALLOW_CONFIG_FILE_OVERWRITE"))
    {
    update_file_=(bool)getInt("ALLOW_CONFIG_FILE_OVERWRITE");
    }
  }

//----------------------------------------------------------------------
// If auto-comments are not enabled, then the destructor does nothing.
// Otherwise, the destructor re-writes the entire config file based on
// auto-comment data collected during the Config objects existence.
// User comments starting with % are left in place.
//----------------------------------------------------------------------

Config::~Config()
  {
  if (!auto_comments_ || !update_file_) return;
  cerr << "Updating " << filename_ << endl;
  writeAutoCommentFile(filename_ + ".tmp"); 
  string command("/bin/mv -f ");
  command += filename_ + ".tmp " + filename_;
  if (system(command.c_str()) != 0)
    {
    throw ErrorMessage("Config: System command: " + command + ", failed");
    }
  }

//-----------
// Operators
//-----------

//----------------------------------------------------------------------
// operator[]
//
// Calls num() to look up and return the number (possibly with units)
// associated with the key string s.
// If no match is found, an exception is thrown.
// Only const access is allowed; The config list can only be setup by
// a constructor.
// If a string member is desired, use str().
//----------------------------------------------------------------------

UnitVar Config::operator[](const string& s)  
  {
  return(num(s));
  }

//----------------
// Predicates
//----------------

//--------------------------------------------------------------------
// keywordExists(keyword)
//
// Returns true if the specified keyword is in this Config object.
// Returns false otherwise.
// The special character '*' present in the argument string
// is treated as a wildcard and will match any substring.
// Currently, the wildcard is only recognized as a trailing suffix; Any
// characters following it will cause an exception.
//--------------------------------------------------------------------

bool Config::keywordExists(const string& s) const // no exceptions
  {
  // Look for wildcard character '*' in argument
  string::size_type n = s.find('*');
  if (n != string::npos)
    {
    string s1 = s.substr(0,n);
    strings::const_iterator p;
    for (p = string_list_.begin(); p != string_list_.end(); ++p)
      {
      string s2 =p->first.substr(0,n);
      if (s1==s2)
        {
        break;
        }
      }
    if (p == string_list_.end())
      {  // no match found, check number list
      numbers::const_iterator p2;
      for (p2 = number_list_.begin(); p2 != number_list_.end(); ++p2)
        {
	string s2 =p2->first.substr(0,n);
        if (s1==s2)
          {
          break;
          }
        }
      if (p2 == number_list_.end())
        {  // no match in number list either
        return(false);
        }
      }

    if (n == s.size()-1)
      {  // nothing left after the wildcard, so this is a match
      return(true);
      }
    else
      {  // extra characters after wildcard - not supported 
      string emsg("Config::keywordExists: ");
      emsg += "Found extra characters after wildcard\n";
      emsg += "File: ";
      emsg += filename_;
      emsg += " Keyword: ";
      emsg += s;
      ErrorMessage e(emsg);
      e.throwMe();
      }
    }

  // Check regular keywords (no wildcard character present)

  strings::const_iterator p = string_list_.find(s);
  if (p != string_list_.end())
    {  // found a match
    return(true);
    }
  else
    {  // no match found, check number list
    numbers::const_iterator p2 = number_list_.find(s);
    if (p2 != number_list_.end())
      {  // found a match
      return(true);
      }
    }
  return(false);  // no match found
  }

//--------------------------------------------------------------------
// numberKeywordExists(keyword)
//
// Returns true if the specified keyword is in this Config object.
// Returns false otherwise.
//--------------------------------------------------------------------

bool Config::numberKeywordExists(const string& s) const // no exceptions
  {
  numbers::const_iterator p2 = number_list_.find(s);
  return(p2 != number_list_.end());
  }

//--------------------------------------------------------------------
// stringKeywordExists(keyword)
//
// Returns true if the specified keyword is in this Config object.
// Returns false otherwise.
//--------------------------------------------------------------------

bool Config::stringKeywordExists(const string& s) const // no exceptions
  {
  strings::const_iterator p2 = string_list_.find(s);
  return(p2 != string_list_.end());
  }

//--------------------------------------------------------------------
// isDeprecated(keyword)
//
// Returns true if the specified keyword is in the deprecated keyword
// list.  Deprecated keywords are identified in the master config file.
// Returns false otherwise.
//--------------------------------------------------------------------

bool Config::isDeprecated(const string& s) const
  {
  map<string,list<string> >::const_iterator p = deprecated_keywords_.find(s);
  return(p != deprecated_keywords_.end());
  }

//--------------------------------------------------------------------
// isDisabled(keyword)
//
// Returns true if the specified keyword is in the disabled keyword
// list.  Disabled keywords are identified in the master config file.
// Returns false otherwise.
//--------------------------------------------------------------------

bool Config::isDisabled(const string& s) const
  {
  map<string,list<string> >::const_iterator p = disabled_keywords_.find(s);
  return(p != disabled_keywords_.end());
  }

//----------------
// File I/O
//----------------

//--------------------------------------------------------------------------
// readFile(filename)
//
// Read in a configuration file containing keyword,value pairs (1 per line).
// The '%' character starts a comment line in the file.
// The '#' character starts an auto-comment line in the file.
// Blank lines are ignored.
//--------------------------------------------------------------------------

void Config::readFile(const string& filename)
  {
  ifstream filestream(filename.c_str());
  if (!filestream)
    {
    ErrorMessage e("Can't open configuration file: " + filename);
    e.throwMe();
    }
  string line;
  list<string> user_comment_lines;
  while (getline(filestream,line))
    {  // process lines
    ISTRINGSTREAM linestream(line.c_str());
    string keyword;
    UnitVar num_value;
    string str_value;
    if (linestream >> keyword)
      {  // process one keyword,value pair per line
      streampos sp = linestream.tellg();  // remember position of value

      //-------------------------------------------------------
      // Look for special comment that identifies a master
      // config file.  Its presence enables auto-comments.
      //-------------------------------------------------------

      if (keyword.substr(0,2) == "#!")
        {  // the next space delimited string is the master config file name
        if (!(linestream >> master_cfg_))
          {
          ErrorMessage e("Config: Can't read master config file name");
          e.throwMe();
          }
        auto_comments_ = true;
        readSpecialKeywords();
        }

      //-------------------------------------------------------
      // User comments are stored and put at the top of the
      // file or right before the associated keyword
      // when auto-comments are enabled.
      //-------------------------------------------------------

      if (keyword[0] == '%')
        {
        if (number_list_.size() == 0 && string_list_.size() == 0)
          {  // no keyword value pairs read yet
          user_header_comments_.push_back(line);
          }
        else
          {
          user_comment_lines.push_back(line);
          }
        continue;
        }

      //-------------------------------------------------------
      // Prior auto-comments are ingored.  They are regenerated
      // after each usage when auto-comments are enabled.
      //-------------------------------------------------------

      if (keyword[0] == '#') continue;

      //-------------------------------------------------------
      // Read in a keyword value pair. 
      //-------------------------------------------------------

      if (linestream >> num_value)
        {  // successfully read in a keyword,UnitVar pair
        // Link any user comments preceeding this keyword
        user_comments_[keyword] = user_comment_lines;
        user_comment_lines.clear();

        char c;
        linestream >> c;
        if (!linestream || isspace(c))
          {
          // There is no more contiguous non-space input so this should
          // be a number value.  Otherwise it should be a string value.
          number_list_[keyword] = num_value;
          line_list_[keyword] = line;
          line_list_copy_[keyword] = line;
          continue;  // next line
          }
        }

      // need to reset stream before next attempt
      linestream.clear();
      linestream.seekg(sp);
      if (linestream >> str_value)
        {  // successfully read in a keyword,string pair
        // Link any user comments preceeding this keyword
        user_comments_[keyword] = user_comment_lines;
        user_comment_lines.clear();

        string_list_[keyword] = str_value;
        line_list_[keyword] = line;
        line_list_copy_[keyword] = line;
        }
      else
        {
         ErrorMessage e("Config encountered keyword: " + keyword +
          " with no value in " + filename);
	 e.throwMe();
        }
      }
    }
  }

//--------------------------------------------------------------------------
// writeAutoCommentFile(filename)
//
// Read in the original configuration file, and then write out a new version
// to filename with auto-comments inserted appropriately based on the master
// config file.
// The '%' character starts a comment line in the file.
// The '#' character starts a previous auto-comment line and is ignored.
// Blank lines are ignored.
//--------------------------------------------------------------------------

void Config::writeAutoCommentFile(const string& filename)
  {
  DebugInfo dbg("Config::writeAutoCommentFile");
  
  ifstream input_stream(master_cfg_.c_str());
  if (!input_stream)
    {
    ErrorMessage e("Can't open master configuration file: " + master_cfg_);
    e.throwMe();
    }
  ofstream output_stream(filename.c_str());
  if (!output_stream)
    {
    ErrorMessage e("Can't open output configuration file: " + filename);
    e.throwMe();
    }

  // Write out master config file specification at the beginning.
  output_stream << "#! " << master_cfg_ << endl << endl;

  // Write out any user comments at the beginning of the file
  writeLines(output_stream,user_header_comments_);
  output_stream << endl;

  //----------------------------------------------------------------------
  // Read master config file and write out config file based on it.
  // Each block is read completely, and then if any of its keywords
  // were present in this Config object, then the block auto-comments
  // are written out followed by the keyword value lines from this
  // Config object (in the same order as they appear in the master config
  // file.
  //----------------------------------------------------------------------

  string line;
  string last_keyword;
  list<string> lines;
  list<string> keywords;
  modeE mode = none;
  bool block_enabled = false;
  bool value_header = false;

  while (getline(input_stream,line))
    {  // process lines
    ISTRINGSTREAM linestream(line.c_str());
    if (line == "") continue;
    string keyword;
    if (linestream >> keyword)
      {  // examine keywords and process their associated blocks
      if (keyword[0] == '%') continue;

      if (keyword == "!block")
        {
        mode = block_mode;
        if (block_enabled)
          {  // prior block had at least one keyword present, so write it out
          output_stream << "#-----------------------------------------------------------------------------" << endl;
          writeLines(output_stream,lines);
          output_stream << "#-----------------------------------------------------------------------------" << endl;
          output_stream << endl;  // one blank line after auto-comments
          expandKeywordList(keywords);  // expand out any wildcards
          writeKeywordLines(output_stream,keywords);
          output_stream << endl;  // one blank line after keyword block
          lines.clear();
          keywords.clear();
          block_enabled = false;
          }
        else
          {  // If no keywords present, don't write out this block's comments
          lines.clear();
          keywords.clear();
          }
        }
      else if (keyword == "!keyword" || keyword == "!deprecated_keyword" ||
               keyword == "!disabled_keyword")
        {
        mode = keyword_mode;
        }
      else if (keyword == "!value")
        {
        mode = value_mode;
        }
      else if (mode == block_mode)
        {
        lines.push_back("# " + line);
        }
      else if (mode == keyword_mode)
        {
        lines.push_back("# " + line);
        keywords.push_back(keyword);
        if (keywordExists(keyword)) block_enabled = true;
        last_keyword = keyword;
        mode = keyword_cont;
        value_header = true;
        }
      else if (mode == value_mode)
        {
        if (value_header)
          {
          lines.push_back("# Valid values for: " + last_keyword);
          value_header = false;
          }
        lines.push_back("# " + line);
        mode = value_cont;
        }
      else if (mode == keyword_cont)
        {  // keyword continuation lines
        lines.push_back("# " + line);
        }
      else if (mode == value_cont)
        {  // value continuation lines
        lines.push_back("# " + line);
        }
      }
    }

  //-------------------------------
  // Output last block if needed
  //-------------------------------

  if (block_enabled)
    {  // prior block had at least one used keyword, so write it out
    output_stream << "#-----------------------------------------------------------------------------" << endl;
    writeLines(output_stream,lines);
    output_stream << "#-----------------------------------------------------------------------------" << endl;
    output_stream << endl;  // one blank line after auto-comments
    writeKeywordLines(output_stream,keywords);
    output_stream << endl;  // one blank line after keyword block
    }

  //----------------------------------------------------------------------
  // Output any orphaned keywords (ie., not listed in master config file)
  //----------------------------------------------------------------------

  if (line_list_.size() > 0)
    {  // some keywords aren't in the master config file
    output_stream << "#-----------------------------------------------------------------------------" << endl;
    output_stream << "# Orphaned keywords:" << endl;
    output_stream << "# Need to add these to master config file!" << endl;
    output_stream << "#-----------------------------------------------------------------------------" << endl;
    output_stream << endl;  // one blank line after auto-comments
    if(dbg.level>0){
      cerr << "Warning: The following orphaned keywords were found in: "
      << filename_ << endl;
      cerr << "  ";
    }
    for (strings::const_iterator p = line_list_.begin();
         p != line_list_.end(); ++p)
      {
      output_stream << p->second << endl;
      if(dbg.level > 0 ) cerr << p->first << " " << endl;
      }
    output_stream << endl;  // one blank line after keyword block
    }

  // Write out keyword usage by this caller.
  output_stream << endl;
  output_stream << "# Keywords accessed by: " << caller_id_ << endl;
  writeCallerKeywords(output_stream);
  line_list_ = line_list_copy_;
  }

//----------------
// Access Methods
//----------------

//------------------------------------------------------------------
// str(s)
//
// Look up and return the string associated with the key string s.
// If no match is found, an exception is thrown.
//------------------------------------------------------------------

string Config::str(const string& s) 
  {
    static char cptr[100];
  string s2; // used by interactive mode
  map<string,list<string> >::const_iterator p1 = disabled_keywords_.find(s);
  if (p1 != disabled_keywords_.end())
    {
    ErrorMessage e("Config keyword (" + s +
      ") is disabled\n" + lineMerge(p1->second));
    e.throwMe();
    }

  strings::const_iterator p = string_list_.find(s);
  if (p == string_list_.end())
    {  // no match found, throw an exception
      // check to see if user mixed up type
      numbers::const_iterator p2 = number_list_.find(s);
      if (p2 != number_list_.end())
	{  // there is a match in the number list
	  ErrorMessage e("Config keyword " + s +
			 " matches a number, not a string");
	  e.throwMe();
	}
      else if(interactive_mode_on_){
	cerr << "Config keyword " << s <<
	  " is not in the configuration list." << endl;
	cerr << "Enter string value for this keyword:";
	cin.getline(cptr,100);
	s2=cptr;
        cerr << endl;
	addKeyword(s,s2);
        return(s2);
      }
      else{
	ErrorMessage e("Config keyword " + s +
		       " is not in the configuration list");
	e.throwMe();
      }
    }
  caller_keywords_.insert(s);
  return(p->second);
  }

//------------------------------------------------------------------
// num(s)
//
// Look up and return the number (possibly with units)
// associated with the key string s.
// If no match is found, an exception is thrown.
//------------------------------------------------------------------

UnitVar Config::num(const string& s)
  {
  // variables  used in interactive mode
  double value;
  static char cptr[100];
  Uvar u;
  string ustr;

  map<string,list<string> >::const_iterator p1 = disabled_keywords_.find(s);
  if (p1 != disabled_keywords_.end())
    {
    ErrorMessage e("Config keyword " + s +
      " is disabled" + lineMerge(p1->second));
    e.throwMe();
    }

  numbers::const_iterator p = number_list_.find(s);
  if (p == number_list_.end())
    {  // no match found, throw an exception
      // check to see if user mixed up type
      strings::const_iterator p2 = string_list_.find(s);
      if (p2 != string_list_.end())
	{  // there is a match in the string list
	  ErrorMessage e("Config keyword " + s +
			 " matches a string, not a number");
	  e.throwMe();
	}
      else if(interactive_mode_on_){
	cerr << "Config keyword " << s <<
	  " is not in the configuration list." << endl;
	cerr << "Enter numeric value for this keyword (w/o units):";
	cin.getline(cptr,100);
        value=strtod(cptr,NULL);
	cerr << endl;
        cerr << "Enter unit string:(return if no units)";
	cin.getline(cptr,100);
        ustr=cptr;
        cerr << endl;
        u=coerce_base(value,ustr);
	addKeyword(s,u,ustr);
        return(u);
      }
      
      else{
	ErrorMessage e("Config keyword " + s +
		       " is not in the configuration list");
	e.throwMe();
      }
    }
  caller_keywords_.insert(s);
  return(p->second);
  }

//------------------------------------------------------------------
// getInt(s)
//
// Look up and return the integer 
// associated with the key string s.
// If no match is found, an exception is thrown.
//------------------------------------------------------------------

int Config::getInt(const string& s) 
  {
  UnitVar x=num(s);
  int i;
  try
    {
    i=convert_to_int(x);
    }
  catch(Unit::UnitError)
    {
    OSTRINGSTREAM os;
    os << "Config::getInt unable to convert UnitVar (" << x
      << ") indexed by keyword " << s << " to integer";
    ErrorMessage e(toStr(os));
    e.throwMe();
    }
  return(i);
  }

Uvar Config::getTime( const string& keyword)
  {
  Uvar value;
  if (numberKeywordExists(keyword))
    {
    value = num(keyword);
    }
  else
    {
    string s = str(keyword);
    if (isUtc(s))
      {
      Time t(s);
      value = t;
      }
    else
      {
      ErrorMessage e("Invalid keyword, string value pair: (" +
        keyword + " " + s + ")");
      e.throwMe();
      }
    }
  return(value);
  }

//------------------------------------------------------------------------
// convertWestLongitude(str)
//
// Convert the value of the specified keyword from East longitude to
// West longitude.  This assumes that the keyword has an angle value,
// otherwise a unit error will occur.
//------------------------------------------------------------------------

Uvar Config::convertWestLongitude(const string& s)
  {
  Uvar value(2*pi,"rad");
  value-=num(s);
  return(value);
  }

//------------------------------------------------------------------
// getDouble(s)
//
// Look up and return the double 
// associated with the key string s.
// If no match is found, an exception is thrown.
//------------------------------------------------------------------

double Config::getDouble(const string& s)
  {
  UnitVar x=num(s);
  double d;
  try
    {
    d=convert_to_double(x);
    }
  catch(Unit::UnitError)
    {
    OSTRINGSTREAM os;
    os << "Config::getDouble unable to convert UnitVar (" << x
      << ") indexed by keyword " << s << " to double";
    ErrorMessage e(toStr(os));
    e.throwMe();
    }
  return(d);
  }

//------------------------------------------------------------------------
// removeKeyword(keyword)
//
// Remove the indicated keyword (and its associated value) from
// this config object.
// User comments are left in place and will slide to the next keyword
// the next time the config file is read.
//------------------------------------------------------------------------

void Config::removeKeyword(const string& keyword)
  {
  string_list_.erase(keyword);
  number_list_.erase(keyword);
  line_list_.erase(keyword);
  line_list_copy_.erase(keyword);
  }

//------------------------------------------------------------------------
// addKeyword(keyword,value)
// addKeyword(keyword,value,unit_str)
//
// Insert a new keyword value pair into this config object.
// The value can be either a string, or a numeric value (Uvar).
// If a numeric value is supplied, then a unit string must also be
// supplied (3rd argument).  The Uvar will be written out in these units.
// The new keyword/value pair will be written out to the config file
// if auto-comments are enabled.  Otherwise, they will only exist in
// this object.  The master config file determines where the new keyword
// will show up.
//------------------------------------------------------------------------

void Config::addKeyword(const string& keyword, const string& value)
  {
  string_list_[keyword] = value;
  line_list_[keyword] = keyword + " " + value;
  line_list_copy_[keyword] = keyword + " " + value;
  }

void Config::addKeyword(const string& keyword, const Uvar& value,
  const string& ustr)
  {
  number_list_[keyword] = value;
  line_list_[keyword] = keyword + " " + toStr(get_in_units(value,ustr))
    +" "+ ustr;
  line_list_copy_[keyword] = keyword + " " + toStr(get_in_units(value,ustr))
    + " "+ustr;
  }

//---------------------------------------------------------------------------
// suffixSet(keyword,subset)
//
// Find all the keywords that match the input keyword except for a
// suffix string.  Place these suffixes and their values into the subset
// numbers map supplied as an argument.
//---------------------------------------------------------------------------

void Config::suffixSet(const string& keyword, numbers& subset) const
  {
  for (numbers::const_iterator p = number_list_.begin();
       p != number_list_.end(); ++p)
    {
    if (keyword == p->first.substr(0,keyword.size()))
      {
      string suffix = p->first.substr(keyword.size(),
        p->first.size()-keyword.size());
      subset[suffix] = p->second;
      }
    }
  }

//---------------------------------------------------------------------------
// suffixSet(keyword,subset)
//
// Find all the keywords that match the input keyword except for a
// suffix string.  Place these suffixes and their values into the subset
// strings map supplied as an argument.
//---------------------------------------------------------------------------

void Config::suffixSet(const string& keyword, strings& subset) const
  {
  for (strings::const_iterator p = string_list_.begin();
       p != string_list_.end(); ++p)
    {
    if (keyword == p->first.substr(0,keyword.size()))
      {
      string suffix = p->first.substr(keyword.size(),
        p->first.size()-keyword.size());
      subset[suffix] = p->second;
      }
    }
  }

//---------------------------------------------------------------------------
// suffixSet(keyword,keyword_list,suffix_list)
//
// Find all the keywords (number and string)
// that match the input keyword except for a suffix string.
// Append these keywords to the supplied string list.
// Append the suffixes to the supplied suffix list.
//---------------------------------------------------------------------------

void Config::suffixSet(const string& keyword, list<string>& keyword_list,
  list<string>& suffix_list) const
  {
  for (strings::const_iterator p = string_list_.begin();
       p != string_list_.end(); ++p)
    {
    if (keyword == p->first.substr(0,keyword.size()))
      {
      string suffix = p->first.substr(keyword.size(),
        p->first.size()-keyword.size());
      keyword_list.push_back(p->first);
      suffix_list.push_back(suffix);
      }
    }
  for (numbers::const_iterator p = number_list_.begin();
       p != number_list_.end(); ++p)
    {
    if (keyword == p->first.substr(0,keyword.size()))
      {
      string suffix = p->first.substr(keyword.size(),
        p->first.size()-keyword.size());
      keyword_list.push_back(p->first);
      suffix_list.push_back(suffix);
      }
    }
  }

//----------------------------------------------------------------------
// s = filename
//
// Return name of the config file used to construct this Config object.
//----------------------------------------------------------------------

string Config::filename() const
  {
  return(filename_);
  }

//---------------------------
// Other support methods
//---------------------------

//----------------------------------------------------------------------
// suffixSort(suffix_list,keyword_list)
//
// Sorts a list of strings by their suffixes.
// Keywords with the same suffix are held in their original order.
// Suffixes are sorted into the order of the suffix list.
//----------------------------------------------------------------------

void Config::suffixSort(const list<string>& suffix_list,
  list<string>& keyword_list)
  {
  list<string> new_keyword_list;
  for (list<string>::const_iterator psuffix = suffix_list.begin();
       psuffix != suffix_list.end(); ++psuffix)
    {
    list<string>::iterator p = keyword_list.begin();
    list<string>::iterator tmp_p;
    while (p != keyword_list.end())
      {
      if (*psuffix == p->substr(p->size()-psuffix->size(),psuffix->size()))
        {  // this keyword has the current suffix
        // move the matching element to the end of the new keyword list
        tmp_p = p;
        p++;  // bump the loop pointer before this element gets removed
        new_keyword_list.splice(new_keyword_list.end(),keyword_list,tmp_p);
        }
      else
        {
        p++;
        }
      }
    }
  keyword_list.merge(new_keyword_list);  // add sorted keywords back to list
  }

//---------------------------
// Debugging support
//---------------------------

void Config::numbersShow(const string& name, numbers& nn)
  {
  cout << "Map: " << name << endl;
  for (Config::numbers::const_iterator p = nn.begin(); p != nn.end(); ++p)
    {
    cout << "  " << p->first << "," << p->second << endl;
    }
  }
 
void Config::showStringList(const string& name, const list<string>& slist)
  {
  cout << name << ": list size = " << slist.size() << endl;
  for (list<string>::const_iterator p = slist.begin(); p != slist.end(); ++p)
    {
    cout << "  " << (*p) << endl;
    }
  }
 
//---------------------------
// Private support methods
//---------------------------

//--------------------------------------------------------------------
// expandKeywordList(keywords)
//
// Locates any keywords with wildcard characters, and expands
// the wildcard.  The result is that the keyword with the special
// wildcard character is replaced with all of the matching keywords
// in this Config object.
// The special character '*' is the wildcard character.
// Currently, the wildcard is only recognized as a trailing suffix.
//--------------------------------------------------------------------

void Config::expandKeywordList(list<string>& keywords) const
  {
  list<string> suffix_list,t_suffix_list;
  list<string> keyword_list;  // accumulate matched keywords to add
  list<string>::iterator p = keywords.begin();
  while (p != keywords.end())
    {
    // Look for wildcard character '*' in keyword
    string::size_type n = p->find('*');
    if (n != string::npos)
      {
      string s1 = p->substr(0,n);
      t_suffix_list.clear();
      suffixSet(s1,keyword_list,t_suffix_list);  // accumulate matched keywords
      if (keyword_list.size() > 0 && suffix_list.size() == 0)
        {  // save the first non-empty suffix set for later sorting
        suffix_list = t_suffix_list;
        }
      p = keywords.erase(p);  // remove the wildcarded keyword
      }
    else
      {
      p++;
      }
    }

  // Reorder list of keywords by suffix
  suffixSort(suffix_list,keyword_list);

  // Add the matched keywords to the input keyword list
  keywords.merge(keyword_list);

  }

//--------------------------------------------------------------------------
// readSpecialKeywords()
//
// Scan the master config file for disabled and deprecated keywords and
// build internal lists to check against.
//--------------------------------------------------------------------------

void Config::readSpecialKeywords()
  {
  ifstream input_stream(master_cfg_.c_str());
  if (!input_stream)
    {
    ErrorMessage e("Can't open master configuration file: " + master_cfg_);
    e.throwMe();
    }

  //----------------------------------------------------------------------
  // Read master config file and write out config file based on it.
  // Each block is read completely, and then if any of its keywords
  // were present in this Config object, then the block auto-comments
  // are written out followed by the keyword value lines from this
  // Config object (in the same order as they appear in the master config
  // file.
  //----------------------------------------------------------------------

  string line;
  string last_keyword;
  list<string> lines;
  modeE mode = none;
  bool in_deprecated_keyword = false;
  bool in_disabled_keyword = false;

  while (getline(input_stream,line))
    {  // process lines
    ISTRINGSTREAM linestream(line.c_str());
    if (line == "") continue;
    string keyword;
    if (linestream >> keyword)
      {  // scan for deprecated and disabled keyword commands
      if (keyword[0] == '%') continue;

      if (keyword[0] == '!' && in_disabled_keyword)
        {  // finished reading disabled keyword description
        disabled_keywords_[last_keyword] = lines;
        lines.clear();
        in_disabled_keyword = false;
        mode = none;
        }
      else if (keyword[0] == '!' && in_deprecated_keyword)
        {  // finished reading deprecated keyword description
        deprecated_keywords_[last_keyword] = lines;
        lines.clear();
        in_deprecated_keyword = false;
        mode = none;
        }

      if (keyword == "!disabled_keyword")
        {
        mode = keyword_mode;
        in_disabled_keyword = true;
        }
      else if (keyword == "!deprecated_keyword")
        {
        mode = keyword_mode;
        in_deprecated_keyword = true;
        }
      else if (mode == keyword_mode)
        {
        lines.push_back(line);
        last_keyword = keyword;
        mode = keyword_cont;
        }
      else if (mode == keyword_cont)
        {  // keyword continuation lines
        lines.push_back(line);
        }
      }
    }
  }

//--------------------------------------------------------------
// writeLines(output_stream,lines)
//
// Write out a list of lines to the indicated output stream
//--------------------------------------------------------------

void Config::writeLines(ofstream& output_stream, const list<string>& lines)
  {
  for (list<string>::const_iterator p = lines.begin(); p != lines.end(); ++p)
    {
    output_stream << *p << endl;
    }
  }

//--------------------------------------------------------------
// writeKeywordLines(output_stream,keywords)
//
// Write out the stored keyword value lines from the config file
// used to initialize this Config object for the keyords supplied
// in the keywords argument.
// Also write out any user comment lines that are linked to the
// keywords (right before the keyword).
// Written keywords are removed from the line list.
//--------------------------------------------------------------

void Config::writeKeywordLines(ofstream& output_stream,
  const list<string>& keywords)
  {
  for (list<string>::const_iterator p = keywords.begin();
       p != keywords.end(); ++p)
    {
    if (keywordExists(*p))
      {
      writeLines(output_stream,user_comments_[*p]); 
      output_stream << line_list_[*p] << endl;
      line_list_.erase(*p);
      if (isDeprecated(*p) && DebugInfo::allWarnings)
        {
        cerr << "Warning: Found deprecated keyword (" << *p << ") in "
          << filename_ << endl;
        }
      if (isDisabled(*p) && DebugInfo::allWarnings)
        {
        cerr << "Warning: Found disabled keyword (" << *p << ") in "
          << filename_ << endl;
        }
      }
    }
  }

//--------------------------------------------------------------
// writeCallerKeywords(output_stream)
//
// Write out the stored caller keyword list for this config file.
//--------------------------------------------------------------

void Config::writeCallerKeywords(ostream& output_stream) const
  {
  for (set<string>::const_iterator p = caller_keywords_.begin();
       p != caller_keywords_.end(); ++p)
    {
    output_stream << "#   " << *p << endl;
    }
  }

//--------------------------------------------------------------
// str = lineMerge(string_list)
//
//
// Merge all the strings in a list of strings into one big string
// with newlines inserted between each string element.
//--------------------------------------------------------------

string Config::lineMerge(const list<string>& lines) const
  {
  string s;
  for (list<string>::const_iterator p = lines.begin();
       p != lines.end(); ++p)
    {
    s += *p + '\n';
    }
  return(s);
  }

