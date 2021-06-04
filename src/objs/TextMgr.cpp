#include <iostream>
#include "TextMgr.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::endl;

//---------------------------
// Methods for TextMgr
//---------------------------

//--------------
// Constructors
//--------------

TextMgr::TextMgr(const string& output_filename)
  : output_filename_(output_filename), special_str_("@@"),
    comment_sync_str_("Lineout mode"), mode_(source)
  {
  // Delete any previously existing file
  ofstream output_stream(output_filename_.c_str());
  if (!output_stream)
    {
    throw ErrorMessage("TextMgr::lineOut: Error opening " +
      output_filename_ + " for output");
    }
  }

//--------
// Setup
//--------

void TextMgr::setSpecialStr(const string& special_str)
  {
  special_str_ = special_str;
  }

void TextMgr::setCommentSyncStr(const string& comment_sync_str)
  {
  comment_sync_str_ = comment_sync_str;
  }

//--------
// I/O
//--------

//------------------------------------------------------------------------- 
// lineOut(input_filename,id_str)
//
// Scan the indicated input file looking for a c-style comment block 
// identified by the comment_sync_str, and followed by the indicated
// id_str.  If found, output all the lines in the remainder of the
// comment block (not including the final line with the end of comment).
// When special string expressions are found followed by (name), then
// substitute for the special_str(name) with the value associated
// with name (see addVar).
//------------------------------------------------------------------------- 

void TextMgr::lineOut(const string& input_filename, const string& id_str)
  throw(ErrorMessage)
  {
  if (id_str.size() == 0)
    {
    throw ErrorMessage("Need an id_str for TextMgr::lineOut");
    }

  //---------------------------------------------------------------------- 
  // Setup to read lines from the input file and output to the output file
  //---------------------------------------------------------------------- 

  ifstream input_stream(input_filename.c_str());
  if (!input_stream)
    {
    throw ErrorMessage("TextMgr::lineOut: Error opening " +
      input_filename + " for input");
    }

  ofstream output_stream(output_filename_.c_str(),std::ios::app);
  if (!output_stream)
    {
    throw ErrorMessage("TextMgr::lineOut: Error opening " +
      output_filename_ + " for output");
    }

  // Process one line at a time.
  string line;
  unsigned int line_number = 0;
  while (getline(input_stream,line))
    {
    line_number++;
    string::size_type len = line.size();
    if (len > 2 && line[0] == '/' && line[1] == '/')
      {  // skip c++ comments
      continue;
      }

    //-------------------------------------------------------------------
    // When in source mode, look for the indicated comment block.
    // When in text mode, output lines, handling special strings, until
    // the end of the comment block.
    //-------------------------------------------------------------------

    if (mode_ == source)
      {  // look for c-style comment start
      string::size_type n = line.find("/*");
      if (n == string::npos)
        {  // no comment found, so skip this source line
        continue;
        }
      else
        {  // check for the id string after the comment start
        string rest_of_line = line.substr(n+2,len-n-2);
        string::size_type i = rest_of_line.find(id_str);
        if (i == string::npos)
          {  // incorrect id str, so skip this comment line
          continue;
          }
        else
          {
          mode_ = text;
          }
        }
      }
    else if (mode_ == text)
      {  // xfer lines, looking for special string and end of comment block
      string::size_type n = line.find("*/");
      if (n != string::npos)
        {  // found the end of comment block
        mode_ = source;  // back to source mode
        continue;        // keep scanning for another matching comment block
        }

      string out_line;  // holds characters to be output

      while (1)
        {  // scan repeatedly for special string
        string::size_type n = line.find(special_str_);
        string::size_type remaining_len = line.size();
        if (n == string::npos)
          {  // no more special strings, xfer rest of line and finish
          out_line += line;
          break;
          }

        //--------------------------------------
        // Found a special string, so handle it.
        //--------------------------------------

        // xfer up to the special string
        out_line += line.substr(0,n);

        // toss characters already xfer'd
        line = line.substr(n,remaining_len-n);
        remaining_len = line.size();

        // look for properly formed special_str( .. )
        string::size_type n_left_paren = line.find("(");
        string::size_type n_right_paren = line.find(")");
        if (n_left_paren != special_str_.size() ||
            n_right_paren <= n_left_paren+1 ||
            n_right_paren == string::npos)
          {  // didn't find () after special_str_
          OSTRINGSTREAM ost;
          ost << "TextMgr::lineOut: Invalid special string on line "
            << line_number;
          throw ErrorMessage(ost.str());
          }

        // We have a properly formed special string.
        string::size_type arg_len = n_right_paren - n_left_paren - 1;
        string arg_str = line.substr(n_left_paren+1, arg_len);
          
        // xfer special string argument translation
        out_line += getVar(arg_str);

        // toss characters already xfer'd
        line = line.substr(n_right_paren+1, remaining_len-n_right_paren-1);
        }

      // Write the output line to the output file.
      output_stream << out_line << endl;
      }
    }
  }

//------------------
// Variable handling
//------------------

//---------------------------------------------------------------------------
// addVar(varname,varstring)
//
// Add the variable called varname to this TextMgr with the associated
// varstring to be printed whenever this varname is called for by a special
// string embedded in the text being output.
//---------------------------------------------------------------------------

void TextMgr::addVar(const string& varname, const string& varstring) throw()
  {
  variables_[varname] = varstring;
  }

//---------------------------------------------------------------------------
// varstring = getVar(varname)
//
// Return the varstring associated with varname in this TextMgr.
// Throws an exception if varname is unknown.
//---------------------------------------------------------------------------

string TextMgr::getVar(const string& varname) throw(ErrorMessage)
  {
  map<string,string>::const_iterator p = variables_.find(varname);
  if (p == variables_.end())
    {
    throw ErrorMessage("TextMgr::getVar: No variable called " + varname);
    }
  return(p->second);
  }

