//----------------------------------------------------------------------------
// Units.cpp
//
// This file contains method definitions for the classes supporting UnitVar
// and supporting function definitions.  Template definitions need to be
// included in the .h file, so they are put in Units_impl2.h.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_c[] =
  "@(#) $Id: TUnits.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <iostream>
#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <strstream>
#include <fstream>
#include <string>
#include "Units.h"
#include "Io.h"
#include "Utils.h"
#include "Constants.h"

using std::list;
using std::map;
using std::string;
using std::istream;
using std::ostream;
using std::ifstream;
using std::min;

//-------------------
// Class Unit Methods
//-------------------

//-----------------
// Construction
//-----------------

Unit::Unit() { }
Unit::Unit(UTI utp) : utp_(utp) { }

//-------
// Setup
//-------

//----------------------------------------------------------------------
// specifyUnitConversions
// This static member function provides the public unit table setup
// interface.  When called the first time, it allocates and initializes
// the unit tables using data in the specified file.  Subsequent calls
// add data to the unit tables as if the files were merged together.
//----------------------------------------------------------------------

void Unit::specifyUnitConversions(const string& filename)
  throw(Unit::UnitError)
  {
  enum modeE {none,category,derived};
  modeE mode = none;

  string cur_cat_name,cur_derived_name;
  int category_index = -1;

  ifstream filestream(filename.c_str());
  if (!filestream)
    {
    throw Unit::UnitError("Error opening unit spec file: " + filename);
    }

  string line;
  while (getline(filestream,line))
    {  // process lines
    ISTRINGSTREAM linestream(line.c_str());
    string word;
    while (linestream >> word)
      {  // process words within a line
      if (word[0] == '%') break;  // comment char dumps rest of this line
      if (word == "category")
        {  // process category definition and switch to category mode
        if (!(linestream >> cur_cat_name))
          {
          throw Unit::UnitError("Missing category name");
          }
        
        mode = category;
        map<string,int>::iterator p
          = category_table().find(cur_cat_name);
        if (p != category_table().end())
          {  // category already defined, so use its index (add to it)
          category_index = p->second;
          }
        else
          {  // No match found, so create a new category type
          category_index = category_table().size();
          category_table()[cur_cat_name] = category_index;
          }
        }
      else if (word == "derived")
        {  // process derived definition and switch to derived mode
        // derived category name is read but not used (yet)
        if (!(linestream >> cur_derived_name))
          {
          throw Unit::UnitError("Missing derived category name");
          }
        mode = derived;
        }
      else if (mode == category)
        {  // process category unit definition
        double conv_factor;
        if (!(linestream >> conv_factor))
          {
          throw Unit::UnitError("Missing unit conversion factor");
          }
        conversion c;
        c.conv_factor = conv_factor;
        c.category_index = category_index;
        c.category_name = cur_cat_name;
        unit_table()[word] = c;
        }
      else if (mode == derived)
        {  // process derived category unit definition
        string rest_of_line;
        if (!(getline(linestream,rest_of_line)))
          {
          throw Unit::UnitError("Missing derived unit conversion");
          }
        unit_tokens tokens;
//        rest_of_line += " ";  // sidesteps a bug in std streams.
        ISTRINGSTREAM ist(rest_of_line.c_str());
        Units::streamUnitListsIn(ist,tokens.top_tokens,tokens.bot_tokens);
        unit_derived()[word] = tokens;
        }
      else
        {
        throw Unit::UnitError("Specified unit conversion without category");
        }
      } 
    }
  
  }

// Predicates

//-------------------------------------------------------------------------
// categoryMatch
//
// Category matches lie in the same unit category (eg., length) but are not
// necessarily exact matches (eg., m and km).
//-------------------------------------------------------------------------

bool Unit::categoryMatch(const Unit& u) const throw()
  {
  return(utp_->second.category_index == u.utp_->second.category_index);
  }

//---------------------------------
// == and != require exact matches
//---------------------------------

bool Unit::operator==(const Unit& u) const throw()
  {
  return(utp_ == u.utp_);
  }

bool Unit::operator!=(const Unit& u) const throw()
  {
  return(utp_ != u.utp_);
  }

//--------------------------------------------------------------------
// isUnit(ustr)
// Static member function returns true if ustr is a valid string
// representation of a single unit term.
//--------------------------------------------------------------------

bool Unit::isUnit(const string& ustr) throw()
  {
  if (ustr == "1") return(true);  // 1 is a special case

  Unit::UTI utp = Unit::unit_table().find(ustr);
  if (utp != Unit::unit_table().end())
    {  // found a matching unit name in the unit table
    return(true);
    }

  //--------------------------------------------------------------------
  // Didn't find the term in the unit table, so check the derived unit
  // table.
  //--------------------------------------------------------------------

  Unit::UDI udp = Unit::unit_derived().find(ustr);
  if (udp != Unit::unit_derived().end())
    {  // found a matching unit name in the derived unit table
    return(true);
    }

  //--------------------------------------------------------------------
  // Didn't find the term in the derived unit table either, so it isn't
  // a known unit term.
  //--------------------------------------------------------------------

  return(false);
  }

//--------------------
// Class Units Methods
//--------------------

// Construction/Destruction

Units::Units() { }

//-------------------------------------------------------------------------
// Construct from string
// This constructor converts a string into a Units object.
// An exception is thrown if the unit string is unrecognized, or formatted
// incorrectly.
//-------------------------------------------------------------------------

Units::Units(const string& ustr) throw(Unit::UnitError)
  {
  string ustr1 = ustr + " ";
  ISTRINGSTREAM ist(ustr1.c_str());
  streamUnitsIn(ist);
  }

//----------------------------------------------------------------------------
// Construct from file
// This constructor reads and converts a character string from the indicated
// file.  By default, the number of characters is stored first as an
// unsigned char, followed by the characters themselves.
// Like above, an exception is thrown for bad unit strings.
// See read().
//----------------------------------------------------------------------------

Units::Units(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  read(file,usize);
  }

//---------
// Testing
//---------

//----------------------------------------------------------
// selfTest
//
// Excercize the UnitVar type and verify proper operation.
//----------------------------------------------------------

bool Units::selfTest(const string& filename)
  {
  typedef UnitVar<double> Uvar;
  Unit::specifyUnitConversions(filename);

  Uvar cc("speed_light");
  string s = "2.997924580e8 m/s";  // check reading of exponential notation
  ISTRINGSTREAM ist1(s.c_str());
  ist1 >> cc;
  if (cc != Uvar(2.997924580e8,"m/s")) return(false);
  // check builtin constant
  if (cc != Uvar("speed_light")) return(false);
  
  Uvar frequency("frequency");
  frequency = Uvar(14,"GHz");
  if (!frequency.hasUnits()) return(false);
  if (fabs(frequency.getInUnits("Hz") - 14e9) > 1e-5) return(false);
  Uvar vel("vel",4,"m/s");
  Uvar mass("mass",3,"kg");
  Uvar ke("ke");
  ke = mass*pow(vel,2);
  Uvar kematch(3,"J");
  ke += Uvar(3,"J");  // check for proper matching of derived units
  Uvar a;
  if (a.coerce("") != 0) return(false);
  if (a.hasUnits()) return(false);
  Uvar b(3.5);
  if (b.coerce("km") != 3.5) return(false);
  if (!b.hasUnits()) return(false);
  // Can't mix different UnitVar types yet.
  // Need to figure out how to generalize the template definitions, or supply
  // specializations.
  //  UnitVar<float> c(7.25,"m");
  Uvar c(7.25,"km");
  if (c.coerce("cm") != 725000) return(false);
  if (c.getInUnits("mm") != 7250000) return(false);
  if (!(b+c == Uvar(10.75,"km"))) return(false);
  Uvar t("t");
  t.stringIn("240 s");
  if (t.coerce("min") != 4) return(false);
  bool errflag = false;
  try
    {
    a = c + t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  // Operator tests
  float g = 5.2;
  Uvar g1 = g*t;
  Uvar g2 = t*g;
  g1 *= g;
  g1 *= t;
  g1 = g/t;
  g2 = t/g;
  g1 /= g;
  g1 /= t;

  try
    {
    g1 = g+t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  try
    {
    g1 = g-t;
    }
  catch(Unit::UnitError& e)
    {
    if (e.error_type == Unit::mismatch) errflag = true;
    }
  if (!errflag) return(false);

  Uvar theta(30,"deg");
  double d = sin(theta);
  if (fabs(d - 0.5) > 1e-10) return(false);
  Uvar theta2 = asin(Uvar(d));
  if (abs(theta2 - theta) > 1e-10) return(false);

  theta = Uvar(60,"deg");
  d = cos(theta);
  if (fabs(d - 0.5) > 1e-10) return(false);
  theta2 = acos(Uvar(d));
  if (abs(theta2 - theta) > 1e-10) return(false);

  d = tan(theta);
  if (fabs(d - sqrt(3.0)) > 1e-10) return(false);
  theta2 = atan(Uvar(d));
  if (abs(theta2 - theta) > 1e-10) return(false);

  double aa = exp(Uvar(3.5));
  double bb = log(aa);
  if (abs(bb - 3.5) > 1e-10) return(false);

  string ustr = "4 m m/(s s)";
  ISTRINGSTREAM ist(ustr.c_str());
  ist >> a;
  if (a != Uvar(4,"m m/(s s)")) return(false);

  b = sqrt(a);
  c = Uvar(2,"m / s");
  if (b != c) return(false);

  a = Uvar(3.5,"km"); 
  b = Uvar(3600,"m");
  if (a > b) return(false);
  if (a >= b) return(false);
  if (a == b) return(false);
  if (!(a != b)) return(false);

  return(true); 
  }

// Predicates

//--------------------------------------------------------------
// hasUnits
// Returns true if there are any numerator or denominator units.
//--------------------------------------------------------------

bool Units::hasUnits() const throw()
  {
  return(units_top_.size() != 0 || units_bot_.size() != 0);
  }

// I/O

//-------------------------------------------------------------------------
// write(file,usize)
//
// Convert this Units object to a unit string, and write the string to the
// indicated file.
//
// file = FileMgr object to write to
// usize = the number of characters to write for the unit string.
//   If usize is 0, then no units are written.
//   If usize is -1, then the actual size of the unit string is written.
//      The length of the unit string is written as a unsigned char, followed
//      by the string itself (no trailing 0).  If the unit string exceeds 255
//      chars, then an exception is thrown.
//      (this is the default value if usize is not explicitely passed)
//   If usize is positive, then precisely usize characters are written
//      after the size field, and no trailing 0.  If the unit string is
//      larger than usize chars, then an exception is thrown.
//      Excess characters are padded with spaces.
//-------------------------------------------------------------------------

void Units::write(const FileMgr& file, int usize) const
  throw(Unit::UnitError, FileMgr::IoError)
  {
  if (usize == 0) return;

  string ustr;
  toString(ustr);
  if (ustr.length() > 255)
    {
    throw
      Unit::UnitError("More than 255 characters in units to write to file " +
      file.name());
    }
  unsigned char len = ustr.length();

  if (usize == -1)
    {
    file.write(len);
    file.write(ustr);
    }
  else if (len > usize)
    {
    OSTRINGSTREAM os;
    os << "More than " << usize << " characters in units to write to file "
      << file.name();
    throw Unit::UnitError(os.str());
    }
  else if (len == usize)
    {
    file.write(len);
    file.write(ustr);
    }
  else
    {
    string s(usize-len,' ');
    file.write(len);
    file.write(ustr);
    file.write(s);
    }

  }

//-------------------------------------------------------------------------
// read(file,usize)
//
// Read a unit string from the file.  usize determines how the string is
// stored just like write() above.
//-------------------------------------------------------------------------

void Units::read(const FileMgr& file, int usize)
  throw(Unit::UnitError, FileMgr::IoError)
  {
  if (usize == 0) return;

  if (usize == -1)
    {
    unsigned char len;
    file.read(len);
    string ustr(' ',len);
    file.read(ustr);
//    ustr += " ";  // sidestep stream problem repositioning at end of stream
    ISTRINGSTREAM ist(ustr.c_str());
    streamUnitsIn(ist);
    }
  else if (usize > 0)
    {
    string ustr(' ',usize);
    file.read(ustr);
//    ustr += " ";  // sidestep stream problem repositioning at end of stream
    ISTRINGSTREAM ist(ustr.c_str());
    streamUnitsIn(ist);
    }
  else
    {
    throw Unit::UnitError("Bad usize parameter passed to Units::read");
    }

  }

//------------------------
// Unit conversion support
//------------------------

//-------------------------------------------------------------------------
// matchConversion(u2)
//
// Compute the unit conversion factor which needs to multiply the argument
// units (u2) to make them match *this's units (u_).  If all unit factors
// can't be matched, then a UnitMismatch exception is thrown.
//-------------------------------------------------------------------------

double Units::matchConversion(const Units& u2) const throw(Unit::UnitError)
  {
  // For each u_ top/bot unit, find a category match in u2 and accumulate
  // the unit conversions (top/bot separately).

  double conv_factor_top = 1.0;

  list<Unit> unit_match_list = u2.units_top_;  // make copy
  list<Unit>::iterator u2p;
  for (list<Unit>::const_iterator u1p=units_top_.begin();
       u1p != units_top_.end(); ++u1p)
    {
    if (factorMatch(*u1p,unit_match_list,conv_factor_top,u2p))
      {
      unit_match_list.erase(u2p);
      }
    else
      {
      // No match found, so this is a unit mismatch because terms require
      // a category match for every unit factor.
      throw Unit::UnitError(Unit::mismatch);
      }
    }

  // If there are any unmatched units remaining in the match list, then we
  // have a mismatch because every unit in u2 should be matched in *this.
  if (unit_match_list.size() != 0) throw Unit::UnitError(Unit::mismatch);

  // Again for the denominators

  unit_match_list = u2.units_bot_;  // make copy
  double conv_factor_bot = 1.0;
  for (list<Unit>::const_iterator u1p=units_bot_.begin();
       u1p != units_bot_.end(); ++u1p)
    {
    if (factorMatch(*u1p,unit_match_list,conv_factor_bot,u2p))
      {
      unit_match_list.erase(u2p);
      }
    else
      {
      // No match found, so this is a unit mismatch because terms require
      // a category match for every unit factor.
      throw Unit::UnitError(Unit::mismatch);
      }
    }
  if (unit_match_list.size() != 0) throw Unit::UnitError(Unit::mismatch);
  
  return(conv_factor_top/conv_factor_bot);
  }

//--------------------------------------------------------------------------
// factorMerge
//
// Merge the unit terms of u2 into *this's units.
// The terms are multiplied together if mtype is product, and divided
// together if mtype is quotient.
// The return value is the accumulated conversion factor that was applied
// to make u2 match *this when category matches occurred.
//--------------------------------------------------------------------------

double Units::factorMerge(const Units& u2, mergeE mtype)
  throw()
  {
  double conv_factor_top = 1.0;
  double conv_factor_bot = 1.0;

  // Make copies which are stripped of matches below.
  list<Unit> u2_copy_top;
  list<Unit> u2_copy_bot;

  if (mtype == product)
    {
    u2_copy_top = u2.units_top_;
    u2_copy_bot = u2.units_bot_;
    }
  else
    {
    u2_copy_top = u2.units_bot_;
    u2_copy_bot = u2.units_top_;
    }

  // Process u1 (ie., *this) numerator.

  list<Unit>::iterator u1p = units_top_.begin();
  list<Unit>::iterator u2p;
  while (u1p != units_top_.end())
    {
    if (factorMatch(*u1p,u2_copy_bot,conv_factor_bot,u2p))
      {  // found a match, so cancel these units
      u1p = units_top_.erase(u1p);
      u2p = u2_copy_bot.erase(u2p);
      }
    else if (factorMatch(*u1p,u2_copy_top,conv_factor_top,u2p))
      {  // check for numerator match (no cancellation) and move on
      *u2p = *u1p;  // make the copy match for proper merging later
      ++u1p;
      }
    else
      {
      ++u1p;
      }
    }

  // Process u1 (ie., *this) denominator.

  u1p = units_bot_.begin();
  while (u1p != units_bot_.end())
    {
    if (factorMatch(*u1p,u2_copy_top,conv_factor_top,u2p))
      {  // found a match, so cancel these units
      u1p = units_bot_.erase(u1p);
      u2p = u2_copy_top.erase(u2p);
      }
    else if (factorMatch(*u1p,u2_copy_bot,conv_factor_bot,u2p))
      {  // check for denominator match (no cancellation) and move on
      *u2p = *u1p;  // make the copy match for proper merging later
      ++u1p;
      }
    else
      {
      ++u1p;
      }
    }
  
  // Merge the remaining units in u2 into the units of u1.
  units_top_.splice(units_top_.end(),u2_copy_top);
  units_bot_.splice(units_bot_.end(),u2_copy_bot);

  return(conv_factor_top/conv_factor_bot);
  }

//---------------
// Other services
//---------------

//-------------------
// clear()
// Removes any units
//-------------------

void Units::clear() throw()
  {
  units_top_.clear();
  units_bot_.clear();
  }

//------------------------------------------------------
// invert()
// Swaps the numerator units with the denominator units.
//------------------------------------------------------

void Units::invert() throw()
  {
  std::swap(units_top_,units_bot_);
  }

//-------------------------------------------------------------------
// combine()
// Merges the units (num and den) of supplied Units with *this units.
//-------------------------------------------------------------------

void Units::combine(Units ucopy) throw()
  {
  units_top_.splice(units_top_.end(),ucopy.units_top_);
  units_bot_.splice(units_bot_.end(),ucopy.units_bot_);
  }

//---------------------------------------------------------
// toString
// Convert this Units object into a string representation.
// (the opposite of toUnits)
//---------------------------------------------------------

void Units::toString(string& ustr) const throw()
  {
  string s;
  if (units_top_.size() == 0 && units_bot_.size() == 0)
    {
    return;
    }
  if (units_top_.size() == 0)
    {
    s += "1";
    }
  else if (units_top_.size() == 1)
    {
    s += (units_top_.front()).utp_->first;
    }
  else
    {
    list<Unit>::const_iterator p = units_top_.begin();
    for (unsigned int i=0; i < units_top_.size()-1; ++i)
      {
      s += p->utp_->first + " ";
      ++p;
      }
    s += p->utp_->first;
    }

  if (units_bot_.size() == 0)
    {
    ustr = s;
    return;
    }
  else if (units_bot_.size() == 1)
    {
    s += "/" +  (units_bot_.front()).utp_->first;
    }
  else
    {
    s += "/(";
    list<Unit>::const_iterator p = units_bot_.begin();
    for (unsigned int i=0; i < units_bot_.size()-1; ++i)
      {
      s += p->utp_->first + " ";
      ++p;
      }
    s += p->utp_->first + ")";
    }

  ustr = s;
  }

/****
This stream version doesn't work even though it should.  There appears to be
some errors in the strstream library implementation on cygwin/egcs.  Also,
the str() method of a stringstream doesn't return a string like it should.
Instead, it returns a char*.

void Units::toString(string& ustr) const
  {
  OSTRINGSTREAM s;
  s.flush();
  cout << "toString:";
//  cout << s.str() << " ";
  if (units_top_.size() == 0 && units_bot_.size() == 0)
    {
    return;
    }
  if (units_top_.size() == 0)
    {
    s << "1";
    }
  else if (units_top_.size() == 1)
    {
    s << units_top_.front();
    s.flush();
    cout << (units_top_.front()).utp_->first;
//    cout << s.str() << " ";
    }
  else
    {
    list<Unit>::const_iterator p = units_top_.begin();
    for (unsigned int i=0; i < units_top_.size()-1; ++i)
      {
      s << *p << " ";
      ++p;
      }
    s << *p;
    }

  if (units_bot_.size() == 0)
    {
    ustr = s.str();
    return;
    }
  else if (units_bot_.size() == 1)
    {
    s << "/" << units_bot_.front();
    s.flush();
    cout << (units_bot_.front()).utp_->first;
//    cout << s.str() << " ";
    }
  else
    {
    s << "/(";
    list<Unit>::const_iterator p = units_bot_.begin();
    for (unsigned int i=0; i < units_bot_.size()-1; ++i)
      {
      s << *p << " ";
      ++p;
      }
    s << *p << ")";
    }

//  cout << s.str() << " ";
  if (s.str()) ustr = s.str();
  cout << ustr << " " << s.str() << ":toString" << endl;
  cout << ustr.size() << endl;
  exit(1);
  }
*****/

//---------------------------------------------------------
// toUnits()
// Convert the supplied string into *this Units.
// (the opposite of toString)
// Throws an exception if the string is not valid units.
//---------------------------------------------------------

void Units::toUnits(const string& ustr) throw(Unit::UnitError)
  {
  string ustr1 = ustr + " ";
  ISTRINGSTREAM ist(ustr1.c_str());
  streamUnitsIn(ist);
  }

//-------------------------------------------------
// reduce(n)
//
// Reduce the units in *this by the n'th power.
// Thus, "m m/(s s)" is reduced to "m/s" if n = 2.
// Empty units mean no action.
// n = 0 or 1 means no action.
//-------------------------------------------------

void Units::reduce(unsigned int n) throw(Unit::UnitError)
  {
  if (!hasUnits()) return;
  reduceList(units_top_, n);
  reduceList(units_bot_, n);
  }

//----------------------------------------------------------------------
// ***** This method is currently not being used.
// tokenize()
// Parse's strings into subwords suitable for unit analysis.
// Declared static in Units.h because it does not need to access a Units
// object.
// **** Needs to handle () just like stream output.
//----------------------------------------------------------------------

void Units::tokenize(const string& ustr, list<string>& tokens)
  throw(Unit::UnitError)
  {
  ISTRINGSTREAM ist(ustr.c_str());
  string token;
  char ch = 0;
  while (ist >> ch)
    {
    if (ch == '/')
      {
      tokens.push_back("/");
      }
    else if (ch == '1')
      {
      tokens.push_back("1");
      }
    else if (isalpha(ch))
      {
      ist.putback(ch);
      ist >> token;
      tokens.push_back(token);
      }
    else
      {
      throw Unit::UnitError("Bad unit string: " + ustr);
      }
    }
  }

//----------------------------
// Protected internal methods
//----------------------------

//--------------------------------------------------
// reduceList(ulist,n)
//
// Reduce a unit list by the order n (see reduce())
//--------------------------------------------------

void Units::reduceList(list<Unit>& ulist, unsigned int n) throw(Unit::UnitError)
  {
  //---------------------------------------
  // Reduce successively by dumping matches
  //---------------------------------------

  for (unsigned int i=2; i <= n; ++i)
    {
    for (list<Unit>::iterator p = ulist.begin();
         p != ulist.end(); ++p)
      {  // look for a match to each unit term
      Unit u = *p;  // unit term to match
      ++p;          // advance to start of search list
      list<Unit>::iterator pmatch = find(p,ulist.end(),u);
      --p;          // backup to original position
      if (pmatch == ulist.end())
        {  // no match - this is an error
        throw Unit::UnitError("UnitVar reduction error; No match for " +
          p->utp_->first, Unit::reduction_error);
        }
      else
        {  // found match, dump the match
        ulist.erase(pmatch);
        }
      }
    }
      
  }

//-----------------------------------
// tokensToUnits
//
// Convert token lists into unit list.
//-----------------------------------

// Need to rework some of this so that stream input is used for all processing
// so that I don't need to put the () processing into tokenize.
// tokensToUnits should accept numerator and denominator strings separately,
// and it should do the recursive call via toUnits.  The derived unit
// mapping should probably be strings instead of token lists.

void Units::tokensToUnits(list<string>& top_tokens, list<string>& bot_tokens,
  list<Unit>& utop, list<Unit>& ubot) throw(Unit::UnitError)
  {
  // Process numerator tokens, expanding out derived units, and adding
  // the Unit terms to utop and ubot.

  for (list<string>::iterator p=top_tokens.begin();
       p != top_tokens.end(); ++p)
    {
    if (*p == "1") continue;  // skip over numerator '1'
    Unit::UTI utp = Unit::unit_table().find(*p);
    if (utp != Unit::unit_table().end())
      {  // found a matching unit name in the unit table
      utop.push_back(Unit(utp));
      }
    else
      { // No fundamental unit match, so look for a derived unit match.
      Unit::UDI udp = Unit::unit_derived().find(*p);
      if (udp != Unit::unit_derived().end())
        {  // found matching derived unit: recursively process it
        tokensToUnits(udp->second.top_tokens,udp->second.bot_tokens,utop,ubot);
        }
      else
        {  // No matching derived unit either, so this is an unknown unit.
        throw Unit::UnitError("Bad unit string: " + *p);
        }
      }
    }
      
  // Process denominator tokens, expanding out derived units, and adding
  // the Unit terms to utop and ubot in inverted sense (ie., numerator of
  // denominator units go in the denominator, and denominator of denominator
  // units go in the numerator!)

  for (list<string>::iterator p=bot_tokens.begin();
       p != bot_tokens.end(); ++p)
    {
    Unit::UTI utp = Unit::unit_table().find(*p);
    if (utp != Unit::unit_table().end())
      {  // found a matching unit name in the unit table
      ubot.push_back(Unit(utp));
      }
    else
      { // No fundamental unit match, so look for a derived unit match.
      Unit::UDI udp = Unit::unit_derived().find(*p);
      if (udp != Unit::unit_derived().end())
        {  // found matching derived unit: recursively process it
        // Note reversal of utop,ubot order because these are denominator
        // units.
        tokensToUnits(udp->second.top_tokens,udp->second.bot_tokens,ubot,utop);
        }
      else
        {  // No matching derived unit either, so this is an unknown unit.
        throw Unit::UnitError("Bad unit string: " + *p);
        }
      }
    }

  }

//-------------------------------------------------------
// *****This method not used anymore.
// unitSplit
// Split token list into numerator and denominator parts.
//-------------------------------------------------------

void Units::unitSplit(list<string>& tokens,
  list<string>& toptokens, list<string>& bottokens)
  throw()
  {
  list<string>::iterator pdiv = find(tokens.begin(),tokens.end(),"/");
  if (pdiv == tokens.end())
    {  // no division symbol, so everything is in the numerator
    toptokens = tokens;
    return;
    }
  for (list<string>::iterator p = tokens.begin(); p != pdiv; ++p)
    {
    toptokens.push_back(*p);
    }
  ++pdiv;  // move to start of denominator list
  for (list<string>::iterator p = pdiv; p != tokens.end(); ++p)
    {
    if (*p == "(" || *p == ")") continue;  // skip any enclosing ()
    bottokens.push_back(*p);
    }
  }
  
//-------------------------------------------------------------------
// factorMatch
// Look for a match for u1 in the list unit_match_list.  When found,
// accumulate the conversion, and return an iterator to the matching
// element.  Return true if a match found, false otherwise.
//-------------------------------------------------------------------

bool Units::factorMatch(const Unit& u1,
  list<Unit>& unit_match_list, double& conv_factor,
  list<Unit>::iterator& u2p)
  throw()
  {
  u2p=find_if(unit_match_list.begin(),unit_match_list.end(),CategoryEqual(u1));

  // If no match found, return without doing anything.
  if (u2p == unit_match_list.end()) return(false);

  // If the units within the category don't match, then a conversion
  // needs to be applied to the 2nd term (ie., *u2p) to match the
  // 1st term (u1).
  if (*u2p != u1)
    {
    conv_factor *= u2p->utp_->second.conv_factor /
                       u1.utp_->second.conv_factor;
    }
  return(true);
  }

//---------------------------------------------------------------------------
// streamUnitTermsIn(s,uterms_list)
//
// Stream input of Unit terms into a list of strings.
// Each string in the list is a potential individual unit term.
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
// This function only handles simple lists of terms (eg., m m s s) but not
// numerator/denominator combinations which are handled by the << operator
// (with help from this protected function).
// The '/' and ')' characters are recognized as special, however,
// so that spaces are not required between the last numerator term and the '/',
// or between the last denominator term and the ')'.
// The '1' character is also recognized as special (for cases where there
// are no numerator units), and is a legitimate unit term.
// The last character read is returned because the std strstream appears to
// have a bug that prevents the last character on a stream from being putback.
// Thus, it needs to be checked separately.
//---------------------------------------------------------------------------

char Units::streamUnitTermsIn(istream& s, list<string>& uterms_list)
  throw()
  {
  string word;
  string ustr;  // accumulates the valid unit string
  string unit_term;

  //----------------------------------------------------------------------
  // Pull in whitespace delimited words and check each to see if it might
  // contribute to a valid unit string.
  //----------------------------------------------------------------------

  while (1)
    {
    streampos sp = s.tellg();  // remember current stream position
    if (!(s >> word)) break;   // no more input

    string::size_type n1 = word.find('/');
    string::size_type n2 = word.find(')');
    if (n1 == string::npos && n2 == string::npos && Unit::isUnit(word))
      {  // simple unit term
      uterms_list.push_back(word);
      continue;
      }
    else if (n1 == string::npos && n2 == string::npos)
      {  // invalid unit term; this ends further unit input
      s.seekg(sp);
      break;  // done reading stream
      }
    else
      {  // Found a '/' or ')' so we have a possible unit term; xxx/.. or xxx)
      string::size_type n = min(n1,n2); // select the first occuring
      unit_term = word.substr(0,n);  // extract xxx before the first occuring
      if (!Unit::isUnit(unit_term))
        {  // invalid unit term; this ends further unit input
        if (unit_term.size() == 0)
          {  // the '/' or ')' is by itself, so reposition on it
          s.seekg(n - word.size(), std::ios::cur);
          }
        else
          {  // reposition at start of invalid unit term
          s.seekg(sp);
          }
        }
      else
        {  // unit term before '/' or ')' is legitimate
        uterms_list.push_back(unit_term);  // final unit term
        // reposition on the '/' or ')'
        s.seekg(n - word.size(), std::ios::cur);
        // We need to clear the flags here just in case the stream ended right
        // after the last unit term which would leave the eof flag set, and
        // cause subsequent processing (of the denominator) to fail in an
        // obscure way!
        s.clear(std::ios::goodbit);
        }
      break;  // done reading stream
      }
    }

  if (uterms_list.empty())
    {  // no valid unit terms and not eof, so set failed
    s.clear(std::ios::failbit);
    }

  // Return last character read from stream
  return(word[word.size()-1]);

  }
      
//---------------------------------------------------------------------------
// streamUnitListsIn(s,top_tokens,bot_tokens)
//
// Stream input of Unit lists.
// The numerator list is returned in top_tokens, and the denominator list
// is returned in bot_tokens.  The basic structure has to be n1/d1 or
// n1 n2 n3/(d1 d2 d3).
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
//---------------------------------------------------------------------------

void Units::streamUnitListsIn(istream& s, list<string>& top_tokens,
  list<string>& bot_tokens) throw()
  {
  top_tokens.clear();
  bot_tokens.clear();

  // First read in the numerator terms
  char last_c = streamUnitTermsIn(s,top_tokens);    

  if (s)
    {  // There is more input to process after the numerator terms
    // Check for a '/' which introduces denominator terms.
    // More than one denominator terms will be enclosed by ().
    char c = 0;
    streampos sp = s.tellg();  // remember position of '/' (if present)
    s.get(c);
//    if (!(s >> c)) return;  // Couldn't actually read anything!
    if (!s) return;
    if (c == '/')
      {  // we have a denominator
      c = ' ';
      while (isspace(c))
        {  // skip white space after '/'
        if (!(s >> c))
          {  // the '/' is the last character, so no units left
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        }

      if (c == '(')
        {  // we have a multi-term denominator inside ()
        last_c = streamUnitTermsIn(s,bot_tokens);    
        if (!s)
          {  // no more input or invalid input
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        s >> c;
        if ((s && c != ')') || (!s && last_c != ')'))
          {
          // no more input or there is more input, but it isn't the needed ')'
          s.seekg(sp);                 // reposition at '/'
          return;
          }
        }
      else
        {  // a simple denominator term
        s.seekg(sp+1);                 // reposition after '/'
        streamUnitTermsIn(s,bot_tokens);    
        }
      }
    }
  
  }

//---------------------------------------------------------------------------
// streamUnitsIn(s)
//
// Stream input of Units (into *this).
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
//---------------------------------------------------------------------------

void Units::streamUnitsIn(istream& s) throw()
  {
  list<string> top_tokens;
  list<string> bot_tokens;

  // Read in the numerator and denominator unit lists
  streamUnitListsIn(s,top_tokens,bot_tokens);

  if (top_tokens.empty())
    {  // no valid unit string and not eof, so set failed
    s.clear(std::ios::failbit);
    }
  else
    {  // we have some legitimate unit string, so convert it
    // Here is where *this is actually set.
    tokensToUnits(top_tokens,bot_tokens,units_top_,units_bot_);
    }
  }

//------------------------------------------------------------------------
// Initialize unit tables
// These static member functions provide construct on first use semantics.
//------------------------------------------------------------------------

Unit::UNIT_TABLE& Unit::unit_table() throw()
  {
  static Unit::UNIT_TABLE* p = new Unit::UNIT_TABLE;
  return(*p);
  }

Unit::CATEGORY_TABLE& Unit::category_table() throw()
  {
  static Unit::CATEGORY_TABLE* p = new Unit::CATEGORY_TABLE;
  return(*p);
  }

Unit::UNIT_DERIVED& Unit::unit_derived() throw()
  {
  static Unit::UNIT_DERIVED* p = new Unit::UNIT_DERIVED;
  return(*p);
  }

Unit& Unit::suppressed_unit() throw()
  {
  static Unit* p = new Unit;
  return(*p);
  }

//----------------------
// Class UnitVar Methods
//----------------------

// These are all in Units_impl2.h because g++ doesn't support the
// export keyword.

//----------------------------------------//
// I/O Supporting functions and operators //
//----------------------------------------//

//----------------------
// Stream output of Unit
//----------------------

ostream& operator<<(ostream& s, const Unit& u)
  {
  return s << u.utp_->first;
  }

//--------------------------------------------------------------
// Stream output of Units (string representation is written out)
//--------------------------------------------------------------

ostream& operator<<(ostream& s, const Units& u)
  {
  string ustr;
  u.toString(ustr);
  s << ustr;
  return s;
  }

//---------------------------------------------------------------------------
// Stream input of Units.
// Input stream is read until the first unrecognized or invalid element
// is encountered. Check the state of the stream to see if input was
// terminated by the end of stream or something else.
// If no valid unit string elements are found, but there are some characters
// before the end of the stream, then the state is set to FAIL.
//---------------------------------------------------------------------------

istream& operator>>(istream& s, Units& u)
  {
  u.streamUnitsIn(s);
  return(s);
  }

