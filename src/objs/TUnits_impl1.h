//----------------------------------------------------------------------------
// Units_impl1.h
//
// This file contains declarations for the classes supporting UnitVar.
// They need to be included before the regular UnitVar interface delcarations,
// but they are not considered part of the public interface, so they are
// separated into this file.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_units_impl1_h[] =
  "@(#) $Id: TUnits_impl1.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

//---------------------
// Forward declarations
//---------------------

class Unit;
class Units;

using std::string;
using std::list;
using std::map;
using std::istream;
using std::ostream;

typedef std::istrstream ISTRINGSTREAM;
typedef std::ostrstream OSTRINGSTREAM;

//----------------------
// Supporting structures
//----------------------

struct conversion
  {
  double conv_factor;
  int category_index;
  string category_name;
  };

struct unit_tokens
  {
  list<string> top_tokens;
  list<string> bot_tokens;
  };

//-------------------
// Class declarations
//-------------------

class Unit
  {
  public:

  typedef map<string,conversion> UNIT_TABLE;
  typedef map<string,int> CATEGORY_TABLE;
  typedef map<string,unit_tokens> UNIT_DERIVED;
  typedef UNIT_TABLE::iterator UTI;
  typedef UNIT_DERIVED::iterator UDI;

  // Friend classes that need access to Unit's representation
  friend Units;

  //----------------
  // Error Handling
  //----------------

  enum errorE {unspecified, read_error, write_error, internal_error,
    mismatch, formation, reduction_error};
  class UnitError;

  //-------------------------
  // Constructors/Destructors
  //-------------------------

  Unit();
  Unit(UTI utp);

  //------
  // Setup
  //------

  static void specifyUnitConversions(const string& filename)
    throw(UnitError);

  //-----------
  // Predicates
  //-----------

  // Category matches lie in the same unit category (eg., length) but are not
  // necessarily exact matches (eg., m and km).
  bool categoryMatch(const Unit& u) const throw();
  bool operator==(const Unit& u) const throw();
  bool operator!=(const Unit& u) const throw();
  static bool isUnit(const string& ustr) throw();

  // I/O
  friend ostream& operator<<(ostream& s, const Unit& u);

  private:

  // Internal representation
  UTI utp_;

  // Static unit conversion tables (filled by unit_setup).
  // The unit table is a mapping between a string indicating the unit name
  // and a pair consisting of an enum (categoryE) and a double
  // (multiplicative conversion factor to the base unit of the category).
  // The unit derived table is a mapping between a string (the derived unit
  // name) and a list of strings (the names of the units that make up the
  // derived unit which may themselves be derived units).
  // The following static member functions provide construct on first use
  // semantics for the unit tables.  The unit tables are allocated and
  // initialized the first time they are used.  When the program ends, the
  // unit tables are abandoned on the heap.  This approach avoids potential
  // problems with static initialization order.

  static UNIT_TABLE& unit_table() throw();
  static CATEGORY_TABLE& category_table() throw();
  static UNIT_DERIVED& unit_derived() throw();
  static Unit& suppressed_unit() throw();

  };

class Units
  {
  public:

  enum mergeE {product, quotient};

  // Constructors/Destructors
  Units();
  Units(const string& ustr) throw(Unit::UnitError);
//  Units(istream& is) throw(Unit::UnitError);
  Units(const FileMgr& file, int usize = -1)
    throw(Unit::UnitError, FileMgr::IoError);

  // Testing
  static bool selfTest(const string& filename);

  // Predicates
  bool hasUnits() const throw();

  // I/O
  void write(const FileMgr& file, int usize) const
    throw(Unit::UnitError, FileMgr::IoError);
  void read(const FileMgr& file, int usize)
    throw(Unit::UnitError, FileMgr::IoError);
  friend ostream& operator<<(ostream& s, const Units& u);
  friend istream& operator>>(istream& s, Units& u);

  // Unit conversion support
  double matchConversion(const Units& u2) const throw(Unit::UnitError);
  double factorMerge(const Units& u2, mergeE mtype) throw();

  // Other services
  void clear() throw();
  void invert() throw();
  void combine(Units u) throw();
  void toString(string& ustr) const throw();
  void toUnits(const string& ustr) throw(Unit::UnitError);
  void reduce(unsigned int n) throw(Unit::UnitError);

  static void tokenize(const string& ustr, list<string>& tokens)
    throw(Unit::UnitError);
  static void streamUnitListsIn(istream& s, list<string>& top_tokens,
    list<string>& bot_tokens) throw();

  private:

  // Internal methods
  void reduceList(list<Unit>& ulist, unsigned int n)
    throw(Unit::UnitError);
  static void
    tokensToUnits(list<string>& top_tokens, list<string>& bot_tokens,
      list<Unit>& utop, list<Unit>& ubot)
    throw(Unit::UnitError);
  static void unitSplit(list<string>& tokens,
    list<string>& toptokens, list<string>& bottokens)
    throw();
  static bool factorMatch(const Unit& u1,
    list<Unit>& unit_match_list, double& conv_factor,
    list<Unit>::iterator& u2p) throw();
  static char streamUnitTermsIn(istream& s, list<string>& uterms_list)
    throw();
  void streamUnitsIn(istream& s) throw();

  // Internal representation
  list<Unit> units_top_;
  list<Unit> units_bot_;
  };

