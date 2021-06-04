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
  "@(#) $Id: Units_impl1.h,v 11.5 2011/09/16 00:03:30 richw Exp $";

//---------------------
// Forward declarations
//---------------------

class Unit;
class Units;
struct conversion;
struct unit_tokens;
class UnitTable;
class UnitDerivedTable;

using std::string;
using std::list;
using std::map;
using std::istream;
using std::ostream;

//-------------------
// Class declarations
//-------------------

class Unit
  {
  public:

  typedef UnitTable UNIT_TABLE;
  typedef UnitDerivedTable UNIT_DERIVED;
//  typedef map<string,conversion> UNIT_TABLE;
//  typedef map<string,int> CATEGORY_TABLE;
//  typedef map<string,unit_tokens> UNIT_DERIVED;
  typedef map<string,conversion>::iterator UTI;
  typedef map<string,unit_tokens>::iterator UDI;

  // Friend classes that need access to Unit's representation
  friend class Units;

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

  //-----------
  // Predicates
  //-----------

  // Category matches lie in the same unit category (eg., length) but are not
  // necessarily exact matches (eg., m and km).
  bool categoryMatch(const Unit& u) const ; // no exceptions
  bool operator==(const Unit& u) const ; // no exceptions
  bool operator!=(const Unit& u) const ; // no exceptions
  static bool isUnit(const string& ustr) ; // no exceptions

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

  static UNIT_TABLE& unit_table(); // no exceptions
//  static CATEGORY_TABLE& category_table(); // no exceptions
  static UNIT_DERIVED& unit_derived(); // no exceptions
  static Unit& suppressed_unit(); // no exceptions

  };

class Units
  {
  public:

  enum mergeE {product, quotient};
  enum term_handlingE {allow_any_terms, known_terms_only};
  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  // Constructors/Destructors
  Units();
  Units(const string& ustr) ;
//  Units(istream& is) ;
  Units(const FileMgr& file, int usize = -1);


  // Testing
  static bool selfTest(const string& filename);

  // Predicates
  bool hasUnits() const ; // no exceptions
  bool hasNonSuppressedUnits() const;
  bool operator!=(const Units& u) const ; // no exceptions

  // I/O
  void write(const FileMgr& file, int usize) const;
  void read(const FileMgr& file, int usize);
  friend ostream& operator<<(ostream& s, const Units& u);
  friend istream& operator>>(istream& s, Units& u);
  void show() const;

  // Unit conversion support
  double matchConversion(const Units& u2) const ;
  double factorMerge(const Units& u2, mergeE mtype) ; // no exceptions
  Units baseUnits() const;

  // Other services
  void clear(); // no exceptions
  void invert(); // no exceptions
  void combine(Units u); // no exceptions
  void toString(string& ustr) const; // no exceptions
  string toString() const;
  void toUnits(const string& ustr);
  void reduce(unsigned int n);
  void removeSuppressedUnit();

  static void tokenize(const string& ustr, list<string>& tokens);
  static void streamUnitListsIn(istream& s, term_handlingE term_handling,
    list<string>& top_tokens, list<string>& bot_tokens); // no exceptions

  private:

  // Internal methods
  void reduceList(list<Unit>& ulist, unsigned int n);
  static void
    tokensToUnits(list<string>& top_tokens, list<string>& bot_tokens,
      list<Unit>& utop, list<Unit>& ubot);
  static void unitSplit(list<string>& tokens,
    list<string>& toptokens, list<string>& bottokens); // no exceptions
  static bool factorMatch(const Unit& u1,
    list<Unit>& unit_match_list, double& conv_factor,
    list<Unit>::iterator& u2p); // no exceptions
  // no exceptions
  static char streamUnitTermsIn(istream& s, term_handlingE term_handling,
    list<string>& uterms_list);
  void streamUnitsIn(istream& s); // no exceptions
  static void removeSuppressedUnit(list<Unit>& ulist);

  // Internal representation
  list<Unit> units_top_;
  list<Unit> units_bot_;
  };

//----------------------
// Supporting structures
//----------------------

struct conversion
  {
  double conv_factor;
  int category_index;
  string category_name;
  Unit base_unit;
  };

struct unit_tokens
  {
  list<string> top_tokens;
  list<string> bot_tokens;
  };

//-----------------------------------------------------------------------
// Class UnitTable
//
// Class UnitTable derives from a mapping between unit strings and the
// associated conversion.  It adds a special constructor which loads
// static data into the map.  It also contains the static data to be
// loaded into the unit tables.  The actual data lists are in Units.cpp
//-----------------------------------------------------------------------

class UnitTable : public map<string,conversion>
  {
  public:

  typedef map<string,conversion>::iterator UTI;

  // Constructor
  UnitTable();

  private:

  const static unsigned int N_units_ = 43;
  const static double conv_factor_[];
  const static char* unit_name_[];
  const static int category_index_[];
  const static char* category_name_[];
  };

//-----------------------------------------------------------------------
// Class UnitDerivedTable
//
// Class UnitDerivedTable derives from a mapping between unit strings and the
// associated string representing the unit in fundamental units.
// It adds a special constructor which loads static data into the map.
// It also contains the static data to be loaded.
// The actual data lists are in Units.cpp
//-----------------------------------------------------------------------

class UnitDerivedTable : public map<string,unit_tokens>
  {
  public:

  typedef map<string,unit_tokens>::iterator UDI;
  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  // Constructor
  UnitDerivedTable();

  private:

  const static unsigned int N_derived_ = 53;
  const static char* derived_unit_name_[];
  const static char* derived_unit_string_[];
  };
