//==============================================================================
// Plot.h
//
// This file contains the Plot class declaration and some supporting class
// declarations (line, sym).
// The Plot class provides a convenient interface to external plotting
// packages (eg., xmgr).
// This header comment summarizes the interface.
// For details about a specific function,
// look at the header comment in Plot.cpp.
//
// Currently, xmgr is the only external package supported.
//
// Interface summary:
//
// class Plot;
//   Defines a single plot (which may have many points, curves etc).
//
// Plot Methods:
//
//   Construction:
//
//   Plot();  Only the default constructor is supplied
//
//   Get/Set methods:
//
//   setTool(tool_name);  Specify the external plotting tool (xmgr)
//
//   Data set control:
//
//   addY(y,line,sym);      Add a single vector y
//   addY(y,ustr,line,sym); Add a single vector y and force units ustr
//   addXY(x,y,line,sym);   Add data set (x,y)
//   addXY(x,ustrx,y,line,sym);
//   addXY(x,y,ustry,line,sym);
//   addXY(x,ustrx,y,ustry,line,sym);
//   addLegend(lstr);
//
//   Graph settings:
//
//   setLimits();  set axes to autoscale
//   setLimits(xmin,xmax,ymin,ymax);  set plot limits for axes
//   setXtick(xmajor,xminor);
//   setXtick(xmajor);  xminor defaults to xmajor/2
//   setYtick(ymajor,yminor);
//   setYtick(ymajor);  yminor defaults to ymajor/2
//   setTitle(title_str);
//   setSubTitle(subtitle_str);
//   setXlabel(xlabel_str);
//   setXlabelandUnits(xlabel_str);
//   setYlabel(ylabel_str);
//   setYlabelandUnits(ylabel_str);
//   setLegendPosition(x,y);
//   setLegendTextSize(s); 
//
//
// class line,sym;
//   Defines the characteristics of lines and symbols on a plot.
//   These classes only provide constructors.
//
//  line constructors:
//    line();
//    line(type_str);
//    line(type_str,width);
//    line(type_str,color_str);
//    line(type_str,color_str,width);
//
//  sym constructors:
//    sym();
//    sym(type_str);
//    sym(type_str,size);
//    sym(type_str,color_str);
//    sym(type_str,color_str,size);
//    sym(type_str,color_str,size,fill);
//
// line types: none, solid, dash, dot, dot-dash, long-dash
// sym types:  none, dot, circle, square, diamond, up-triangle, left-triangle,
//             right-triangle, down-triangle
// colors: black, red, green, blue, yellow
//
//
//  Contour plot specific tools
// void addXYZ(const Uvec& x, const string& ustrx, 
//      const Uvec& y, const string& ustry, 
//      const Umat& z, const string& ustrz);
//    add data set of z=z(x,y) where x ans y are 1D Uvar array
// void addXYZ(const Umat& x, const string& ustrx, 
//      const Umat& y, const string& ustry, 
//      const Umat& z, const string& ustrz);
//    add data set of z=z(x,y) where x ans y are 2d non-uniform  Uvar array
//void addRectangle(const Uvec& x, const string& ustrx, 
//	    const Uvec& y, const string& ustry);
//     add rectangle in contour plot
//void addMark(const Uvec& x, const string& ustrx,
//       const Uvec& y,const string& ustry);
//     add a white mark of point (markersize 4 of matlab) 
//     at specified locations 
//    
//==============================================================================

#ifndef Plot_H
#define Plot_H

#include <string>
#include <vector>
#include <list>
#include <stdarg.h>
#include <sstream>
#include "Array.h"
#include "Units.h"
#include "Error.h"

//---------------------
// Forward declarations
//---------------------

using std::list;
using std::string;
using std::vector;

class line;
class sym;

// Supporting class legend defined here to avoid compiler error.

class legend
  {
  public:

  legend();
  void clear();

  // public variables

  list<string> strlist;
  bool fill_frame;
  bool frame;
  Uvar xpos;
  Uvar ypos;
  Uvar text_scale;
  };

//-------------------------------
// Plot making support
//-------------------------------

class Plot
  {
  public:

  typedef std::istringstream ISTRINGSTREAM;
  typedef std::ostringstream OSTRINGSTREAM;

  
  typedef list<Dvec>::const_iterator CI;
  typedef list<Dmat>::const_iterator DMI;
  typedef list<line>::const_iterator LI;
  typedef list<sym>::const_iterator SI;
  typedef list<string>::const_iterator STI;

  //--------------
  // Constructors
  //--------------

  Plot() throw();

  //-----------
  // Operators
  //-----------

  //-------------------
  // High level control
  //-------------------

  void setTool(const string& tool_name) throw(ErrorMessage);
  void clear();

  //-------------------
  // Data set control
  //-------------------

  void addY(const Uvec& y, const line& ll, const sym& ss)
    throw(ErrorMessage);
  void addY(const Uvec& y, const string& ustr,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addY(const vector<Uvar>& y, const string& ustr,
	    const line& ll, const sym& ss) throw(ErrorMessage);
  void addY(const vector<double>& y,
	    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Uvec& x, const Uvec& y,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Uvec& x, const string& ustrx,
    const Uvec& y,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Uvec& x,
    const Uvec& y, const string& ustry,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Uvec& x, const string& ustrx,
    const Uvec& y, const string& ustry,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Dvec& x, const Dvec& y, 
	     const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const vector<Uvar>& x, const string& ustrx,
	     const vector<Uvar>& y, const string& ustry,
	     const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const vector<double>& x,
	     const vector<double>& y,
	     const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Umat& x, const string& ustrx,
    const Umat& y, const string& ustry,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addXY(const Uvec& x, const string& ustrx,
    const Umat& y, const string& ustry,
    const line& ll, const sym& ss) throw(ErrorMessage);
  void addLegend(string s);

  //-----------------------
  //For 2Dcontour or image  plot
  //---------------------
  void addXYZ(const Uvec& x, const string& ustrx, 
	      const Uvec& y, const string& ustry, 
	      const Uvec& z, const string& ustrz);
  //generate color image: z becomes color value

  void addXYZ(const Uvec& x, const string& ustrx, 
	      const Uvec& y, const string& ustry, 
	      const Umat& z2, const string& ustrz);
  void addXYZ(const Dvec& x,  const Dvec& y, const Dmat& z2);
  void addXYZ(const Umat& x2, const string& ustrx, 
	      const Umat& y2, const string& ustry, 
	      const Umat& z2, const string& ustrz);

  void addRectangle(const Uvec& x, const string& ustrx, 
		    const Uvec& y, const string& ustry);
  void addMark(const Uvec& x, const string& ustrx,
	       const Uvec& y,const string& ustry);
  void addMark(const vector<Uvar>& x, const string& ustrx,
	       const vector<Uvar>& y, const string& ustry);

  //----------------
  // Graph settings
  //----------------

  void setLimits() throw();
  void setLimits(const Uvar& xmin, const Uvar& xmax,
                 const Uvar& ymin, const Uvar& ymax)
    throw(ErrorMessage);
  void setXtick(const Uvar& xmajor, const Uvar& xminor);
  void setXtick(const Uvar& xmajor);
  void setYtick(const Uvar& ymajor, const Uvar& yminor);
  void setYtick(const Uvar& ymajor);
  void setTitle(const std::string& title);
  void setSubTitle(const std::string& subtitle);
  void setXlabel(const std::string& label);
  void setXlabelandUnits(const std::string& label);
  void setYlabel(const std::string& label);
  void setYlabelandUnits(const std::string& label);
  void setLegendPosition(Uvar x, Uvar y);
  void setLegendTextSize(Uvar s);

  //--------------------
  // Display and output
  //--------------------

  void setFile(const string& filename);
  void setCmdFile(const string& cmd_filename);
  void show(const std::string& target_device) throw(ErrorMessage);
  void show2Dcontourplot(const std::string& target_device);
  void show2Dimage(const std::string& target_device);
  private:

  // Internal representation
  string cmd_filename_;
  string filename_;
  string tool_name_;
  string x_units_;
  string y_units_;
  string z_units_;
  list<Dvec> x_,y_,z_;//1D array of x , y , z
  list<Dmat> x2_,y2_,z2_;//2D array of x, y, z
 
  list<Dvec> rectanglex_,rectangley_;
  list <Dvec> markx_, marky_;
  Uvar xmin_,xmax_,ymin_,ymax_,zmin_,zmax_;
  Uvar xmajor_,xminor_,ymajor_,yminor_,zmajor_,zminor_;
  bool autox_,autoy_,autoz_;
  list<line> l_;
  list<sym> s_;
  string title_,subtitle_,xlabel_,ylabel_,zlabel_;
  legend leg_;
  
  };

//-------------------
// Supporting classes
//-------------------

class line
  {
  public:

  //--------------
  // Constructors
  //--------------

  line() throw();
  line(const string& str1);
  line(const string& str1, const string& str2);
  line(const string& str1, const float& wid);
  line(const string& str1, const string& str2, const float& wid);

  void setParam(const string& str) throw(ErrorMessage);

  //-------------
  // Public data
  //-------------
 
  string type;
  int inttype;
  string color;
  int intcolor;
  float width; 
  };

class sym
  {
  public:

  //--------------
  // Constructors
  //--------------

  sym() throw();
  sym(const string& str1);
  sym(const string& str1, const string& str2);
  sym(const string& str1, const float& wid);
  sym(const string& str1, const string& str2, const float& wid);
  sym(const string& str1, const string& str2, const float& wid,
      const int& fil);

  void setParam(const string& str) throw(ErrorMessage);

  //-------------
  // Public data
  //-------------
 
  string type;
  int inttype;
  string color;
  int intcolor;
  float size; 
  int fill;
  };


#endif




