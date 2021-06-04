//----------------------------------------------------------------------------
// Plot.cpp
//
// This file contains method definitions for the Plot classes.
// These classes use an external plotting package to make plots.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_plot_c[] =
  "@(#) $Id: Plot.cpp,v 11.7 2017/04/07 22:14:16 richw Exp $";

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "Plot.h"
#include "Io.h"
#include "Array.h"
#include "Error.h"

using std::cout;
using std::endl;
using std::flush;
using std::setprecision;

//--------------------
// Class Plot methods
//--------------------

//--------------
// Construction
//--------------

Plot::Plot() throw()
   : cmd_filename_("tmp_file.xmgr"), filename_("xmgr_plot.eps"),
     tool_name_("xmgrace"), x_units_(""), y_units_(""),z_units_(""),
     xmin_(0), xmax_(1), ymin_(0), ymax_(1),zmin_(0),zmax_(1),
     xmajor_(1), xminor_(0.5),
     ymajor_(1), yminor_(0.5), 
     zmajor_(1), zminor_(0.5),
     autox_(true), autoy_(true), autoz_(true) { }

//-------------------
// High level control
//-------------------

void Plot::setTool(const string& tool_name) throw(ErrorMessage)
  {
  if (tool_name == "xmgr")
    {  // only accept valid tool names
    tool_name_ = tool_name;
    cmd_filename_="tmp_file.xmgr"; 
    }
  else if(tool_name=="xmgrace")
    {
      tool_name_=tool_name;
      cmd_filename_="tmp_file.xmgr"; 
    }
  else if (tool_name =="matlab")
    { // let's test mlab
    tool_name_ = tool_name;
    cmd_filename_="tmp_file.m";
    if(filename_=="xmgr_plot.eps")
      {
	filename_ = "matlab_plot.fig";
	//this is to avoid possible conflicts.  For example, user might
	//set the matlab file name first and then choose the tool later.
	//in that case, the file name should keep the one set by setFile()
      }
    }
  else
    {
    ErrorMessage e("Invalid tool name: " + tool_name + " in Plot");
    e.throwMe();
    }
  }

//-------------------------------------------------------------------------
// clear()
//
// Clears all data structures except the tool name so that this Plot object
// can be used again for a new plot.
//-------------------------------------------------------------------------

void Plot::clear()
  {
  cmd_filename_ = "";
  filename_ = "";
  autox_ = true;
  autoy_ = true;
  x_units_ = "";
  y_units_ = "";
  z_units_="";
  xmin_ = 0;
  xmax_ = 1;
  ymin_ = 0;
  ymax_ = 1;
  zmin_ = 0;
  zmax_ = 1;
  xmajor_ = 1;
  xminor_ = 0.5;
  ymajor_ = 1;
  yminor_ = 0.5;
  zmajor_ = 1;
  zminor_ = 0.5;
  x_.clear();
  y_.clear();
  z_.clear();
  rectanglex_.clear();
  rectangley_.clear();
  l_.clear();
  s_.clear();
  title_ = "";
  subtitle_ = "";
  xlabel_ = "";
  ylabel_ = "";
  leg_.clear();
  }

//------------------
// Data set control
//------------------

//-----------------------------------
// addY
//
// Add a single vector to this plot.
//-----------------------------------

void Plot::addY(const Uvec& y, const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addY received zero length vector: " + y.name());
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec x("x",y.size());
  for (unsigned int i=0; i < x.size(); ++i) x(i) = i;
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(y(0));

  // Add a new set at the end of the lists (x is indices)
  x_.push_back(x);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);
  }

void Plot::addY(const Uvec& y, const string& ustr,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addY received zero length vector: " + y.name());
    }

  // Enforce specified unit conformity and xfer to Dvec storage.
  Dvec x("x",y.size());
  for (unsigned int i=0; i < x.size(); ++i) x(i) = i;
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(ustr);

  // Add a new set at the end of the lists (x is indices)
  x_.push_back(x);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);

  y_units_ = ustr;
  }

void Plot::addY(const vector<Uvar>& y, const string& ustr, 
		 const line& ll, const sym& ss) throw(ErrorMessage)
  {
    if(y.size()==0)
      {
	throw ErrorMessage("Plot::addY receive zero length Y vector container");
      }
    // Enforce specified unit conformity and xfer to Dvec storage.
    Dvec x("x",y.size());
    for (unsigned int i=0; i < x.size(); ++i) x(i) = i;
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y[i].getInUnits(ustr);

    // Add a new set at the end of the lists (x is indices)
    x_.push_back(x);
    y_.push_back(newy);
    l_.push_back(ll);
    s_.push_back(ss);

    y_units_ = ustr;
  }


void Plot::addY(const vector<double>& y, 
		 const line& ll, const sym& ss) throw(ErrorMessage)
  {
    if(y.size()==0)
      {
	throw ErrorMessage("Plot::addY receive zero length Y vector container");
      }
    // Enforce specified unit conformity and xfer to Dvec storage.
    Dvec x("x",y.size());
    for (unsigned int i=0; i < x.size(); ++i) x(i) = i;
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y[i];

    // Add a new set at the end of the lists (x is indices)
    x_.push_back(x);
    y_.push_back(newy);
    l_.push_back(ll);
    s_.push_back(ss);
  }

void Plot::addXY(const Uvec& x, const Uvec& y,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length yvector: " + y.name());
    }
  else if (x.size() != y.size())
    {
    OSTRINGSTREAM os;
    os << "Plot::addXY received mismatched vectors: " <<
      x.name() << "(" << x.size() << ") and " <<
      y.name() << "(" << y.size() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec newx("newx",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(x(0));
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(y(0));

  x_.push_back(newx);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);
  }

void Plot::addXY(const Uvec& x, const string& ustrx,
		 const Uvec& y,
		 const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length yvector: " + y.name());
    }
  else if (x.size() != y.size())
    {
    OSTRINGSTREAM os;
    os << "Plot::addXY received mismatched vectors: " <<
      x.name() << "(" << x.size() << ") and " <<
      y.name() << "(" << y.size() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec newx("newx",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(ustrx);
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(y(0));

  x_.push_back(newx);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);

  x_units_ = ustrx;
  }

void Plot::addXY(const Uvec& x, const string& ustrx,
  const Uvec& y, const string& ustry,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length yvector: " + y.name());
    }
  else if (x.size() != y.size())
    {
    OSTRINGSTREAM os;
    os << "Plot::addXY received mismatched vectors: " <<
      x.name() << "(" << x.size() << ") and " <<
      y.name() << "(" << y.size() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec newx("newx",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(ustrx);
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(ustry);

  x_.push_back(newx);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);

  x_units_ = ustrx;
  y_units_ = ustry;
  }

void Plot::addXY(const Uvec& x,
  const Uvec& y, const string& ustry,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXY received zero length yvector: " + y.name());
    }
  else if (x.size() != y.size())
    {
    OSTRINGSTREAM os;
    os << "Plot::addXY received mismatched vectors: " <<
      x.name() << "(" << x.size() << ") and " <<
      y.name() << "(" << y.size() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec newx("newx",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(x(0));
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(ustry);

  x_.push_back(newx);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);

  y_units_ = ustry;
  }



void Plot::addXY(const Dvec& x, const Dvec& y, 
		 const line& ll, const sym& ss) throw(ErrorMessage)
  {
    //check size first
    if (x.size() == 0)
      {
	throw ErrorMessage("Plot::addXY received zero length xvector: " + x.name());
      }
    if (y.size() == 0)
      {
	throw ErrorMessage("Plot::addXY received zero length yvector: " + y.name());
      }
    else if (x.size() != y.size())
      {
	OSTRINGSTREAM os;
	os << "Plot::addXY received mismatched vectors: " <<
	  x.name() << "(" << x.size() << ") and " <<
	  y.name() << "(" << y.size() << ")";
	throw ErrorMessage(toStr(os));
    }
    //store data into private variables
    x_.push_back(x);
    y_.push_back(y);
    l_.push_back(ll);
    s_.push_back(ss);
    x_units_="";
    y_units_="";
  }

void Plot::addXY(const vector<Uvar>& x, const string& ustrx,
	     const vector<Uvar>& y, const string& ustry,
	     const line& ll, const sym& ss) throw(ErrorMessage)
  {
    //check size first
    if (x.size() == 0)
      {
	throw ErrorMessage("Plot::addXY received zero length xvector: " );
      }
    if (y.size() == 0)
      {
	throw ErrorMessage("Plot::addXY received zero length yvector: " );
      }
    else if (x.size() != y.size())
      {
	OSTRINGSTREAM os;
	os << "Plot::addXY received mismatched vectors: " <<
	  "x size"  << "(" << x.size() << ") and " <<
	  "y size" << "(" << y.size() << ")";
	throw ErrorMessage(toStr(os));
    }

    Dvec newx("x",x.size());
    for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x[i].getInUnits(ustrx);
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y[i].getInUnits(ustry);
    //store data into private variables
    x_.push_back(newx);
    y_.push_back(newy);
    l_.push_back(ll);
    s_.push_back(ss);
    x_units_=ustrx;
    y_units_=ustry;
    
  }

 void Plot::addXY(const vector<double>& x,
	     const vector<double>& y,
	     const line& ll, const sym& ss) throw(ErrorMessage)
{
  //check size first
  if (x.size() == 0)
    {
      throw ErrorMessage("Plot::addXY received zero length xvector: " );
    }
  if (y.size() == 0)
    {
      throw ErrorMessage("Plot::addXY received zero length yvector: " );
    }
  else if (x.size() != y.size())
    {
      OSTRINGSTREAM os;
      os << "Plot::addXY received mismatched vectors: " <<
	"x size"  << "(" << x.size() << ") and " <<
	"y size" << "(" << y.size() << ")";
      throw ErrorMessage(toStr(os));
    }
  
  Dvec newx("x",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x[i];
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y[i];
  //store data into private variables
  x_.push_back(newx);
  y_.push_back(newy);
  l_.push_back(ll);
  s_.push_back(ss);
  
}
//---------------------------------------------------------------------------
// addXY(x,ustrx,y,ustry,line,sym)
//
// When the input data arrays are matrices (with matching dimensions), then
// treat each row of input Array2D's x and y as separate sets on the graph.
// This allows for easy drawing of many line segments or other repeated
// symbol sets.
//---------------------------------------------------------------------------

void Plot::addXY(const Umat& x, const string& ustrx,
  const Umat& y, const string& ustry,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.rows() != y.rows() || x.cols() != y.cols())
    {
    throw ErrorMessage("Plot::addXY received mismatched matrices: " +
      x.name() + " and " + y.name());
    }
  for(unsigned int i=0; i < x.rows(); ++i)
    {  // loop over sets (one per row)
    addXY(x.getRow(i),ustrx,y.getRow(i),ustry,ll,ss);
    }

  x_units_ = ustrx;
  y_units_ = ustry;
  }

//---------------------------------------------------------------------------
// addXY(x,ustrx,y,ustry,line,sym)
//
// When the first input data array is a vector (1D) and the second
// input data array is a matrix, then the input vector is used for each
// row or column of the matrix.  If the matrix has only one dimension matching
// the vector, then that determines if rows or columns are used.  If the
// matrix is square and matches the vector x, then rows are used by default.
// (ie., each row becomes one set).
// This allows for easy drawing of many line segments or other repeated
// symbol sets.
//---------------------------------------------------------------------------

void Plot::addXY(const Uvec& x, const string& ustrx,
  const Umat& y, const string& ustry,
  const line& ll, const sym& ss)
  throw(ErrorMessage)
  {
  if (x.size() == y.cols())
    {
    for(unsigned int i=0; i < y.rows(); ++i)
      {  // loop over sets (one per row)
      addXY(x,ustrx,y.getRow(i),ustry,ll,ss);
      }
    }
  else if (x.size() == y.rows())
    {
    for(unsigned int i=0; i < y.cols(); ++i)
      {  // loop over sets (one per column)
      addXY(x,ustrx,y.getCol(i),ustry,ll,ss);
      }
    }
  else
    {
    throw ErrorMessage("Plot::addXY received mismatched matrices: " +
      x.name() + " and " + y.name());
    }

  x_units_ = ustrx;
  y_units_ = ustry;
  }

//--------------------
// Add a legend
//--------------------
void Plot::addLegend(const string s)
{
  leg_.strlist.push_back(s);
}

//--------------------
// Graph settings
//--------------------

void Plot::setLimits() throw()
  {
  autox_ = true;
  autoy_ = true;
  //cout << "setLimits()" << endl;
  //cout << "  autox_ = " << autox_ << endl;
  //cout << "  autoy_ = " << autoy_ << endl;
  }

void Plot::setLimits(const Uvar& xmin, const Uvar& xmax,
                     const Uvar& ymin, const Uvar& ymax)
  throw(ErrorMessage)
  {
  xmin_ = xmin;
  xmax_ = xmax;
  ymin_ = ymin;
  ymax_ = ymax;
  autox_ = false;
  autoy_ = false;
  }

void Plot::setXtick(const Uvar& xmajor, const Uvar& xminor)
  {
  xmajor_ = xmajor;
  xminor_ = xminor;
  }

void Plot::setXtick(const Uvar& xmajor)
  {
  setXtick(xmajor,xmajor/2.0);
  }

void Plot::setYtick(const Uvar& ymajor, const Uvar& yminor)
  {
  ymajor_ = ymajor;
  yminor_ = yminor;
  }

void Plot::setYtick(const Uvar& ymajor)
  {
  setYtick(ymajor,ymajor/2.0);
  }

void Plot::setTitle(const std::string& title)
  {
  title_ = title;
  }

void Plot::setSubTitle(const std::string& subtitle)
  {
  subtitle_ = subtitle;
  }

void Plot::setXlabel(const std::string& xlabel)
  {
  xlabel_ = xlabel;
  }

void Plot::setXlabelandUnits(const std::string& xlabel)
  {
  xlabel_ = xlabel + " (" + x_units_ + ")";
  }

void Plot::setYlabel(const std::string& ylabel)
  {
  ylabel_ = ylabel;
  }

void Plot::setYlabelandUnits(const std::string& ylabel)
  {
  ylabel_ = ylabel + " (" + y_units_ + ")";
  }

void Plot::setLegendPosition(Uvar x, Uvar y)
  {
  leg_.xpos=x;
  leg_.ypos=y;
  }

void Plot::setLegendTextSize(Uvar s)
  {
  leg_.text_scale=s;
  }

//------------------------------------------------------------------
//addXYZ(x,ustrx,y,ustry,z,ustrz)
// x: Uvec with n elements and unit of ustrx
// y: Uvec with n elements and unit of ustry
// z: Uvec with n elements and unit of ustrz
// This method is used to generate color coded image with z as color value
//-------------------------------------------------------------------
void Plot::addXYZ(const Uvec& x, const string& ustrx, 
	      const Uvec& y, const string& ustry, 
	      const Uvec& z, const string& ustrz)
  {
    if (x.size() == 0)
      {
	throw ErrorMessage("Plot::addXYZ received zero length xvector: " + x.name());
      }
    if (y.size() == 0)
      {
	throw ErrorMessage("Plot::addXYZ received zero length yvector: " + y.name());
    }

    if (z.size() == 0)
      {
	throw ErrorMessage("Plot::addXYZ received zero length zvector: " + z.name());
    }
    
    if (x.size() != y.size() || y.size() != z.size() || z.size() != x.size())
      {
	throw ErrorMessage("Plot::addXYZ received different sizes of x, y, z");
      }
    // Enforce unit conformity and xfer to Dvec storage.
    Dvec newx("newx",x.size());
    for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(ustrx);
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(ustry);
    Dvec newz("newz",z.size());
    for (unsigned int i=0; i < newz.size(); ++i) newz(i) = z(i).getInUnits(ustrz);
    
    x_.push_back(newx);
    y_.push_back(newy);
    z_.push_back(newz);
    
    x_units_ = ustrx;
    y_units_ = ustry;
    z_units_ = ustrz;
  }

//------------------------------------------------------------------
//addXYZ(x,ustrx,y,ustry,z2,ustrz)
// x: Uvec with n elements and unit of ustrx
// y: Uvec with n elements and unit of ustry
// z2: Umat with n x n elements and unit of ustrz
// As of now, x-dimension(number of elements) should be same as y-dimension
// and z should be a square matrix of (x-dimension) x (y-dimension)
// 
// This method works only with matlab tool. If xmgr were chosen, it would throw
// an errormessage.
//
//----------------------------------------------------------------
void Plot::addXYZ(const Uvec& x, const string& ustrx, 
		  const Uvec& y, const string& ustry, 
		  const Umat& z2, const string& ustrz)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXYZ received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXYZ received zero length yvector: " + y.name());
    }
  
  if(x.size() !=z2.rows())
    {
    OSTRINGSTREAM os;
    os<< "Plot::addXYZ received mismatched vectors: "<<
      x.name()<< "(" << x.size() << ") and "<<
      z2.name()<< "(" << z2.rows() << ")";
    throw ErrorMessage(toStr(os));
    }
  else if (y.size() !=z2.cols())
    {
    OSTRINGSTREAM os;
    os<< "Plot::addXYZ received mismatched vectors: "<<
      y.name()<< "(" << y.size() << ") and "<<
      z2.name()<< "(" << z2.cols() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
  Dvec newx("newx",x.size());
  for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(ustrx);
  Dvec newy("newy",y.size());
  for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(ustry);
  Dmat newz("newz",x.size(),y.size());
 
  for (unsigned int i = 0; i < newx.size();++i)
  for (unsigned int j = 0; j < newy.size();++j)
    {
    newz(i,j)= z2(i,j).getInUnits(ustrz);
    }
  x_.push_back(newx);
  y_.push_back(newy);
  z2_.push_back(newz);
  
  x_units_ = ustrx;
  y_units_ = ustry;
  z_units_ = ustrz;
  }

//-----------------------------------------------
//2D contour plot for double input
//-----------------------------------------------

void Plot::addXYZ(const Dvec& x, const Dvec& y, const Dmat& z2)
  {
  if (x.size() == 0)
    {
    throw ErrorMessage("Plot::addXYZ received zero length xvector: " + x.name());
    }
  if (y.size() == 0)
    {
    throw ErrorMessage("Plot::addXYZ received zero length yvector: " + y.name());
    }
  
  if(x.size() !=z2.rows())
    {
    OSTRINGSTREAM os;
    os<< "Plot::addXYZ received mismatched vectors: "<<
      x.name()<< "(" << x.size() << ") and "<<
      z2.name()<< "(" << z2.rows() << ")";
    throw ErrorMessage(toStr(os));
    }
  else if (y.size() !=z2.cols())
    {
    OSTRINGSTREAM os;
    os<< "Plot::addXYZ received mismatched vectors: "<<
      y.name()<< "(" << y.size() << ") and "<<
      z2.name()<< "(" << z2.cols() << ")";
    throw ErrorMessage(toStr(os));
    }

  // Enforce unit conformity and xfer to Dvec storage.
 
  x_.push_back(x);
  y_.push_back(y);
  z2_.push_back(z2);
  
  x_units_ = " ";
  y_units_ = " ";
  z_units_ = " ";
  }


//------------------------------------------------------------------
//addXYZ(x2,ustrx,y2,ustry,z2,ustrz)
// x2: Umat with n x n elements and unit of ustrx
// y2: Umat with n x n elements and unit of ustry
// z2: Umat with n x n elements and unit of ustrz
// all three x,y,z have the same matrix size: n x n
// 
// This method works only with matlab tool. 
//If xmgr were chosen, it would throw
// an errormessage.
//
//----------------------------------------------------------------
void Plot::addXYZ(const Umat& x2, const string& ustrx, 
		  const Umat& y2, const string& ustry, 
		  const Umat& z2, const string& ustrz)
  {
  if (x2.rows() != y2.rows() 
      || y2.rows() != z2.rows() 
      || z2.rows() != x2.rows())
    {
    throw ErrorMessage("Plot::addXYZ mismatched x,y,z  matrix sizes");
    }
  if (x2.cols() != y2.cols() 
      || y2.cols() != z2.cols() 
      || z2.cols() != x2.cols())
    {
    throw ErrorMessage("Plot::addXYZ : mismatched x,y,z matrix sizes");
    }
  if (x2.cols() != x2.rows() 
      || y2.cols() != y2.rows() 
      || z2.cols()!=z2.rows())
    {
    throw ErrorMessage("Plot::addXYZ : x,y,z are not square matrix");
    }
  
 
  // Enforce unit conformity and xfer to Dvec storage.
  unsigned int Ndat = x2.cols();
  Dmat newx("newx",Ndat,Ndat);
  Dmat newy("newy",Ndat,Ndat);
  Dmat newz("newz",Ndat,Ndat);
  for (unsigned int i=0; i < Ndat; ++i)
  for (unsigned int j = 0 ; j < Ndat;++j)
    {
      newx(i,j) = x2(i,j).getInUnits(ustrx);
      newy(i,j) = y2(i,j).getInUnits(ustry);
      newz(i,j) = z2(i,j).getInUnits(ustrz);
    }
  x2_.push_back(newx);
  y2_.push_back(newy);
  z2_.push_back(newz);
  x_units_ = ustrx;
  y_units_ = ustry;
  z_units_ = ustrz;
  }



//------------------------------------------------------------------
//addRectangle(x,ustrx,y,ustry)
// x: Uvec with 4 elements and unit of ustrx
// y: Uvec with 4 elements and unit of ustry
// this method will add an rectangle in 2D contourplot
// 
// This method works only with matlab tool. If xmgr were chosen, it would throw
// an errormessage.
//
//----------------------------------------------------------------
void Plot::addRectangle(const Uvec& x, const string& ustrx, 
			const Uvec& y, const string& ustry)
  {
    if (x.size() !=4)
      {
	ErrorMessage("Plot::addRectangle: only four x coordinate values are required to define a rectangle").throwMe();
      }
    if(y.size() !=4)
      {
	ErrorMessage("Plot::addRectangle: only four y coordinate values are required to define a rectangle").throwMe();
      }
    if(ustrx !=x_units_)
      {
	ErrorMessage("Plot::addRectangle: x coordinate unit mismatch to that used for data");
      }
    if(ustry !=y_units_)
      {
	ErrorMessage("Plot::addRectangle: y coordinate unit mismatch to that used for data");
      }
    //Enforce unit conformity and xfer to Dvec storage.
    Dvec newx("newx",x.size());
    for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(x_units_);
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(y_units_);

    rectanglex_.push_back(newx);
    rectangley_.push_back(newy);

  }


//------------------------------------------------------------------
//addMark(const Uvec& x, const Uvec& y): add mark(point) at (x,y)
//                                         in contour plot
// x: x position of mark
// y: y position of mark
// this method will add marks in 2D contour plot
// 
// This method works only with matlab tool. If xmgr were chosen, it would throw
// an errormessage.
//
//----------------------------------------------------------------
void Plot::addMark(const Uvec& x,const string& ustrx,
		   const Uvec& y,const string& ustry)
  {
    if (x.size() != y.size())
      {
	ErrorMessage("Plot::addMark: size of x is different from size of y").throwMe();
      }
   
    if(ustrx !=x_units_)
      {
	ErrorMessage("Plot::addMark: x coordinate unit mismatch to that used for data");
      }
    if(ustry !=y_units_)
      {
	ErrorMessage("Plot::addMark: y coordinate unit mismatch to that used for data");
      }
    //Enforce unit conformity and xfer to Dvec storage.
    Dvec newx("newx",x.size());
    for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x(i).getInUnits(x_units_);
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y(i).getInUnits(y_units_);

    markx_.push_back(newx);
    marky_.push_back(newy);

  }


void Plot::addMark(const vector<Uvar>& x,const string& ustrx,
		   const vector<Uvar>& y,const string& ustry)
  {
    if (x.size() != y.size())
      {
	ErrorMessage("Plot::addMark: size of x is different from size of y").throwMe();
      }
   
    if(ustrx !=x_units_)
      {
	ErrorMessage("Plot::addMark: x coordinate unit mismatch to that used for data");
      }
    if(ustry !=y_units_)
      {
	ErrorMessage("Plot::addMark: y coordinate unit mismatch to that used for data");
      }

    //Enforce unit conformity and xfer to Dvec storage.
    Dvec newx("newx",x.size());
    for (unsigned int i=0; i < newx.size(); ++i) newx(i) = x[i].getInUnits(x_units_);
    Dvec newy("newy",y.size());
    for (unsigned int i=0; i < newy.size(); ++i) newy(i) = y[i].getInUnits(y_units_);

    markx_.push_back(newx);
    marky_.push_back(newy);
  }





//--------------------
// Display and output
//--------------------

//---------------------------------------------------------
// setFile(filename)
//
// Set the name of the file used for postscript output
//---------------------------------------------------------

void Plot::setFile(const string& filename)
  {
  filename_ = filename;
  }

//---------------------------------------------------------
// setCmdFile(cmd_filename)
//
// Set the name of the file used to hold the grahic commands
// that the external tool uses to make the plot.
//---------------------------------------------------------

void Plot::setCmdFile(const string& cmd_filename)
  {
  cmd_filename_ = cmd_filename;
  }

//----------------------------------------------------------------
// show
//
// Fill the command file with the graphic commands that make
// the plot specified in this object, and then use a system call
// to invoke the external plotting tool to make the plot.
//----------------------------------------------------------------

void Plot::show(const string& target_device) throw(ErrorMessage)
  {

/*
  cout << "Plot:" << endl;
  cout << "  x_units_ = " << x_units_ << endl;
  cout << "  y_units_ = " << y_units_ << endl;
  cout << "  autox_ = " << autox_ << endl;
  cout << "  autoy_ = " << autoy_ << endl;
  cout << "  xmin_.hasUnits = " << xmin_.hasUnits() << endl;
  cout << "  ymin_.hasUnits = " << ymin_.hasUnits() << endl;
  cout << "  xmin_ = " << xmin_ << endl;
  cout << "  xmax_ = " << xmax_ << endl;
  cout << "  ymin_ = " << ymin_ << endl;
  cout << "  ymax_ = " << ymax_ << endl;
  cout << "  xmajor_ = " << xmajor_ << endl;
  cout << "  xminor_ = " << xminor_ << endl;
  cout << "  ymajor_ = " << ymajor_ << endl;
  cout << "  yminor_ = " << yminor_ << endl;
*/
  // setup units
  double xmin,xmax,xmajor,xminor;
  double ymin,ymax,ymajor,yminor;

  if (autox_ && autoy_ && tool_name_ == "xmgrace")
    {  // grace needs autoscaling and autoticking supplied externally!
    // Autoscale
    xmin = x_.front()(0);
    xmax = xmin;
    ymin = y_.front()(0);
    ymax = ymin;
    Plot::CI py = y_.begin();
    for (Plot::CI px = x_.begin();
	 px != x_.end() && py != y_.end(); ++px, ++py)
      {
      for (unsigned int i=0; i < (*px).size(); ++i)
	{
        if ((*px)(i) < xmin) xmin = (*px)(i);
        if ((*px)(i) > xmax) xmax = (*px)(i);
        if ((*py)(i) < ymin) ymin = (*py)(i);
        if ((*py)(i) > ymax) ymax = (*py)(i);
	}
      }
    }

  if (x_units_ != "" && !autox_)
    {
    if (!xmin_.hasUnits())
      {
      ErrorMessage("Plot.cpp::show: Need to setLimits with matching units").throwMe();
      }
    xmin = xmin_.getInUnits(x_units_);
    xmax = xmax_.getInUnits(x_units_);
    xmajor = xmajor_.getInUnits(x_units_);
    xminor = xminor_.getInUnits(x_units_);
    }
  else if (autox_ && tool_name_ == "xmgrace")
    {
    // Auto-tick X using already computed xmin,xmax
    //cout << "Auto-tick X" << endl;
    double tick = (xmax - xmin)/6;
    if (tick > 0)
      {
      int pp = (int)(log(tick)/log(10.0)); 
      if (pp < 0) pp--;
      xmajor = tick/pow(10.0,pp);
      if (xmajor > 0.0)
        {
        xmajor += 0.5;
        }
      else
        {
        xmajor -= 0.5;
        }
      xmajor = (int)xmajor;
      xmajor *= pow(10.0,pp);
      xminor = xmajor/2.0;
      }
    else
      {
      xmajor = 1;
      xminor = 0.5;
      }
    }
  else
    {
    xmin = xmin_.getValue();
    xmax = xmax_.getValue();
    xmajor = xmajor_.getValue();
    xminor = xminor_.getValue();
    }

  if (y_units_ != "" && !autoy_)
    {
    if (!ymin_.hasUnits())
      {
      ErrorMessage("Plot.cpp::show: Need to setLimits with matching units").throwMe();
      }
    ymin = ymin_.getInUnits(y_units_);
    ymax = ymax_.getInUnits(y_units_);
    ymajor = ymajor_.getInUnits(y_units_);
    yminor = yminor_.getInUnits(y_units_);
    }
  else if (autoy_ && tool_name_ == "xmgrace")
    {
    // Auto-tick Y
    //cout << "Auto-tick Y" << endl;
    double tick = (ymax - ymin)/6;
    if (tick > 0)
      {
      int pp = (int)(log(tick)/log(10.0)); 
      if (pp <= 0) pp--;
      ymajor = tick/pow(10.0,pp);
      if (ymajor > 0.0)
        {
        ymajor += 0.5;
        }
      else
        {
        ymajor -= 0.5;
        }
      ymajor = (int)ymajor;
      ymajor *= pow(10.0,pp);
      yminor = ymajor/2.0;
      }
    else
      {
      ymajor = 1;
      yminor = 0.5;
      }
    }
  else
    {
    ymin = ymin_.getValue();
    ymax = ymax_.getValue();
    ymajor = ymajor_.getValue();
    yminor = yminor_.getValue();
    }

  std::ofstream g_commands(cmd_filename_.c_str());
  string tmp_shell_filename ="tmp_file.sh";//matlab only
  std::ofstream s_commands(tmp_shell_filename.c_str());//matlab only    
  
  if (tool_name_ =="xmgr" || tool_name_=="xmgrace")
    {
    if (tool_name_ == "xmgrace")
      {
      g_commands << "@version 50114" << endl;
      }
    g_commands << "@with g0" << endl << "@g0 on" << endl;
    if (!autox_ || tool_name_ == "xmgrace")
      {  // auto-scaling works for xmgr, but not for grace
      g_commands << "@world xmin " <<
	xmin << endl;
      g_commands << "@world xmax " <<
	xmax << endl;
      g_commands << "@xaxis tick major " <<
	xmajor << endl;
      if (tool_name_ == "xmgr")
        {
        g_commands << "@xaxis tick minor " <<
	  xminor << endl;
        }
      else if (xminor != 0.0)
        {
        g_commands << "@xaxis tick minor ticks " <<
	  xmajor/xminor - 1 << endl;
        }
      else
        {
        g_commands << "@xaxis tick minor ticks 1" << endl;
        }
      }
    if (!autoy_ || tool_name_ == "xmgrace")
      {  // auto-scaling works for xmgr, but not for grace
      g_commands << "@world ymin " <<
	ymin << endl;
      g_commands << "@world ymax " <<
	ymax << endl;
      g_commands << "@yaxis tick major " <<
	ymajor << endl;
      if (tool_name_ == "xmgr")
        {
        g_commands << "@yaxis tick minor " <<
	  yminor << endl;
        }
      else if (yminor != 0.0)
        {
        g_commands << "@yaxis tick minor ticks " <<
	  ymajor/yminor - 1 << endl;
        }
      else
        {
        g_commands << "@yaxis tick minor ticks 1" << endl;
        }
      }

    g_commands << "@xaxis tick on" << endl;
    g_commands << "@xaxis tick major grid on" << endl;
    g_commands << "@xaxis tick minor grid on" << endl;
    g_commands << "@xaxis tick major linestyle 2" << endl;
    g_commands << "@xaxis tick minor linestyle 2" << endl;
    
    g_commands << "@yaxis tick on" << endl;
    g_commands << "@yaxis tick major grid on" << endl;
    g_commands << "@yaxis tick minor grid on" << endl;
    g_commands << "@yaxis tick major linestyle 2" << endl;
    g_commands << "@yaxis tick minor linestyle 2" << endl;
    
    if (title_.size() > 0)
      {
	g_commands << "@title \"" << title_ << "\"" << endl;
      }
    if (subtitle_.size() > 0)
      {
	g_commands << "@subtitle \"" << subtitle_ << "\"" << endl;
      }
    if (xlabel_.size() > 0)
      {
	g_commands << "@xaxis label \"" << xlabel_ << "\"" << endl;
      }
    if (ylabel_.size() > 0)
      {
	g_commands << "@yaxis label \"" << ylabel_ << "\"" << endl;
      }
    
    // Write out line and symbol settings
    
    int i = 0;
    for (Plot::LI pl = l_.begin(); pl != l_.end(); ++pl, ++i)
      {
	g_commands << "@s" << i << " linestyle " << pl->inttype << endl;
	g_commands << "@s" << i << " color " << pl->intcolor << endl;
	g_commands << "@s" << i << " linewidth " << pl->width << endl;
      }
    
    i = 0;
    for (Plot::SI ps = s_.begin(); ps != s_.end(); ++ps, ++i)
      {
	g_commands << "@s" << i << " symbol " << ps->inttype << endl;
	g_commands << "@s" << i << " symbol color " << ps->intcolor << endl;
	g_commands << "@s" << i << " symbol size " << ps->size << endl;
	g_commands << "@s" << i << " symbol fill " << ps->fill << endl;
      }
    
    // Write out legend settings
    // These need to come after the line and symbol settings, otherwise
    // Grace will complain about unallocated sets.

    if(!leg_.strlist.empty())
      {
      g_commands << "@ legend on" << endl;
      g_commands << "@ legend loctype view" << endl;
      g_commands << "@ legend layout 0" << endl;
      g_commands << "@ legend vgap 2" << endl;
      g_commands << "@ legend hgap 1" << endl;
      g_commands << "@ legend length " << leg_.strlist.size() << endl;
      string fill_frame_str="off";
      string frame_str="off";
      if(leg_.frame)
        {
	frame_str="on";
	if(leg_.fill_frame) fill_frame_str="on";
        }
      g_commands << "@ legend box "<< frame_str << endl;
      g_commands << "@ legend box fill  "<< fill_frame_str << endl;
      g_commands << "@ legend box fill with color "<< endl;
      g_commands << "@ legend box fill color 0"<< endl;
      g_commands << "@ legend box fill pattern 1"<< endl;
      g_commands << "@ legend box color 1"<< endl;
      g_commands << "@ legend box linewidth 1"<< endl;
      g_commands << "@ legend box linestyle 1"<< endl;
      g_commands << "@ legend x1 "<< leg_.xpos << endl;
      g_commands << "@ legend y1 "<< leg_.ypos << endl;
      g_commands << "@ legend font 4 "<< endl;
      g_commands << "@ legend char size "<< leg_.text_scale << endl;
      g_commands << "@ legend linewidth 1"<< endl;
      g_commands << "@ legend linestyle 1"<< endl;
      }
    i = 0;  
    for (Plot::STI pstr = leg_.strlist.begin(); pstr != leg_.strlist.end(); 
	 ++pstr, ++i)
      {
      g_commands << "@ legend string " << i << "" << '"' << *pstr <<'"'
        << endl;
      }

    // Write each data set to file with & in between each set.
    Plot::CI py = y_.begin();
    for (Plot::CI px = x_.begin();
	 px != x_.end() && py != y_.end(); ++px, ++py)
      {
      for (unsigned int i=0; i < (*px).size(); ++i)
	{
	g_commands << setprecision(16) << (*px)(i) << " " << (*py)(i) << endl;
	}
      g_commands << "&" << endl;
      }
    
    if (autox_ && tool_name_ != "xmgrace")
      {
	g_commands << "@autoscale xaxes" << endl;
      }
    if (autoy_ && tool_name_ != "xmgrace")
      {
	g_commands << "@autoscale yaxes" << endl;
      }
    
    g_commands << flush;
    }
  else if (tool_name_ =="matlab")
    {  
    //create binary data file
    string tmp_binary= cmd_filename_+".dat";
    FileMgr fp(tmp_binary,"wb");

    unsigned int Nset = 0;
    Plot::CI py = y_.begin();
    for (Plot::CI px = x_.begin();
	 px != x_.end() && py != y_.end(); ++px, ++py)
      {
      Nset += 1;
      }
    //cout<<"data size"<<Ndata<<endl;
    fp.write(Nset);

    py = y_.begin();
    for (Plot::CI px = x_.begin(); px != x_.end() && py !=y_.end() ; ++px,++py)
      {
	unsigned int Ndata = (*px).size();
	fp.write(Ndata);
	for (unsigned int i=0; i < (*px).size(); ++i)
	  {
	    fp.write((*px)(i));	 
	  }
	for (unsigned int i=0; i < (*py).size(); ++i)
	  {
	    fp.write((*py)(i));	 
	  }
      }
    fp.close();
    
    g_commands<<"fid=fopen('"<<tmp_binary<<"');"<<endl;//file open
    g_commands<<"Nset = fread(fid,1,'uint32') ; "<<endl;//read number of data
    g_commands<<" for i = 1:Nset "<<endl;
    g_commands<<"Ndata = fread(fid,1,'uint32');"<<endl;
    g_commands<<"x = fread(fid,Ndata,'double'); "<<endl;
    g_commands<<"y=fread(fid,Ndata,'double'); "<<endl;
    g_commands<<"plot(x,y,'-r'); "<<endl;
    g_commands<<"hold on"<<endl;
    g_commands<<"end "<<endl;
  
    if(title_.size() > 0)
      {
      g_commands<<"title('"<<title_<<"');"<<endl;  
      }
    if (xlabel_.size() > 0)
      {
      g_commands<<"xlabel('"<<xlabel_<<"');"<<endl;   
      }
       
    if(ylabel_.size() >0 )
      {
      g_commands<<"ylabel('"<<ylabel_<<"');"<<endl;   
      }

    if(target_device=="x")
      {//no matlab filename is specified
	g_commands<<"uiwait;"<<endl;//display and wait
      }
    else if (target_device=="matlab_fig")
      {
	g_commands<<"saveas(gcf,'"<<filename_<<"');"<<endl;//display and wait
      }
    else if (target_device=="matlab_eps")
      {
	g_commands<<"saveas(gcf,'"<<filename_<<"','psc2');"<<endl;
      }
    else
      {
	ErrorMessage("Plot.cpp::show: Invalid target device(x or matlab-fig)").throwMe();
      }
    
    g_commands<<"quit;"<<endl;
    g_commands << flush;

    string cmd_filename;
    unsigned int Nstring = cmd_filename_.size();
    cmd_filename.resize(Nstring-2);
    for (unsigned int i = 0; i<Nstring-2;i++)
      {
      cmd_filename += toStr(cmd_filename_[i]);
      }       
    s_commands<<"#! /bin/sh "<<endl;
    s_commands<<" matlab -nodesktop -nosplash -r "<<cmd_filename <<endl;
    s_commands<<flush;
    //---------------------------------
    //system command matlab -nodesktop -nosplash -r tm_file.m
    // did not work.  As a result, we decide to sidestep this problem
    // by making a shell command file named tmp.sh.  
    // content of tmp.sh
    //  #! /bin/sh
    //  matlab -nodesktop -nosplash -r cmd_filename
    //--------------------------------
    }
  else
    {
      throw ErrorMessage("Invalid tool name "+tool_name_+"in Plot");
    }
  
  string buf;

  if (tool_name_ == "xmgr")
    {
    if (target_device == "x")
      {
      buf = "xmgr -maxsets 5000 " + cmd_filename_;
      }
    else if (target_device == "eps-landscape")
      {
      buf = "xmgr -maxsets 5000 -noask -device 2 -hardcopy -eps -printfile " +
        filename_ + " " + cmd_filename_;
      }
    else if (target_device == "eps-portrait")
      {
      buf = "xmgr -maxsets 5000 -noask -device 1 -hardcopy -eps -printfile " +
        filename_ + " " + cmd_filename_;
      }
    else
      {
      throw ErrorMessage("Invalid target device: " + target_device);
      }
    }
  else if (tool_name_ == "xmgrace")
    {
      if (target_device == "x")
	{
	buf = "xmgrace -free " + cmd_filename_;
	}
      else if (target_device == "eps-landscape")
	{
	buf = "grace  -noask -hdevice EPS -hardcopy  -printfile " +
	  filename_ + " " + cmd_filename_;
	}
      else
	{
	  throw ErrorMessage("Invalid target device: " + target_device);
	}
    }
  else if(tool_name_=="matlab")
    {
    buf ="sh " + tmp_shell_filename;
    }
  else
    {
    throw ErrorMessage("Invalid tool name: " + tool_name_ + " in Plot");
    }

  if (system(buf.c_str()) != 0)
    {
    throw ErrorMessage("System command: " + buf + ", failed");
    }

  }

//---------------------------
//Show 2D contour plot
//--------------------------
void Plot::show2Dcontourplot(const std::string& target_device)
  {
  std::ofstream g_commands(cmd_filename_.c_str());
  string tmp_shell_filename ="tmp_file.sh";//matlab only
  std::ofstream s_commands(tmp_shell_filename.c_str());//matlab only    
  
  if (tool_name_ !="matlab")
    {
    ErrorMessage("Plot::show2Dcontourplot: 2D contour plot is availalbe only for matlab tool").throwMe();
    }
 
  //create binary data file
  string tmp_binary= cmd_filename_+".dat";
  FileMgr fp(tmp_binary,"wb");
  
  //--------------------------------------------------------
  //there are two different types of data for 2d plot:
  // (x,y,z2) where x and y are 1D vector and z is 2D matrix
  // (x2,y2,z2) where x2, y2, and z are all 2D matrix 
  // and only "Single set" of data should be present before 
  // making output file
  // Nset : count data set for (x,y,z2)
  // N2Dset: count data set for (x2,y2,z2)
  // Either one of  Nset and  N2Dset should be 1 and the other 0 
  //-------------------------------------------------------
  
  //add(x,y,z2):x and y are 1D array and z is 2D matrix of
  // (x.size(),y.size())
  unsigned int N1set = 0;
  Plot::CI py = y_.begin();
  Plot::DMI pz2 = z2_.begin();//set to the beginning
  
  for (Plot::CI px = x_.begin();
       px != x_.end() && py != y_.end() && pz2 != z2_.end() ;
       ++px, ++py,++pz2)
    {
      N1set += 1;
    }
  
  unsigned int N2set = 0;
  Plot::DMI py2 = y2_.begin();
  pz2= z2_.begin();
  for (Plot::DMI px2 = x2_.begin();
       px2 != x2_.end() && py2 != y2_.end() && pz2 != z2_.end() ;
       ++px2,++py2, ++pz2)
    {
      N2set += 1;
    }
  
  
  
  //count number of data set
  
  if (N1set == 0 &&  N2set == 0)
    {
      ErrorMessage("Plot::show2Dcontourplot: no data set").throwMe();
    }
  if (N1set !=0 && N2set !=0)
    {
      ErrorMessage("Plot::show2Dcontourplot: both types of  data set  (x,y,z2) and (x2,y2,z2) are present").throwMe();
    }
  
  if (N1set != 0)
    {
      fp.write(N1set);//number of data set
      py = y_.begin();
      pz2=z2_.begin();
      for (Plot::CI  px = x_.begin();
	   px != x_.end() && py != y_.end() && pz2 != z2_.end();
	   ++px, ++py,++pz2)
	{
	  unsigned int Nxdata = (*px).size();
	  fp.write(Nxdata);
	  unsigned int Nydata=(*py).size();
	  fp.write(Nydata);
	  for (unsigned int i=0; i < (*px).size(); ++i)
	    {
	      fp.write((*px)(i));	 
	    }
	  for (unsigned int i=0; i < (*py).size(); ++i)
	    {
	      fp.write((*py)(i));	 
	    }
	  for (unsigned int i = 0; i < (*px).size();++i)
	  for (unsigned int j = 0; j < (*py).size();++j)
	    {
	      fp.write((*pz2)(i,j));
	    }
	}
    }
  else if (N2set != 0)
    {
      fp.write(N2set);//number of data set
      py2 = y2_.begin();
      pz2 = z2_.begin();
      for (Plot::DMI px2 = x2_.begin();
	   px2 != x2_.end() && py2 != y2_.end() && pz2 != z2_.end();
	   ++px2, ++py2,++pz2)
	{
	  fp.write( (*px2).cols());
	  for (unsigned int i = 0; i < (*px2).cols();++i)
	  for (unsigned int j = 0; j < (*py2).cols();++j)
	    {
	      fp.write((*px2)(i,j));
	      fp.write((*py2)(i,j));
	      fp.write((*pz2)(i,j));
	    }
	}
    }
  else
    {
      ErrorMessage("Plot::show2Dcontourplot: no data present").throwMe();
    }
  //--------------------
  //write number of rectangles and their four cooridinate system
  //-------------------
  
  unsigned int Nrectangle = 0;//reset set value
  Plot::CI ry = rectangley_.begin();
  for (Plot::CI rx = rectanglex_.begin();
       rx != rectanglex_.end() && ry != rectangley_.end(); ++rx, ++ry)
    {
      Nrectangle += 1;
    }
  fp.write(Nrectangle);
  
  ry=rectangley_.begin();
  for (Plot::CI rx = rectanglex_.begin();
       rx != rectanglex_.end() && ry != rectangley_.end(); ++rx, ++ry)
    {
      for (unsigned int i =0; i < (*rx).size();++i)
	{
	  fp.write((*rx)(i));
	}
      for (unsigned int j = 0; j < (*ry).size();++j)
	{
	  fp.write((*ry)(j));
	}
    }
  
  //-------------------------------------------
  //add marks at (markx,marky)
  //-------------------------------------------    
  unsigned int Nmark  = 0;
  Plot::CI my=marky_.begin();
  for (Plot::CI mx=markx_.begin(); mx != markx_.end() && my !=marky_.end();++mx,++my)
    {
      Nmark += 1;
    }
  fp.write(Nmark);
  
  my=marky_.begin();
  for (Plot::CI mx = markx_.begin();
       mx != markx_.end() && my != marky_.end(); ++mx, ++my)
    {
      unsigned int Nxdata = (*mx).size();
      fp.write(Nxdata);
      for (unsigned int i =0; i < (*mx).size();++i)
	{
	  fp.write((*mx)(i));
	}
      for (unsigned int j = 0; j < (*my).size();++j)
	{
	  fp.write((*my)(j));
	}
    }
  
  fp.close();//close binary data file
  
  g_commands<<"fid=fopen('"<<tmp_binary<<"');"<<endl;//file open
  if (N1set != 0)
    {
      g_commands<<"N1set = fread(fid,1,'uint32');"<<endl;//read number of data
      g_commands<<"for i_set= 1: N1set "<<endl;
      g_commands<<"Nxdata= fread(fid,1,'uint32') ;"<<endl;
      g_commands<<"Nydata= fread(fid,1,'uint32') ;"<<endl;	
      g_commands<<" x = fread(fid,Nxdata,'double');"<<endl;
      g_commands<<" y = fread(fid,Nydata,'double');"<<endl;
      g_commands<<" for i=1:Nxdata "<<endl;
      g_commands<<" for j=1:Nydata "<<endl;
      g_commands<<" z(i,j) = fread(fid,1,'double');"<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" contourf(x,y,z',20) ; "<<endl;
      g_commands<<"hold on;"<<endl;
      g_commands<<"axis auto; "<<endl;
      g_commands<<" end "<<endl;
      g_commands<<"colorbar; "<<endl;
    }
  else if (N2set != 0)
    {
      g_commands<<"N2set = fread(fid,1,'uint32') ; "<<endl;//read number of data
      g_commands<<"for i_set= 1: N2set "<<endl;
      g_commands<<"Nxdata= fread(fid,1,'uint32') ;"<<endl;
      g_commands<<" for i=1:Nxdata "<<endl;
      g_commands<<" for j=1:Nxdata "<<endl;
      g_commands<<" x(i,j) = fread(fid,1,'double');"<<endl;
      g_commands<<" y(i,j) = fread(fid,1,'double');"<<endl;
      g_commands<<" z(i,j) = fread(fid,1,'double');"<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" contourf(x,y,z,20) ; "<<endl;
      g_commands<<"hold on;"<<endl;
      g_commands<<"axis auto; "<<endl;
      g_commands<<" end "<<endl;
      g_commands<<"colorbar; "<<endl;
    }
  
  //-------------------------
  //Read rectangle cooridinates
  //-------------------------
  g_commands<<"NsetRec=fread(fid,1,'uint32') ;"<<endl;
  g_commands<<" for i = 1:NsetRec "<<endl;
  //if NsetRec = 0,this part would not be executed
  g_commands<<" xcorner=fread(fid,4,'double');"<<endl;
  g_commands<<" ycorner=fread(fid,4,'double');"<<endl;
  g_commands<<" minx = min(xcorner); "<<endl;
  g_commands<<" maxx = max(xcorner); "<<endl;
  g_commands<<" miny = min(ycorner); "<<endl;
  g_commands<<" maxy = max(ycorner); "<<endl;
  g_commands<<" rectangle('Position',[minx,miny,maxx-minx,maxy-miny]);"<<endl;
  g_commands<<"  end           "<<endl;//end of for i = 1:Nset
  
  
  //---------------------------
  //read mark positions
  //----------------------------
  g_commands<<"NsetRec=fread(fid,1,'uint32'); "<<endl;
  g_commands<<"for i = 1: NsetRec "<<endl;
  g_commands<<"Nxdata=fread(fid,1,'uint32'); "<<endl;
  g_commands<<"x =fread(fid,Nxdata,'double');"<<endl;
  g_commands<<"y=fread(fid,Nxdata,'double');"<<endl;
  g_commands<<"plot(x,y,'.w','MarkerSize',8)"<<endl;
  g_commands<<"end "<<endl;//end of plotting marks 
  
  
  
  
  if(title_.size() > 0)
    {
      g_commands<<"title('"<<title_<<"');"<<endl;  
    }
  if (xlabel_.size() > 0)
    {
      g_commands<<"xlabel('"<<xlabel_<<"');"<<endl;   
    }       
  if(ylabel_.size() >0 )
    {
      g_commands<<"ylabel('"<<ylabel_<<"');"<<endl;   
    }
  if(target_device=="x")
    {//no matlab filename is specified
      g_commands<<"uiwait;"<<endl;//display and wait
    }
  else if (target_device=="matlab_fig")
    {
      g_commands<<"saveas(gcf,'"<<filename_<<"');"<<endl;//display and wait
    }
  else if (target_device=="matlab_eps")
    {
      g_commands<<"saveas(gcf,'"<<filename_<<"','psc2');"<<endl;
    }
  else
    {
      ErrorMessage("Plot.cpp::show2Dcontourplot: Invalid target device(x or matlab-fig)").throwMe();
    }
  g_commands<<"quit;"<<endl;
  g_commands << flush;
  
  string cmd_filename;
  unsigned int Nstring = cmd_filename_.size();
  cmd_filename.resize(Nstring-2);
  for (unsigned int i = 0; i<Nstring-2;i++)
    {
      cmd_filename += toStr(cmd_filename_[i]);
    }       
  s_commands<<"#! /bin/sh "<<endl;
  s_commands<<" matlab -nodesktop -nosplash -r "<<cmd_filename <<endl;
  s_commands<<flush;
  //---------------------------------
  //system command matlab -nodesktop -nosplash -r tm_file.m
  // did not work.  As a result, we decide to sidestep this problem
  // by making a shell command file named tmp.sh.  
  // content of tmp.sh
  //  #! /bin/sh
  //  matlab -nodesktop -nosplash -r cmd_filename
  //--------------------------------
  
  string buf;
  if(target_device !="x")
    {
      ErrorMessage("Plot::show2Dcontourplot: x is the only target_device avail");
    }
  buf ="sh " + tmp_shell_filename;
  if (system(buf.c_str()) != 0)
    {
      throw ErrorMessage("System command: " + buf + ", failed");
    }
  
  }


//---------------------------
//Show 2D color image
//--------------------------
void Plot::show2Dimage(const std::string& target_device)
  {
  std::ofstream g_commands(cmd_filename_.c_str());
  string tmp_shell_filename ="tmp_file.sh";//matlab only
  std::ofstream s_commands(tmp_shell_filename.c_str());//matlab only    
  
  if (tool_name_ !="matlab")
    {
    ErrorMessage("Plot::show2Dimage: 2D color image plot is availalbe only for matlab tool").throwMe();
    }
 
  //create binary data file
  string tmp_binary= cmd_filename_+".dat";
  FileMgr fp(tmp_binary,"wb");
  
  //--------------------------------------------------------
  //there are three different types of data for color image plot:
  //(x,y,z) where all of them are 1D array
  //(x,y,z2) where x and y are 1D vector and z is 2D matrix
  //(x2,y2,z2) where x2, y2, and z are all 2D matrix 
  // N0set : count data set for (x,y, z);
  // N1set : count data set for (x,y,z2)
  // N2set: count data set for (x2,y2,z2)
  // only one of  N0set ,N1set and  N2set should be 1 and the others 0 
  //-------------------------------------------------------
  
 
  //data set counter
  unsigned int N0set = 0;
  unsigned int N1set = 0;
  unsigned int N2set = 0;

  Plot::CI py = y_.begin();
  Plot::CI pz = z_.begin();
  Plot::DMI py2 = y2_.begin(); 
  Plot::DMI pz2 = z2_.begin();
  

  //count number of data set (x,y,z)
  for (Plot::CI px = x_.begin();
       px != x_.end() && py != y_.end() && pz != z_.end();
       ++px, ++py,++pz)
    {
      N0set += 1;
    }
  //return iterators
  py=y_.begin();
  pz=z_.begin();

  //count number of data set (x,y,z2)
  for (Plot::CI px = x_.begin();
       px != x_.end() && py != y_.end() && pz2 != z2_.end() ;
       ++px, ++py,++pz2)
    {
      N1set += 1;
    }
  py=y_.begin();
  pz2=z2_.begin();
 
  //count number of data set (x2,y2,z2)
  for (Plot::DMI px2 = x2_.begin();
       px2 != x2_.end() && py2 != y2_.end() && pz2 != z2_.end() ;
       ++px2,++py2, ++pz2)
    {
      N2set += 1;
    }
  py2=y2_.begin();
  pz2=z2_.begin();
  
  
  //count number of data set
  //make sure only a sinlge type of data set is present
  //otherwise it will throw an exception
  if (N0set==0 && N1set == 0 &&  N2set == 0)
    {
      ErrorMessage("Plot::show2Dimage: no data set").throwMe();
    }
  if (N0set != 0 && N1set !=0 && N2set !=0)
    {
      ErrorMessage("Plot::show2Dimage: all three types of  data set of (x,y,z),(x,y,z2) and (x2,y2,z2) are present").throwMe();
    }
  if ( (N0set != 0 && N1set != 0)
       ||(N1set !=0 && N2set !=0)
       ||(N2set !=0 && N0set !=0))
    {
      ErrorMessage("Plot::show2Dimage::at least two different types of data are present among (x,y,z),(x,y,z2), and (x2,y2,z2)").throwMe();
    }  
  if (N0set != 0)
    {
      //count number of data elements
      py = y_.begin();
      pz = z_.begin();
      unsigned int i_count = 0;
      for (Plot::CI  px = x_.begin();
	   px != x_.end() && py != y_.end() && pz != z_.end();
	   ++px, ++py,++pz)
	{
	  i_count += (*px).size();
	}
      fp.write(i_count);//save total number of data
      
      //write down a complete list of each component (x,y,z)
      for (Plot::CI px = x_.begin();px != x_.end();++px)
	{
	  for (unsigned int i = 0; i < (*px).size();++i)
	    {
	      fp.write((*px)(i));
	    }
	}

      for (py = y_.begin();py != y_.end();++py)
	{
	  for (unsigned int i = 0; i < (*py).size();++i)
	    {
	      fp.write((*py)(i));
	    }
	}
      py=y_.begin();

      for (pz = z_.begin();pz != z_.end();++pz)
	{
	  for (unsigned int i = 0; i < (*pz).size();++i)
	    {
	      fp.write((*pz)(i));
	    }
	}
      pz=z_.begin();
    }
  else if (N1set != 0)
    {
      fp.write(N1set);//number of data set
      py = y_.begin();
      pz2=z2_.begin();
      for (Plot::CI  px = x_.begin();
	   px != x_.end() && py != y_.end() && pz2 != z2_.end();
	   ++px, ++py,++pz2)
	{
	  unsigned int Nxdata = (*px).size();
	  fp.write(Nxdata);
	  unsigned int Nydata=(*py).size();
	  fp.write(Nydata);
	  for (unsigned int i=0; i < (*px).size(); ++i)
	    {
	      fp.write((*px)(i));	 
	    }
	  for (unsigned int i=0; i < (*py).size(); ++i)
	    {
	      fp.write((*py)(i));	 
	    }
	  for (unsigned int i = 0; i < (*px).size();++i)
	  for (unsigned int j = 0; j < (*py).size();++j)
	    {
	      fp.write((*pz2)(i,j));
	    }
	}
      py=y_.begin();
      pz2=z2_.begin();
    }
  else if (N2set != 0)
    {

      //count number of data elements
      py2 = y2_.begin();
      pz2 = z2_.begin();
      unsigned int i_count = 0;
      for (Plot::DMI  px2 = x2_.begin();
	   px2 != x2_.end() && py2 != y2_.end() && pz2 != z2_.end();
	   ++px2, ++py2,++pz2)
	{
	  i_count += (*px2).rows() * (*px2).cols();
	}
      fp.write(i_count);//save total number of data
    
      for (Plot::DMI px2 = x2_.begin(); px2 != x2_.end();++px2)
	{
	  
	  for (unsigned int i = 0; i < (*px2).rows();++i)
	  for (unsigned int j = 0; j < (*px2).cols();++j)
	    {
	      fp.write((*px2)(i,j));
	    }
	}   
      for ( py2 = y2_.begin(); py2 != y2_.end();++py2)
	{
	  
	  for (unsigned int i = 0; i < (*py2).rows();++i)
	  for (unsigned int j = 0; j < (*py2).cols();++j)
	    {
	      fp.write((*py2)(i,j));
	    }
	}
      py2=y2_.begin();
      for ( pz2 = z2_.begin(); pz2 != z2_.end();++pz2)
	{
	  
	  for (unsigned int i = 0; i < (*pz2).rows();++i)
	  for (unsigned int j = 0; j < (*pz2).cols();++j)
	    {
	      fp.write((*pz2)(i,j));
	    }
	}
      pz2=z2_.begin();
    }
  else
    {
      ErrorMessage("Plot::show2Dcontourplot: no data present").throwMe();
    }
 
  fp.close();//close binary data file
  
  g_commands<<"fid=fopen('"<<tmp_binary<<"');"<<endl;//file open
 
  if (N0set != 0)
    {
      g_commands<<"Ndata = fread(fid,1,'uint32');"<<endl;//read number of data
      g_commands<<" x = fread(fid,Ndata,'double');"<<endl;
      g_commands<<" y = fread(fid,Ndata,'double');"<<endl;
      g_commands<<" z = fread(fid,Ndata,'double');"<<endl;
      g_commands<<" imagesc(x,y,z) ; "<<endl;
      g_commands<<"hold on;"<<endl;
      g_commands<<"axis xy; "<<endl;
      g_commands<<"colorbar; "<<endl;
    }
  else if (N1set != 0)
    {
      g_commands<<"N1set = fread(fid,1,'uint32');"<<endl;//read number of data
      g_commands<<"for i_set= 1: N1set "<<endl;
      g_commands<<"Nxdata= fread(fid,1,'uint32') ;"<<endl;
      g_commands<<"Nydata= fread(fid,1,'uint32') ;"<<endl;	
      g_commands<<" x = fread(fid,Nxdata,'double');"<<endl;
      g_commands<<" y = fread(fid,Nydata,'double');"<<endl;
      g_commands<<" for i=1:Nxdata "<<endl;
      g_commands<<" for j=1:Nydata "<<endl;
      g_commands<<" z(i,j) = fread(fid,1,'double');"<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" end "<<endl;
      g_commands<<" imagesc(x,y,z') ; "<<endl;
      g_commands<<"axis auto; "<<endl;
      g_commands<<"axis xy; "<<endl;      
      g_commands<<"hold on;"<<endl;
      g_commands<<" end "<<endl;
      g_commands<<"colorbar; "<<endl;
    }
  else if (N2set != 0)
    {
      g_commands<<"N2set = fread(fid,1,'uint32');"<<endl;//read number of data
      g_commands<<" x= fread(fid,N2set,'double');"<<endl;
      g_commands<<" y = fread(fid,N2set,'double');"<<endl;
      g_commands<<" z = fread(fid,N2set,'double');"<<endl;
      g_commands<<" imagesc(x,y,z) ; "<<endl;
      g_commands<<"hold on;"<<endl;
      g_commands<<"axis auto; "<<endl;
      g_commands<<"axis xy; "<<endl;
      g_commands<<"colorbar; "<<endl;
    }
  
  if(title_.size() > 0)
    {
      g_commands<<"title('"<<title_<<"');"<<endl;  
    }
  if (xlabel_.size() > 0)
    {
      g_commands<<"xlabel('"<<xlabel_<<"');"<<endl;   
    }       
  if(ylabel_.size() >0 )
    {
      g_commands<<"ylabel('"<<ylabel_<<"');"<<endl;   
    }
  if(target_device=="x")
    {//no matlab filename is specified
      g_commands<<"uiwait;"<<endl;//display and wait
    }
  else if (target_device=="matlab_fig")
    {
      g_commands<<"saveas(gcf,'"<<filename_<<"');"<<endl;//display and wait
    }
  else if (target_device=="matlab_eps")
    {
      g_commands<<"saveas(gcf,'"<<filename_<<"','psc2');"<<endl;
    }
  else
    {
      ErrorMessage("Plot.cpp::show2Dimage: Invalid target device(x or matlab-fig)").throwMe();
    }
 
  g_commands<<"quit;"<<endl;
  g_commands << flush;
  
  string cmd_filename;
  unsigned int Nstring = cmd_filename_.size();
  cmd_filename.resize(Nstring-2);
  for (unsigned int i = 0; i<Nstring-2;i++)
    {
      cmd_filename += toStr(cmd_filename_[i]);
    }       
  s_commands<<"#! /bin/sh "<<endl;
  s_commands<<" matlab -nodesktop -nosplash -r "<<cmd_filename <<endl;
  s_commands<<flush;
  //---------------------------------
  //system command matlab -nodesktop -nosplash -r tm_file.m
  // did not work.  As a result, we decide to sidestep this problem
  // by making a shell command file named tmp.sh.  
  // content of tmp.sh
  //  #! /bin/sh
  //  matlab -nodesktop -nosplash -r cmd_filename
  //--------------------------------
  
  string buf;
  if(target_device !="x")
    {
      ErrorMessage("Plot::show2Dcontourplot: x is the only target_device avail");
    }
  buf ="sh " + tmp_shell_filename;
  if (system(buf.c_str()) != 0)
    {
      throw ErrorMessage("System command: " + buf + ", failed");
    }
  
  }

//--------------------
// Class line methods
//--------------------

//--------------
// Construction
//--------------

line::line() throw()
  : type("solid"), inttype(1), color("black"), intcolor(1), width(1.0) { }

line::line(const string& str1)
  : type("solid"), inttype(1), color("black"), intcolor(1), width(1.0)
  {
  setParam(str1);
  }

line::line(const string& str1, const float& wid)
  : type("solid"), inttype(1), color("black"), intcolor(1), width(wid)
  {
  setParam(str1);
  }

line::line(const string& str1, const string& str2)
  : type("solid"), inttype(1), color("black"), intcolor(1), width(1.0)
  {
  setParam(str1);
  setParam(str2);
  }

line::line(const string& str1, const string& str2, const float& wid)
  : type("solid"), inttype(1), color("black"), intcolor(1), width(wid)
  {
  setParam(str1);
  setParam(str2);
  }

void line::setParam(const string& str) throw(ErrorMessage)
  {
  if (str == "black")
    {
    color = str;
    intcolor = 1;
    }
  else if (str == "red")
    {
    color = str;
    intcolor = 2;
    }
  else if (str == "green")
    {
    color = str;
    intcolor = 3;
    }
  else if (str == "blue")
    {
    color = str;
    intcolor = 4;
    }
  else if (str == "yellow")
    {
    color = str;
    intcolor = 5;
    }
  else if (str == "brown")
    {
    color = str;
    intcolor=6;
    }
  else if (str == "gray")
    {
    color = str;
    intcolor=7;
    }
  else if (str == "violet")
    {
    color = str;
    intcolor=8;
    }
  else if (str == "cyan")
    {
    color = str;
    intcolor=9;
    }
  else if (str == "magenta")
    {
    color = str;
    intcolor=10;
    }
  else if (str == "orange")
    {
    color = str;
    intcolor=11;
    }
  else if (str == "none")
    {
    type = str;
    inttype = 0;
    }
  else if (str == "solid")
    {
    type = str;
    inttype = 1;
    }
  else if (str == "dash")
    {
    type = str;
    inttype = 2;
    }
  else if (str == "dot")
    {
    type = str;
    inttype = 3;
    }
  else if (str == "dot-dash")
    {
    type = str;
    inttype = 3;
    }
  else if (str == "long-dash")
    {
    type = str;
    inttype = 4;
    }
  else
    {
    throw ErrorMessage("Invalid line parameter: " + str);
    }
  }

//--------------------
// Class sym methods
//--------------------

//--------------
// Construction
//--------------

sym::sym() throw()
  : type("circle"), inttype(1), color("black"), intcolor(1), size(1.0), fill(0)
  { }

sym::sym(const string& str1)
  : type("circle"), inttype(1), color("black"), intcolor(1), size(1.0), fill(0)
  {
  setParam(str1);
  }

sym::sym(const string& str1, const float& siz)
  : type("circle"), inttype(1), color("black"), intcolor(1), size(siz), fill(0)
  {
  setParam(str1);
  }

sym::sym(const string& str1, const string& str2)
  : type("circle"), inttype(1), color("black"), intcolor(1), size(1.0), fill(0)
  {
  setParam(str1);
  setParam(str2);
  }

sym::sym(const string& str1, const string& str2, const float& siz)
  : type("circle"), inttype(1), color("black"), intcolor(1), size(siz), fill(0)
  {
  setParam(str1);
  setParam(str2);
  }

sym::sym(const string& str1, const string& str2, const float& siz,
  const int& fil)
  : type("circle"), inttype(1), color("black"), intcolor(1), size(siz),
    fill(fil)
  {
  setParam(str1);
  setParam(str2);
  }

// Symbol settings for Grace
void sym::setParam(const string& str) throw(ErrorMessage)
  {
  if (str == "black")
    {
    color = str;
    intcolor = 1;
    }
  else if (str == "red")
    {
    color = str;
    intcolor = 2;
    }
  else if (str == "green")
    {
    color = str;
    intcolor = 3;
    }
  else if (str == "blue")
    {
    color = str;
    intcolor = 4;
    }
  else if (str == "yellow")
    {
    color = str;
    intcolor = 5;
    }
  else if (str == "brown")
    {
    color = str;
    intcolor=6;
    }
  else if (str == "gray")
    {
    color = str;
    intcolor=7;
    }
  else if (str == "violet")
    {
    color = str;
    intcolor=8;
    }
  else if (str == "cyan")
    {
    color = str;
    intcolor=9;
    }
  else if (str == "magenta")
    {
    color = str;
    intcolor=10;
    }
  else if (str == "orange")
    {
    color = str;
    intcolor=11;
    }
  else if (str == "none")
    {
    type = str;
    inttype = 0;
    }
  else if (str == "plus")
    {
    type = str;
    inttype = 8;
    }
  else if (str == "circle")
    {
    type = str;
    inttype = 1;
    }
  else if (str == "square")
    {
    type = str;
    inttype = 2;
    }
  else if (str == "diamond")
    {
    type = str;
    inttype = 3;
    }
  else if (str == "up-triangle")
    {
    type = str;
    inttype = 4;
    }
  else if (str == "left-triangle")
    {
    type = str;
    inttype = 5;
    }
  else if (str == "right-triangle")
    {
    type = str;
    inttype = 7;
    }
  else if (str == "down-triangle")
    {
    type = str;
    inttype = 6;
    }
  else
    {
    ErrorMessage e("Invalid sym parameter: " + str);
    e.throwMe();
    }
  }

/****   old xmgr symbol settings
void sym::setParam(const string& str) throw(ErrorMessage)
  {
  if (str == "black")
    {
    color = str;
    intcolor = 1;
    }
  else if (str == "red")
    {
    color = str;
    intcolor = 2;
    }
  else if (str == "green")
    {
    color = str;
    intcolor = 3;
    }
  else if (str == "blue")
    {
    color = str;
    intcolor = 4;
    }
  else if (str == "yellow")
    {
    color = str;
    intcolor = 5;
    }
  else if (str == "brown")
    {
    color = str;
    intcolor=6;
    }
  else if (str == "gray")
    {
    color = str;
    intcolor=7;
    }
  else if (str == "violet")
    {
    color = str;
    intcolor=8;
    }
  else if (str == "cyan")
    {
    color = str;
    intcolor=9;
    }
  else if (str == "magenta")
    {
    color = str;
    intcolor=10;
    }
  else if (str == "orange")
    {
    color = str;
    intcolor=11;
    }
  else if (str == "none")
    {
    type = str;
    inttype = 0;
    }
  else if (str == "dot")
    {
    type = str;
    inttype = 1;
    }
  else if (str == "circle")
    {
    type = str;
    inttype = 2;
    }
  else if (str == "square")
    {
    type = str;
    inttype = 3;
    }
  else if (str == "diamond")
    {
    type = str;
    inttype = 4;
    }
  else if (str == "up-triangle")
    {
    type = str;
    inttype = 5;
    }
  else if (str == "left-triangle")
    {
    type = str;
    inttype = 6;
    }
  else if (str == "right-triangle")
    {
    type = str;
    inttype = 7;
    }
  else if (str == "down-triangle")
    {
    type = str;
    inttype = 8;
    }
  else
    {
    ErrorMessage e("Invalid sym parameter: " + str);
    e.throwMe();
    }
  }
***/

//----------------------
// Class legend methods
//----------------------

// Constructor
legend::legend()
  : fill_frame(true),frame(true),xpos(0.8,""),ypos(0.8,""),text_scale(1,"")
  { 
  return;
  }

// Clear
void legend::clear()
  {
  strlist.clear();
  }
