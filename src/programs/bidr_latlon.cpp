//----------------------------------------------------------------------
// bidr_latlon:  Given a standard latitude and longitude in a BIDR file,
//               display the corresponding line and sample numbers.
//
// bidr_pixel:   Given line and sample numbers in a BIDR file, display
//               the corresponding standard latitude and longitude.
//
// If the executable is invoked as "bidr_latlon," it performs the
// conversion to lat/lon.  Otherwise, it performs the conversion to
// line/sample.
//
// Both conversions rely on the oblique cylindrical projection parameters
// stored in the PDS label of the BIDR file.  The results for the conversions
// have inherent inaccuracies because of the limited precision used in the
// label.  The lat/lon outputs are good to about 7 significant digits,
// while the line/sample outputs are good only to about 4 significant
// digits.
//----------------------------------------------------------------------

#include <stdlib.h>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include <string.h>
#include <getopt.h>
#include "BIDRFile.h"
#include "Utils.h"

static const char rcs_id_bidr_latlon_c[] =
  "@(#) $Id: bidr_latlon.cpp,v 11.6 2017/04/11 20:35:14 cveerama Exp $";

#define PRECISION_MIN 1
#define PRECISION_MAX 9

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::fixed;
using std::setprecision;

void usage(char *cmd, bool flag);

int main(int argc, char* argv[])
{
  try
  {

    //------------------------
    // Parse the command line
    //------------------------  

    double line;
    double sample;
    double slat;
    double slon;
    double lat;
    double lon;
    int c;
    int precision = 6;
    int option_index;
    bool UNITS_DEGREES = true; // false => radians
    bool POSITIVE_WEST = false; // true => positive-east longitudes
    bool PIXEL_TO_LATLON = true; // false => reverse conversion
    char *cmd = argv[0];
    char *bidr_input = NULL;
    extern char *optarg;
    extern int optind;
    string msg;
    struct option long_options [] =
    {
      {"pre", 1, 0, 10},
      {"pw", 0, 0, 11},
      {"rad", 0, 0, 12},
      {0, 0, 0, 0}
    };

    // HACK that allows negative operands
    // if the command line argument list is too short for any options
    // then it does not look for any
    if(argc>4){
      while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
	{
	  switch (c)
	    {
	    case 10:
	      precision = atoi(optarg);
	      break;
	    case 11:
	      POSITIVE_WEST = true;
	      break;
	    case 12:
	      UNITS_DEGREES = false;
	      break;
	    default:
	      ;
	      break;
	    }
	}
    }
    string s = cmd;
    size_t index = s.find("bidr_latlon", 0);
    if (index == string::npos)
    {
      PIXEL_TO_LATLON = false;
    }

    if (argc - optind != 3)
    {
      usage(cmd, PIXEL_TO_LATLON);
    }
    bidr_input = argv[optind++];
    if (PIXEL_TO_LATLON)
    {
      line = strtod(argv[optind++], (char **) NULL);
      sample = strtod(argv[optind++], (char **) NULL);
    }
    else
    {
      slat = strtod(argv[optind++], (char **) NULL);
      slon = strtod(argv[optind++], (char **) NULL);
    }

    // Check options and arguments
    if (!fileExists(bidr_input))
    {
      msg = "Input file " + string(bidr_input) + " not found or unreadable";
      ErrorMessage e(msg);
      e.throwMe();
    }
    BIDRFile bf_in(bidr_input,"r");
    bf_in.readHeader();

    if (PIXEL_TO_LATLON)
    {
      if (line < 1.0 || line > (double) bf_in.numLons())
      {
	msg = "Invalid line " + toStr(line) + ".  Line must be between" +
	  " 1 and " + toStr(bf_in.numLons()) + ".";
	ErrorMessage e(msg);
	e.throwMe();
      }
      if (sample < 1.0 || sample > (double) bf_in.numLats())
      {
	msg = "Invalid sample " + toStr(sample) + ".  Sample must be between" +
	  " 1 and " + toStr(bf_in.numLats()) + ".";
	ErrorMessage e(msg);
	e.throwMe();
      }
    }
    if (precision < PRECISION_MIN || precision > PRECISION_MAX)
    {
      msg = "Invalid precision " + toStr(precision) + ".  Precision must be" +
        " between " + toStr(PRECISION_MIN) + " and " + toStr(PRECISION_MAX) + ".";
      ErrorMessage e(msg);
      e.throwMe();
    }

    if (PIXEL_TO_LATLON)
    {
      // Calculate the oblique cylindrical latitude and longitude in radians at
      // the selected point.
      lat = bf_in.OCLatFromGrid(sample);
      lon = bf_in.OCLonFromGrid(line);

      // Convert the latitude and longitude to standard coordinates.
      slat = bf_in.standardLatInRad(lon, lat);
      slon = bf_in.standardLonInRad(lon, lat);
      if (POSITIVE_WEST)
      {
	slon = 2*pi - slon;
      }
      if (UNITS_DEGREES)
      {
	slat *= radtodeg;
	slon *= radtodeg;
      }
      cout << setprecision(precision) << std::fixed << slat << " " << slon << endl;
    }
    else
    {
      // Convert the latitude and longitude to oblique cylindrical coordinates.
      // First make sure that the values are in radians and that the longitude
      // is in positive-east format.
      if (UNITS_DEGREES)
      {
	slon *= degtorad;
	slat *= degtorad;
      }
      if (POSITIVE_WEST)
      {
	while (slon < 0) { slon += 2*pi; }
	while (slon >= 2*pi) { slon -= 2*pi; }
	slon = 2*pi - slon;
      }
      lat = bf_in.latInRad(slon, slat);
      lon = bf_in.lonInRad(slon, slat);

      // Calculate the (fractional) line and sample at those coordinates.
      int ppd=bf_in.pixelsPerDegree();
      line = bf_in.OCLonToGrid(lon);
      int absmaxnumline=360*ppd;
      while(line>absmaxnumline) line-=absmaxnumline;

      sample = bf_in.OCLatToGrid(lat);
      cout << setprecision(precision) << std::fixed << line << " " << sample << endl;
    }
  }
  catch(ErrorMessage& e)
  {
    cerr << "Error:  " << e.msg << endl;
  }
  catch(std::exception& ex)
  {
    cerr << typeid(ex).name() << ": " << ex.what() << endl;
  }
  catch(...)
  {
    cerr << "Exception caught" << endl;
  }
}

void usage(char *cmd, bool pixel_to_latlon)
{
  if (pixel_to_latlon)
  {
    cerr << "Usage:  " << cmd << " [--pre <precision>] [--pw] [--rad] " <<
      "<BIDR file> <line> <sample>" << endl;
    cerr << "outputs the standard latitude and longitude at the given line" <<
      " and sample," << endl;
    cerr << "where" << endl;
    cerr << "  --pre specifies the number of digits after the decimal point " <<
      "(default 6)" << endl;
    cerr << "  --pw outputs the longitude in positive-west coordinates " <<
      "(default positive-east)" << endl;
    cerr << "  --rad outputs both values in radians (default degrees)" << endl;
    cerr << "Lines and samples are numbered starting at 1, with fractional values allowed." << endl;
  }
  else
  {
    cerr << "Usage:  " << cmd << " [--pre <precision>] [--pw] [--rad] " <<
      "<BIDR file> <latitude> <longitude>" << endl;
    cerr << "outputs the line and sample at the given standard latitude" <<
      " and longitude," << endl;
    cerr << "where" << endl;
    cerr << "  --pre specifies the number of digits after the decimal point " <<
      "(default 6)" << endl;
    cerr << "  --pw reads the longitude in positive-west coordinates " <<
      "(default positive-east)" << endl;
    cerr << "  --rad reads both values in radians (default degrees)" << endl;
    cerr << "Lines and samples are numbered starting at 1." << endl;
  }
  exit(1);
}
