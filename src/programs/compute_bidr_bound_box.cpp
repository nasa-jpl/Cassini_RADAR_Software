//----------------------------------------------------------------------
//  Compute latitude min/max and eastermnost and westernmost longitude
//----------------------------------------------------------------------
double ANGDIF(double A,double B){
  double d=A-B;
  if(d>180) d=d-360;
  if(d<-180) d=d+360;
  return(d);
}

#include <stdlib.h>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include <getopt.h>
#include "BIDRFile.h"
#include "Utils.h"

static const char rcs_id_compute_bidr_bound_box_c[] =
  "@(#) $Id: compute_bidr_bound_box.cpp,v 11.6 2017/04/07 22:09:41 cveerama Exp $";

#define PRECISION_MIN 1
#define PRECISION_MAX 9

using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::fixed;
using std::setprecision;



int main(int argc, char* argv[])
{
  try
  {

    //------------------------
    // Parse the command line
    //------------------------  

 
    double slat;
    double slon;
    double lat;
    double lon;
    int c;
    int precision = 8;
    int option_index;
    bool UNITS_DEGREES = true; // false => radians
    bool POSITIVE_WEST = true; 
    bool PIXEL_TO_LATLON = true; // false => reverse conversion
    char *cmd = argv[0];
    char *bidr_input = NULL;
    int clidx=1;
    string msg;

    if(argc!=2){
      cerr << "Usage: compute_bidr_bound_box bidrfile" << endl;
      exit(1);
    } 
    bidr_input = argv[clidx++];


    // Check options and arguments
    if (!fileExists(bidr_input))
    {
      msg = "Input file " + string(bidr_input) + " not found or unreadable";
      ErrorMessage e(msg);
      e.throwMe();
    }
    BIDRFile bf_in(bidr_input,"r");
    bf_in.readHeader();

    int nsamp=0;
    double minlat=90.0, maxlat=-90.0, elon=0, wlon=0, angbet=0;
    for(int line=1;line<=bf_in.numLons();line++){
      for(int sample=1;sample<=bf_in.numLats();sample++){
      // Calculate the oblique cylindrical latitude and longitude in radians at
      // the selected point.
      lat = bf_in.OCLatFromGrid(sample);
      lon = bf_in.OCLonFromGrid(line);

      // Convert the latitude and longitude to standard coordinates.
      slat = bf_in.standardLatInRad(lon, lat);
      slon = bf_in.standardLonInRad(lon, lat);
      slon = 2*pi - slon;
      slat *= radtodeg;
      slon *= radtodeg;
     
      //--------- find min max lat
      if(slat<minlat) minlat=slat;
      if(slat>maxlat) maxlat=slat;

      //--------- find easternmost and westernmost lon 
      // set elon first time through  
      //#define DEBUG
      if( nsamp== 0 ) elon = slon;

      // second time through set wlon (skip samples until slon changes)
      // also set angbet
      else if( nsamp== 1 ){
	double dlon = ANGDIF(slon,elon);
	if( dlon> 0 ){
	  wlon = slon;
          angbet=ANGDIF(wlon,elon);
#ifdef DEBUG
	  cerr << elon << " " << wlon << " " << angbet << " " << nsamp << "DLON=" << dlon << endl;
#endif
	}
	else if( dlon !=0 ){
	  wlon = elon;
	  elon = slon;
          angbet=ANGDIF(wlon,elon);
#ifdef DEBUG
	  cerr << elon << " " << wlon << " " << angbet << " " << nsamp << "DLON=" << dlon << endl;
#endif
	}
        else{
	  nsamp=0; // if the first two long are exactly equal wait for next one
	}    
      }
      // after first two samples update elon, wlon, and angbet
      // by expanding the the minimal amount to include each new point
      else{
	double eang=ANGDIF(slon,elon);
        if(eang>angbet){
	  wlon=slon;
	  angbet=eang;
	}
	if(eang<0 && (angbet-eang)<360){
	  elon=slon;
	  angbet=angbet-eang;
	}
#ifdef DEBUG
	  cerr << elon << " " << wlon << " " << angbet << " " << nsamp << "EANG=" << eang << endl;
#endif
      }

        nsamp++; 
      }

    }
    if(angbet>355){
      wlon=360;
      elon=0;
    }
    cout << setprecision(precision) << std::fixed << minlat << " " << maxlat << " " << elon << " " << wlon << endl;
    return(0);
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


