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

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include <getopt.h>
#include "BIDRFile.h"
#include "Utils.h"

static const char rcs_id_bidr_latlon_c[] =
  "@(#) $Id: bidr_to_outline.cpp,v 11.1 2012/03/23 18:33:01 bstiles Exp $";

#define PRECISION_MIN 1
#define PRECISION_MAX 9
#define BAD_VALUE_HEX 0xff7ffffb



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
    int v=BAD_VALUE_HEX;
    float* p=(float*)&v;
    float BAD_VALUE =*p;

    //------------------------
    // Parse the command line
    //------------------------  

    double line;
    double sample;
    double slat1,slat2;
    double slon1,slon2;
    double lat;
    double lon;
    int c;
    int precision = 9;
    bool UNITS_DEGREES = true; // false => radians
    bool POSITIVE_WEST = true; // true => positive-east longitudes
    bool PIXEL_TO_LATLON = true; // false => reverse conversion
    char *cmd = argv[0];
    char *bidr_input = NULL;
    string msg;

    
    int clidx=1;
    if(argc!=3){
      fprintf(stderr,"Usage %s  bidr_input_file  out_text_file\n",cmd);
      exit(1);
    }
    bidr_input=argv[clidx++];
    char* outfile=argv[clidx++];
     
    if (!fileExists(bidr_input))
    {
      msg = "Input file " + string(bidr_input) + " not found or unreadable";
      ErrorMessage e(msg);
      e.throwMe();
    }
    BIDRFile bf_in(bidr_input,"r");
    bf_in.readHeader();

    // open input file and header output file
    FileMgr ifm(bidr_input,"r");


    // read label from input file
    PDSLabel label(ifm);

    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< cmd << ": input file " << bidr_input
           << " does not appear to be PDS." << endl;
      exit(1);
    }


    // skip to correct position in file
    float value=0;

    int headrecs=label.labelRecords();
    int nrows=label.fileRecords()-headrecs;
    int offset =headrecs*label.recordLength();


    nrows=nrows-headrecs;
    int ncols=label.recordLength()/sizeof(float); // HACK assumes 4 byte records
    float *y = (float*)malloc(ncols*sizeof(float));

    ifm.close();

    FILE* ifp=fopen(bidr_input,"r");
    FILE* ofp=fopen(outfile,"w");

    if(ifp==NULL){
      fprintf(stderr,"Error opening BIDR input file %s\n",bidr_input);
      exit(1);
    }

    if(ofp==NULL){
      fprintf(stderr,"Error creating output text file %s\n",outfile);
      exit(1);
    }
    

    for(int c=0;c<nrows;c+=10)
    {
      fseek(ifp,(unsigned int)(offset+c*ncols*sizeof(float)),0);
      fread(&y[0],sizeof(float),ncols,ifp);
      int j=0;
      while(j<ncols && y[j]==BAD_VALUE) j++;
      if(j==ncols) continue;

      int j2=ncols-1;
      while(j2>j && y[j2]==BAD_VALUE) j2--;
  
      // Calculate the oblique cylindrical latitude and longitude in radians at
      // the selected point.
      lat = bf_in.OCLatFromGrid(j+1);
      lon = bf_in.OCLonFromGrid(c+1);

      // Convert the latitude and longitude to standard coordinates.
      slat1 = bf_in.standardLatInRad(lon, lat);
      slon1 = bf_in.standardLonInRad(lon, lat);


     // Calculate the oblique cylindrical latitude and longitude in radians at
      // the selected point.
      lat = bf_in.OCLatFromGrid(j2+1);
      lon = bf_in.OCLonFromGrid(c+1);

      // Convert the latitude and longitude to standard coordinates.
      slat2 = bf_in.standardLatInRad(lon, lat);
      slon2 = bf_in.standardLonInRad(lon, lat);
      if (POSITIVE_WEST)
      {
	slon1 = 2*pi - slon1;
	slon2 = 2*pi - slon2;
        
      }
      if (UNITS_DEGREES)
      {
	slat1 *= radtodeg;
	slon1 *= radtodeg;
	slat2 *= radtodeg;
	slon2 *= radtodeg;
      }

      // lon1 lat1 lon2 lat2 line_number
      fprintf(ofp,"%7.3f %7.3f %7.3f %7.3f %d\n",slon1,slat1,slon2,slat2,c+1);
    } // end of loop over rows
    free(y);
    fclose(ifp);
    fclose(ofp);
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


