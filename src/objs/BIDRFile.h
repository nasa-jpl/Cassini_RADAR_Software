#ifndef BIDRFILE_H
#define BIDRFILE_H

static const char rcs_id_bidrfile_h[] =
  "@(#) $Id: BIDRFile.h,v 11.5 2011/09/16 00:03:29 richw Exp $";

#include <map>
#include "BIDR.h"
#include "Constants.h"
#include "Io.h"
#include "PDSLabel.h"
#include "Plot.h"
#include "SARFunctions.h"

#define BAD_VALUE_HEX 0xff7ffffb
// Default resolution for "averaged" beam mask product
#define DEFAULT_NEW_BEAMMASK_RESOLUTION 128
// Marker for "--ppd" option in process_beammask
#define PPD_OPTION_MARKER 12

//enum BIDRTypeE {INC, BEAMMASK, S0_UNCORR, S0_CORR, LAT, LON, START_BURST_NUM, END_BURST_NUM, NUM_LOOKS, S0_CORR_DB, INVALID};
enum BIDROptionErrorE { RESOLUTION_INVALID = -10, RESOLUTION_TOO_HIGH };

class BIDRFile : public FileMgr
{
 public:
   BIDRFile();
   BIDRFile(const string& filename, const string& mode);
   ~BIDRFile();
   void slideShow(int lons_per_plot, int lon_skip, bool plot_average);
   string& typeString();
   /*** These may not be needed
   float minLat();
   float maxLat();
   float minLon();
   float maxLon();
   int numLons();
   int numLats();
   ***/

   void readLine();
   void readHeader();
   void fixPole(double tca);
   void skipInvalidLines(bool b){skip_invalid_=b;}
   void skipLowValuedLines(bool b){skip_low_=b;}
   void useDecibels(bool b){use_dB_=b;}
   void setRange(float minlon,float maxlon){min_lon_=minlon; max_lon_=maxlon;}
   BIDRTypeE getBIDRType();
   int getScaleFactor(int new_resolution);

   float BAD_VALUE;
   int numLons(){return num_lons_;}
   int numLats(){return num_lats_;}
   int pixelsPerDegree(){return pixelsperdegree;}
   int labelLength(){return label_length_;}
   BIDRTypeE type(){return bte_;}
   string productID(){return product_id_;}
   string label(){return label_;}
   KeyValSeq labelList(){return label_list_;}
   int *invalidValue(){return (int *) &invalid_val_;}
   double standardLatInRad(double lonrad, double latrad);
   double standardLonInRad(double lonrad, double latrad);
   double OCLatFromGrid(double sample);
   double OCLonFromGrid(double line);
   double latInRad(double slonrad, double slatrad);
   double lonInRad(double slonrad, double slatrad);
   double OCLatToGrid(double latrad);
   double OCLonToGrid(double lonrad);

 private:

   void setDataValid();
   bool header_handled_;
   BIDRTypeE bte_;
   string type_string_;
   string y_units_;
   string product_id_;
   float min_lat_;
   int pixelsperdegree;
   float res_;
   float max_lat_;
   float start_lon_;
   float low_val_;
   float invalid_val_;
   bool lon_inc_;
   double d_min_lat_;
   double d_start_lon_;
   double d_res_;
   double pole_latitude_;
   double pole_longitude_;
   double pole_rotation_;
   double rotation_matrix_[3][3];
   double rotation_matrix_iau[3][3];
   double rotation_matrix_6para[3][3];
   
   int num_lats_;
   int num_lons_;
   Dvec x_;
   Dvec y_;
   Dvec y2_;
   Plot plot_;
   int current_lon_idx_;
   bool skip_invalid_;
   bool skip_low_;
   bool use_dB_;
   bool data_valid_;
   float min_lon_;
   float max_lon_;
   string label_;
   int label_length_;
   KeyValSeq label_list_;
};

class BIDRFileUpdates
{
 public:
   BIDRFileUpdates(char *input_file, char *output_file,
     int flagged_ppd_new, int ppd_new);
   ~BIDRFileUpdates();
   int average();
   int flip();
   int readInputToBuffer();
   void averageBuffer(unsigned char **cbuf);
   void averageBuffer(int **ibuf);
   void averageBuffer(float **fbuf);
   void flipBuffer(char **cbuf);
   void flipBuffer(int **ibuf);
   void flipBuffer(float **fbuf);
   void setProductID();
   void updateLabel(int pixelsize=1,bool flipped=false);
   void updateFlippedLabel(int pixelsize=1);
   int writeBufferToOutput();
 private:
   BIDRFile bf_in_;
   BIDRFile bf_out_;
   char **buf_;
   int line_start_;
   int line_end_;
   int pixel_start_;
   int pixel_end_;
   int orig_line_start_;
   int orig_line_end_;
   int orig_pixel_start_;
   int orig_pixel_end_;
   int scale_factor_;
   int buf_numlats_;
   int buf_numlons_;
   int checksum_;
   int ppd_;
   int truncated_;
   string product_id_;
   KeyValSeq label_list_;
};
#endif
