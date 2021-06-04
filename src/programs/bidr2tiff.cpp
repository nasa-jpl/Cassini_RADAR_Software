// bidr2tiff
//
// Convert a BIDR file into a tiff image file.
//

#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <getopt.h>
#include <string.h>
#include <tiff.h>
#include <tiffio.h>
#include "Io.h"
#include "PDSLabel.h"
#include "BIDR.h"

using std::cerr;
using std::cout;
using std::endl;


using std::set_unexpected;
using std::terminate;

// Define command line option flags
#define OPTSTRING "hdni:o:g:s:l:r:"
extern int optind;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  terminate();
  }


int main(int argc, char* argv[])

{


set_unexpected(myUnexpected);
string input;

try
  {

  //------------------------
  // Parse the command line
  //------------------------

  char* infile = NULL;
  char* outfile = NULL;
  char* gamma_str = NULL;
  char* startline_str = NULL;
  char* numline_str = NULL;
  char* reduce_str = NULL;
  bool use_log_scaling = false;
  bool zero_collapse = false;

  while (1)
    {
    int c = getopt(argc,argv,OPTSTRING);
    if (c == 'h')
      {
      cout << "Option -h prints this message" << endl;
      cout << "Option -g specifies the gamma to apply, default 1.0 "
        << endl;
      cout << "Option -d if present specifies log scaling (dB), (-g ignored)"
        << endl;
      cout << "Option -n if negative values should be collapsed to zero"
        << endl;
      cout << "Option -s specifies the starting line number"
        << endl;
      cout << "Option -l specifies the length (number of lines to output)"
        << endl;
      cout << "Option -r specifies the reduction factor"
        << endl;
      cout << "Option -i specifies the input BIDR filename" << endl;
      cout << "Option -o specifies the output jpeg filename" << endl;
      }
    else if (c == 'g')
      {
      gamma_str = optarg; 
      }
    else if (c == 'd')
      {
      use_log_scaling = true;
      }
    else if (c == 'n')
      {
      zero_collapse = true;
      }
    else if (c=='o')
      {
      outfile = optarg;
      }
    else if (c == 'i')
      {
      infile = optarg;
      }
    else if (c == 's')
      {
      startline_str = optarg;
      }
    else if (c == 'l')
      {
      numline_str = optarg;
      }
    else if (c == 'r')
      {
      reduce_str = optarg;
      }
    else if (c == -1) break;
    }

  if (!infile || !outfile)
    {
    cout <<
      "Usage: bidr2tiff -h -i BIDR_Filename -o tiff_filename";
    cout << endl;
    return(-1);
    }

  float gamma = 1.0; // default value
  if (gamma_str)
    {
    gamma = atof(gamma_str);
    }

  int reduction = 1; // default value
  if (reduce_str)
    {
    reduction = atoi(reduce_str);
    }

  if (use_log_scaling)
    {
    cout << "using log scaling (proportional to dB) so gamma is not used"
      << endl;
    }
  else
    {
    cout << "gamma = " << gamma << endl;
    }

  int val=BAD_VALUE_HEX;
  float* ptr=(float*) &val;
  float bad_value=*ptr;

  // open input file and header output file
  FileMgr ifm(infile,"r");
  TIFF* ofm = TIFFOpen(outfile,"wb");

  // read label from input file
  PDSLabel label(ifm);

  // if file is not PDS
  if(label.label() == EMPTY_STRING){
    cerr << "Fatal error: " << infile 
      << " does not appear to be PDS." << endl;
    exit(1);
  }

  float value=0;
  int ncols=label.recordLength()/sizeof(int);
  int nrows=label.fileRecords()-label.labelRecords();
  cout << "rows = " << nrows << endl;
  cout << "cols = " << ncols << endl;

  int startline = 1; // default value
  if (startline_str)
    {
    startline = atoi(startline_str);
    }
  if (startline < 1) startline = 1;
  if (startline >= nrows) startline = nrows;

  int numline = nrows; // default value
  if (numline_str)
    {
    numline = atoi(numline_str);
    }
  if (numline < 1) numline = 1;
  if (numline >= nrows) numline = nrows;

  cout << "Outputting from row " << startline << " covering " << numline 
    << " rows " << endl;

  // offset over the PDS label
  int pos = ifm.getPosition();
  //cout << "pos = " << pos << endl;
  ifm.setPosition(pos + label.labelLength());

  // Pass through data collecting scaling info
  double maxvalue = 0;
  double minvalue = 1e6;
  double sumvalue = 0;
  double sumvalue2 = 0;
  double sumgvalue = 0;
  double sumgvalue2 = 0;
  double gvalue = 0;
  int nbad = 0;
  int first_nonempty_row_index = 0;
  int last_nonempty_row_index = 0;
  bool watch_start_empty = true;
  for (int i=0; i < nrows; ++i)
    {
    bool empty_row = true;
    for (int j=0; j < ncols; ++j)
      {
      ifm.read(value);
      if (value == bad_value)
        {
        nbad++;
        }
      else if (finite(value))
        {
        empty_row = false;
        if (zero_collapse & value < 0.0) value = 0.0;
        sumvalue += value;
        sumvalue2 += value*value;
        if (value > maxvalue) maxvalue = value;
        if (value < minvalue) minvalue = value;
        gvalue = pow(value,gamma);
        sumgvalue += gvalue;
        sumgvalue2 += gvalue*gvalue;
        }
      }
    if (empty_row && watch_start_empty)
      {
      first_nonempty_row_index = i+1;
      }
    else if (!empty_row)
      {
      watch_start_empty = false;
      last_nonempty_row_index = i;
      }
    }
  double meanvalue = sumvalue/(nrows*ncols - nbad);
  double mom2 = sumvalue2/(nrows*ncols - nbad);
  double sdev = sqrt(mom2 - meanvalue*meanvalue);
  double meangvalue = sumgvalue/(nrows*ncols - nbad);
  double gmom2 = sumgvalue2/(nrows*ncols - nbad);
  double gsdev = sqrt(gmom2 - meangvalue*meangvalue);
//  double gmaxvalue = pow(3.0*sdev+meanvalue,gamma);
  double gmaxvalue = meangvalue + 3.0*gsdev;
  double gminvalue = pow(minvalue,gamma); // a^b is monotonic in a
  double gm = 255.0/(gmaxvalue - gminvalue);
  double gb = -gm*gminvalue;
  double smaxvalue = log(maxvalue);
  double sminvalue = log(minvalue);
  double sm = 255.0/(smaxvalue - sminvalue);
  double sb = -sm*sminvalue;
  cout << "number of invalid pixels = " << nbad << endl;
  cout << "minvalue = " << minvalue << endl;
  cout << "maxvalue = " << maxvalue << endl;
  cout << "meanvalue = " << meanvalue << endl;
  //cout << "mom2 = " << mom2 << endl;
  cout << "sdev = " << sdev << endl;
  cout << "gminvalue = " << gminvalue << endl;
  cout << "gmaxvalue = " << gmaxvalue << endl;
  cout << "meangvalue = " << meangvalue << endl;
  cout << "gsdev = " << gsdev << endl;
  cout << "gm,gb = " << gm << ", " << gb << endl;
  //cout << "sminvalue = " << sminvalue << endl;
  //cout << "smaxvalue = " << smaxvalue << endl;
  //cout << "sm,sb = " << sm << ", " << sb << endl;

  ifm.setPosition(pos + label.labelLength());

  // Compute range of rows to output
  int istart,iend;
  if (first_nonempty_row_index > startline)
    {
    istart = first_nonempty_row_index;
    }
  else
    {
    istart = startline;
    }
  iend = istart + numline - 1;
  if (iend > last_nonempty_row_index) iend = last_nonempty_row_index;
  int Nrows_valid = iend - istart + 1;
//  int Nrows_valid = last_nonempty_row_index - first_nonempty_row_index + 1;

  // Apply subsampling factor
  if (reduction != 1)
    {
    cout << "Applying subsampling factor of: " << reduction << endl;
    }

  // Tiff setup

  uint16 config = PLANARCONFIG_CONTIG;
  uint32 rowsperstrip = (uint32) -1;
  uint32 rows = Nrows_valid/reduction;
  uint32 cols = ncols/reduction;

  TIFFSetField(ofm, TIFFTAG_IMAGELENGTH, rows);
  TIFFSetField(ofm, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(ofm, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(ofm, TIFFTAG_PLANARCONFIG, config);

  TIFFSetField(ofm, TIFFTAG_IMAGEWIDTH, cols);
  TIFFSetField(ofm, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(ofm, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

  tsize_t scanline = TIFFScanlineSize(ofm);
  if ((uint32)scanline != cols)
    {
    printf("scanline = %d, cols = %d\n",(int)scanline,(int)cols);
    }
  unsigned char* buf = (unsigned char *)_TIFFmalloc(scanline);
  rowsperstrip = 1;
  TIFFSetField(ofm, TIFFTAG_ROWSPERSTRIP,
      TIFFDefaultStripSize(ofm, rowsperstrip));

  for (int i=0; i < nrows; ++i)
    {
    cout << "row = " << i << endl;
    bool empty_row = true;
    for (int j=0; j < ncols; ++j)
      { // one row at a time
      ifm.read(value);
      if (j % reduction == 0)
        {
        int jj = j/reduction;
        if (value == bad_value)
          {
          // replace with white to save ink
          buf[jj] = 255;
          }
        else if (!finite(value))
          {
          // replace numerical screw-up's with black 
          buf[jj] = 0;
          }
        else
          {
          empty_row = false;
          if (use_log_scaling)
            {
            value = sm*log(value) + sb;
            }
          else
            {
            value = gm*pow(value,gamma) + gb;
            }
          if (value > 255.0) value = 255.0;
          if (value < 0.0) value = 0.0;
          buf[jj] = (unsigned char)value;
          }
        }
      }
//    if (i >= first_nonempty_row_index && i <= last_nonempty_row_index)
    if (i >= istart && i <= iend && i % reduction == 0)
      {
      if (TIFFWriteScanline(ofm, buf, i, 0) < 0) break;
      }
    }
  (void) TIFFClose(ofm);

  // close input file
  ifm.close();
  }
catch(...)
  {
  cerr << "Exception caught" << endl;
  }
}
