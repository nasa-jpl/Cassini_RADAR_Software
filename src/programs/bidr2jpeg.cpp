// bidr2jpeg
//
// Convert a BIDR file into a jpeg image file.
//

#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <getopt.h>
#include <string.h>
#include <jpeglib.h>
#include "Io.h"
#include "PDSLabel.h"
#include "BIDR.h"
#include "Error.h"

using std::cerr;
using std::cout;
using std::endl;


using std::set_unexpected;
using std::terminate;

// Define command line option flags
#define OPTSTRING "hdneq:i:o:g:s:l:r:"
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
  char* quality_str = NULL;
  char* gamma_str = NULL;
  char* startline_str = NULL;
  char* numline_str = NULL;
  char* reduce_str = NULL;
  bool use_log_scaling = false;
  bool zero_collapse = false;
  bool empty_output = false;

  while (1)
    {
    int c = getopt(argc,argv,OPTSTRING);
    if (c == 'h')
      {
      cout << "Option -h prints this message" << endl;
      cout << "Option -q specifies the jpeg quality [0,100] default 75 "
        << endl;
      cout << "Option -g specifies the gamma to apply, default 1.0 "
        << endl;
      cout << "Option -d if present specifies log scaling (dB), (-g ignored)"
        << endl;
      cout << "Option -n if negative values should be collapsed to zero"
        << endl;
      cout << "Option -e if all empty rows should be output"
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
    else if (c == 'q')
      {
      quality_str = optarg; 
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
    else if (c == 'e')
      {
      empty_output = true;
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
      "Usage: bidr2jpeg -h -i BIDR_Filename -o jpeg_filename";
    cout << endl;
    return(-1);
    }

  int quality = 75; // default value
  if (quality_str)
    {
    quality = abs(atoi(quality_str));
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

  cout << "MAKE SURE YOU HAVE WRITE PERMISSION INTO DESTINATION DIRECTORY!"
    << endl << "  if not then this program will segment fault" << endl;
  cout << "quality = " << quality << endl;
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

  // Get actual file size
  long filebytes = 0;
  FILE* fileptr = fopen(infile,"rb");
  if (!fileptr)
    {
    cout << "Error opening " << infile << endl;
    }
  else
    {
    fseek(fileptr,0,SEEK_END);
    filebytes = ftell(fileptr);
    fclose(fileptr);
    }

  // open input file and header output file
  FileMgr ifm(infile,"r");
  cout << "Input file: " << infile << endl;
  FILE *ofm = fopen(outfile,"wb");
  cout << "Output file: " << outfile << endl;

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

  if (filebytes/(4*ncols) < nrows)
    {
    nrows = filebytes/(4*ncols);
    cout << "partial file, actual rows = " << nrows << endl;
    }

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
  if (first_nonempty_row_index > startline && !empty_output)
    {
    istart = first_nonempty_row_index;
    }
  else
    {
    istart = startline;
    }
  iend = istart + numline - 1;
  if (iend > last_nonempty_row_index && !empty_output)
    {
    iend = last_nonempty_row_index;
    }
  int Nrows_valid = iend - istart + 1;
//  int Nrows_valid = last_nonempty_row_index - first_nonempty_row_index + 1;

  // Apply average down factor
  if (reduction != 1)
    {
    cout << "Applying average down factor of: " << reduction << endl;
    }

  // Initialize the JPEG decompression object with default error handling.
  // Structure declarations moved here to avoid some kind of memory corruption
  // problem when they were at the top.  Can't see what the problem was.
  unsigned char* buffer[1];
  buffer[0] = (unsigned char*)malloc(3*(ncols/reduction));
  float* linebuf = (float *)malloc(sizeof(float)*ncols/reduction);
  int* nlinebuf = (int *)malloc(sizeof(int)*ncols/reduction);
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  // Specify data destination for compression
  jpeg_stdio_dest(&cinfo, ofm);

  // set basic image parameters
  cinfo.image_width = ncols/reduction;
  cinfo.image_height = Nrows_valid/reduction;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
      
  // set default compression parameters
  jpeg_set_defaults(&cinfo);

  // set user supplied quality parameter
  jpeg_set_quality(&cinfo,quality,TRUE);

  // Start compressor
  (void) jpeg_start_compress(&cinfo,TRUE);

  // Initialize line buffers
  for (int jj=0; jj < ncols/reduction; ++jj)
    {
    linebuf[jj] = bad_value;
    nlinebuf[jj] = 0;
    }

  for (int i=0; i < nrows; ++i)
    {
    for (int j=0; j < ncols; ++j)
      { // one row at a time
      ifm.read(value);
      int jj = j/reduction;
      if (finite(value) && value != bad_value)
        {  // sum together rows and cols within reduction factor
        if (nlinebuf[jj] == 0)
          {  // first entry - overwrite bad_value
          linebuf[jj] = value;
          }
        else
          {
          linebuf[jj] += value;
          }
        nlinebuf[jj] += 1;
        }
      }
    if (i >= istart && i <= iend && (i+1) % reduction == 0)
      {  // output current line buffer and setup for the next
      for (int jj=0; jj < ncols/reduction; ++jj)
        {
        if (linebuf[jj] == bad_value)
          {
          // replace with white to save ink
          value = 255;
          }
        else if (use_log_scaling)
          {
          value = sm*log(linebuf[jj]/nlinebuf[jj]) + sb;
          }
        else
          {
          value = gm*pow(linebuf[jj]/nlinebuf[jj],gamma) + gb;
          }
        if (value > 255.0) value = 255.0;
        if (value < 0.0) value = 0.0;
        buffer[0][3*jj] = (unsigned char)value;
        buffer[0][3*jj+1] = (unsigned char)value;
        buffer[0][3*jj+2] = (unsigned char)value;
        linebuf[jj] = bad_value;
        nlinebuf[jj] = 0;
        }
      (void)jpeg_write_scanlines(&cinfo, buffer, 1);
      }
    }

  (void) jpeg_finish_compress(&cinfo);
  fclose(ofm);
  jpeg_destroy_compress(&cinfo);

  // free temporary buffers
  free(buffer[0]);
  free(linebuf);
  free(nlinebuf);

  // close input file
  ifm.close();
  }

catch(ErrorMessage& e)
  {
  cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
  cerr << "Exception caught" << endl;
  }

}
