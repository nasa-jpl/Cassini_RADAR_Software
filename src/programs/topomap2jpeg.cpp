// topomap2jpeg
//
// Convert a Topomap file into a jpeg image file.
//

#include <math.h>
#include <stdlib.h>
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
//#define OPTSTRING "hdneq:i:o:g:s:l:r:"
//extern int optind;

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

  
  int clidx=1;
  if(argc!=6 && argc!=8 && (argc<10 || argc> 14)){
    fprintf(stderr,"topomap2jpeg topomapfile bidrmapfile numrows numcols outfile [mindB] [maxdB] [minht] [maxht] [nvertgrid] [nhorizgrid] [nleggrid] [ispolarmap]\n");
    exit(1);
  }
  char * tfile=argv[clidx++];
  char * bfile=argv[clidx++];
  int nrows=atoi(argv[clidx++]);
  int ncols=atoi(argv[clidx++]);
  int nrowout=nrows*21/20;
  char * outfile=argv[clidx++];
  float mindb=-25;
  float maxdb=7;
  float minht=-1500;
  float maxht=1000;
  int nvertgrid=1;
  int nhorizgrid=1;
  int nleggrid=1;
  int ispolar=0;
  if(argc>=8){
    mindb=atof(argv[clidx++]);
    maxdb=atof(argv[clidx++]);
  }
  if(argc>=10){
    minht=atof(argv[clidx++]);
    maxht=atof(argv[clidx++]);
  }
  if(argc>=11) nvertgrid=atoi(argv[clidx++]);
  if(argc>=12) nhorizgrid=atoi(argv[clidx++]);
  if(argc>=13) nleggrid=atoi(argv[clidx++]);
  if(argc>=14) ispolar=atoi(argv[clidx++]);
  if(ispolar) nrowout=nrows*22/20;
  int vstep=nrows/nvertgrid;
  int hstep=ncols/nhorizgrid;
  int lstep=ncols/nleggrid;

  int quality = 90; // default value

  cout << "MAKE SURE YOU HAVE WRITE PERMISSION INTO DESTINATION DIRECTORY!"
    << endl << "  if not then this program will segment fault" << endl;
  cout << "quality = " << quality << endl;

  int val=BAD_VALUE_HEX;
  float* ptr=(float*) &val;
  float bad_value=*ptr;

  // Get actual file size
  int filebytes = 0;
  FILE* fileptr = fopen(tfile,"rb");
  if (!fileptr)
    {
    cout << "Error opening " << tfile << endl;
    }
  else
    {
    fseek(fileptr,0,SEEK_END);
    filebytes = ftell(fileptr);
    fclose(fileptr);
    }

  // open files
  FileMgr tfm(tfile,"r");
  FileMgr bfm(bfile,"r");
  cout << "Topomap file: " << tfile << endl;
  cout << "BIDR mosaic file: " << bfile << endl;
  FILE *ofm = fopen(outfile,"wb");
  cout << "Output file: " << outfile << endl;

  float value=0,h=0,s0=0;
  int cmapidx, red,green, blue;
  float rmap[64]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,32,48,64,80,96,112,128,144,160,176,192,208,224,240,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,224,208,192,176,160,144};
  float gmap[64]={0,0,0,0,0,0,0,0,16,32,48,64,80,96,112,128,144,160,176,192,208,224,240,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,224,208,192,176,160,144,128,112,96,80,64,48,32,16,0,0,0,0,0,0,0,0};
  float bmap[64]={144,160,176,192,208,224,240,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,224,208,192,176,160,144,128,112,96,80,64,48,32,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  cout << "rows = " << nrows << endl;
  cout << "cols = " << ncols << endl;

  if (filebytes/(4*ncols) < nrows)
    {
    nrows = filebytes/(4*ncols);
    cout << "Fatal error: partial file, actual rows = " << nrows << endl;
    exit(1);
    }

  int startline = 1; // default value
  /**** comment out for now
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
  */



  // Initialize the JPEG decompression object with default error handling.
  // Structure declarations moved here to avoid some kind of memory corruption
  // problem when they were at the top.  Can't see what the problem was.
  unsigned char* buffer[1];
  buffer[0] = (unsigned char*)malloc(3*(ncols));
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  // Specify data destination for compression
  jpeg_stdio_dest(&cinfo, ofm);

  // set basic image parameters
  cinfo.image_width = ncols;
  cinfo.image_height = nrowout;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
      
  // set default compression parameters
  jpeg_set_defaults(&cinfo);

  // set user supplied quality parameter
  jpeg_set_quality(&cinfo,quality,TRUE);

  // Start compressor
  (void) jpeg_start_compress(&cinfo,TRUE);

 
  for (int i=0; i < nrows; ++i)
    {
    for (int j=0; j < ncols; ++j)      { // one row at a time
     
      float di=i-(nrows-1)/2.0;
      float dj=j-(ncols-1)/2.0;
      float r=(int)sqrt(di*di+dj*dj);
      float ang=atan2(di/r,dj/r);
      if(ang<0)ang=ang+2*pi;

      tfm.read(h);
      bfm.read(s0);
      // for polar case white out everything below min lat
      
      if(ispolar && r>nrows/2){
	 buffer[0][3*j] = 255;
	 buffer[0][3*j+1] = 255;
	 buffer[0][3*j+2] = 255;
      }
      else{
	int s0bad=0;
	int hbad=0;
	if(s0<-999999 || !finite(s0)){
	  value=1.0;
	  s0bad=1;
	}
	else{
	  if(s0<0) s0=0.0001;
	  value=10*log10(s0);
	  value=(value-mindb)/(maxdb-mindb);
	  if(value>1) value=1;
	  if(value<0) value=0;
	}
	if(h<-999999) hbad=1;

	if(hbad){
	  red=(int)(255*value);
	  green=(int)(255*value);
	  blue=(int)(255*value);
	}
	else{
	  cmapidx=(int)(64*(h-minht)/(maxht-minht));
	  if(cmapidx<0) cmapidx=0;
	  if(cmapidx>63) cmapidx=63;
	  red=(int)(rmap[cmapidx]);
	  green=(int)(gmap[cmapidx]);
	  blue=(int)(bmap[cmapidx]);
	}
	buffer[0][3*j] = (unsigned char)red;
	buffer[0][3*j+1] = (unsigned char)green;
	buffer[0][3*j+2] = (unsigned char)blue;
  
	// Add gridlines
	if(!ispolar){
	  if(i<=1 || i >= nrows-2 || j<=1 || j>=ncols-2 || j%hstep<=1 || i%vstep<=1){
	    buffer[0][3*j] = 0;
	    buffer[0][3*j+1] = 0;
	    buffer[0][3*j+2] = 0;
	  }
	}
	// add polar gridlines
	else{
	  // radial lines
	  for(int ir=0;ir<=nvertgrid;ir++){
	    float dr=r-(nrows-1)*ir/(2.0*nvertgrid);
	    if(fabs(dr)<2){
	      buffer[0][3*j] = 0;
	      buffer[0][3*j+1] = 0;
	      buffer[0][3*j+2] = 0;
	    }
	  }
	  for(int ia=0;ia<=nhorizgrid;ia++){
	    float dang=ang-(2.0*pi*ia)/nhorizgrid;
	    if(fabs(dang)*r<2){
	      buffer[0][3*j] = 0;
	      buffer[0][3*j+1] = 0;
	      buffer[0][3*j+2] = 0;
	    }
	  }
	} // end polar gridlines case 
      } // end valid region case
    }
    (void)jpeg_write_scanlines(&cinfo, buffer, 1);
    }
    
  for (int i=nrows; i < nrowout; ++i)
    {
    for (int j=0; j < ncols; ++j)      { // one row at a time
      cmapidx=(int)(64*(j)/float(ncols));
        if(cmapidx<0) cmapidx=0;
        if(cmapidx>63) cmapidx=63;
        red=(int)(rmap[cmapidx]);
        green=(int)(gmap[cmapidx]);
        blue=(int)(bmap[cmapidx]);
      
      buffer[0][3*j] = (unsigned char)red;
      buffer[0][3*j+1] = (unsigned char)green;
      buffer[0][3*j+2] = (unsigned char)blue;
      if(i<=nrows+1 || i >= nrowout-2 || j<=1 || j>=ncols-2 || j%lstep<=1){
	buffer[0][3*j] = 0;
	buffer[0][3*j+1] = 0;
	buffer[0][3*j+2] = 0;
      }
      // leave space for polar legend
      if(ispolar && i < (nrows+nrowout)/2 ){
	buffer[0][3*j] = 255;
	buffer[0][3*j+1] = 255;
	buffer[0][3*j+2] = 255;
      }
      if(ispolar && i >= (nrows+nrowout)/2 && i <= (nrows+nrowout)/2 + 1 ){
	buffer[0][3*j] = 255;
	buffer[0][3*j+1] = 255;
	buffer[0][3*j+2] = 255;
      }
   }
 
  (void)jpeg_write_scanlines(&cinfo, buffer, 1);
    }
  (void) jpeg_finish_compress(&cinfo);
  fclose(ofm);
  jpeg_destroy_compress(&cinfo);

  // free temporary buffers
  free(buffer[0]);

  // close input file
  tfm.close();
  bfm.close();
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
