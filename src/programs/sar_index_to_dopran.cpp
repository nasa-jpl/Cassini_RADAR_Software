#include<stdio.h>
#include<stdlib.h>
#include"Utils.h"
#include<iostream>
#include"Error.h"
#include"L1I.h"

using std::cerr;
using std::endl;
using std::cout;

// this routine is used to compute the doppler and range used by the
// sar processor for given range and azimuth indices output by ptang
// it takes four command line parameters
// azimuth_index range_index doppler_file range_file

int main(int argc, char* argv[]){



  // allocate arrays
  float fdop[L1I_DOPPLER_WIDTH][L1I_RANGE_WIDTH];
  float ran[L1I_DOPPLER_WIDTH][L1I_RANGE_WIDTH];

  try{


  // parse command line
  int clidx=1;
  if(argc!=5){
    cerr<<"Usage:" << argv[0] 
	<< " azimuth_index range_index doppler_file range_file"
	<< endl;
    exit(1);
  }

  float azimuth_idx= atof(argv[clidx++]);
  float range_idx = atof(argv[clidx++]);
  string azimuth_file=argv[clidx++];
  string range_file=argv[clidx++];


  // open inputs file
  FILE* afp=fopen(azimuth_file,"r");
  FILE* rfp=fopen(range_file,"r");

 

  
  // read 
  for(int i=0;i<L1I_DOPPLER_WIDTH;i++){
    fread((void*)&ran[i][0],sizeof(float),L1I_RANGE_WIDTH,rfp);
  }  
  for(int i=0;i<L1I_DOPPLER_WIDTH;i++){
    fread((void*)&fdop[i][0],sizeof(float),L1I_RANGE_WIDTH,afp);
  }  

  // close input files
  fclose(afp);
  fclose(rfp);


  // interpolate range and doppler arrays
  float range;
  float doppler;

  int ri1,ri2,di1,di2;
  float cr1,cr2,cd1,cd2;
  
  // ptang indices start at 1 (center of first pixel)  
  // convert to starting edge of first pixel  = 0 notation
  range_idx--;
  azimuth_idx--;

  // get linear and bilinear interpolation coefficients and indices
  ri1=(int)floor(range_idx);
  ri2=ri1+1;

  di1=(int)floor(azimuth_idx);
  di2=di1+1;
  
  cr1=ri2-range_idx;
  cr2=1-cr1;

  cd1=di2-azimuth_idx;
  cd2=1-cd1;

  
  range=cr1*cd1*ran[di1][ri1]+cd1*cr2*ran[di1][ri2]
    +cr1*cd2*ran[di2][ri1]+cr2*cd2*ran[di2][ri2];

  doppler=cr1*cd1*fdop[di1][ri1]+cd1*cr2*fdop[di1][ri2]+cr1*cd2*fdop[di2][ri1]+
    cr2*cd2*fdop[di2][ri2];


  // output range and doppler
  cout.precision(20);
  cout << "range=" << range <<" km, doppler =" << doppler << " Hz" << endl;
  }
  catch( ErrorMessage e){
    cerr << e.msg << endl;
    exit(2);
  }
}
