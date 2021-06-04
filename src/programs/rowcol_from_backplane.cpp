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

  try{


  // parse command line
  int clidx=1;
  if(argc!=6){
    cerr<<"Usage:" << argv[0] 
	<< " row col nrows ncols infile"
	<< endl;
    exit(1);
  }

  double row= strtod(argv[clidx++],NULL);
  double col =strtod(argv[clidx++],NULL);
  int nrows =atoi(argv[clidx++]);
  int ncols =atoi(argv[clidx++]);
  string infile=argv[clidx++];



  // ptang indices start at 1 (center of first pixel)  
  // convert to starting edge of first pixel  = 0 notation
  row--;
  col--;

  if( row < 0 || row >= nrows || col < 0 || col >= ncols ){
    cerr << argv[0] <<": Bad row or column" << endl;
  }
  // open inputs file
  FileMgr f(infile,"r");
 
  float v00,v01,v10,v11;
  f.setPosition(4*ncols*int(row)+4*int(col));
  f.read(v00);
  f.read(v01);
  f.setPosition(4*ncols*int(row+1)+4*int(col));
  f.read(v10);
  f.read(v11);
  



  // interpolate backplane
  double value;
  

  int ri1,ri2,di1,di2;
  double cr1,cr2,cd1,cd2;
  


  // get linear and bilinear interpolation coefficients and indices
  ri1=(int)floor(row);
  ri2=ri1+1;

  di1=(int)floor(col);
  di2=di1+1;
  
  cr1=ri2-row;
  cr2=1-cr1;

  cd1=di2-col;
  cd2=1-cd1;

  
  value=cr1*cd1*(double)v00+cd2*cr1*(double)v01
    +cr2*cd1*(double)v10+cr2*cd2*(double)v11;


  // output value
  cout.precision(20);
  cout << value << endl;
  }
  catch( ErrorMessage e){
    cerr << e.msg << endl;
    exit(2);
  }
}
