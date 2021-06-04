#include<stdio.h>
#include<stdlib.h>
#include"Utils.h"
#include<iostream>
#include"Error.h"
#include"L1I.h"
#include"SimpleArray.h"

using std::cerr;
using std::endl;
using std::cout;

// this routine is used to compute the doppler and range used by the
// sar processor for given range and azimuth indices output by ptang
// it takes four command line parameters
// xvalue yvalue rows cols ximagefile yimagefile
 

int main(int argc, char* argv[]){



  try{


  // parse command line
  int clidx=1;
  if(argc!=7){
    cerr<<"Usage:" << argv[0] 
	<< "xvalue yvalue rows cols ximagefile yimagefile "
	<< endl;
    exit(1);
  }

  float x= atof(argv[clidx++]);
  float y= atof(argv[clidx++]);
  int rows= atoi(argv[clidx++]);
  int cols= atoi(argv[clidx++]);
  string xfile=argv[clidx++];
  string yfile=argv[clidx++];



  // allocate arrays
  float** xs=(float**)make_array(sizeof(float),2,rows,cols);
  float** ys=(float**)make_array(sizeof(float),2,rows,cols);
 


  // open inputs file
  FILE* xfp=fopen(xfile,"r");
  FILE* yfp=fopen(yfile,"r");

 

  
  // read 
  for(int i=0;i<rows;i++){
    fread((void*)&xs[i][0],sizeof(float),cols,xfp);
  }  
  for(int i=0;i<rows;i++){
    fread((void*)&ys[i][0],sizeof(float),cols,yfp);
  }  

  // close input files
  fclose(xfp);
  fclose(yfp);


  // find best pixel
 int ibest=0,jbest=0;
 float mindist=HUGE_VAL;
 float dist = 0; 
 float xdist,ydist;
 float dxr=xs[1][0]-xs[0][0];
 float dxc=xs[0][1]-xs[0][0];
 float dyr=ys[1][0]-ys[0][0];
 float dyc=ys[0][1]-ys[0][0];
 float dx=sqrt(dxr*dxr+dxc*dxc);
 float dy=sqrt(dyr*dyr+dyc*dyc);
 for(int i=0;i<rows;i++){
   for(int j=0;j<cols;j++){
     xdist = (x -xs[i][j])/dx;
     xdist*=xdist;
     ydist = (y -ys[i][j])/dy;
     ydist*=ydist;
     dist = xdist + ydist;
     if(dist < mindist){
       ibest=i;
       jbest=j;
       mindist=dist;
     }
   }
 }
 
 // output pixel i,j, mindist
 cout.precision(20);
 cout << ibest+1 << " " << jbest+1 << " " << mindist  << endl;

 free_array((void*)xs,2,rows,cols);
 free_array((void*)ys,2,rows,cols);
  }
  catch( ErrorMessage e){
    cerr << e.msg << endl;
    exit(2);
  }
}
