#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include "Error.h"
#include <string.h>
#include "PDSLabel.h"
#include "Utils.h"
#include "Constants.h"

using std::cerr;
using std::cout;
using std::endl;

#define BAD_VALUE_HEX 0xff7ffffb

using std::set_unexpected;
using std::terminate;

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
    int v=BAD_VALUE_HEX;
    float* p=(float*)&v;
    float BAD_VALUE =*p;

    //------------------------------
    // Parse Command Line
    //------------------------------
    int clidx=0;
    string command=argv[clidx++];
    if(argc!=8 && argc!=9){
      cerr << "Usage: " << command << " infile beamfile outfile mindist maxdist numdiststeps closestapproachrow [ucharvalued]"
      << endl;
      exit(1);
    }
    
    string infile=argv[clidx++];
    string beamfile=argv[clidx++];
    string outfile=argv[clidx++];
    float mindist=atof(argv[clidx++]);
    float maxdist=atof(argv[clidx++]);
    int N=atoi(argv[clidx++]);
    int zerorow=atoi(argv[clidx++]);
    bool intval=false;
    if(argc==9){
      intval=(bool)atoi(argv[clidx++]);
    }

    // open files
    FileMgr ifm(infile,"r");
    FileMgr bfm(beamfile,"r");


    // read label from input files
    PDSLabel inlabel(ifm);
    PDSLabel blabel(bfm);
    
    // compute file size info 
    int nrows=inlabel.fileRecords()-1;
    int ncols=inlabel.recordLength()/sizeof(float);
    int bheadsz=blabel.labelRecords();
    if(intval){
      nrows=inlabel.fileRecords()-bheadsz;
      ncols=inlabel.recordLength()/sizeof(unsigned char);
    }


    // if file is not PDS
    if(inlabel.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << infile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }




    // close input file
    ifm.close();
    bfm.close(); 

    // reopen files as standard FILE* files for fast read/write
    FILE* ifp=fopen(infile,"r");
    FILE* bfp=fopen(beamfile,"r");
    FILE* ofp=fopen(outfile,"w");


    
    // skip headers in BIDR (assumes 1 line in header)
    // for float valued BIDRs 
    // uchar valued cases have header size of beammask backplane
    if(!intval){
      fseek(ifp,ncols*sizeof(float),SEEK_SET);
    }
    else{
      fseek(ifp,ncols*sizeof(unsigned char)*bheadsz,SEEK_SET);
    }
    fseek(bfp,ncols*bheadsz*sizeof(char),SEEK_SET);



    // allocate buffers
    float dist_res=(maxdist-mindist)/N;
    float* minv=(float*)malloc(sizeof(float)*N*5);
    float* maxv=(float*)malloc(sizeof(float)*N*5);
    double* meanv=(double*)malloc(sizeof(double)*N*5);
    double* stdv=(double*)malloc(sizeof(double)*N*5);
    int* num=(int*)malloc(sizeof(int)*N*5);
  
    for(int i=0;i<N*5;i++){
      minv[i]=1e99;
      maxv[i]=-1e99;
      meanv[i]=0;
      stdv[i]=0;
      num[i]=0;
    }

    float* x=(float*)malloc(sizeof(float)*ncols);
    unsigned char* ix=(unsigned char*)malloc(sizeof(unsigned char)*ncols);
    char * bc=(char*)malloc(sizeof(char)*ncols);


    // loop over rows 
    for(int r=0;r<nrows;r++){
  
      // read row of data from each input file
      if(!intval){
	if((int)fread(&x[0],sizeof(float),ncols,ifp)!=ncols){
	  ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
			 +infile);
	  e.throwMe();
	}
      }
      else{
      // read row of integer data from each input file
	if((int)fread(&ix[0],sizeof(unsigned char),ncols,ifp)!=ncols){
	  ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
			 +infile);
	  e.throwMe();
	}
	for(int i=0;i<ncols;i++) x[i]=(float)ix[i];
      }
   
      float d=(r-zerorow)*2575*2*pi/360.0/256.0;
      printf("dist = %g\n",d);
      // read row of data from each input file
      if((int)fread(&bc[0],sizeof(char),ncols,bfp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
		       +beamfile);
	e.throwMe();
      }


      

      // loop over columns
      int bcvals[]={0,1,2,0,3,0,0,0,4,0,0,0,0,0,0,0,5};
      for(int c=0;c<ncols;c++){
        int b=0;
        if(bc[c]<=16) b=bcvals[(int)bc[c]];
        if(x[c]==BAD_VALUE || b==0 
	   || d<mindist || d>maxdist) continue;
	else{
	  int i=(int)(floor((d-mindist)/dist_res));
          int i2=(b-1)*N+i;
          meanv[i2]+=x[c];
          stdv[i2]+=x[c]*x[c];
          num[i2]++;
          if(x[c]>maxv[i2])maxv[i2]=x[c];
          if(x[c]<minv[i2])minv[i2]=x[c];
	}
      } // end column loop

      printf("Completed row %d.\n",r);
     } // end row loop

    // write output file
    fprintf(ofp,"# atdist,beam,mean,mean-2*std,mean+2*std,min,max,std,num\n");
    for(int i=0;i<N;i++){
      float d=(mindist+i*dist_res+0.5*dist_res);
      for(int b=0;b<5;b++){
	int c=N*b+i;
	if(num[c]>=2){
	  meanv[c]/=num[c];
	  stdv[c]=sqrt((stdv[c]-num[c]*meanv[c]*meanv[c])/(num[c]-1));
	  fprintf(ofp,"%g %d %g %g %g %g %g %g %d\n",
		  d,b+1,meanv[c],meanv[c]-2*stdv[c],meanv[c]+2*stdv[c],minv[c],maxv[c],stdv[c],num[c]);
	}
      }
    }

    // close files
    fclose(ifp);
    fclose(ofp);
    fclose(bfp);

    // free arrays
    free(x);
    free(ix);
    free(bc);
    free(num);
    free(meanv);
    free(maxv);
    free(minv);
    free(stdv);

  }
catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
    cerr << "Non Error Message Exception caught" << endl;
  }
}
