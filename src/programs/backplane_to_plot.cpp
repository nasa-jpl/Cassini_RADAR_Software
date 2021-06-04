#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include "Error.h"
#include <string.h>
#include "PDSLabel.h"
#include "Utils.h"

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
    if(argc!=5 && argc!=6){
      cerr << "Usage: " << command << " infile beamfile outfile beam [ucharvalued]"
      << endl;
      exit(1);
    }
    
    string infile=argv[clidx++];
    string beamfile=argv[clidx++];
    string outfile=argv[clidx++];
    int beam=atoi(argv[clidx++]);
    bool intval=false;
    if(argc==6){
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
    float* minv=(float*)malloc(sizeof(float)*nrows);
    float* maxv=(float*)malloc(sizeof(float)*nrows);
    double* meanv=(double*)malloc(sizeof(double)*nrows);
    double* stdv=(double*)malloc(sizeof(double)*nrows);
    int* num=(int*)malloc(sizeof(int)*nrows);
  
    for(int i=0;i<nrows;i++){
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
        if(x[c]==BAD_VALUE || b!=beam) continue;
	else{
          meanv[r]+=x[c];
          stdv[r]+=x[c]*x[c];
          num[r]++;
          if(x[c]>maxv[r])maxv[r]=x[c];
          if(x[c]<minv[r])minv[r]=x[c];
	}
      } // end column loop


      // write output file
      if(num[r]>=2){
	meanv[r]/=num[r];
	stdv[r]=sqrt((stdv[r]-num[r]*meanv[r]*meanv[r])/(num[r]-1));
	fprintf(ofp,"%d %g %g %g %g %g %g %d\n",
		r,meanv[r],meanv[r]-stdv[r],meanv[r]+stdv[r],minv[r],maxv[r],stdv[r],num[r]);
      }
      printf("Completed row %d.\n",r);
     } // end row loop



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
