#include <iomanip>
#include <stdlib.h>
#include <sstream>
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
    if(argc!=5 && argc!=4 ){
      cerr << "Usage: " << command << " neigborhood_size infile kpfile [meanfile]"
      << endl;
      exit(1);
    }
    
    int n=atoi(argv[clidx++]);
    if(n%2==0 || n<1){
      cerr << command << 
	": Neighborhood_size should be an odd positive integer" << endl;
      exit(1);
    }
    
    string infile=argv[clidx++];
    string kpfile=argv[clidx++];
    bool compute_mean = (argc==5);
    string meanfile;
    if(compute_mean) meanfile=argv[clidx++];

    // open input file and header output file
    FileMgr ifm(infile,"r");


    // read label from input file
    PDSLabel label(ifm);
    
    // compute file size info (ASSUMES FLOATING POINT FORMAT)
    int nrows=label.fileRecords()-1;
    int ncols=label.recordLength()/sizeof(float);


    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << infile 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }

    // write labels to output files
    FileMgr kpfm(kpfile,"w");
    kpfm.write(label.label());
    if(compute_mean){
          FileMgr meanfm(meanfile,"w");
          meanfm.write(label.label());
          meanfm.close();
    }

    // close files
    kpfm.close();
    ifm.close();

    // reopen files as standard FILE* files for fast read/write
    FILE* ifp=fopen(infile,"r");
    FILE* kpfp=fopen(kpfile,"r+");
    FILE* meanfp=NULL;
    if(compute_mean) meanfp=fopen(meanfile,"r+");

    
    // skip headers in BIDR (assumes 1 line in header)
    fseek(ifp,ncols*sizeof(float),SEEK_SET);
    fseek(kpfp,ncols*sizeof(float),SEEK_SET);
    if(compute_mean) fseek(meanfp,ncols*sizeof(float),SEEK_SET);

    // allocate result buffer
    float* kp=(float*)malloc(sizeof(float)*ncols);
    float* mean=(float*)malloc(sizeof(float)*ncols);

    // allocate ring buffer for reading file
    float* value=(float*)malloc(sizeof(float)*ncols*n);
    

    // read first n-1 rows
    if((int)fread(&value[0],sizeof(float),(ncols*(n-1)),ifp)!=(ncols*(n-1))){
      ErrorMessage e("Error Reading Rows 1 to "+toStr(n-1)+" from file "
		     +infile);
      e.throwMe();
    }

    // write first n/2 bad_value rows 
    for(int c=0;c<ncols;c++) kp[c]=BAD_VALUE;
    for(int i=0;i<n/2;i++){
      if((int)fwrite(&kp[0],sizeof(float),ncols,kpfp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(i+1)+" to file "
		     +kpfile);
	e.throwMe();
      }
      if(compute_mean){
        // kp (all BAD_VALUES) is written to mean file 
	if((int)fwrite(&kp[0],sizeof(float),ncols,meanfp)!=ncols){
	  ErrorMessage e("Error writing Row "+toStr(i+1)+" to file "
			 +meanfile);
	  e.throwMe();
	}
      }
    }


    // loop over rows with full neighborhood r= n/2 to nrows-(n/2)-1
    for(int r=(n/2);r<nrows-(n/2);r++){
      int rowmod=(r+n/2-1)%n;
      
      // read row from file
      if((int)fread(&value[ncols*rowmod],sizeof(float),ncols,ifp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
		       +infile);
	e.throwMe();
      }
      
      float kp_norm = 1.0/((n*n)-1);
      float mean_norm = 1.0/(n*n);

      // loop over columns
      for(int c=0;c<ncols;c++){

	// compute row neighborhood bounds
	int startcol= c-n/2;
	int endcol = c+n/2;
        
        // set edge columns in which full neighborhood is not available
	// to BAD_VAL
	if(startcol<0 || endcol>=ncols){
	  mean[c]=BAD_VALUE;
	  kp[c]=BAD_VALUE;
	}
        // otherwise compute mean and kp from neighborhood
        else{
          bool badval_found=false;
	  mean[c]=0;
	  kp[c]=0;
          for(int i=0;i<n;i++){
	    int offset=i*ncols;
	    for(int j=startcol;j<=endcol;j++){
              float val = value[offset+j];
              if(val==BAD_VALUE){
		badval_found=true;
	      }
	      mean[c]+=val;
	    } // end neighborhood column loop
            if(badval_found) break;
	  } // end neighborhood row loop
	  if(badval_found){
	    mean[c]=BAD_VALUE;
	    kp[c]=BAD_VALUE;
	  }
	  else{
	    mean[c]*=mean_norm;
	    for(int i=0;i<n;i++){
	      int offset=i*ncols;
	      for(int j=startcol;j<=endcol;j++){
		float val = (value[offset+j]-mean[c]);
		kp[c]+=val*val;
	      } // end neighborhood column loop
	    } // end neighborhood row loop
            kp[c]*=kp_norm;
            kp[c]=sqrt(kp[c]);
            kp[c]=kp[c]/mean[c];
	  }
	} // end else (FULL NEIGHBORHOOD AVAILABLE case)
      } // end column loop

      // write output rows

      if((int)fwrite(&kp[0],sizeof(float),ncols,kpfp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(r+1)+" to file "
		     +kpfile);
	e.throwMe();
      }
      if(compute_mean){
        // kp (all BAD_VALUES) is written to mean file 
	if((int)fwrite(&mean[0],sizeof(float),ncols,meanfp)!=ncols){
	  ErrorMessage e("Error writing Row "+toStr(r+1)+" to file "
			 +meanfile);
	  e.throwMe();
	}
      }

      cout << "finished row " << r <<endl;
    } // end row loop

    // write last n/2 bad_value rows 
    for(int c=0;c<ncols;c++) kp[c]=BAD_VALUE;
    for(int i=0;i<n/2;i++){
      if((int)fwrite(&kp[0],sizeof(float),ncols,kpfp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(nrows-(n/2)+i)+" to file "
		     +kpfile);
	e.throwMe();
      }
      if(compute_mean){
        // kp (all BAD_VALUES) is written to mean file 
	if((int)fwrite(&kp[0],sizeof(float),ncols,meanfp)!=ncols){
	  ErrorMessage e("Error writing Row "+toStr(nrows-(n/2)+i)+" to file "
			 +meanfile);
	  e.throwMe();
	}
      }
    }

    // close files
    fclose(ifp);
    fclose(kpfp);
    if(meanfp)fclose(meanfp);

    // free arrays
    free(value);
    free(kp);
    free(mean);

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
