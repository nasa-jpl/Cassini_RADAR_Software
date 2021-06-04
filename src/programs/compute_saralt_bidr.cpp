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

#define MAX_NUM_SABS 60000
int main(int argc, char* argv[])

{
set_unexpected(myUnexpected);
float ht[MAX_NUM_SABS];
try
  {
    int v=BAD_VALUE_HEX;
    float* p=(float*)&v;
    float BAD_VALUE =*p;
    for(int i=0;i<MAX_NUM_SABS;i++)ht[i]=-1000;

    //------------------------------
    // Parse Command Line
    //------------------------------
    int clidx=0;
    string command=argv[clidx++];
    if(argc!=6 ){
      cerr << "Usage: " << command << " filebase s0ext estrdanal_file htoffset htscale"
      << endl;
      exit(1);
    }

    
    string base = argv[clidx++];
    string s0ext = argv[clidx++];
    string eafile = argv[clidx++];
    float htoff=atof(argv[clidx++]);
    float htsc=atof(argv[clidx++]);
    string s0file = base + "." + s0ext;
    string startfile = base + ".start_burst_num";
    string endfile = base + ".end_burst_num";
    string outfile = base + ".saralt";

    // open files
    FileMgr s0fm(s0file,"r");
    FileMgr ofm(outfile,"w");

    PDSLabel label(s0fm);      

    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << s0file 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }


  
    int nrows=label.fileRecords()-1;
    int ncols=label.recordLength()/sizeof(float);
    ofm.write(label.label());
    ofm.close();
    s0fm.close();

    //read estrdanal file
   FileMgr efm(eafile,"r");
   string str;
   str.resize(200);
   efm.readLine(str);
   efm.readLine(str);
    while(!efm.eof()){
      str=efm.parseNextString();
      str=efm.parseNextString();
      int sabc=atoi(str.c_str());
      for(int i=3;i<=10;i++) str=efm.parseNextString();
      ht[sabc]=atof(str.c_str());
      for(int i=11;i<=15;i++) str=efm.parseNextString();      
    }  
    efm.close();

    // reopen files as standard FILE* files for fast read/write
    FILE* s0fp=fopen(s0file,"r");
    FILE* sfp=fopen(startfile,"r");

    int use_start_only=0;
    FILE* efp;
    try{
      efp=fopen(endfile,"r");
    }
    catch(ErrorMessage e){
      efp=fopen(startfile,"r");
      fprintf(stderr,"endfile not found use start_sab as end_sab\n");
    }
    FILE* ofp=fopen(outfile,"r+");


    
    // skip headers in BIDR (assumes 1 line in header)
    fseek(s0fp,ncols*sizeof(float),SEEK_SET);
    fseek(sfp,ncols*sizeof(float),SEEK_SET);
    fseek(efp,ncols*sizeof(float),SEEK_SET);
    fseek(ofp,ncols*sizeof(float),SEEK_SET);


    // allocate buffers
    float* s0=(float*)malloc(sizeof(float)*ncols);
    float* y=(float*)malloc(sizeof(float)*ncols);
    int* sab1=(int*)malloc(sizeof(int)*ncols);
    int* sab2=(int*)malloc(sizeof(int)*ncols);


    // loop over rows 
    for(int r=0;r<nrows;r++){
  
      // read row of data from each input file
      if((int)fread(&s0[0],sizeof(float),ncols,s0fp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "
		       +s0file);
	e.throwMe();
      }

      if((int)fread(&sab1[0],sizeof(int),ncols,sfp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "+
		       startfile);
	e.throwMe();
      }

     if((int)fread(&sab2[0],sizeof(int),ncols,efp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(r+1)+" from file "+
		       endfile);
	e.throwMe();
      }

      

      // loop over columns
      for(int c=0;c<ncols;c++){
        int sab;
        if(sab1[c]<0 || sab2[c]<0) y[c]=BAD_VALUE;
        else{
	  if((sab2[c]-sab1[c])%2==0) sab=(sab1[c]+sab2[c])/2;
	  else sab=(sab1[c]+sab2[c]-5)/2;
	  if(ht[sab]>-500) y[c]=htsc*(ht[sab]-htoff);
	  else y[c]=s0[c];
	}
      } // end column loop

      // write output rows

      if((int)fwrite(&y[0],sizeof(float),ncols,ofp)!=ncols){
	ErrorMessage e("Error writing Row "+toStr(r+1)+" to file "
		     +outfile);
	e.throwMe();
      }

      cout << "finished row " << r <<endl;
    } // end row loop


    // close files
    fclose(ofp);
    fclose(sfp);
    fclose(efp);
    fclose(s0fp);

    // free arrays
    free(y);
    free(sab1);
    free(sab2);
    free(s0);


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
