#include "BIDRFile.h"
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include "Error.h"
#include "PDSLabel.h"

using std::cerr;
using std::endl;
using std::cout;

int main(int argc, char* argv[])
{
  try{
  int clidx=0;
  string command=argv[clidx++];
  if(argc!=5){
    cerr << "Usage:"<< command << " bidr_infile bidr_outfile factor use_db"
	 << endl;
    exit(1);
  }
  string bf1_str=argv[clidx++];
  string of_str=argv[clidx++];
  int factor=atoi(argv[clidx++]);
  bool use_db=(bool)atoi(argv[clidx++]);

  BIDRFile bf1(bf1_str,"r");
  FileMgr of(of_str,"w");

  // read input headers output first header and header of output file
  bf1.readHeader();

  // for now no label is output to output file

  // check for error conditions
  // for now only floating point types are down averaged
  
  BIDRTypeE typ1;
  typ1=bf1.type();
 
  ErrorMessage e(command + " bad BIDR type");

  switch(typ1){
    case INC:
    case S0_UNCORR:
    case S0_CORR:
    case LAT:
    case LON:
      break;
  default:
    e.throwMe();
  }

  
  // Main loop
  int num_lons=bf1.numLons();
  int num_lats=bf1.numLats();
  float* input_array=(float*)malloc(sizeof(float)*num_lats*factor);  
  float* output_array=(float*)malloc(sizeof(float)*num_lats/factor);  

  float bad=bf1.BAD_VALUE;

  for(int i=0;i<num_lons/factor;i++){
    // read 
    for(int j=0; j< num_lats*factor; j++){
      bf1.read(input_array[j]);
    }
    // down average
    for(int j=0; j< num_lats/factor;j++){
      output_array[j]=0.0;
      int n=0;
      for(int k=0;k<factor;k++){
	for(int m=0;m<factor;m++){
	  int offset=k*num_lats+j*factor+m;
	  float val=input_array[offset];
	  if(val!=bad){
	    n++;
	    output_array[j]+=val;
	  }
	}
      }
      if(n){
	output_array[j]/=n;
	if(use_db) output_array[j]=10*log10(output_array[j]);
      }
      else output_array[j]=bad;
    }
    //write
    for(int j=0; j< num_lats/factor;j++){
      of.write(output_array[j]);
    }
  }
  free(output_array);
  free(input_array);
  cout << "Number of columns is " << num_lats/factor << endl;
  }

  catch( ErrorMessage e){
    cerr << e.msg;
    exit(1);
  }
}
