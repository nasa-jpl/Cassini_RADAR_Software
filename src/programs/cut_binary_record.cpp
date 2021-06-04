#include<stdio.h>
#include<stdlib.h>
#include"Utils.h"
#include<iostream>
#include"Error.h"

using std::cerr;
using std::endl;


// this routine is used to clip a single binary record from a file
// it takes four command line parameters
// record_size, record_num (starting from 1), in_file, and out_file
// a final optional header size parameter is also allowed

int main(int argc, char* argv[]){

  try{
  // parse command line
  int clidx=1;
  if(argc!=5 && argc!=6){
    cerr<<"Usage:" << argv[0] 
	<< " record_size record_num(first=1) in_file out_file [header_size]"
	<< endl;
    exit(1);
  }

  unsigned int record_size = atoi(argv[clidx++]);
  unsigned int record_num = atoi(argv[clidx++]);
  string in_file=argv[clidx++];
  string out_file=argv[clidx++];
  int header_size=0;
  if(argc==6){
    header_size=atoi(argv[clidx++]);
  }

  // open input file
  FILE* ifp=fopen(in_file,"r");

  // open output file
  FILE* ofp=fopen(out_file,"w");

  // skip to start of record
  unsigned int offset=header_size+(record_num-1)*record_size;
  fseek(ifp,offset,0);
  
  // read
  char* bfr=(char*)malloc(record_size);
  fread(&bfr[0],1,record_size,ifp);
  
  // write
  fwrite(&bfr[0],1,record_size,ofp);

  // close files
  fclose(ifp);
  fclose(ofp);

  free(bfr);
  }
  catch( ErrorMessage e){
    cerr << e.msg << endl;
    exit(2);
  }
}
