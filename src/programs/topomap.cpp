//----------------------------------------------
//  topomap
// Routine for producing a floating point topo map from sar data
// routine reads an intermed file and an aux file
// produced by sar_topo to get lines, sample, height, dheight etc
//
// This routine uses a bidr file to transform line/sample to lat lon
// The topomap.cpp routine gets lat and lons from the intermed file.


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "BIDR.h"
#include "Utils.h"
#include "Config.h"
#include "TopoMap.h"
#include "SimpleArray.h"

using std::cerr;
using std::endl;
int readfilelist(string lname, int& nfiles, char**& fnames, int*& typs, double*& times){
  FILE* lfp=fopen(lname,"r");
  char line[200];
  fnames=(char**)malloc(sizeof(char*)*1000);
  typs=(int*)malloc(sizeof(int)*1000);
  times=(double*)malloc(sizeof(double)*1000);
  nfiles=0;
  while(!feof(lfp)){
    fgets(line,200,lfp);
    if(feof(lfp)) break;
    if(strcmp(line,"\n")==0) break;
    nfiles++;
    if(nfiles>1000){
      fprintf(stderr,"Exceeded maximum number of input files (1000)\n");
      exit(1);
    }
    char* str=strtok(line," \t");
    fnames[nfiles-1]=(char*)malloc(sizeof(char)*200);
    strcpy(fnames[nfiles-1],str);
    str=strtok(NULL," \t");
    typs[nfiles-1]=atoi(str); // 0=SARTOPO CSV ,1=ALTIMETRY CSV, 2=BIDR
    if(typs[nfiles-1]==2){
    str=strtok(NULL," \t");
    times[nfiles-1]=atof(str); // TIME in seconds since J2000.0
    }
  }
  fclose(lfp);
  return(1);
}
int main(int argc, char* argv[]){



  if(argc!=4 && argc !=5){
    fprintf(stderr,"Usage: topomap_from_bidr cfgfile filelist outmap [heightonlyfile]\n");
    exit(0);
  }
  // start try
  try{
  int clidx=1;
  const char* cfgfile = argv[clidx++];
  string filelist = argv[clidx++];
  string outmap = argv[clidx++];
  string honlyfile;
  bool height_only=false;
  if(argc==5){
    height_only=true;
    honlyfile=argv[clidx++];
  }
  Uvar::setMode("no_unit_support");
  Config cfg("topomap_from_bidr",cfgfile);

  fprintf(stderr,"Config file %s read\n",cfgfile);
 
  // construct topomap input map
  TopoMap tmap;
  tmap.config(cfg,"NONE");

  fprintf(stderr,"TopoMap Object constructed\n");
  
 
  // read list of files
  int filec;
  char** filev;
  int* filetyp;
  double* ftimes;
  readfilelist(filelist,filec,filev,filetyp,ftimes);  
  for(int i=0;i<filec;i++){
    fprintf(stderr,"Processing file %s ...\n",filev[i]);
    tmap.addFile(filev[i],filetyp[i],ftimes[i]);
  } // end i loop


  // output map
  tmap.write(outmap);
  bool spmode=filetyp[0]==2 || filetyp[0]==3; // use a special mode that does not filter and write out a filenum backplane when using bidrs
  if(height_only) tmap.writeHeightsOnly(honlyfile,spmode);
  
 

  } // end try

  catch( ErrorMessage e){
    cerr << "Error: " << e.msg << endl;
  }
 }
