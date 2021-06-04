#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define MAXCOUNT 1000000 
int main(int argc, char* argv[]){
  const char* command=argv[0];
  bool unzipping=false;
  if(strstr(command,"zero_run_float_unzip")!=NULL){
    unzipping=true;
  }
  char* rawfile=NULL;
  char* zipfile=NULL;

  FILE* ifp;
  FILE* ofp;
  float val=0;
  int count=0;
  int clidx=1;
  float* zeroes=(float*)calloc(MAXCOUNT,sizeof(float));
  if(unzipping){
    if(argc!=3){
      fprintf(stderr,"Usage: zero_run_float_unzip zipfile uncompressedfile\n");
      exit(1);
    }
    else{
      zipfile=argv[clidx++];
      rawfile=argv[clidx++];
      ifp=fopen(zipfile,"r");
      ofp=fopen(rawfile,"w");
      while(!feof(ifp)){
	fread(&val,sizeof(float),1,ifp);
        if(feof(ifp)) break;
        if(val==0){
	  fread(&count,sizeof(int),1,ifp);
          while(count>0){
	    if(count<MAXCOUNT) fwrite(&zeroes[0],sizeof(float),count,ofp);
            else fwrite(&zeroes[0],sizeof(float),MAXCOUNT,ofp);
            count=count-MAXCOUNT;
	  }
	}
        else{
	  fwrite(&val,sizeof(float),1,ofp);
	}
      } // end while loop
    } // end of valid command line
  } // end of unzipping case

  // zipping case
  else{

    if(argc!=3){
      fprintf(stderr,"Usage: zero_run_float_zip rawfile zipfile\n");
      exit(1);
    }
    else{
      rawfile=argv[clidx++];
      zipfile=argv[clidx++];
      ifp=fopen(rawfile,"r");
      ofp=fopen(zipfile,"w");
      while(!feof(ifp)){
	fread(&val,sizeof(float),1,ifp);
        if(feof(ifp)) break;
	fwrite(&val,sizeof(float),1,ofp);
        count=0;
        while(val==0){
          count++;
	  fread(&val,sizeof(int),1,ifp);
          if(feof(ifp)) break;
	}
	if(count>0)fwrite(&count,sizeof(int),1,ofp);
	if(!feof(ifp) && count>0 ) fwrite(&val,sizeof(int),1,ofp);
      } // end while loop
    }
  } // end of zipping case
  fclose(ifp);
  fclose(ofp);
  return(0);
}
