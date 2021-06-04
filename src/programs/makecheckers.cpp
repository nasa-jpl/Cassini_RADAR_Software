#include<stdio.h>
#include<stdlib.h>

int main(int argc, char* argv[]){

  if(argc!=6){
    fprintf(stderr,"Usage: makecheckers latbins lonbins latsquares lonsquares outfile\n");
    exit(1);
  }

  int cli=1;
  int latbins=atoi(argv[cli++]);
  int lonbins=atoi(argv[cli++]);
  int latsq=atoi(argv[cli++]);
  int lonsq=atoi(argv[cli++]);
  char* outfile=argv[cli++];
  FILE* ofp=fopen(outfile,"w");
  if(ofp==NULL){
    fprintf(stderr,"Cannot create outfile file %s\n",outfile);
    exit(1);
  }
  int lats=latbins/latsq;
  int lons=lonbins/lonsq;
  for (int i=0;i<latbins;i++){
    for (int j=0;j<lonbins;j++){
      float val=0.0;
      /** if within a white square then set value to 100 **/
      bool oddlat=(bool)((i/lats)%2);
      bool oddlon=(bool)((j/lons)%2);
      if((oddlon && oddlat)  || (!oddlat && !oddlon)){
	val =100;
      }

      if(fwrite(&val,sizeof(float),1,ofp)!=1){
	fprintf(stderr,"Cannot write to outfile file %s\n",outfile);
	exit(1);
      }
    }
  }
  fclose(ofp);
}
