#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char* argv[]){

  int clidx=1;
  if(argc!=3){
    fprintf(stderr,"Usage: maker_bounds_str oldstring newstring\n");
    exit(1);
  }
  char* str1=argv[clidx++];
  char* str2=argv[clidx++];
  int len = strlen(str1);
  int len2 = strlen(str2);
  int startoff=0;
  for(int c=0;c<len;c++){
    if((str1[c]>='0' && str1[c]<='9') || str1[c]=='.' || str1[c]=='-'){
      if(startoff<len2){
	str1[c]=str2[startoff];
	startoff++;
      }
      else{
	str1[c]='0';
	startoff++;
      }
    }
    else if(startoff == 0){
      if(str1[c]==' ') str1[c]='@'; // replace spaces to alow top level script to work
      continue;
    }
    else{
      if(str1[c]==' ') str1[c]='@'; // replace spaces to alow top level script to work
      continue;
    }
    
  }
  printf("%s",str1);
  return(0);
}

