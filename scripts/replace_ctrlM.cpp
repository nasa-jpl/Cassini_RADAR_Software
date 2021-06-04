#include<stdio.h>
// Routine to remove \r and replace with \n throughout an ASCII file
int main(){
  while(1){
    
    char c;
    c=fgetc(stdin);
    if(feof(stdin)) break;
    if(c=='\r') c='\n';
    printf("%c",c);
  }
  return(0);
}
