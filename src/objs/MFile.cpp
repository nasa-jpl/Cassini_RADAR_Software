#include <sstream>
#include "MFile.h"
#include "Error.h"
#include "Frame.h"
#include "Utils.h"

// Methods for MFile

MFile::MFile()
  : fp(NULL)
{
};

MFile::~MFile()
{
  close();
}

void
MFile::close(){
  if(fp){
    fclose(fp);
    fp=NULL;
  }
}

void 
MFile::flush()
{
  dieOnUnopenedFile("flush");
  fflush(fp);
}

void
MFile::open(const char* s)
{
  close();
  fp=fopen(s,"w");
  if(!fp){
    ostringstream os;
    os << "MFile cannot create file " << s;
    ErrorMessage e(toStr(os));
    e.throwMe();
  }
}

void
MFile::comment(const string& s)
{
  dieOnUnopenedFile("comment");
  fprintf(fp,"%% %s\n",s.c_str()); 
}

void 
MFile::set(const string & left, float v, int d1, int d2)
{
  leftHandSide(left,d1,d2);
  fprintf(fp,"= %15.10g;\n",v);
}

void 
MFile::set(const string & left, float* v, int len, int d1, int d2)
{
  leftHandSide(left,d1,d2);
  char c;
  if(d1==-1){
    c=';';
  }
  else c = ',';
  fprintf(fp,"= [");
  for(int i=0;i<len;i++){
    fprintf(fp,"%15.10g",v[i]);
    if(i<len-1) fprintf(fp,"%c",c);
  }
  fprintf(fp,"];\n");
}

void 
MFile::set(const string & left, complex<float>* v, int len, int d1, 
	   int d2)
{
  leftHandSide(left,d1,d2);
  char c;
  if(d1==-1){
    c=';';
  }
  else c = ',';
  fprintf(fp,"= [");
  for(int i=0;i<len;i++){
    fprintf(fp,"%15.10g+%15.10gi",real(v[i]),imag(v[i]));
    if(i<len-1) fprintf(fp,"%c",c);
  }
  fprintf(fp,"];\n");
}

void
MFile::leftHandSide(const string& left, int d1, int d2){
  dieOnUnopenedFile("leftHandSide");

  fprintf(fp,"%s",left.c_str());

  if(d1 || d2) fprintf(fp,"(");

  if(d1>0) fprintf(fp,"%d",d1);
  else if(d1==-1) fprintf(fp,":");

  if(d1 && d2) fprintf(fp,",");

  if(d2>0) fprintf(fp,"%d",d2);
  else if(d2==-1) fprintf(fp,":");
  
  if(d1 || d2) fprintf(fp,")");

}

void MFile::dieOnUnopenedFile(const string& funcname){
  if(!fp){
    ErrorMessage e("MFILE::"+funcname+" File unopened");
    e.throwMe();
  }
}
