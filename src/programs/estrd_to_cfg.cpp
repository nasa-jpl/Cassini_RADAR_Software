#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Utils.h"
#include<string.h>

using std::cerr;
using std::cout;
using std::endl;



int main(int argc, char* argv[])
{
  // parse command line
  int clidx=0;
  string command=argv[clidx++];
  if(argc!=6 && argc!=8){
    cerr << "Usage: " << command << " in_estrdfile out_cfgfile range_poly_order  dop_poly_order separate_beams? [start_time] [end_time]" 
	 << endl;
    exit(1);
  }
    
  string infile=argv[clidx++];
  string outfile=argv[clidx++];

  int rpoly_order=atoi(argv[clidx++]);
  int dpoly_order=atoi(argv[clidx++]);
  bool sep=(bool)atoi(argv[clidx++]);
  float start_time=-100000,end_time=100000;
  if(argc==8){
    start_time=atof(argv[clidx++]);
    end_time=atof(argv[clidx++]);
  }

  double range_mean=0, range_std=0, dop_mean=0, dop_std=0, time_mean=0, time_std=0;
  
  // assume 20000 is enough points for now
  int NMAX=20000;

  int numbeams=1;
  if(sep) numbeams=5;

  Dmat t("t",NMAX,numbeams);
  Dmat r("r",NMAX,numbeams);
  Dmat d("d",NMAX,numbeams);
  Dmat rp("p",rpoly_order+1,numbeams);
  Dmat dp("p",dpoly_order+1,numbeams);


  // open input file
  FILE* ifp=fopen(infile,"r");
  

  // loop through lines and accumulate stats
  // stats are NOT computed separately for different beams
  
  char line[200];
  int num[]={0,0,0,0,0};
  while(!feof(ifp)){
    fgets(line,200,ifp);
    if(feof(ifp)) break;

    // read time from closest approach in s
    char* str=strtok(line," \t");
    if(str==NULL) break;
    double t0=strtod(str,NULL); 

    // skip sab number 
    str=strtok(NULL," \t"); 

    // read beam number
    str=strtok(NULL," \t");
    int b=atoi(str)-1;
 
    // eliminate separation by beam if sep=false
    if(!sep) b=0;

    // skip center indices
    str=strtok(NULL," \t");
    str=strtok(NULL," \t");

    // read estimated range correction in km
    str=strtok(NULL," \t");  
    double r0=strtod(str,NULL); 


    // read estimated doppler correction in Hz
    str=strtok(NULL," \t");  
    double d0=strtod(str,NULL);
    if(argc!=8 || (t0>start_time && t0<end_time)){
      // fill arrays
      t(num[b],b)=t0;
      r(num[b],b)=r0;
      d(num[b],b)=d0;
      
      // accumulate stats
      time_mean+=t0; 
      range_mean+=r0; 
      dop_mean+=d0; 
      time_std+=t0*t0; 
      range_std+=r0*r0; 
      dop_std+=d0*d0; 
      num[b]++;
    }
  }

  fclose(ifp);
  // complete stat computation (using data from all 5 beams)
  int n=num[0]+num[1]+num[2]+num[3]+num[4];
  time_mean/=n;
  range_mean/=n;
  dop_mean/=n;
  time_std-=time_mean*time_mean*n;
  time_std=sqrt(time_std/(n-1));

  range_std-=range_mean*range_mean*n;
  range_std=sqrt(range_std/(n-1));

  dop_std-=dop_mean*dop_mean*n;
  dop_std=sqrt(dop_std/(n-1));
  

  //--------------------------
  // compute polynomial fit 
  //--------------------------
 
  // normalize data to mean=0 var=1
  for(int b=0;b<numbeams;b++){
    for(int i=0;i<num[b];i++){
      t(i,b)-=time_mean;
      t(i,b)/=time_std;
      r(i,b)-=range_mean;
      r(i,b)/=range_std;
      d(i,b)-=dop_mean;
      d(i,b)/=dop_std;
    }
  }

  for(int b=0;b<numbeams;b++){
    if(num[b]==0) continue;
    Dvec x=t.getPartialCol(b,0,num[b]-1);
    Dvec y=r.getPartialCol(b,0,num[b]-1);
    Dvec tmp("tmp",rpoly_order+1);
    polyfit(tmp,y,x);
    for(int j=0;j<rpoly_order+1;j++){
      rp(j,b)=tmp(j);
    }
    y=d.getPartialCol(b,0,num[b]-1);
    tmp.resize(dpoly_order+1);
    polyfit(tmp,y,x);
    for(int j=0;j<dpoly_order+1;j++){
      dp(j,b)=tmp(j);
    }
  }  

  // output config file
  FILE* ofp=fopen(outfile,"w");

  // Output means, standard deviations, orders  and SET_DOPPLER_CENTROID
  fprintf(ofp,"SET_DOPPLER_CENTROID 1\n");
  fprintf(ofp,"DOPPLER_CENTROID_POLYNOMIAL_ORDER %d\n",dpoly_order);
  fprintf(ofp,"RANGE_CENTROID_POLYNOMIAL_ORDER %d\n",rpoly_order);
  fprintf(ofp,"ESTRD_TIME_MEAN %6.10g s\n",time_mean);
  fprintf(ofp,"ESTRD_TIME_STD %6.10g s\n",time_std);
  fprintf(ofp,"DOPPLER_CENTROID_MEAN %6.10g Hz\n",dop_mean);
  fprintf(ofp,"DOPPLER_CENTROID_STD %6.10g Hz\n",dop_std);
  fprintf(ofp,"RANGE_CENTROID_MEAN %6.10g km\n",range_mean);
  fprintf(ofp,"RANGE_CENTROID_STD %6.10g km\n",range_std);


  // output doppler coeffs
  for(int b=0;b<5;b++){
    int k=b;
    if(!sep) k=0;
    if(num[k]==0) continue;

    for(int i=0;i<dpoly_order+1;i++){
      fprintf(ofp,"DOPPLER_CENTROID_BEAM_%d_POLYNOMIAL_COEFF_%d %3.20g\n",b+1,i,	     
	      dp(i,k));
    }
  }

  // output range coeffs 
  for(int b=0;b<5;b++){
    int k=b;
    if(!sep) k=0;
    if(num[k]==0) continue;

    for(int i=0;i<rpoly_order+1;i++){
      fprintf(ofp,"RANGE_CENTROID_BEAM_%d_POLYNOMIAL_COEFF_%d %3.20g\n",b+1,i,
	      rp(i,k));
    }
  }


  // write xmgr fitted data output to stdout
  // for now do -25 to +25 minutes
  // t, r_b1, d_b1, r_b2 d_b2 ... d_b5
  for(int b=0;b<numbeams;b++){

    for(int i=0;i<num[b];i++){ 
      double t0=t(i,b);
      cout << (t0*time_std+time_mean)/60.0 << " ";
      double r0=0, d0=0;
      if(num[b]>0){
	for(int j=rpoly_order;j>0;j--){
	  r0+=rp(j,b);
	  r0*=t0;
	}
	r0+=rp(0,b);
	
	for(int j=dpoly_order;j>0;j--){
	  d0+=dp(j,b);
	  d0*=t0;
	}
	d0+=dp(0,b);
	r0=r0*range_std+range_mean;
	d0=d0*dop_std+dop_mean;
      }
      else{
	r0=0;
        d0=0;
      }
      cout << r0 <<" " << r(i,b)*range_std+range_mean << " " << d0 << " " <<d(i,b)*dop_std+dop_mean << endl;
    }
    cout << "&" << endl;
  }
  return(0);
}
