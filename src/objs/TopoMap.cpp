#include <stdlib.h>
#include <string.h>
#include"TopoMap.h"
#include"Utils.h"
#include"SimpleArray.h"
#include"Constants.h"
#include"BIDRFile.h"
#include<iostream>

using std::cerr;
using std::endl;

TopoMap::TopoMap()
  : polar_proj(false), wraparound(false), num(NULL)
{}

TopoMap::~TopoMap()
{
  if(num!=NULL){
    free_array(num,2,nLats,nLons);
    free_array(height,2,nLats,nLons);
  }
}


// routine to initialize paarmeters called by read and config
int 

TopoMap::init(){
  latStep=(maxLat-minLat)/(nLats-1); 
  lonStep=(maxLon-minLon)/nLons;

  // set wraparound indicator
  if(maxLon-minLon>2*pi-0.01*lonStep) wraparound=true;
  else wraparound=false;

  //sanity check
  if(maxLon-minLon>2*pi+0.01*lonStep){
    fprintf(stderr,"TopoMap::config Longitude range spans more than 360 degrees\n");
    exit(1);
  }

  // special polar processing North mode
  if(maxLat>89.5*degtorad && minLat>30*degtorad && maxLon-minLon>=359.9*degtorad){
    polar_proj=true;
    if(nLons!=nLats){
      fprintf(stderr,"Trying to use Polar Projection need Nlons=Nlats.\n");
	exit(1);
    }
    // set lat step and lon step to 2*(90-minlat)/nlats
    latStep=2*(maxLat-minLat)/(nLats-1); // actually vertical step
    lonStep=2*(maxLat-minLat)/(nLats-1); // actually horizontal step
      minPLat=minLat-maxLat;
      maxPLat=maxLat-minLat;
      minPLon=minLat-maxLat;
      maxPLon=maxLat-minLat;
  }
  // special polar processing mode South pole
  else if(maxLat<-30*degtorad && minLat<3-89.5*degtorad && maxLon-minLon>=359.9*degtorad){
    polar_proj=true;
    if(nLons!=nLats){
      fprintf(stderr,"Trying to use Polar Projection need Nlons=Nlats.\n");
	exit(1);
    }
    // set lat step and lon step to 2*(90-minlat)/nlats
    latStep=2*(maxLat-minLat)/(nLats-1); // actually vertical step
    lonStep=2*(maxLat-minLat)/(nLats-1); // actually horizontal step
      minPLat=minLat-maxLat;
      maxPLat=maxLat-minLat;
      minPLon=minLat-maxLat;
      maxPLon=maxLat-minLat;
  }
  else{
    minPLat=minLat;
    maxPLat=maxLat;
    minPLon=minLon;
    maxPLon=maxLon;
  }
  return(1);
}

int
TopoMap::config(Config& cfg, string file){

  if(file!="None" && file!="none" && file!="NONE"){
    read(file);
  }
  else{
    nLons=cfg.getInt("TOPOMAP_NUM_LONS");
    nLats=cfg.getInt("TOPOMAP_NUM_LATS");
    dhthresh=(float)cfg.getDouble("TOPOMAP_DH_THRESH");
    maxNeighborDist=(float)cfg.getDouble("TOPOMAP_MAX_NEIGHBOR_DIST");
    gridLonStep=(float)cfg.getDouble("TOPOMAP_GRID_LON_STEP");
    gridLatStep=(float)cfg.getDouble("TOPOMAP_GRID_LAT_STEP");
    gridThickness=cfg.getInt("TOPOMAP_GRID_THICKNESS");
    minLat=(float)cfg.getDouble("TOPOMAP_GRID_MIN_LAT");
    maxLat=(float)cfg.getDouble("TOPOMAP_GRID_MAX_LAT");
    minLon=(float)cfg.getDouble("TOPOMAP_GRID_MIN_LON");
    maxLon=(float)cfg.getDouble("TOPOMAP_GRID_MAX_LON");
    allocateArrays();  
    init(); 
  }

  return(1);
}

// Lon input is East Longitude !!!
float TopoMap::getHeightInKm(double lat, double lon, bool& inbounds){
       if(!polar_proj){
	  while(lon< minLon) lon+=2*pi;
	  while(lon>= maxLon) lon-=2*pi;
	}
	else projectPolar(lat,lon);
	int i=int((lat-minPLat)/latStep);
	int j=int((lon-minPLon)/lonStep);
        // North pole HACK (will allow values exactly at LatMax to be accessed)
        // Fragile will also access up to latStep beyond upper latitude 
	// boundary
        if(i==nLats) i=nLats-1;
        if(i<0 || j<0 || i>nLats || j>=nLons || num[i][j]==0){
          inbounds=false;
	  return(0.0);          
	}

        inbounds=true;
	float h=height[i][j]/num[i][j];
        return(h/1000.0); // returns values in km
}


// Lon input is East Longitude !!!
float TopoMap::getInterpolatedHeightInKm(double lat, double lon, bool& inbounds){
       if(!polar_proj){
	  while(lon< minLon) lon+=2*pi;
	  while(lon>= maxLon) lon-=2*pi;
	}
       else{
	 // projectPolar(lat,lon);
	 fprintf(stderr,"TopoMap:getInterpolatedHeightInKm files for polar projections!\n");
	 exit(1);
       }
       int i=int(floor((lat-minPLat)/latStep));
       int j=int(floor((lon-minPLon)/lonStep));
       int i2=i+1;
       int j2=j+1;
       float ifloat=(lat-minPLat)/latStep;
       float jfloat=(lon-minPLon)/lonStep;
       float c11=(i2-ifloat)*(j2-jfloat);
       float c12=(i2-ifloat)*(jfloat-j);
       float c21=(ifloat-i)*(j2-jfloat);
       float c22=(ifloat-i)*(jfloat-j);


       if(i<0 || j2<0 || i>=nLats-1 || j>nLons || (j>=nLons-1 && !wraparound) || (j==-1 && !wraparound)){
          inbounds=false;
	  return(0.0);          
       }
        
       if(j==nLons-1 && wraparound) j2=0;
       if(j==nLons && wraparound){
	 j2=1;
	 j=0;
       }
       if(j==-1 && wraparound) j=nLons-1;

       if(num[i][j]==0 || num[i][j2]==0 || num[i2][j]==0 || num[i2][j2]==0){
          inbounds=false;
	  return(0.0);          
        }
        inbounds=true;
	float h11=height[i][j]/num[i][j];
	float h12=height[i][j2]/num[i][j2];
	float h21=height[i2][j]/num[i2][j];
	float h22=height[i2][j2]/num[i2][j2];

        float h=c11*h11+c12*h12+c21*h21+c22*h22;
        return(h/1000.0); // returns values in km
}


float TopoMap::getNearbyHeightInKm(double lat, double lon, bool& inbounds, float radbound_radians){
       if(!polar_proj){
	  while(lon< minLon) lon+=2*pi;
	  while(lon>= maxLon) lon-=2*pi;
       }
       else{
	 inbounds=false;
	 fprintf(stderr,"getNearbyHeightInKm doesn't work for Polar projections\n");
	 exit(1);
       }
       int radbound=int(radbound_radians/latStep +0.5);
       int i=int((lat-minPLat)/latStep);
       int j=int((lon-minPLon)/lonStep);
       // North pole HACK (will allow values exactly at LatMax to be accessed)
       // Fragile will also access up to latStep beyond upper latitude 
       // boundary
       if(i==nLats) i=nLats-1;

       if(i<0 || j<0 || i>nLats || j>=nLons){
	 inbounds=false;
	  return(0.0);          
       }

       int nfound=0;
       int rad=1;
       int lonrad=(int)((rad/cos(lat))*(latStep/lonStep)+0.5);
       	
       if(num[i][j]>0) nfound=1;
       float h=height[i][j]/num[i][j];
       int totnlons=(int)(2*pi/lonStep +0.5);
       if(lonrad>totnlons/2) lonrad=totnlons/2;
       while(nfound==0 && rad<radbound){
	 h=0;
	 for(int i2=i-rad;i2<i+rad;i2++){
	   for(int j2=j-lonrad;j2<j+lonrad;j2++){
	     int j3=j2;
	     if(j2<0) j3=j2+totnlons;
             if(j2>=nLons) j3=j2-totnlons;
	     if(j3<0 || j3>=nLons || i2>=nLats || i2<0) continue;
	     if(num[i2][j3]==0) continue;
	     h+=height[i2][j3]/num[i2][j3];
	     nfound++;
	   }
	 }
         h=h/nfound;
         rad++;
	 lonrad=(int)(rad/cos(lat)+0.5);
       }
       if(rad>=radbound){
	 inbounds=false;
	 return(0);
       }
       inbounds=true;

       return(h/1000.0); // returns values in km
}


int TopoMap::projectPolar(double& lat,double& lon){

  // North pole case
  if(minLat > 0 ){
    double latold=lat;
    double lonold=lon;
    double r=pi/2-latold;
    double ang=lonold-3*pi/2;
    // Makes prime meridian extend toward top of image with clockwise = west
    while(ang<0)ang+=2*pi;
    while(ang>2*pi)ang-=2*pi;
    lon=r*cos(ang);
    lat=r*sin(ang);
    return(1);
  }
  // South pole case
  else if(maxLat < 0){
    double latold=lat;
    double lonold=lon;
    double r=latold+pi/2;
    double ang=pi/2-lonold;
    // Makes prime meridian extend toward top of image with counterclockwise = west
    while(ang<0)ang+=2*pi;
    while(ang>2*pi)ang-=2*pi;
    lon=r*cos(ang);
    lat=r*sin(ang);
    return(1);  
  }
  
  else{
    fprintf(stderr,"Bad polar projection; neither North nor South ??\n");
    exit(1);
  }
  
}


// takes radians 
int TopoMap::unprojectPolar(double& lat,double& lon){
  // North pole case
  if(minLat>0){
    double latold=lat;
    double lonold=lon;
    double  r=sqrt(latold*latold+lonold*lonold);
    lat=pi/2-r;
    if(r==0){
      lon=0;
    }
    else{
      double ang=atan2(latold/r,lonold/r);
      lon=ang+3*pi/2;
    }  
  }

  // South pole case
  if(maxLat<0){
    double latold=lat;
    double lonold=lon;
    double  r=sqrt(latold*latold+lonold*lonold);
    lat=r-pi/2;
    if(r==0){
      lon=0;
    }
    else{
      double ang=atan2(latold/r,lonold/r);
      lonold=pi/2-ang;
    }  
  }
  else{
    fprintf(stderr,"Bad polar projection; neither North nor South ??\n");
    exit(1);
  }
  return(1);
}



int TopoMap::read(string file){
  FILE* fp=fopen(file,"r");
  if(fread(&nLats,sizeof(int),1,fp)!=1 ||
     fread(&nLons,sizeof(int),1,fp)!=1 ||
     fread(&dhthresh,sizeof(float),1,fp)!=1 ||
     fread(&minLat,sizeof(float),1,fp)!=1 ||
     fread(&maxLat,sizeof(float),1,fp)!=1 ||
     fread(&latStep,sizeof(float),1,fp)!=1 ||
     fread(&minLon,sizeof(float),1,fp)!=1 ||
     fread(&maxLon,sizeof(float),1,fp)!=1 ||
     fread(&lonStep,sizeof(float),1,fp)!=1 ||
     fread(&polar_proj,sizeof(int),1,fp)!=1 ||
     fread(&maxNeighborDist,sizeof(float),1,fp)!=1 ||
     fread(&gridLonStep,sizeof(float),1,fp)!=1 ||
     fread(&gridLatStep,sizeof(float),1,fp)!=1 ||
     fread(&gridThickness,sizeof(int),1,fp)!=1){
    cerr << "TopoMap::read Error reading from " << file << endl;
    exit(1);
  }
  allocateArrays();
  if( ! read_array(fp,height,sizeof(float),2,nLats,nLons) ||
      ! read_array(fp,num,sizeof(int),2,nLats,nLons)){
    cerr << "TopoMap::read Error reading from " << file << endl;
    exit(1);

  }
  fclose(fp);

  init();
  return(1);
}

int TopoMap::allocateArrays(){
  height=(float**)make_array(sizeof(float),2,nLats,nLons);
  num=(int**)make_array(sizeof(int),2,nLats,nLons);
  if(!height || !num){
    fprintf(stderr,"TopoMap: error cannot allocate Arrays\n");
    exit(1);
  }
  for(int i=0;i<nLats;i++){
    for(int j=0;j<nLons;j++){
      height[i][j]=0;
      num[i][j]=0;
    }
  }
  return(1);
}


int TopoMap::write(string file){
  FILE* fp=fopen(file,"w");
  if(fwrite(&nLats,sizeof(int),1,fp)!=1 ||
     fwrite(&nLons,sizeof(int),1,fp)!=1 ||
     fwrite(&dhthresh,sizeof(float),1,fp)!=1 ||
     fwrite(&minLat,sizeof(float),1,fp)!=1 ||
     fwrite(&maxLat,sizeof(float),1,fp)!=1 ||
     fwrite(&latStep,sizeof(float),1,fp)!=1 ||
     fwrite(&minLon,sizeof(float),1,fp)!=1 ||
     fwrite(&maxLon,sizeof(float),1,fp)!=1 ||
     fwrite(&lonStep,sizeof(float),1,fp)!=1 ||
     fwrite(&polar_proj,sizeof(int),1,fp)!=1 ||
     fwrite(&maxNeighborDist,sizeof(float),1,fp)!=1 ||
     fwrite(&gridLonStep,sizeof(float),1,fp)!=1 ||
     fwrite(&gridLatStep,sizeof(float),1,fp)!=1 ||
     fwrite(&gridThickness,sizeof(int),1,fp)!=1){
    cerr << "TopoMap::write Error writing to " << file << endl;
    exit(1);
  }
  if( ! write_array(fp,height,sizeof(float),2,nLats,nLons) ||
      ! write_array(fp,num,sizeof(int),2,nLats,nLons)){
    cerr << "TopoMap::write Error writing to " << file << endl;
    exit(1);
  }
  fclose(fp);
  return(1);
}

int
TopoMap::writeHeightsOnly(string file, bool spmode = false){
  // if spmode it true then it outputs another file with num values in it and does no filtering on the data
  FILE* fp=fopen(file,"w");

  if( spmode ){
    FILE* fp2 = fopen(file+".num","w");
    for(int i=nLats-1;i>=0;i--){
      if( fwrite(&height[i][0],sizeof(float),nLons,fp) != nLons){
	cerr << "TopoMap::writeHeightsOnly Error writing to " << file << endl;
	exit(1);
      }
      if( fwrite(&num[i][0],sizeof(int),nLons,fp2) !=nLons){
	cerr << "TopoMap::writeHeightsOnly Error writing to " << file+".num" << endl;
	exit(1);
      }
    }
    fclose(fp);    
    fclose(fp2);    
    return(1);
  }
  float** h=(float**)make_array(sizeof(float),2,nLats,nLons);

  float** hfilted=(float**)make_array(sizeof(float),2,nLats,nLons);
  if(!h || !hfilted){
    fprintf(stderr,"TopoMap::Error making temporary array\n");
    exit(1);
  }
  for(int i=0;i<nLats;i++){
    for(int j=0;j<nLons;j++){
      if(num[i][j]!=0){
	h[nLats-i-1][j]=height[i][j]/num[i][j];
      }
      else h[nLats-i-1][j]=-1000000;
    }
  }

  // nearest neighbor filter  (make this a separate routine later)
  for(int i=0;i<nLats;i++){
    for(int j=0;j<nLons;j++){
      if(h[i][j]<-999999){
	float mind=500;
        float bestfilenum=100000;
        int mini2=0;
        int minj2=0;
        int imaxdist=(int)ceil(maxNeighborDist/(2575*latStep));
        int rad=2*imaxdist;
        float trad=2575;
        double lat1=(minPLat+latStep*i);
        double lon1=(minPLon+lonStep*j);
	double r=lat1*lat1+lon1*lon1;
	double jacob=sin(r)/r;

	int jmaxdist=(int)ceil(maxNeighborDist/(2575*latStep*cos(lat1)));
        if(polar_proj){
	  jmaxdist=(int)ceil(maxNeighborDist/(2575*latStep*jacob));
	  imaxdist=jmaxdist;
          rad=2*imaxdist;
	}
        int jrad=2*jmaxdist;
        if(jrad>nLons/2 +1) jrad=nLons/2 +1;
        double x1,x2,y1,y2,z1,z2;
	if(polar_proj){
	  x1=trad*lon1*jacob;
	  y1=trad*lat1*jacob;
	  z1=0;
	}
        else{
	  x1=trad*cos(lat1)*sin(lon1);
	  y1=trad*cos(lat1)*cos(lon1);
	  z1=trad*sin(lat1);
	}
	for(int i2=i-rad;i2<=i+rad;i2++){
	for(int j2=j-jrad;j2<=j+jrad;j2++){
          int j3=j2;
 	  if(i2<0 || i2>=nLats) continue;
          if(j2<0) j3=j2+nLons;
          else if(j2>=nLons) j3=j2-nLons;
	  if(h[i2][j3]<-999999) continue;
 
          double lat2=(minPLat+latStep*i2);

          double lon2=(minPLon+lonStep*j3);


          if(polar_proj){
	    x2=trad*lon2*jacob;
	    y2=trad*lat2*jacob;
	    z2=0;
	  }
          else{
	    x2=trad*cos(lat2)*sin(lon2);
	    y2=trad*cos(lat2)*cos(lon2);
	    z2=trad*sin(lat2);
	  }

	  double dz=z2-z1;
          double dy=y2-y1;
          double dx=x2-x1;
          double d=sqrt(dx*dx+dy*dy+dz*dz);
          // this assures that we use the preferred image (lowest number)
          if(d<mind) 
	    { 
	      mind=d;
	      mini2=i2;
	      minj2=j3;
	    }
	} // end j2 loop
  
	} // end i2 loop
	if(mind<maxNeighborDist) hfilted[i][j]=h[mini2][minj2];
	else hfilted[i][j]=-1000000;
      } // h<-999999 case
      else hfilted[i][j]=h[i][j];
      } // end j loop
    } // end i loop

  // setGridLines(hfilted);
  if( ! write_array(fp,hfilted,sizeof(float),2,nLats,nLons)){
    cerr << "TopoMap::writeHeightsOnly Error writing to " << file << endl;
    exit(1);
  }
  fclose(fp);
  free_array(h,2,nLats,nLons);
  free_array(hfilted,2,nLats,nLons);
  return(1);
}


/*** Incomplete version
int
TopoMap::writeHeightsOnly(string file){
  FILE* fp=fopen(file,"w");
  if( ! write_array(fp,height,sizeof(float),2,nLats,nLons)){
    cerr << "TopoMap::writeHeightsOnly Error writing to " << file << endl;
    exit(1);
  }
  fclose(fp);
  return(1);
}
***/
int
TopoMap::setGridLines(float** h){
  if(gridThickness==0) return(1);
  int nlon_lines=(int)ceil((maxLon-minLon)/gridLonStep) -1;
  int nlat_lines=(int)ceil((maxLat-minLat)/gridLonStep) -1;
  int lrad = gridThickness/2;
  
  // make longitude lines
  for(int ilonlines=1;ilonlines<=nlon_lines;ilonlines++){
    float lon=minLon+ilonlines*gridLonStep;
    int j=int((lon-minLon)/lonStep);
    for(int j2=j-lrad;j2<=j+lrad;j2++){
      for(int i=0;i<nLats;i++){
	h[i][j2]=-2;
      }
    }
  }
  // make latititude_lines
  for(int ilatlines=1;ilatlines<=nlat_lines;ilatlines++){
    float lat=minLat+ilatlines*gridLatStep;
    int i=int((lat-minLat)/latStep);
    for(int i2=i-lrad;i2<=i+lrad;i2++){
      for(int j=0;j<nLons;j++){
	h[i2][j]=-2;
      }
    }
  }
  return(1);
}
int
TopoMap::addFile(char* fname, int typ, double tca){
  static int fileno = 1;
  if(typ==0){   // SARTopo CSV case
    char line[200];
    FILE* ifp=fopen(fname,"r");
    if(ifp==NULL){
	fprintf(stderr,"Error reading file %s\n",fname);
        fprintf(stderr,"Skipping it\n");
    }
    else{   
      while(!feof(ifp)){
	fgets(line,200,ifp);
	if(feof(ifp)) break;
	char* str=strtok(line,",");//wlon
	lon=-atof(str)*degtorad;
        str=strtok(NULL,",");//lat
        lat=atof(str)*degtorad;
        str=strtok(NULL,",");//inc
        str=strtok(NULL,",");//croswid
        double wid1=atof(str);
        str=strtok(NULL,",");//alongwid
        double wid2=atof(str);
        str=strtok(NULL,",");//h
        double h=atof(str);
        str=strtok(NULL,",");//herr
        str=strtok(NULL,",");//flag
        int flag=atoi(str);
	if(lat==90) lat=89.99999;

        //if(flag) continue;
        if(lat<minLat || lat > maxLat) continue;
        if(!polar_proj){
	  while(lon< minLon) lon+=2*pi;
	  while(lon>= maxLon) lon-=2*pi;
	}
	else projectPolar(lat,lon);
	int i=int((lat-minPLat)/latStep);
	int j=int((lon-minPLon)/lonStep);
	if(i<0 || j<0 || i>=nLats || j>=nLons) continue; // continue on out of bounds
	
	height[i][j]=h;
        num[i][j]=1;
      }
      fclose(ifp);        
    }
  }
  else if(typ==1){   // Altimeter CSV case
    char line[200];
    FILE* ifp=fopen(fname,"r");
      if(ifp==NULL){
	fprintf(stderr,"Error reading file %s\n",fname);
	fprintf(stderr,"Skipping it\n");
      }
      else{
	while(!feof(ifp)){
	  fgets(line,200,ifp);
	  if(feof(ifp)) break;
	  char* str=strtok(line,",");// sab
	  str=strtok(NULL,","); // t
	  str=strtok(NULL,","); // range
          double range=atof(str);
	  str=strtok(NULL,","); // wlon
 	  lon=-atof(str)*degtorad;
	  str=strtok(NULL,","); // lat
	  lat=atof(str)*degtorad;
	  str=strtok(NULL,",");//h_thres
	  str=strtok(NULL,",");//h_mle
	  str=strtok(NULL,",");//h_cent
	  str=strtok(NULL,",");//h_cent_corr
	  double h=atof(str)-2575000;
	  str=strtok(NULL,",");//h_thre_corr
	  str=strtok(NULL,",");//h_span
          double h_span=atof(str);
	  str=strtok(NULL,",");//h_skew
	  str=strtok(NULL,",");//inc
	  double inc=atof(str);
	  if(lat==90) lat=89.99999;

	  if(inc>0.4 || h_span>300 || range>20000000 || lat<minLat || lat>maxLat) continue;
          if(!polar_proj){
	    while(lon< minLon) lon+=2*pi;
	    while(lon>= maxLon) lon-=2*pi;
	  }
	  else projectPolar(lat,lon);
	  int i=int((lat-minPLat)/latStep);
	  int j=int((lon-minPLon)/lonStep);
	  if(i<0 || j<0 || i>=nLats || j>=nLons) continue; // continue on out of bounds

	  height[i][j]=h;
	  num[i][j]=1;
	}
	fclose(ifp);
      }
  }
  else if(typ==2 || typ==3){   // BIDR case   2= IAU pole  3= new pole
    // Read in BIDR file header
    BIDRFile bf(fname,"r");
    bf.readHeader();
    bf.rewind();
    PDSLabel label(bf);
  
 
    // if file is not PDS
    if(label.label() == EMPTY_STRING){
      fprintf(stderr,"%s is not a bidr file\n",fname);
      exit(1);
    }


    // skip to correct position in file
    int nrows=label.fileRecords();
    int headrecs=label.labelRecords();
    nrows=nrows-headrecs;
    int ncols=label.recordLength()/sizeof(float);
 
    if(typ==2) bf.fixPole(tca);  // IAU pole case  
    // reopen files as standard FILE* files for fast read/write
    FILE* ifp=fopen(fname,"r");
 

    
    // skip headers in BIDR 
    fseek(ifp,headrecs*ncols*sizeof(float),SEEK_SET);



    // allocate buffers
    float* s0=(float*)malloc(sizeof(float)*ncols);

    // Accumulate map
    for(int j=0;j<nrows;j++){
      if(j%100==0)fprintf(stderr,"%s: %d out %d rows mapped.\n",fname,j,nrows);
      // read row of data from each input file
      if((int)fread(&s0[0],sizeof(float),ncols,ifp)!=ncols){
	ErrorMessage e("Error Reading Row "+toStr(j+1)+" from file "
		       +fname);
	e.throwMe();
      }
      for(int i=0;i<ncols;i++){
	// compute lat and lon
	double slat,slon;
	
	// Calculate the oblique cylindrical latitude and longitude in radians at
	// the selected point.
	
	double lat = bf.OCLatFromGrid(i+1);
	double lon = bf.OCLonFromGrid(j+1);
	
	// Convert the latitude and longitude to standard coordinates.
	slat = bf.standardLatInRad(lon, lat);
	slon = bf.standardLonInRad(lon, lat);
      
	// convert to positive west
	slon = 2*pi - slon;

	if(s0[i]<-1000000) continue;
	float val=s0[i];

	add(slat,slon,val,0.0001,false,fileno); // only keep first entry do not average!
      } // end i loop
    } // end j loop
  }
  // Triaxial body file case
  // Format is one line with four real values
  // x_rad y_rad (prime meridian radius) z_rad (polar_radius)
  // all radii in meters.
  // Optional comment lines (starting with #) are ignored
  else if (typ==4){
        char line[200];
        double rad=0,xrad=0,yrad=0,zrad=0;
	FILE* ifp=fopen(fname,"r");
	if(ifp==NULL){
	  fprintf(stderr,"Error reading file %s\n",fname);
	  fprintf(stderr,"Skipping it\n");
	}
	else{    
	  while(!feof(ifp)){
	    fgets(line,200,ifp);
            if(line[0]=='#')
	      continue; // skip comments
	    if(feof(ifp)) break;
	    char* str=strtok(line,", \t");
            if(str[0]=='#')
	      continue; // skip comments
	    xrad=atof(str)*1000; // heights are stored in m
                                 // radii in triaxial fiels are in km
	    str=strtok(NULL,", \t");
	    yrad=atof(str)*1000;
	    str=strtok(NULL,", \t");
	    zrad=atof(str)*1000;
	    rad=sqrt(xrad*xrad+yrad*yrad+zrad*zrad)/sqrt(3);
            break;
 	  }
	  if(xrad==0 || yrad==0 || zrad==0 || rad==0){
	    fprintf(stderr,"Bad file read for triaxial file %s\n");
	    fprintf(stderr,"Skipping ..\n");
	  }
	  else{
	    for(int i=0;i<nLats;i++){
	      for(int j=0;j<nLons;j++){
		double lat=minPLat+i*latStep+0.5*latStep;
		double lon=minPLon+j*lonStep+0.5*lonStep;
		if(polar_proj) unprojectPolar(lat,lon);
		double theta=pi/2-lat;
		// computed from
		//  (x/xrad)^2+(y/yrad)^2+(z/zrad)^2=1
		//  z/r=cos(theta)
		//  x/r=cos(lon)sin(theta)
		//  y/r=sin(lon)sin(theta)
		//  h=r-rad;
		double zoverr=cos(theta);
		double yoverr=sin(theta)*sin(lon);
		double xoverr=sin(theta)*cos(lon);
		zoverr=zoverr/zrad;
		yoverr=yoverr/yrad;
		xoverr=xoverr/xrad;
		zoverr*=zoverr;
		yoverr*=yoverr;
		xoverr*=xoverr;
		double r=sqrt(1/(xoverr+yoverr+zoverr));
		height[i][j]=(r-rad);
		num[i][j]=1; 
	      } // end nlons loop	
	    } // end nlats loop
	  } // triaxial good file case
	} // triaxial file found case
  } // end triaxial case
  fileno++;
  return(1);
}

int
TopoMap::add(double slat, double slon, float h, float dh, bool average=true, int fileno=1){

  // invalid data trap  turned off for now
  // if(dh>dhthresh || h<-2 || h>2) return(0);

  // convert slat slon to proper format and compute indices
  lat=slat;
  lon=-slon; // Map indices used positive East Lon so that when North is up
             // East is right;

  if(!polar_proj){
    while(lon< minLon) lon+=2*pi;
    while(lon>= maxLon) lon-=2*pi;
  }

  // North polar projection case
  if(polar_proj & minLat > 0){
    if(lat<minLat) return(0);
    projectPolar(lat,lon);
  }

  // South polar projection case
  if(polar_proj & maxLat <0 ){
    if(lat> maxLat) return(0);
    projectPolar(lat,lon);
  }

  int i=int((lat-minPLat)/latStep);
  int j=int((lon-minPLon)/lonStep);


  // out of bounds trap
  
  if(i<0 || i>nLats) return(0);
  if(i==nLats) i=nLats-1;
  if(j<0 || j>=nLons) return(0); 

  if(average){
  height[i][j]+=h;
  num[i][j]++;
  }
  else if(num[i][j]==0){
    num[i][j]=fileno;
    height[i][j]=h;
  }
  return(1);  
}
