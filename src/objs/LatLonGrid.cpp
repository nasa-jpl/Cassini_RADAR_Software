#include"LatLonGrid.h"
#include<math.h>
#include"Error.h"
#include"Constants.h"
#include"DebugInfo.h"



using std::endl;
using std::cerr;

LatLonGrid::LatLonGrid()
  :num_lats_(0),num_lons_(0)
{
}

LatLonGrid::LatLonGrid(int pixels_per_degree, double min_lat_deg, double max_lat_deg, double start_lon_deg, 
		       double lon_width_deg, bool lon_increasing, bool lon360, int lon_buffer_in_deg)
{
  // compute grid size

  if(min_lat_deg>=max_lat_deg || min_lat_deg<-90 || max_lat_deg>90){
    ErrorMessage e("Bad latitude bounds in LatLonGrid constructor.");
    e.throwMe();
  }
  if(pixels_per_degree<=0){
    ErrorMessage e("Bad pixelsperdegree in LatLonGrid constructor.");
    e.throwMe();
  }
  num_lats_=(int)((max_lat_deg-min_lat_deg)*pixels_per_degree+0.5);

  total_num_lons_=(int)(ceil(lon_width_deg*pixels_per_degree));

  int bufsz=lon_buffer_in_deg*pixels_per_degree;
  
  // make buffer size 45 degrees if pixels_per_degree is less than 256
  // helps with HISAR processing
  // but causes turn processing to bomb
  // if(pixels_per_degree<129) bufsz=45*pixels_per_degree; // HACK

  // for 16 or fewer pixels per degree use one big buffer
  if(pixels_per_degree <=16 || lon360){
    bufsz=360*pixels_per_degree;
  }
  // set buffer size optimially for small regions
  int smaller_bufs[7]={8,6,5,4,3,2,1}; // values are factors of 360
  int idx=-1;
  int done=0;
  while(bufsz > total_num_lons_ && !done && pixels_per_degree> 16 & !lon360){
    idx++;
    if(idx==7) done=1;
    else if (smaller_bufs[idx]*pixels_per_degree> total_num_lons_){
      bufsz=smaller_bufs[idx]*pixels_per_degree;
    }
    else done=1;
  }
  num_lons_=bufsz;
  num_abs_lons_=360*pixels_per_degree;

  // compute start of latitude range (edge of first pixel)
  min_lat_=min_lat_deg*degtorad; // convert to radians

  res_=degtorad/(double)pixels_per_degree; // compute resolution in radians

  // computing starting valid longitude range
  abs_start_valid_i_=(int)floor(start_lon_deg*pixels_per_degree);
  lon_inc_=lon_increasing;

  start_lon_=start_lon_deg*degtorad;
  setValidRange();
}

// sets abs_start_valid_i and abs_end_valid_i so that they are both greater
// than or equal to zero and one of them is less than num_abs_lons_
// Needs initial value for abs_start_valid_i_ and value for lon_inc_
void LatLonGrid::setValidRange()
{
  if(lon_inc_)
    abs_end_valid_i_=abs_start_valid_i_+num_lons_-1;
  else
    abs_end_valid_i_=abs_start_valid_i_-num_lons_+1;

  while(abs_end_valid_i_ < 0 || abs_start_valid_i_ <0){
    abs_end_valid_i_+=num_abs_lons_;
    abs_start_valid_i_+=num_abs_lons_;
  }
  while(abs_end_valid_i_ >= num_abs_lons_ && 
	abs_start_valid_i_ >= num_abs_lons_)
    {
      abs_end_valid_i_-=num_abs_lons_;
      abs_start_valid_i_-=num_abs_lons_;
    }
}

void LatLonGrid::getBoundary(const double lonlat_bounds[4], int ij_bounds[4])
{
  float lonmin, lonmax, latmin,latmax;
  int imin,imax,jmin,jmax;

  // unpack input bounds
  lonmin=lonlat_bounds[0];
  lonmax=lonlat_bounds[1];
  latmin=lonlat_bounds[2];
  latmax=lonlat_bounds[3];

  // quantize
  imin=(int)floor(lonmin/res_);
  imax=(int)ceil(lonmax/res_);
  jmin=(int)floor((latmin-min_lat_)/res_);
  jmax=(int)ceil((latmax-min_lat_)/res_);

  // check for bad bounds
  if(jmin>jmax){

    ErrorMessage e("LatLonGrid::getBoundary Bad lat bounds latmin>latmax.");
    e.throwMe();
  }

  if(imin>imax){
    ErrorMessage e("LatLonGrid::getBoundary Bad lon bounds lonmin>lonmax.");
    e.throwMe();
  }
  // make sure all indices are positive
  while(imin<0){
    imin+=num_abs_lons_;
    imax+=num_abs_lons_;
  }
  if(jmin < 0 || jmax > num_lats_-1){
    /*** This was core dumping ......
    ErrorMessage e("LatLonGrid::getBoundary latitude out of bounds.");
    e.throwMe();
    ***/
    // And I decided a warning was preferable.
    if(DebugInfo::allWarnings){
      cerr << "Warning LatLonGrid::getBoundary latitude out of bounds." << endl;
      cerr << "Restricting LatLon Boundaries of burst to grid bounadries." << endl;
      cerr << "You will be warned again if this impacts performance." << endl;
    }
    if(jmin < 0) jmin=0;
    if(jmax > num_lats_-1) jmax=num_lats_-1;
  }

  // pack quantized bounds
  ij_bounds[0]=imin;
  ij_bounds[1]=imax;
  ij_bounds[2]=jmin;
  ij_bounds[3]=jmax;
}
 
int LatLonGrid::relLonIndex(int abs_i){
  return(abs_i%num_lons_);
}



double LatLonGrid::latInRad(int j)
{
  double lat=min_lat_+(j+0.5)*res_;
  return(lat);
}

double LatLonGrid::lonInRad(int abs_i)
{
  double lon=(abs_i+0.5)*res_;
  return(lon);
}

int LatLonGrid::absStartValidIndex(){
  return(abs_start_valid_i_);
}

int LatLonGrid::absEndValidIndex(){
  return(abs_end_valid_i_);
}

int LatLonGrid::absIndexToLineNumber(int abs_i){
  double lon=lonInRad(abs_i);
  double lon_offset;
  if(lon_inc_) lon_offset=lon-start_lon_;
  else lon_offset= start_lon_-lon;
  while(lon_offset<0) lon_offset+=2*pi;
  lon_offset/=res_;
  return((int)lon_offset); 
}

double LatLonGrid::lineNumberToLonInRad(int lnum)
{
  double lon=start_lon_;
  if(lon_inc_) lon+=(lnum+0.5)*res_;
  else lon-=(lnum+0.5)*res_;
  while(lon>2*pi)lon-=2*pi;
  while(lon<0) lon+=2*pi;
  return(lon);
}


void LatLonGrid::advance()
{
  if(lon_inc_){
    abs_start_valid_i_++;
  }
  else{ 
    abs_start_valid_i_--;
  }
  setValidRange();
}


int LatLonGrid::firstValidLongitudeIndex()
{
  return(abs_start_valid_i_%num_lons_);
}

int LatLonGrid::lastValidLongitudeIndex()
{
  return(abs_end_valid_i_%num_lons_);
}

bool LatLonGrid::checkValidLongitude(int abs_i){
  bool retval=true;
  int mini,maxi;
  if(lon_inc_) {
    mini= abs_start_valid_i_;
    maxi= abs_end_valid_i_;
  }
  else{
    maxi= abs_start_valid_i_;
    mini= abs_end_valid_i_;
  }
  while(abs_i < mini ) abs_i+=num_abs_lons_;
  while(abs_i > maxi ) abs_i-=num_abs_lons_;

  if(abs_i < mini ) {
    retval=false;
  }
  return(retval);
}
