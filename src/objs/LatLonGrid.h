#ifndef LATLONGRID_H
#define LATLONGRID_H


class LatLonGrid{

 public:
  LatLonGrid();
  LatLonGrid(int pixels_per_degree, double min_lat_deg, double max_lat_deg,
	     double start_lon_deg, double lon_width_deg, bool lon_increasing,
	     bool lon360, int lon_buffer_in_deg);
  void getBoundary(const double lonlat_bounds[4], int ij_bounds[4]);
  double latInRad(int j);
  double lonInRad(int abs_i);
  double lineNumberToLonInRad(int lnum);
  int absIndexToLineNumber(int abs_i);
  void advance();
  bool checkValidLongitude(int abs_i);
  int absStartValidIndex();
  int absEndValidIndex();
  int firstValidLongitudeIndex();
  int lastValidLongitudeIndex();
  int relLonIndex(int abs_i);
  int numLats(){return(num_lats_);}
  int numLons(){return(num_lons_);}
  int totalNumLons(){return(total_num_lons_);}
  double minLat(){return(min_lat_);}
  double startLon(){return(start_lon_);}
  bool lonIncreasing(){return(lon_inc_);}
 private:
  void setValidRange();
  int abs_start_valid_i_;
  int abs_end_valid_i_;
  int num_lats_;
  int num_lons_;
  int num_abs_lons_;
  double min_lat_;
  double res_;
  bool lon_inc_;
  double start_lon_;
  int total_num_lons_;
};
#endif
