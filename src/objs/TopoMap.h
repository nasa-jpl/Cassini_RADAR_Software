# ifndef  TOPOMAP_H
#define TOPOMAP_H
#include"Config.h"
#include<string>

class TopoMap{

 public:
  TopoMap();
  int config(Config& cfg, string file);
  int read(string file);
  int write(string file);
  int writeHeightsOnly(string file, bool spmode);
  int add(double slat, double slon, float h, float dh, bool average, int fileno);
  int addFile(char* fname,int typ,double tca);
  float getHeightInKm(double lat, double lon, bool& inbounds);
  float getInterpolatedHeightInKm(double lat, double lon, bool& inbounds);
  float getNearbyHeightInKm(double lat, double lon, bool& inbounds, float radbound_radians);

  ~TopoMap();

 private:
  int  init();  
  bool polar_proj;
  bool wraparound;
  int** num;
  float** height;
  double lon;
  double lat;
  int nLats;
  int nLons;
  float minLon;
  float minLat;
  float maxLon;
  float maxLat;
  float minPLon;
  float minPLat;
  float maxPLon;
  float maxPLat;
  float latStep;
  float lonStep;
  float dhthresh;
  float maxNeighborDist;
  float gridLonStep;
  float gridLatStep;
  int gridThickness;
  int projectPolar(double& slat,double& slon);
  int unprojectPolar(double& plat,double& plon);
  int allocateArrays();
  int setGridLines(float** h);

};


#endif
