//---------------------------------------------
// Estimates j2000 coordinates from TBF lat/lons
//---------------------------------------------
// the is the code for both get_tiepoint_info
// and get_tiepoint_info_from_lbdr
// If it is the command is get_tiepoint_info it reads
// from an sbdr file if it is get_tiepoint_info_from_lbdr
// then it reads from a lbdr file
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"
#include "L1B.h"
#include "config_keywords.h"
#include "DebugInfo.h"
#include "Beam.h"
#include "Plot.h"
#include "Constants.h"
#include "TargetGeom.h"
#include "SARFunctions.h"
#include "DebugInfo.h"
#include "BurstData.h"
#include "BIDRFile.h"
#include "TopoMap.h"
#include "Array.h"

using std::cout;
using std::cerr;
using std::endl;
using std::set_unexpected;
using std::terminate;
using std::isfinite;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  terminate();
  }




int main(int argc, char* argv[])

{

set_unexpected(myUnexpected);
typedef list<Uvec>::const_iterator UVECI;
 
try
  {

  //------------------------
  // Parse the command line
  //------------------------  

    const char* command = argv[0];
    if (argc != 7 && argc !=9)
      {
        if(strcmp(argv[0],"get_tiepoint_info_from_lbdr")==0)
	  cerr << "Usage: " << command << " cfg_file lbdrfile bidrfile(should be beammask) topomapfile(or none) line sample [lat] [wlon]" << endl;
	else
	  cerr << "Usage: " << command << " cfg_file sbdrfile bidrfile(should be beammask) topomapfile(or none) line sample [lat] [wlon]" << endl;
	exit(-1);
      }

    int clidx = 1;

    
    const char* config_file = argv[clidx++];
    const char* bodp_file = argv[clidx++];
    const char* bidr_file = argv[clidx++];
    const char* tmap_file = argv[clidx++];

    double line = atof(argv[clidx++]);
    double sample = atof(argv[clidx++]);
    int iline=(int)floor(line+0.5);
    int isample=(int)floor(sample+0.5);

    bool optchecklatlon=false;
    double clat,clon;
    double tol=degtorad/512;
    double lontol=tol/cos(clat);
    if(argc==9){
      clat=atof(argv[clidx++])*degtorad;
      clon=2*pi-atof(argv[clidx++])*degtorad;
      optchecklatlon=true;
    }

    //------------------------------
    // Load configuration parameters
    //------------------------------
    set_unit_option("none");
  

    //----------------------
    //config setup
    //----------------------
    Config cfg(config_file);

    //-----------------------------------------------
    // Load spice files & intermediate geometry file if there is any
    //------------------------------------------------
    Frame::config(cfg,true);
    BurstData::config(cfg);



    //-----------------
    //Debug info setup
    //------------------
    DebugInfo::config(cfg);
    DebugInfo dbg("main");
    dbg.file<<"latlon_to_j2000 process debug level is "<<dbg.level<<endl;

    Uvar carrier_freq=cfg["carrier_frequency"];

    //---------------------------
    // set up frames
    //----------------------------
    Frame j2000("J2000","Titan");
    Frame tbf("IAU_Titan","Titan");


    //--------------------------
    // Set up beam and beam frames usng CassiniSim
    //--------------------------
    CassiniSim cassini_sim; // Object for handling beam and target geometry
                          // during a Cassini data take
 
     vector<Beam> beam;       // Beam objects-- trying to get rid of dependency on CassiniSim

     cassini_sim.config(cfg);
    // Setup beams
    bool read_beam_pattern=true;

    beam.resize(5);
    float maxgain[5];  
    for (unsigned int i=1; i <= 5; ++i)
      {
	beam[i-1] = Beam(i,cfg,read_beam_pattern);
	Uvar elev=0;
	Uvar azim=beam[i-1].maxGainLineToAzim(elev);
	maxgain[i-1]=beam[i-1].getMaxGain();
      }


    string msg;
    //--------------------------
    // Convert line sample to latlon check lat lon values if check values available
    //---------------------------
    if (!fileExists(bidr_file))
    {
      msg = "Input file " + string(bidr_file) + " not found or unreadable";
      ErrorMessage e(msg);
      e.throwMe();
    }
    BIDRFile bf_in(bidr_file,"r");
    bf_in.readHeader();

    if (line < 1.0 || line > (double) bf_in.numLons())
      {
	msg = "Invalid line " + toStr(line) + ".  Line must be between" +
	  " 1 and " + toStr(bf_in.numLons()) + ".";
	ErrorMessage e(msg);
	e.throwMe();
      }
    if (sample < 1.0 || sample > (double) bf_in.numLats())
      {
	msg = "Invalid sample " + toStr(sample) + ".  Sample must be between" +
	  " 1 and " + toStr(bf_in.numLats()) + ".";
	ErrorMessage e(msg);
	e.throwMe();
      }
    
    // get beam number from beammask file 
    char maskvalue;
    unsigned int offset=sizeof(char)*(isample-1)+(iline-1)*sizeof(char)*bf_in.numLats();
    offset+=bf_in.labelLength();
    bf_in.setPosition(offset);
    bf_in.read(maskvalue);
    int beam_number=0;
    switch ((int)maskvalue){
    case 1:
      beam_number=1;
      break;
    case 2:
      beam_number=2;
      break;
    case 4:
      beam_number=3;
      break;
    case 8:
      beam_number=4;
      break;
    case 16:
      beam_number=5;
      break;
    default:
      beam_number=0;
      break;
    }

    if(beam_number==0){
      char trymask=0;
      int c = 0;
      while(trymask==0){
	  c=c+1;
	  bf_in.setPosition(offset+c);
          bf_in.read(trymask);
          if(trymask==0){
            c=-c;
          }
	  bf_in.setPosition(offset+c);
          bf_in.read(trymask);
          if(trymask==0){
            c=-c;
	  }
      }
      printf("INVALID_TIEPOINT %d pixels from edge",c);
      exit(0);
    }
    // Calculate the oblique cylindrical latitude and longitude in radians at
    // the selected point.
    double lat = bf_in.OCLatFromGrid(sample);
    double lon = bf_in.OCLonFromGrid(line);
    
    // Convert the latitude and longitude to standard coordinates.
    double slat = bf_in.standardLatInRad(lon, lat);
    double slon = bf_in.standardLonInRad(lon, lat);
  
    //Check lat lon if check values available
    if(optchecklatlon){
      double dlat=fabs(slat-clat);
      double dlon=fabs(slon-clon);
      if(dlat>tol || dlon>lontol){
	fprintf(stderr,"%s:Latlon tolerance not met Err=(dlat,dlon)=(%g,%g)deg Tolerance=(%g,%g)deg\n",argv[0],dlat*radtodeg,dlon*radtodeg,tol*radtodeg,tol*radtodeg);
        exit(1);
      }
    }

    char* actstr="active";
    char* passstr="passive";
    char* tstr=passstr;
    if(strcmp(argv[0],"get_tiepoint_info_lbdr")==0){
      tstr=actstr;
    }

    //------------------------------------
    // Setup LBDR object
    //------------------------------------
    L1B l1b(bodp_file,"r",tstr);
    unsigned int Nrecord = l1b.recordCount();
    

    //------------------------------------
    // Get height from Topomap
    //------------------------------------
    
    TopoMap tmap;
    double height=0;
    if(strcasecmp("none",tmap_file)!=0){
      tmap.config(cfg,tmap_file);
      bool inbounds=false;
      float radbound=0.038835; //  radians of latitude upper bound distance to look for SARTopo (100 km) 
      height=tmap.getNearbyHeightInKm(slat,slon,inbounds,radbound);
      if(!inbounds){
	fprintf(stderr,"Warning %s:lat,wlon=(%g,%g) is not in topomap bounds\n",argv[0],slat*radtodeg,360-slon*radtodeg);
	height=0.0;
      }
    }

    PositionVector pos_j2000;

    // find sab with highest relative gain at  Latlon and compute all values

    double bestgain=0; 
    double x,y,z,xc,yc,zc,vx,vy,vz;
    double range,fdop;
    double lamb;
    double best_ground_impact_time_in_s;
    int sabnum;
    for(int c=0;c<(int)Nrecord;c++){

      l1b.readParameter("rc_bw"); // skip nonSAR bursts
      l1b.readParameter("sab_counter");
      if(l1b.rc_bw<Uvar(400000,"Hz")||l1b.rc_bw>Uvar(2000000,"Hz"))
	{ 
          //fprintf(stderr,"sab %d  rc_bw %g\n",l1b.sab_counter,l1b.rc_bw.getInUnits("Hz"));
	  l1b.skipRecord();
	  continue;
	}     

      l1b.readParameter("beam_number");
      if(l1b.beam_number!=beam_number)
	{ 
	  l1b.skipRecord();
	  continue;
	}     
  
      l1b.readParameter("ctrx");
      l1b.readParameter("pri");
      l1b.readParameter("t");
      l1b.readParameter("pul");
      l1b.readParameter("rwd");
      l1b.readParameter("act_centroid_lat");
      l1b.readParameter("act_centroid_lon");
      l1b.readParameter("fast_csf");
      l1b.readParameter("slow_cfs");
      l1b.readParameter("csq");

      // use time halfway between transmit window center and receive window
      // center to get nominal_doppler_centroid and nominal_delay
      Uvar receive_window=l1b.ctrx*l1b.pri;
      Uvar transmit_window=l1b.pul*l1b.pri;
      Uvar transmit_center=transmit_window/2;
      Uvar receive_center=l1b.rwd+receive_window/2;
      Uvar time_offset=(transmit_center+receive_center)/2;

 
      // Choice of ground impact time has a very small impact on range
      // Even if it is off several pri ( 1 ms) if will be at most 3 m off.
      Time ground_impact_time=l1b.t+time_offset;
      StateVector scstate, scstate_tbf;
      double ground_impact_time_in_s;
      ground_impact_time.getEt(ground_impact_time_in_s);

      DirectionVector dc(tbf,ground_impact_time,0,0,0);
      DirectionVector dc2(tbf,ground_impact_time,0,0,0);
      dc.setPlanetodetic(l1b.act_centroid_lat,l1b.act_centroid_lon);
      dc.representIn(j2000);
      dc2.setPlanetodetic(slat,slon);
      dc2.representIn(j2000);

 
      PositionVector pos(dc);
      pos=pos*2575;
       
      PositionVector pos2(dc2);
      pos2=pos2*2575;
      PositionVector pos2_j2000=pos2.representIn(j2000);
      j2000.ephemeris(scstate,cassini_spice_id,ground_impact_time_in_s);
      tbf.ephemeris(scstate_tbf,cassini_spice_id,ground_impact_time_in_s);
      PositionVector pos_j2000=scstate.position();
      PositionVector look=pos2_j2000-pos_j2000;
    
      DirectionVector ulook(look);

      DirectionVector lookinbeam=ulook.representIn(Frame(beam_frame_spice_id[l1b.beam_number-1],cassini_spice_id));
      
      Uvar azim,elev;
      lookinbeam.getAzimuthElevation(azim,elev);


      // sanity check to avoid Beam pattern bug for ridiculously large angles
      if(fabs(azim)>Uvar(0.1,"rad") || fabs(elev) > Uvar(0.1,"rad") ){
	//fprintf(stderr,"sab %d pos2 %g %g %g pos %g %g %g look %g %g %g azim %g elev %g\n",l1b.sab_counter,pos2_j2000.km(PositionVector::X),pos2_j2000.km(PositionVector::Y),pos2_j2000.km(PositionVector::Z),pos_j2000.km(PositionVector::X),pos_j2000.km(PositionVector::Y),pos_j2000.km(PositionVector::Z),look.km(PositionVector::X),look.km(PositionVector::Y),look.km(PositionVector::Z),azim.getInUnits("deg"),elev.getInUnits("deg"));
	l1b.skipRecord();
	continue;
      }
      double gain=beam[l1b.beam_number-1].bilinear(azim,elev);
      gain=gain/maxgain[l1b.beam_number-1];

      //fprintf(stderr,"sab %d pos2 %g %g %g pos %g %g %g look %g %g %g azim %g elev %g gain %g\n",l1b.sab_counter,pos2_j2000.km(PositionVector::X),pos2_j2000.km(PositionVector::Y),pos2_j2000.km(PositionVector::Z),pos_j2000.km(PositionVector::X),pos_j2000.km(PositionVector::Y),pos_j2000.km(PositionVector::Z),look.km(PositionVector::X),look.km(PositionVector::Y),look.km(PositionVector::Z),azim.getInUnits("deg"),elev.getInUnits("deg"),gain);

      if(gain>bestgain){
	bestgain=gain;
	best_ground_impact_time_in_s = ground_impact_time_in_s;

	FloatVector vel_j2000=scstate.velocity();
	x=pos_j2000[PositionVector::X].getInUnits("km");
	y=pos_j2000[PositionVector::Y].getInUnits("km");
	z=pos_j2000[PositionVector::Z].getInUnits("km");
	
	vx=vel_j2000[FloatVector::X].getInUnits("km/s");
	vy=vel_j2000[FloatVector::Y].getInUnits("km/s");
	vz=vel_j2000[FloatVector::Z].getInUnits("km/s");
	  
	  
	xc=dc[DirectionVector::X];
	yc=dc[DirectionVector::Y];
	zc=dc[DirectionVector::Z];


	Uvar chirp_freq=carrier_freq;
	chirp_freq+=l1b.fast_csf; // add chirp start freq
	chirp_freq+=l1b.slow_cfs*(l1b.csq-1)/2; // add half of chirp bandwidth
	Uvar lambda_chirp=speed_light/chirp_freq; // convert to wavelength
	lamb=lambda_chirp.getInUnits("km");   
	range=(look.magnitude()).getInUnits("km");
	look=look.representIn(tbf);
	DirectionVector lookdir(look);
	FloatVector vel_tbf=scstate_tbf.velocity();
	fdop=(2.0*dot(vel_tbf,lookdir)/lamb).getInUnits("Hz");      
        sabnum=l1b.sab_counter;
      } // end if closest to latlon so far
     
     l1b.skipRecord();
    } // end loop over records

    printf("%23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.9f %23.15f %23.9f %d\n",range,fdop,x,y,z,vx,vy,vz,xc,yc,zc,best_ground_impact_time_in_s,lamb,height,sabnum); 
  return(0);
  }
 
 
catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
    cerr << "Non Error Message Exception caught" << endl;
  }

 return(0);
 
}

