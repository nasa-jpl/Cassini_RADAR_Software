//-----------------------------
// Routine to create an extra BIDR file with the same format as the primary sigma0 BIDR
// but with values chosen from the appropriate SAB either from the SBDR or from
// a ASCII file with comma or space separated columns in which one column is sab_number and the other is the image value
// Bryan Stiles Aug 31 2012
// If the column numbers for the ASCII file are not specified assumes it is an SBDR file


#include <string.h>
#include <iostream>
#include <stdlib.h>
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

#define BAD_VALUE_HEX 0xff7ffffb

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

    // set up BAD_VALUE
    int v=BAD_VALUE_HEX;
    float* p=(float*)&v;
    float BAD_VALUE =*p;

    //------------------------
    // Parse the command line
    //------------------------  

    const char* command = argv[0];
    if (argc != 8 && argc !=13)
      {
	  cerr << "Usage: " << command << " cfgfile sabordered_file(sbdr or ASCII) paramname(only used for sbdr case)"<< endl; 
	  cerr << "      startnum_bidrfile endnum_bidrfile s0_bidrfile(or just header) outbidrfile [ASCII_sabcol] [ASCII_valcol] [ASCII_latcol] [ASCII_wloncol] [ASCII_beamcol]" << endl;
	exit(-1);
      }

    int clidx = 1;

    const char* cfgfile = argv[clidx++];    
    const char* bodp_file = argv[clidx++];
    const char* param_name = argv[clidx++];
    const char* bidr1_file = argv[clidx++];
    const char* bidr2_file = argv[clidx++];
    const char* s0_bidr_file = argv[clidx++];
    const char* outfile = argv[clidx++];
    int sabcol=0,valcol=0,loncol=0,latcol=0,beamcol=0;
    bool use_sbdr=true;
    if(argc==13){
      sabcol= atoi(argv[clidx++]);
      valcol = atoi(argv[clidx++]);
      latcol = atoi(argv[clidx++]);
      loncol = atoi(argv[clidx++]);
      beamcol = atoi(argv[clidx++]);
      use_sbdr=false;
    }

    //------------------------------
    // set unit options (none is used to speed things up)
    //------------------------------
    set_unit_option("none");
    Config cfg("make_extra_backplane",cfgfile);

    //---------------------------------
    // Load spice files & intermediate geometry file
    //---------------------------------
    Frame::setGeometryMode(Frame::READ_FILE);
    bool needs_ckernel=true;
    Frame::config(cfg,needs_ckernel); 
    BurstData::config(cfg);

    //-------------------
    // Open StartNum BIDR file
    //---------------------- 
    string msg;
    if (!fileExists(bidr1_file))
    {
      msg = "Input file " + string(bidr1_file) + " not found or unreadable";
      ErrorMessage e(msg);
      e.throwMe();
    }
    BIDRFile bf_in(bidr1_file,"r");
    bf_in.readHeader();


    //-- open s0 BIDR file
    FileMgr sfm(s0_bidr_file,"r");

    //-- read s0 BIDR header
    PDSLabel label(sfm);


    // open output file
    FileMgr ofm(outfile,"w");

    //-- write s0 BIDR header
    if(label.label() == EMPTY_STRING){
      cerr << "Fatal error:"<< command << ": input file " << s0_bidr_file 
	   << " does not appear to be PDS." << endl;
      exit(1);
    }

    ofm.write(label.label());

    // close s0 BIDR input file
    sfm.close();


    // allocate sabcounter map array    
    int Nlats=bf_in.numLats();
    int Nlons=bf_in.numLons();
    int** minsabmap=(int**)make_array(sizeof(int),2,Nlons,Nlats);
    int** maxsabmap=(int**)make_array(sizeof(int),2,Nlons,Nlats);


     
    // get minsabmap
    int s;
    for(int i=0;i<bf_in.numLons();i++){
      for(int j=0;j<bf_in.numLats();j++){
	bf_in.read(s);
        minsabmap[i][j]=s; 
      }
    }

    // get maxsabmap
    BIDRFile bf_in2(bidr2_file,"r");
    bf_in2.readHeader();
    for(int i=0;i<bf_in2.numLons();i++){
      for(int j=0;j<bf_in2.numLats();j++){
	bf_in2.read(s);
        maxsabmap[i][j]=s; 
      }
    }

    // allocate list of vals,lats,lon, and beam numbers indexed by SAB
    int MAXNSABS = 100000;
    int* beam=(int*)calloc(MAXNSABS,sizeof(int)); // initialized to zero
    float* val=(float*)malloc(sizeof(float)*MAXNSABS);
    float* lat=(float*)malloc(sizeof(float)*MAXNSABS);
    float* lon=(float*)malloc(sizeof(float)*MAXNSABS);
    int minsab=100000;
    int maxsab=0;

    //----------------
    // read in list of vals,lats,lons,and beamnumbers
    //----------------
    if(use_sbdr){  

      //------------------------------------
      // initialize input file
      //------------------------------------
      L1B l1b(bodp_file,"r","active");
      l1b.readHeader();
      
      // Configure SarProcessorParams
      // includes Control Booleans, Processing constants and Calibration
      // coefficients
      SARProcParams spp(cfg,l1b.maxSabCounter()); 
      
      // Construct OblCyl projection/
      StateVector s=spp.getClosestApproachState();
      OblCylProj proj(s);
      
      while(!l1b.eof()){    

	// Reading individual parameters instead of whole record 
	// for faster (but more fragile) access 
	// 
	l1b.readParameter("sclk"); // needed for PDS label start/stop timestamps
	l1b.readParameter("brst"); // needed for PDS label start/stop timestamps
	l1b.readParameter("t"); // needed for PDS label start/stop timestamps
	l1b.readParameter("beam_number"); // needed for single beam mode check
	l1b.readParameter("adc"); // needed for isSAR check
	l1b.readParameter("csr"); // needed for isCAl check in badSARData
	l1b.readParameter("engineer_qual_flag"); // needed for badSARData check
	l1b.readParameter("science_qual_flag"); // needed for badSARData check
	l1b.readParameter("act_centroid_lon"); // needed to compute lat bounds
	l1b.readParameter("act_centroid_lat"); // needed to compute lat bounds
	l1b.readParameter("record_id"); // needed for spp::inRegion
	l1b.readParameter("time_from_closest_approach"); // needed for spp.inRegion
	l1b.readParameter("time_from_epoch"); // needed for spp.inRegion
	l1b.readParameter("sab_counter"); // needed for quality override check
	l1b.readParameter("dtn"); // needed for PDS label product ID
	l1b.readParameter("num_pulses_received"); // needed for check in skipBurst
        l1b.readParameter(param_name);
	// Nominal skipping of bursts (also implemented in BIDR::BIDR)
	if(spp.skipBurst(proj,l1b)){
	  l1b.skipRecord();
	  continue;
	}
        
	int offset=l1b.sab_counter;
        if(offset>=MAXNSABS){
	  fprintf(stderr,"Found sab_counter %d >= MAXNSABS(%d)\n",offset,MAXNSABS);
	  exit(1);
	}
	if(offset>maxsab) maxsab=offset;
        if(offset<minsab) minsab=offset;
        beam[offset]=l1b.beam_number;
        val[offset]=l1b.getParam(param_name);
        lat[offset]=get_in_base_units(l1b.act_centroid_lat); // radians
        lon[offset]=get_in_base_units(l1b.act_centroid_lon); // radians positive east

        l1b.skipRecord();
      } // end loop through l1b
    } // end SBDR file case

    // ASCII column file case
    else{
      // determine last column needed
      float fcols[200];
      int icols[200];
      int maxcol=sabcol;
      if(valcol>maxcol)maxcol=valcol;
      if(latcol>maxcol)maxcol=latcol;
      if(loncol>maxcol)maxcol=loncol;
      if(beamcol>maxcol)maxcol=beamcol;
      if(maxcol>200){
	fprintf(stderr,"Max column entry %d exceeds 200\n",maxcol);
	exit(1);
      }
      // open ASCII input file
      FILE* ifp=fopen(bodp_file,"r");
      if(ifp==NULL){
	fprintf(stderr,"Error opening file %s\n",bodp_file);
	exit(1);
      }
      // loop over lines in file parsing each line
      char linestr[2000];
      while(1){
	fgets(linestr,2000,ifp);
	if(feof(ifp))break;
        char* str=strtok(linestr," \t,");
        fcols[0]=atof(str);
        icols[0]=atoi(str);
	for(int c=1;c<maxcol;c++){
	  str=strtok(NULL," \t,");
	  fcols[c]=atof(str);
	  icols[c]=atoi(str);
	}
	int offset=icols[sabcol-1];
        if(offset>=MAXNSABS){
	  fprintf(stderr,"Found sab_counter %d >= MAXNSABS(%d)\n",offset,MAXNSABS);
	  exit(1);
	}
	if(offset>maxsab) maxsab=offset;
        if(offset<minsab) minsab=offset;
	beam[offset]=icols[beamcol-1];
	val[offset]=fcols[valcol-1];
	lat[offset]=fcols[latcol-1]*degtorad;
	lon[offset]=2*pi-fcols[loncol-1]*degtorad; 
      } // end parsing loop
      fclose(ifp);
    } // end ASCII column input file case

    // account for special case of angle parameter (for now this only works for act_azimuth_angle act_inc
    bool isang=false;
    if(strstr(param_name,"lon")!=NULL || strstr(param_name,"angle")!=NULL) isang=true;
    //---- average and write out vals indexed by minsabmap through maxsabmap
    for(int i=0;i<Nlons;i++){
      for(int j=0;j<Nlats;j++){
	float outval=BAD_VALUE;
 	if(minsabmap[i][j]>0){
	  int sab1=minsabmap[i][j];
          int sab2=maxsabmap[i][j];
	  int n=0;
          outval=0;
          int b=beam[sab1];
	  float val1=val[sab1];

          for(int s=sab1;s<=sab2;s++){
	    if(beam[s]==b){
              if(isang){
		float dang=val1-val[s];
		if(dang>180) val[s]+=360;
		else if(dang<-180) val[s]-=360;
	      }
	      outval+=val[s];
	      n++;
	    }
	  }
	  outval/=n;
          if(isang){
	    while(outval>360)outval-=360;
	    while(outval<0)outval+=360;
	  }
	}
	ofm.write(outval);
      }
    }
    ofm.close();
    free_array(minsabmap,2,Nlons,Nlats);
    free_array(maxsabmap,2,Nlons,Nlats);
    free(beam);
    free(lat);
    free(val);
    free(lon);
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

