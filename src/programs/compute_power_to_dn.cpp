//--------------------------------------------------
//compute power to dn when sc is pointing cold space
//--------------------------------------------------
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
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
#include "RADFunctions.h"
#include "Radiometer.h"
using std::cout;
using std::cerr;
using std::endl;
using std::set_unexpected;
using std::terminate;

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
    if (argc != 4)
      {
	cerr << "Usage: " << command << " config_file unit mode plot_mode(x or prefix_str)" << endl;
	exit(-1);
      }

    int clidx = 1;
    const char* config_file = argv[clidx++];
    char* units_mode = argv[clidx++];
    char* plot_mode=argv[clidx++];

    //----------------------
    //No range flag to be used while accessing L1BP
    //---------------------
    bool range_flag = false;
    double rangeLow=0.0;//no use of range flag
    double rangeHi=0.0;//no use of range flag

    //------------------------------
    // Load configuration parameters
    //------------------------------
    set_unit_option(units_mode);
  
    //-----------------------------------------------------------
    // Set plot options 
    //-----------------------------------------------------------
    cout << "Building plots" << endl;
    string plot_target;
    string plot_mode_str = plot_mode;
    string prefix_str ="";
    if (plot_mode_str=="x")
      {
	plot_target = "x";
      }
    else
      {
	plot_target = "eps-landscape";
	prefix_str = plot_mode;
      }
    
    //----------------------
    //config setup
    //----------------------
    Config cfg(config_file);

    //-----------------------------------------------
    // Load spice files & intermediate geometry file if there is any
    //------------------------------------------------
    Frame::config(cfg);
    BurstData::config(cfg);

    //-------------------------------------
    //First assign color to beam
    //-------------------------------------
    Array1D<string.h>  color_table("beam color",5);//beam color
    color_table(0) = "black";
    color_table(1) = "red";
    color_table(2) = "green";
    color_table(3) = "blue";
    color_table(4) = "magenta";



    //---------------------------
    //targets in solar system
    //----------------------------
    unsigned int Ntarget = 1;
    string target= cfg.str("target");//assume single target  
    Frame j2000c("J2000","Cassini");
   
    Array1D<string.h> target_names("target names");
    Array1D<Uvar> target_radii("target radii");
    Array1D<Frame> target_frames("target frames");
    Array1D<SpiceInt> target_ids("target ids");
    if(cfg.keywordExists("Number_of_targets"))
      {//multiple target
	Ntarget=cfg.getInt("Number_of_targets");
	target_names.resize(Ntarget);
	target_radii.resize(Ntarget);
	target_frames.resize(Ntarget);
	target_ids.resize(Ntarget);
	//-----------------------------------------
	//default target, save it to first element
	//---------------------------------------
	target_names(0)=cfg.str("target");
	for (unsigned int i_target=1;i_target<Ntarget;i_target++)
	  target_names(i_target)=cfg.str("target_"+toStr(i_target));
	for(unsigned int i=0;i<Ntarget;++i)
	  {
	    TargetGeom tg;
	    tg.setTarget(target_names(i));
	    target_radii(i)=tg.radius();
	    cout<<"target radius "<<target_names(i)<<" "<<target_radii(i)<<endl;   
	    target_frames(i)=Frame("J2000",target_names(i));
	    spice_target_id(target_names(i),target_ids(i));
	  }
      }
    else
      {//single default target
	Ntarget=1;
	target_names.resize(Ntarget);
	target_radii.resize(Ntarget);
	target_frames.resize(Ntarget);
	target_ids.resize(Ntarget);
	target_names(0)=cfg.str("target");//default target
	for(unsigned int i=0;i<Ntarget;++i)
	  {
	    TargetGeom tg;
	    tg.setTarget(target_names(i));
	    target_radii(i)=tg.radius();
	    cout<<"target radius "<<target_names(i)<<" "<<target_radii(i)<<endl;   
	    target_frames(i)=Frame("J2000",target_names(i));
	    spice_target_id(target_names(i),target_ids(i));
	  }
      }
    cout<<"Number of targets listed in config file "<<Ntarget<<endl;
   

  

    //------------------------------------------
    //Construct Beams and Frames
    //-------------------------------------------
    vector<Beam> Bvector;
    vector<Frame> Fvector;
    Uvec elev_beam3dB_oneway("beam 3db oneway along elev ",5);
    Uvec azi_beam3dB_oneway("beam 3db oneway along azi",5);
    Bvector.clear();
    Fvector.clear();
    //make all 5 beams and store
    for (unsigned int i = 0; i < 5;i++)
      {
	Beam beam(i+1,cfg,true);
	Bvector.push_back(beam); 
	Frame frame("CASSINI_RADAR_"+toStr(i+1),"Cassini");
	Fvector.push_back(frame);
	elev_beam3dB_oneway(i)=beam.getElevationWidthOneWay();
	azi_beam3dB_oneway(i)=beam.getAzimuthWidthOneWay();
      }
 

   
    //------------------------------------
    // Set SBDR file name
    //------------------------------------
    //    BurstData* data;
    string l1b_p_filename = cfg.str("L1B_P_filename");
    string l1b_a_filename = cfg.str("L1B_A_filename");
    L1B sbdr(l1b_p_filename,"r+","passive");
    L1B lbdr(l1b_a_filename,"r+","active");
    unsigned int Nrecord = sbdr.recordCount();
    sbdr.mapRecords();
    lbdr.mapRecords();
  

    //----------------------------------------------------------------
    //(0) Setting up resulting variables which are being read from L1BP
    //-----------------------------------------------------------------
    Uvec sab_counter("sab counter",Nrecord);
    Array1D<unsigned int> quality_flag("quality flag",Nrecord);
    Array1D<unsigned int>  beam_number("beam number",Nrecord);
    Array1D<Time> time("time",Nrecord);
    Array1D<unsigned int> rmode("",Nrecord);

    
    unsigned int Nparam=6;
    double *return_doubles; 
    return_doubles =new double[Nparam] ;
    char* param_list[Nparam];
    param_list[0]="sab_counter";
    param_list[1]="sclk";
    param_list[2]="brst";
    param_list[3]="beam_number";
    param_list[4]="r_mode";
    param_list[5]="engineer_qual_flag";
    unsigned int i_count = 0;
    double sclk_double, brst_double;
   
     while(!sbdr.eof()){
       if (!sbdr.returnParams(return_doubles,param_list,
				 Nparam,range_flag,rangeLow,
				rangeHi))
	{
	  sab_counter(i_count)= round_double(*return_doubles);
	  sclk_double=*(return_doubles+1);
	  brst_double=*(return_doubles+2);
	  beam_number(i_count)=(unsigned int) round_double(*(return_doubles+3));
	  rmode(i_count)=(unsigned int) round_double(*(return_doubles+4));
	  time(i_count).setSclk("Cassini",(unsigned long) round_double(sclk_double));
	  time(i_count)+= Uvar(brst_double,"s");
	  quality_flag(i_count)=(unsigned int) round_double(*(return_doubles+5));
	  i_count++;
	}
     }
     if(i_count !=Nrecord) ErrorMessage("rad_proc::expected number of records is not read in").throwMe();
     cout<<"Finished extracting radiometer relevant paramenters from L1BP"<<endl;     
    
     /*
     Plot a;
     Uvec tmp_y("",Nrecord);
     tmp_y=rmode;
     a.addXY(sab_counter,"",tmp_y," ",line("solid","red",1),sym("none"));
     a.setTitle("radar mode");
     a.setXlabelandUnits("sab number");
     a.setYlabelandUnits("r mode");
     a.show("x");
     */
	
     
     //----------------------------------------------------------------
     //(3) compute the angle  between limb direction and boresight direction for solar targets
     // in units of beam oneway 3dB
     // positive: boresight off the target
     // negative: boresight on the target
     //----------------------------------------------------------------
     Array2D<double> solar_target_distances("min target distances",
					    Nrecord,
					    Ntarget);//records * targets
     solar_target_distances=0.0;
     
     //recycled variables for speed up
     DirectionVector boresight,target_dir;
     PositionVector target_position;
     Uvar distance, limb_angle, boresight_angle;
     double dummy_a, dummy_b;
     for(unsigned int i_rec=0;i_rec<Nrecord;i_rec++)
       {
	 if(quality_flag(i_rec)==1 ||
	    quality_flag(i_rec)==2 ||
	    beam_number(i_rec)==0) continue;
	 boresight=DirectionVector("bore",Fvector[beam_number(i_rec)-1],time(i_rec),0,0,1);
	 for(unsigned int i_target=0;i_target<Ntarget;++i_target)
	   {
	     target_position=PositionVector("",target_frames(i_target),time(i_rec),0,0,0).representIn(j2000c);
	     target_dir = target_position;
	     distance = target_position.magnitude();
	     limb_angle =asin(target_radii(i_target)/distance);
	     boresight_angle = target_dir.angle(boresight);
	     dummy_a =((boresight_angle-limb_angle)/elev_beam3dB_oneway(beam_number(i_rec)-1)).getInUnits("");
	     dummy_b=((boresight_angle-limb_angle)/azi_beam3dB_oneway(beam_number(i_rec)-1)).getInUnits("");
	     //take smaller value
	     solar_target_distances(i_rec,i_target)=dummy_a;
	     if(dummy_b < solar_target_distances(i_rec,i_target)) solar_target_distances(i_rec,i_target)=dummy_b;
	   }
       }

     
       //debug::plot beam distance
     /*
       if(target!="")
	 {
	  Plot a;
	  Uvec distance_to_target(" ",Nrecord);
	  for(unsigned int i=0;i<Ntarget;++i)
	    {
	      distance_to_target=0.0;
	      for(unsigned int j=0;j<Nrecord;++j)
		distance_to_target(j)=solar_target_distances(j,i);
	      a.addY(distance_to_target,line("solid",color_table(i),1),sym("none"));
	    }
	  a.setTitle(" beam distance to specified targets in config");
	  a.setFile(prefix_str +"solar_targets.eps");
	  a.setCmdFile(prefix_str+ "solar_targets.xmgr");
	  a.show(plot_target);
	 }

     */

       //-------------------------------------------------------------------
       //(4) determine the  nearest target among solar and microwave targets
       //----------------------------------------------------------------------
       Array1D<int>  target_record("sab target",Nrecord);//each beam's target
       Array1D<double> smallest_beam_distance("target whose beam distance is smallest",Nrecord);
       Array1D<double> largest_beam_distance("",Nrecord);
       //assume default target, each time it computes
       //target id which has the nearest range
       for (unsigned int i_rec=0;i_rec<Nrecord;i_rec++)
	 {
	   //select min target distance among target names
	   // variable: solar_target_distances
	   unsigned int index_min_distance = target_ids(0);
	   double min_distance = solar_target_distances(i_rec,0);

	   unsigned int index_max_distance = target_ids(0);
	   double max_distance = solar_target_distances(i_rec,0);

	   
	   //for solar targets
	   for (unsigned int i_target=0;i_target<Ntarget;++i_target)
	     {
	       if(min_distance>solar_target_distances(i_rec,i_target))
		 {
		   index_min_distance = target_ids(i_target);
		   min_distance = solar_target_distances(i_rec,i_target);
		 }

	       if(max_distance<solar_target_distances(i_rec,i_target))
		 {
		   index_max_distance = target_ids(i_target);
		   max_distance = solar_target_distances(i_rec,i_target);
		 }
	     }
	   target_record(i_rec)=index_min_distance;
	   smallest_beam_distance(i_rec)= min_distance;
  	   largest_beam_distance(i_rec)= max_distance;
	 }//loop over record













       //when pointing cold space, write out into text file
       Fvec echo("");
       for(unsigned int i=0;i<Nrecord;++i){
	 if(rmode(i)>3 && rmode(i)<8) continue;;//not active mode
	 if(!lbdr.foundSABRecord((unsigned long) round_double( sab_counter(i).getInUnits("")))) continue;//do not care
	 lbdr.readRecord();
	 double adc= lbdr.adc.getInUnits("Hz");
	 double att;
	 unsigned int beam_id= lbdr.beam_number -1;
	 if(beam_id <2) att=lbdr.at1.getInUnits("");
	 else if(beam_id==2) att=lbdr.at3.getInUnits("");
	 else att=lbdr.at4.getInUnits("");
	 
	 Uvar rcv;
	 if(adc==2.0e6)
	   rcv= sarh_rcv_bandwidth;
	 else if(adc==1.0e6)
	   rcv= sarl_rcv_bandwidth;
	 else if(adc==250e3)
	   rcv= altl_rcv_bandwidth;	  
	 else if (adc==10e6)
	   rcv= alth_rcv_bandwidth;
	 else{
	   cout<<"Adc "<< adc<<endl;
	   ErrorMessage("invalid adc ").throwMe();
	 }
	 //------------
	 //lbdr length
	 //--------------
	 unsigned int Nlength=lbdr.Nradar_data;
	 if(Nlength==0) continue;
	 echo.resize(Nlength);
	 for(unsigned int k=0;k<Nlength;++k)
	   echo(k)=lbdr.radar_data(k);
	 if(lbdr.baq_mode==3)
	   echo /= lbdr.pul;
	 else{
	   echo -= echo.mean();
	 }

	 if(lbdr.csr==1){
	   cout<<round_double( sab_counter(i).getInUnits(""))<<" "<<lbdr.beam_number<<" "" "<<smallest_beam_distance(i)<<" "<< att<<" "<< rcv.getInUnits("KHz")<<" "<< echo.mean()<<" "<< echo.rms(echo.size())<<endl;
	 }
       }//loop over the record

     
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


 
