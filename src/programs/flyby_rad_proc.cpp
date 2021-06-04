//------------------------------
// radiometer processor
//  this program is to process radiometer data during 
//  Titan flyby
// Date: April 13, 2006
// Written by Y. Gim
//-------------------------------
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
#include "DebugInfo.h"
#include "BurstData.h"
#include "Utils.h"
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

    //----------------------------------------------
    //No range flag to be used while accessing L1BP
    //-----------------------------------------------
    bool range_flag = false;
    double rangeLow=0.0;//no use of range flag
    double rangeHi=0.0;//no use of range flag

    //---------------------------------------------------
    // Load configuration parameters
    //----------------------------------------------------
    set_unit_option(units_mode);
  
    //-----------------------------------------------------------
    // Set plot options 
    //-----------------------------------------------------------
    cout << "Building plots" << endl;
    string plot_target;
    string plot_mode_str = plot_mode;
    string prefix_str ="";
    string plot_tool ="xmgr";
    if (plot_mode_str=="x")
      {
	plot_target = "x";
      }
    else
      {
	plot_target = "eps-landscape";
	prefix_str = plot_mode;
      }

    //------------------------------------------
    //config setup
    //-----------------------------------------------
    Config cfg(config_file);


    //-----------------------------------------------
    // Load spice files & intermediate geometry file if there is any
    //------------------------------------------------
    Frame::config(cfg);
    BurstData::config(cfg);

    //------------------------------------------------------
    //First assign color to beam
    //-------------------------------------------------------
    Array1D<string.h>  color_table("beam color",5);//beam color
    color_table(0) = "black";
    color_table(1) = "red";
    color_table(2) = "green";
    color_table(3) = "blue";
    color_table(4) = "magenta";


    //--------------------------------------------------
    //Debug info setup
    //----------------------------------------------------
    DebugInfo::config(cfg);
    DebugInfo dbg("main");
    dbg.file<<"radiometer process degub level is "<<dbg.level<<endl;


    //---------------
    //conversion factors
    //-----------------
    double d2r= pi/180.0;
    double r2d= 180.0/pi;

    //------------------------------------------
    //Construct Beams and Frames
    //-------------------------------------------
    vector<Beam> Bvector;
    vector<Frame> Fvector;
    Uvec elev_beam3dB_oneway("beam 3db oneway along elev ",5);
    Uvec azi_beam3dB_oneway("beam 3db oneway along azi",5);
    Bvector.clear();
    Fvector.clear();
    vector<double> integrated_beam_gain;
   
    //make all 5 beams and store
    for (unsigned int i = 0; i < 5;i++)
      {
	Beam beam(i+1,cfg,true);
	Bvector.push_back(beam); 
	Frame frame("CASSINI_RADAR_"+toStr(i+1),"Cassini");
	Fvector.push_back(frame);
	elev_beam3dB_oneway(i)=beam.getElevationWidthOneWay();
	azi_beam3dB_oneway(i)=beam.getAzimuthWidthOneWay();

	double  directivity_dB, solid_angle, sum;
	beam.integrateBeamPattern(solid_angle, directivity_dB,sum);
	cout<<"solid angle "<< solid_angle/(4.0*pi)<<endl;
	cout<<"directivity "<< directivity_dB<<endl;
	cout<<"beam integration "<< sum<<endl;
	integrated_beam_gain.push_back(sum);
      }

    

    //---------------------------
    //targets in solar system: usually 1
    //----------------------------
    unsigned int Ntarget = 1;
    string target= cfg.str("target");//assume single target  
    Frame j2000c("J2000","Cassini");
   
    //------------------------------------------
    //target names, radii, frame, spice id
    //--------------------------------------------
    Array1D<string.h> target_name("target names",1);
    Array1D<Uvar> target_radius("target radii");
    Array1D<Frame> target_frame("target frames");
    Array1D<SpiceInt> target_id("target ids");
    Array1D<Uvar> target_temperature("target temperaturea");
    //flyby: just one target
    Ntarget=1;
    target_name.resize(Ntarget);
    target_radius.resize(Ntarget);
    target_frame.resize(Ntarget);
    target_id.resize(Ntarget);
    target_name(0)=cfg.str("target");//default target
    target_temperature(0)=cfg["target_temperature"];
    for(unsigned int i=0;i<Ntarget;++i)
      {
	TargetGeom tg;
	tg.setTarget(target_name(i));
	target_radius(i)=tg.radius();
	cout<<"target radius "<<target_name(i)<<" "<<target_radius(i)<<endl;   
	target_frame(i)=Frame("J2000",target_name(i));
	spice_target_id(target_name(i),target_id(i));
      }
      
   

    //---------------------------------------------------------
    // Fixed parameters from MJ's IDL program
    // 
    //----------------------------------------------------------
    Uvar Trec=cfg["receiver_temperature"];
    Uvar radiometer_cal=1.0/cfg["radiometer_cal_gain"];
    Uvar Tcoldsky=cfg["cold_sky_Temperature"];//Ku band 
    Uvar Troom=Uvar(300,"K");
    Uvar norm_gain= (Troom + Trec)/cfg["norm_rl_count"];
    unsigned int Npolyfit=(unsigned int) cfg.getInt("cold_sky_cal_polynomial_fit_order");//not sure this is useful

    cout<<"Receiver temp "<< Trec<<endl;
    cout<<"radiometer cal  "<< radiometer_cal<<endl;
    cout<<"cold sky "<< Tcoldsky<<endl;
    cout<<"norm gain "<< norm_gain<<endl;
    //-----------------------------------------------------------------
    // Parameters from config file
    //-------------------------------------------------------------------
    Uvar min_beam_distance_in_units_of_beam3dB=cfg["min_beam_distance"];
    double min_beam_distance =  min_beam_distance_in_units_of_beam3dB.getInUnits("");
    Uvar avg_time=cfg["cold_sky_average_time"];
    Uvar max_cold_sky_cal_separation=cfg["maximum_separation_bet_cold_sky_meas"];
   
  

    //----------------------------------------------------------------------
    //If nearest cold sky measurement occur before or later than 
    //max_cold_sky_cal_separation
    // System gain will be computed using resistive load measurement
    //----------------------------------------------------------------------
    cout<<"receiver temperature "<<Trec<<endl;
    cout<<"min beam distance between boresight and limb direction "<<min_beam_distance_in_units_of_beam3dB<<endl;
    cout<<"cold sky temperature "<<Tcoldsky<<endl;
    cout<<"cold sky average time "<< avg_time<<endl;
    cout<<"Max cold sky cal separation "<< max_cold_sky_cal_separation<<endl;
  

   
 

   
    //------------------------------------
    // Set SBDR file name
    //------------------------------------
    string l1b_p_filename = cfg.str("L1B_P_filename");
    string l1b_a_filename = cfg.str("L1B_A_filename");
    L1B data(l1b_p_filename,"r+","passive");
    unsigned int Nrecord = data.recordCount();
  
  

    //----------------------------------------------------------------
    //(0) Setting up resulting variables which are being read from L1BP
    //-----------------------------------------------------------------
    Uvec sab_counter("sab counter",Nrecord);
    Array1D<unsigned int> quality_flag("quality flag",Nrecord);
    Uvec pass_geom_time_offset("",Nrecord);
    Uvec transmit_time_offset("",Nrecord);
    Uvec ncnt_rl("ncnt_rl",Nrecord);
    Uvec ncnt_nd("ncnt_nd",Nrecord);
    Uvec ncnt_radio("ncnt_radio",Nrecord);
    Array1D<unsigned int>  beam_number("beam number",Nrecord);
    Array1D<Time> time("time",Nrecord);
    Array1D<Uvar> rlotmp("rlotmp",Nrecord);
    Array1D<Uvar> lnatmp("lna tmp",Nrecord);
    Array1D<Uvar> adctmp("adc tmp",Nrecord);
    vector<Uvec>  wgbtmp; wgbtmp.clear();//lower waveguide temperature
    Array1D<Uvar> fwdtmp("fwdtmp",Nrecord);
    vector<Uvec>  rfes_beam_tmp; rfes_beam_tmp.clear();//beam temp inside RFES_FEE
    Array1D<Uvar> scwgtmp("scwb tmp",Nrecord);
    Array1D<Uvar> hgatmp("hga tmp",Nrecord),feedtmp("",Nrecord);
  
    //---------------------------------------------------
    //(1) Extracted parameters from L1B
    // sab_counter, norm_cnt_rl, norm_cnt_nd, norm_cnt_radio
    // qaulity_flag, sclk, brst, beam_number, rlotmp, lnatmp
    // adctmp, wgb1t1, wgb3t1, wgb3t2, wgb3t3, wgb5t1
    //----------------------------------------------------
    unsigned int Nparam=27;
    double *return_doubles; 
    return_doubles =new double[Nparam] ;
    char* param_list[Nparam];
    param_list[0]="sab_counter";
    param_list[1]="norm_cnt_rl";
    param_list[2]="norm_cnt_nd";
    param_list[3]="norm_cnt_radio";
    param_list[4]="engineer_qual_flag";
    param_list[5]="sclk";
    param_list[6]="brst";
    param_list[7]="beam_number";
    param_list[8]="rlotmp";
    param_list[9]="lnatmp";
    param_list[10]="adctmp";
    param_list[11]="wgb1t1";
    param_list[12]="wgb3t1";
    param_list[13]="wgb3t2";
    param_list[14]="wgb3t3";
    param_list[15]="wgb5t1";
    param_list[16]="fwdtmp";
    param_list[17]="be1tmp";
    param_list[18]="be2tmp";
    param_list[19]="be3tmp";
    param_list[20]="be4tmp";
    param_list[21]="be5tmp";
    param_list[22]="scwg_tmp";
    param_list[23]="hga_tmp";
    param_list[24]="feed_tmp";
    param_list[25]="pass_geom_time_offset";
    param_list[26]="transmit_time_offset";

   
    unsigned int i_count = 0;
    double sclk_double, brst_double;
    Uvec wgb_dummy("dummy",5);
    Uvec beam_tmp_dummy("dummy",5);

     while(!data.eof()){
       if (!data.returnParams(return_doubles,param_list,
				 Nparam,range_flag,rangeLow,
				rangeHi))
	{
	  sab_counter(i_count)= round_double(*return_doubles);
	  ncnt_rl(i_count) = Uvar(*(return_doubles+1),"1/s");
	  ncnt_nd(i_count) =Uvar(*(return_doubles+2),"1/s");
	  ncnt_radio(i_count)=Uvar(*(return_doubles+3),"1/s");
	  quality_flag(i_count)= (unsigned int) round_double(*(return_doubles+4));
	  sclk_double=*(return_doubles+5);
	  brst_double=*(return_doubles+6);
	  beam_number(i_count)=(unsigned int) round_double(*(return_doubles+7));
	  time(i_count).setSclk("Cassini",(unsigned long) round_double(sclk_double));
	  time(i_count)+= Uvar(brst_double,"s");
	  rlotmp(i_count) = Uvar( *(return_doubles+8),"K");
	  lnatmp(i_count)=Uvar( *(return_doubles+9),"K");
	  adctmp(i_count)=Uvar( *(return_doubles+10),"K");
	  wgb_dummy=Uvar(0,"K");//reset to 0 K
	  for(unsigned int i=0;i<5;++i)
	    {
	      wgb_dummy(i)=Uvar( *(return_doubles+11+i),"K");
	    }
	  wgbtmp.push_back(wgb_dummy);
	  fwdtmp(i_count) = *(return_doubles+16);
	  
	  beam_tmp_dummy=Uvar(0,"K");
	  for(unsigned int i=0;i<5;++i)
	    {
	      beam_tmp_dummy(i)= *(return_doubles+17+i);
	    }
	  rfes_beam_tmp.push_back(beam_tmp_dummy);//save into vector array
	  scwgtmp(i_count) = Uvar(*(return_doubles+22),"K");
	  hgatmp(i_count) =Uvar(*(return_doubles+23),"K");
	  feedtmp(i_count)=Uvar( *(return_doubles+24),"K");
	  pass_geom_time_offset(i_count)=Uvar ( *(return_doubles+25),"s");
	  transmit_time_offset(i_count)=Uvar( *(return_doubles+26),"s");
	  time(i_count)+= transmit_time_offset(i_count);
	  time(i_count) += pass_geom_time_offset(i_count);
	  i_count++;
	}
     }
     if(i_count !=Nrecord) ErrorMessage("rad_proc::expected number of records is not read in").throwMe();
     cout<<"Finished extracting radiometer relevant paramenters from L1BP"<<endl;     
     
     //-----------------------------------
     //plot data :normialized RL , ND, and Radiometer
     //---------------------------------
     if (target!=""){
       Plot a;
       a.addXY(sab_counter,"",ncnt_rl,"1/s",line("solid","red",1),sym("none"));
       a.setTitle("normalized RL count vs sab counter");
       a.setXlabelandUnits("sab number");
       a.setYlabelandUnits("normalized RL count");
       a.setFile(prefix_str + "norm_RL.eps");
       a.setCmdFile(prefix_str + "norm_RL.xmgr");
       a.show(plot_target);
       
       Plot b;
       b.addXY(sab_counter,"",ncnt_nd,"1/s",line("solid","red",1),sym("none"));
       b.setTitle("normalized ND count vs sab counter");
       b.setXlabelandUnits("sab number");
       b.setYlabelandUnits("normalized ND count");
       b.setFile(prefix_str + "norm_ND.eps");
       b.setCmdFile(prefix_str + "norm_ND.xmgr");
       b.show(plot_target);
       
       Plot c;
       c.addXY(sab_counter,"",ncnt_radio,"1/s",line("solid","red",1),sym("none"));
       c.setTitle("normalized Radiometer count vs sab counter");
       c.setXlabelandUnits("sab number");
       c.setYlabelandUnits("normalized radio count");
       c.setFile(prefix_str + "norm_Radio.eps");
       c.setCmdFile(prefix_str + "norm_Radio.xmgr");
       c.show(plot_target);


       Plot e;
       e.addXY(sab_counter,"",feedtmp,"K",line("solid","black",1),sym("none"));
       e.addLegend("feed tmp");
       e.addXY(sab_counter,"",scwgtmp,"K",line("solid","red",1),sym("none"));
       e.addLegend("swg tmp");
       e.addXY(sab_counter,"",hgatmp,"K",line("solid","green",1),sym("none"));
       e.addLegend("hga tmp");
       e.setTitle("feed swg hga  temperature vs sab ");
       e.setXlabelandUnits("sab number");
       e.setYlabelandUnits("temperature");
       e.setFile("feed_swg_hga.eps");
       e.setCmdFile("feed_swg_hga.xmgr");
       e.show(plot_target);
     }
     
    
     
     //----------------------------------------------------------------
     //(3) compute the angle  between limb direction and boresight
     // direction for solar targets
     // in units of beam oneway 3dB
     // positive: boresight off the target
     // negative: boresight on the target
     //----------------------------------------------------------------
     Array2D<double> solar_target_distances("min target distances",
					    Nrecord,
					    Ntarget);//records * targets

     Umat azim_distance("azim distance", Nrecord, Ntarget);//records * targets
     Umat elev_distance("elev distance",Nrecord,Ntarget);//records * targets
     solar_target_distances=0.0;
     azim_distance=Uvar(0,"rad");
     elev_distance=Uvar(0,"rad");

     //recycled variables for speed up
     DirectionVector boresight,target_dir;
     PositionVector target_position;
     Uvar distance, limb_angle, boresight_angle;
     double dummy_a, dummy_b;
     for(unsigned int i_rec=0;i_rec<Nrecord;i_rec++){
       //------------------------------
       //Is this a valid  record?
       //-------------------------------
       if(bitget(quality_flag(i_rec),0,0)==1 ||
	  bitget(quality_flag(i_rec),1,1)==1 ||
	  beam_number(i_rec)==0) continue;
      
       //--------------------
       //boresight direction
       //---------------------
       boresight=DirectionVector("bore",Fvector[beam_number(i_rec)-1],
				 time(i_rec),0,0,1);

       //-------------------------------------------------
       //for each target, compute limb angle and beam angle
       //--------------------------------------------------
       for(unsigned int i_target=0;i_target<Ntarget;++i_target){
	 target_position=PositionVector("",target_frame(i_target),
					time(i_rec),0,0,0).representIn(j2000c);
	 target_dir = target_position;
	 target_dir.representIn(Fvector[beam_number(i_rec)-1]);
	 target_dir.getAzimuthElevation(azim_distance(i_rec,i_target),
					elev_distance(i_rec,i_target));
	 distance = target_position.magnitude();
	 limb_angle =asin(target_radius(i_target)/distance);
	 
	 //-------------------------------------------------------
	 //special care when boresight is lined up with target direction
	 //-----------------------------------------------------------
	 double bore_target=dot(boresight, target_dir);
	 if(bore_target>=1.0) boresight_angle=Uvar(0,"rad");
	 else
	   boresight_angle = Uvar(acos(bore_target),"rad");

	 dummy_a =((boresight_angle-limb_angle)
		   /elev_beam3dB_oneway(beam_number(i_rec)-1)).getInUnits("");
	 dummy_b=((boresight_angle-limb_angle)
		  /azi_beam3dB_oneway(beam_number(i_rec)-1)).getInUnits("");

	 //-----------------------------------
	 //take smaller value of beam distances
	 //------------------------------------
	 solar_target_distances(i_rec,i_target)=dummy_a;
	 if(dummy_b < solar_target_distances(i_rec,i_target))
	   solar_target_distances(i_rec,i_target)=dummy_b;
	 
	 //-----------------------------------------
	 //when something is not going well(level==0)
	 //------------------------------------------
	 if(dbg.level==0)
	   cout<<"target distance "<<target_name(i_target)<<" "<< solar_target_distances(i_rec,i_target)<<endl;
       }//loop over targets
     }//end of all records
     cout<<"Finish calculating the beam angles of solar targets "<<endl;
       
     //----------------------------------------------
     //debug::plot beam distances for various targets
     //---------------------------------------------
     if(dbg.level==1){
       Plot a;
       Uvec distance_to_target(" ",Nrecord);
       for(unsigned int i=0;i<Ntarget;++i)
	 {
	   distance_to_target=0.0;
	   for(unsigned int j=0;j<Nrecord;++j)
	     distance_to_target(j)=solar_target_distances(j,i);
	   a.addXY(sab_counter,"",distance_to_target,line("solid",color_table(i),1),sym("none"));
	   a.addLegend(target_name(i));//add target names
	 }
       
       a.setTitle(" beam distance to specified targets in config");
       a.setFile(prefix_str +"solar_targets.eps");
       a.setCmdFile(prefix_str+ "solar_targets.xmgr");
       a.show(plot_target);
       
       Plot b;
       b.setTitle("target in beam frame:r-azim, b-elev");
       //for(unsigned int i=0;i<Ntarget;++i)
       //{
       b.addXY(sab_counter,"",azim_distance.getCol(0),"deg",line("none"),sym(color_table(0),"circle",1));
       b.addXY(sab_counter,"",elev_distance.getCol(0),"deg",line("none"),sym(color_table(0),"circle",1));
       // }
       b.setFile(prefix_str +"target_angle.eps");
       b.setCmdFile(prefix_str+ "target_angle.xmgr");
       b.setXlabelandUnits("sab number");
       b.setYlabelandUnits("deg");
       b.setFile(prefix_str +"target_angle.eps");
       b.setCmdFile(prefix_str+ "target_angle.xmgr");
       b.show(plot_target);
     }


     //-------------------------------------------------------------------
     //(4) determine the  nearest target among solar targets
     //----------------------------------------------------------------------
     Array1D<int>  target_record("sab target",Nrecord);//each beam's target
     Array1D<double> smallest_beam_distance("target whose beam distance is smallest",Nrecord);
     Array1D<Uvar> nearest_target_temperature("nearest target temperature",Nrecord);
     
     //----------------------------------------------
     //assume default target, each time it computes
     //target id which has the nearest range
     //-----------------------------------------------
     for (unsigned int i_rec=0;i_rec<Nrecord;i_rec++)
       {
	 //select min target distance among target names
	 // variable: solar_target_distances
	 unsigned int index_min_distance = target_id(0);
	 double min_distance = solar_target_distances(i_rec,0);
	 Uvar target_tmp=target_temperature(0);
	 //for solar targets
	 for (unsigned int i_target=0;i_target<Ntarget;++i_target)
	   {
	     if(min_distance>solar_target_distances(i_rec,i_target))
	       {
		 index_min_distance = target_id(i_target);
		 min_distance = solar_target_distances(i_rec,i_target);
		 target_tmp=target_temperature(i_target);
	       }
	   }
	 target_record(i_rec)=index_min_distance;
	 smallest_beam_distance(i_rec)= min_distance;
	 nearest_target_temperature(i_rec)=target_tmp;
	 
       }//loop over record
     
     //-----------------------------------------
     //debug: plot min distance
     //-----------------------------------------
     if(dbg.level==1)
       {
         Plot a;
         Uvec y("",smallest_beam_distance.size());
         for(unsigned int i=0;i<y.size();++i)y(i)=smallest_beam_distance(i);
         a.addY(y,line("solid","red",1),sym("none"));
         a.setTitle("min  beam distance to solar targets");
	 a.setFile(prefix_str +"min_beam_dist.eps");
	 a.setCmdFile(prefix_str+ "min_beam_dist.xmgr");	 
	 a.show(plot_target);
       }       

    

    

     //----------------------------------
     //MAIN:   process each beam separately
     //----------------------------------
     unsigned int i_counter=0;
     string output_filename=prefix_str+".txt";
     ofstream output_file(output_filename.c_str());
     if(!output_file.is_open()) ErrorMessage("Could not create a file: output_filename").throwMe();

     //----------------
     //Declare container
     //---------------
     map<Uvar, Uvar> sab_vs_antenna_temperature;
     map<Uvar, Uvar> sab_vs_system_gain;
     map<Uvar, Uvar> sab_vs_receiver_temperature;
     map<Uvar, Uvar> sab_vs_antenna_temperature_stddev;


     for (unsigned int i_beam=0; i_beam < 5;++i_beam)
       {
	 unsigned int beam_id = i_beam + 1;
	 string prefix_str1=prefix_str+"beam"+toStr(beam_id)+"_";
	 cout<<"Processing beam: "<<beam_id<<endl;
	 
	 //------------------------------
	 //collect beam_id specific data 
	 //-------------------------------
	 vector<Uvar> beam_sab_counter;
	 vector<Uvar> beam_rel_time;
	 vector<Uvar> beam_ncnt_rl;
	 vector<Uvar> beam_ncnt_nd;
	 vector<Uvar> beam_ncnt_radio;
	 vector<Uvar> beam_distance_from_limb_to_boresight;
	 vector<Uvar> beam_rlotmp;
	 vector<Uvar> beam_lnatmp;
	 vector<Uvar> beam_adctmp;
	 vector<Uvar> beam_wgbtmp;
	 vector<int> beam_target;
	 vector<Uvar> beam_fwdtmp;
	 vector<Uvar> beam_rfes_tmp;
	 vector<Uvar> beam_feed_tmp;
	 vector<Uvar> beam_scwg_tmp;
	 vector<Uvar> beam_hga_tmp;
	 vector<Uvar> beam_target_temperature;
	 //---------------------
	 //clear the container
	 //----------------------
	 beam_sab_counter.clear();
	 beam_rel_time.clear();
	 beam_ncnt_rl.clear();
	 beam_ncnt_nd.clear();
	 beam_ncnt_radio.clear();
	 beam_distance_from_limb_to_boresight.clear();
	 beam_rlotmp.clear();
	 beam_lnatmp.clear();
	 beam_adctmp.clear();
	 beam_wgbtmp.clear();
	 beam_target.clear();
	 beam_fwdtmp.clear();
	 beam_rfes_tmp.clear();
	 beam_feed_tmp.clear();
	 beam_hga_tmp.clear();
	 beam_scwg_tmp.clear();
	
	 
	 
	 Uvec dummy("dummy",5), dummy2("",5); 
	 unsigned int bit_reading;
	 for (unsigned int i_rec=0;i_rec<Nrecord;i_rec++)
	   {
	     //----------------------------------------------
	     //check the validity of each record
	     //select data specific to the designated beam
	     //-----------------------------------------------
	     bit_reading = quality_flag(i_rec);
	     if(bitget(bit_reading,0,0)==1 ||
		bitget(bit_reading,1,1)== 1 ||
		beam_number(i_rec) != beam_id) continue;

	     beam_sab_counter.push_back(sab_counter(i_rec));
	     beam_rel_time.push_back(time(i_rec)-time(0));
	     beam_ncnt_rl.push_back(ncnt_rl(i_rec));
	     beam_ncnt_nd.push_back(ncnt_nd(i_rec));
	     beam_ncnt_radio.push_back(ncnt_radio(i_rec));
	     beam_rlotmp.push_back(rlotmp(i_rec));
	     beam_lnatmp.push_back(lnatmp(i_rec));
	     beam_adctmp.push_back(adctmp(i_rec));
	     beam_target.push_back(target_record(i_rec));//target name(integer index)
	     beam_fwdtmp.push_back(fwdtmp(i_rec));
	     dummy=rfes_beam_tmp[i_rec];
	     beam_rfes_tmp.push_back(dummy(i_beam));
	     dummy2 = wgbtmp[i_rec];
	     if(beam_id ==1 || beam_id==2) beam_wgbtmp.push_back(dummy2(0));
	     if(beam_id==3) beam_wgbtmp.push_back(dummy2(2));
	     if(beam_id==4 || beam_id==5) beam_wgbtmp.push_back(dummy2(4));
	     beam_distance_from_limb_to_boresight.push_back(smallest_beam_distance(i_rec)); 
	     beam_feed_tmp.push_back(feedtmp(i_rec));
	     beam_hga_tmp.push_back(hgatmp(i_rec));
	     beam_scwg_tmp.push_back(scwgtmp(i_rec));
	     beam_target_temperature.push_back(nearest_target_temperature(i_rec));
	   }
	 
	 //------------------------------
	 //check the size of containter
	 //--------------------------------
	 unsigned int Nbeam_record=beam_sab_counter.size();
	 cout<<Nbeam_record<<"  samples were chosen for beam number "
	     <<beam_id<<endl;
	 if(beam_sab_counter.size()==0) continue;//go over next beam

	 //------------------------------------------------
	 //make sure all the container have the same size
	 //--------------------------------------------------
	 if(Nbeam_record!=beam_rel_time.size() ||
	    Nbeam_record!=beam_ncnt_rl.size()  ||
	    Nbeam_record!=beam_ncnt_nd.size()  ||
	    Nbeam_record!=beam_ncnt_radio.size() ||
	    Nbeam_record!=beam_distance_from_limb_to_boresight.size() ||
	    Nbeam_record!=beam_rlotmp.size()   ||
	    Nbeam_record != beam_lnatmp.size() ||
	    Nbeam_record != beam_adctmp.size() ||
	    Nbeam_record != beam_target.size() ||
	    Nbeam_record != beam_wgbtmp.size() ||
	    Nbeam_record != beam_hga_tmp.size()||
	    Nbeam_record != beam_target_temperature.size())
	   
	   
	   {
	     cout<<beam_rel_time.size()<<" "<<beam_ncnt_rl.size()<<" "<<beam_ncnt_nd.size()<<
	       " "<<beam_ncnt_radio.size()<<" "<<beam_distance_from_limb_to_boresight.size()<< " " 
		 <<beam_rlotmp.size()<<" "<< beam_hga_tmp.size()<<" "<<beam_target_temperature.size()<<endl;
	     ErrorMessage("Container size mismatch").throwMe();
	   }
	 
	 //------------------------------------------------
	 //plots of each beam's radio count, RL count, and beam distance
	 //---------------------------------------------
	 if(dbg.level==1)
	   {
	     //---------------
	     //normalized count
	     //---------------
	     Plot a;
	     a.addXY(beam_sab_counter,"",beam_ncnt_radio,"1/s",line("solid","red",1),sym("none"));
	     a.setTitle("normalized count for beam vs sab counter:beam" +toStr(beam_id));
	     a.setXlabelandUnits("sab number");
	     a.setYlabelandUnits(" norm_Radio");
	     a.setFile(prefix_str1 +"norm_Radio.eps");
	     a.setCmdFile(prefix_str1+ "norm_Radio.xmgr");
	     a.show(plot_target);
	     
	     //--------------------
	     //lna adc temperature 
	     //--------------------     
	     Plot b;
	     b.addXY(beam_sab_counter,"",beam_rlotmp,"K",line("solid","black",1),sym("none"));
	     b.addLegend("rlo temp");
	     b.addXY(beam_sab_counter,"",beam_lnatmp,"K",line("solid","red",1),sym("none"));
	     b.addLegend("lna temp");
	     b.addXY(beam_sab_counter,"",beam_adctmp,"K",line("solid","blue",1),sym("none"));
	     b.addLegend("adc temp");
	     b.setTitle("RL(BL), LNATMP(R) , adc(B),waveguide temperature vs sab counter:beam" +toStr(beam_id));
	     b.setXlabelandUnits("sab number");
	     b.setYlabelandUnits("temperature");
	     b.setFile(prefix_str1 +"RL_LNA_ADC_temp.eps");
	     b.setCmdFile(prefix_str1 +"RL_LNA_ADC_temp.xmgr");
	     b.show(plot_target);
	     
	     //----------------------
	     //waveguid temperature
	     //---------------------
	     Plot b1;
	     b1.addXY(beam_sab_counter,"",beam_wgbtmp,"K",line("dash","red",1),sym("none"));
	     b1.setTitle("Waveguid temp:");
	     b1.setXlabelandUnits("sab number");
	     b1.setYlabelandUnits("beam waveguide temperature");
	     b1.setFile(prefix_str1 +"WGBT.eps");
	     b1.setCmdFile(prefix_str1 +"WGBT.xmgr");
	     b1.show(plot_target);
	     
	     //---------
	     //hgatmp
	     //-----------
	     Plot b2;
	     b2.addXY(beam_sab_counter,"",beam_hga_tmp,"K",line("dash","red",1),sym("none"));
	     b2.setTitle("HGA temp");
	     b2.setXlabelandUnits("sab number");
	     b2.setYlabelandUnits("beam hga temperature");
	     b2.setFile(prefix_str1 +"hga.eps");
	     b2.setCmdFile(prefix_str1 +"hga.xmgr");
	     b2.show(plot_target);
	     
	     //----------------
	     //normal RL       
	     //-----------------    
	     Plot c;
	     c.addXY(beam_sab_counter,"",beam_ncnt_rl,"1/s",line("solid","red",1),sym("none"));
	     c.setTitle("normalized RL count");
	     c.setXlabelandUnits("sab number");
	     c.setYlabelandUnits(" norm_RL");
	     c.setFile(prefix_str1 +"norm_RL.eps");
	     c.setCmdFile(prefix_str1 +"norm_RL.xmgr");      
	     c.show(plot_target); 
	     
	     //-----------
	     //distance
	     //--------------
	     Plot d;
	     d.addXY(beam_sab_counter,"", beam_distance_from_limb_to_boresight," ",line("solid","red",1),sym("none"));
	     d.setTitle("beam distance in units of elev 3dB oneway:beam"+toStr(beam_id));
	     d.setXlabelandUnits("sab number");
	     d.setYlabelandUnits(" distance");
	     d.setFile(prefix_str1 +"distance.eps");
	     d.setCmdFile(prefix_str1 +"distance.xmgr");
	     d.show(plot_target); 
	   }
	 
  	 


	 //------------------------
	 //Compute time-averaged beam_ncnt_rl
	 //------------------------
	 cout<<"Compute time-averaged RL counts "<<endl;
	 vector<Uvar> avg_beam_ncnt_rl;
	 vector<Uvar> beam_relative_gain;

	 beam_relative_gain.clear();
	 avg_beam_ncnt_rl.clear();
	 compute_time_averaged_rl(avg_time, beam_ncnt_rl,beam_rel_time, avg_beam_ncnt_rl);
	 if(avg_beam_ncnt_rl.size()!=Nbeam_record){
	   cout<<"Input rl size "<< beam_ncnt_rl.size()<<endl;
	   cout<<"Output avg rl size "<< avg_beam_ncnt_rl.size()<<endl;
	   ErrorMessage("input rl and output avg rl size mismatch").throwMe();
	 }

	 //-----------------------------------
	 //compute relative gain
	 //-----------------------------------
	 cout<<"Compute relative gain "<<endl;
	 for (unsigned int i_rec=0;i_rec<Nbeam_record;i_rec++){
	   beam_relative_gain.push_back(norm_gain*avg_beam_ncnt_rl[i_rec]/(beam_rlotmp[i_rec]+Trec));
	 }
 
	 //--------------
	 //Must be shown
	 //-------------
	 cout<<"beam_relative gain "<< beam_relative_gain[0]<<endl;
	 Plot Rel_Gain;
	 Rel_Gain.addXY( beam_rel_time,"min",beam_relative_gain," ",line("solid","red",1),sym("none"));
	 Rel_Gain.setTitle("Relative gain  close to 1 ");
	 Rel_Gain.setXlabelandUnits("time in min");
	 Rel_Gain.setYlabelandUnits("relative gain");
	 Rel_Gain.setFile(prefix_str1 +"rel_gain.eps");
	 Rel_Gain.setCmdFile(prefix_str1 +"rel_gain.xmgr");
	 Rel_Gain.show(plot_target); 

	
	   
	 //--------------------
	 //Uncalibrated Ta
	 //--------------------
	 vector<Uvar> uncal_Ta;
	 uncal_Ta.clear();
	 for (unsigned int i_rec=0;i_rec<Nbeam_record;i_rec++){
	   uncal_Ta.push_back(radiometer_cal*(beam_ncnt_radio[i_rec]-avg_beam_ncnt_rl[i_rec])*Troom/((beam_rlotmp[i_rec]-Tcoldsky)*beam_relative_gain[i_rec]));
	 }
	 //-------------------
	 //uncal Ta
	 //---------------------
	 Plot Uncal_Ta;
	 Uncal_Ta.addXY(beam_rel_time,"min", uncal_Ta,"K",line("solid","red",1),sym("none"));
	 Uncal_Ta.setTitle("Ta without cold sky reference");
	 Uncal_Ta.setXlabelandUnits("time in min");
	 Uncal_Ta.setYlabelandUnits("Ta before cd ref");
	 Uncal_Ta.setFile(prefix_str1 +"uncal_Ta.eps");
	 Uncal_Ta.setCmdFile(prefix_str1 +"uncal_Ta.xmgr");
	 Uncal_Ta.show(plot_target); 


	 //---------------------------------------------------
	 //determine cold sky pointing  and target pointing
	 //-------------------------------------------------
	 vector<bool> is_beam_coldsky;
	 is_beam_coldsky.clear();
	 for(unsigned int i_rec=0;i_rec<Nbeam_record;++i_rec)
	   {
	     if(beam_distance_from_limb_to_boresight[i_rec]<min_beam_distance)
	       {
		 is_beam_coldsky.push_back(false);
	       }
	     else
	       {
		 is_beam_coldsky.push_back(true);
	       }       
	   }
	 if(is_beam_coldsky.size()!=Nbeam_record)
	   ErrorMessage("rad_proc:: size of is_beam_coldsky does not match to total number of record, Nbeam_record ").throwMe();
	 

	 //------------------------------------------
	 //compute calibration start and end points
	 //-------------------------------------------
	 vector<Uvar> zero_Ta;
	 vector<Uvar> zero_Ta_cal_time;
	 vector<Uvar> cal_start, cal_end, Time_cal_start, Time_cal_end, Sab_cal_start,Sab_cal_end;
	 for(unsigned int i_rec=0;i_rec<Nbeam_record-1;i_rec++)
	   {
	     if(!is_beam_coldsky[i_rec] && is_beam_coldsky[i_rec+1])
	       {
		 cal_start.push_back(i_rec+1);  
		 Time_cal_start.push_back(beam_rel_time[i_rec+1]);
		 Sab_cal_start.push_back(beam_sab_counter[i_rec+1]);
	       }
	     if(is_beam_coldsky[i_rec] && !is_beam_coldsky[i_rec+1]) 
	       {
		 cal_end.push_back(i_rec);
		 Time_cal_end.push_back(beam_rel_time[i_rec]);  
		 Sab_cal_end.push_back(beam_sab_counter[i_rec]);
	       }
	   } 
	 
	 unsigned int Ncal_start = cal_start.size();
	 cout<<"Cal start points :  "<< Ncal_start<<endl;
	 if(Ncal_start != Time_cal_start.size()||
	    Ncal_start != Sab_cal_start.size()) { 	 
	   cout<<"sizes "<<Time_cal_start.size()
	       <<" "<<Sab_cal_start.size()<<" "<<endl;
	   ErrorMessage("Ncal start container size mismatch").throwMe();
	 }
	 
	 unsigned int Ncal_end = cal_end.size();
	 cout<<"Cal end points:  "<< Ncal_end<<endl;
	 if(Ncal_end != Time_cal_end.size() ||
	    Ncal_end != Sab_cal_end.size()) {
	   cout<<"sizes "<<" "<<Time_cal_end.size()
	       <<" "<<Sab_cal_end.size()<<" "<<endl;
	   ErrorMessage("Ncal end container size mismatch");
	 }
	 
	 if(Ncal_start ==0 && Ncal_end==0) continue;//go to next beam
	 //----------------------------
	 //Find cold sky reference: farthest off-titan points for reference
	 //---------------------------	

	 //-------------------------------------------------------------------
	 //collect cold sky reference points before first "Scan End"
	 //----------------------------------------------------------------------
	 while(Sab_cal_end.size()!=0){
	   unsigned int i_rec_start= (unsigned int) cal_end[0].getInUnits("");//get index: first "scan end"
	   unsigned int i_begin=0;//set to a record 1 min before the i_rec_start
	   for(unsigned int i=i_rec_start;i!=0;i--){//search backward
	     if((beam_rel_time[i_rec_start]-beam_rel_time[i])>Uvar(60,"s")){
	       i_begin=i;
	       break;
	     }
	   }
	   //scan up
	   for(unsigned int i=i_begin;i<i_rec_start;++i){//search backward
	     if(is_beam_coldsky[i]){ 
	       //valid record for cold sky calibration
	       zero_Ta.push_back(uncal_Ta[i]);
	       zero_Ta_cal_time.push_back(beam_rel_time[i]);
	     }
	   }
	   break;//collect cold sky reference points before first "Scan end"
	 }

	 //--------------------------------------------------------------------------
	 //Collect cold sky reference between first "Scan end" and last "Scan start"
	 //-----------------------------------------------------------------------------
	 for(unsigned int i_s=0;i_s<Ncal_start-1;++i_s){
	   //---------------
	   //Find scan end: except first and last segment
	   // segment start with : scan_start --- scan_end
	   // Find scan end of the current segment starting from scan start
	   //---------------
	   unsigned int i_rec_start= (unsigned int) cal_start[i_s].getInUnits("");//get index
	   unsigned int i_rec_end=0;
	   bool found_sab_end=false;
	   for(unsigned int i_e=0;i_e<Ncal_end;++i_e){
	     
	     if(Sab_cal_end[i_e]> Sab_cal_start[i_s] &&  Sab_cal_end[i_e]< Sab_cal_start[i_s+1]){
	       i_rec_end= (unsigned int) cal_end[i_e].getInUnits("");//get index
	       found_sab_end=true;
	       break;
	     }
	   }
	   //---------------------------------------------------------
	   //Check the existence of sab_end: Sab end needs to exist
	   //----------------------------------------------------------
	   if(!found_sab_end){
	     cout<<"Could not find end of the current scan: scan number "<<Sab_cal_start[i_s]<<" "<<Sab_cal_start[i_s+1]<<endl;
	     for(unsigned int i_e=0;i_e<Ncal_end;++i_e) cout<<Sab_cal_end[i_e]<<" ";
	     cout<<endl;
	     continue;//go to next one
	   }
	   else{
	     if(i_rec_end< i_rec_start) ErrorMessage("Scan end occurs earlier than scan start").throwMe();
	   }  
	   
	   //------------------------
	   //When scan end is detected
	   // find farthest point
	   // record number:i_s, i_e
	   //-----------------------
	   cout<<"Scan start and end record id "<< i_rec_start<<" "<<i_rec_end<<endl;
	   vector<Uvar> scan_Ta;//Ta while pointing cold sky
	   vector<Uvar> scan_beam_distance;//beam distance while pointing cold sky
	   vector<Uvar> scan_time;//scan time
	   vector<unsigned int> farthest_distance_index;
	   
	   for(unsigned int i=i_rec_start;i<i_rec_end;++i){
	     scan_Ta.push_back(uncal_Ta[i]);
	     scan_beam_distance.push_back(beam_distance_from_limb_to_boresight[i]);
	     scan_time.push_back(beam_rel_time[i]);
	   }
	   Uvar max_distance, min_distance;
	   min_max(min_distance,max_distance,scan_beam_distance);//assign min and max beam distance
	   find_greater_than_target_value(scan_beam_distance, (min_distance+max_distance)/2, farthest_distance_index);//find index
	   if(farthest_distance_index.size() ==0 && scan_beam_distance.size()!=0)
	     ErrorMessage("Failed to find an index value corresponding max value in the vector array").throwMe();

	  
	   //cout<<"farthest distance index size "<< farthest_distance_index.size()<<endl;
	   for(unsigned int i=0;i<farthest_distance_index.size();++i){
	     zero_Ta.push_back(scan_Ta[farthest_distance_index[i]]);//uncal Ta when farthest away from Target
	     zero_Ta_cal_time.push_back(scan_time[farthest_distance_index[i]]);// time when farthest away from Target
	   }
	 }//loop over scan start

	 //---------------------------------------
	 //Collect cold sky reference After last "Scan start"
	 //----------------------------------------
	 while(Sab_cal_start.size()!=0){
	   unsigned int i_rec_start= (unsigned int) cal_start[cal_start.size()-1].getInUnits("");//get index: last "scan start"
	   for(unsigned int i=i_rec_start;i<Nbeam_record;++i){
	     if(is_beam_coldsky[i] &&(beam_rel_time[i]-beam_rel_time[i_rec_start])<Uvar(60,"s")){
	       //valid record for cold sky calibration
	       zero_Ta.push_back(uncal_Ta[i]);
	       zero_Ta_cal_time.push_back(beam_rel_time[i]);
	     }
	   }
	   break;//collect cold sky reference points after last "scan start"
	 }




	 //-----------------------------------------------------
	 //Apply median filter within 1 min of calibration time
	 //------------------------------------------------------
	 cout<<"apply median filter to an array whose size is "<< zero_Ta.size()<<endl;
	 vector<Uvar> median_zero_Ta;
	 median_zero_Ta.resize(zero_Ta.size());//make the same size vector
	 
	 for(unsigned int i=0;i<zero_Ta.size();++i){
	   //cout<<"i: "<< i<<endl;
	   
	   vector<Uvar> tmp_ta;
	   
	   unsigned int j=i;
	   //backward search
	   while(j >=0){
	     if((zero_Ta_cal_time[i]-zero_Ta_cal_time[j])<Uvar(1*60,"s"))
	       tmp_ta.push_back(zero_Ta[j]);
	     
	     if((zero_Ta_cal_time[i]-zero_Ta_cal_time[j])>Uvar(1*60,"s"))
	       break;

	     if(j==0) break;
	     else j--;
	   }

	   //foward search
	   j=i+1;
	   while(j<zero_Ta.size()){
	     if((zero_Ta_cal_time[j]-zero_Ta_cal_time[i])<Uvar(1*60,"s"))
	       tmp_ta.push_back(zero_Ta[j]);
	     
	     if((zero_Ta_cal_time[j]-zero_Ta_cal_time[i])>Uvar(1*60,"s"))
	       break;

	     if(j==(zero_Ta.size()-1)) break;
	     else j++;
	   }

	   //for(unsigned int j=0;j<zero_Ta.size();++j){
	   //if(fabs((zero_Ta_cal_time[i]-zero_Ta_cal_time[j])) > Uvar(1*60,"s")) continue;
	   //if((zero_Ta_cal_time[j]-zero_Ta_cal_time[i]) > Uvar(1*60,"s")) break;
	   //
	   //if(fabs(zero_Ta_cal_time[i]-zero_Ta_cal_time[j])<Uvar(1*60,"s")){
	   //  tmp_ta.push_back(zero_Ta[j]);
	   //}
	   //}

	   //---------------------------
	   //Sorting for median filter
	   //--------------------------
	   if(tmp_ta.size()==1) {
	     median_zero_Ta[i]= zero_Ta[i];
	     continue;//skip this record
	   }
	   sort(tmp_ta.begin(),tmp_ta.end());
	   median_zero_Ta[i]= tmp_ta[(unsigned int)tmp_ta.size()/2];
	 }
	 

	 //------------------------------------------------------
	 //Need to add  effect of target with a finite temperature
	 // reference point is set to 0, although there is a finite
	 // contribution fromt target
	 // time: ero_Ta_cal_time=beam_relative_time
	 // absolute time= time(0) + zero_Ta_cal_time
	 // zero_Ta = zero_Ta - Tcold_sky -Target(T=T-Tcoldsky)
	 // target_name,Ntarget, target_temperature, target_frame
	 //---------------------------------------------------
	 
	 cout<<"compute Beam integration over the sample "<<endl;
	 cout<<"cal points "<< median_zero_Ta.size()<<endl;
	 vector<Uvar> beam_integration_power;
	 beam_integration_power.resize(median_zero_Ta.size());

	 //for each calibration point looking at cold sky
	 for(unsigned int i_cal=0;i_cal<median_zero_Ta.size();++i_cal){
	   Time t= time(0);
	   t += zero_Ta_cal_time[i_cal];
	   Uvar sum("",0.0,"K");

	   if(i_cal != 0){
	     //compute beam integration every 30 second
	     if( (zero_Ta_cal_time[i_cal]- zero_Ta_cal_time[i_cal-1]) <  Uvar(30,"s")){
	       beam_integration_power[i_cal]= beam_integration_power[i_cal-1];
	       median_zero_Ta[i_cal]= median_zero_Ta[i_cal-1];
	       continue;//go to next cal point
	     }
	   }

	   //for each target: mostly 1 target
	   for(unsigned int j_target=0;j_target<Ntarget;++j_target){
	     
	     //IAU frame
	     Frame ftarget("IAU_"+target_name(j_target), target_name(j_target));//target
	     StateVector sc;
	     ftarget.ephemeris(sc,"Cassini",t,"NONE");//state
	     
	     double distance=sc.position().magnitude().getInUnits("km");//direction to position
	     double r_a=target_radius(j_target).getInUnits("km");
	     double max_theta= asin(r_a/distance);//in beam frame
	     max_theta= max_theta *r2d;// in degrees


	     //beam frame
	     Frame fbeam=Fvector[i_beam];
	     DirectionVector bore("",fbeam,t,0,0,1);
	     //construct a new z-axis
	     DirectionVector beam_z= -sc.position();//from spacecraft to target
	     beam_z.representIn(fbeam);//from sc to target in beam frame

	     //construct beam_y using beam_z and beam's x-axis
	     DirectionVector look_x=DirectionVector("",fbeam,t,1,0,0);//beam x
	     DirectionVector beam_y=cross(beam_z, look_x);

	     //construct beam_x using beam_y and beam_z
	     DirectionVector beam_x=cross(beam_y,beam_z);


	     //new beam frame
	     Frame newbeam(sc.position(),beam_x,beam_y,beam_z);//NO need,  just need TR matrix
	     DirectionVector new_bore("",newbeam,t,0,0,1);
	     cout<<"angle between bore and target center "<< new_bore.angle(bore).getInUnits("rad")*r2d<<endl;
	     

	     //represent new xyz in old beam frame
	     beam_x.representIn(fbeam);
	     beam_y.representIn(fbeam);
	     beam_z.representIn(fbeam);
	     
	     //debug
	     //cout<<"new beam x "<< beam_x<<endl;
	     //cout<<"new beam y "<< beam_y<<endl;
	     //cout<<"new beam z "<< beam_z<<endl;

	     //matrix from new beam to old beam
	     double TR_from_newbeam_to_oldbeam[3][3];

	     TR_from_newbeam_to_oldbeam[0][0]=beam_x[DirectionVector::X];
	     TR_from_newbeam_to_oldbeam[1][0]=beam_x[DirectionVector::Y];
	     TR_from_newbeam_to_oldbeam[2][0]=beam_x[DirectionVector::Z];
	     
	     TR_from_newbeam_to_oldbeam[0][1]=beam_y[DirectionVector::X];
	     TR_from_newbeam_to_oldbeam[1][1]=beam_y[DirectionVector::Y];
	     TR_from_newbeam_to_oldbeam[2][1]=beam_y[DirectionVector::Z];
	     
	     TR_from_newbeam_to_oldbeam[0][2]=beam_z[DirectionVector::X];
	     TR_from_newbeam_to_oldbeam[1][2]=beam_z[DirectionVector::Y];
	     TR_from_newbeam_to_oldbeam[2][2]=beam_z[DirectionVector::Z];


	     unsigned int Num_theta_step= 1000;
	     double theta_start=0.0;//in deg
	     double theta_end= max_theta;// in deg
	     double theta_step= (theta_end - theta_start)/float(Num_theta_step);// 1000 elev steps

	     unsigned int Num_phi_step=360;
	     double phi_start=-180.0;
	     double phi_end=180.0;
	     double phi_step= (phi_end-phi_start)/float(Num_phi_step);//assume 1 deg in phi angle
	     

	     double theta, phi;
	     double look_new_beam[3];
	     double look_old_beam[3];
	     double elev, azi,gain, solid_angle;
	     double solid_angle_integral;
	     Uvar beam_integral;
	     solid_angle_integral=0.0;
	     beam_integral=Uvar(0.0,"K");
 
	     for(unsigned int j=0;j<Num_theta_step;++j)
	       for(unsigned int k=0;k<Num_phi_step;++k){
		 theta= theta_start + float(j)*theta_step;//deg
		 phi= phi_start + float(k)*phi_step;//deg
		 
		 theta= theta * d2r;//radian
		 phi= phi*d2r;//radian


		 //solid angle
		 solid_angle= sin(theta)*(theta_step*d2r)*(phi_step*d2r);
		 solid_angle_integral += solid_angle;

		 //look direction in target frame      
		 look_new_beam[0]= sin(theta)*cos(phi);//ftarget x
		 look_new_beam[1]= sin(theta)*sin(phi);//ftarget y
		 look_new_beam[2]= cos(theta);//ftarget z

		 //look direction in beam frame
		 for(unsigned int l=0;l<3;++l){
		   look_old_beam[l]=0.0;
		   for(unsigned int m=0;m<3;++m)
		     look_old_beam[l] += TR_from_newbeam_to_oldbeam[l][m]*look_new_beam[m];
		 }

		 //azi elev
		 beam_azim_elev(look_old_beam[0], look_old_beam[1], look_old_beam[2], azi,elev);

		 //debug
		 

		 //gain
		 gain=Bvector[i_beam].bilinear(azi,elev);
		 beam_integral=beam_integral  + gain*(target_temperature(j_target)-Tcoldsky)*solid_angle/integrated_beam_gain[i_beam];
	       }


	     double true_value_of_solid_angle= 2.0*pi*(1.0-cos(max_theta*d2r));
	     if( fabs(solid_angle_integral-true_value_of_solid_angle)/true_value_of_solid_angle*100.0>1.0)
	       ErrorMessage("mj_rad_proc::Solid angle integration error while beam integration over Target").throwMe();
	    
	     cout<<"beam integral over target( Ttarget- Tcoldsky) : "<< beam_integral<<endl;
	     sum= sum + beam_integral;
	   }//for each target
	   median_zero_Ta[i_cal]= median_zero_Ta[i_cal]  - sum;
	   beam_integration_power[i_cal]=sum;
	   cout<<"cal index, sum "<< i_cal<<" "<< "of "<< median_zero_Ta.size()<<" "<<sum<<endl;
	 }//for each cal point


	 //-------------------------
	 //show this integrated beam power
	 //-----------------------
	 Plot beam_power;
	 beam_power.addXY(zero_Ta_cal_time,"min",beam_integration_power,"K",line("none"),sym("circle","red",1));
	 beam_power.setTitle("integrated beam power vs time");
	 beam_power.setFile(prefix_str1 +"beam_integral.eps");
	 beam_power.setCmdFile(prefix_str1 +"beam_integra.xmgr");
	 beam_power.show(plot_target);




	 //------------------------------------------
	 //Use polynomial fitting of median zero Ta
	 //-----------------------------------------
	 cout<<"Poly fit "<<endl;
	
	 Dvec input_time("",median_zero_Ta.size());
	 Dvec input_ta("",median_zero_Ta.size());
	 
	 Dvec poly_fit_coeff("",Npolyfit+1);
	 Uvar med_time=zero_Ta_cal_time[(unsigned int) zero_Ta_cal_time.size()/2];
	 for(unsigned int i=0;i<zero_Ta_cal_time.size();++i){
	   input_time(i)= (zero_Ta_cal_time[i]-med_time).getInUnits("s")/60.0;//use min as base
	   input_ta(i)= median_zero_Ta[i].getInUnits("K");
	 }
	 polyfit(poly_fit_coeff,input_ta,input_time);

	
	 //-----------------------------------------------------
	 //Display uncal Ta and time when farthest from Target
	 //------------------------------------------------------
	 if(zero_Ta.size()!=zero_Ta_cal_time.size())
	   ErrorMessage("zero Ta and time containers size mismatch").throwMe();
	 if(zero_Ta.size()==0){
	   ErrorMessage("there is no cold sky pointing: may need to adjust cold sky condition in the config file").throwMe();
	 }
	 else{
	   vector<Uvar> fit_ta;
	   fit_ta.resize(zero_Ta_cal_time.size());
	   for(unsigned int i=0;i<zero_Ta_cal_time.size();++i){
	     double y=poly_fit_coeff(0);
	     double x= (zero_Ta_cal_time[i]-med_time).getInUnits("s")/60.0;//use min as base
	     for(unsigned int j=0;j<Npolyfit;++j){
	       y += poly_fit_coeff(j+1)*pow(x,j+1);
	     }
	     fit_ta[i]=Uvar(y,"K");
	   }
	   
	   Plot Ta_cal;
	   Ta_cal.addXY(zero_Ta_cal_time,"min", median_zero_Ta,"K",line("none"),sym("circle","red",1));
	   Ta_cal.addXY(zero_Ta_cal_time,"min", fit_ta,"K",line("solid","red",1),sym("none"));
	   Ta_cal.setTitle("zero Ta when pointing is farthest from Target");
	   Ta_cal.setXlabelandUnits("time in min");
	   Ta_cal.setYlabelandUnits("cold sky reference Ta");
	   Ta_cal.setFile(prefix_str1 +"zero_temp.eps");
	   Ta_cal.setCmdFile(prefix_str1 +"zero_temp.xmgr");
	   Ta_cal.show(plot_target); 
	 }

	 


	 //-------------------
	 //Rescale uncal_Ta
	 //--------------------
	 cout<<"Rescale uncal Ta "<<endl;
	 for (unsigned int i=0;i<Nbeam_record;++i){
	   Uvar zero;
	   //-----
	   //time
	   //----------
	   if(beam_rel_time[i]<=zero_Ta_cal_time[0])
	     zero= median_zero_Ta[0];
	   else if(beam_rel_time[i] >= zero_Ta_cal_time[zero_Ta_cal_time.size()-1])
	     zero=median_zero_Ta[median_zero_Ta.size()-1];
	   else{
	     //-------------------------
	     //Use polynomial fit
	     //-------------------------
	     double x= (beam_rel_time[i]-med_time).getInUnits("s")/60.0;//use min as base
	     double y=poly_fit_coeff(0);
	    
	     for(unsigned int j=0;j<Npolyfit;++j){
	       y += poly_fit_coeff(j+1)*pow(x,j+1);
	     }
	     zero= Uvar(y,"K");
	   }//end of else
	   uncal_Ta[i]= uncal_Ta[i]-zero;
	 }
	 

	 //-----------------------
	 //antenna temperature, gain
	 //-----------------------
	 cout<<"compute antenna temperature "<<endl;
	 vector<Uvar> beam_ant_temperature, beam_sys_gain, receiver_temp;
	
	 for (unsigned int i=0;i<Nbeam_record;++i){
	   beam_ant_temperature.push_back(uncal_Ta[i]+Tcoldsky);//add Tcoldsky
	   beam_sys_gain.push_back(radiometer_cal*Troom/((beam_rlotmp[i]-Tcoldsky)*beam_relative_gain[i]));
	   receiver_temp.push_back(beam_sys_gain.back()*beam_ncnt_radio[i]- beam_ant_temperature.back());
	 }
	 Uvar a0,a1,a2,a_avg,a_sqr_diff;
	 vector<Uvar> beam_ant_temperature_std;
	 beam_ant_temperature_std.resize(Nbeam_record);
	 for(unsigned int i=1;i<Nbeam_record-1;++i){
	   a0= beam_ant_temperature[i-1];
	   a1= beam_ant_temperature[i];
	   a2= beam_ant_temperature[i+1];
	   a_avg= (a0+a1+a2)/3.0;
	   a_sqr_diff= pow(a0-a_avg,2) + pow(a1-a_avg,2)+pow(a2-a_avg,2);
	   beam_ant_temperature_std[i]= sqrt(a_sqr_diff)/sqrt(2.0);
	 }
	 beam_ant_temperature_std[0]=beam_ant_temperature_std[1];
	 beam_ant_temperature_std[Nbeam_record-1]= beam_ant_temperature_std[Nbeam_record-2];
	



	 //------------------
	 //Save output
	 //-----------------
	 for(unsigned int i=0;i<Nbeam_record;++i){
	   sab_vs_antenna_temperature[beam_sab_counter[i]]= beam_ant_temperature[i];
	   sab_vs_system_gain[beam_sab_counter[i]]= beam_sys_gain[i];
	   sab_vs_receiver_temperature[beam_sab_counter[i]]= receiver_temp[i];
	   sab_vs_antenna_temperature_stddev[beam_sab_counter[i]]=beam_ant_temperature_std[i];

	   output_file<<i_counter<<" "<<beam_sys_gain[i].getInUnits("K s")<<" "<<beam_ant_temperature[i].getInUnits("K")<<" "<<receiver_temp[i].getInUnits("K")<<" "<< beam_ant_temperature_std[i].getInUnits("K")<<endl;
	   i_counter++;
	 }

	 

	 //--------------
	 //Must see plot
	 //--------------
	 cout<<"Make plots "<<endl;
	 Plot p_ant;
	 p_ant.addXY(beam_rel_time,"min", beam_ant_temperature,"K",line("solid","red",1),sym("none"));
	 p_ant.setTitle("beam antenna temperature");
	 p_ant.setXlabelandUnits("time in min");
	 p_ant.setYlabelandUnits("Ta");
	 p_ant.setFile(prefix_str1 +"ant_temp.eps");
	 p_ant.setCmdFile(prefix_str1 +"ant_temp.xmgr");
	 p_ant.show(plot_target); 

	

	 Plot p_ant_std;
	 p_ant_std.addXY(beam_rel_time,"min", beam_ant_temperature_std,"K",line("solid","red",1),sym("none"));
	 p_ant_std.setTitle("beam antenna temperature  STDDEV");
	 p_ant_std.setXlabelandUnits("time in min");
	 p_ant_std.setYlabelandUnits("Ta STDDEV");
	 p_ant_std.setFile(prefix_str1 +"ant_std_temp.eps");
	 p_ant_std.setCmdFile(prefix_str1 +"ant_std_temp.xmgr");
	 p_ant_std.show(plot_target); 
	



	 Plot p_gain;
	 p_gain.addXY(beam_rel_time,"min", beam_sys_gain,"K s",line("solid","red",1),sym("none"));
	 p_gain.setTitle("beam system gain");
	 p_gain.setXlabelandUnits("time in min");
	 p_gain.setYlabelandUnits("Gain K-s");
	 p_gain.setFile(prefix_str1 +"gain.eps");
	 p_gain.setCmdFile(prefix_str1 +"gain.xmgr");
	 p_gain.show(plot_target); 

	 

	 Plot p_rec;
	 p_rec.addXY(beam_rel_time,"min", receiver_temp,"K",line("solid","red",1),sym("none"));
	 p_rec.setTitle("receiver temperature");
	 p_rec.setXlabelandUnits("time in min");
	 p_rec.setYlabelandUnits("receiver temperature(K)");
	 p_rec.setFile(prefix_str1 +"rec_temp.eps");
	 p_rec.setCmdFile(prefix_str1 +"rec_temp.xmgr");
	 p_rec.show(plot_target); 



   
       }//for each beam


     



     //-----------------
     //Write outputs back to LBDR(if exists) and SBDR
     // L1BP Data
     //----------------
     //-----------------------------------
     //Start with SBDR:  set file pointer right after header
     //-----------------------------------     
     data.gotoFirstRecord();//SBDR
     
     //default value setting
     Uvar no_avail_cal = Uvar(0,"K");
     Uvar sab_antenna_brightness_temp=Uvar(0,"K");
     Uvar sab_antenna_brightness_temp_stddev=Uvar(0,"K");
     Uvar sab_system_gain=Uvar(0,"K s");
     Uvar sab_receiver_temp=Uvar(0,"K");
     
     //iterator to each map variable
     map<Uvar,Uvar>::const_iterator brightness_temp_pointer;
     map<Uvar,Uvar>::const_iterator system_gain_pointer;
     map<Uvar,Uvar>::const_iterator system_noise_temp_pointer;
     map<Uvar,Uvar>::const_iterator brightness_temp_stddev_pointer;
     Uvar  sab_counter_check;
     bool sys_gain_set, ant_temp_set, system_noise_temp_set,abt_std_set;
     unsigned long rad_science_flag;
     cout<<"write radiometer data into SBDR "<<endl;
     for(unsigned int i_rec=0;i_rec<Nrecord;++i_rec)
       {
	 //-------------------
	 //default setting: no value is set
	 //--------------------
	 sys_gain_set = false;
	 ant_temp_set = false;
	 system_noise_temp_set = false;
	 abt_std_set = false;
	 
	 //---------------------------
	 //start with default value
	 //----------------------------
	 sab_antenna_brightness_temp=Uvar(0,"K");
	 sab_antenna_brightness_temp_stddev=Uvar(0,"K");
	 sab_system_gain=Uvar(0,"K s");
	 sab_receiver_temp=Uvar(0,"K");
	 
	 //------------------------------------------------------------------------
	 //check sab counter match: sab_counter(i_rec) should match to sab_counter
	 //-----------------------------------------------------------------------
	 sab_counter_check=data.getSabCounter();
	 if(sab_counter(i_rec) !=sab_counter_check) {
	   ErrorMessage("rad_proc.cpp::sab counter mismatch for input and output").throwMe();  } 
	 
	 //---------------------------------------------------
	 //check whether brightness temperature is calibrated
	 //if not, write default value
	 //sab_vs_antenna_temperature
	 //sab_vs_system_gain[beam_sab_counter
	 //sab_vs_receiver_temperature[beam_sab_counter
	 //sab_vs_antenna_temperature_stddev
	 //----------------------------------------------------
	 brightness_temp_pointer = sab_vs_antenna_temperature.find(sab_counter(i_rec));
	 if(brightness_temp_pointer!=sab_vs_antenna_temperature.end())
	   {
	     ant_temp_set = true;
	     sab_antenna_brightness_temp = (brightness_temp_pointer)->second;
	   }

	 //----------------------
	 //check std dev of antenna temperature
	 //----------------------
	 brightness_temp_stddev_pointer=sab_vs_antenna_temperature_stddev.find(sab_counter(i_rec));
	 if(brightness_temp_stddev_pointer!=sab_vs_antenna_temperature_stddev.end()){
	   abt_std_set=true;
	   sab_antenna_brightness_temp_stddev=(brightness_temp_stddev_pointer)->second;
	 }


	 //----------------------------
	 //systerm gain pointer
	 //-----------------------------
	 system_gain_pointer = sab_vs_system_gain.find(sab_counter(i_rec));
	 if(system_gain_pointer != sab_vs_system_gain.end())
	   {
	     sys_gain_set = true;
	     sab_system_gain = (system_gain_pointer)->second;
	   }

	 //---------------------
	 //system receiver temperature
	 //--------------------------
	 system_noise_temp_pointer = sab_vs_receiver_temperature.find(sab_counter(i_rec));
	 if(system_noise_temp_pointer!=sab_vs_receiver_temperature.end())
	   {
	     system_noise_temp_set = true;
	     sab_receiver_temp =(system_noise_temp_pointer)->second;
	   }
	 
	 if( sys_gain_set = true &&
	     ant_temp_set ==true &&
	     system_noise_temp_set == true&&
	     abt_std_set == true){
	   data.writeParameter("antenna_brightness_temp", sab_antenna_brightness_temp,"K");
	   data.writeParameter("abt_std",sab_antenna_brightness_temp_stddev,"K");
	   data.writeParameter("system_gain",sab_system_gain,"K s");
	   data.writeParameter("system_noise_temp",sab_receiver_temp,"K");
	 }
	 
	 
	 //set science flag
	 data.readParameter("science_qual_flag");//read first
	 rad_science_flag = data.science_qual_flag;//retrieve value
	 if( sys_gain_set = true &&
	     ant_temp_set ==true &&
	     system_noise_temp_set == true&&
	     abt_std_set == true)
	   {
	     bitset(rad_science_flag,4,0);
	   }
	 else
	   {
	     bitset(rad_science_flag,4,1);
	   }
	 data.writeParameter("science_qual_flag",rad_science_flag);
	 data.skipRecord();//go to next record
       }
     

     //------------------------------------------
     // LBDR file update if there is any LBDR file
     //check there is an LBDR file
     //---------------------------------------------
     cout<<"write radiometer data into LBDR "<<endl;
     ifstream active_file(l1b_a_filename.c_str());
     if(!active_file.is_open()){
       cout<<"No active record is avail "<<endl;
       exit(0);//stop here
     }
     active_file.close();
     
     L1B data_a(l1b_a_filename,"r+","active");
     data_a.readHeader();
     unsigned int Nactive_record = data_a.recordCount();
     data_a.gotoFirstRecord();
     
     for(unsigned int i=0; i<Nactive_record;++i)
       {
	 //-------------------
	 //default setting: no value is set
	 //--------------------
	 sys_gain_set = false;
	 ant_temp_set = false;
	 system_noise_temp_set = false;
	 abt_std_set = false;
	 
	 //---------------------------
	 //start with default value
	 //----------------------------
	 sab_antenna_brightness_temp=Uvar(0,"K");
	 sab_antenna_brightness_temp_stddev=Uvar(0,"K");
	 sab_system_gain=Uvar(0,"K s");
	 sab_receiver_temp=Uvar(0,"K");
	 
	 //------------------------------------------------------------------------
	 //check sab counter match: sab_counter(i_rec) should match to sab_counter
	 //-----------------------------------------------------------------------
	 sab_counter_check=data_a.getSabCounter();
	
	 //---------------------------------------------------
	 //check whether brightness temperature is calibrated
	 //if not, write default value
	 //sab_vs_antenna_temperature
	 //sab_vs_system_gain[beam_sab_counter
	 //sab_vs_receiver_temperature[beam_sab_counter
	 //sab_vs_antenna_temperature_stddev
	 //----------------------------------------------------
	 brightness_temp_pointer = sab_vs_antenna_temperature.find(sab_counter_check);
	 if(brightness_temp_pointer!=sab_vs_antenna_temperature.end())
	   {
	     ant_temp_set = true;
	     sab_antenna_brightness_temp = (brightness_temp_pointer)->second;
	   }

	 //----------------------
	 //check std dev of antenna temperature
	 //----------------------
	 brightness_temp_stddev_pointer=sab_vs_antenna_temperature_stddev.find(sab_counter_check);
	 if(brightness_temp_stddev_pointer!=sab_vs_antenna_temperature_stddev.end()){
	   abt_std_set=true;
	   sab_antenna_brightness_temp_stddev=(brightness_temp_stddev_pointer)->second;
	 }
	 
	 //----------------------------
	 //systerm gain pointer
	 //-----------------------------
	 system_gain_pointer = sab_vs_system_gain.find(sab_counter_check);
	 if(system_gain_pointer != sab_vs_system_gain.end())
	   {
	     sys_gain_set = true;
	     sab_system_gain = (system_gain_pointer)->second;
	   }

	 //---------------------
	 //system receiver temperature
	 //--------------------------
	 system_noise_temp_pointer = sab_vs_receiver_temperature.find(sab_counter_check);
	 if(system_noise_temp_pointer!=sab_vs_receiver_temperature.end())
	   {
	     system_noise_temp_set = true;
	     sab_receiver_temp =(system_noise_temp_pointer)->second;
	   }
	 
	 if( sys_gain_set = true &&
	     ant_temp_set ==true &&
	     system_noise_temp_set == true&&
	     abt_std_set == true){
	   data_a.writeParameter("antenna_brightness_temp", sab_antenna_brightness_temp,"K");
	   data_a.writeParameter("abt_std",sab_antenna_brightness_temp_stddev,"K");
	   data_a.writeParameter("system_gain",sab_system_gain,"K s");
	   data_a.writeParameter("system_noise_temp",sab_receiver_temp,"K");
	 }
	 
	 
	 //set science flag
	 data_a.readParameter("science_qual_flag");//first read before setting value
	 rad_science_flag = data_a.science_qual_flag;//read first
	 if( sys_gain_set = true &&
	     ant_temp_set ==true &&
	     system_noise_temp_set == true&&
	     abt_std_set == true)
	   {
	     bitset(rad_science_flag,4,0);
	   }
	 else
	   {
	     bitset(rad_science_flag,4,1);
	   }
	 data_a.writeParameter("science_qual_flag",rad_science_flag);
	 data_a.skipRecord();//go to next record
       }       

    

  }//end of try

	 
 
 
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


