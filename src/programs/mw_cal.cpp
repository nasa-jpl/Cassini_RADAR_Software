//------------------------------
// radiometer processor
//-------------------------------
#include <string.h>
#include <iostream>
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
    string target= cfg.str("target");//assume single target 
    if(strcasecmp(target.c_str(),"source")==0) target="none";
    Frame fj2000e("J2000","Earth");
   
  
    //---------------------------------------------------------
    //read parameters from config file
    //         : system temperature
    //         : front_end_loss, min_beam_distance
    //         : cold sky temperature
    //         : system gain
    //
    //         system gain: 0.106 counts/ms/K --> 106 counts/s/K 
    //        : cold sky average time
    //----------------------------------------------------------
    /*
    Uvar Tcfg_sys=cfg["Tsys"];
    Uvar front_end_loss_dB=cfg["frontend_loss_dB"];
    double front_end_loss_figure = pow(10, front_end_loss_dB.getInUnits("")/10.0);
    Uvar min_beam_distance_in_units_of_beam3dB=cfg["min_beam_distance"];
    double min_beam_distance =  min_beam_distance_in_units_of_beam3dB.getInUnits("");
    Uvar Tcoldsky=cfg["cold_sky_Temperature"]; 
    Uvar Gcfg_sys=cfg["radiometer_system_gain"];
    Uvar avg_time=cfg["cold_sky_average_time"];

    cout<<"system noise temperature "<<Tcfg_sys<<endl;
    cout<<"front_end_losss_in_dB "<<front_end_loss_dB<<endl;
    cout<<"min beam distance between boresight and limb direction "<<min_beam_distance_in_units_of_beam3dB<<endl;
    cout<<"cold sky temperature "<<Tcoldsky<<endl;
    cout<<"const system gain "<<Gcfg_sys<<endl;
    cout<<"cold sky average time "<< avg_time<<endl;
    */
    //------------------------------------------------------
    //Read in Microwave source names
    // and their RADEC coordinate values from config file
    //--------------------------------------------------------
    unsigned int Nsource=(unsigned int) cfg.getInt("NUM_SOURCES");
    Ivec source_ids("source id",Nsource);
    Uvec source_decs("MW_source_dec",Nsource);
    Uvec source_ras("MW_source_ra",Nsource);
    Array1D<string.h> source_names("source names",Nsource);
    for(unsigned int i=0;i<Nsource;i++)
      {
	source_ids(i)=cfg.getInt("SOURCE_ID_"+toStr(i+1));
	source_decs(i)=cfg["SOURCE_DEC_"+toStr(i+1)];
	source_ras(i)=cfg["SOURCE_RA_"+toStr(i+1)];
	source_names(i)=cfg.str("SOURCE_NAME_"+toStr(i+1));
	//cout<<source_names(i)<<endl;
	//cout<<source_ras(i)<<" "<<source_decs(i)<<endl;
      }
    if(source_names.size() !=Nsource) ErrorMessage("Source name list size and source size does not match").throwMe();


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
    BurstData* data;
    string l1b_filename = cfg.str("L1B_P_filename");
    data = new L1B(l1b_filename,"rb","passive");
    unsigned int Nrecord = data->recordCount();
  
 
       
    //----------------------------------------------------------------
    //(0) Setting up resulting variables which are being read from L1BP
    //-----------------------------------------------------------------
    Uvec sab_counter("sab counter",Nrecord);
    Array1D<unsigned int> quality_flag("quality flag",Nrecord);
    Uvec ncnt_rl("ncnt_rl",Nrecord);
    Uvec ncnt_nd("ncnt_nd",Nrecord);
    Uvec ncnt_radio("ncnt_radio",Nrecord);
    Array1D<unsigned int>  beam_number("beam number",Nrecord);
    Array1D<Time> time("time",Nrecord);
    Array1D<Uvar> rlotmp("rlotmp",Nrecord);
    Array1D<Uvar> lnatmp("lna tmp",Nrecord);
    Array1D<Uvar> adctmp("adc tmp",Nrecord);
    vector<Uvec>  wgbtmp; wgbtmp.clear();//lower waveguide temperature
   

    //---------------------------------------------------
    //(1) Extracted parameters from L1B
    // sab_counter, norm_cnt_rl, norm_cnt_nd, norm_cnt_radio
    // qaulity_flag, sclk, brst, beam_number, rlotmp, lnatmp
    // adctmp, wgb1t1, wgb3t1, wgb3t2, wgb3t3, wgb5t1
    //----------------------------------------------------
    unsigned int Nparam=16;
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

    unsigned int i_count = 0;
    double sclk_double, brst_double;
    Uvec wgb_dummy("dummy",5);
    
     while(! data->eof()){
       if (!(*data).returnParams(return_doubles,param_list,
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
	  i_count++;
	}
     }
     cout<<"Finished extracting radiometer relevant paramenters from L1BP"<<endl;     
     
   
     //--------------------
     //(2) cold sky determination and microwave source identification
     //--------------------
     vector<bool> is_beam_coldsky;
     vector<string.h> sab_source_name;
     is_beam_coldsky.clear();
     sab_source_name.clear();
    
     Uvec  mw_source_distances("", Nsource);
     unsigned int beam_id =3;//beam 3 only calibration
   
     //collect beam_id specific data 
     vector<Uvar> beam_sab_counter;
     vector<Uvar> beam_rel_time;
     vector<Uvar> beam_ncnt_rl;
     vector<Uvar> beam_ncnt_nd;
     vector<Uvar> beam_ncnt_radio;
     vector<Uvar> beam_rlotmp;
     vector<Uvar> beam_lnatmp;
     vector<Uvar> beam_adctmp;
     vector<Uvar> beam_wgbtmp;
     vector<Uvar> signal_count;
     vector<Uvar> signal_sab;
     vector<Uvar> signal_time;
     vector<Uvar> beam_distance;
     vector<Uvar> signal_id;
     Uvec source_azim("",Nsource);
     Uvec source_elev("",Nsource);

     vector<Uvar>  azim_index;
     vector<Uvar> elev_index;

     //clear the container
     beam_sab_counter.clear();
     beam_rel_time.clear();
     beam_ncnt_rl.clear();
     beam_ncnt_nd.clear();
     beam_ncnt_radio.clear();
     beam_rlotmp.clear();
     beam_lnatmp.clear();
     beam_adctmp.clear();
     beam_wgbtmp.clear();
     signal_count.clear();
     signal_sab.clear();
     beam_distance.clear();
     azim_index.clear();
     elev_index.clear();
     signal_id.clear();
     

     Uvec dummy("dummy",5);  
     string mw_source_name;
     bool coldsky_indicator;
     int mw_source_id;
     Uvar mw_azi, mw_elev;
       
     for(unsigned int i_rec=0;i_rec<Nrecord;i_rec++)
       {
	 if( bitget(quality_flag(i_rec),0,0)==1 ||
	     bitget(quality_flag(i_rec),1,1)==2 || 
	     beam_number(i_rec)!=beam_id ) continue;

	 beam_sab_counter.push_back(sab_counter(i_rec));
	 beam_rel_time.push_back(time(i_rec)-time(0));
	 beam_ncnt_rl.push_back(ncnt_rl(i_rec));
	 beam_ncnt_nd.push_back(ncnt_nd(i_rec));
	 beam_ncnt_radio.push_back(ncnt_radio(i_rec));
	 beam_rlotmp.push_back(rlotmp(i_rec));
	 beam_lnatmp.push_back(lnatmp(i_rec));
	 beam_adctmp.push_back(adctmp(i_rec));
	
	 if(beam_id ==1 || beam_id==2) beam_wgbtmp.push_back(dummy(0));
	 if(beam_id==3) beam_wgbtmp.push_back(dummy(2));
	 if(beam_id==4 || beam_id==5) beam_wgbtmp.push_back(dummy(4));


	
	 //-----------------------------------------
	 // Compute Microwave sources' position in
	 // beam3 frame (azi, elev)
	 //-------------------------------------------
	 for(unsigned int i_mw=0;i_mw<Nsource;i_mw++)
	   { 
	     DirectionVector mw_source_dir("mc_source",fj2000e,
					   time(i_rec),0,0,1);//fj2000e: J2000 centered at earth
	     mw_source_dir.setRADEC(source_ras(i_mw),source_decs(i_mw));
	     mw_source_dir.representIn(Fvector[beam_number(i_rec)-1]);
	     mw_source_dir.getAzimuthElevation(mw_azi,mw_elev);
	     source_azim(i_mw) = mw_azi;
	     source_elev(i_mw)= mw_elev;
	     mw_source_distances(i_mw)= sqrt(mw_azi*mw_azi + mw_elev*mw_elev);
	     azim_index.push_back(mw_azi);
	     elev_index.push_back(mw_elev);
	   }//loop over each source
	 
	 //take the source name
	 //determine beam's pointing
	 coldsky_indicator= true;

	 for(unsigned int i_mw=0;i_mw<Nsource;++i_mw)
	   {
	     if(source_azim(i_mw) > -azi_beam3dB_oneway(beam_id-1)*2
		&& source_azim(i_mw) < azi_beam3dB_oneway(beam_id-1)*2
		&& source_elev(i_mw) > -elev_beam3dB_oneway(beam_id-1)*2
		&& source_elev(i_mw) < elev_beam3dB_oneway(beam_id-1)*2)
	       {
		 coldsky_indicator= false;
		 mw_source_id=i_mw;
		 mw_source_name=source_names(i_mw);
		 //cout<<"target "<< mw_source_name<<endl;
	       }
	   }
	 is_beam_coldsky.push_back(coldsky_indicator);
	 sab_source_name.push_back(mw_source_name);
	 if(!coldsky_indicator)
	   {
	     signal_sab.push_back(sab_counter(i_rec));
	     signal_count.push_back(ncnt_radio(i_rec));
	     signal_time.push_back((time(i_rec)-time(0)));
	     signal_id.push_back(mw_source_id);
	   }
       }

     //-------------------------------
     //(3) Plot each microwave source's  azi and elev
     //-------------------------------
     cout<<"sab beam counter "<<beam_sab_counter.size()<<endl;
     cout<<"index size "<<azim_index.size()<<" "<<elev_index.size()<<endl;
     for(unsigned int i_mw = 0; i_mw<Nsource;i_mw++)
       {
	 Plot a;
	 vector<Uvar> x; vector<Uvar> y;
	 x.clear(); y.clear();
	 if(azim_index.size()!= elev_index.size()) ErrorMessage("index mismatch").throwMe();
	 for(unsigned int j=0;j<azim_index.size();++j)
	   {
	     if(j%Nsource == i_mw)
	       {
		 x.push_back(azim_index[j]);
		 y.push_back(elev_index[j]);
	       }
	   }
	 cout<<"size "<<x.size()<<" "<<y.size()<<endl;
	 a.addXY(x,"deg",y,"deg",line("solid","red",2),sym("none"));
	 a.setFile(prefix_str +toStr(i_mw+1)+ "mw_pointing.eps");
	 a.setCmdFile(prefix_str +toStr(i_mw+1)+ "mw_pointing.xmgr");
	 a.setTitle(source_names(i_mw));
	 a.show(plot_target);
       }
    
     unsigned int Nbeam_record=beam_sab_counter.size();
     cout<<Nbeam_record<<"  samples were chosen for beam number "
	 <<beam_id<<endl;
     cout<<signal_sab.size()<< "targeted sabs"<<endl;
     //make sure all the container have the same size
     if(Nbeam_record!=beam_rel_time.size() ||
	Nbeam_record!=beam_ncnt_rl.size() ||
	Nbeam_record!=beam_ncnt_nd.size() ||
	Nbeam_record!=beam_ncnt_radio.size() ||
	Nbeam_record!=beam_rlotmp.size()||
	Nbeam_record != beam_lnatmp.size()||
	Nbeam_record != beam_adctmp.size()||
	Nbeam_record != beam_wgbtmp.size() ||
	Nbeam_record != is_beam_coldsky.size()||
	Nbeam_record != sab_source_name.size())
       {
	 cout<<Nbeam_record<<" "<<beam_rel_time.size()
	     <<" "<<beam_ncnt_rl.size() 
	     <<" "<<beam_ncnt_nd.size()
	     <<" "<<beam_ncnt_radio.size() 
	     <<" "<<beam_rlotmp.size()
	     <<" "<< beam_lnatmp.size()
	     <<" "<< beam_adctmp.size()
	     <<" "<< beam_wgbtmp.size() 
	     <<" "<< is_beam_coldsky.size()
	     <<" "<< sab_source_name.size()<<endl;
	 ErrorMessage("Container size mismatch").throwMe();
       }

   
     //debug
     //--------------------------
     //(4) normalized count with target id
     //---------------------------
     //plot sab normalized radio count
     if(target!="")
       {
	 Plot a;
	 a.addXY(beam_sab_counter," ",beam_ncnt_radio,"1/ms",line("solid","yellow",1),sym("none"));
	 a.addLegend("norm count");
	
	 for(unsigned int i_mw =0;i_mw<Nsource;++i_mw)
	   {
	     vector<Uvar> x,y;
	     x.clear();y.clear();
	     for(unsigned int i=0;i<signal_sab.size();++i)
	       {
		 if(signal_id[i]==i_mw)
		   {
		     x.push_back(signal_sab[i]);
		     y.push_back(signal_count[i]);
		   }
	       }
	     if(x.size()!=0 && y.size()!=0)
	       {
		 a.addLegend(source_names(i_mw)); 
		 a.addXY(x," ",y,"1/ms",line("none"),sym("circle",color_table(i_mw),.25));
	       }
	   }
	 a.setTitle("normalized radio count (counts/ms)");
	 a.setFile(prefix_str + "norm_radio_targets.eps");
	 a.setCmdFile(prefix_str + "norm_radio_targets.xmgr");
	 a.show(plot_target);
       }

     //--------------------------------------
     /*  stop here for  microwave calibration
    //--------------------------------------------
     //scan segmentation
     vector<unsigned int> cal_start,cal_end;
     vector<Uvar> Gsys_cal_start, Gsys_cal_end;
     vector<Uvar> Time_cal_start,Time_cal_end;
     vector<Uvar> Sab_cal_start, Sab_cal_end;
     cal_start.clear();
     cal_end.clear();
     Gsys_cal_start.clear();
     Gsys_cal_end.clear();
     Time_cal_start.clear();
     Time_cal_end.clear();
     Sab_cal_start.clear();
     Sab_cal_end.clear();
     for(unsigned int i=0;i<Nbeam_record-1;i++)
       {
	 
	 if(!is_beam_coldsky[i] && is_beam_coldsky[i+1])
	   {
	     cal_start.push_back(i+1);  
	     Gsys_cal_start.push_back(CalibrateGsysUsingConstTsys
				      (i+1, Tcfg_sys, 
				       avg_time,  Tcoldsky,
				       beam_rel_time,
				       beam_ncnt_radio,
				       is_beam_coldsky));
	     Time_cal_start.push_back(beam_rel_time[i+1]);
	     Sab_cal_start.push_back(beam_sab_counter[i+1]);
	     
	   }
	 if(is_beam_coldsky[i] && !is_beam_coldsky[i+1]) 
	   {
	     cal_end.push_back(i);
	     Gsys_cal_end.push_back(CalibrateGsysUsingConstTsys
				    (i, Tcfg_sys, avg_time,  Tcoldsky,
				     beam_rel_time,beam_ncnt_radio,is_beam_coldsky));
	     Time_cal_end.push_back(beam_rel_time[i]);  
	     Sab_cal_end.push_back(beam_sab_counter[i]);
	   }
       } 
	   
     unsigned int Ncal_start = cal_start.size();
     if(Ncal_start != Gsys_cal_start.size() ||
	Ncal_start != Time_cal_start.size()||
	Ncal_start != Sab_cal_start.size()) { 	 
       cout<<"sizes "<<Gsys_cal_start.size()<<" "<<Time_cal_start.size()
	   <<" "<<Sab_cal_start.size()<<" "<<endl;
       ErrorMessage("Ncal start container size mismatch").throwMe();
     }
     
     unsigned int Ncal_end = cal_end.size();
     if(Ncal_end != Gsys_cal_end.size() ||
	Ncal_end != Time_cal_end.size() ||
	Ncal_end != Sab_cal_end.size()) {
       cout<<"sizes "<<Gsys_cal_end.size()<<" "<<Time_cal_end.size()
	   <<" "<<Sab_cal_end.size()<<" "<<endl;
       ErrorMessage("Ncal end container size mismatch");}
     
     
     cout<<"cal start size "<<Ncal_start<<" cal end size "<<Ncal_end<<endl;
     //-----------------------------
     //Plot scan segmented system gain
     //---------------------------
     if(Ncal_end !=0 && Ncal_start !=0)
       {
	 Plot a;
	 a.addXY(Sab_cal_start,"",Gsys_cal_start,"1/(s K)",line("solid","black",2),sym("none"));
	 a.addXY(Sab_cal_end,"",Gsys_cal_end,"1/(s K)",line("solid","red",2),sym("none"));
	 a.setTitle("black:cal start, red: cal end vs sab counter:beam"+toStr(beam_id));
	 a.setXlabelandUnits("sab number");
	 a.setYlabelandUnits("scan edge cal. temperature");
	 a.setFile(prefix_str + "scan_edge.eps");
	 a.setCmdFile(prefix_str + "scan_edge.xmgr");
	 a.show(plot_target);
       }

     //-----------------------------
     //(10): Now compute antennam brightness temperature
     //calibrate at scan edges
     // cold_sky_counts/s = G (Tsys + Tsky)
     // RL_counts/s = G (Tsys + Trl)
     // Y = cold_sky_counts/ RL_counts
     //  Ts = (Tsky - Y * Trl)/(Y-1)
     //  G (counts/s/T) = RL_counts/s /(Tsys+Trl)
     //
     //------------------------------
     
     Uvec antenna_temperature("antenna",Nbeam_record);
     Uvec system_gain("system gain",Nbeam_record);
     Uvec sab_counter_index("sab index",Nbeam_record);
     antenna_temperature=Uvar(0,"K");
     system_gain=Uvar(0,"1/(s K)");
	   
     Uvar Gsys,Tbright, Ta,alpha;
     for (unsigned int i_beam_rec=0; i_beam_rec<Nbeam_record;i_beam_rec++)
       {
	 //reset recycled variables
	 sab_counter_index(i_beam_rec)=beam_sab_counter[i_beam_rec];
	 Gsys = Uvar(0,"1/(s K)");//reset
	 Tbright=Uvar(0,"K");//reset
	 Ta=Uvar(0,"K");
	 alpha = 0;
	 
	 if(is_beam_coldsky[i_beam_rec])
	   {
	     //----------------------
	     //cold sky pointing
	     //---------------------
	     Gsys=CalibrateGsysUsingConstTsys
	       (i_beam_rec,Tcfg_sys, avg_time,  Tcoldsky,
		beam_rel_time,  beam_ncnt_radio,  is_beam_coldsky);
	     
	     Ta = beam_ncnt_radio[i_beam_rec]/Gsys - Tcfg_sys;
	     system_gain(i_beam_rec) = Gsys;
	     antenna_temperature(i_beam_rec)=Ta;
	   }
	             
	 else if(cal_end.size()!= 0 && cal_end.size()!=0)
	   {
	     //-----------------
	     //(solar) target pointing
	     //CASE-A: there are cold sky calibration points avail
	     //------------------
	     unsigned int last_cal_end ;
	     bool last_cal_end_set = false;
	     for(unsigned int i_cal=0;i_cal<Ncal_end-1;i_cal++)
	       {
		 if(i_beam_rec >cal_end[i_cal] && i_beam_rec <cal_end[i_cal+1])
		   {
		     last_cal_end = i_cal;
		     last_cal_end_set = true;
		     break;
		   } 
	       }
	     if(last_cal_end_set==false && i_beam_rec>cal_end[Ncal_end-1])
	       {//take care of last "cal_end"
		 //this happens last pointing stays inside target 
		 last_cal_end = Ncal_end -1;
		 last_cal_end_set=true;
	       }
	     
	     unsigned int next_cal_start;
	     bool next_cal_start_set = false;
	     if(i_beam_rec < cal_start[0] )
	       {//take care of first "cal_start"
		 //when scan start even before calibration is done
		 next_cal_start = 0;
		 next_cal_start_set= true; 
	       }
	     else
	       {
		 for(unsigned int i_cal=1;i_cal<Ncal_start;i_cal++)
		   {
		     if(i_beam_rec>cal_start[i_cal-1] && i_beam_rec<cal_start[i_cal])
		       {
			 next_cal_start =i_cal;
			 next_cal_start_set = true;
			 break;
		       }
		   }
	       }
	     
	     if(last_cal_end_set && next_cal_start_set)
	       {
		 //-------------------------------------------
		 //if both calibration end and start are set
		 //normal case: calibration data avail at both ends
		 // cal_end index should be smaller than cal_start index: always true!!!1
		 //------------------------------------------------
		 if(!(cal_end[last_cal_end]<i_beam_rec &&
		      i_beam_rec<cal_start[next_cal_start])) 
		   ErrorMessage("beam record is not lar").throwMe();
		 alpha = Gsys_cal_start[next_cal_start]-Gsys_cal_end[last_cal_end];
		 alpha /= Time_cal_start[next_cal_start] - Time_cal_end[last_cal_end];
		 Gsys= alpha *(beam_rel_time[i_beam_rec] - Time_cal_end[last_cal_end]);
		 Gsys += Gsys_cal_end[last_cal_end];
		 system_gain(i_beam_rec) = Gsys;	  		  
		 Ta = beam_ncnt_radio[i_beam_rec]/Gsys- Tcfg_sys;
		 antenna_temperature(i_beam_rec)=Ta ;
	       }
	     else if( (last_cal_end_set && !next_cal_start_set))
	       {//use last gain as constant
		 cout<< "sab counter "<<sab_counter_index(i_beam_rec)<<endl;
		 cout<<"only cal end set "<<endl;
		 Gsys=Gsys_cal_start[last_cal_end];
		 Ta =  beam_ncnt_radio[i_beam_rec]/Gsys- Tcfg_sys;
		 system_gain(i_beam_rec)=Gsys;
		 antenna_temperature(i_beam_rec)=Ta;
		 cout<<Gsys<<" "<<beam_ncnt_radio[i_beam_rec]<<" "<<Ta<<endl;
	       } 
	     else if (!last_cal_end_set && next_cal_start_set)
	       {//use next cal start gain as constant
		 cout<< "sab counter "<<sab_counter_index(i_beam_rec)<<endl;
		 cout<<"only cal start is set "<<endl;
		 Gsys=Gsys_cal_start[next_cal_start];
		 Ta =  beam_ncnt_radio[i_beam_rec]/Gsys- Tcfg_sys;
		 system_gain(i_beam_rec)=Gsys;
		 antenna_temperature(i_beam_rec)=Ta;
		 cout<<Gsys<<" "<<beam_ncnt_radio[i_beam_rec]<<" "<<Ta<<endl;
	       }
	     else
	       {//no calibration
		 ErrorMessage("Should not happen: No calibration data avail").throwMe();
	       }
	   }//pointing target with calibration avail
	 else
	   {
	     //should not happen in calibration mode
	     ErrorMessage("should not happen, no cal data at all").throwMe();
	   }
       }
     
     if(antenna_temperature.size()!=0)
       {
	 Plot a;
	 
	 a.addXY(sab_counter_index,"",antenna_temperature,"K",line("solid","red",1),sym("none"));
	 a.setTitle("antenna temperature vs sab counter:beam" +toStr(beam_id));
	 a.setXlabelandUnits("sab number");
	 a.setYlabelandUnits("antenna  temperature");
	 a.setFile(prefix_str +"antenna.eps");
	 a.setCmdFile(prefix_str +"antenna.xmgr");
	 a.show(plot_target);
	       
	 Plot b;
	 b.addXY(sab_counter_index,"",system_gain,"1/(s K)",line("solid","red",1),sym("none"));
	 b.setTitle("system gain:beam" +toStr(beam_id));
	 b.setXlabelandUnits("sab number");
	 b.setYlabelandUnits("calibrated system gain");
	 b.setFile(prefix_str +"cal_system_gain.eps");
	 b.setCmdFile(prefix_str +"cal_system_gain.xmgr");
	 b.show(plot_target);
       }
  
     */


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


 
