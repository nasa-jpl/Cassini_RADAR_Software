#include<time.h>
#include<stdio.h>

// POINT TARGET SIM TIMING
double run_time=0.0;
double pts_config_time=0.0;
double pts_run_time=0.0;
double pts_addecho_time=0.0;
double pts_addecho_getPTE_time=0.0;
double pts_getPTAN_time=0.0;
double pts_range_time=0.0;
double pts_ideal_track_time=0.0;
double pts_scale_time=0.0;
double pts_fdop_time=0.0;
double pts_locate_time=0.0;
double pts_l1b_write_time=0.0;
double pts_update_radar_params_time=0.0;
double pts_set_midpulse_time=0.0;

double frame_ephemeris_time=0.0;
double pv_framerotate_time=0.0;


//#define PROFILING_ON
#ifdef  PROFILING_ON
double get_time(){return((double)clock()/(double)CLOCKS_PER_SEC);}

void report_profile()
{

  double run_time_accounted_for=pts_run_time+pts_config_time;
  fprintf(stderr,"\n\nTOP LEVEL PROFILE\n");
  fprintf(stderr,"total_run_time = %g (%g percent)\n",run_time,100.0);
  fprintf(stderr,"PTS::config time = %g (%g percent)\n",pts_config_time,
	  100.*pts_config_time/run_time);  
  fprintf(stderr,"PTS::run time = %g (%g percent)\n",pts_run_time,
	  100.*pts_run_time/run_time);
  fprintf(stderr,"---------------------------------\n");
  fprintf(stderr,"Total run time accounted for = %g (%g percent)\n\n"
	  ,run_time_accounted_for,
	  100.*run_time_accounted_for/run_time);

  double pts_run_time_accounted_for=pts_addecho_time+pts_range_time+
    pts_ideal_track_time+pts_scale_time+pts_fdop_time+pts_locate_time +
    pts_l1b_write_time+pts_update_radar_params_time+pts_set_midpulse_time;
  fprintf(stderr,"\n\nLEVEL 2: PTS::RUN PROFILE\n");
  fprintf(stderr,"PTS::run time = %g (%g percent)\n",pts_run_time,
	  100.*pts_run_time/run_time);
  fprintf(stderr,"PTS::addecho time = %g (%g percent)\n",pts_addecho_time,
	  100.0*pts_addecho_time/run_time);
  fprintf(stderr,"PTS::idealtrack time = %g (%g percent)\n",
	  pts_ideal_track_time,
	  100.0*pts_ideal_track_time/run_time);
  fprintf(stderr,"scale time = %g (%g percent)\n",pts_scale_time,
	  100.0*pts_scale_time/run_time);
  fprintf(stderr,"fdop time = %g (%g percent)\n",pts_fdop_time,
	  100.0*pts_fdop_time/run_time);
  fprintf(stderr,"range time = %g (%g percent)\n",pts_range_time,
	  100.0*pts_range_time/run_time);
  fprintf(stderr,"L1B::locate time = %g (%g percent)\n",pts_locate_time,
	  100.0*pts_locate_time/run_time);
  fprintf(stderr,"L1B::write time = %g (%g percent)\n",pts_l1b_write_time,
	  100.0*pts_l1b_write_time/run_time);
  fprintf(stderr,"PTS::updateRadarParams time = %g (%g percent)\n",
	  pts_update_radar_params_time,
	  100.0*pts_update_radar_params_time/run_time);
  fprintf(stderr,"PTS::set_midpulse_time time = %g (%g percent)\n",
	  pts_set_midpulse_time,
	  100.0*pts_set_midpulse_time/run_time);
  fprintf(stderr,"---------------------------------\n");
  fprintf(stderr,"PTS::RUN time accounted for =  %g (%g percent)\n",
	  pts_run_time_accounted_for, 100.0*pts_run_time_accounted_for
	  /run_time);

  double pts_addecho_time_accounted_for=pts_addecho_getPTE_time+
  pts_getPTAN_time;

  fprintf(stderr,"\n\nLEVEL 3: PTS::ADDECHO PROFILE\n");
  fprintf(stderr,"PTS::addecho time = %g (%g percent)\n",pts_addecho_time,
	  100.0*pts_addecho_time/run_time);
  fprintf(stderr,"PTS::addecho->getPTE time = %g (%g percent)\n",
	  pts_addecho_getPTE_time,
	  100.0*pts_addecho_getPTE_time/run_time);
  fprintf(stderr,"PTS::addecho->getPTAN time = %g (%g percent)\n",
	  pts_getPTAN_time,
	  100.0*pts_getPTAN_time/run_time);
  fprintf(stderr,"---------------------------------\n");
  fprintf(stderr,"PTS::ADDECHO time accounted for =  %g (%g percent)\n",
	  pts_addecho_time_accounted_for, 100.0*pts_addecho_time_accounted_for
	  /run_time);


  fprintf(stderr,"\n\nINTERMEDIATE LEVEL GEOMETRY ROUTINES\n");
  fprintf(stderr,"Frame::ephemeris time = %g (%g percent)\n",
	  frame_ephemeris_time,
	  100.0*frame_ephemeris_time/run_time);
  fprintf(stderr,"PositionVector::frameRotate time = %g (%g percent)\n",
	  pv_framerotate_time,
	  100.0*pv_framerotate_time/run_time);
}

#else
double get_time(){return(0.0);}
void report_profile(){}
#endif



