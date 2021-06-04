#ifndef PROFILING_H
#define PROFILING_H
#include<time.h>
// POINT TARGET SIM TIMING
extern double run_time;
extern double pts_config_time;
extern double pts_run_time;
extern double pts_addecho_time;
extern double pts_addecho_getPTE_time;
extern double pts_getPTAN_time;
extern double pts_range_time;
extern double pts_ideal_track_time;
extern double pts_scale_time;
extern double pts_fdop_time;
extern double pts_locate_time;
extern double pts_l1b_write_time;
extern double pts_update_radar_params_time;
extern double pts_set_midpulse_time;
extern double frame_ephemeris_time;
extern double pv_framerotate_time;

double get_time();
void report_profile();
#endif



