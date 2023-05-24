%%% load .mat file with the "current" and "final" spin model structures
%%% see get_residual_for_case.m for details of the parameters in each
%%% structure
load titan_spin_models_20190330.mat

%%% List of tiepoint numbers that failed Bruce Bills' quality control check
bruce_outliers=[211,233,277,278,347,480,481,723,844,845,1127,2270, ...
    2272,2303,2305,2306,2312,2314,2315,2317,2320,2321,243,244,245, ...
    246,247,248,249,250,251,252,257,258,259,352,353,354,355,2412];  

%%% compute all information for "good" tiepoints with the "current" Titan
%%% spin model
[resid.current,tiepoints.current_good]=get_residual_for_case(current,0,bruce_outliers)

%%% compute all information for "bad" tiepoints with the "current" Titan
%%% spin model
usebadtiesonly=1;
[bad.current,tiepoints.current_bad]=get_residual_for_case(current,usebadtiesonly,bruce_outliers)


%%%% The "final" titan spin model is the same as "dynamic_spin_and_pole" 
%%%% Tiepoints information computed with this spin model can be found in
%%%% the tiepoint files with "dsap" in the name
dynamic_spin_and_pole=final;

%%%% compute all information for "good" tiepoints with the "dynamic_spin_and_pole" Titan spin model
[resid.dynamic_spin_and_pole,tiepoints.dynamic_spin_and_pole_good]=get_residual_for_case(dynamic_spin_and_pole,0,bruce_outliers)

%%%% compute all information for "bad" tiepoints with the "dynamic_spin_and_pole" Titan spin model
usebadtiesonly=1;
[bad.dynamic_spin_and_pole,tiepoints.dynamic_spin_and_pole_bad]=get_residual_for_case(dynamic_spin_and_pole,usebadtiesonly,bruce_outliers)

%%% Uncomment this line if you want to save tiepoint info in a .mat file
%%%save tiepoints.mat tiepoints dynamic_spin_and_pole current

%%% List of output filenames
%%% Revision 20230523 fixes a bug in the reported look direction and wavelength
 outfile_gm = 'tiepoints_good_measured_20230523.csv'
 outfile_bm = 'tiepoints_bad_measured_20230523.csv'
 outfile_gc = 'tiepoints_good_current_pos_20230523.csv'
 outfile_bc = 'tiepoints_bad_current_pos_20230523.csv'
 outfile_gd = 'tiepoints_good_dsap_pos_20230523.csv'
 outfile_bd = 'tiepoints_bad_dsap_pos_20230523.csv'  
 outfile_gpls = 'tiepoints_good_dsap_onsphere_20230523.csv'
 outfile_bpls = 'tiepoints_bad_dsap_onsphere_20230523.csv'
 outfile_gplsc = 'tiepoints_good_current_onsphere_20230523.csv'
 outfile_bplsc = 'tiepoints_bad_current_onsphere_20230523.csv' 

 %%% Write _measured_ files which contain fundamental measured info for
 %%% each tiepoint
 outfields={'height','time_1','doppler_1','range_1','scpos_inert_1','scvel_inert_1','lookvec_inert_1','wavelength_1', 'time_2','doppler_2','range_2','scpos_inert_2','scvel_inert_2','lookvec_inert_2','wavelength_2'};
fmt='%28s,%7.3f,%14.3f,%11.2f,%11.2f,%11.2f,%11.2f,%11.2f,%10.6f,%10.6f,%10.6f,%11.2f,%11.2f,%11.2f,%16.11g,%14.3f,%11.2f,%11.2f,%11.2f,%11.2f,%11.2f,%10.6f,%10.6f,%10.6f,%11.2f,%11.2f,%11.2f,%16.11g\n';
write_tie_file(outfile_gm,tiepoints.current_good,outfields,fmt); % good tiepoint file
write_tie_file(outfile_bm,tiepoints.current_bad,outfields,fmt);% bad tiepoint file


%%% Write _current_ files which contain tiepoint positions computed using
%%% the "current" spin model which was used for the latest PDS archived
%%% version of the Cassini SAR data. However locations differ because the
%%% archived data were geolocated assuming all pixels were on the 2575-km
%%% reference sphere. Using the SARTopo heights instead achieves much
%%% smaller( better) residual position between different observations of the
%%% same tiepoint.
outfields={'height','time_1','pos_inert_1','lat_1','wlon_1','time_2','pos_inert_2','lat_2','wlon_2','residual_dist'};
fmt='%28s,%7.3f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%11.3f\n';
write_tie_file(outfile_gc,tiepoints.current_good,outfields,fmt); % good tiepoint file
write_tie_file(outfile_bc,tiepoints.current_bad,outfields,fmt);% bad tiepoint file

%%% Write _dsap_ files which contain tiepoint positions computed using
%%% the "dynamic_spin_and_pole" spin model which is a better fit that
%%% allows spin rate and pole to vary with time and thus has smaller
%%% residual position differences. 
write_tie_file(outfile_gd,tiepoints.dynamic_spin_and_pole_good,outfields,fmt); % good tiepoint file
write_tie_file(outfile_bd,tiepoints.dynamic_spin_and_pole_bad,outfields,fmt); % bad tiepoint file

%%%% Write _dsap_onsphere_ files which contain tiepoint locations computed
%%%% using the "dynamic_spin_and_pole" spin model but assuming every
%%%% tiepoint is on the 2575.0 sphere.
outfields={'height','time_1','pos_inert_zeroh_1','lat_zeroh_1','wlon_zeroh_1','time_2','pos_inert_zeroh_2','lat_zeroh_2','wlon_zeroh_2','residual_dist_zeroh'};
fmt='%28s,%7.3f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%11.3f\n';
write_tie_file(outfile_gpls,tiepoints.dynamic_spin_and_pole_good,outfields,fmt)
write_tie_file(outfile_bpls,tiepoints.dynamic_spin_and_pole_bad,outfields,fmt)

%%%% Write _current_onsphere_ files which contain tiepoint locations computed
%%%% using the "current" spin model but assuming every
%%%% tiepoint is on the 2575.0 sphere. This yields positions that are
%%%% identical to the geolocations of tiepoints in the PDS archive. 
outfields={'height','time_1','pos_inert_zeroh_1','lat_zeroh_1','wlon_zeroh_1','time_2','pos_inert_zeroh_2','lat_zeroh_2','wlon_zeroh_2','residual_dist_zeroh'};
fmt='%28s,%7.3f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%14.3f,%11.2f,%11.2f,%11.2f,%8.4f,%9.4f,%11.3f\n';
write_tie_file(outfile_gplsc,tiepoints.current_good,outfields,fmt)
write_tie_file(outfile_bplsc,tiepoints.current_bad,outfields,fmt)
