Cassini  RADAR Tiepoint Description

Bryan Stiles, Randy Kirk, Alex Hayes, and Bruce Bills

This document describes the tiepoint files in the Cassini RADAR software github archive https://github.com/nasa-jpl/Cassini_RADAR_Software/ in the TiePoints subdirectory.   The data are a set of pairs of observations of the same surface features on Titan at different times in Cassini RADAR SAR images. Each such pair is referred to as a “tiepoint”. Errors in the rotational model of Titan lead to differences in the apparent latitudes and longitudes of the same feature as observed by Cassini RADAR at different times. The tiepoints archived here were used by the Cassini RADAR Team to determine the rotation axis and pole of Titan used to locate the Cassini SAR images [1,2]. 

The following 5 tiepoint files are included:

1. tiepoints_good_measured_20230523.csv: quantities measured from the radar such as Doppler and range rather than quantities like latitude and longitude that depend upon the Titan rotation model and Titan reference surface. 
2. tiepoints_good_current_onsphere_20230523.csv: latitudes and longitudes computed for each tiepoint as computed using the current Titan rotation model and Titan surface shape, sphere of radius 2575.0 that were used to locate the Cassini SAR images that are archived in PDS. The root mean square distance between observations is 3.1 km for all 2421 good tiepoints.
3. tiepoints_good_current_pos_20230523.csv: latitudes and longitudes computed using the current Titan rotation model and SARTopo [3] surface heights. The root mean square distance between observations is 2.4 km for all 2421 good tiepoints, significantly smaller than when assuming tiepoints are on the reference sphere.
4. tiepoints_good_dsap_onsphere_20230523.csv: latitudes and longitudes computed with a dynamically varying Titan spin rate and axis model fit using physical constraints as described in [1] and assuming tiepoints are on the 2575.0 radius Titan reference sphere. The root mean square distance between observations is 2.6 km for all 2421 good tiepoints, which demonstrates that getting the surface height correct is more important for reducing the residual geolocation error than further optimizing the spin model.
5. tiepoints_good_dsap_pos_20230523.csv: latitudes and longitudes computed using the best fit dynamically varying Titan spin rate and axis model [1] and SARTopo surface heights. The root mean square distance between observations is 1.6 km for all 2421 good tiepoints, a significantly lower residual error than the spin model used for the archived data set.

In addition to these 5 files which document high quality tiepoints that passed quality control, there are five other files (e.g. tiepoints_bad_measured_20230523.csv) that document other tiepoints which were excluded when quality control was applied. There are 2,427 high quality tiepoints and 56 tiepoints with bad quality. Bad quality tiepoints were not used to fit Titan spin models. The high quality tiepoints have residual apparent dislocation between observations in the archived SAR imagery between 0 and 14 km. The bad tiepoints have dislocations between 0 and 2,343 km. The bad tiepoints with the smallest dislocations are cases that made use of lower resolution SAR imagery, observations with very small temporal differences, or observations that were far removed from the nearest available SARTopo surface height information.

The tiepoint files are comma separated files in which each row is a single tiepoint and each column is a particular quantity. Different files have different columns. The headings for each column that appears in one or more of the files are defined as follows:

* Tiepoint_ID:  A unique identifier for each tiepoint. The format of the identifier varies but all identifiers container a string that indicate the affiliation or name of the investigator who first identified the tiepoint and the Titan flyby and segment numbers of the two observations in chronological order.
* height: the nearest available SARTopo height to the location in the archived SAR imagery in km above the 2575.0 km radius reference sphere. We used the  SARTopo obtained during the second observation because SARTopo was deemed to be more error prone earlier in the mission. 
* time_1: time of first observation in seconds since J2000 (Jan 1, 2000, 12:00:00 TT).
* doppler_1: Doppler of first observation of tiepoint location in Hz, a measured quantity based upon the observed frequency shift of the radar signal. An image (backplane) of Doppler can be optionally produced by the Cassini SAR processor that is coregistered with the SAR imagery.
* range_1: Range in km of first observation from Cassini antenna to tiepoint location in km, a measured quantity based upon the observed round trip time of the radar signal. An image (backplane) of Range can be optionally produced by the Cassini SAR processor that is coregistered with the SAR imagery.
* sc_pos_inert_1_{x,y,z}:  Titan-centered Cartesian coordinates in km of the inertial (J2000) position of the Cassini spacecraft at the time of the first observation. 
* sc_vel_inert_1_{x,y,z}:  Titan-centered Cartesian coordinates in km/s of the inertial (J2000) velocity of the Cassini spacecraft at the time of the first observation. 
* lookvec_inert_1{x,y,z}: The vector from the spacecraft to the peak of the antenna pattern on the ground used to determine the gross pointing direction in order to choose among the two possible solutions for the tiepoint location. The vector is in inertial coordinates. Units are km. 
* wavelength_1: Wavelength of transmitted signal for the first observation. This value is used to convert Doppler to the vertex angle of the Doppler cone. 
* pos_inert_1_{x,y,z}: Titan-centered Cartesian coordinates of the tiepoint from the Cassini Spacecraft at the time of the first observation.  The tiepoint position in inertial space is computed from height, range, Doppler, and spacecraft position and velocity, by intersecting three surfaces: 1)  a sphere of radius 2575.0 km + the measured SARTopo height centered on the center of Titan, 2) a cone with its point at the spacecraft location oriented along the spacecraft velocity vector and its vertex angle such that points on its surface have the measured Doppler, and 3) a sphere centered on the spacecraft position with a radius of the measured range. Intersecting these spheres yields two solutions, we choose the one which is closest to the antenna pointing direction.
* lat_1: latitude in degrees of Titan body fixed (TBF) position for the first observation of the tiepoint computed by multiplying the inertial coordinate by the rotation matrix from inertial to TBF at time_1. Different rotation models are used for different files. Files with “current” in the name use the rotation model of the Cassini SAR PDS archived data set which is a synchronous constant spin rate (Titan tidally locked with Saturn), constant spin pole model. Files with “dsap” in the name use a titan rotation model with a varying spin rate and pole that better fits the tiepoint data. The “dsap” (dynamic spin and pole) model  is still synchronous on long time scales. The subSaturn point on Titan varies by a couple degrees due to Titan’s orbit not being circular.
* wlon_1: positive west longitude in degrees of TBF position for the first observation of the tiepoint.
* pos_inert_zeroh_1_{x,y,z}: same as pos_inert_1_{x,y,z} except that the tiepoint location is computed by assuming it is on the 2575-km radius reference sphere. This value is found in the files with “onsphere” in the name.
* lat_zeroh_1: same as lat_1 except that the tiepoint location is computed by assuming it is on the 2575-km radius reference sphere. This value is found in the files with “onsphere” in the name.
* wlon_zeroh_1: same as wlon_1 except that the tiepoint location is computed by assuming it is on the 2575-km radius reference sphere. This value is found in the files with “onsphere” in the name.
* time_2, doppler_2 ,range_2, sc_pos_inert_2_{x,y,z}, sc_vel_inert_2_{x,y,z}, lookvec_inert_2_{x,y,z}, wavelength_2, pos_inert_2_{x,y,z},lat_2,wlon_2, pos_inert_zeroh_2_{x,y,z}, lat_zeroh_2, wlon_zeroh_2: Same as above except for the second observation of the tiepoint rather than the first.
* residual_dist: apparent distance in km between TBF positions for the two observations of the tiepoint.

In additions to tiepoint files, this directory also contains other files including MATLAB code and data files used to compute the tiepoint archive, to convert from inertial to titan body fixed coordinates for various spin models, and to plot tiepoints. The following is a description of each of those files. Researchers can use this code to generate new tiepoints and compute residual tiepoint location errors for a given candidate Titan spin model. In this way, they can utilize Cassini data to refine the estimates of Titan’s pole and spin rate.

* TiePointREADME.docx:  MS Word version of this file.
* TiePointREADME.txt: ASCII version of this file.
* jplusgscornell_2016_newsartopo_randop.tab:  This is the raw tiepoint file used to generate the tiepoint archive. Each row of this file was generated by twice running the get_tiepoint_info executable compiled from src/programs/get_tiepoint_info.cpp. The executable needs to be run separately for each observation of a tiepoint. This file was made by concatenating data from multiple runs of get_tiepoint_info and handediting it to remove changes in the format such as the addition of ground impact time and SAB number. Since the output of get_tiepoint_info is in inertial coordinates, the values in this file are insensitive to the particular version of the Titan spin model used to generate the tiepoint data. The usage of get_tiepoint_info is: 
o get_tiepoint_info cfg_file lbdrfile bidrfile(should be beammask) topomapfile line sample
* titan_spin_models_20190330.mat:  MATLAB data file containing the parameters for the “current” and “dsap” Titan spin models.
* tqual_2016.txt: ASCII file containing the quality value of each tiepoint. Between 0 and 0.5 is good. Everything else is bad.
* make_tiepoint1_images.m: example MATLAB code to make images of each observation of a tiepoint.
* TiepointImage_tas01t23s01_jpl2_p1_TA.png: Image of tiepoint tas01t23s01_jpl2_p1 as seen in the TA flyby generated by make_tiepoint1_images.m.
* TiepointImage_tas01t23s01_jpl2_p1_T23.png: Image of tiepoint tas01t23s01_jpl2_p1 as seen in the T23 flyby generated by make_tiepoint1_images.m.
* read_bidr.m: MATLAB function used to read a BIDR image
* make_tiepoint_archive_csv.m:  Top level MATLAB script used to generate the tiepoint archive.
* get_residual_for_case.m: MATLAB function which generates a tiepoint information structure given a spin model, and quality control options.
* get_j2000_pos.m:  MATLAB function which compute tiepoint inertial coordinates given  spacecraft position and velocity, Titan’s spin rate, radar wavelength, doppler, range, the center of the radar footprint in inertial coordinates and the local radius of Titan (including surface height estimate).
* ellsphcone_intersect.m: MATLAB function called by get_j2000_pos.m find the intersection of two spheres and a cone.
* pos_to_lonlat.m: MATLAB function to compute longitude and latitude from a Titan body fixed position vector.
* read_dopran_old.m: MATLAB function to read raw tiepoint info file.
* write_tie_file.m: MATLAB function used to write lines to tiepoint archive files.

Acknowledgements:
Tiepoint location (line and sample) within Cassini RADAR images were determined manually in a labor intensive manner by Alexander Hayes at Cornell University, Randolph Kirk and Ella Lee at USGS, and Bryan Stiles at JPL.  Tiepoint quality control was performed by Bryan Stiles and Bruce Bills at JPL. Titan spin models from these tiepoints were estimated by Bruce Bills and Bryan Stiles.

References:

[1] Bills, Bruce, Bryan W. Stiles, and Alexander Hayes. "Constraints on Titan's rotation from Cassini mission radar data." AAS/Division of Dynamical Astronomy Meeting# 46. Vol. 46. 2015.

[2] Bryan W Stiles, Randolph L Kirk, Ralph D Lorenz, Scott Hensley, Ella Lee, Steven J Ostro, Michael D Allison, Philip S Callahan, Yonggyu Gim, Luciano Iess, Paolo Perci del Marmo, Gary Hamilton, William TK Johnson, Richard D West, and the Cassini RADAR Team, “Determining Titan's spin state from Cassini RADAR images,” The Astronomical Journal, Volume 135, No. 5, pp 1669-1680  March, 2008.

[3] Bryan W Stiles, Scott Hensley, Yonggyu Gim, David M Bates, Randolph L Kirk, Alex Hayes, Jani Radebaugh, Ralph D Lorenz, Karl L Mitchell, Philip S Callahan, Howard Zebker, William TK Johnson, Stephen D Wall, Jonathan I Lunine, Charles A Wood, Michael Janssen, Frederic Pelletier, Richard D West, Chandini Veeramacheneni, snd the Cassini RADAR Team, “Determining Titan surface topography from Cassini SAR data,” Icarus, Volume 202, Issue 2, August 2009, Pages 584-598, https://doi.org/10.1016/j.icarus.2009.03.032






