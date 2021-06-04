List of executables:

sar_proc: Cassini SAR processor, operates on a LBDR echo data file to generate SAR imagery. 
 Usage: sar_proc config_file 
 Several of the more useful keywords in the config file are described as follows:
 spice_directory: path to NAIF files
 ckernel_directory: path to NAIF ckernel files
 L1B_A_filename: filename of LBDR (if LBDR is broken into two files the code looks for files with the same name
                with a .2 appended.)
 ras_filename: Filename of intermediate file generated from raw data needed by noise subtraction
 t_feed_filename:
 t_hga_filename:
 t_scwg_filename: Three ancillary files used in calibrating LBDR data
 intermediate_geometry_file: Geometry data file needed by sar_proc created by generate_geometry_file
 min_oneway_gaindB_wrt_peak: Gain cutoff for excluding data.
 signal_to_amb_ratiodB: Signal to ambiguity cutoff for excluding data.
 noise_equivalent_sigma0dB: Noise equivalent sigma0 threshold for excluding data.
 burst_ordered_ieblist: Filename of intermediate file generated from raw data needed by noise subtraction
 PRODUCER_FULL_NAME: Name of data producer (e.g. whoever runs the SAR processor and wants to take credit)
 SAR_PROC_SEGMENT_TYPE:  Specifies method for defining limit of data to process in this segment 
      e.g. TIME_FROM_CLOSEST_APPROACH, TIME_FROM_EPOCH, TIME_FROM_TRIGGER, ABSOLUTE_TIME, FULL, LATLON,
      or LONLAT
 SAR_PROC_START_TIME: Start time for segment when a temporal method is used.
 SAR_PROC_END_TIME: End time for segment when a temporal method is used.
 SINGLE_BEAM_MODE: 0 indicates process all beams, 1 indicates only a single beam is processed
 SARPROC_BEAM_NUMBER: Number of beam to process when in single beam mode
 sar_annotation_file: Filename of text output file that provides information about run
 FULL_USABILITY_CALC_ENABLED: 1 is nominal method for selecting valid data bounds in each single look image using
               gain threshold, signal to ambigutiy threshold and noise equivalent sigma0 threshold
               0: turns off the nominal method
 GAIN_CUTOFF_ENABLED: 1 is an alternative method for selecting valid data bounds that only use antenna gain
              0: turns off alternative method
 GAIN_CUT_POLYNOMIAL_ORDER: Polynomial order of Gain cutoff threshold as a function of time (usually 0)
 GAIN_CUTOFF_BEAM_1_POLYNOMIAL_COEFF_0:  Constant coefficient of Gain cutoff polynomial for Beam 1
 BIDR_PIXELS_PER_DEGREE: Resoution of output images.
 SAR_PROC_NUMBER_OF_LOOKS_REQUIRED: Only pixels with this many looks or more included in SAR imagery output
 BIDR_FILENAME_PREFIX: Prefix for output SAR imagery files. A suffix is added for each backplane.
 (For more informtion on keywrods see '#' comments in example configuration file) 

preproc: preprocessor which generates LBDR file from EGSE file
decode_ccsd: decoder that generates EGSE file from raw telemetry SFDU file.
generate_geometry_file: routine to generate intermediate geometry file needed to make sar_proc computationally
                        efficient
average_down_bidr: Average down a BIDR to create a lower resolution SAR image
subsample_bidr: Better but more complicated routine to do the same thing as average_down_bidr and handle different backplane types appropriately. 
bidr2jpeg: Make a BIDR image into a JPEG
bidr_latlon: Take a row and column and return the latitude and longitude in a BIDR.
bidr_rowcol: Take latitude and longitude and return the row and column in a BIDR.
burst_data_extract: extract desired parameter from a LBDR or SBDR file
clip_header: Remove attached header from a BIDR, LBDR or SBDR file
flip_bidr: Routine to correct the format of a BIDR if it was processed from data in which the oblique cylindrical longitude is decreasing.
get_tiepoint_info: Get information about a tiepoint location from an SBDR file that is useful for usings those tiepoints to determine
                   Titan's rotation model
make_extra_backplane: Make a BIDR backplane at the real aperture resolution from any parameter in an SBDR file
sar_topo: Run the SARTopo processor without the corrections for attitude error.
topomap: Make a topomap from SARTopo data.
topomap2jpeg: Convert the topomap to a JPEG
zero_run_float_zip: Simple but useful method for compressing large single look SAR data files.
