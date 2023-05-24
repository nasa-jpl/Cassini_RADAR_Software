This file describes how to build the Cassini RADAR software
and run an example data case with the SAR processor.
Note:: This code has only been tested on Linux systems.
It was most recently compiled on CentOS Linux 7 (Core)
using gcc version 4.8.5 20150623 (Red Hat 4.8.5-39) 
 
Step 1: Compile the cspice libraries used by the software.
(This is the version which was used by the Cassini RADAR code without any commit info.
More up-to-date versions can be found at https://naif.jpl.nasa.gov )

>> cd cspice
>> mkdir exe
>> mkdir lib
>> makeall.csh

Step 2: From the top directory compile the rest of the Cassini software
>> make tree
>> make

If step 2  fails on one of the executables, you can still compile the 
rest one at a time by

>> cd src/programs
>> make your_favorite_executable (e.g. >> make sar_proc )

See src/programs/README.TXT for full list of executables with brief description. 
The command line usage for each executable can be determined by running it without any command line parameter for example:

>> sar_proc

yields


>>Usage: ../src/programs/sar_proc config_file[single_burst_sab]

So, the sar_proc executable takes a configuration filename and optionally a SAB (burst identifier) if ones wants to run the 
processor on a single burst.


Step3: Acquire Cassini data.  


To run the SAR processor you will need NAIF geometry data available from naif.jpl.nasa.gov and LBDR, beam pattern, and 
configuration file data available from pds.jpl.nasa.gov. Sufficient data from T108 that can be used as a test is in the 
repository in the inputs directory. At present the SAR processor can only be run with noise subtraction turned off unless raw telemetry is available. 
Noise subtraction requires raw data that may be found at:https://zenodo.org/record/5636891#.YYArPyV7ltI 


Step 4: Create configuration files. The archived Cassini configuration file needs to be modified to find the input data on 
your computer. Three modified configuration files are provided to show how this can be done.
The original T108 primary segment configuration file is inputs/cfg/PCF_0265_V02.CFG.
There are five other config files in the test-ref directory. They have input filenames with relative paths and can be 
used to run routines from the test subdirectory. Use "diff" to see how they were modified from the archived config file:

sar_proc_no_rawdata.cfg:  Configuration file for running SAR processor without raw data or noise subtraction.
generate_geometry_file.cfg:  Configuration file to create geometry file needed by SAR processor
sar_proc.cfg: Configuration file to run SAR processor with noise subtraction when raw data is available
preproc.cfg: Configuration file used by preproc to recreate LBDR from EGSE file and make other intermediate files needed for
             noise subtraction
decode_ccsd.cfg: Configuration file used by decode_ccsd to generate EGSE file from raw data (SFDU) file.

 

Step 5: Run the example case from the LBDR data without raw data (and without noise subtraction)

Copy test-ref/sar_proc_no_rawdata.cfg to the test directory
Copy test-ref/generate_geometry_file.cfg to the test directory
Unzip inputs/data/LBDR_14_D265_P1_V02.ZIP to get LBDR_14_D265_P1_V02.TAB
Copy or link LBDR_14_D265_P1_V02.TAB to test-ref/t108_reproc1.lbdr
Unzip inputs/data/LBDR_14_D265_P2_V02.ZIP to get LBDR_14_D265_P2_V02.TAB
Copy or link LBDR_14_D265_P2_V02.TAB to test-ref/t108_reproc1.lbdr.2

Run the following in the test-ref directory after setting your PATH to include src/programs
>> mkdir no_rawdata
>> generate_geometry_file generate_geometry_file.cfg
>> sar_proc sar_proc_no_rawdata.cfg 

This will generate a bunch of SAR data files with names like t108_reproc1_corramb_bidr.*  
The file  no_rawdata/t108_reproc1_corramb_bidr.s0_uncorr should be identical to 
test-ref/no_rawdata/t108_reproc1_corramb_bidr.s0_uncorr


Step 6: Run the example case from raw data with noise subtraction
(You can skip running generate_geometry_file if you already did that in step 5.)
Copy test-ref/sar_proc.cfg to the test directory
Copy test-ref/generate_geometry_file.cfg to the test directory
Copy test-ref/preproc.cfg to the test directory
Copy test-ref/decode_ccsd.cfg to the test directory

Run the following in the test-ref directory after setting your PATH to include src/programs
>> decode_ccsd decode_ccsd.cfg none
>> preproc preproc.cfg none
>> generate_geometry_file generate_geometry_file.cfg
>> sar_proc sar_proc.cfg 

This will generate a bunch of SAR data files with names like t108_reproc1_corramb_bidr.*  
The file  t108_reproc1_corramb_bidr.s0nsqt_uncorr should be identical to 
test-ref/t108_reproc1_corramb_bidr.s0nsqt_uncorr

 
For more information contact Bryan.W.Stiles@jpl.nasa.gov.
