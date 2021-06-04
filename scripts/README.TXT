Brief description of files in this directory:
To get usage statement run routines without command line parameters, or read comments in file
All files in this directory are human readable.

cutrecs: executable Perl script to cut out fixed length records from a file
fix_bounds.csh: executable csh script to fix bounds in BIDR files using a precomputed text file containing 
                the correct bounds for a list of BIDRs.

fndr: executable csh script to print the full path name of the directory under the
      search_tree that matches the target_dir.

getbounds.csh: executable csh script that computes the correct bounding box for a BIDR file and outputs the 
               correct values followed by the possibly incorrect values in the header of the file

get_correct_bounds.csh: executable csh script that uses getbound.csh to loop through all BIDR files 
               in a specific directory structure and write their correct bounds to file usable 
               by fix_bound.csh

grabhead: executable perl script to copy a fixed length header from any file

inclist: executable csh script to generate a list of include files needed by a C++ source file

make_volume_tree: executable csh script used to make a Planetary Data System (PDS) Volume directory structure

pds_split3: executable csh script used to break up a large PDS volume into three gzipped tar files under the 2 GB limit

replace_ctrlM.cpp:  C++ source code for simple routine to remove ^M end-of-line characters created by Windows systems
		    To compile type "g++ -o replace_ctrlM replace_ctrlM.cpp"

replace_header: executable csh script to replace the header of a data file

substring.csh: executable csh script to replace a pattern of characters in a file with a different pattern of characters
               For example the command "substring.csh XXX YYY  *.txt" replaces "XXX" with "YYY" everywhere "XXX"
               occurs in every file that ends in ".txt"
