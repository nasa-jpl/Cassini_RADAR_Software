// SAR stereo processor which generates backplanes needed for SAR stereo processing
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"
#include "L1I.h"
#include "L1B.h"
#include "config_keywords.h"
#include "DebugInfo.h"
#include "SARFunctions.h"
#include "SARProcParams.h"
#include "RangeDopplerFile.h"
#include "BIDR.h"
#include "RasList.h"

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
 string input;
try
  {

  //------------------------
  // Parse the command line
  //------------------------  

    const char* command = argv[0];
    if (argc != 2 && argc!=3 )
      {
	cerr << "Usage: " << command << " config_file" 
	     << "[single_burst_sab]" << endl;
	exit(-1);
      }


    int clidx = 1;
    const char* config_file = argv[clidx++];
    bool single_burst_mode = false;
    int single_burst_sab=0;
    if(argc==3){
      single_burst_mode=true;
      single_burst_sab=atoi(argv[clidx++]);
      cerr << "Operating sar_proc in single burst mode ......." << endl;
    }

    //------------------------------
    // Load configuration parameters
    //------------------------------
    Uvar::setMode("no_unit_support");
    Config cfg("sar_stereo",config_file);

    //------------------------------------------------------------------
    // Load static map for DebugInfo class which stores which routines 
    // (methods) have
    // debug on for what level with what debug parameters, etc
    //------------------------------------------------------------------
    DebugInfo::config(cfg); //INCOMPLETE

    // set up debugging info for main routine
    DebugInfo dbg("main");
 
    dbg.file << "sar_proc: main routine debug level is "<< dbg.level << endl;


    //---------------------------------
    // Load spice files & intermediate geometry file
    //---------------------------------
    Frame::setGeometryMode(Frame::READ_FILE);
    bool needs_ckernel=true;
    Frame::config(cfg,needs_ckernel); 
    BurstData::config(cfg);

    //------------------------------------
    // initialize input file
    //------------------------------------
    string l1b_filename = cfg.str(L1B_ACTIVE_MODE_FILENAME);
    L1B l1b(l1b_filename,"r+","active");
    l1b.readHeader();

    // Configure SarProcessorParams
    // includes Control Booleans, Processing constants and Calibration
    // coefficients
    SARProcParams spp(cfg,l1b.maxSabCounter()); 

    RasList* raslist=NULL;
    
    if(spp.simBaqForNoise || spp.noise_subtraction_on){
      string ras_file = cfg.str("ras_filename");
      raslist= new RasList(ras_file,"r");
    }

    // Initialize Output Classes

    string l1i_filename = cfg.str("L1I_filename");
    bool write_l1i = true;
    if(l1i_filename == "none" || l1i_filename == "NONE" || l1i_filename == "None"){
      l1i_filename="dummy.l1i";
      write_l1i=false;
    }
    L1I l1i(l1i_filename,"wb");
    l1i.config(cfg);

    // check for special error condition (trying to write ambiguity ratio to X)
    // when ambiguity computation is disabled
    if(l1i.replace_X_param==L1I::AMBIG_RATIO 
       && !spp.full_useable_calc_enabled){
      cerr << "Warning replace_X_param is AMBIG_RATIO but ambiguities" << endl;
      cerr << "         not computed. Reverting to standard X value." << endl;
      l1i.replace_X_param=L1I::NONE;
      l1i.replace_X_mode=false;
    }

    
    // writer L1I header
    if(write_l1i) l1i.writeHeader(l1b); 

    // get BIDR filename prefix (each type has adifferent extension)
    // (eventually filenames need to follow PDS convention)
    // I may do this with a script.
    string bidr_prefix = cfg.str("BIDR_FILENAME_PREFIX");
    BIDR bidr(bidr_prefix,cfg,l1b,spp,"w",STEREO);
 
    // Class for manipulating BIDR files 
    // handles beam mask, range, doppler, and backscatter
    // images simultaneously
    // set up BIDR projection and write BIDR headers

    bidr.writeHeaders(); 


    // open annotation file
    string afile=cfg.str("sar_annotation_file");

    ofstream afs;
    afs.open(afile.c_str());
    if(afs==NULL){
      ErrorMessage e("sar_proc: Unable to create Annotation File "+afile);
      e.throwMe();
    }
    
    afs.precision(5);
    afs << "!! Cassini RADAR Sar Stereo Processor Version 0.1" << endl;
    afs << "!! Burst by Burst Annotation File" << endl << endl;
 
    // Output single beam mode comment 
    if(spp.single_beam_mode){
      afs << "!! Proceesing in Single Beam mode." << endl;
      afs << "!! Only Beam #" << spp.beam_number << " bursts are processed."
	  << endl << endl;
    }

    DebugInfo::report(afs); // documents special debugging modes if enabled

    if(dbg.level > 0) 
      dbg.file << "Initializaton Complete. Starting main loop ..." << endl;
     


    // Output Labels for standard info line
    cout << "SAB_COUNTER  tfc lat lon beam s0 X num_pulses_received sampmin sampmax linemin linemax"
	 << endl;
    cout << "--------------------------------------------------------------------------------------"
         << endl;

    // Main SAR burst by burst loop
    while(!l1b.eof()){
      l1b.readParameter("sab_counter");

      // Nominal skipping of bursts (also implemented in BIDR::BIDR)
      // to determine size of swath
      if(spp.skipBurst(l1b.sab_counter)){
	l1b.skipRecord(); // this is done to avoid reading the whole record
                          // unless the burst is processed
	continue;
      }

      // Special Skipping of Bursts performed in single burst mode 
      // 
      if(single_burst_mode && l1b.sab_counter!=(unsigned int)single_burst_sab){
	l1b.skipRecord(); // this is done to avoid reading the whole record
                          // unless the burst is processed
	continue;
      }

      // zero L1I arrays
      l1i.zeroArrays();

      l1b.copyTo(l1i); // reads L1B record and copies to L1I
                       // needs to call L1I.resetRecord() first
                       // set l1i.rms_radar_data
      
      // determine which bursts to skip using precomputed map

  
      // Compute burst_no for debugging purposes                 
      int burst_no=(int)(l1i.record_id%1000000);
      
      afs << endl << "Starting Burst Number " << burst_no << " ..." << endl;

  
      // get beam number
      int bn=l1i.beam_number;

      // computes all burst parameters and converts to doubles
      // includes anti-aliasing filter start frequency and 
      // nominal delay and doppler centroid
      // and windowing parameters

      afs << " Beam number is " << bn << endl;
      l1i.computeBurstParams(spp,afs,raslist); // use L1I class to contain all values
      
      // If no pulse found write out a bad record and continue
      if(l1i.num_pulses_received==0) 
	{
          cerr << "Warning: SAB number " << l1i.sab_counter 
	       << "sar_proc says no pulses received." << endl
               << "even though preprocessor found pulses!!!!" << endl;

 	  afs << "Error: sar_proc says no pulses received; cannot process this burst" << endl
	     << "even though preprocessor found pulses!!!!" << endl;
	  continue;
	}

      // Process burst
      l1i.rangeCompressHensley(spp,afs); 

      l1i.azimuthCompress(spp,afs);   


      // Initialize standard info output line
      // L1I::calibrate and BIDR::stripProcessAndWrite populate the line
      spp.cout_line="";


      // spp.calibration_on==false case is handled within the L1I::calibrate 
      // routine
      // Useable extent is computed within L1I::calibrate  
      l1i.calibrate(spp,afs);      


      // estimate noise subtraction parameters if noise subtraction enabled
      if(spp.noise_subtraction_on)
	spp.estimateNoiseSubtractionParameters(l1i,raslist);

      
      afs << "End Burst " << burst_no << endl << endl;

      
      l1i.initLatLonTransformation(spp,bidr.getProjection());
        // sets up a simple conversion
	// from (fdop,range) to projected (OblCyl) (lat,lon) coordinates
        // that only applies to current burst

      if(write_l1i) l1i.writeRecord(); // output to L1I file
      
      if(spp.strip_processor_on){

	// Inserts calibrated burst and 
	// normalizes and writes
	// a line of the BIDR if it is ready
	// this routine also handles computation of and
	// writes to auxiliary BIDR files
	bidr.stripProcessAndWrite(l1i,spp,single_burst_mode); 
      }

      // output standard info line
      cout << spp.cout_line << endl;
    }

   // Output last portions of BIDR images.
   if(spp.strip_processor_on) bidr.flush(spp);  

   // Update the PDS label in each output product
   if(spp.strip_processor_on) bidr.rewriteHeaders();

   // output burst skip report to annotation file
   spp.reportSkippedBursts(afs);


   if(raslist) delete raslist;
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



