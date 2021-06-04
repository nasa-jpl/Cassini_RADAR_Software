// SAR processor
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

#define OUTPUT_ALL_DATA 1
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
    if (argc != 2 && argc!=3 && argc!=4)
      {
	cerr << "Usage: " << command << " config_file" 
	     << "[input_intermed_filename]" << "noise_correction_factor" << endl;
	exit(-1);
      }


    int clidx = 1;
    const char* config_file = argv[clidx++];
    bool single_burst_mode = false;
    int single_burst_sab=0;
    string intermed_file;
    bool read_intermed=false;
    float noise_correction_factor=0.5;// default value (works much better than 1)
    if(argc>=3){
      intermed_file = argv[clidx++];
      read_intermed=true;
    }
    if(argc==4){
      noise_correction_factor=atof(argv[clidx++]);
    }

    //------------------------------
    // Load configuration parameters
    //------------------------------
    Uvar::setMode("no_unit_support");
    Config cfg("sar_topo",config_file);


    // hard code ambiguity ratio threshold for now
    float ambrat_thresh=0.1;

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

    string topopre=cfg.str("SAR_TOPO_PREFIX");
    float dhthresh=0.075;
    string auxfile= topopre + "_aux.dat";
    string heightfile= topopre + "_height.txt";
    string csvfile[5];
    csvfile[0]=topopre+"_1.csv";
    csvfile[1]=topopre+"_2.csv";
    csvfile[2]=topopre+"_3.csv";
    csvfile[3]=topopre+"_4.csv";
    csvfile[4]=topopre+"_5.csv";
    if(!read_intermed) intermed_file= topopre + "_intermed.dat";
    FILE* afp=fopen(auxfile,"w");
    FILE* hfp=fopen(heightfile,"w");
    FILE* csvfp[5];
    for(int d=0;d<5;d++) csvfp[d]=fopen(csvfile[d],"w");
    FILE* mfp=NULL;
    if(read_intermed) mfp=fopen(intermed_file,"r");
    else mfp=fopen(intermed_file,"w");

    // assign height loop parameters
    float min_trial_height=(float)cfg.getDouble("SAR_TOPO_MIN_HEIGHT");
    int num_heights=cfg.getInt("SAR_TOPO_NUM_HEIGHTS");
    float height_step=(float)cfg.getDouble("SAR_TOPO_HEIGHT_STEP");
    float oneway_gain_cutoff=(float)cfg.getDouble("SAR_TOPO_ONEWAY_GAIN_CUTOFF");
    int runave_wid=cfg.getInt("SAR_TOPO_RUN_AVERAGE_RADIUS");
    float azwidper=(float)cfg.getDouble("SAR_TOPO_AZIMUTH_WIDTH_PERCENT");
    int height_fit_order=1;

    // allocate intermediate parameter arrays

    
    // MAX_BIDR_LINES is defined in BIDR
    
    float*** dsigma0=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** rms_ds0=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** noisefloor1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** noisefloor2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** power1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** power2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** ambrat1=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** ambrat2=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** dsigma0std=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** dsigma0bias=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float** heightbias=(float**)make_array(sizeof(float),2,5,MAX_BIDR_LINES);
    float** height=(float**)make_array(sizeof(float),2,5,MAX_BIDR_LINES);
    float** dheight=(float**)make_array(sizeof(float),2,5,MAX_BIDR_LINES);
    float** dhdfnf=(float**)make_array(sizeof(float),2,5,MAX_BIDR_LINES);
    bool** monotonic=(bool**)make_array(sizeof(bool),2,5,MAX_BIDR_LINES);
    int** cmin=(int**)make_array(sizeof(int),2,5,MAX_BIDR_LINES);
    int** cmin2=(int**)make_array(sizeof(int),2,5,MAX_BIDR_LINES);
    int** nzerocross=(int**)make_array(sizeof(int),2,5,MAX_BIDR_LINES);
    float*** sumsigma0=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** sumsigma0andnoise=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);
    float*** wnumlook=(float***)make_array(sizeof(float),3,4,MAX_BIDR_LINES,num_heights);

    float*** heightpoly=(float***)make_array(sizeof(float),3,5,MAX_BIDR_LINES,
				   height_fit_order+1);

    int** overlap_col=(int**)make_array(sizeof(int),2,4,MAX_BIDR_LINES);

    int** overlap_width=(int**)make_array(sizeof(int),2,4,MAX_BIDR_LINES);
    float** inc1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    float** inc2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    float** lat=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
    float** lon=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);

    
  

    if(!read_intermed){    
  
    for(int c=-1;c<num_heights;c++){
      float trial_height=min_trial_height+c*height_step;
      if(c==-1) trial_height=0; // special case
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
      if(c==-1)spp.topoInit(oneway_gain_cutoff,azwidper); // sets parameters to configuration for topo processing
	else spp.topoInit(100.0,azwidper); // sets parameters to configuration for topo processing
	// ignoring config file
	
	RasList* raslist;
	string ras_file = cfg.str("ras_filename");
	raslist= new RasList(ras_file,"r");

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
	string bidr_prefix = cfg.str("BIDR_FILENAME_PREFIX")+"_"+toStr(trial_height);
	if(c==-1) bidr_prefix+="_init";
      
	BIDR bidr(bidr_prefix,cfg,l1b,spp,"w",TOPO);
 
	// set ranges of pixels to use for topo after first run through
	if(c!=-1){
	  bidr.setOverlap(overlap_col,overlap_width);
	}

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

	  if(dbg.level > 0) 
	    dbg.file << "Processing " << l1b.sab_counter << " ..." << endl;

	  // zero L1I arrays
	  l1i.zeroArrays();
	  if(dbg.level > 0) 
	    dbg.file << "     L1I arrays zeroed" << endl;

	  l1b.copyTo(l1i); // reads L1B record and copies to L1I

	  if(dbg.level > 0) 
	    dbg.file << "     Copied Fields from L1B" << endl;
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
	if(dbg.level > 0) 
	  dbg.file << "     Computed Burst Parameters" << endl;
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
	if(dbg.level > 0) 
	  dbg.file << "     Range Compressed" << endl;

	l1i.azimuthCompress(spp,afs);   
	if(dbg.level > 0) 
	  dbg.file << "     Azimuth Compressed" << endl;
	
	// Initialize standard info output line
	// L1I::calibrate and BIDR::stripProcessAndWrite populate the line
	spp.cout_line="";


        // Used new calibration routine that handles ambiguities
        // Sept 26 2012
        bool init_mode=(c==-1);
        l1i.sartopo_calibrate(spp,afs,trial_height,init_mode);
	//l1i.calibrate(spp,afs,trial_height);
 	if(dbg.level > 0) 
	  dbg.file << "     Calibrated" << endl;
	// estimate noise subtraction parameters if noise subtraction enabled
	if(spp.noise_subtraction_on)
	  spp.estimateNoiseSubtractionParameters(l1i,raslist);
	if(dbg.level > 0) 
	  dbg.file << "     Noise subtracted" << endl;
      
	afs << "End Burst " << burst_no << endl << endl;

      
	l1i.initLatLonTransformation(spp,bidr.getProjection());
	if(dbg.level > 0) 
	  dbg.file << "     LatLon Transformation Initialized" << endl;
        // sets up a simple conversion
	// from (fdop,range) to projected (OblCyl) (lat,lon) coordinates
        // that only applies to current burst

	if(write_l1i) l1i.writeRecord(); // output to L1I file
      
	
	bidr.topoProcessAndWrite(l1i,spp,single_burst_mode,trial_height); 
	if(dbg.level > 0) 
	  dbg.file << "     BIDR::topoProcessAndWrite complete" << endl;

	// output standard info line
	cout << spp.cout_line << endl;
      }
      
      // Output last portions of BIDR images.
      bidr.flush(spp);  

      // Update the PDS label in each output product
      bidr.rewriteHeaders();

      // output burst skip report to annotation file
      spp.reportSkippedBursts(afs);


      if(raslist) delete raslist;

      // Copy Topo arrays (Taking Running Average)
      if(c==-1){
        
	for(int i=0; i<4;i++){
	  int last_good_j=-1;
          int last_good_col=-1;
	  for(int j=0;j<MAX_BIDR_LINES;j++){
	    if(bidr.overlap_col[i][j]){
	      // fill in any gaps
	      if(last_good_col<0) last_good_col=bidr.overlap_col[i][j];
	      float slope=(float)(bidr.overlap_col[i][j]-last_good_col)/
		(float)(j-last_good_j);
	      for(int j2=last_good_j+1;j2<j;j2++){
		overlap_col[i][j2]=(int)(last_good_col+slope*(j2-last_good_j)
					+0.5);
		overlap_width[i][j2]=-1;
	      }
	      last_good_col=bidr.overlap_col[i][j];
	      last_good_j=j;
	    }
	    overlap_col[i][j]=bidr.overlap_col[i][j];
	    overlap_width[i][j]=bidr.overlap_width[i][j];
	  }
	  for(int j2=last_good_j+1;j2<MAX_BIDR_LINES;j2++){
	    overlap_col[i][j2]=last_good_col;
	    overlap_width[i][j2]=-1;
	  }  
	}
	bidr.smoothOverlap(overlap_col,overlap_width);

	if( !write_array(mfp,overlap_col,sizeof(int),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,overlap_width,sizeof(int),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.topo_inc1,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.topo_inc2,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.topo_lat,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.topo_lon,sizeof(float),2,4,MAX_BIDR_LINES)
	    ){
	  cerr << "Error writing intermediate file " << intermed_file << endl;
	  exit(1);
	}

#define DEBUGSARTOPO
#ifdef DEBUGSARTOPO
        for(int j=0;j<MAX_BIDR_LINES;j++){
          fprintf(stderr,"%d ",j);
	  for(int i=0;i<4;i++){
	    fprintf(stderr,"%d %d %d %d ",overlap_col[i][j],overlap_width[i][j],
		   bidr.overlap_col[i][j],bidr.overlap_width[i][j]);
	  }
	  fprintf(stderr,"\n");
	}
#endif

      } // end c==-1
      else{
 	if( !write_array(mfp,bidr.dsigma0,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.sumsigma0,sizeof(float),2,4,MAX_BIDR_LINES) ||	    
	    !write_array(mfp,bidr.dsigma0std,sizeof(float),2,4,MAX_BIDR_LINES) ||	    
	    !write_array(mfp,bidr.dsigma0bias,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.topo_wnl,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.mid_ds0,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !write_array(mfp,bidr.slope_ds0,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.rms_ds0,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.power1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.power2,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.noisefloor1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.noisefloor2,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.ambrat1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !write_array(mfp,bidr.ambrat2,sizeof(float),2,4,MAX_BIDR_LINES) 
	    ){
	  cerr << "Error writing intermediate file " << intermed_file << endl;
	  exit(1);
	}
	for(int i=0; i<4;i++){
	  for(int j=0;j<MAX_BIDR_LINES;j++){
	    if(j==900 && i==0 ){
	      int debugpoint=0;
	      debugpoint++;
	    }
	    dsigma0[i][j][c]=0;
	    rms_ds0[i][j][c]=0;
	    noisefloor1[i][j][c]=0;
	    noisefloor2[i][j][c]=0;
	    power1[i][j][c]=0;
	    power2[i][j][c]=0;
	    ambrat1[i][j][c]=0;
	    ambrat2[i][j][c]=0;

	    dsigma0std[i][j][c]=0;
	    dsigma0bias[i][j][c]=0;
	    sumsigma0[i][j][c]=0;
	    sumsigma0andnoise[i][j][c]=0;
	    wnumlook[i][j][c]=0;
	    int nvalid=0;
	    float corrstd=0;
	    
	    for(int r=-runave_wid;r<=runave_wid;r++){
	      if(j+r>0 && j+r<MAX_BIDR_LINES){
                // BWS Oct 30 2011 exclude BAD_VALUES but not valid negatives
		if(bidr.sumsigma0[i][j+r]>-1e10){
		  dsigma0[i][j][c]+=bidr.dsigma0[i][j+r];
		  sumsigma0[i][j][c]+=bidr.sumsigma0[i][j+r];
                  wnumlook[i][j][c]+=bidr.topo_wnl[i][j+r];
		  nvalid++;
		  dsigma0bias[i][j][c]+=bidr.dsigma0bias[i][j+r];
		  rms_ds0[i][j][c]+=bidr.rms_ds0[i][j+r]*bidr.rms_ds0[i][j+r];
		  power1[i][j][c]+=bidr.power1[i][j+r];
		  power2[i][j][c]+=bidr.power2[i][j+r];
		  noisefloor1[i][j][c]+=bidr.noisefloor1[i][j+r];
		  noisefloor2[i][j][c]+=bidr.noisefloor2[i][j+r];
		  ambrat1[i][j][c]+=bidr.ambrat1[i][j+r];
		  ambrat2[i][j][c]+=bidr.ambrat2[i][j+r];
		}
	      }
	    } // end for r
	    if(nvalid>0){
	      dsigma0[i][j][c]/=nvalid;
	      sumsigma0[i][j][c]/=nvalid;
	      dsigma0bias[i][j][c]/=nvalid;
	      rms_ds0[i][j][c]=sqrt(rms_ds0[i][j][c]/nvalid);
	      power1[i][j][c]/=nvalid;
	      power2[i][j][c]/=nvalid;
	      noisefloor1[i][j][c]/=nvalid;
	      noisefloor2[i][j][c]/=nvalid;
              sumsigma0andnoise[i][j][c]=sumsigma0[i][j][c]+noise_correction_factor*noisefloor1[i][j][c]+noise_correction_factor*noisefloor2[i][j][c];
              // HACK undo noise subtraction
              dsigma0[i][j][c]=dsigma0[i][j][c]+noise_correction_factor*noisefloor1[i][j][c]-noise_correction_factor*noisefloor2[i][j][c];
              // END HACK
	      ambrat1[i][j][c]/=nvalid;
	      ambrat2[i][j][c]/=nvalid;
	      float s01 = (sumsigma0[i][j][c]+dsigma0[i][j][c])/2;
	      float s02 = (sumsigma0[i][j][c]-dsigma0[i][j][c])/2;
	      float SNR1=s01/noisefloor1[i][j][c];
	      float SNR2=s02/noisefloor2[i][j][c];
	      float var1=s01*s01*(1+2/SNR1+1/SNR1/SNR1)/wnumlook[i][j][c];
	      float var2=s02*s02*(1+2/SNR2+1/SNR2/SNR2)/wnumlook[i][j][c];
	      dsigma0std[i][j][c]=sqrt(var1+var2);

	    }
	  } // end for j
	} // end for i
      } // end c!=-1
    } // end trial height loop (for c)
    } // end !read_intermed loop
    else{
      float** ds0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** sums0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** ds0std=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** ds0bias=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** unused=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** wnl=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** rmsds0=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** p1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** p2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** nf1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** nf2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** ar1=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      float** ar2=(float**)make_array(sizeof(float),2,4,MAX_BIDR_LINES);
      for(int c=-1;c<num_heights;c++){
       // Copy Topo arrays (Taking Running Average)
      if(c==-1){
	if( !read_array(mfp,overlap_col,sizeof(int),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,overlap_width,sizeof(int),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,inc1,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,inc2,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,lat,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,lon,sizeof(float),2,4,MAX_BIDR_LINES)
	    ){
	  cerr << "Error reading intermediate file " << intermed_file << endl;
	  exit(1);
	}
      }
      else{
 	if( !read_array(mfp,ds0,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,sums0,sizeof(float),2,4,MAX_BIDR_LINES) ||	    
	    !read_array(mfp,ds0std,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,ds0bias,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,wnl,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,unused,sizeof(float),2,4,MAX_BIDR_LINES)||
	    !read_array(mfp,unused,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,rmsds0,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,p1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,p2,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,nf1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,nf2,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,ar1,sizeof(float),2,4,MAX_BIDR_LINES) ||
	    !read_array(mfp,ar2,sizeof(float),2,4,MAX_BIDR_LINES)
	    ){
	  cerr << "Error reading intermediate file " << intermed_file << endl;
	  exit(1);
	}
	for(int i=0; i<4;i++){
	  for(int j=0;j<MAX_BIDR_LINES;j++){
	    if(j==900 && i==0 ){
	      int debugpoint=0;
	      debugpoint++;
	    }
	    dsigma0[i][j][c]=0;
	    dsigma0std[i][j][c]=0;
	    dsigma0bias[i][j][c]=0;
	    sumsigma0[i][j][c]=0;
	    sumsigma0andnoise[i][j][c]=0;
            wnumlook[i][j][c]=0;
	    rms_ds0[i][j][c]=0;
	    noisefloor1[i][j][c]=0;
	    noisefloor2[i][j][c]=0;
	    power1[i][j][c]=0;
	    power2[i][j][c]=0;
	    ambrat1[i][j][c]=0;
	    ambrat2[i][j][c]=0;


	    int nvalid=0;
	    float corrstd=0;
	    
	    for(int r=-runave_wid;r<=runave_wid;r++){
	      if(j+r>0 && j+r<MAX_BIDR_LINES){
                // BWS Oct 30 2011 exclude BAD_VALUES but not valid negatives
		if(sums0[i][j+r]>-1e10){
		  dsigma0[i][j][c]+=ds0[i][j+r];
		  sumsigma0[i][j][c]+=sums0[i][j+r];
		  wnumlook[i][j][c]+=wnl[i][j+r];
		nvalid++;
		dsigma0bias[i][j][c]+=ds0bias[i][j+r];
                rms_ds0[i][j][c]+=rmsds0[i][j+r];
		power1[i][j][c]+=p1[i][j+r];
		power2[i][j][c]+=p2[i][j+r];
		noisefloor1[i][j][c]+=nf1[i][j+r];
		  noisefloor2[i][j][c]+=nf2[i][j+r];
		  ambrat1[i][j][c]+=ar1[i][j+r];
		  ambrat2[i][j][c]+=ar2[i][j+r];
		}
	      }
	    } // end for r
	    if(nvalid>0){
	      dsigma0[i][j][c]/=nvalid;
	      sumsigma0[i][j][c]/=nvalid;
	      dsigma0std[i][j][c]/=nvalid*sqrt((float)nvalid);
	      dsigma0bias[i][j][c]/=nvalid;
	      rms_ds0[i][j][c]=sqrt(rms_ds0[i][j][c]/nvalid);
	      power1[i][j][c]/=nvalid;
	      power2[i][j][c]/=nvalid;
	      noisefloor1[i][j][c]/=nvalid;
	      noisefloor2[i][j][c]/=nvalid;
              sumsigma0andnoise[i][j][c]=sumsigma0[i][j][c]+noise_correction_factor*noisefloor1[i][j][c]+noise_correction_factor*noisefloor2[i][j][c];
              // HACK undo noise subtraction
              dsigma0[i][j][c]=dsigma0[i][j][c]+noise_correction_factor*noisefloor1[i][j][c]-noise_correction_factor*noisefloor2[i][j][c];
              // END HACK
	      ambrat1[i][j][c]/=nvalid;
	      ambrat2[i][j][c]/=nvalid;
	      float s01 = (sumsigma0[i][j][c]+dsigma0[i][j][c])/2;
	      float s02 = (sumsigma0[i][j][c]-dsigma0[i][j][c])/2;
	      float SNR1=s01/noisefloor1[i][j][c];
	      float SNR2=s02/noisefloor2[i][j][c];
	      float var1=s01*s01*(1+2/SNR1+1/SNR1/SNR1)/wnumlook[i][j][c];
	      float var2=s02*s02*(1+2/SNR2+1/SNR2/SNR2)/wnumlook[i][j][c];
	      dsigma0std[i][j][c]=sqrt(var1+var2);
	    }
	  } // end for j
	} // end for i
      } // end c!=-1
    } // end trial height loop (for c)
      free_array(ds0,2,4,MAX_BIDR_LINES);
      free_array(sums0,2,4,MAX_BIDR_LINES);
      free_array(ds0std,2,4,MAX_BIDR_LINES);
      free_array(ds0bias,2,4,MAX_BIDR_LINES);
      free_array(unused,2,4,MAX_BIDR_LINES);
      free_array(wnl,2,4,MAX_BIDR_LINES);
      free_array(rmsds0,2,4,MAX_BIDR_LINES);
      free_array(p1,2,4,MAX_BIDR_LINES);
      free_array(p2,2,4,MAX_BIDR_LINES);
      free_array(nf1,2,4,MAX_BIDR_LINES);
      free_array(nf2,2,4,MAX_BIDR_LINES);
      free_array(ar1,2,4,MAX_BIDR_LINES);
      free_array(ar2,2,4,MAX_BIDR_LINES);
    } // end if read_intermed
    fclose(mfp);

    // do nothing if num_heights<0 (a debugging case)
    if(num_heights>0){
    for(int i=0; i<5;i++){
	for(int j=0;j<MAX_BIDR_LINES;j++){
	  // Fit polynominals of dsigma0 to height

         if(j==900 && i==0 ){
	    int debugpoint=0;
	    debugpoint++;
	  }
        

	 // Originally a line was valid only if it had sumsigma0!=0
         // for all heights
         // This threw out lines where certain heights eliminate overlap
         // For the original SARTopo processing this never happened
         // but the latest version eliminates region with high ambiguity
         // and/or very low gain. To accomodate this I have modified the
         // code to allow processing when so heights are bad so long as
         // a zero crossing in d=dsigma0/sumsigma0 exists and d is montonic
         // BWS Oct 8 2012
	  bool valid_line=true;
          int num_valid_heights=0;
          double d=0,dold=0;
          bool d_increasing=false;
	  nzerocross[i][j]=0;
          int zero_cross_c[20];	  
          monotonic[i][j]=true;
          bool zero_double_cross=false;
          if(i<4){
	    for(int c=0;c<num_heights;c++){
	      dold=d;

                // Using sumsigma0andnoise instead od sumsigma0 to resolve negative sumsigma0 errors
                // should not change crossover point for positive sumsigma0 cases
                // March 7 2013 - BWS
	      if(sumsigma0andnoise[i][j][c]>0){
		num_valid_heights++;

 
		d=dsigma0[i][j][c]/(sumsigma0andnoise[i][j][c]);

                // extrapolate to large negative heights
                if(num_valid_heights==1){
		  for(int cc=c-1;cc>=0;cc--){ 
		    dsigma0[i][j][cc]=dsigma0[i][j][c];
		    sumsigma0andnoise[i][j][cc]=sumsigma0andnoise[i][j][c];
		  }

		}
		if(num_valid_heights>1){
		  // check for zero_crossing
		  if(d*dold<0){
		    zero_cross_c[nzerocross[i][j]]=c;
		    if(fabs(dold)<fabs(d))zero_cross_c[nzerocross[i][j]]=c-1;
		    nzerocross[i][j]++;
		  }
		  
                  // check for monotonicity
		  if(num_valid_heights>2){
		    if(d_increasing!=(d>=dold)){
		      monotonic[i][j]=false;
		      break;
		    }
		  }
		  d_increasing=(d>=dold);		    		  
		} // end if first valid height found
	      } // end if sumsigma0 > 0
	      else if(num_valid_heights>0){
                // extrapolate dsigma0 and sumsigma0 to large positive heights
		for(int cc=c;cc<num_heights;cc++){ 
		  dsigma0[i][j][cc]=dsigma0[i][j][c-1];
		  sumsigma0andnoise[i][j][cc]=sumsigma0andnoise[i][j][c-1];
		}
		break;
	      }
            } // end height loop
	    if(num_valid_heights<5 || nzerocross[i][j]==0 ){
	      valid_line=false;
	    }
	  } // end if i<4 case

          // i==4 case needs both 2/3 and 3/4 overlap
	  else{
	    for(int c=0;c<num_heights;c++){
	      dold=d;
	      if(sumsigma0andnoise[1][j][c]>0 && sumsigma0andnoise[2][j][c]>0){
		num_valid_heights++;


 
		d=dsigma0[1][j][c]/(sumsigma0andnoise[1][j][c]);
		d+=dsigma0[2][j][c]/(sumsigma0andnoise[2][j][c]);

                // extrapolate to large negative heights
                // should not be needed already done for i==1 and i==2
                //if(num_valid_heights==1){
		//  for(int cc=c-1;cc>=0;cc--){ 
		//    dsigma0[1][j][cc]=dsigma0[1][j][c];
		//    dsigma0[2][j][cc]=dsigma0[2][j][c];
		//    sumsigma0andnoise[1][j][cc]=sumsigma0andnoise[1][j][c];
		//    sumsigma0andnoise[2][j][cc]=sumsigma0andnoise[2][j][c];
		//  }
		//}
		if(num_valid_heights>1){
		  // check for zero_crossing
		  if(d*dold<0){
		    zero_cross_c[nzerocross[i][j]]=c;
                    if(fabs(dold)<fabs(d)) zero_cross_c[nzerocross[i][j]]=c-1;
		    nzerocross[i][j]++;
		  }
		  
                  // check for monotonicity
		  if(num_valid_heights>2){
		    if(d_increasing!=(d>=dold)){
		      monotonic[i][j]=false;
		    }
		  }
		  d_increasing=(d>=dold);		    		  
		} // end if num_valid_heights>1
	      } // end if sumsigma0s > 0
	      else if(num_valid_heights>0){
                // extrapolate dsigma0 and sumsigma0 to large positive heights
                // should not be needed already done for i=1 and i=2
		//for(int cc=c;cc<num_heights;cc++){ 
		//  dsigma0[1][j][cc]=dsigma0[1][j][c-1];
		//  dsigma0[2][j][cc]=dsigma0[2][j][c-1];
		//  sumsigma0andnoise[1][j][cc]=sumsigma0andnoise[1][j][c-1];
		//  sumsigma0andnoise[2][j][cc]=sumsigma0andnoise[2][j][c-1];
		//}
                
		break;
	      }
	    } // end height loop
	    if(num_valid_heights<5 || nzerocross[i][j]==0){
	      valid_line=false;
	    }
	  } // end i==4 case (beam24)
	    
          if(valid_line){
	    // Estimate height and std of height
            Dvec x("",height_fit_order+2);
            Dvec y("",height_fit_order+2);
            Dvec p("",height_fit_order+1); // BWS Oct 25 2012 fixed BUG
	                                   // p was too large, so a quadratic
	                                   // fit was performed...

	    cmin2[i][j]=0;
            int d=0;
            float minval2=100000000000;
            for(int c=0;c<num_heights;c++){
	      float val,val2;
               // modifying val to avoid negative sumsigma0 problems
                // should not change crossover point for posiitve sumsigma0 cases
                // March 7 2013 - BWS
	      if(i<4){
		val=fabs(dsigma0[i][j][c]/(sumsigma0andnoise[i][j][c]));
		val2=rms_ds0[i][j][c];
	      }
	      else{
		val=fabs(dsigma0[1][j][c]/(sumsigma0andnoise[1][j][c])
			 +dsigma0[2][j][c]/(sumsigma0andnoise[2][j][c]));
		val2=rms_ds0[1][j][c]+rms_ds0[1][j][c];
	      }
	      if(val2<minval2 ){
		minval2=val2;
		cmin2[i][j]=c;
	      }
	    }
	    cmin[i][j]=0;
            int minval=100000000;

            // find zero_crossing point of ds0/sums0 closest to rms_ds0 minimum
	    for(int c=0;c<nzerocross[i][j];c++){
	      int val=abs(zero_cross_c[c]-cmin2[i][j]);
	      if(val<minval){
		cmin[i][j]=zero_cross_c[c];
		minval=val;
	      }
	    }

	    d=cmin[i][j]-(height_fit_order+1)/2;
            if(d<0)d=0;
	    if(d>num_heights-height_fit_order-2){
	      d=num_heights-height_fit_order-2;
	    }
            for(int c=d;c<d+height_fit_order+2;c++){
              if(i<4){
                // modifying d to avoid negative sumsigma0 problems
                // should not change crossover point for posiitve sumsigma0 cases
                // March 7 2013 - BWS
		x(c-d)=dsigma0[i][j][c]/(sumsigma0andnoise[i][j][c]);
	      }
	      else{
		x(c-d)=dsigma0[1][j][c]/(sumsigma0andnoise[1][j][c]);
		x(c-d)+=dsigma0[2][j][c]/(sumsigma0andnoise[2][j][c]);
	      }
	      y(c-d)=min_trial_height+c*height_step;
	    }
	  
	    try{
	      polyfit(p,y,x);
	      height[i][j]=p(0);
	      float dhneg=0;
	      float dhpos=0;
              float dhmin=0;
	      float dhbneg=0;
	      float dhbpos=0;
              // store coefficients into heightpoly
              for(int c=0;c<height_fit_order+1;c++){
		heightpoly[i][j][c]=p(c);
	      }

	      if(i<4){
		float dx=dsigma0std[i][j][cmin[i][j]]; // standard deviation in sigma0 difference between beams
		float dxmin=dsigma0[i][j][cmin[i][j]]/sumsigma0andnoise[i][j][cmin[i][j]]; // diff/sum criteria at minimum value
		                                                                   // best candidate height

                // dxbneg and dxbpos values used to compute heighbias quantity
		float dxbneg=dxmin-(dsigma0[i][j][cmin[i][j]]-dsigma0bias[i][j][cmin[i][j]])/(sumsigma0andnoise[i][j][cmin[i][j]]+dsigma0[i][j][cmin[i][j]]);
		float dxbpos=(dsigma0[i][j][cmin[i][j]]+dsigma0bias[i][j][cmin[i][j]])/(sumsigma0andnoise[i][j][cmin[i][j]]-dsigma0[i][j][cmin[i][j]])-dxmin;
                
		for(int k=height_fit_order;k>0;k--){
		  dhneg+=p(k);
		  dhneg*=-dx;
		  dhpos+=p(k);
		  dhpos*=dx;

		  dhbneg+=p(k);
		  dhbneg*=-dxbneg;
		  dhbpos+=p(k);
		  dhbpos*=dxbpos;
		  
		  dhmin+=p(k);
		  dhmin*=dxmin;
		}
		dhmin+=p(0);
		dheight[i][j]=(fabs(dhpos)+fabs(dhneg))/2.0;
		heightbias[i][j]=(fabs(dhbpos)+fabs(dhbneg))/2.0;
		dheight[i][j]+=fabs(dhmin-min_trial_height-cmin[i][j]*height_step);

                // compute derivative of height with respect to multiplicative error in noise floor.
		float dfnf=0.1; // small multiplicative delta on noisefloor
		float ds0fnf=dsigma0[i][j][cmin[i][j]]+noisefloor1[i][j][cmin[i][j]]-noisefloor2[i][j][cmin[i][j]]-(1+dfnf)*noisefloor1[i][j][cmin[i][j]]+(1+dfnf)*noisefloor2[i][j][cmin[i][j]];
		float sums0fnf=sumsigma0[i][j][cmin[i][j]]+noisefloor1[i][j][cmin[i][j]]+noisefloor2[i][j][cmin[i][j]]-(1+dfnf)*noisefloor1[i][j][cmin[i][j]]-(1+dfnf)*noisefloor2[i][j][cmin[i][j]];

		float xfnf=(ds0fnf/sums0fnf);
                float hfnf=0;                
		for(int k=height_fit_order;k>0;k--){
		  hfnf+=p(k);
		  hfnf*=xfnf;
		}
		hfnf+=p(0);
                dhdfnf[i][j]=(hfnf-min_trial_height-cmin[i][j]*height_step)/dfnf;
	      }


              // special case for beam2/4 (i=4)
	      else{
		dhmin+=p(0);
		dheight[i][j]=sqrt((dheight[1][j]*dheight[1][j]
				    +dheight[2][j]*dheight[2][j])/2);
		heightbias[i][j]=0; // 1/2 and 2/3 biases should cancel out

                // compute derivative of height with respect to multiplicative error in noise floor.
                dhdfnf[i][j]=dhdfnf[1][j]-dhdfnf[2][j];
	      }
	    }
            catch(ErrorMessage e){
	      cerr << "Warning Bad Height to dsigma0 fit for line " << j 
		   << " overlap beam " << i+1 <<" and "<<i+2 << endl;
	      height[i][j]=-1000;
	      dheight[i][j]=0;
              heightbias[i][j]=0;	      
	    }
	  } // end if valid_line
          else{
	    height[i][j]=-1000;
	    dheight[i][j]=0;
	    heightbias[i][j]=0;	      
	  }
	} // end for j
      } // end for i

      if(dbg.level > 0) 
	dbg.file << "Height computed ..." << endl;
	
      // Output arrays
      float dinc=0;
      for(int j=0;j<MAX_BIDR_LINES;j++){
	if(inc1[0][j]>0 && inc1[3][j]>0){
	  dinc+= inc1[3][j]-inc1[0][j];
	}
      }
      int bad_beam=3;
      if(dinc<0) bad_beam=0;
      for(int i=0;i<5;i++){ 
	for(int c=0;c<MAX_BIDR_LINES;c++){
	  if(height[i][c]>-999){
            int mask=0;
            if(dheight[i][c]>dhthresh) mask+=16;
            int outcol, outwid;
            float outinc, outlat,outlon,ar;
            if(i<4){
	      if(inc1[i][c]*radtodeg<10.0) mask+=1;
	      if(lat[i][c]<-100000) mask+=2;
	      if(i==1){  
		if(overlap_width[1][c]>4.0*overlap_width[2][c]) mask+=4; 
	      }
	      if(i==2){  
		if(overlap_width[2][c]>4.0*overlap_width[1][c]) mask+=4; 
	      }
	      if(i==bad_beam) mask+=8;
	      if(heightbias[i][c]>dhthresh && (i==1 || i==2)) mask+=32;


 	      outcol=overlap_col[i][c];
              outwid=overlap_width[i][c];
              outinc=radtodeg*(inc1[i][c]+inc2[i][c])/2;
              outlat=lat[i][c]*radtodeg;
              outlon=360-lon[i][c]*radtodeg;
	      ar=MAX(ambrat1[i][c][cmin[i][c]],ambrat2[i][c][cmin[i][c]]);
	    }
            // beam 2/4 special case
            else{
	      if(inc1[1][c]*radtodeg<10.0 || inc1[2][c]*radtodeg<10.0)
		mask+=1;
	      if(lat[1][c]<-100000 || lat[2][c]<-100000) mask+=2;
	      outcol=(overlap_col[1][c]+overlap_col[2][c])/2;
              outwid=abs(overlap_col[1][c]-overlap_col[2][c]);
              outinc=radtodeg*(inc1[1][c]+inc2[2][c])/2;
              outlat=radtodeg*(lat[1][c]+lat[2][c])/2;
              float lon1=lon[1][c];
              float lon2=lon[2][c];
              if(lon1<lon2-pi) lon1+=2*pi;
              if(lon1>lon2+pi) lon1-=2*pi;
              outlon=(lon1+lon2)/2;
              if(outlon<0) outlon+=2*pi;
              if(outlon>2*pi) outlon-=2*pi;
              outlon=360-radtodeg*outlon;

	      float ar1=MAX(ambrat1[1][c][cmin[i][c]],ambrat2[1][c][cmin[i][c]]);
	      float ar2=MAX(ambrat1[2][c][cmin[i][c]],ambrat2[2][c][cmin[i][c]]);             ar=MAX(ar1,ar2);
	    }

	      if(!monotonic[i][c]) mask+=64;
              if(nzerocross[i][c]>1) mask+=128;
              if(nzerocross[i][c]==0) mask+=256;
              if(abs(cmin[i][c]-cmin2[i][c])>1) mask+=512;

              // flag ambiguity ratio > 20%
              
	      if(ar>0.2)mask+=1024;
	      if(fabs(dhdfnf[i][c])>10) mask+=2048; //noisefloor needs to be known too precisely 100% change -> 10 km height change.
 
	    fprintf(hfp,"%d %g %g %d %d %d %g\n",c+1,height[i][c],
		    dheight[i][c],outcol,mask,outwid,
		    heightbias[i][c]);
            if(mask==0 || OUTPUT_ALL_DATA){
              float along_width=(runave_wid*2+1)*0.175;
	      fprintf(csvfp[i],"%7.3f,%6.3f,%5.2f,%5.2f,%5.2f,%d,%d,%d,%d,%d %7.2f\n",
		      outlon,outlat,outinc,outwid*0.175,
		      along_width,int(floor(1000*height[i][c]+0.5)),
		      int(floor(1000*dheight[i][c]+0.5)),
		      mask,c,outcol,1000*dhdfnf[i][c]);
		      
	    } // end loop to output only mask ==0 data unless OUTPUT_ALL_DATA is set
	  } //end if statement height > -999 (km) sanity check to skip invalid lines
     
	} // end loop over BIDR lines
	fprintf(hfp,"&\n");
      } // end loop over profiles 0-4
      fclose(hfp);
      for(int d=0;d<5;d++)fclose(csvfp[d]);

      if(dbg.level > 0) 
	dbg.file << "Height file written ..." << endl; 
      int nlines=MAX_BIDR_LINES;
      if(fwrite(&nlines,sizeof(int),1,afp)!=1 ||
	 fwrite(&num_heights,sizeof(int),1,afp)!=1 ||
	 fwrite(&height_fit_order,sizeof(int),1,afp)!=1){
	cerr << "Error writing header to Topo intermediate data file " << auxfile 
	     << endl;
	exit(1);
      }

      if(dbg.level > 0) 
	dbg.file << "Auxiliary file header written ..." << endl; 
      
      if(!write_array(afp,dsigma0,sizeof(float),3,4,MAX_BIDR_LINES,num_heights) ||
	 !write_array(afp,dsigma0std,sizeof(float),3,4,MAX_BIDR_LINES,num_heights) ||
	 !write_array(afp,sumsigma0,sizeof(float),3,4,MAX_BIDR_LINES,num_heights) ||
	 !write_array(afp,heightpoly,sizeof(float),3,5,MAX_BIDR_LINES,height_fit_order+1) ||
	 !write_array(afp,overlap_col,sizeof(int),2,4,MAX_BIDR_LINES) ||
	 !write_array(afp,overlap_width,sizeof(int),2,4,MAX_BIDR_LINES) || 	 
	 !write_array(afp,height,sizeof(float),2,5,MAX_BIDR_LINES) ||
	 !write_array(afp,dheight,sizeof(float),2,5,MAX_BIDR_LINES) ||
	 !write_array(afp,heightbias,sizeof(float),2,5,MAX_BIDR_LINES)){
	cerr << "Error writing arrays to Topo auxiliary data file " << auxfile 
	     << endl;
	exit(1);
      }

      fclose(afp);

      if(dbg.level > 0) 
	dbg.file << "Auxiliary file data written ..." << endl; 
    } // end num_heights > 0 case
      free_array(dsigma0,3,4,MAX_BIDR_LINES,num_heights);
      free_array(dsigma0std,3,4,MAX_BIDR_LINES,num_heights);
      free_array(dsigma0bias,3,4,MAX_BIDR_LINES,num_heights);
      free_array(rms_ds0,3,4,MAX_BIDR_LINES,num_heights);
      free_array(power1,3,4,MAX_BIDR_LINES,num_heights);
      free_array(power2,3,4,MAX_BIDR_LINES,num_heights);
      free_array(noisefloor1,3,4,MAX_BIDR_LINES,num_heights);
      free_array(noisefloor2,3,4,MAX_BIDR_LINES,num_heights);
      free_array(ambrat1,3,4,MAX_BIDR_LINES,num_heights);
      free_array(ambrat2,3,4,MAX_BIDR_LINES,num_heights);
      free_array(sumsigma0,3,4,MAX_BIDR_LINES,num_heights);
      free_array(sumsigma0andnoise,3,4,MAX_BIDR_LINES,num_heights);
      free_array(wnumlook,3,4,MAX_BIDR_LINES,num_heights);
      free_array(overlap_col,2,4,MAX_BIDR_LINES);
      free_array(overlap_width,2,4,MAX_BIDR_LINES);
      free_array(inc1,2,4,MAX_BIDR_LINES);
      free_array(inc2,2,4,MAX_BIDR_LINES);
      free_array(lat,2,4,MAX_BIDR_LINES);
      free_array(lon,2,4,MAX_BIDR_LINES);
      free_array(heightpoly,3,5,MAX_BIDR_LINES,height_fit_order+1);
      free_array(height,2,5,MAX_BIDR_LINES);
      free_array(dheight,2,5,MAX_BIDR_LINES);
      free_array(dhdfnf,2,5,MAX_BIDR_LINES);
      free_array(heightbias,2,5,MAX_BIDR_LINES);
      free_array(monotonic,2,5,MAX_BIDR_LINES);
      free_array(nzerocross,2,5,MAX_BIDR_LINES);
      free_array(cmin,2,5,MAX_BIDR_LINES);
      free_array(cmin2,2,5,MAX_BIDR_LINES);



      if(dbg.level > 0) 
	dbg.file << "Arrays freed ..." << endl;       
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



