#include <string.h>
#include <strings.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <getopt.h>
#include "Plot.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"
#include "BurstData.h"
#include "L1B.h"
#include "config_keywords.h"

//! Define command line option flags
#define OPTSTRING "hxs:c:o:O:r:t:P:p:T:t:"
extern int optind;
using std::cout;
using std::cerr;
using std::endl;
using std::set_unexpected;
using std::terminate;


void myUnexpected() throw()
{
  cout << "Unexpected exception" << std::endl;
  terminate();
}

int main(int argc, char* argv[])
  
{
  typedef std::istringstream ISTRINGSTREAM;
  
  std::set_unexpected(myUnexpected);
  
  try{
    
    //------------------------
    // Parse the command line
    //------------------------
    
    char* cfg_filename = NULL;
    char* output_prefix = NULL;
    char* output_path = NULL;
    
    double rangeLow;
    double rangeHi;
    bool rangeFlag = false;
    
    bool xtarget = false;
    char* typeString = NULL;
    char* titleString = " ";

    int subSample = 1;
    unsigned int nPlot=1,iPlot; // default 1 plot, iplot is index
    unsigned int vPlot[100]; // vector of plots per page
    unsigned int iAbscissa=0;
    bool multipleAbscissa = false;
	
    while (1)
      {
	int c = getopt(argc,argv,OPTSTRING);
	if (c == 'h')
	  {
	    cout << "Option -h prints this message" << endl;
	    cout << "Option -c <filename> identifies the config file to read from" << endl;
	    cout << "Option -t <L1A,L1P> identifies the data type" << endl;
	    cout << "Option -o <prefix> identifies the output prefix for plot files" << endl;
	    cout << "Option -O <path> identifies the output path for plot files" << endl;
	    cout << "Option -s <subsambple factor>" << endl;
	    cout <<
	      "Option -r <low:high> identifies the optional parameter range for the 1st param" << endl;
	    cout << "Option -x means send plots to the current X-terminal" << endl;
	    cout << "Option <-P n1:n2:n3...nm> Multiple plots, 1 abscissa for m plots with n1,n2,..nm ordinates on each page. Paramlist is: x y1_n1 y2_n1...yn1_n1..... y1_nm..ynm_nm " << endl;
	    cout << endl;
	    cout << "Option -p <n1:n2:n3...nm> Multiple plots with multiple abscissae and ordinates, with m plots with n1,n2,..nm ordinates on each page. paramlist is: x1 y1_n1 y2_n1...yn1_n1.... x_nm y1_nm..ynm_nm " << endl;
	    cout << endl;
	    cout << "Option -T <Appended Title String>" << endl;
	    cout << endl;
	    cout << endl;
	    cout << "paramlist is a list of field names including:\n" << endl;
	    cout << "be1tmp be2tmp be3tmp be4tmp be5tmp tadcal1 tadcal2 tadcal3 tadcal4 wgb1t1 wgb3t1 wgb3t2 wgb3t3 wgb5t1 fwdtmp nsdtmp dcgttm cucttm fguttm usotmp cputmp memtmp mratmp mruttm scwg_tmp feed_tmp hga_tmp adctmp sadctmp epctmp ep1ttm p_stmp p_sttm esstmp pcutmp twttmp tw1ttm rlotmp lnatmp diptmp evdtmp ecltmp tadcal5 pcu5v_72 pcu5v_74 pcu15v_76 pcu15v_78 pcu12v_80 pcu9v_86 pcu9v_88 pcu5i_73 pcuti_75 pcu15i_77 pcu15i_79 pcu9i_87 pcu12i_81 pcu9i_89 frwdpw nsdcur pcucur catcur hpapsm p_smon cpbnkv essvlt svlsta tadcal5 tadcl7rc_bw adc fast_csf pllmon slow_cfs scpr rel_time sclk header_tfi slow_tfi fast_tfi csd chirp_length bpd pri rwd hip cip rip ctrx tro brst tn baq_mode rad fin sin cnt_rl at1_db at3_db at4_db sync sclk record_id header_tnc header_typ header_tca header_tcb header_end pwri vicc vimc tail_len tail_id sab_counter sab_len fswm fswc ctbc ctps " << endl;
	  }
	else if (c == 'c')
	  cfg_filename = optarg;
	else if (c == 't')
	  typeString=optarg;
	else if (c == 'T')
	  titleString=optarg;
	else if (c == 'o')
	  output_prefix = optarg;
	else if (c == 'O')
	  output_path = optarg;
	else if (c == 'x')
	  xtarget = true;
	else if (c == 'r')
	  {
	    rangeFlag = true;
	    char c1 = ' ';
	    ISTRINGSTREAM is(optarg);
	    is >> rangeLow;
	    is >> c1;
	    is >> rangeHi;
	    if (c1 != ':' )
	      {
		cout << "Error determining parameter range: " << optarg << endl;
		return(-1);
	      }
	  }
	else if (c == 's')
	  {
	    ISTRINGSTREAM is(optarg);
	    is >> subSample;
	  }
	else if ( c == 'P' || c == 'p' )
	  {
	    if ( c == 'p' ) {
	      multipleAbscissa= true;
	    }
	    char c2 = ' ';
	    ISTRINGSTREAM is(optarg);
	    is >> nPlot;
	    for(iPlot=0;iPlot<nPlot;iPlot++){
	      is >> c2;
	      if ( c2 != ':' ) {
		cout << "Error parsing Page Option" << optarg << endl;
	      }
	      is >> vPlot[iPlot]; //read in plots per page
	    }
	  }
	else if (c == -1) break;
      }
    
    // check to make sure all necessary command line info is available
    if (cfg_filename == NULL || typeString == NULL )
      {
	cerr << "Usage: monitor -h -x -c cfgfile -r lo:hi -t data_type paramlist";
	cerr << endl;
	return(-1);
      }
    
    //------------------------------
    // Load configuration parameters
    //------------------------------
    
    string o_prefix_str;
    if (output_prefix) o_prefix_str = output_prefix;
    string o_path_str;
    if (output_path) o_path_str = output_path;

    Config cfg(cfg_filename);
    
    Frame::config(cfg);
    BurstData::config(cfg);
    
    // Set up a data object according to data type argument
    string types[]={"L1BP","L1BA"};
    string filetypes[]={"passive","active"};
    string cfg_tag[]={L1B_PASSIVE_MODE_FILENAME, 
		      L1B_ACTIVE_MODE_FILENAME};
    
    //    int nTag =sizeof(cfg_tag)/sizeof(cfg_tag[0]);
    //    int nTypes = sizeof(types)/sizeof(types[0]);
    int nFileTypes = sizeof(filetypes)/sizeof(filetypes[0]);
    


    BurstData* data;

    unsigned int iRow=0;
    unsigned int nRecords;
    bool data_type_found=false;
    for(int c=0;c<nFileTypes;c++){
      if(strcasecmp(typeString,types[c].c_str())==0){
	string filename=cfg.str(cfg_tag[c]);
	if(c>1){ // L1I 
	  //  data=new L1I(filename,"rb",filetypes[c]);
	}
	else{ // L1B
	  data=new L1B(filename,"rb",filetypes[c]);
	}
	data_type_found=true;	
	break;
      }
    }
    if(!data_type_found){
      cerr << "monitor: Invalid data type " << typeString << endl;
      cerr << "valid data types are:\n";
      for (int c=0;c<nFileTypes;c++){
	cerr << types[c] << endl;
      }
      return(-1);
    }
    
    // Read Header
    
    data->readHeader();
    // Create Arrays for Data
    double *prow; // pointer to table row. double is coded and working.
    UnitVar *prow2; // using a unit variable instead of a double (debugging)
    unsigned int nFields;
    nFields=argc-optind;
    prow = new double [nFields];
    prow2= new UnitVar[nFields];
    unsigned int iCol;
    
    cout << "Counting Good Data:";

    /*TEST
    for (iCol=0;iCol < nFields;iCol++){
      cout << *(argv+optind+iCol) << endl;
      cout << isParam(*(argv+optind+iCol)) << endl;
    }
    exit(0);
    *///TEST


    while(! data->eof()){
      if (!(*data).returnParams(prow,argv+optind,1,rangeFlag,rangeLow,
				rangeHi))
	{
	  iRow++;
	}  
    }
    
    char xString[170]="";
    char  xunit[10]="";
    
    char yString[170]="";
    char yunit[10]="";
    
    
    
    nRecords = iRow;
    cout << nRecords << " Records" << endl;
    if ( !nRecords ) { cout << "No Elements to Plot" << endl; exit(0);}
    
    //   new Umat(filetypes[c],nRecords,nFields);
    Umat returnData("returnData",1+(nRecords-1)/subSample,nFields);
    Uvec returnIndex("index",1+(nRecords-1)/subSample);
    
    
    //    returnData = new Umat("returnData",nRecords,nFields);

    data_type_found=false;
    for(int c=0;c<nFileTypes;c++){
      if(strcasecmp(typeString,types[c].c_str())==0){
	string filename=cfg.str(cfg_tag[c]);
	if(c>1){ // L1I 
	  //  data=new L1I(filename,"rb",filetypes[c]);
	}
	else{ // L1B
	  data=new L1B(filename,"rb",filetypes[c]);
	}
	data_type_found=true;	
	break;
      }
    }
    if(!data_type_found){
      cerr << "monitor: Invalid data type " << typeString << endl;
      cerr << "valid data types are:\n";
      for (int c=0;c<nFileTypes;c++){
	cerr << types[c] << endl;
      }
      return(-1);
    }
    
    // Read Header
    
    data->readHeader();
    iRow=0;

    cout << "Retrieving Data" ;
    if ( subSample != 1 ){
      cout << " with sub-sampling factor:" << subSample;
    }
    cout << endl;

    while(! data->eof()){
      //      data->extractParams(argv+optind,argc-optind,rangeFlag,rangeLow,rangeHi);
      if (!(*data).returnParams(prow,argv+optind,nFields,rangeFlag,rangeLow,
				rangeHi))
	{
	  if ( (iRow%subSample) == 0 ) { 
	    for (iCol=0;iCol<nFields;iCol++){
	      returnData(iRow/subSample,iCol)=*(prow+iCol);
	    }
	    returnIndex(iRow/subSample)=iRow/subSample;
	  }
	  iRow++;
	}
      
    }
    
    
    //-----------------------------------------------------------
    // Make plots
    //-----------------------------------------------------------
    /* to do:
       (1) change returnParams from double array to UNIT array
       (2) count number of different units = nUnit
       (3) make Plot class an array of [nUnit]
       (4) plot same unit stuff on same plot.
    */
    


    cout << "Building plots..." << endl;

    string target;
    if (xtarget)
      {
	target = "x";
      }
    else
      {
	target = "eps-landscape";
      }

    // Define Plot Colors For Multiple fields
    string symColor[] ={"black","red","blue","green",
			       "cyan","orange",};
    string symSymbol[]={"circle","up-triangle","square",
			"diamond","down-triangle"};
    int nColor=sizeof(symColor)/sizeof(symColor[0]);
    int nSymbol=sizeof(symSymbol)/sizeof(symSymbol[0]);
    // Define Title Strings
    char TitleString[170]="";
    strcpy(TitleString,cfg_filename);
    strcat(TitleString,":");
    strcat(TitleString,typeString);
    strcat(TitleString," ");
    strcat(TitleString,titleString);
    strcat(TitleString,"#");

    int b = strlen(TitleString);

    int pOff=1;
    if ( nPlot == 1 ) {
      vPlot[0]=nFields-1;
    }
    //Plot plot[nPlot];
    vector<Plot> plot;
    plot.resize(nPlot);
    /* Algorithm is as follows:
       If there is one field then:
       .     plot Field vs. Index*SubSample (loop over plots is moot)

       If there is more than one field :
       .     Loop over plot number (to # of Plots)
       .       Plot X and Add Y_1 to Y_n where

    */
    for(iPlot=0;iPlot<nPlot;iPlot++){
      cout << "Plotting :"<< *(argv+optind+iAbscissa) << endl;

      if ( nFields == 1) {   // plot single field vs. index
	plot[iPlot].addXY(returnIndex*subSample,"Index",returnData.getCol(0),
			  "Y",line("none"),sym("circle","black"));
	plot[iPlot].setYlabel(*(argv+optind));
	plot[iPlot].setXlabel("Index");
      }
      else { // plot multiple fields vs. 1st field
	for (iCol = pOff; iCol<(pOff+vPlot[iPlot]);iCol++){
	  cout << "                   " << *(argv+optind+iCol) << endl;
	  plot[iPlot].addXY(returnData.getCol(0+iAbscissa),"X",
			    returnData.getCol(iCol),"Y",
			    line("none"),
			    sym(//use mod (%) to wrap and run color and symbol 
				symSymbol[((iCol/nColor)%nSymbol
					   +iCol%nColor)%nSymbol  ],
				symColor[iCol%nColor]
				));
	  plot[iPlot].addLegend(*(argv+optind+iCol));
	}
	
	//	cout << "x1 " << xString << endl;
	strcpy(xString,*(argv+optind+iAbscissa));
	//cout << "x2 " << xString << endl;
	strcat(xString," (");
	//cout << "x3 " << xString << endl;
	strcat(xString,xunit);
	//cout << "x4 " << xString << endl;
	strcat(xString,")");
	//cout << "x5 " << xString << endl;
	plot[iPlot].setXlabel(xString);
	
	//cout << "y1 " << yString << endl;
	strcpy(yString,"(");
	//cout << "y1 " << yString << endl;
	strcat(yString,yunit);
	//cout << "y1 " << yString << endl;
	strcat(yString,")");
	//cout << "y1 " << yString << endl;
	plot[iPlot].setYlabel(yString);
	
      }
      TitleString[b]=char(48+(iPlot+1)/10);
      TitleString[b+1]=char(48+(iPlot+1)%10);
      plot[iPlot].setTitle( TitleString );
      plot[iPlot].setFile(o_path_str+o_prefix_str
			  + char(48+(iPlot+1)/10)+char(48+(iPlot+1)%10)
			  +".eps");
      plot[iPlot].show(target);
      pOff = pOff + vPlot[iPlot];
      if ( multipleAbscissa )
	{
	  iAbscissa = pOff++; // correct x & shift y columns by one. 
	}
    }// iPlot loop
    exit(0);
    
}

catch(ErrorMessage& e)
  {
  cerr << "Error: " << e.msg << std::endl;
  }
catch(...)
  {
  cerr << "Exception caught" << std::endl;
  }

return(0);
}







