
// ----------------------------------------------------------
//		IO: CLASSES / routines / definitions 
// ----------------------------------------------------------

#include <stdlib.h>
#include <string.h>
#include <string>
#include <fstream>
#include <iostream>

#include "IOinput.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::terminate;


//----------------------------------------------------
// --------- PEF input file definition
//----------------------------------------------------

/* **************************************
   during initialize:
	- set all parameters to NA 
****************************************** */

   void PEFCUT_FILES::initialize ( int cmddefinition )
	{ 

	int loop ;		
				
		for (loop = 0 ; loop < MAXNUMBEROFRECORDS ; loop++ ) {
         		//strcpy (sclk_time[loop] , "NA" );
			sclk[loop] = 0 ;
			command[loop] = COMMANDERROR ;
			action[loop] = NODATA ;
		} //end for loop

		cmdtype = cmddefinition;
		numberofcommands = 0;
		strcpy (current_sclk, "NA" );
        	strcpy (current_command, "NA" );
		strcpy (current_action, "NA" );
	}


// ------------------------------------------

/* **************************************
   during readin:
	- set all parameters to NA 
****************************************** */

  int PEFCUT_FILES::readin ( const char* filename )
						// readin file...
	{
	char inputstring [ MAXFILESTRINGSIZE ];
	char stringtoparse [ MAXFILESTRINGSIZE ];

	char timestringend[] = ":" ;
	char commandstringend[] = "," ;
	char actionstringend[] = ";" ;
	char cmddelimiter [] = "CMD";

	char timestring [ MAXFILESTRINGSIZE ] ;
	char commandstring [ MAXFILESTRINGSIZE ] ;

	char cmd81vs6 [ 40 ] ;	// either: CMD,6CHG_SC_TM_IMM or CMD,81

	char tempinputstring [ MAXFILESTRINGSIZE ] ;
	char actionstring [ MAXFILESTRINGSIZE ] ;
	int cmdlinenumber = 0;


	ifstream infile ( filename );

 	if ( infile.good () ) {
	   // should get header and place it in a report
	   // skip for now 2 lines
           	infile.getline ( inputstring , MAXFILESTRINGSIZE ) ;
           	infile.getline ( inputstring , MAXFILESTRINGSIZE ) ;
	} else {
	  cout << "ERROR (IOinput): readin of PEF file failed.  Verify Path and name.  Exiting pef_read.... " << endl;
	  return (RETURN_ERROR);
	  //exit(0);
	}

if (cmdtype == 6 ) strcpy ( cmd81vs6 , "CMD,6CHG_SC_TM_IMM") ;
if (cmdtype == 81 ) strcpy ( cmd81vs6 , "CMD,81") ;

//cout << cmdtype << ")" << cmd81vs6 << ":"<< strlen (cmd81vs6) ;

 	while ( infile.good ()  ){

// strpbrk: returns from delimiter on including delimiter (first character based..do not use)
// strtok: return upto delimiter, NOT including delimiter (returns to itself (char1))
// strstr: retrurn after delimiter, including delimeter (returns only to return)

infile.getline ( inputstring , MAXFILESTRINGSIZE ) ;
strcpy ( tempinputstring, &inputstring[37]);

if ( strncmp ( tempinputstring, cmd81vs6, strlen (cmd81vs6) ) == 0 ) {

	if (strlen (inputstring) > 10 ) { //there is a string no blank last line..process
		strcpy ( timestring , inputstring ) ;
	
	// get time before :
		strtok ( timestring , timestringend ) ;
	// get command after 'CMD,' prep for action string
        // first check delimiter exist, else this will create string with no end
           if (strstr (inputstring , cmddelimiter ) != NULL ) { 
		strcpy ( commandstring, strstr ( inputstring , cmddelimiter ) );
		strcpy ( commandstring, &commandstring[strlen(cmddelimiter)+1] ) ;
			strcpy ( stringtoparse , commandstring ) ;
		strtok ( commandstring , commandstringend ) ;

                // first check delimiter exist, else this will create string with no end
                if (strstr (stringtoparse , cmddelimiter ) != NULL ) {
			strcpy ( actionstring, (strstr ( stringtoparse , cmddelimiter )) );
			strcpy ( actionstring, &actionstring[strlen(cmddelimiter)+1] ) ;
			strtok ( actionstring , actionstringend ) ;

//cout  << cmdlinenumber<< timestring << "::_::" << commandstring << "::_::" << actionstring << "::_::" << "::::" ;

         strcpy (current_sclk , timestring );
         strcpy (current_command , commandstring );
         strcpy (current_action , actionstring );

         sclk[cmdlinenumber] = (double) atof (current_sclk) ;
         command[cmdlinenumber] = currentcommandnumber () ;
         action[cmdlinenumber] = currentactionnumber () ;

	 //cout <<  (unsigned int)sclk[cmdlinenumber] << "::_::" << command[cmdlinenumber] << "::_::" << action[cmdlinenumber] << "::_::" << "that's it" << endl;

	 cmdlinenumber ++; 

                } else cout << "Warning; stringtoparse-CMDDELIMITER check not found, skip this line: " << stringtoparse << endl;
            } else cout << "Warning; inputstring-CMDDELIMITER check not found, skip this line: " << inputstring << endl;

	 } else  {	 // end if strlen
		
		cout << "end of PEFCUT_FILES file found: " << cmdlinenumber << endl ;
	}

//cout << sclk[0] << "::_::" << command[0] << "::_::" << action[0] << "::_::" << "re" << endl;
} // --- end input string IF
       } // ----- END OF WHILE BASEFILE GOOD LOOP 

 if (cmdlinenumber >= MAXNUMBEROFRECORDS ) {
	cout << "********* PREINPUT ERROR *************" << endl;
	cout << "preprocess is set up (.h file) for a maximum of " << endl;
	cout << MAXNUMBEROFRECORDS << " HHHrecords.  The file: (" << filename << ")" << endl;
	cout << "has at least " << cmdlinenumber << " records." << endl;
	cout << "The preprocess .h file neeed to be modified to handle this large file!" << endl;
	cout << "********* PREINPUT ERROR *************" << endl;
}

	//cout << "END readin PEFCUT_FILES: " <<  filename  << endl ;

	numberofcommands = cmdlinenumber;
	return (cmdlinenumber) ;
	}

   // end of PEFCUT_FILES definition

//----------------------------------------------------


/* **************************************
   during get(time, command, action):
	- return a string at that record number
****************************************** */

   int PEFCUT_FILES::getnumberofcommands (  )
	{ 
	return (numberofcommands);
	}

   void PEFCUT_FILES::printcommands ( void )
	{ int loop= 0;

	cout << "\n ******************** GOING TO PRINT ********************" << endl;

		for (loop = 0; loop < numberofcommands ; loop++){
			cout << "command "<<loop<<"="<< sclk[loop]<<"=";
			cout << command[loop]<<":"<<action[loop]<<endl;
		}
	cout << "\n ******************** EXIT PRINTING ********************" << endl;

	}

   double PEFCUT_FILES::gettime ( int record )
	{ 
	return (sclk[record]);
	}


   int PEFCUT_FILES::getcommand ( int record )
	{ 
	return ( command[record] );
	}


   int PEFCUT_FILES::getaction ( int record )
	{ 
	return ( action [record] );
	}

// ------------------------------------------
// ------------------------------------------


/* ********************************************* */
/* only used internally when readin file process */

   int PEFCUT_FILES::currentcommandnumber ( void )
	{ 

	int commands_enumed = COMMANDERROR;

	if ( strncmp ("81PS_DSS", current_command, 7 ) == 0 ) commands_enumed = PS_DSS ;
	if ( strncmp ("81PS_RFES", current_command, 8 ) == 0 ) commands_enumed = PS_RFES ;
	if ( strncmp ("81PS_ESS", current_command, 7 ) == 0 ) commands_enumed = PS_ESS ;
	if ( strncmp ("81IEB_TRIGGER", current_command, 13 ) == 0 ) commands_enumed = IEB_TRIGGER ;
	if ( strncmp ("81IEB_HALT", current_command, 10 ) == 0 ) commands_enumed = IEB_HALT ;
	if ( strncmp ("81RT_RESET", current_command, 10 ) == 0 ) commands_enumed = RT_RESET ;

	if ( strncmp ("6CHG_SC_TM_IMM", current_command, 14 ) == 0 ) commands_enumed = CHG_SC_TM_IMM ;
	if ( strncmp ("6POL_SCI_CNTL", current_command, 13 ) == 0 ) commands_enumed = POL_SCI_CNTL ;
	if ( strncmp ("6MOD_POL_TBL", current_command, 12 ) == 0 ) commands_enumed = MOD_POL_TBL ;

	return (commands_enumed);
	}

// ------------------------------------------

   int PEFCUT_FILES::currentactionnumber ( void )
	{ 
	int actions_enumed = NODATA;

	if ( strncmp ("OFF", current_action, 3 ) == 0 ) actions_enumed = OFF ;
	if ( strncmp ("ON", current_action, 2 ) == 0 ) actions_enumed = ON ;
	if ( strncmp ("NORM_IEB", current_action, 8 ) == 0 ) actions_enumed = NORM_IEB ;

	if ( strncmp ("RESET", current_action, 5 ) == 0 ) actions_enumed = RESET ;
	if ( strncmp ("RELEASE", current_action, 7 ) == 0 ) actions_enumed = RELEASE ;

	if ( strncmp ("SAF142", current_action, 8 ) == 0 ) actions_enumed = LOWDATA ;
	if ( strncmp ("SAF248", current_action, 8 ) == 0 ) actions_enumed = HIGHDATA ;
	if ( strncmp ("S_N_ER_3", current_action, 8 ) == 0 ) actions_enumed = NODATA ;
	if ( strncmp ("S_N_ER_5A", current_action, 9 ) == 0 ) actions_enumed = ULTRALOWDATA ;
	if ( strncmp ("S_N_ER_8", current_action, 8 ) == 0 ) actions_enumed = HIGHDATA ;

	return (actions_enumed);
	}


// ------------------------------------------
// ---------- END OF PEF CLASS ----------------
// ------------------------------------------


// ------------------------------------------
// --------- NOW FOR SMT CLASS ------------
// ------------------------------------------

/* **************************************
   during initialize:
	- set all parameters to NA or 0
****************************************** */

   void SMT_FILE::initialize ( void )
	{ 
	int loop ;				
		for (loop = 0 ; loop < MAXIEBS ; loop++ ) {
			bits[loop] = 0 ;
			starttime[loop] = "" ;
			endtime[loop]  = "";
		} //end for loop
	allocations = 0;
	}


// ------------------------------------------

/* **************************************
   during readin:
	- set all parameters to NA 
****************************************** */

  int SMT_FILE::readin ( const char * filename  )
						// readin file...
	{
	char inputstring [ MAXFILESTRINGSIZE ];
	//char starttime_string [12];
	//char endtime_string [12];
	char bits_string [12];
	char temp_string [12];
	string year = "0000";
	string start_day = "000";
	string doy = "000";
	string hh = "00";
	string mm = "00";
	string ss = "00.000";

	ifstream infile ( filename );


 	if ( infile.good () ) {
	  ;
	} else {
	 
	  cout << "ERROR (IOinput): readin of SMT file failed.  Verify Path and name.  Exiting smt_read.... " << endl;
	  return (RETURN_ERROR);
	  //exit(0);
	}
	
 	while ( infile.good () ) {		// START WHILE while loop (file is good)
	  				 	// should get header and place it in a report
	 infile.getline ( inputstring , MAXFILESTRINGSIZE ) ;

/* **************** START: FOUND THE DATA VOLUME SECTION IN SMT.REPORT ******************* */
	 if ( strncmp ( inputstring, " DATA VOLUME REPORT", 19 ) == 0 ) { 	
								// found a DATA VOLUE to start

	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be  ----------
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // this should have RADAR in column #X
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be(units)(MB)
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be ----------
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be DATA!!

//cout << "Smt line:"<<inputstring << endl;

	   while ( strncmp ( inputstring, "--------------", 12 ) != 0) {
						// START WHILE while loop (file is in DATA VOLUME)

		 if ( strncmp ( inputstring, " OBSERVATION", 12 ) == 0 ) { 	
								// found a OBSERVATION to start
								// SEARCH FOR RADAR DATA VOLUME ALLOCATION


                    strncpy ( bits_string, &inputstring[100], 6);

		    bits[allocations] = (double) ( atof (bits_string) * 1000 * 1000 ) ;

		    //cout << starttime[allocations] << "::" << endtime[allocations] << "::" << bits[allocations] << endl;

								// REGISTER DATA IN CLASS IF VOLUME IS > 0;
		    if ( bits[allocations] > 0 ) {		// RADAR ALLOCATION!!
					

		      //strncpy ( starttime_string, &inputstring[29], 10);
		      //strncpy ( endtime_string, &inputstring[40], 10);

			strncpy ( temp_string, &inputstring[29], 10);
			strtok ( temp_string, " ");
			doy = temp_string;
// ---------------------
// added for year roll over with SMT file
// seen if DOY wrap....
// added 11/11/2004 (post build 5)
int int_doy, int_start_day, int_year;

		int_doy = atoi ( doy.c_str() );
		int_start_day = atoi ( start_day.c_str() ) ;
		int_year = atoi (year.c_str());

		if (int_doy < int_start_day ) int_year++;

		sprintf (temp_string, "%i", int_year);
		year = temp_string;
// ------------------

			strncpy ( temp_string, &inputstring[33], 10);
			strtok ( temp_string, ":");
			hh = temp_string;

			strncpy ( temp_string, &inputstring[36], 10);
			strtok ( temp_string, " ");
			mm = temp_string;

			starttime[allocations] = year+"-"+doy+"T"+hh+":"+mm+":"+ss;

			strncpy ( temp_string, &inputstring[40], 10);
			strtok ( temp_string, " ");
			doy = temp_string;

			strncpy ( temp_string, &inputstring[44], 10);
			strtok ( temp_string, ":");
			hh = temp_string;

			strncpy ( temp_string, &inputstring[47], 10);
			strtok ( temp_string, " ");
			mm = temp_string;

			endtime[allocations] = year+"-"+doy+"T"+hh+":"+mm+":"+ss;
		
			allocations ++ ;
		    }						// END RADAR ALLOCATION AND STORAGE		
		} 						// END OBSERVATION SERACH

 		infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // get next line

	    }						// END WHILE leaving DATA VOLUME START section

	 
	  }						// END IF that found a DATA VOLUME START section



	 if ( strncmp ( inputstring, " TELEMETRY MODE REPORT", 19 ) == 0 ) { 	
								// found a TELEMETRY MODE (just to get Year value)

	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be  ----------
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // this should have SCET TELEMETRY MODE AND REQUEST
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be ----------
	  infile.getline ( inputstring , MAXFILESTRINGSIZE ) ; // should be DATA!!

	  strncpy (temp_string, &inputstring[6],5);
	  strtok (temp_string, "T");
	  start_day = temp_string;
	//cout << "HHHHHHHHHHHHHHH tempsting DAY: " << temp_string<< endl;
	  strtok (inputstring, "-");
	  year = inputstring;
	 }

	}	// END WHILE while loop (file no longer good)
/* **************** END: FOUND THE DATA VOLUME SECTION IN SMT.REPORT ******************* */
	if ( atoi(year.c_str()) < 2004 || atoi(year.c_str()) > 2030 ){
	    cout << "ERROR: did not find year. Did SMT format change?" << endl;
	    cout << "ERROR: year expected between (2004-2030).  WAS:" << year << endl;
	}

	cout <<  "END readin SMT_FILE: " <<  filename << endl ;
	cout << "\tTOTAL Data Volume allocation found are: " << allocations << endl ;


	if (allocations == 0 ) {			// ERROR no allocations found in SMT
		cout << "\n*********SMT ERROR (readin): no allocations found in SMT file ***** " << endl;
	}


	return (allocations) ;
	}			// END of readin method for SMT_FILE



// ------------------------------------------


   void SMT_FILE::printcommands ( void )
	{ int loop = 0;

	cout << "\n ******************** GOING TO PRINT ********************" << endl;

		cout << "SMT: number of Data Volume allocations = " << allocations << endl;
		for (loop = 0; loop < allocations ; loop++){
			cout << "SMT("<<loop<<") = "<< bits[loop] ;
			cout << " from: "<< starttime[loop] ;
			cout << " to " << endtime[loop]<<endl;

		}
	cout << "\n ******************** EXIT PRINTING ********************" << endl;

	}



// ------------------------------------------

   double SMT_FILE::getbits ( int record )
	{ 
	if (record <= allocations ) return (bits[record]);
			else return (0);
	}


   string SMT_FILE::getstarttime ( int record )
	{ 
	if (record <= allocations ) return ( starttime[record] );
			else return ("");

	}


   string SMT_FILE::getendtime ( int record )
	{ 
	if (record <= allocations ) return ( endtime [record] );
			else return (0);

	}


   int SMT_FILE::getallocations ( void )
	{ 
	return (allocations );
	}

// ------------------------------------------



// ------------------------------------------
// ---------- END OF SMT CLASS ----------------
// ------------------------------------------

// ------------------------------------------
// --------- NOW FOR PDT_CONFIG_FILE CLASS ------------
// ------------------------------------------

//----------------------------------------------------
// --------- SMT input file definition
//----------------------------------------------------

/* **************************************
   during initialize:
	- set all parameters to NA or 0
****************************************** */

   void PDT_CONFIG_FILE::initialize ( void )
	{ 
/*	int loop ;				
		for (loop = 0 ; loop < MAXIEBS ; loop++ ) {
			bits[loop] = 0 ;
			starttime[loop] = 0 ;
			endtime[loop] = 0 ;
		} //end for loop
	allocations = 0;
*/
	filenames_found = 0 ;
	}


// ------------------------------------------

/* **************************************
   during readin:
	- set all parameters to NA 
****************************************** */

  int PDT_CONFIG_FILE::readin ( const char * filename  )
						// readin file...
	{
	char inputstring [ MAXFILESTRINGSIZE ];
	int keywordloop;
	int keyword_found = -1 ;
	int got_keywords = RETURN_ERROR ;
	int keyword_set[7] = {0,0,0,0,0,0,0} ;
	char keywords[7][30] = { 	"Spacecraft_Clock_clock_file", 	"Spacecraft_Clock_sclk_file",
					"Ephemerides_spk_files",	"Ephemerides_naif_leapseconds",
					"Ephemerides_naif_constants",   "I_KERNELS_frame_kernel",
					"I_KERNELS_radar_kernel"  } ;

        char stringtoparse1 [ MAXFILESTRINGSIZE ];
        char stringtoparse2 [ MAXFILESTRINGSIZE ];
	// three ephemeris file names are listed on 1 line on PDT.cfg
	int ephemeris_string_lenght = 0;
	bool sk_found = false;
	bool plteph_found = false;
	bool se_sat_found = false;
	bool ephemeris_used = false ; // used as warning message if 'extra' ephermis files are in pdt.cfg

	ifstream infile ( filename );

 	if ( infile.good () ) {
	  ;
	} else {
	  cout << "ERROR (IOinput): readin of PDT file failed.  Verify Path and name.  Exiting.... " << endl;
	  exit(0);
	}
	
 	while ( infile.good () ) {		// START WHILE while loop (file is good)
	  				 	// should get header and place it in a report
	 infile.getline ( inputstring , MAXFILESTRINGSIZE ) ;

/* **************** START: ANY NON-COMMENT LINE IN PDT CONFIG FILE **************************** */
	 if ( strncmp ( inputstring, "#", 1 ) != 0 ) { 	
								// found a NON-COMMENT LINE to start
//cout << "PDT (non-comment) line:"<<inputstring << endl;

   	for ( keywordloop = 0 ; keywordloop < 7 ; keywordloop++ ){
	 	if ( strncmp ( inputstring, keywords[keywordloop] , strlen (keywords[keywordloop]) ) == 0 ) { 
			strcpy ( stringtoparse1, ( strpbrk ( inputstring , " " )) );// move to next ']'
			keyword_found = keywordloop ;
		}	
	}		// end of keyword search (FOR LOOP)


switch (keyword_found){
case (0):
	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
	}
	keyword_set[keyword_found] = RETURN_OK ;
	PDT_sclk_file = stringtoparse1;
	break;

case (1):
	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
	}
	keyword_set[keyword_found] = RETURN_OK ;
	PDT_sclkscet_file = stringtoparse1;
	break;

case (2):					// special:  will have multiple files.....
 	while (strstr (stringtoparse1, "/") != NULL ) {			// gets the last file after <>/
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
	}
	//PDT_solarsystem_ephemeris_file = stringtoparse1;

	strcpy ( stringtoparse1, ( strpbrk ( inputstring , "{" )) );// move to next '{'
	ephemeris_string_lenght = strlen (stringtoparse1);

	// should file 3 items in this loop
	while ((!sk_found || !plteph_found || !se_sat_found) && ephemeris_string_lenght > 4 ) { // start of while loop

	  //strcpy ( stringtoparse1, ( strpbrk ( inputstring , "{" )) );// move to next '{'
	strcpy ( stringtoparse1, ( strpbrk ( stringtoparse1 , " " )) );// move to next space
	strcpy (stringtoparse1, &stringtoparse1[1] );
	// Have string of file names
	strcpy (stringtoparse2, stringtoparse1 );		// for later...second word
	strtok ( stringtoparse2 , " " ) ;				// read up to next space
	strtok ( stringtoparse2 , "}" ) ;				// or read up to end bracket

	ephemeris_used = false; // used to identify if filename is not used (extra files in cfg file!)

	if (strstr (stringtoparse2, "_SK_") != NULL){ // found SK in file
	  PDT_sk_ephemeris_file = stringtoparse2;
	  ephemeris_used = true;
	  //cout << "in SK...." << endl;
	  if (sk_found) 
	    cout << "ERROR:(IOinput): second SK ephemeris file found in PDT_CFG." << PDT_sk_ephemeris_file << " & " << stringtoparse2 <<endl ; 
	  else sk_found = true;
	} else if ( (strstr (stringtoparse2, "_PLTEPH") != NULL) || (strstr (stringtoparse2, "_PE_") != NULL) ){ // found PLTEPH in file
 	  PDT_solarsystem_ephemeris_file = stringtoparse2;
	  ephemeris_used = true;
	  //cout << "in PLTEPH...." << endl;
	  if (plteph_found) 
	    cout << "ERROR:(IOinput): second PLTEPH ephemeris file found in PDT_CFG." << PDT_solarsystem_ephemeris_file << " & " << stringtoparse2 <<endl ; 
	  else plteph_found = true;
	} else if (strstr (stringtoparse2, "_SE_") != NULL){ // found SE_SAT in file
	  PDT_planet_ephemeris_file= stringtoparse2;
	  ephemeris_used = true;
	  //cout << "in SOLAR(SE)...." << endl;
	  if (se_sat_found) 
	    cout << "ERROR:(IOinput): second SE_SAT ephemeris file found in PDT_CFG." << PDT_planet_ephemeris_file << " & " << stringtoparse2 <<endl ; 
	  else se_sat_found = true;
	}

	// debug
//	cout << "sk/plt/sol: " <<  PDT_sk_ephemeris_file <<":"<<PDT_planet_ephemeris_file<<":"<<PDT_solarsystem_ephemeris_file<< endl;
//	cout << "true Flags:"  << sk_found <<":"<< plteph_found<<":"<<se_sat_found << endl;
//	cout << "OLDstr: " << stringtoparse1 << endl;

	if (!ephemeris_used){
	  cout << "Warning(IOinput): Extra ephemeris filename: ("<< stringtoparse2<<")  File not used.  User should verify." << endl;
	} 
	
	
  //      ephemeris_string_lenght = strlen (stringtoparse1);
//	cout << "early string lenght = " << ephemeris_string_lenght;

	strcpy ( stringtoparse1, ( strpbrk ( stringtoparse1 , " }" )) );// move to next space
	ephemeris_string_lenght = strlen (stringtoparse1);
	
	// debug
//	cout << "NEWstr: " << stringtoparse1 << endl;
//	cout << "StringLEN2 = " << ephemeris_string_lenght << endl;

	} // end of string_lenght and found(flag) while loop
	
	if ( sk_found && plteph_found && se_sat_found ) // if all three found the PDT.cfg is OK for eph
	  keyword_set[keyword_found] = RETURN_OK;

	break;

case (3):
 	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
	}
 	keyword_set[keyword_found] = RETURN_OK ;
	PDT_time_kernel_file = stringtoparse1;
	//cout << "IOinput: case3: time kernel found";
	break;

case (4):
/*	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
        	cout << "case4(loop): "<< stringtoparse1<<endl;
	}
*/
      while (strstr (stringtoparse1, "/") != NULL ) {
                strcpy (stringtoparse1, (strstr (stringtoparse1, "/pck/")));
                strcpy (stringtoparse1, &stringtoparse1[5] );
		strcpy (stringtoparse2, stringtoparse1);//copy 1 to 2
		strtok ( stringtoparse2 , " }" ) ;                               // read up to next space
		if (strstr (stringtoparse2, "ock") == NULL) { // pck has normal and 'rock' version, we mant to skip rock and get normal
			// this means stringtoparse2 has non-rack file. We are done.
                        if (strstr (stringtoparse1, "/pck/") != NULL ){// there are still more files
                		strcpy (stringtoparse1, (strstr (stringtoparse1, "/pck/")));
				cout << "Warning(IOinput): There are Extra cpk files ("<< stringtoparse1<<") NOT used.  User should verify." << endl;
			}
                	strcpy (stringtoparse1, stringtoparse2);//place answer back in #1
		} else { // wrong file...keep looking
          		cout << "Warning(IOinput): Extra cpk file: ("<< stringtoparse2<<")  File not used.  User should verify." << endl;
		}
	}

	keyword_set[keyword_found] = RETURN_OK ;
	PDT_planetary_constants_file = stringtoparse1;
        //cout << "IOinput: case4: planetary file found"<< PDT_planetary_constants_file;
	break;

 case (5): // frame kernel
	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
		//cout << "FRAME(loop) = " << stringtoparse1 << endl;
		if (strstr (stringtoparse1,"cas_v") != NULL) {
			//cout << "got it" << stringtoparse1 << endl;
			strcpy (stringtoparse2, stringtoparse1);
			strtok (stringtoparse2, " ");
			//cout << "parse2 = " << stringtoparse2 << endl;
		} 
	}
	if (strstr (stringtoparse1, "cas_stat") != NULL)
		cout << "Warning(IOinput): Extra frame file: ("<< stringtoparse1<<")  File not used.  User should verify." << endl;
	if (strstr (stringtoparse2, "cas_v") != NULL){
		keyword_set[keyword_found] = RETURN_OK ;
		PDT_frame_kernel_file = stringtoparse2;
		//cout << "FRAME = " << PDT_frame_kernel_file << endl;
	}
	break;

 case (6): // instrument kernel
	while (strstr (stringtoparse1, "/") != NULL ) {
		strcpy (stringtoparse1, (strstr (stringtoparse1, "/")));
		strcpy (stringtoparse1, &stringtoparse1[1] );
	}
	keyword_set[keyword_found] = RETURN_OK ;
	PDT_instrument_kernel_file = stringtoparse1;
	//cout << "INSTR = " << PDT_instrument_kernel_file << endl;
	break;

default:
	break;
}		// end of SWITCH statement on keyword_found

keyword_found = -1;			// clear keyword found for next loop


	 
	  }						// END IF that found a NON-COMMENT LINE section

	}	// END WHILE while loop (file no longer good)
 
got_keywords = RETURN_OK ;
for ( keywordloop = 0 ; keywordloop < 5 ; keywordloop++){
	if ( keyword_set[keywordloop] != RETURN_OK ) got_keywords = RETURN_ERROR ;
}

	return (got_keywords) ;
	}			// END of readin method for PDT_CONFIG_FILE



// ------------------------------------------


   void PDT_CONFIG_FILE::printcommands ( void )
	{
	cout << "\n ******************** GOING TO PDT_CONFIG PRINT ********************" << endl;

		cout <<  "\tPDT_sclk_file: " << get_sclk ( ) << endl ;
		cout <<  "\tPDT_sclkscet_file: " << get_sclkscet ( ) << endl ;
		cout <<  "\tPDT_sk_ephemeris_file: " << get_sk_ephemeris ( ) << endl ;
		cout <<  "\tPDT_planet_ephemeris_file: " << get_planet_ephemeris ( ) << endl ;
		cout <<  "\tPDT_solarsystem_ephemeris_file: " << get_solarsystem_ephemeris ( ) << endl ;
		cout <<  "\tPDT_time_kernel_file: " << get_time_kernel ( ) << endl ;
		cout <<  "\tPDT_planetary_constants_file: " << get_planetary_constants ( ) << endl ;

	cout << "\n ******************** EXIT PDT_CONFIG PRINTING ********************" << endl;
	}


// ------------------------------------------

   string PDT_CONFIG_FILE::get_sclk ( )
	{ 
	return ( "sclk/"+PDT_sclk_file );
	}

   string PDT_CONFIG_FILE::get_sclkscet ( )
	{ 
	return ( "sclk/"+PDT_sclkscet_file );
	}

   string PDT_CONFIG_FILE::get_sk_ephemeris ( )
	{ 
	return ( "sp/"+PDT_sk_ephemeris_file );
	}

   string PDT_CONFIG_FILE::get_planet_ephemeris ( )
	{ 
	return ( "sp/"+PDT_planet_ephemeris_file );
	}

   string PDT_CONFIG_FILE::get_solarsystem_ephemeris ( )
	{ 
	return ( "sp/"+PDT_solarsystem_ephemeris_file );
	}

   string PDT_CONFIG_FILE::get_time_kernel ( )
	{ 
	return ( "lsk/"+PDT_time_kernel_file );
	}

   string PDT_CONFIG_FILE::get_planetary_constants ( )
	{ 
	return ( "pck/"+PDT_planetary_constants_file );
	}

   string PDT_CONFIG_FILE::get_frame_kernel ( )
	{ 
	return ( "fk/"+PDT_frame_kernel_file );
	}

   string PDT_CONFIG_FILE::get_instrument_kernel ( )
	{ 
	return ( "ik/"+PDT_instrument_kernel_file );
	}


// ------------------------------------------
// ---------- END OF PDT_CONFIG_FILE CLASS ----------------
// ------------------------------------------

//----------------------------------------------------
//----------------------------------------------------
//----------------------------------------------------

