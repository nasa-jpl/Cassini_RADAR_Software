
//			PRECLASSES.HEADER
// ----------------------------------------------------------
//			PRECLASSES.h
// 	Date: 				January 23, 2003
//	Author:				Gary Hamiton
//	Subsystem: Cassini Radar Mapping Sequencing Software (RMSS)
//
// ----------------------------------------------------------

#ifndef __PRECLASSES_H
#define __PRECLASSES_H

#include <string>

// -----------FUNCTION PROTOTYPES -----------
//int OutputReport ( char * ) ;
//
#define TNC_MODE	98
#define POWER_MODE	99

#define RETURN_ERROR 		9999
#define RETURN_OK 		1

#define MAXFILESTRINGSIZE 	300  	// must be multiple of 5
#define MAXNUMBEROFRECORDS 	400	// has to cover both CDS (6) and RADAR (81)
#define MAXIEBS		 	15
#define MAXPOWER		40
	#define NOPOWER		0
	#define LOWPOWER	55
	#define HIGHPOWER	86

#define MAXDATA		 	200	// must handle ALL 6CHG CDS commands
	#define NODATA		0	// any unknown mode
	#define ULTRALOWDATA	7600	// S_N_ER_5A? S_N_ER_3??
	#define LOWDATA		30400	// SAF142
	#define HIGHDATA	364800	// SAF248B, S_N_ER_8

// 81COMMAND SEARCHES
	#define COMMANDERROR		810
	#define PS_DSS			811
	#define PS_RFES			812
	#define PS_ESS			813
	#define IEB_TRIGGER		814
	#define IEB_HALT		815
	#define RT_RESET		816

// 6COMMAND SEARCHES
	#define CHG_SC_TM_IMM		650
	#define POL_SCI_CNTL		651
	#define MOD_POL_TBL		652


// 81ACTION SEARCHES
	#define ACTIONERROR		710
	#define OFF			711
	#define ON			712
	#define NORM_IEB		713
//HANDLE POWER GLITCH
	#define RESET			714	// the reset command, RESET set power to 0
	#define RELEASE			715	// the reset command, RELEASE set power to 55


// 6ACTION SEARCHES
	#define SAF142			601		// need to verify text
	#define SAF248			602		// need to verify text
	#define S_N_ER_3		610
	#define S_N_ER_5A		611
	#define S_N_ER_8		612		// need to verify text
	#define SAMEMODE		699		// to reduce duplicate CDS mode changes


// FILE IO
	#define FCLOSE			1001		// need to verify text
	#define FAPPEND			1002


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

using std::string;
using std::cerr;
using std::terminate;

typedef struct { 
		double time		;
		int level		;
} ALLOCATION ;

typedef struct { 
		double time		;
		int duration 		;
		double volume 		;
		int idap 		;

} IEB_STATS ;


// --------------- for PEF input file
class PEFCUT_FILES
	{	
	private: 
		double sclk[MAXNUMBEROFRECORDS]	;
		int command[MAXNUMBEROFRECORDS]	;
		int action[MAXNUMBEROFRECORDS]	;

		int numberofcommands		;
		int cmdtype			;
		char current_sclk[MAXFILESTRINGSIZE]		;
		char current_command[MAXFILESTRINGSIZE]		;
		char current_action[MAXFILESTRINGSIZE]		;

	public:		
		PEFCUT_FILES ( void )  {cout<<"ERROR in PEF INIT: filetype not specified"<<endl;}
		PEFCUT_FILES ( int cmdtype )  {  initialize ( cmdtype ) ; }
		void initialize ( int  ) ;
		int readin ( const char * ) ;

		int getnumberofcommands ( void );
		double gettime (int );
		int getcommand (int );
		int getaction (int );

		int currentcommandnumber ( void );
		int currentactionnumber ( void );

		void printcommands ( void );

	};
// --------------- for PEF input file
class SMT_FILE
	{	
	private: 
		double bits[MAXIEBS]	;
		string starttime[MAXIEBS]	;
		string endtime[MAXIEBS]	;
		int allocations ;

		
	public:		
		SMT_FILE ( void )  { initialize ( ) ;}
		void initialize ( void  ) ;
		int readin ( const char *  ) ;

		double getbits ( int );
		string getstarttime (int );
		string getendtime (int );
                int getallocations ( void  ) ;

		void printcommands ( void );


	};

// --------------- for PDT Config input file
class PDT_CONFIG_FILE
	{	
	private: 
		string PDT_sclk_file ;
		string PDT_sclkscet_file ;
		string PDT_sk_ephemeris_file ;
		string PDT_planet_ephemeris_file ;
		string PDT_solarsystem_ephemeris_file ;
		string PDT_time_kernel_file ;
		string PDT_planetary_constants_file ;
		string PDT_frame_kernel_file ;
		string PDT_instrument_kernel_file ;
		int filenames_found ;
		
	public:		
		PDT_CONFIG_FILE ( void )  { initialize ( ) ;}
		void initialize ( void  ) ;
		int readin ( const char *  ) ;

		string get_sclk ( );
		string get_sclkscet ( );
		string get_sk_ephemeris ( );
		string get_planet_ephemeris ( );
		string get_solarsystem_ephemeris ( );
		string get_time_kernel ( );
		string get_planetary_constants ( );
		string get_frame_kernel ( );
		string get_instrument_kernel ( );
		void printcommands ( void );


	};

// --------------- for PDT_CONFIG_FILE input file




#endif  /* __PRECLASSES_H defined */
// --------------------------------------------------------------
//	END OF __PRECLASSES_H HEADER FILE
// --------------------------------------------------------------
