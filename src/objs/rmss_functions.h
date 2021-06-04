
//			RMSS.HEADER
// ----------------------------------------------------------
//			rmss.h
// 	Date: 				March 04, 2002 
//	Author:				Gary Hamiton
//	Subsystem: Cassini Radar Mapping Sequencing Software (RMSS)
//
// ----------------------------------------------------------

#ifndef __RMSS_H
#define __RMSS_H

	// -----------------------------------------------
	//	Version Control
	// -----------------------------------------------

#define 	RMSS_H_VERSION_MM 	02
#define 	RMSS_H_VERSION_DD 	24
#define 	RMSS_H_VERSION_YYYY 	2003

	// -----------------------------------------------
	//	BASE INPUT Definitions
	// -----------------------------------------------

#define MAX_INPUT_SIZE			100
#define	INPUT_SECONDS_STRING_SIZE	10

#define NUMBER_OF_BASE_PARAMETERS 9 	// elements + space for filename
#define MAXSIZE_OF_BASE_PARAMETER 30

#define NUMBER_OF_ENG_PARAMETERS 23	// elements + space for filename
#define MAXSIZE_OF_ENG_PARAMETER 30

#define NUMBER_OF_SCI_PARAMETERS 42 	// elements + space for filename
#define MAXSIZE_OF_SCI_PARAMETER 30

#define NUMBER_OF_SCIMODE_PARAMETERS 9 	// elements + space for filename
#define SCIMODE_INTFIELDS_START 2 	// first 2 element are characters, then the rest are integers


// POWER LEVELS ---> moved to IOinput.h
//#define OFF_POWER_LEVEL 0
//#define LOW_POWER_LEVEL 55
//#define HIGH_POWER_LEVEL 83

//BASE_PARAMETERS 
	#define b_BASE_FILE	0
	#define b_COMMENT	1
  	#define b_TARGET	2
	#define b_PASS		3
	#define b_RUN		4
	#define b_TRAJ_FILE	5
	#define b_ENG_FILE	6
	#define b_SCI_FILE	7
	#define b_SCIMODE_DIR	8


//ENG_PARAMETERS
	#define e_ENGFILE			0
	#define e_COMMENT			1
	#define e_TRIGGER_TIME			2
  	#define e_IEB_DURATION			3
	#define e_DATA_VOLUME			4
	#define e_DATA_RATE_START		5
	#define e_DATA_RATE_END			13
	#define e_POWER_ALLOCATION_START	14
	#define e_POWER_ALLOCATION_END		22



//SCI_PARAMETERS
	#define s_SCIFILE			0
	#define s_COMMENT			1
	#define s_SCI_MODE_START		2
	#define s_SCI_MODE_END			41


//SCI_PARAMETERS
	#define m_SCIMODEFILE			0
	#define m_COMMENT			1
	#define m_MOD				2
	#define m_CSR				3
	#define m_ADC				4
	#define m_RCV				5
	#define m_TRO				6
	#define m_BAQ				7
	#define m_BEM				8




	// -----------------------------------------------
	//	Global Definitions
	// -----------------------------------------------


#define 	DNE		-99

#define 	PRINTSTRINGSIZE		200
#define 	DIR_STRING_SIZE		180



#define 	READ			1001
#define		WRITE			1002
#define 	TRUE			1
#define 	ENABLE			1
#define		DISABLE			0
#define		FALSE			0

#define 	OPEN			10
#define 	APPEND			15
#define		CLOSE			20

//#define 	NOT_VALID		0x5a
//#define 	VALID			0xa5
//#define 	RETURN_FAULT		0x5a
//#define 	RETURN_OK		0xa5


#define		MAX_IEB_TIME		0xffff
// (SFTI is 16 bits = 65535 seconds)

	// -----------------------------------------------
	//	IEB Parameters
	// -----------------------------------------------

#define 	SLOW_MAXSIZE 		200
#define 	SLOW_INSTRUCTIONS	19
#define 	SLOW_WORDS 		10
#define	 	SLOW_TYPE		3

#define 	FAST_MAXSIZE 		500
#define 	FAST_INSTRUCTIONS	19
#define 	FAST_WORDS 		7
#define 	FAST_TYPE		2

#define 	PWR_MAXSIZE 		20
#define 	PWR_INSTRUCTIONS	3
#define 	PWR_WORDS 		2
#define 	PWR_TYPE		0

#define 	TNC_MAXSIZE 		20
#define 	TNC_INSTRUCTIONS	6
#define 	TNC_WORDS 		5
#define 	TNC_TYPE		1

#define 	EDS_TYPE		5
#define 	ENG_TYPE		6
#define 	TRAJ_TYPE		7



	// ---------  DEFAULT PARAMETERS 
#define 	DEFAULT_ND_IP		2
#define 	DEFAULT_RL_IP		5
#define 	ND_MAX			15
#define 	RL_MAX			15
#define 	ATTN_MAX		65
#define 	BEAM_MASK1_2		3
#define 	BEAM_MASK3		4
#define 	BEAM_MASK4_5		24


#define 	DEFAULT_CSD		5
#define 	ATTN_MAX		65
#define 	BEAM_MASK1_2		3









/* POWER STATES */
/*
#define  enum POWER_STATES 
[WU1-WU2,WU2-WU3,WU3-WU4,WU4-RAI
RAI-WU4,WU4-WU3,WU3-WU2,WU2-WU1,WU1-WU0
WU0-WU1]
*/
/*  put general structures here*/


// Parameter passing strucrure.............
typedef struct {
	int request		;
	int time		; 
	char answer	[50]	;
} PASSING_STRUCTURE ;



// main IEB strucrure.............
typedef struct {

/* SLOW PARAMETERS: 160 bits or ten 16-bit words */
	 int STFI			; 
	 	int TFI_h		;
		int TFI_m		; 
		int TFI_s		;
	 int STYP			;
	 int DTN			;
	 int SIN			;
	 int MOD			; 
	 	char MOD_char[6]	;
	 int CSR			; 
	 	char CSR_char[6]	;
	 int ADC			; 
	 	int ADC_khz		;
	 int RCV			; 
	 	int RCV_khz		;
	 int TRO			; 
	 	int TRO_pri		;
	 int BAQ			; 
	 	char BAQ_char[10]	;
	 int BEM			;
	 int AT1			; 
		 int AT1_db		; 
	 int AT3			; 
		 int AT3_db		;
	 int AT4			; 
		 int AT4_db		;
	 int RIP			; 
	 	int RIP_msec		;
	 int RAD			;
	 int CSD			; 
	 	int CSD_nsec		;
	 int CSQ			; 
	 	int CSQ_usec		;
	 int CFS			; 
	 	int CFS_khz		;


/* FAST PARAMETERS: 112 bits or seven 16-bit words */
 int FTFI				;
	 int FTYP			;
	 int FIN			;
	 int BII			;
	 int PUL			; 
	 	int B_on		;
	 int BPD			;
	 int RWD			; 
	 	int RWD_msec	;
	 int PRI			;
	 	float XMTR_percent	; 
	 	float XMON_percent	;
	 int CSF			;
	 	int CSF_khz		; 
	 	int CSF_center_freq	;


/* POWER PARAMETERS: 32 bits or two 16-bit words */
	 int PTFI			;
	 int PTYP			;
	 int PMD			; 
	 	char PMD_char[10]	;

/* TNC PARAMETERS: 80 bits or five 16-bit words */
	 int TTFI			;
	 int TTYP			;
	 int TNC			;
	 int TCA			;
	 int TCB			;
	 int TCC			;


/* EDS INFORMATION  240 bits or Fifteen 16-bit words */	

	 int RECF			; 
	 int DEL1			; 
	 int DEL2			; 
	 int EDS_spare1			;
	 int DELB			;
	 int SCL1			;
	 int SCL2			;
	 int NCOF			;
	 int UPCN			;
	 int ESIN			;
	 int EFIN			;
	 int EDSM			;
	 int ATSM			;
	 int ATSA			;
	 int RATN			;
	 int EDS_DATA_VOLUME		;


/* ENGINEERING/SCIENCE PARAMETERS  240 bits or Seventeen 16-bit words*/
	 int VALIDTYPE		;

	 int AUTO_RAD_ENABLE	;
	 int AUTO_GAIN_ENABLE	; 
	 int CTU_TAIL_ENABLE	;
	 int ENG_TAIL_ENABLE	;
	 int PICK_UP_MODE		;   
	 int ENG_FLAG3		; 
	 int ENG_FLAG2		; 
	 int ENG_FLAG1		;
 
	 int ND_IP			;
	 int RL_IP			; 
	 int AUTO_ND_MAX		;
	 int AUTO_RL_MAX		;

	 int PRODUCTION_RATE	;
	 int SDB_WORDS		;

	 int SABS_PRODUCED 	;
	 int SAB_SIZE		;

	 int CUM_DATA_VOLUME	;
 



/* TRAJECTORY INFORMATION  */
	 
	 int TARGET_RANGE	;
	 int TARGET_ANGLE	;
	 int TRAJ_spare8	; 
	 int TRAJ_spare7	; 
	 int TRAJ_spare6	; 
	 int TRAJ_spare5	; 
	 int TRAJ_spare4	;
	 int TRAJ_spare3	; 
	 int TRAJ_spare2	;  
	 int TRAJ_spare1	;



	} ARRAY_STRUCTURE ;




/* FUNCTION PROTOTYPES */
/* void to_report (*char);*/
 int     to_report ( char* )   ;

/* *******************************************	*/
/* *******************************************	*/
/* *******************************************	*/


/* to go to get_constants and config file */

/* really: slow+fast+tnc+power (count start 0-739 or 740 commands)*/
#define MAX_ARRAY_SIZE                  739
#define MAX_TIME_SECONDS                64800     /* 18 hours = 64800 */
#define MIN_TIME_INTERVAL               2       /* every 2 seconds = new inst */

#define IEB_ON_TIME	               	12      /* every 2 seconds = new inst */
#define MAX_IDEAL_ARRAY_SIZE 		19020 /* if bigger c crashes  */
/* 19020 ok */
/* 19050 core dump*/


#define ALL_TYPES                       10
#define NUMBER_OF_TYPES                 7
#define TOSCREEN                        110
#define TOWTKFILE                       111
#define TOWTKXCLFILE                    112





/*  putGLOBAL structure  here*/

//typedef struct {

/* SLOW PARAMETERS: 160 bits or ten 16-bit words */
//	 char HOME_DIR[50]	;


//	} GLOBAL_STRUCTURE ;



// 	Standard function required to be called on in new versions
//	of gcc (g++)						


//#include <string>
//#include <fstream>
//#include <iostream>

using std::cout;
using std::endl;
using std::ifstream;

using std::string;
//using std::cerr;
//using std::terminate;


#endif  /* __RMSS_H defined */


// --------------------------------------------------------------
//	END OF RMSS HEADER FILE
// --------------------------------------------------------------

