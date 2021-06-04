/* 
 * definitions.h
 *
 * Generic defines and typedefs that are used across several of the 
 * processes 
 *
 */

#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

/* added 2/1/2003 by GH */

#define TRUE	1
#define FALSE	0

/** definitions file for all A', A and B software modules **/


/* internal process control and definitions */
#define MAX_SEQ_LEN 2048 
#define RTIU_COL_LEN 30
#define SEQ_STRUCT_LEN 384
#define SEQ_STRUCT_LEN16 192  /* in 16-bit words */
#define LIMIT_LEN 2048  /* in bytes */
#define DSS_YMIN_IDX 128
#define RFES_YMIN_IDX 129
#define ESS_YMIN_IDX 130
#define DSS_YMAX_IDX 384
#define RFES_YMAX_IDX 385
#define ESS_YMAX_IDX 386
#define DSS_RMIN_IDX 640
#define RFES_RMIN_IDX 641
#define ESS_RMIN_IDX 642
#define DSS_RMAX_IDX 896
#define RFES_RMAX_IDX 897
#define ESS_RMAX_IDX 898
#define DCPWR_MAX_REQ 30
#define FLAGS_BUFSIZE 1024
#define RTIU_COL_CMD 196
#define RTIU_TIME_OFFSET 10
#define DSS_RFESSIM_START 64
#define NUM_RFESSIM_BYTES 12
#define HEAPSIZE 1048576  /* 1MB allocated for memload */
#define SLOW_LEN 10 /* 16 bit */
#define FAST_LEN 7  /* 16 bit */
#define POWER_LEN 2 /* 16 bit */
#define TEST_LEN 5  /* 16 bit */
#define CMDQ_LEN 1024
#define DISPQ_LEN 128
#define ERR_ARCQ_LEN 128
#define DCPWRQ_LEN 128
#define MAX_PORTS 5
#define SEND_MSG_OPEN 0
#define SEND_MSG_READ 1
#define SEND_MSG_WRIT 2
#define SEND_MSG_CLOS 3
#define RECV_MSG_OPEN 4
#define RECV_MSG_READ 5
#define RECV_MSG_WRIT 6
#define RECV_MSG_CLOS 7
#define QUIT 127
#define RTIU_TYPE 0
#define DSSSE_SEQ_TYPE 1
#define DSSSE_IMM_TYPE 11
#define RFESSE_TYPE 2
#define COPY 0
#define INS 1
#define APP 2
#define DEL 3

#define CONNECT_TO_SERVER 0
#define SEND_TO_SERVER 1
#define CLOSE_SERVER_CONNECTION 2
#define MKEY1 ((key_t) 1234L)
#define MKEY2 ((key_t) 2345L)
#define MKEY3 ((key_t) 4321L)
#define MKEY4 ((key_t) 4331L)
#define MKEY5 ((key_t) 4341L)
#define MKEY6 ((key_t) 4351L)
#define MKEY7 ((key_t) 4361L)
#define MKEY8 ((key_t) 4371L)
#define MKEY9 ((key_t) 4381L)

#define PERMS 0666

/* Input Message Types for Test Status GUI */

#define RTC_MSG 3L
#define RTCD_MSG 4L
#define IPT_MSG 5L
#define CST_DISP_MSG 6L
#define OST_DISP_MSG 8L
#define GPIB_DISP_MSG 10L
#define DCPWR_MSG 11L
#define TELEM_DISP_MSG 12L
#define RTCD_ERR_MSG 13L

#define TIMERVAL 999000   /* usec */
#define MAP_VME 34


#define FLAGS_SHMKEY ((key_t) 7799)
#define FLAGS_SHM_OPEN 80
#define FLAGS_SHM_CLOS 81

#define RTTS_FLG_SEMKEY ((key_t) 5794)
#define RSVE_FLG_SEMKEY ((key_t) 6795)
#define IOST_FLG_SEMKEY ((key_t) 6796)
#define EOST_FLG_SEMKEY ((key_t) 6797)
#define CST_FLG_SEMKEY ((key_t) 6798)
#define RSYNC_FLG_SEMKEY ((key_t) 6799)
#define DCPWR_FLG_SEMKEY ((key_t) 6794)
#define LIMFLG1_SEMKEY ((key_t) 7786)
#define LIMFLG2_SEMKEY ((key_t) 7787)


#define IMM_SHMKEY ((key_t) 7890)
#define IMM_SEMKEY1 ((key_t) 7891)
#define IMM_SEMKEY2 ((key_t) 7892)
#define IMM_SHM_OPEN 0
#define IMM_SHM_CLOS 1
#define IMM_SHM_READ 2

#define SEQ_SHMKEY ((key_t) 7893)
#define SEQ_SEMKEY1 ((key_t) 7894)
#define SEQ_SEMKEY2 ((key_t) 7895)
#define SEQ_SHM_OPEN 10
#define SEQ_SHM_CLOS 11
#define SEQ_SHM_READ 12

#define SCT_SHMKEY ((key_t) 7777)
#define SCT_SEMKEY ((key_t) 7778)
#define SCT_SHM_OPEN 20
#define SCT_SHM_CLOS 21
#define SCT_SHM_READ 22
#define SCT_SHM_WRIT 23

#define ARC_SHMKEY ((key_t) 7780)
#define ARC_SEMKEY ((key_t) 7781)
#define ARC_SHM_OPEN 30
#define ARC_SHM_CLOS 31
#define ARC_SHM_READ 32
#define ARC_SHM_WRIT 33

#define LIM_SHMKEY ((key_t) 7783)
#define LIM_SEMKEY ((key_t) 7784)
#define LIM_SHM_OPEN 40
#define LIM_SHM_CLOS 41
#define LIM_SHM_READ 42
#define LIM_SHM_WRIT 43

#define TCST_SHMKEY ((key_t) 7788)
#define TCST_SEMKEY ((key_t) 7789)
#define TCST_SHM_OPEN 50
#define TCST_SHM_CLOS 51
#define TCST_SHM_READ 52
#define TCST_SHM_WRIT 53

#define RTTS_SHMKEY ((key_t) 7791)
#define RTTS_SEMKEY ((key_t) 7792)
#define RTTS_SHM_OPEN 60
#define RTTS_SHM_CLOS 61
#define RTTS_SHM_READ 62
#define RTTS_SHM_WRIT 63

#define SETTS_SHMKEY ((key_t) 7796)
#define SETTS_SEMKEY ((key_t) 7797)
#define SETTS_SHM_OPEN 70
#define SETTS_SHM_CLOS 71
#define SETTS_SHM_READ 72
#define SETTS_SHM_WRIT 73

#define RS_SHMKEY ((key_t) 7800)
#define RS_SEMKEY ((key_t) 7801)
#define RS_SHM_OPEN 84
#define RS_SHM_CLOS 85
#define RS_SHM_READ 86
#define RS_SHM_WRIT 87

#define IPT_MSG_OPEN  90
#define IPT_MSG_CLOS  91
#define DISPQ_MSG_OPEN 92
#define DISPQ_MSG_CLOS 93
#define ARCQ_MSG_OPEN 94
#define ARCQ_MSG_CLOS 95
#define GPIB_DATAQ_MSG_OPEN 96
#define GPIB_DATAQ_MSG_CLOS 97
#define GPIB_CMDQ_MSG_OPEN 98
#define GPIB_CMDQ_MSG_CLOS 99
#define ENGQ_MSG_OPEN 100
#define ENGQ_MSG_CLOS 101
#define HDRFTRQ_MSG_OPEN 102
#define HDRFTRQ_MSG_CLOS 103

#define SET_RSVE 74
#define RESET_RSVE 75
#define SET_RTTS 76
#define RESET_RTTS 77

#define IRIG_OPEN 0
#define IRIG_READ 1

#define GPIB_INIT 0
#define GPIB_READ_DSS 1
#define GPIB_READ_RFES 2
#define GPIB_READ_ESS 3

#define BIGCOUNT 10000 

/** Message type to be used for Test message across a hardened socket **/
#define MSG_TEST_ONLY 		150
  
/** message type definitions for Telemetry Processor interface **/
#define MSG_HOUSEKEEP 		161
#define MSG_ENGINEER 		162
#define MSG_HFT 		163
#define MSG_CTU_EXTRA 		164  /* Obsolete. Will be removed */
#define MSG_SCIENCE 		165
#define MSG_SCIENCE_PLBK 	166
#define MSG_HOUSEKEEP_PLBK 	167
#define MSG_ENGINEER_PLBK 	168
#define MSG_HFT_PLBK 		169
#define MSG_SLOW_ERR		241
#define MSG_FAST_ERR		242
#define MSG_TNCI_ERR		243
#define MSG_PWRI_ERR		244
#define MSG_BIU_DISCRETE 	245
#define MSG_SLOW_ERR_PLBK	246
#define MSG_FAST_ERR_PLBK	247
#define MSG_TNCI_ERR_PLBK	248
#define MSG_PWRI_ERR_PLBK	249

/** archive message queue type codes **/
#define HKP_ARC 		200
#define MRO_ARC 		201
#define SCI_ARC 		202
#define ENG_ARC 		203
#define HDRFTR_ARC 		204
#define RDR_LIM_ERR_ARC 	205
#define RDR_STAT_ERR_ARC 	206
#define CCSDS_RTIU_ARC  	207
#define OST_ARC 		208
#define CMD_STRT_ARC 		209
#define CCSDS_ARC 		210
#define RDR_STAT_MSG_ARC	211
#define RDR_LIM_MSG_ARC		212

#define GENERAL_TEXT		17
#define TIME_TAGGED_TEXT	18

#define RTIU_PACKET_ARC		23
#define DSS_SE_ARC		24
#define TRIG_ILX_ARC		25

#define COLLECT_ARC_DATA	220
#define PLAYBACK_ARC_DATA	221

#define MODE_COLLECTION		301
#define MODE_PLAYBACK		302

/** engineering queue data types and lengths **/
#define ENG_FULL_DATA_MSG 300
#define ENGQ_FULL_LEN 196  /** 96 words data plus 2 wd SCT  **/
#define ENG_SUBCOM_DATA_MSG 301
#define ENGQ_SUBCOM_LEN 18 /** 6 words data + 1 wd SCT & subcom  **/
#define SUBCOM_OFFSET 30  /**  byte offset to sct in sci footer **/
#define SAB_HEADER_MSG		302
#define TRIG_TIME_MSG		303
/* pointers into flags shared memory buffer */
#define LTM_A 0
#define LTM_B 1
#define RSVE 2
#define RTTS 3
#define IOST 4
#define EOST 5
#define CST 6
#define RSYNC 7
#define DCPWR_READ 8

/* IEB Instruction Type Codes */
#define SF_TYPE 3
#define FF_TYPE 2
#define TC_TYPE 1
#define PR_TYPE 0

/* Immediate Instruction Type Codes */
#define IMM_PR_TYPE 0
#define IMM_SF_TYPE 1  /* combo slow-fast instruction */
#define IMM_TC_TYPE 2


/* CDS/Radar Commands (described in 3484-CRADAR-039 Muh Yang 3/8/94) */
#define IEBHALT 0x0f
#define IEBTRIG 0x0c
#define LOADCMD 6
#define ILX 9
#define MRO 3
/*** #define MRO 1   OLD OPCODE  Did you Remember to rebuild edit_work?? ***/
#define CHECKSUM 5
#define MEMLOAD 0x20 
#define ALFTERM 0x23
#define ALFSKIP 0x25
#define ALF_START_SEQ 21627 
#define BITE 0x0a
#define RAMEXEC 0x12 /** NEW version !!!! also need to recompile edit_work **/
/** #define RAMEXEC 0x08 **/ /** OLD VERSION using for Jim's tests 04/94 **/
#define BIUDISC 99

/* RTIU Remote Setup and Commands (RTC Tokens) */
#define RUN_PAUSE 20
#define SET_DATA_RETURN 21
#define SET_DATA_COLLECTION 22
#define SET_INIT_COLLECTION 221
#define SET_SEC_COLLECTION 222
#define SET_TEST_COLLECTION 223
#define SEND_MODE 23
#define SEND_CDSCMD 24
#define HANDSHAKE 25
#define SCT_BROADCAST 26
#define SET_SCT 27
#define SET_RTADDR 28
#define RTIU_SET 255

/** Misc RTC and set-up commands (RTC Tokens) **/
#define ARCCTL 30
#define LIMTABLE 31
#define RFESSIM_SELECT 0
#define TSSIM_SELECT 1
#define DCPWR_NOSELECT 110
#define EMERGENCY_OFF 111
#define SEQUENCED_OFF 112
#define DSS_DCPWR 113
#define RFES_DCPWR 114
#define ESS_DCPWR 115
#define GO_IMMEDIATE 11
#define GO 12
#define HALT 13
#define SET_IEBACTIVE 14
#define RESET_IEBACTIVE 15
#define WU1 16

#define MNE_OPEN      57
#define MNE_CLOSE     58

/* RTIU Function Opcodes */
#define OP_CDSXMIT 0xfc0a
#define OP_RUNPAUSE 0xfa14
#define OP_SETDRET 0xfa0b
#define OP_SETDCOL 0xfa07
#define OP_SENDMODE 0xfa10
#define OP_HANDSHAKE 0xfc0b
#define OP_SCTBROAD 0xfc07
#define OP_SETSCT 0xfc05
#define OP_SETRTA 0xfc03

/** Canned Data Collection Schedule values **/
#define BIU_INIT_SCHED_WORDS 4
#define BIU_INIT_SCHED_RATE_MSBS 0x0	/** 64 decimal **/
#define BIU_INIT_SCHED_RATE_LSBS 0x40
#define TLM_INIT_SCHED_WORDS 0x0
#define TLM_INIT_SCHED_RATE_MSBS 0x0	/** 0 decimal **/
#define TLM_INIT_SCHED_RATE_LSBS 0x0
#define HKP_INIT_SCHED_WORDS 85
#define HKP_INIT_SCHED_RATE_MSBS 0	/** 500 decimal **/
#define HKP_INIT_SCHED_RATE_LSBS 0x01f4
#define BIU_SEC_SCHED_WORDS 4
#define BIU_SEC_SCHED_RATE_MSBS 0x0	/** 64 decimal **/
#define BIU_SEC_SCHED_RATE_LSBS 0x40
#define TLM_SEC_SCHED_WORDS 475
#define TLM_SEC_SCHED_RATE_MSBS 0x5	/** 364800 decimal **/
#define TLM_SEC_SCHED_RATE_LSBS 0x9100
#define HKP_SEC_SCHED_WORDS 85
#define HKP_SEC_SCHED_RATE_MSBS 0	/** 500 decimal **/
#define HKP_SEC_SCHED_RATE_LSBS 0x01f4
#define BIU_TEST_SCHED_WORDS 4
#define BIU_TEST_SCHED_RATE_MSBS 0x0	/** 64 decimal **/
#define BIU_TEST_SCHED_RATE_LSBS 0x40
#define TLM_TEST_SCHED_WORDS 475
#define TLM_TEST_SCHED_RATE_MSBS 0x0	/** 30000 decimal **/
#define TLM_TEST_SCHED_RATE_LSBS 0x76C0
#define HKP_TEST_SCHED_WORDS 85
#define HKP_TEST_SCHED_RATE_MSBS 0	/** 500 decimal **/
#define HKP_TEST_SCHED_RATE_LSBS 0x01f4

/** State Check Error Codes **/
#define IEBLOAD_ERROR1 1
#define IEBLOAD_ERROR2 2
#define IEBTRIG_ERROR1 3
#define IEBTRIG_ERROR2 4
#define IEBTRIG_ERROR3 5
#define IEBHALT_ERROR1 6
#define ILX_ERROR1 7

/** more misc (for c2) **/
#define CCSDS_HDR_LEN 6
#define SAB_TYPE 12
#define HKP_TYPE1 7
#define HKP_TYPE2 0x38
#define MRO_TYPE1 0x15
#define MRO_TYPE2 0x2a
#define BIU_TYPE 0xaa   /** not a ccsds code **/
#define MSGQS_OPEN 100
#define CONNECT_TO_RTIU 101
#define CONNECT_TO_SCIPROC 102
#define CONNECT_TO_TELEMPROC 103
#define CLOSE_RTIU_CONNECTION 104
#define CLOSE_TELEMPROC_CONNECTION 105
#define CLOSE_SCIPROC_CONNECTION 106
#define SEND_TO_SCIPROC 107
#define SEND_TO_TELEMPROC 108
#define ARCHIVE_NEXT_N 110
#define ARCHIVE_EVERY_NTH 111
#define ARCHIVE_ALL 112
#define HKP_BYTES 316
#define HKP_SELECT 200
#define ENG_SELECT 201
#define SCI_SELECT 202
#define PN_BYTE_CODE1 0x77 
#define PN_BYTE_CODE2 0x74
#define PN_BYTE_CODE3 0x6b
#define PN_BYTE_CODE4 0x6a
#define SAB_PN_CODE1 0x7774
#define SAB_PN_CODE2 0x6b6a
#define SAB_HEADER_BYTES 80
#define SAB_FOOTER_BYTES 44
#define PN_WORD_SEARCH 0
#define PN_BIT_SEARCH 1
#define MAX_TAIL_TXFR 7000
#define RECV_FROM_RTIU 120
#define RTIU_WRITES_TO_BCE 5555
#define SCI_PORT 3333
#define TLM_PORT 3222

/** Unit Conversions **/
#define RCV_ATTN_CONV 4.0
#define RAD_INTEG_CONV 5.0
#define CSD_CONV 133.333
#define CFS_CONV 1.788
#define CSF_CONV 457.764
#define EFIFO_DELAY_CONV 1.067
#define EFIFO_SCALE_CONV 0.125
#define NCO_FREQ_CONV 1.02997
#define NCO_FREQ_OFFSET 0.098304
#define UPCONV_GAIN_CONV 5.0



/** Here is the socket address structure **/
typedef struct sockaddr_un {
	short sun_family;
	unsigned short sun_path[108];
} SOCKADDR;

/** Here is the message queue stucture **/
typedef struct {
	long mtype;
	unsigned char msgbuf[33000];
} MSGBUF;

#ifndef DESIGN_TIME                                      /* NOTE: The Precompiler variable DESIGN_TIME is
                                                                  defined ONLY for the UIM/X interpreter. 
                                                                  This precompiler directive precludes the
                                                                  UIM/X interpreter from trying to resolve 
                                                                  bit assignments on short int's. It doesn't
                                                                  seem to handle them well. */
/** Here is the RTIU data collection structure **/
typedef struct {
	unsigned  short		num_words	:16;
	unsigned  short		data_rate_msbs	:16;
	unsigned  short		data_rate_lsbs	:16;
} RTIU_COLLECTION;
#endif    /* DESIGN_TIME */

/** Here is the Archive Control data collection structure **/
typedef struct {
	unsigned 	cst_set			: 8;
	unsigned 	cst_spare1		: 8;
	unsigned 	cst_spare2		:16;

	unsigned 	ost_set			: 8;
	unsigned 	ost_spare1		: 8;
	unsigned 	ost_spare2		:16;

	unsigned 	radar_limit_set		: 8;
	unsigned 	radar_limit_spare1	: 8;
	unsigned 	radar_limit_spare2	:16;

	unsigned 	radar_state_set		: 8;
	unsigned 	radar_state_spare1	: 8;
	unsigned 	radar_state_spare2	:16;

	unsigned 	ccsds_rtiu_set		: 8;
	unsigned 	ccsds_rtiu_spare1	: 8;
	unsigned 	ccsds_rtiu_spare2	:16;

	unsigned 	hkp_sab_set		: 8;
	unsigned 	hkp_sab_option		: 8;
	unsigned 	hkp_sab_n		:16;

	unsigned 	mro_set			: 8;
	unsigned 	mro_spare1		: 8;
	unsigned 	mro_spare2		:16;

	unsigned 	sci_sab_set		: 8;
	unsigned 	sci_sab_option		: 8;
	unsigned 	sci_sab_n		:16;

	unsigned 	eng_sab_set		: 8;
	unsigned 	eng_sab_option		: 8;
	unsigned 	eng_sab_n		:16;

	unsigned 	ccsds_set		: 8;
	unsigned 	ccsds_option		: 8;
	unsigned 	ccsds_n			:16;

	unsigned 	hdrftr_set		: 8;
	unsigned 	hdrftr_option		: 8;
	unsigned 	hdrftr_n		:16;

	int	spare_a				: 32;
	int	spare_b				: 32;
	int	spare_c				: 32;
	int	archiver_mode			: 32;
	int	arch_new			: 32;
} ARC;

typedef struct {
	unsigned	ps_addr		: 2,
			ps_mode		: 2,
			ess_addr	: 2,
			ess_mode	: 2,
			fgu_addr	: 1,
			fgu_mode	: 2,
			hpa_addr	: 3,
			hpa_mode	: 2,
			dcg_addr	: 1,
			dcg_mode	: 2,
			cuca_addr	: 1,
			cuca_mode	: 2,
			mr_addr		: 3,
			mr_mode		: 2,
			fee_addr	: 4,
			fee_mode1	: 1,	/** 32 bit boundary **/
			fee_mode2	: 1,
			ant_addr	: 3,
			ant_mode	: 2,
					:26,	/** 32 bit boundary **/
					:23,
					: 1,
			rfessync	: 1,
			ats_addr	: 3,
			ats_mode	: 2,
					: 2;	/** 32 bit boundary **/
}RSIM; /** 96 bits **/

typedef struct {
	unsigned	rec_attn	: 1,
			rec_filter	: 1,
			efifo_delay_1a	: 19,
			efifo_delay_2aa	: 11,	/** 32 bit boundary **/
			efifo_delay_2ab	: 8,
			efifo_delay_b	: 8,
			efifo_scale_1	: 3,
			efifo_scale_2	: 3,
			nco_frequencya	: 10,	/** 32 bit boundary **/
			nco_frequencyb	: 2,
			upconv_noise	: 1,
			upconv_gain	: 4,
			nco_freq_msb	: 1,
			tdas_header_sin : 8,
			tdas_header_fin : 8,
			eds_mode	: 2,	/** 32 bit boundary **/
			ats_mode	: 2, 
			ats_addr	: 3, 
			csync_enable	: 1;	
} TSIM;  /** echo delay simulator -- 96 bits **/

typedef union {
   RSIM rsim;
   TSIM tsim;
} SIM_OVERLAY;  /** 96 bits **/

/** Here is the DSS SE control structure **/
typedef struct {
	SIM_OVERLAY ov;
	unsigned				:32, 
						:16,
						:16,
						:16,
						:16,
						:16,
						:16,
						:16,
						: 8,
						: 1,
						: 1,
						: 1,
						: 1,
						: 1,
				tssim_state	: 1,
				rfessim_state	: 1,
				antsim_state	: 1;

} DSS_SE_CONTROL;  /** 32 bytes **/



/** Here is the test sequence format **/

/***
     These fields must be on a 32-bit boundary, a C restriction.
     Therefore some fields are broken down and require special
     processing.
***/
typedef struct {
	unsigned 	
			sf_time			:16,  /*  1-16  */
						: 6,  /* 17-22  */
			sf_cmd_type		: 2,  /* 23-24  */
			sf_dtn			: 8,  /* 25-32  */

			sf_cmd_number		: 8,  /*  1-8  */
			sf_sci_meas_mode	: 4,  /*  9-12   */
			sf_cal_source		: 4,  /* 13-16  */
			sf_ADC_sample_rate	: 2,  /* 17-18  */
			sf_rcvr_bw		: 2,  /* 19-20  */
			sf_tx_burst_rcv_offset	: 4,  /* 21-24  */
			sf_baq_mode		: 3,  /* 25-27  */
			sf_beam_mask		: 5,  /* 28-32  */

						: 4,  /*  1-4  */
			sf_rcv_attn_12a		: 5,  /*  5-9  */
			sf_rcv_attn_12b		: 5,  /* 10-14  */
			sf_rcv_attn_12c		: 2,  /* 15-16   */
						: 4,  /* 17-20  */
			sf_rcv_attn_3a		: 5,  /* 21-25  */
			sf_rcv_attn_3b		: 5,  /* 26-30  */
			sf_rcv_attn_3c		: 2,  /* 31-32  */

			sf_rcv_attn_45a		: 5,  /* 1-5  */
			sf_rcv_attn_45b		: 5,  /* 6-10  */
			sf_rcv_attn_45c		: 2,  /* 11-12  */
			sf_rad_integ_period	: 4,  /* 13-16  */
			sf_rad_window_cnt	: 8,  /* 17-24 */
			sf_chirp_step_dur	: 8,  /* 25-32 */

						: 4,  /*  1-4 */
			sf_chirp_step_qty	:12,  /*  5-16 */
			sf_chirp_freq_step	:16,  /* 17-32 */


			ff_time			:16,  /*  1-16 */
						: 6,  /* 17-22 */
			ff_cmd_type		: 2,  /* 23-24 */
			ff_cmd_number		: 8,  /* 25-32 */

			ff_bursts_in_cmd	: 8,  /* 1-8 */
			ff_pulses_per_burst	: 8,  /* 9-16 */
						: 4,  /* 17-20 */
			ff_burst_period		:12,  /* 21-32 */

						: 6,  /* 1-6 */
			ff_rcv_window_delay	:10,  /* 7-16 */
						: 6,  /* 17-22 */
			ff_pulse_rep_interval	:10,  /* 23-32 */

			ff_chirp_start_freq	:16,  /* 1-16 */
						:16,  /* 17-32 */


			pr_time			:16,  /*  1-16 */
						: 6,  /* 17-22 */
			pr_cmd_type		: 2,  /* 23-24 */
			pr_power_mode		: 8,  /* 25-32 */

			tc_time			:16,  /*  1-16 */
						: 6,  /* 17-22 */
			tc_ieb_type		: 2,  /* 23-24 */
			tc_mode			: 8,  /* 25-32 */

			tc_parameter_a		:16,  /*  1-16  */
			tc_parameter_b		:16,  /* 17-32 */
			tc_parameter_c		:16,  /*  1-16 */

						:16,  /* 17-32 */

/*** the following fields are used by the command editor s/w only!!  ***/

			sf_num			:16,  /*  1-16 */
			ff_num			:16,  /* 17-32 */
			pr_num			:16,  /* 1-16 */
			tc_num			:16,  /* 17-32 */
/** cmd_mask definition:  4lsbs ==> T&C Pwr Fast Slow **/
			cmd_mask		:16,  /*  1-16 */
			radar_lock		:16;  /* 17-32 */

} RADAR_CONTROL;  /** 32 16-bit words total **/

#ifndef DESIGN_TIME
typedef union {
   struct {
	RADAR_CONTROL radar;
	DSS_SE_CONTROL dss_se;
	RTIU_COLLECTION	rtiu_seq[RTIU_COL_LEN];
	unsigned 	rtiu_index		: 16;
	unsigned 	rtiu_backptr		: 16;
	unsigned 	rtiu_set		: 8;
	unsigned 	num_config		: 16;
	unsigned char	rpb_spare[101];
   } rpb;
unsigned short us_buffer[SEQ_STRUCT_LEN16];
unsigned char buffer[SEQ_STRUCT_LEN];
} SEQUENCE;
#endif   /* DESIGN_TIME */

/** Here is the BCE to RTIU message format **/
typedef struct {
	unsigned
		length			:16,  /*  1-16 */
		fcode			:16,  /* 17-32 */
		status			:16,  /*  1-16 */
		sequence		:16,  /* 17-32 */
		value			:16,  /*  1-16 */
		time_a			:16,  /* 17-32 */
		time_b			:16,  /*  1-16 */
					:13,  /* 17-29 */
		rti_counter		: 3;  /* 30-32 */
	unsigned short   data[128];
} RTIU;

#ifndef DESIGN_TIME
/** Here is the Spacecraft Time data format **/
typedef struct {
	unsigned short	msbs		:16;
	unsigned short	lsbs		:16;
	unsigned short			:13,
			rti_counter	: 3;
} SCT;
#endif   /* DESIGN_TIME */

/** Here is the Radar State shared memory format **/
typedef struct {
	unsigned 	dcreg_pwr	:16;  /** 11 lsbs are valid **/
	unsigned 	fsw_mode	:8;
	unsigned 	ieb_halt	:8;
	unsigned 	ieb_trig	:8;
	unsigned 	ieb_load	:8;
	unsigned 	ilx		:8;
} RADAR_STATE;

/** Here is the SAB Header (Radar State) shared memory data format **/
typedef struct {
	unsigned 	sync		:32; 
	unsigned 	sct		:32;
	unsigned 	pickup_rate	:4,
			bst		:12;
	unsigned char	sf[20];
	unsigned char	ff[14];
	unsigned char	tc[10];
	unsigned char	pr[4];
	unsigned 	sc_cmd		:16;
	unsigned 	sc_msg		:16;
	unsigned 	tail_len	:16;
	unsigned 	tail_id		:16;
	unsigned 	sab_ctr		:16;
	unsigned 	sab_len		:16;
	unsigned 	fsw_err		:8;
	unsigned 	fsw_err_code	:8;
	unsigned 	ctr_bii		:8;
	unsigned 	ctu_rx		:8;
	unsigned 	ctu_dcreg	:11,
			ctu_beam	:5;
	unsigned	instr_readback	:17,
			spare1		:1,
			sab_type	:1,
			spare2		:1,
			run_mode	:1,
			ieb_halt	:1,
			ieb_trig	:1,
			ieb_load	:1,
			ilx		:1,
			last_ieb	:1,
			ctu_irc		:1,
			spare3		:5;
} SAB_HEADER;


/** Here is the Housekeeping Frame format with RTC header info for 
    the Telemetry Processor TCP socket and Archive Processor queue **/
typedef struct {
	unsigned long	type;
	unsigned long	length;		/** 228 bytes **/
	SAB_HEADER	sab_header;
	unsigned short	eng_data[96];
	unsigned short	eng_footer[22];
} HOUSEKEEP_FRAME;


/** Here is the Engineering Hdr/Data/Ftr/Tail Frame format with RTC header info 
    for the Telemetry Processor TCP socket and Archive Processor queue **/
typedef struct {
   unsigned long	type;
   unsigned long	length;			/* min 228 bytes, max 268 bytes;	*/
						/* depending on length of tail data	*/

   unsigned short	sab_header  [40];	/* See Table7-1 in the HLD		*/
   unsigned short	eng_data    [96];	/* See Table7-3 in the HLD		*/
   unsigned short	eng_footer  [22];	/* See Table7-4 in the HLD		*/
   unsigned short	tail_data   [40];	/* Format is TBD			*/
} ENGINEER_FRAME;

/** Here is the Science Hdr/Ftr/Tail Frame format with RTC header info 
	for the Telemetry Processor TCP socket **/
typedef struct {
	unsigned long	type;
	unsigned long	length;		/** variable **/
	SAB_HEADER	sab_header;
	unsigned short	sci_footer[22];
	unsigned short	tail_data[32768];
} HFT_FRAME;

/** Here is the Science/Engineering Hdr/Ftr Frame format with RTC header info 
	for the Archive Processor HDR/FTR Archive Buffer **/
typedef struct {
	unsigned long	type;
	unsigned long	length;		/** variable **/
	SAB_HEADER	sab_header;
	unsigned short	eng_footer[22];
} HF_FRAME;

/** Here is the continuation of science Tail data (CTU waveform)
	for the Telemetry Processor TCP socket
	   	   This Structure is Obsolete. It will be removed **/
typedef struct {
	unsigned long	type;
	unsigned long	length;	
	unsigned short	xfer_number;
	unsigned short	tail_data[4086];
} TAIL_FRAME;



/** Here is the SAB Frame format for the Science & Eng Processor 
	TCP socket and Archive Processor SAB Archive Buffer **/
typedef struct {
	unsigned long	type;
	unsigned long	length;	
	SAB_HEADER	sab_header;
	unsigned short	sab_dft[32790]; /** 16KW + 22W + 16KW **/
} SAB_FRAME;

/** Here is the Eng Ftr (Rad data) Frame format for the Science Processor 
	TCP socket **/
typedef struct {
	unsigned short	eng_footer[22];
} ENGFTR_FRAME;


/** Here is the Telemetry Processor BIU Discrete Status Frame format **/
typedef struct {
	unsigned long	type;
	unsigned long	length;
	unsigned short	BIU_status[4];  /** see 3-271 pg 74 **/
} BIU_STATUS_FRAME;


/** note: CCSDS_HEADER forces word boundaries to avoid fills by compiler **/
typedef struct {
	unsigned 	ver 	: 3,
			type	: 1,
			hflag	: 1,
			id	: 5,
			pkttyp	: 6;
	unsigned 	segflg	: 2,
			seq	:14;  /** word boundary **/
	unsigned 	len	:16;
	unsigned 	err	:8;
	unsigned 	sct1	:8;   /** word boundary **/
	unsigned 	sct2	:16;
	unsigned 	sct3	:8;
	unsigned 	timer	:8;   /** word boundary **/
} CCSDS_HEADER;


/** Here is the CCSDS packet format, including RTIU-to-BCE header **/
typedef struct {
	unsigned 
                        length                  :16,  /*  1-16 */
                        fcode                   :16,  /* 17-32 */
                        status                  :16,  /*  1-16 */
                        rt_addr                 : 8,  /* 17-24 */
                        rt_subaddr              : 8,  /* 25-32 */
                        value                   :16,  /*  1-16 */
                        time_a                  :16,  /* 17-32 */
                        time_b                  :16,  /*  1-16 */
                                                :13,  /* 17-29 */
                        rti_counter             : 3;  /* 30-32 */
        unsigned short  status_1553B[18];

	CCSDS_HEADER	ccsds_header;
	unsigned char	data[1078];  /** may need to change later **/
} TELEMETRY; 


/** Here is the memload PSDL buffer space **/
typedef union {
   unsigned char cval[64];
   unsigned short sval[32];
} PSDL;

/** Here is the memload ALF data format **/
typedef struct {
	unsigned 	op_code		:8;
	unsigned 	seq_num		:8;
	unsigned 	seg_num		:8;
	unsigned 	tot_num		:8;
	unsigned short	data[16];
	unsigned 	even_chksum	:16;
	unsigned 	odd_chksum	:16;
} ALF;

	/* Slow Field State Error */
typedef struct {
   unsigned long	type;
   unsigned long	length;			/* 48 bytes */
						/* Section 12 of the HLD defines the Slow Field */

   unsigned short	mask      [10];		/* Bits in 'mask' that are non-zero indicate	*/
						/* which fields are currently in error		*/

   unsigned short	expected  [10];		/* Fields in 'expected' are only meaningful 	*/
						/* where the cooresponding field is non-zero	*/
						/* in 'mask' above. 'expected' shows the 	*/
						/* expected value for the fields that are 	*/
						/* currently in error				*/
} SLOW_ERR_STRUCT;


	/* Fast Field State Error */
typedef struct {
   unsigned long	type;
   unsigned long	length;			/* 36 bytes */
						/* Section 13 of the HLD defines the Slow Field */

   unsigned short	mask      [7];		/* Bits in 'mask' that are non-zero indicate	*/
						/* which fields are currently in error		*/

   unsigned short	expected  [7];		/* Fields in 'expected' are only meaningful 	*/
						/* where the cooresponding field is non-zero	*/
						/* in 'mask' above. 'expected' shows the 	*/
						/* expected value for the fields that are 	*/
						/* currently in error				*/
} FAST_ERR_STRUCT;


	/* T & C Instruction State Error */
typedef struct {
   unsigned long	type;
   unsigned long	length;			/* 28 bytes */
						/* Section 11 of the HLD defines the Slow Field */

   unsigned short	mask      [5];		/* Bits in 'mask' that are non-zero indicate	*/
						/* which fields are currently in error		*/

   unsigned short	expected  [5];		/* Fields in 'expected' are only meaningful 	*/
						/* where the cooresponding field is non-zero	*/
						/* in 'mask' above. 'expected' shows the 	*/
						/* expected value for the fields that are 	*/
						/* currently in error				*/
} TNCI_ERR_STRUCT;


	/* Power Instruction State Error */
typedef struct {
   unsigned long	type;
   unsigned long	length;			/* 16 bytes */
						/* Section 10 of the HLD defines the Slow Field */

   unsigned short	mask      [2];		/* Bits in 'mask' that are non-zero indicate	*/
						/* which fields are currently in error		*/

   unsigned short	expected  [2];		/* Fields in 'expected' are only meaningful 	*/
						/* where the cooresponding field is non-zero	*/
						/* in 'mask' above. 'expected' shows the 	*/
						/* expected value for the fields that are 	*/
						/* currently in error				*/
} PWRI_ERR_STRUCT;


/* Test Status GUI specific typedefs and definitions */

/* error message */
typedef struct {
        int     err_time;               /* 32 bit seconds, displayed as DDDDD/HH:MM:SS */
        char    err_text[256];          /* message ascii text */
} ERR_MESSAGE_DISPLAY;
   
/* time display message */
typedef struct {
        int     sc_time;		/* seconds */
        int     irig_time;		/* seconds */
        int     sc_irig_delta;		/* seconds */
        int     sc_irig_cmp_error;	/* 0 for good, 1 for error */
} TIME_DISPLAY;
   
/* command start message, sent at every command start */
typedef struct {
	int	sc_time;		/* seconds */
	int	config_number;		/* configuration ID */
	int	cmd_start_error;	/* 0=no error, 1=missed a cmd start,
					   2=received unexpected cmd start */
	int	rsync_error;		/* 0=no error, 1=missed a rfes sync,
					   2=received unexpected rfes sync */
} CST_DISPLAY;
   
/* dss status display message, sent every second */
typedef struct {
	int	sc_time;		/* seconds */
	float	dss_voltage;
	float	dss_margin_error;
	float	ess_voltage;
	float	ess_margin_error;
	float	rfes_voltage;
	float	rfes_margin_error;
} DSS_DISPLAY;

/* tsb error message, sent when in error */
typedef struct {
	int		sc_time;	/* seconds */
	int		config_number;	/* configuration ID */
	unsigned short	tsb[6];		/* 96bit register contents */
	unsigned short	gse[6];		/* expected values */
} TSB_DISPLAY;


/* telemetry status display message */ 
typedef struct {
        unsigned	length		:16,
			fcode		:16,
			status		:16,
			sequence	:16,
			value		:16,
			time_a		:16,
			time_b		:16,
					:13,
			rti_counter	: 3;
	unsigned short	num_packets; 
} RTIU_HDR_DISPLAY;
   
/* ccsds header status display message */ 
typedef struct {
        CCSDS_HEADER    ccsds_header;   /* ccsds header record */
        unsigned        short   pass_indicator; /* pass fail indicator, 
                                                   0=all valid packets received,
                                                   1=no packets received, 
                                                   2=got some error packets 
						   3=exceeded error threshold */
        unsigned int  num_packets;    	/* number of CCSDS packets received 
				           int the last sec */
        unsigned int  percent_failed_sec;  /* number of failed packets 
					      in the last sec */
} CCSDS_HDR_DISPLAY;
   
/* Input Message Types for Test Status GUI */ 
#define RTIU_HDR_DISP_MSG       13L
#define TIME_DISP_MSG           14L
#define CCSDS_HDR_DISP_MSG      15L
#define TSB_DISP_MSG            16L
#define CSTRT_DISP_MSG 		17L
#define DSS_DISP_MSG            18L
#define ARCH_CNTL_MSG		19L
#define CMP_ENG_MSG		20L
#define CMP_STATE_MSG		21L
#define TS_MSG			22L   


#endif /* __DEFINITIONS_H */



