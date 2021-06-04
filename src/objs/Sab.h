#ifndef Sab_H
#define Sab_H

#include <vector>

//----------------------------
// Forward declarations
//---------------------------

class Sab;
class IncompleteSab;
class OverfullSab;

#include "Error.h"
#include "Io.h"
#include "Units.h"
#include "Array.h"
#include "Ieb.h"
#include "CAGTable.h"

class Sab
  {
  public:

  //--------------------
  // Typedef's and enums
  //--------------------

  typedef std::vector<unsigned char> cdata;
  typedef std::vector<unsigned char>::const_iterator const_CDI;
  
  enum record_typeE {archive_header = 81, hkp = 200, mro = 201,
    sci = 202, eng = 203, ccsds = 210}; 
  // number of blocks for threshold values
  // unit_file_size sets the number of characters written to a file to
  // store units.
  static const int unit_file_size = 10;

  // Constructors
  // See copy constructor in Sab.cpp
  Sab(FileMgr::orderE mt);

  // Predicates
  bool complete() const;
  bool partial() const;
  bool empty() const;
  bool decoded() const;
  bool baq_decoded() const;
 
  
  //special method to correct sclk by subtracting 1 from Sab's original record
  void subtract1secfromSCLK();

 

  // I/O
  friend std::ostream& operator<<(std::ostream& s, const cdata& c);
  friend std::ostream& operator<<(std::ostream& s,
    const std::vector<unsigned short>& sint);

  // Data handling
  
  void decode() ;
  bool Valid() const;
  void decodeSclk() throw(ErrorMessage);
  void decodeSABcounter() throw(ErrorMessage) ;
 
  void baq_decode() ;
  void addData(const cdata& record) throw(ErrorMessage);
  void addPartialRecord(const cdata& record)throw(ErrorMessage);
  void xfer1_eng_temp_data();
  void xfer2_sub_comm_data();
  int getSabEntry(unsigned int tablenum, FileMgr::orderE machine_type,
    const cdata& data) const throw(ErrorMessage);
  double tadc_conversion(const unsigned short  i ,
			 const unsigned int adc_sub_comm);
  void extract_baq_raw_data(vector<unsigned short>&baq, cdata& rdata) const;
			   
  void get_raw_data(vector<unsigned char>& data);

  //---------------------------------------
  //debugging purpose
  //plotting active data on the screen   
  //----------------------------------------
  void show_data_Onscreen();
  void generatePlot();
  void get_Sab_data_statistics(double& standard_deviation,double& average,Uvar& receiver_bandwidth,double& attenuation);
  void constructActiveRecord(const Time& t, 
			     const unsigned int& burt_number,
			     const unsigned int& beam_number,
			     Ieb& ieb, 
			     const Array1D<float>& radar_data);

  // Other methods
  void clear();

  //{get} Ieb 
  Ieb getIeb() const;

  //{get} error location
  int getErrorInfo();

  //------------------------
  // SAB public data fields
  //------------------------

  // header
  //variables listed in the same order described in the Blue Book
  //do not attemp to change the order 

  unsigned int sync;
  unsigned int sclk;
  unsigned int scpr;
  Uvar brst;
  unsigned int int_brst;//integer value

  Uvar header_tfi;
  unsigned int header_tnc,header_typ,header_tca,header_tcb,header_tcc;

  unsigned int pwri, vicc,vimc,tail_len,tail_id,sab_counter,sab_len;
  unsigned int fswm,fswc,ctbc,ctrx;

  unsigned int ctps,ctbe,ctps_ctbe,header_end;

  //slow field variables
  Uvar slow_tfi;
  unsigned int dtn;
  unsigned int slow_typ;
  unsigned int csr,r_mode,slow_instruction_number,bem,baq_mode;
  Uvar tro;
  int tro_in_units_of_pri;

  Uvar rc_bw,adc;

  double at1,at3,at4;
  unsigned int at1_db,at3_db,at4_db;
  unsigned int at1_db_mask, at3_db_mask,at4_db_mask;
  Uvar rip;

  Uvar csd;
  unsigned int rad;

  unsigned int csq;
  Uvar chirp_length;

  Uvar slow_cfs;

  //fast field variables
  Uvar fast_tfi;
  unsigned int fin,fast_typ,pul,bii;
  Uvar bpd;

  Uvar pri,rwd,fast_csf;

  //sci data field
  double  scatt_dc_offset;
  fdata  decoded_data;
  unsigned int Ndecoded;
  float rms_radar_data;
 
  //footer1 variables    
  LUvec  baq_threshold;
  LUvec  eng_tmp;
  LUmat  sub_comm_data;
  
  unsigned int  cnt_rl,cnt_radio,cnt_nd;
  unsigned int   eout;
  int subr;
  unsigned int space_craft_time;
  Uvar hip,cip;

  unsigned int  iebtth, iebttl;
  unsigned int  bgcalls,delvmn,delvda,delvyr;
     
  //sub_comm_data variables 
  Uvar fwdtmp,be1tmp,be2tmp,be3tmp,be4tmp,be5tmp;
  Uvar diptmp,rlotmp,tadcal1,nsdtmp,lnatmp,evdtmp;
  Uvar mratmp, mruttm,dcgttm,cucttm,twttmp,epctmp;
  Uvar tw1ttm,ep1ttm,p_stmp,p_sttm,fguttm,tadcal4;
  Uvar esstmp,wgb1t1,wgb3t1,wgb3t2,wgb3t3,wgb5t1;
  Uvar pcutmp,adctmp,tadcal2,ecltmp,cputmp,memtmp;
  Uvar sadctmp;

  Uvar tadcal3,frwdpw;
  unsigned int  dcgmon;
  double  lpltlm_db;

  Uvar nsdcur,hpapsm,catcur,p_smon,svlsta,usotmp;
  Uvar cpbnkv, essvlt,tadcal5,pcu5v_72,pcu5i_73,pcu5v_74;
  Uvar pcuti_75,pcu15v_76,pcu15i_77,pcu15v_78,pcu15i_79,pcu12v_80;
 
  Uvar pcu12i_81,pcucur,pllmon,ctu5i;
  unsigned int  tadcal6;
  Uvar pcu9v_86;
  
  Uvar pcu9i_87,pcu9v_88,pcu9i_89;
  Uvar  tadcl7;
  Uvar shpttm;

  private:

  //---------------------
  // Debug access methods
  //---------------------

  void showHeader() const;
  void showFooter() const;
  void showTail() const;
  void showData() const;

  //------------------------
  // Internal representation
  //------------------------
  // Tracking variables
  unsigned int tail_length_;
  unsigned int Nrdat_;
  FileMgr::orderE machine_type_;
  bool decoded_;
  bool baq_decoded_;
  bool is_sab_valid_;

  // Data variables
  cdata header_;
  cdata footer_;
  cdata tail_;
  cdata rdat_;

  //ieb 
  Ieb ieb_;

  //container for RL and ND integration time
  Uvar hip_hold_,cip_hold_;

  //error location
  int error_code_;
 
  //-------------
  //calibration table  variables
  //---------------
  CAGTable cag_table_;


  //-------------------------------------------------------------------
  // ByteIndexTable is a sub-class which aids the static initialization
  // of the header and footer byte index tables.
  //-------------------------------------------------------------------

  class ByteIndexTable : public cdata
    {
    public:
    // Special constructor
    ByteIndexTable(const_CDI p1, const_CDI p2);
    };

 
  //--------------------
  // Constant class data
  //--------------------

  const static unsigned int Nheader_bytes_ = 80;
  const static unsigned int Nfooter_bytes_ = 44;
  const static unsigned int Nb_ = 24;
  const static unsigned int Neng_tmp_=12;
  const static unsigned int Ndat_max_ = 32768;//32*1024
  const static int Nsub_row_ = 15;
  const static int Nsub_col_ = 6;
 
  //-----------------------------------------------------------------
  // Static data Nb_: number of blocks for threshold values
  // Static data tables used to initialize the static _byte_pos, and _bytes
  // Static data Ndmax_: max. number of data (16K x 16)
  // Static Neng_tmp_ : 12 engineering temperature values stored in first 192 bits in footer
  // tables.
  //------------------------------------------------------------------------

  const static unsigned char hb_[];
  const static unsigned char fb_[];
  const static unsigned char hnum_[];
  const static unsigned char fnum_[];

  //--------------------------------------------------------------------
  // The byte_pos tables index the byte position of an entry in the SAB
  // header table (BB, 7-2).  The _bytes tables index the size of an entry
  // in bytes.
  //-------------------------------------------------------------------

  const static ByteIndexTable header_item_byte_pos_;
  const static cdata header_item_bytes_;
  const static ByteIndexTable footer_item_byte_pos_;
  const static cdata footer_item_bytes_;
  };

//-------------------
// Exception classes
//-------------------

class IncompleteSab : public ErrorMessage
  {
  public:
  IncompleteSab() throw() { msg = "Incomplete SAB"; }
  };
  
class OverfullSab : public ErrorMessage
  {
  public:
  OverfullSab() throw() { msg = "Overfull SAB"; }
  };
  
#endif

