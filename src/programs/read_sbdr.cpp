//---------------------------------------------
// Filename: read_lbdr.cpp
//
// Date Written:  Feb. 18, 2004
//
// Programmer: Young Gim at JPL
//
// Description: This program will take inputs of
//               LBDR (long burst-ordered data record) file name
//               and sab_counter(or burst number), and display 
//               information of each LBDR parameter of the sab counter.
//
//  Usage: read_lbdr lbdr_file_name sab_counter
//
//
//---------------------------------------------

#define SYNC_VAL 0x77746B6A
#define UTC_STRING_SIZE 21
#define UTC_STRING_PAD 3
#define TIME_LENGTH_IN_FILE  (UTC_STRING_SIZE+UTC_STRING_PAD+16)
#define Nradar_data_max (32*1024)
#define TARGET_NAME_STRING_SIZE 16
#define TARGET_FRAME_NAME_STRING_SIZE 24
#define SBDR_RECORD_LENGTH 1272
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::cerr;
using std::swap;


void readFloat(float& x, const bool& byte_flip, FILE* fp);
void readDouble(double& x, const bool& byte_flip, FILE* fp);
void readUnsignedInt(unsigned int& i, const bool& byte_flip,FILE* fp);

int main(int argc, char* argv[])
{

  
  //----------------------------
  // Parse the command line
  //------------------------
   
  const char* command = argv[0];
  if (argc != 3)
    {
      cerr << "Usage: " << command << " file_name sab_number " << endl;
      exit(-1);
    }

  int clidx = 1;
  char* lbdr_filename=argv[clidx++];
  char* sab_input_record = argv[clidx++];
  int sab_number = (int) atoi(sab_input_record);
  if(sab_number<0) {
    cerr<<"Input sab number should be a positive integer number"<<endl;
    exit(-1);
  }
  unsigned int input_sab_number = (unsigned int)sab_number;

  //---------------
  //header size and byte_flip flag  declarations
  // Information from SIS
  //-------------
  bool byte_flip= false;
  int fread_return;
  unsigned int Nrecord;
  unsigned int from_top_to_sab_counter=3*sizeof(float);
  from_top_to_sab_counter+= 13*sizeof(unsigned int);
  unsigned int from_sab_counter_to_end;
  //-----------------
  //variables to be read
  //--------------------

  //--------------
  //SAB data
  //-----------------
  unsigned int  sync;
  unsigned int  sclk;  // actual telemetry value
  unsigned int  record_id; 
  float  scpr;
  float  brst;
  float  header_tfi;
  unsigned int  header_tnc,header_typ,header_tca,header_tcb,header_tcc;
  unsigned int  pwri, vicc,vimc,tail_len,tail_id,sab_counter,sab_len;
  unsigned int  fswm,fswc,ctbc,ctrx,ctps,ctbe,ctps_ctbe,header_end;
  float slow_tfi;
  unsigned int dtn;
  unsigned int slow_typ;
  unsigned int csr,r_mode,slow_instruction_number,bem,baq_mode;
  float tro;
  float rc_bw,adc;
  float  at1,at3,at4;
  unsigned int at1_each,at3_each,at4_each;
  float rip,csd;
  unsigned int rad,csq;
  float chirp_length,slow_cfs;
  float fast_tfi;
  unsigned int fin,fast_typ,pul,bii;
  float bpd,pri,rwd,fast_csf;
  unsigned int  iebtth, iebttl;
  unsigned int  bgcalls,delvmn,delvda,delvyr;
  unsigned int  cnt_rl,cnt_radio,cnt_nd;
  unsigned int   eout,subr,space_craft_time;
  float hip,cip;
  float fwdtmp,be1tmp,be2tmp,be3tmp,be4tmp,be5tmp;
  float diptmp,rlotmp,tadcal1,nsdtmp,lnatmp,evdtmp;
  float mratmp, mruttm,dcgttm,cucttm,twttmp,epctmp;
  float tw1ttm,ep1ttm,p_stmp,p_sttm,fguttm,tadcal4;
  float esstmp,wgb1t1,wgb3t1,wgb3t2,wgb3t3,wgb5t1;
  float pcutmp,adctmp,tadcal2,ecltmp,cputmp,memtmp;
  float sadctmp;
  float tadcal3,frwdpw;
  float dcgmon;
  float lpltlm_db;
  float nsdcur,hpapsm,catcur,p_smon,svlsta,usotmp;
  float cpbnkv, essvlt,tadcal5,pcu5v_72,pcu5i_73,pcu5v_74;
  float pcuti_75,pcu15v_76,pcu15i_77,pcu15v_78,pcu15i_79,pcu12v_80;
  float pcu12i_81,pcucur,pllmon,ctu5i;
  float tadcal6;
  float pcu9v_86;
  float pcu9i_87,pcu9v_88,pcu9i_89, tadcl7,shpttm;
  unsigned int num_bursts_in_flight;
 
  //---------------
  //Geometry variable declarations
  // quality flag
  // for burst
  // 0 = good
  // 1 = ckernel error
  // 255 = misc. geometry error
  //----------------
  double t,t_sclk;
  char time_utc[UTC_STRING_SIZE+UTC_STRING_PAD];
  float T_scwg;                  
  float T_feed;
  float T_hga;
  double transmit_time_offset;
  double time_from_closest_approach;
  double time_from_epoch;
  unsigned int quality_flag; 
  char  target_name[TARGET_NAME_STRING_SIZE];
  char tbf_frame_name[TARGET_FRAME_NAME_STRING_SIZE];

  double pole_right_ascension;
  double pole_declination;
  double target_rotation_rate;
  double target_rotation_angle;



  unsigned int  beam_number;   
  double  sc_position_J2000_x, sc_position_J2000_y, sc_position_J2000_z;     
  double sc_position_target_x,sc_position_target_y,sc_position_target_z;      
  double  sc_velocity_J2000_x, sc_velocity_J2000_y, sc_velocity_J2000_z;     
  double sc_velocity_target_x,sc_velocity_target_y,sc_velocity_target_z;      
  double  sc_X_J2000_x, sc_X_J2000_y, sc_X_J2000_z ;
  double  sc_Y_J2000_x, sc_Y_J2000_y, sc_Y_J2000_z ;
  double  sc_Z_J2000_x, sc_Z_J2000_y, sc_Z_J2000_z ;
  double   sc_X_target_x,  sc_X_target_y,  sc_X_target_z;
  double   sc_Y_target_x,  sc_Y_target_y,  sc_Y_target_z;
  double   sc_Z_target_x,  sc_Z_target_y,  sc_Z_target_z;
  double  rot_vel_J2000_x,  rot_vel_J2000_y,  rot_vel_J2000_z;
  double  rot_vel_target_x, rot_vel_target_y, rot_vel_target_z;
  float norm_cnt_rl;
  float  norm_cnt_nd;
  float  norm_cnt_radio;
 
  //---------------------
  //SBDR Science Field data declarations
  //------------------------------
  unsigned int science_qual_flag;
  float system_gain;
  float antenna_brightness_temp;
  float system_noise_temp;
  float abt_std;
  float pass_geom_time_offset;
  float pass_pol_angle;
  float pass_emission_angle;
  float pass_azimuth_angle;
  float pass_centroid_lon;
  float pass_centroid_lat;
  float pass_major_width;
  float pass_minor_width;
  float pass_ellipse_pt1_lon;
  float pass_ellipse_pt2_lon;
  float pass_ellipse_pt3_lon;
  float pass_ellipse_pt4_lon;
  float pass_ellipse_pt1_lat;
  float pass_ellipse_pt2_lat;
  float pass_ellipse_pt3_lat;
  float pass_ellipse_pt4_lat;
  unsigned int num_pulses_received;
  float total_echo_energy;
  float noise_echo_energy;
  float x_factor;
  float sigma0_uncorrected;
  float sigma0_corrected;
  float sigma0_uncorrected_std;
  float altitude_means;
  float altitude_means_std;
  
  float act_geom_time_offset;
  float act_pol_angle;
  float act_incidence_angle;
  float act_azimuth_angle;
  float act_centroid_lon;
  float act_centroid_lat;
  float act_major_width;
  float act_minor_width;
  float act_ellipse_pt1_lon;
  float act_ellipse_pt2_lon;
  float act_ellipse_pt3_lon;
  float act_ellipse_pt4_lon;
  float act_ellipse_pt1_lat;
  float act_ellipse_pt2_lat;
  float act_ellipse_pt3_lat;
  float act_ellipse_pt4_lat;
  float altimeter_profile_range_start;
  float altimeter_profile_range_step;
  unsigned int altimeter_profile_length;
  float  sar_azimuth_res;
  float  sar_range_res;
  float  sar_centroid_bidr_lon;
  float sar_centroid_bidr_lat;
 
  //-----------------
  //Active data variable and array declarations
  //-------------------

  unsigned int Nradar_data;
  float rms_radar_data;
  
  //------------------
  //Open input file
  //-------------------
  FILE  *fp;
  int f_position=0;
  if((fp=fopen(lbdr_filename,"rb"))==NULL){
    cerr<<"Can not open file "<<endl;
    exit(-1);
  }
  
  //-----------------------
  //set file position at the beginning
  //------------------------
  if(fseek(fp,0,0)!=0){
    cerr<<"Can't set file pointer at the beginning of the file"<<endl;
    exit(-1);
  }
 
  //--------------------------------------
  //read total number of active record
  //
  //total number of active record is stored in PDS label
  // starting and ending bytes: 192-196
  //--------------------------------------
  if(fseek(fp,190,0) !=0){
    cerr<<"Can't set file pointer at 192 bytes from beginning of file"<<endl;
    exit(-1);
  }
  char num_of_record_char[5];
  fread_return=fread(num_of_record_char,5,1,fp);
  if(fread_return!=1){
    cerr<<"fail to read total number of active record"<<endl;
    exit(-1);
  }
  for(unsigned int k=0;k<5;++k) cout<<num_of_record_char[k];
  cout<<endl;
  Nrecord= (int) atoi(num_of_record_char); 
  cout<<"total number of active records "<<Nrecord<<endl;
  

  //-----------
  //move file pointer to the beginning of record
  //----------
  if(fseek(fp,2*SBDR_RECORD_LENGTH,0)!=0){
    cerr<<"fail to move file pointer to the beginning of record"<<endl;
    exit(-1);
  }


  //----------------------------------------------------------------
  //Check Endedness
  //-----------------------------------------------------------
  unsigned int sync_no_flip, sync_flip;
  readUnsignedInt(sync_no_flip,false,fp);
  f_position =ftell(fp);
  f_position -= sizeof(unsigned int);
  if(fseek(fp,f_position,0)!=0){
    cerr<<"Can not set file position right before sync "<<endl;
    exit(-1);
  }  
  readUnsignedInt(sync_flip,true,fp);
  if(SYNC_VAL==sync_no_flip) byte_flip= false;
  else if(SYNC_VAL==sync_flip) byte_flip= true;
  else{
    cerr<<"Invalid sync code : can not determine byte-swapping order"<<endl;
    exit(-1);
  }
  
  //--------------  
  //rewind
  //---------------
  if(fseek(fp,2*SBDR_RECORD_LENGTH,0)!=0){
    cerr<<"fail to move file pointer to the beginning of record"<<endl;
    cerr<<"after reading sync code "<<endl;
    exit(-1);
  }

  
  //-------------------------
  // Read header and determine number of records
  // and record length
  //------------------------
  from_sab_counter_to_end= SBDR_RECORD_LENGTH - from_top_to_sab_counter -sizeof(unsigned int);
  cout<<"Num Records  "<<Nrecord<<endl;
  cout<<"Record length "<<SBDR_RECORD_LENGTH<<endl;
  cout<<"from top to sab counter "<< from_top_to_sab_counter<<endl;
  cout<<"size of sab counter "<<sizeof(sab_counter)<<endl;
  cout<<"after sab counter to end "<< from_sab_counter_to_end<<endl;
  

  //---------------
  //Loop over
  //----------------
  unsigned int sab_counter_start = 0;
  unsigned int sab_counter_end = 0;
  bool found_sab_counter= false;
  for(unsigned int i=0;i<Nrecord;++i){
    //--------------------------------
    //read SAB counter only
    // length from top of the record to sab counter: from_top_to_sab_counter
    // length from sab counter to the end: from_sab_counter_to_end
    // record legnth: SBDR_RECORD_LENGTH
    //--------------------------------
    f_position=ftell(fp);
    f_position = f_position+(int) from_top_to_sab_counter;
    if(fseek(fp,f_position,0)!=0){
      cerr<<"Error setting current file position "<<endl;
      exit(-1);
    }
    
    readUnsignedInt(sab_counter,byte_flip,fp);
    
    if (i == 0)
      {
      sab_counter_start = sab_counter;
      }
    else if (i == Nrecord - 1)
      {
      sab_counter_end = sab_counter;
      }
      cout<<"read in sab counter "<<sab_counter<<endl;
    if(sab_counter>input_sab_number) {
      //cout<<"record sab number "<< sab_counter<<endl;
      break;//no need to see other record
    }
    else if(sab_counter< input_sab_number){
      f_position=ftell(fp);
      f_position = f_position +(int) from_sab_counter_to_end;
      if(fseek(fp,f_position,0)!=0){
	cerr<<"Error setting current file position "<<endl;
	exit(-1);
      }
      continue;
    }
    else{}

    found_sab_counter=true;
    cout<<"Display sab information for sab counter number "<<sab_counter<<endl;
    //--------------------
    //If sab counter == input sab number
    // read record after move the pointer to the
    // begining of the record
    //---------------------
    f_position= ftell(fp);
    f_position -= (int) from_top_to_sab_counter + sizeof(sab_counter);
    if(fseek(fp,f_position,0)!=0){
      cerr<<"Error setting current file position "<<endl;
      exit(-1);
    }    
    int start_position=f_position;
    //sync code  
    readUnsignedInt(sync,byte_flip, fp);
    cout<<"sync "<< sync<<endl;
    
    //sclk
    readUnsignedInt(sclk,byte_flip,fp);
    cout<<"sclk "<< sclk<<endl;
    
    //record id
    readUnsignedInt(record_id,byte_flip,fp);
    cout<<"reacord id "<< record_id<<endl;
    
    //scpr
    readFloat(scpr,byte_flip,fp);
    cout<<"scpr "<< scpr<<" Hz"<<endl;
    
    //brst
    readFloat(brst,byte_flip,fp);
    cout<<"brst  "<<brst<<" s"<<endl;
    
    //header_tfi
    readFloat(header_tfi,byte_flip,fp);
    cout<<"header TFI "<< header_tfi<<" s"<<endl;
    
    //header tnc
    readUnsignedInt(header_tnc,byte_flip,fp);
    cout<<"header tnc "<< header_tnc<<endl;
    
    //header typ
    readUnsignedInt(header_typ,byte_flip,fp);
    cout<<"header typ "<<header_typ<<endl;
    
    //header tcn
    readUnsignedInt(header_tca,byte_flip,fp);
    cout<<"header tca "<< header_tca<<endl;
    
    //header tcb
    readUnsignedInt(header_tcb,byte_flip,fp);
    cout<<"header tcb "<<header_tcb<<endl;
    
    //header tcc
    readUnsignedInt(header_tcc,byte_flip,fp);
    cout<<"header tcc "<< header_tcc<<endl;
    
    //pwri
    readUnsignedInt(pwri,byte_flip,fp);
    cout<<"pwri "<<pwri<<endl;
    
    //vicc
    readUnsignedInt(vicc,byte_flip,fp);
    cout<<"vicc "<<vicc<<endl;
    
    //vimc
    readUnsignedInt(vimc,byte_flip,fp);
    cout<<"vimc "<<vimc<<endl;
    
    //tail length
    readUnsignedInt(tail_len,byte_flip,fp);
    cout<<"tail length "<<tail_len<<endl;
    
    //tail_id
    readUnsignedInt(tail_id,byte_flip,fp);
    cout<<"tail id "<<tail_id<<endl;
    
    //sab counter
    readUnsignedInt(sab_counter,byte_flip,fp);
    cout<<"sab counter "<< sab_counter<<endl;
    
    //sab length
    readUnsignedInt(sab_len,byte_flip,fp);
    cout<<"sab length "<<sab_len<<endl;
    
    //fswm
    readUnsignedInt(fswm,byte_flip,fp);
    cout<<"fswm "<< fswm<<endl;
    
    //fswc
    readUnsignedInt(fswc,byte_flip,fp);
    cout<<"fswc "<<fswc<<endl;
    
    //ctbc
    readUnsignedInt(ctbc,byte_flip,fp);
    cout<<"ctbc "<< ctbc<<endl;
    
    //ctrx
    readUnsignedInt(ctrx,byte_flip,fp);
    cout<<"ctrx "<<ctrx<<endl;
    
    //cpts
    readUnsignedInt(ctps,byte_flip,fp);
    cout<<"ctps "<< ctps<<endl;

    //ctbe
    readUnsignedInt(ctbe,byte_flip,fp);
    cout<<"ctbe "<<ctbe<<endl;

    //ctps_ctbe
    readUnsignedInt(ctps_ctbe,byte_flip,fp);
    cout<<"ctps ctbe "<<ctps_ctbe<<endl;

    //header end
    readUnsignedInt(header_end,byte_flip,fp);
    cout<<"header end "<< header_end<<endl;

    //slow tfi
    readFloat(slow_tfi,byte_flip,fp);
    cout<<"slow tfi "<< slow_tfi<<" s"<<endl;
    
    //dtn
    readUnsignedInt(dtn,byte_flip,fp);
    cout<<"dtn "<< dtn<<endl;

    //slow_typ
    readUnsignedInt(slow_typ,byte_flip,fp);
    cout<<"slow typ "<< slow_typ<<endl;

    //csr
    readUnsignedInt(csr,byte_flip,fp);
    cout<<"csr "<<csr<<endl;

    //r_mode
    readUnsignedInt(r_mode,byte_flip,fp);
    cout<<"r_mode "<<r_mode<<endl;

    //slow instruction number
    readUnsignedInt(slow_instruction_number,byte_flip,fp);
    cout<<"slow instruction number "<<slow_instruction_number<<endl;

    //bem
    readUnsignedInt(bem,byte_flip,fp);
    cout<<"bem "<<bem<<endl;
    
    //baq mode
    readUnsignedInt(baq_mode,byte_flip,fp);
    cout<<"baq mode "<<baq_mode<<endl;

    //tro
    readFloat(tro,byte_flip,fp);
    cout<<"tro "<<tro<<" s"<<endl;
    
    //rc_bw
    readFloat(rc_bw,byte_flip,fp);
    cout<<"rc_bw "<<rc_bw<<" Hz"<<endl;

    //adc
    readFloat(adc,byte_flip,fp);
    cout<<"Adc "<< adc<<"  Hz"<<endl;

    //att1
    readFloat(at1,byte_flip,fp);
    cout<<"at1 "<<at1<<" "<<endl;
    
    //at3
    readFloat(at3,byte_flip,fp);
    cout<<"at3 "<<at3<<" "<<endl;
    
    //at4
    readFloat(at4,byte_flip,fp);
    cout<<"at4 "<<at4<<" "<<endl;

    //at1 db
    readUnsignedInt(at1_each,byte_flip,fp);
    cout<<"at1 "<<at1_each<<endl;

    //at3 db
    readUnsignedInt(at3_each,byte_flip,fp);
     cout<<"at3 "<<at3_each<<endl;

    //at4 db
    readUnsignedInt(at4_each,byte_flip,fp);
    cout<<"at4 "<<at4_each<<endl;
  
    //rip
    readFloat(rip,byte_flip,fp);
    cout<<"rip "<<rip<<" s"<<endl;

    //csd
    readFloat(csd,byte_flip,fp);
    cout<<"csd "<<csd<<" s"<<endl;

    //rad
    readUnsignedInt(rad,byte_flip,fp);
    cout<<"rad "<<rad<<endl;

    //csq
    readUnsignedInt(csq,byte_flip,fp);
    cout<<"csq "<<csq<<endl;

    //chirp length
    readFloat(chirp_length,byte_flip,fp);
    cout<<"chirp length "<<chirp_length<<" s"<<endl;
   
    //cfs
    readFloat(slow_cfs,byte_flip,fp);
    cout<<"chirp frequency step "<< slow_cfs<<" Hz"<<endl;


    //fast field variables
    readFloat(fast_tfi,byte_flip,fp);
    cout<<"fast tfi "<< fast_tfi<<endl;
  
    //fin
    readUnsignedInt(fin,byte_flip,fp);
    cout<<"fin "<< fin<<endl;

    //fast typ
    readUnsignedInt(fast_typ,byte_flip,fp);
    cout<<"fast typ "<< fast_typ<<endl;

    //pul
    readUnsignedInt(pul,byte_flip,fp);
    cout<<"number of pulses "<< pul<<endl;

    //bii
    readUnsignedInt(bii,byte_flip,fp);
    cout<<"bii "<<bii<<endl;

    //bpd
    readFloat(bpd,byte_flip,fp);
    cout<<"burst period "<<bpd<<" s"<<endl;

    //pri
    readFloat(pri,byte_flip,fp);
    cout<<"pri "<< pri<<" s"<<endl;
  
    //rwd
    readFloat(rwd,byte_flip,fp);
    cout<<"rwd "<<rwd<<" s"<<endl;

    //fast csf
    readFloat(fast_csf,byte_flip,fp);
    cout<<"fast chirp start frequency "<<fast_csf<<" Hz"<<endl;
  

    //footer variables
    readUnsignedInt(iebtth,byte_flip,fp);
    readUnsignedInt(iebttl,byte_flip,fp);
    readUnsignedInt(bgcalls,byte_flip,fp);
    readUnsignedInt(delvmn,byte_flip,fp);
    readUnsignedInt(delvda,byte_flip,fp);
    readUnsignedInt(delvyr,byte_flip,fp);
    cout<<"iebtth iebtt1 "<< iebtth<<" "<<iebttl<<endl;
    cout<<"delvyr da mn "<< delvyr<<" "<<delvmn<<" "<<delvda<<endl;

    readUnsignedInt(cnt_rl,byte_flip,fp);
    readUnsignedInt(cnt_radio,byte_flip,fp);
    readUnsignedInt(cnt_nd,byte_flip,fp);
    readUnsignedInt(eout,byte_flip,fp);
    readUnsignedInt(subr,byte_flip,fp);
    readUnsignedInt(space_craft_time,byte_flip,fp);
    readFloat(hip,byte_flip,fp);
    readFloat(cip,byte_flip,fp);
    cout<<"cnt_rl radio nd "<< cnt_rl<<" "<<cnt_radio<<" "<<cnt_nd<<endl;
    cout<<"hip cip "<<hip<<" "<<cip<<endl;
       
    
    //sub_comm_data variables 
    
    readFloat(fwdtmp,byte_flip,fp);
    readFloat(    be1tmp,byte_flip,fp);
    readFloat(    be2tmp,byte_flip,fp);
    readFloat(    be3tmp,byte_flip,fp);
    readFloat(    be4tmp,byte_flip,fp);
    readFloat(    be5tmp,byte_flip,fp);
    
    
    readFloat(    diptmp,byte_flip,fp);
    readFloat(    rlotmp,byte_flip,fp);
    readFloat(    tadcal1,byte_flip,fp);
    readFloat(    nsdtmp,byte_flip,fp);
    readFloat(    lnatmp,byte_flip,fp);
    readFloat(    evdtmp,byte_flip,fp);
    
    readFloat(    mratmp,byte_flip,fp);
    readFloat(    mruttm,byte_flip,fp);
    readFloat(    dcgttm,byte_flip,fp);
    readFloat(    cucttm,byte_flip,fp);
    readFloat(    twttmp,byte_flip,fp);
    readFloat(    epctmp,byte_flip,fp); 
    
    readFloat(    tw1ttm,byte_flip,fp);
    readFloat(    ep1ttm,byte_flip,fp);
    readFloat(    p_stmp,byte_flip,fp);
    readFloat(    p_sttm,byte_flip,fp);
    readFloat(    fguttm,byte_flip,fp);
    readFloat(    tadcal4,byte_flip,fp); 
    
    readFloat(    esstmp,byte_flip,fp);
    readFloat(    wgb1t1,byte_flip,fp);
    readFloat(    wgb3t1,byte_flip,fp);
    readFloat(    wgb3t2,byte_flip,fp);
    readFloat(    wgb3t3,byte_flip,fp);
    readFloat(    wgb5t1,byte_flip,fp);
    
    readFloat(    pcutmp,byte_flip,fp);
    readFloat(    adctmp,byte_flip,fp);
    readFloat(    tadcal2,byte_flip,fp);
    readFloat(    ecltmp,byte_flip,fp);
    readFloat(    cputmp,byte_flip,fp);
    readFloat(    memtmp,byte_flip,fp);
    
    readFloat(    sadctmp,byte_flip,fp);
    
    readFloat(    tadcal3,byte_flip,fp);
    readFloat(    frwdpw,byte_flip,fp);
    readFloat(    dcgmon,byte_flip,fp);
    readFloat(    lpltlm_db,byte_flip,fp);
  
    readFloat(    nsdcur,byte_flip,fp);
    readFloat(    hpapsm,byte_flip,fp);
    readFloat(    catcur,byte_flip,fp);
    readFloat(    p_smon,byte_flip,fp);
    readFloat(    svlsta,byte_flip,fp);
    readFloat(    usotmp,byte_flip,fp);
    
    readFloat(cpbnkv,byte_flip,fp);
    readFloat(essvlt,byte_flip,fp);
    readFloat(tadcal5,byte_flip,fp);
    readFloat(pcu5v_72,byte_flip,fp);
    readFloat(pcu5i_73,byte_flip,fp);
    readFloat(pcu5v_74,byte_flip,fp);
    
    readFloat(    pcuti_75,byte_flip,fp);
    readFloat(    pcu15v_76,byte_flip,fp);
    readFloat(    pcu15i_77,byte_flip,fp);
    readFloat(    pcu15v_78,byte_flip,fp);
    readFloat(    pcu15i_79,byte_flip,fp);
    readFloat(    pcu12v_80,byte_flip,fp);
    
    readFloat(    pcu12i_81,byte_flip,fp);
    readFloat(    pcucur,byte_flip,fp);
    readFloat(    pllmon,byte_flip,fp);
    readFloat(    ctu5i,byte_flip,fp);
    readFloat(    tadcal6,byte_flip,fp);
    readFloat(    pcu9v_86,byte_flip,fp);
		  
    readFloat(    pcu9i_87,byte_flip,fp);
    readFloat(    pcu9v_88,byte_flip,fp);
    readFloat(    pcu9i_89,byte_flip,fp);
    
    readFloat(   tadcl7,byte_flip,fp);
		  
    readFloat(    shpttm,byte_flip,fp);

    readUnsignedInt(num_bursts_in_flight,byte_flip,fp);
    cout<<"number of bursts in flight "<<num_bursts_in_flight<<endl;

    //read SAR data
    readUnsignedInt(Nradar_data,byte_flip,fp);
    cout<<"number of data points "<< Nradar_data<<endl;

    readFloat(rms_radar_data,byte_flip,fp);
    cout<<"rms value "<<rms_radar_data<<endl;
	
    //----------------	  
    //read Geometry
    //------------------
    readUnsignedInt(quality_flag,byte_flip,fp);
    readDouble(t_sclk,byte_flip,fp);//sclk
    readDouble(t,byte_flip,fp);//ephemeris time
    //read isoc format first
    fread_return= fread(time_utc,sizeof(time_utc),1,fp);
    if(fread_return!=1){
      cerr<<"UTC time read error"<<endl;
      exit(-1);
    }
    for(unsigned int j=0;j<sizeof(time_utc);++j)
      cout<<time_utc[j];
    cout<<endl;
    
    //next isod format
    fread_return= fread(time_utc,sizeof(time_utc),1,fp);
    if(fread_return!=1){
      cerr<<"UTC time read error"<<endl;
      exit(-1);
    }
    for(unsigned int j=0;j<sizeof(time_utc);++j)
      cout<<time_utc[j];
    cout<<endl;

  
    //transmit time offset
    readDouble(transmit_time_offset,byte_flip,fp);
    cout<<"transmit time offset "<< transmit_time_offset<<endl;

    readDouble(time_from_closest_approach,byte_flip,fp);
    readDouble(time_from_epoch,byte_flip,fp);
   

    
    fread_return=fread(target_name,TARGET_NAME_STRING_SIZE,1,fp);
    if(fread_return!=1){
	cerr<<"TARGET NAME reading error"<<endl;
	exit(-1);
      }
    for(int c=0;c<TARGET_NAME_STRING_SIZE;c++) cout<<target_name[c];
    cout<<endl;

    fread_return=fread(tbf_frame_name,TARGET_FRAME_NAME_STRING_SIZE,1,fp);
      if(fread_return!=1){
	cerr<<"TARGET FRAME NAME reading error"<<endl;
	exit(-1);
      }
    for(int c=0;c<TARGET_FRAME_NAME_STRING_SIZE;++c) cout<<tbf_frame_name[c];
    cout<<endl;


    readDouble(pole_right_ascension,byte_flip,fp);
    readDouble(pole_declination,byte_flip,fp);
    readDouble(target_rotation_rate,byte_flip,fp);
    readDouble(target_rotation_angle,byte_flip,fp);

    cout<<"pole right ascention, dec "<< pole_right_ascension<<" deg "<< pole_declination<<" deg "<<endl;
    cout<<"target rotation rate and angle "<< target_rotation_rate<<" deg/s "<<target_rotation_angle<<" deg "<<endl; 

    readFloat(T_scwg,byte_flip,fp);
    readFloat(T_feed,byte_flip,fp);
    readFloat(T_hga,byte_flip,fp);

    readUnsignedInt(beam_number,byte_flip,fp);
    cout<<"used beam number "<< beam_number<<endl;

    //sc position and velocity in J2000
    readDouble(sc_position_J2000_x,byte_flip,fp);
    readDouble(sc_position_J2000_y,byte_flip,fp);
    readDouble(sc_position_J2000_z,byte_flip,fp);
    readDouble(sc_velocity_J2000_x,byte_flip,fp);
    readDouble(sc_velocity_J2000_y,byte_flip,fp);
    readDouble(sc_velocity_J2000_z,byte_flip,fp);

    //sc  position and velocity in target body fixed frame	      
    readDouble(sc_position_target_x,byte_flip,fp);
    readDouble(sc_position_target_y,byte_flip,fp);
    readDouble(sc_position_target_z,byte_flip,fp);
    readDouble(sc_velocity_target_x,byte_flip,fp);
    readDouble(sc_velocity_target_y,byte_flip,fp);
    readDouble(sc_velocity_target_z,byte_flip,fp);


    cout<<"sc position in target frame "<< sc_position_target_x<<" "<<
      sc_position_target_y<<" "<< sc_position_target_z<<" "<<endl;

    cout<<"sc velocity in target frame "<< sc_velocity_target_x<<" "<<
   sc_velocity_target_y<<" "<< sc_velocity_target_z<<" "<<endl;

    cout<<"sc position in J2000 "<< sc_position_J2000_x<<" "<<
      sc_position_J2000_y<<" "<< sc_position_J2000_z<<" "<<endl;
    

    cout<<"sc velocity  in J2000 "<< sc_velocity_J2000_x<<" "<<
      sc_velocity_J2000_y<<" "<< sc_velocity_J2000_z<<" "<<endl;


  

    //sc XYZ direction vectors in J2000
    readDouble(sc_X_J2000_x,byte_flip,fp);
    readDouble(sc_X_J2000_y,byte_flip,fp);
    readDouble(sc_X_J2000_z,byte_flip,fp);
    readDouble(sc_Y_J2000_x,byte_flip,fp);
    readDouble(sc_Y_J2000_y,byte_flip,fp);
    readDouble(sc_Y_J2000_z,byte_flip,fp);
    readDouble(sc_Z_J2000_x,byte_flip,fp);
    readDouble(sc_Z_J2000_y,byte_flip,fp);
    readDouble(sc_Z_J2000_z,byte_flip,fp);
    cout<<"sc in J2000 X "<< sc_X_J2000_x<<" "<<sc_X_J2000_y<<" "<<sc_X_J2000_x<<endl; 
    cout<<"sc in J2000 Y "<< sc_Y_J2000_x<<" "<<sc_Y_J2000_y<<" "<<sc_Y_J2000_x<<endl; 
    cout<<"sc in J2000 Z"<< sc_Z_J2000_x<<" "<<sc_Z_J2000_y<<" "<<sc_Z_J2000_x<<endl;    

    //sc XYZ direction vectors in target body fixed frame
    readDouble(sc_X_target_x,byte_flip,fp);
    readDouble(sc_X_target_y,byte_flip,fp);
    readDouble(sc_X_target_z,byte_flip,fp);
    readDouble(sc_Y_target_x,byte_flip,fp);
    readDouble(sc_Y_target_y,byte_flip,fp);
    readDouble(sc_Y_target_z,byte_flip,fp);
    readDouble(sc_Z_target_x,byte_flip,fp);
    readDouble(sc_Z_target_y,byte_flip,fp);
    readDouble(sc_Z_target_z,byte_flip,fp);

    cout<<"sc in target X "<< sc_X_target_x<<" "<<sc_X_target_y<<" "<<sc_X_target_x<<endl; 
    cout<<"sc in target Y "<< sc_Y_target_x<<" "<<sc_Y_target_y<<" "<<sc_Y_target_x<<endl; 
    cout<<"sc in target Z"<< sc_Z_target_x<<" "<<sc_Z_target_y<<" "<<sc_Z_target_x<<endl;  


    //sc roational velocity in J2000
    readDouble(rot_vel_J2000_x,byte_flip,fp);
    readDouble(rot_vel_J2000_y,byte_flip,fp);
    readDouble(rot_vel_J2000_z,byte_flip,fp);

    //sc roational velocity in target
    readDouble(rot_vel_target_x,byte_flip,fp);
    readDouble(rot_vel_target_y,byte_flip,fp);
    readDouble(rot_vel_target_z,byte_flip,fp);

    readFloat(norm_cnt_rl,byte_flip,fp);
    readFloat(norm_cnt_nd,byte_flip,fp);
    readFloat(norm_cnt_radio,byte_flip,fp);
  cout<<"norm counts (rl radio nd) "<< norm_cnt_rl<<" "<<norm_cnt_radio<<" "<<norm_cnt_nd<<endl; 
    //read Measurement Geometry
    readUnsignedInt(science_qual_flag,byte_flip,fp);
    readFloat(system_gain,byte_flip,fp);
    readFloat(  antenna_brightness_temp,byte_flip,fp);
    readFloat(  system_noise_temp,byte_flip,fp);
    readFloat(  abt_std,byte_flip,fp);
  
  cout<<"Ta Tsys sys_gain dTa " << antenna_brightness_temp <<" "
    <<system_noise_temp<<" "
    <<system_gain<<" "<<abt_std<<endl; 

    readFloat(  pass_geom_time_offset,byte_flip,fp);
    readFloat(  pass_pol_angle,byte_flip,fp);
    readFloat(  pass_emission_angle,byte_flip,fp);
    readFloat(  pass_azimuth_angle,byte_flip,fp);
    readFloat(  pass_centroid_lon,byte_flip,fp);
    readFloat(  pass_centroid_lat,byte_flip,fp);
    readFloat(  pass_major_width,byte_flip,fp);
    readFloat(  pass_minor_width, byte_flip,fp);
    cout<<"passive pole angle "<< pass_pol_angle<<" deg"<<endl;
    cout<<"pass_emission_angle "<<pass_emission_angle<<" deg"<<endl;
    cout<<"pass major minor width "<< pass_major_width<<" km "<< pass_minor_width<<" km "<<endl;
    cout<<"pass azimuth angle "<< pass_azimuth_angle<<" deg"<<endl;

    cout<<"pass centroid lon and lat "<< pass_centroid_lon<<" "<< pass_centroid_lat<<endl;
    readFloat(  pass_ellipse_pt1_lon,byte_flip,fp);
    readFloat(  pass_ellipse_pt2_lon,byte_flip,fp);
    readFloat(  pass_ellipse_pt3_lon,byte_flip,fp);
    readFloat(  pass_ellipse_pt4_lon,byte_flip,fp);
    

    cout<<"pass_ellipse_pt1_lon "<< pass_ellipse_pt1_lon <<" deg"<<endl;
   
    readFloat(  pass_ellipse_pt1_lat,byte_flip,fp);
    readFloat(  pass_ellipse_pt2_lat,byte_flip,fp);
    readFloat(  pass_ellipse_pt3_lat,byte_flip,fp);
    readFloat(  pass_ellipse_pt4_lat,byte_flip,fp);
    readUnsignedInt(num_pulses_received,byte_flip,fp);
    readFloat(total_echo_energy,byte_flip,fp);
    
    readFloat(  noise_echo_energy,byte_flip,fp);
    readFloat(  x_factor,byte_flip,fp);
    readFloat(  sigma0_uncorrected,byte_flip,fp);
    cout<<"total_echo_energy "<< total_echo_energy << endl;
    cout<<"noise_echo_energy "<< noise_echo_energy << endl;
    cout<<"x_factor "<< x_factor << endl;
    readFloat(  sigma0_corrected,byte_flip,fp);
    readFloat(  sigma0_uncorrected_std,byte_flip,fp);
    cout<<"sigma0 uncorr "<< sigma0_uncorrected<<endl;
    cout<<"sigma0 corr "<< sigma0_corrected<<endl;
    
    readFloat(  altitude_means,byte_flip,fp);
    readFloat(  altitude_means_std,byte_flip,fp);
   

    readFloat(  act_geom_time_offset,byte_flip,fp);
    readFloat(  act_pol_angle,byte_flip,fp);
    readFloat(  act_incidence_angle,byte_flip,fp);
    readFloat(  act_azimuth_angle,byte_flip,fp);  
    cout<<"act azimuth angle "<< act_azimuth_angle<< "deg"<<endl;

    cout<<"act_incidence_angle "<< act_incidence_angle<<" deg"<<endl;    
    readFloat(  act_centroid_lon,byte_flip,fp);
    readFloat(  act_centroid_lat,byte_flip,fp);
    readFloat(  act_major_width,byte_flip,fp);
    readFloat(  act_minor_width,byte_flip,fp);
    cout<<"act centroid lon"<< act_centroid_lon<<" deg"<<endl;
    cout<<"act centroid lat "<< act_centroid_lat<<" deg"<<endl;
    cout<<"act major minor width "<< act_major_width<<" km "<< act_minor_width<<" km "<<endl;

    readFloat(  act_ellipse_pt1_lon,byte_flip,fp);
    readFloat(  act_ellipse_pt2_lon,byte_flip,fp);
    readFloat(  act_ellipse_pt3_lon,byte_flip,fp);
    
    readFloat(  act_ellipse_pt4_lon,byte_flip,fp);
    readFloat(  act_ellipse_pt1_lat,byte_flip,fp); 
    readFloat(  act_ellipse_pt2_lat,byte_flip,fp);
    readFloat(  act_ellipse_pt3_lat,byte_flip,fp);
    readFloat(  act_ellipse_pt4_lat,   byte_flip,fp);
    
    readFloat(  altimeter_profile_range_start,byte_flip,fp);
    readFloat(  altimeter_profile_range_step ,byte_flip,fp);
    readUnsignedInt(altimeter_profile_length,byte_flip,fp);
    readFloat(  sar_azimuth_res,byte_flip,fp);
    readFloat(  sar_range_res, byte_flip,fp);
    readFloat(  sar_centroid_bidr_lon,byte_flip,fp);
    readFloat(  sar_centroid_bidr_lat,byte_flip,fp);
    cout<<"sar resolution azi,range: "<< sar_azimuth_res<<"," <<
      sar_range_res << " km"<<endl;
    cout<<"sar centroid lat,wlon: "<< sar_centroid_bidr_lat<<"," <<
      sar_centroid_bidr_lon << " deg"<<endl;
    int end_position=ftell(fp);
    if( (end_position- start_position)!=SBDR_RECORD_LENGTH){
      cerr<<"Reading error: Wrong size of data length has benn read "<<endl;
      exit(-1);
    }
  }
  if(!found_sab_counter){
    cout<<"No Sab record of "<< input_sab_number<<" in this lbdr file"<<endl;
    cout << "First Sab counter: " << sab_counter_start << endl;
    cout << "Last Sab counter: " << sab_counter_end << endl;
  }

  
  return(0);
  
}




void readFloat(float& x, const bool& byte_flip, FILE* fp)
  {
   
    int fread_return=fread(&x,sizeof(float),1,fp);
    if(fread_return!=1){
      cerr<<"read error"<<endl;
      exit(-1);
    }
    if(byte_flip){
      char* p= (char *) &(x);
      for(unsigned int i=1;i<=sizeof(float)/2;++i){
	swap(*p, *(p+sizeof(float)-(i*2-1)));
	++p;
      }
    }
  }

void readDouble(double& x, const bool& byte_flip, FILE* fp)
  {
   
    int fread_return=fread(&x,sizeof(double),1,fp);
    if(fread_return!=1){
      cerr<<"read error"<<endl;
      exit(-1);
    }
    if(byte_flip){
      char* p= (char *) &(x);
      for(unsigned int i=1;i<=sizeof(double)/2;++i){
	swap(*p, *(p+sizeof(double)-(i*2-1)));
	++p;
      }
    }
  }

void readUnsignedInt(unsigned int& i, const bool& byte_flip,FILE* fp)
  {
    
    int fread_return=fread(&i,sizeof(unsigned int),1,fp);
    if(fread_return!=1){
      cerr<<"read error"<<endl;
      exit(-1);
    }
    if(byte_flip){
      char* p= (char *) &(i);
      for(unsigned int i=1;i<=sizeof(unsigned int)/2;++i){
	swap(*p,*(p+sizeof(unsigned int)-(i*2-1)));
	++p;
      }
    }
  }
