//-------------------------------------------
//Read ccsds packets directly to produce L0 ROME file
//-----------------------------------------------

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "Config.h"
#include "Io.h"
#include "Error.h"
#include "Sab.h"
#include "Config.h"




using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::cerr;
using std::hex;
using std::oct;

void myUnexpected() throw()
  {
  cout << "Unexpected exception" << endl;
  std::terminate();
  }

int  readByteFlipLong( const int& input, const bool& flip)
  {
    if(!flip) return(input);//no need to flip bytes

    int i= input;
    char* p= (char *) &(i);
    for(unsigned int i=1;i<=sizeof(unsigned int)/2;++i){
      swap(*p,*(p+sizeof(unsigned int)-(i*2-1)));
      ++p;
    }
    return(i);
  }


int main(int argc, char* argv[])

{

std::set_unexpected(myUnexpected);

try
{

  //------------------------
  // Parse the command line
  //------------------------

  const char* command = argv[0];
  if (argc != 3)
    {
    cerr << "Usage: " << command << " config_file  units_mode" << endl;
    
    exit(-1);
    }
  

  char* config_file =NULL;
  char* units_mode =NULL;

  //-----------------
  //Parse the input options
  //----------------
  int clidx = 1;
  config_file = argv[clidx++];
  units_mode = argv[clidx++];
 

  //----------------------------
  // Set unit mode
  //----------------------------
  set_unit_option(units_mode);

  //------------------------------
  // Load configuration parameters
  //------------------------------
  Config cfg(config_file);

  
 
 
  //--------------------------------
  //set Level0 filename
  // and read record counts
  //--------------------------------
  bool byte_flip=false;
  string input_filename = cfg.str("ccsds_sci_file");
  FileMgr input_file(input_filename,"rb");
  unsigned short header= 20042;
  unsigned short header2;
  input_file.read(header2);//read as default
  input_file.rewind();
  
  input_file.readKnown(header);
  if(header != header2) byte_flip=true;
  cout<<"machine type "<< input_file.machineType()<<endl;
  cout<<"Ready to read input ccsds file "<< input_filename<<endl;
  if(byte_flip) cout<<"Byte flip: different machine types"<<endl;
  input_file.rewind();
  
  

  //------------------------
  //EGSE file
  //------------------------
  string output_filename=cfg.str("egse_filename");
  int l0_record_type= Sab::archive_header;
  int l0_record_count=0;//temporaray setting
  vector<unsigned char> empty_header; 
  empty_header.resize(40);
  int l0_record_length=84;
  FileMgr output_file(output_filename,"wb");

  
 
  output_file.write(readByteFlipLong(l0_record_type,byte_flip));
  output_file.write(readByteFlipLong(l0_record_length,byte_flip));
  output_file.write(readByteFlipLong(l0_record_count,byte_flip));

  output_file.write(empty_header);
  output_file.write(empty_header);

  //---------------
  //recycled variable
  //--------------
  unsigned int version_number, type, header_flag,system_id_header,system_id_tail, packet_type;
  unsigned int packet_length,segment_flag,cds_error_flag;
  int seq_count,  seq_count_save;
  unsigned int i_count=0;//ccsds counter
  seq_count=0;
  seq_count_save=0;
  unsigned int i_sab=0;//sab counter
 
  //----------------------------
  //Read record length and count
  //-----------------------------
  vector<unsigned char>  scisfdu_header, inspect_header;
  vector<unsigned char> ccsds_header;
  vector<unsigned char> ccsds_body;
  vector<unsigned char> sab_partial;
  vector<unsigned char> l0_container; 
  vector<unsigned char> wtkj;
  vector<unsigned char> sab_full;
  wtkj.push_back('w');
  wtkj.push_back('t');
  wtkj.push_back('k');
  wtkj.push_back('j');
  scisfdu_header.resize(150);//sfdu overhead
  inspect_header.resize(150);
  ccsds_header.resize(12);
  ccsds_body.resize( 950-12);
  sab_partial.resize(2);
  l0_container.clear();
  unsigned int ccsds_length=1100;
  //----------------------
  //count how many ccsds 
  // each ccsds header has 1100 bytes
  //--------------------
  input_file.rewind();
  unsigned int Nccsds=0;
  while(!input_file.eof()){
    input_file.setPosition(input_file.getPosition()+1100);
    Nccsds++;
  }
  if(!input_file.eof()) ErrorMessage("file size is not multipole of 1100 bytes").throwMe();
  cout<<"Total number of packets "<< Nccsds<<endl;
  input_file.rewind();


  //------------------
  //first header
  //------------------
  input_file.read(inspect_header);
  input_file.rewind();

  //----------------
  //Sab
  //----------------
  Sab sab(input_file.machineType());
  sab.clear();
  int sab_type=Sab::sci;
  int sab_length;

  unsigned int count_missing_ccsds=0;
  unsigned int count_dropped_sabs=0;
  unsigned int count1_dropped_sabs=0;
 
  //----------------
  //Read
  //------------------

  for(unsigned int i_ccsds=0; i_ccsds<Nccsds;++i_ccsds){
    input_file.setPosition(i_ccsds*ccsds_length);
    int p_start= input_file.getPosition();
    //cout<<i_ccsds <<endl;

    //-----------------
    //Read scisfdu header: 150 bytes
    //-------------------
    input_file.read(scisfdu_header);

    //Not doing anything with scisfdu header 
    

    //----------------
    //Read ccsds_header: 12 bytes
    //------------------
    input_file.read(ccsds_header);
  
    
   
    //-------------------------------------
    //CCSDS header error checking
    //-------------------------------------
    version_number=bitget(ccsds_header[0],5,7);
    if(version_number!=0){
      cout<<"Wrong version number "<< version_number<<endl;
      continue;//next ccsds packet
    }

    type=bitget(ccsds_header[0],4,4);
    if(type!=0){
      cout<<"Wrong type "<< type<<endl;
      //ErrorMessage("type is not correct").throwMe();
      continue;
    }
    
    header_flag=bitget(ccsds_header[0],3,3);
    if(header_flag!=1){
      cout<<"Wrong header flag "<<endl;
      //ErrorMessage("header flag is wrong ").throwMe();
      continue;
    }
    system_id_header=bitget(ccsds_header[0],0,2);//expect 111
    system_id_tail=bitget(ccsds_header[1],6,7);//expect 01
    if(system_id_header!=7 || system_id_tail!=1){
      cout<<"system id header tail "<< system_id_header<<" "<<system_id_tail<<endl;
      //ErrorMessage("system id error").throwMe();
      continue;
    }
    
    packet_type=bitget(ccsds_header[1],0,5);
    if(packet_type!=12){
      cout<<"read packed type "<< packet_type<<endl;
      cout<<"Wrong packet type "<<endl;
      //ErrorMessage("wrong packet type").throwMe();
      continue;
    }
    
    segment_flag=bitget(ccsds_header[2],6,7);
    if(segment_flag!=3){
      cout<<"Wrong segment flag "<<segment_flag<<endl;
      cout<<"I am gong to quit now, good bye "<<endl;
      // ErrorMessage("wrong segment flag ").throwMe();
      continue;
    }
    

    bool missing_ccsds=false;//assume no missing ccsds
    seq_count= bitget(ccsds_header[2],0,5)*256 + bitget(ccsds_header[3],0,8);
    
    if( ((seq_count/100)*100) == seq_count){
      cout<< "CCSDS counter "<< seq_count<<endl;
    }
   
    if((seq_count -seq_count_save)!=1){
      if(seq_count_save==0){
	//do nothing, normal case
      }
      else if(seq_count==0 && seq_count_save==16383){
	cout<<"seq count roll over"<<endl;
      }
      else{
	missing_ccsds=true;//missing ccsds
	cout<<"Missing CCSDS packets"<<endl;
	count_missing_ccsds += seq_count - seq_count_save - 1;
      }
    }
    //cout<<"ccsds counter "<< seq_count<<endl;
    
    packet_length= bitget(ccsds_header[4],0,7)*256 + bitget(ccsds_header[5],0,7);
    if(packet_length!=943){
      cout<<"wrong package length "<< packet_length<<endl;
      ErrorMessage("packet length error").throwMe();
    }
    
    cds_error_flag=bitget(ccsds_header[6],0,7);
    if(cds_error_flag!=0){
      cout<<"Wrong cds error flag "<<endl;
      ErrorMessage("Error in csd error flag").throwMe();
    }
    //------------------------------
    //End of CCSDS error checking
    //-------------------------------   



    //----------------------------
    //Read ccsds body: 938 bytes : one ccsds size= 150 + 12 + 938=1100
    //--------------------------
    input_file.read(ccsds_body);
    
   
    bool wtkj_flag=false;
    unsigned int wtkj_index=0;
    if(missing_ccsds){

      //------------------
      //Find wtkj
      //--------------------
      for(unsigned int i=0;i<938-4;++i){
	if(ccsds_body[i]=='w' &&
	   ccsds_body[i+1]=='t' &&
	   ccsds_body[i+2]=='k' &&
	   ccsds_body[i+3]=='j'){
	  wtkj_flag=true;
	  wtkj_index=i;
	  break;
	}
      }
      if(wtkj_flag){
	cout<<"CCSDS after gap "<<endl;
	cout<<"Find wtkj in ccsds packet"<<endl;
	cout<<" starting index "<< wtkj_index<<endl;
      }
      
      //do not need to add sab any more
      cout<<"-------------------------"<<endl;
      cout<<"Missing CCSDS packet(s) "<<endl;
      cout<<"Previous seq "<< seq_count_save<<endl;
      cout<<"Current seq "<< seq_count<<endl;
      if(sab.complete()){
	//----------------------------------------------------------
	//write into file if it happens to be that Sab was full before
	// data gap occurs
	//-----------------------------------------------------------
	i_sab++;
	sab.decode();
	cout<<"current sab counter "<< sab.sab_counter<<endl;

	sab.get_raw_data(sab_full);
	if(( sab.sab_len*2+124) != sab_full.size()){
	  cout<<"expected sab length "<< sab.sab_len*2+124;
	  cout<<"actual sab size "<< sab_full.size()<<endl;
	  ErrorMessage("sab length mismatch").throwMe();
	}
	


	sab_length= sab_full.size();
       
	output_file.write(readByteFlipLong(sab_type,byte_flip));//record type
	output_file.write(readByteFlipLong(sab_length,byte_flip));//record length
	output_file.write(sab_full);
	
	l0_container.clear();
	sab.clear();
	//
	//cout<<"completed number of sabs "<<i_sab<<endl;
	//cout<<"sab lenght "<< sab_full.size()<<endl;
	//
      }
      else{
	cout<<"----------------------"<<endl;
	cout<<"Due to data gap, this sab becomes invalid"<<endl;
	cout<<"------------------------"<<endl;
	count_dropped_sabs++;
	if (seq_count - seq_count_save == 2) count1_dropped_sabs++;
	sab.clear();
	l0_container.clear();
      }
    }


    
    for(unsigned int i=0;i<938/2;++i){
      l0_container.push_back(ccsds_body[2*i]);
      l0_container.push_back(ccsds_body[2*i+1]);
      sab_partial[0]=ccsds_body[2*i];
      sab_partial[1]=ccsds_body[2*i+1];
      if(l0_container.size()<4) continue;
      if(l0_container.size()==4){
	if(l0_container!=wtkj){
	  l0_container.erase(l0_container.begin());
	  l0_container.erase(l0_container.begin());
	  sab.clear();
	  continue;//keep reading 
	}
      }
      if(l0_container.size()==4){
	sab.addPartialRecord(l0_container);
      }
      else if(l0_container.size()>4){
	sab.addPartialRecord(sab_partial);
      }
      
      if(sab.complete()){
	i_sab++;
	sab.decode();
	//cout<<"current sab counter "<< sab.sab_counter<<endl;

	sab.get_raw_data(sab_full);
	sab_length= sab_full.size();

	if(( sab.sab_len*2+124) != sab_full.size()){
	  cout<<"expected sab length "<< sab.sab_len*2+124;
	  cout<<"actual sab size "<< sab_full.size()<<endl;
	  ErrorMessage("sab length mismatch").throwMe();
	}

	
	output_file.write(readByteFlipLong(sab_type,byte_flip));//record type
	output_file.write(readByteFlipLong(sab_length,byte_flip));//record length

	output_file.write(sab_full);

	l0_container.clear();
	sab.clear();
	//cout<<"completed number of sabs "<<i_sab<<endl;
	//cout<<"sab lenght "<< sab_full.size()<<endl;
      }
    }
    
       
    
    i_count++;
    int p_end=input_file.getPosition();
    if((p_end-p_start)!=1100)
      ErrorMessage("You read more than 1100 bytes").throwMe();
    
    
    seq_count_save= seq_count;//save it
  }
 

  //---------------
  //write back into l0
  //---------------
  l0_record_count=i_sab;
  output_file.rewind();
  output_file.write(readByteFlipLong(l0_record_type,byte_flip));
  output_file.write(readByteFlipLong(l0_record_length,byte_flip));
  output_file.write(readByteFlipLong(l0_record_count,byte_flip));
  output_file.close();

  cout<<"total ccsds packets decoded "<< i_count<<endl;
  cout<<"total ccsds packets "<< Nccsds<<endl;
  cout<<"total sab "<< i_sab<<endl;
  cout << "missing ccsds records " << count_missing_ccsds << endl;
  cout << "missing SAB records " << count_dropped_sabs << endl;
  cout << "missing SAB records due to one missing ccsds packet "
    << count1_dropped_sabs << endl;

  
 }//end of try

catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
catch(...)
  {
  cerr << "Exception caught" << endl;
  }

return(0);

}


