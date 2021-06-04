#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "Ivd.h"
#include "Array.h"
#include "Config.h"
#include "Io.h"
#include "Units.h"
#include "Error.h"
#include "Config.h"
#include "Frame.h"
#include "Time.h"
#include "Constants.h"
#include "TemplateUtils.h"
#include "SpiceUsr.h"
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;


//------------------------
//Static const member initialization for class Ivd
//------------------------
const unsigned int Ivd::maxrec_ = 20000;
//----------------------------//
//  Methods for Class Ivd     //
//----------------------------//

//--------------
// Constructors
//--------------
Ivd::Ivd() 
  :quats("quaternions"), 
   avvs("avvs"),
   attitude_x("attitude_x"),
   attitude_mz("attitude_mz"),
   et_x("et_x"), 
   et_mz("et_mz"),
   sclkdp_x("sclkdp_x"),
   sclkdp_mz("sclkdp_mz"),
   SCNAME_("CASSINI")
 
  {
    quats.resize(maxrec_,4);
    avvs.resize(maxrec_,3);
    attitude_x.resize(maxrec_,3);
    attitude_mz.resize(maxrec_,3);
    et_x.resize(maxrec_);
    et_mz.resize(maxrec_);
    sclkdp_x.resize(maxrec_);
    sclkdp_mz.resize(maxrec_);
    quats = 0.0;
    avvs = 0.0;
    attitude_x = 0.0;
    attitude_mz = 0.0;
    et_x = Uvar(0.0,"s");
    et_mz = Uvar(0.0,"s");
    sclkdp_x = 0.0;
    sclkdp_mz = 0.0;
    ivd_data_loaded = false;
    N_dat = 0;
    
  }


//----------------------------------------------------------------
// void ReadDataFiles() throw(ErrorMessage);
//----------------------------------------------------------------
void Ivd::ReadDataFile(ifstream& ivd_file) throw(ErrorMessage)
  {
    //local variables
    string str1,str2,str3;
    string input_string;
    string time_string;
    unsigned int str_length;

    //--------------------------
    //Read first 32 lines: header
    //---------------------------
    for (unsigned int i = 0 ; i < 32; i++)
    {
     getline(ivd_file,input_string);
     //     cout<<input_string<<endl;
    }
    
    //---------------------
    //Read start time of plus x
    //---------------------
    ivd_file >>str1;
    ivd_file >>str2;
    str_length = str2.length();  
    time_string.erase();
    for (unsigned int i = 0; i < str_length; i++)
    { 
    if(i>0 && i<(str_length-1)) 
      {
	time_string += str2[i];
	//cout<<"time string "<< time_string<< " "<< str2[i]<<endl;
      }
    }
    start_x.setUtc(time_string);
    time_string.erase(); 
    cout<<"start_x time_string: "<< start_x.utc("ISOD")<<endl;

   
    //-----------------------------
    //Read end time_string of plus x
    //------------------------------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    { 
      if(i>0 && i< (str_length-1)) {
	time_string += str2[i];
	//cout<<"time string "<< time_string<<" "<< str2[i]<<endl;
      }
    }
    end_x.setUtc(time_string);
    time_string.erase();
    cout<<"end_x time_string: "<< end_x.utc("ISOD")<<endl;  
    
    //----------
    //Read head of plus x
    //----------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    { 
      if(i>0 && i<( str_length-1) ) {head_x += str2[i];}
    }
    cout<<"Head_x: "<< head_x<<endl;

    //-------------------------------------------------
    //read base of plus x
    //-------------------------------------------------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    {    
    if(i>0 && i< (str_length-1)) {base_x += str2[i];}
    }
    cout<<"Base_x: "<< base_x<<endl;
     
    //-----------------------------------------------------------
    //Read in time and attitude data until end_time_indicator = end
    //-----------------------------------------------------------
    cout.precision(20);
    unsigned int i_data_count_x = 0;
    Time t = start_x;
        
    for (i_data_count_x = 0; i_data_count_x < maxrec_; ++i_data_count_x)
      {
      double x,y,z;

      getline(ivd_file,str1);

      ivd_file>>str1;
      if (str1 !="Time:") 
        {  // first invalid record ends the x stuff
        break;
        }

      ivd_file>>str2;
      str_length = str2.length();  
      time_string.erase();
      
      for (unsigned int i = 0; i < str_length; i++)
	{ 	
	  if(i>0 && i< (str_length-1)) {time_string += str2[i];}
	}
      //cout << i_data_count_x << ": " << time_string << endl;
      t.setUtc(time_string);
      et_x(i_data_count_x) = t.et(); //save as et to construct quaternions
      sclkdp_x(i_data_count_x) = t.encodedSclk(SCNAME_);
      
      ivd_file>>str3;
      if(str3 != "Position:") 
	throw ErrorMessage("Error in reading Position:");
      
      ivd_file>>x;
      attitude_x(i_data_count_x,0) = x;
      
      ivd_file>>y;
      attitude_x(i_data_count_x,1) = y;
      
      ivd_file>>z;
      attitude_x(i_data_count_x,2) = z;
      }

    if (i_data_count_x >= maxrec_) 
      {
      ErrorMessage e("Increase maxrec_ > " + toStr(maxrec_));
      e.throwMe();
      }

    // transfer total number of data to public variable for later use
    N_dat = i_data_count_x;

   
    //getline(ivd_file,str1);//read in blank line
    time_string.erase();//clear string

    //----------------
    //Read start time of -z
    //---------------
    //ivd_file>>str1; 
    ivd_file>>str2;
    //cout << str1 << " " << str2 << endl;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    { 
      if(i>0 && i< (str_length-1)) {time_string += str2[i];}
    }
    //cout << time_string << endl;
    start_mz.setUtc(time_string);
    time_string.erase(); 
    if (start_mz != start_x) 
      throw ErrorMessage("Start time mismatch:start_x and start_mz");
    //cout<<"start mz  "<<start_mz<<endl;

    //-----------------
    //Read end time of -z
    //-----------------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    { 
    if(i>0 && i< (str_length-1)) {time_string += str2[i];}
    }
    end_mz.setUtc(time_string);
    time_string.erase();
    if(end_mz !=end_x) 
      throw ErrorMessage("End time mismatch: end_x and end_mz");
    //cout<<"end mz "<<end_mz<<endl;

    //------------------------
    //Read head of -z
    //------------------------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    { 
    if(i>0 && i< (str_length-1)) {head_mz += str2[i];}
    }
    //cout<<"Head_mz: "<< head_mz<<endl;
    
    //----------------------
    //Read base of -z
    //---------------------
    ivd_file>>str1; 
    ivd_file>>str2;
    str_length = str2.length();  
    for (unsigned int i = 0; i < str_length; i++)
    {    
    if(i>0 && i< (str_length-1)) {base_mz += str2[i];}
    }
    //cout<<"Base_mz: "<< base_mz<<endl;

    unsigned int i_data_count_mz = 0;
    t = start_mz;

    for (i_data_count_mz = 0; i_data_count_mz < maxrec_; ++i_data_count_mz)
      {
      double x,y,z;
      getline(ivd_file,str1);//read in blank line
      
      ivd_file>>str1;
      if (!ivd_file || str1 !="Time:") 
        {  // first invalid record ends the mz stuff
        break;
        }
      
      ivd_file>>str2;
      str_length = str2.length();  
      time_string.erase();
      for (unsigned int i = 0; i < str_length; i++)
	{ 	
	  if(i>0 && i<( str_length-1)) {time_string += str2[i];}
	}
      //cout << i_data_count_mz << ": " << time_string << endl;
      t.setUtc(time_string);
      sclkdp_mz(i_data_count_mz)=t.encodedSclk(SCNAME_);
      et_mz(i_data_count_mz) = t.et();
      
      ivd_file>>str3;
      if(str3 != "Position:")
	throw ErrorMessage("Error in reading Position:");
      
      ivd_file>>x;
      attitude_mz(i_data_count_mz,0) = x;
      
      ivd_file>>y;
      attitude_mz(i_data_count_mz,1) = y;
      
      ivd_file>>z;
      attitude_mz(i_data_count_mz,2) = z;
      }

    if (i_data_count_mz >= maxrec_) 
      {
      ErrorMessage e("Increase maxrec_ > " + toStr(maxrec_));
      e.throwMe();
      }

    //Compare two data sets and report any discrepancies if there is any
    
    if (N_dat != i_data_count_mz) 
      throw ErrorMessage("Two files have different number of data");
    ivd_data_loaded = true;
    
    /*
    cout<<"Number of record "<<N_dat<<endl;
    for (unsigned int i = 0; i <N_dat;i++)
    {
    cout<<sclkdp_x(i)<<endl;
    cout<<setw(15)<<attitude_x(i,0)<<" "<<setw(15)<<attitude_x(i,1)<<" "<<setw(15)<<attitude_x(i,2)<<endl;
    }

    cout<<"Number of record "<<i_data_count_mz<<endl;
    for (unsigned int i = 0; i <i_data_count_mz;i++)
    {
    cout<<sclkdp_mz(i)<<endl;
    cout<<setw(15)<<attitude_mz(i,0)<<" "<<setw(15)<<attitude_mz(i,1)<<" "<<setw(15)<<attitude_mz(i,2)<<endl;
    }
    */

   
  }




