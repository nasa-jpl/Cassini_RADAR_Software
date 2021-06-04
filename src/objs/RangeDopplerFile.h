//-----------------------------------------------
//RangeDopplerFile.h
// This file contains the RangeDopplerFile class declaration.  It
// defines the file format, and provides methods to read and write
// range_doppler records
//--------------------------------------------------


#ifndef RangeDopplerFile_H
#define RangeDopplerFile_H

//----------------------------
//Forward declaration
//----------------------------
class RangeDopplerFile;

#include <string>
#include <vector>
#include <map>
#include "Array.h"
#include "Error.h"
#include "Units.h"
#include "Time.h"
#include "Io.h"

using std::string;

//using std::endl;
//-------------------------
// Class RasList declaration
//-------------------------
class RangeDopplerFile
  {
  public: 
    //-------------- 
    // Construction
    //-------------- 
    RangeDopplerFile(const string& filename, const string& mode); 

    ~RangeDopplerFile();
    //-----------------------------------------
    //Input method load sab from L0 or echo simulation program
    //-----------------------------------------
   
    void loadBurstData(const Time& burst_time,
		       const Uvar& active_time_offset,
		       const unsigned int& n_sab,
		       const unsigned int beam_id_input,
		       const Uvar& prf_input,
		       const Uvar& lambda_input);
    void loadBurstData(const Time& burst_time,
		       const Uvar& active_time_offset,
		       const unsigned int& n_sab,
		       const unsigned int beam_id_input,
		       const Uvar& prf_input);

    void loadComputedRangeDoppler( const unsigned int& Npoint,
				   const Uvec& range_values,
				   const Dvec& fractional_doppler,
				   const Uvar& bore_range_input,
				   const double& fract0D,
				   const int& bore_doppler_integer_input);

    void loadMeasuredRangeDoppler( const unsigned int& Npoint,
				   const Uvec& range_values,
				   const Dvec& fractional_doppler,
				   const Uvar& bore_range_input,
				   const double& fract0D,
				   const int& bore_doppler_integer_input);
    void loadMeasuredRangeFractDoppler( const unsigned int& Npoint,
					const Uvec& range_values,
					const Dvec& fractional_doppler,
					const Uvar& bore_range_input,
					const double& fract0D);
    void saveLoadedData();

    //------------------------
    //writing int doppler back to file
    //-----------------------
    void writeIntDoppler(const unsigned int& n_sab,
			 const int& int_doppler);

    //------------------------
    //find a  record
    //------------------------  
    bool findRecord(const unsigned int& n_sab);
    bool findRecord(const Time& t);
    bool validGeomData();
    bool validMeasData();
    bool validIntDoppler();
    //--------------------------
    //record manipulation
    //----------------------------
    unsigned int& operator[] (const unsigned int& n);
   
    
    //--------------------------
    //compute record number
    //---------------------------
    unsigned int computeNumberOfRecords();
    void getFirstLastSabNumber(unsigned int& first_sab, unsigned int& last_sab);
  
  

    //------------
    //get beam id of sab
    //-------------
    unsigned int getBeamId(const unsigned int& sab_input);
  
    //--------------------
    //plot echo, matched filter power spectrum
    //---------------------
    void displayRangeDoppler(const unsigned int& n_sab);
  
    //-------------------------------------
    //close the file
    //--------------------------------------
    void close();

    //------------------------------------------
    //public data folder    
    //These parameters are extracted from SAB
    //-------------------------------------------
    Time t;
    Uvar time_offset;
    short unsigned int   data_type;
    unsigned int sab_counter;
    unsigned int beam_id;
    Uvar prf;
    Uvar lambda;
    unsigned int Ndata;
    Uvec range;
    Dvec fract_doppler;
    Uvar bore_range;
    int bore_int_doppler;
    double bore_fract_doppler;
    
    unsigned int Ndata_meas;
    Uvec range_meas;
    Dvec fract_doppler_meas;
    Uvar bore_range_meas;
    int bore_int_doppler_meas;
    double bore_fract_doppler_meas;
    //-----------------------
    //Priviate variables/methods
    //-----------------------
  private:
    FileMgr file_;
    string filename_;
    string mode_;
    unsigned int sync_;
    unsigned int Nmax_;   
   
    //----------------------------------------------------------
    //useful record containter that holds sab, time
    // file position
    //these containers are set onece when rasfile is opened
    //-----------------------------------------------------
    //list of sab number
    vector<unsigned int> sabnumber_list_;
    vector<unsigned int>::const_iterator  sabnumber_list_ptr_;
    //sab number vs file position
    map<unsigned int, int> sabnumber_vs_fileposition_;
    map<unsigned int,int>::const_iterator sabnumber_vs_fileposition_ptr_;
    //list of time for each record
    vector<Time> time_list_;
    vector<Time>::const_iterator time_list_ptr_;
    //time vs sab record
    map<Time, unsigned int> time_vs_sabnumber_;
    map<Time, unsigned int>::const_iterator time_vs_sabnumber_ptr_;
    //sab number vs beam number
    map<unsigned int, unsigned int> sabnumber_vs_beamid_;
    map<unsigned int, unsigned int>::const_iterator sabnumber_vs_beamid_ptr_;

  
    void mapAllRecords();
    void readRecord();
    void writeRecord();
  
  };


#endif



