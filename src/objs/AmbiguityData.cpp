#include <iomanip>
#include "AmbiguityData.h"
#include "TemplateUtils.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::setw;

//---------------------------
// Methods for AmbiuityData
//---------------------------

//--------------
// Constructors
//--------------

AmbiguityData::AmbiguityData() throw(ErrorMessage)
  :max_beam_gain("max_beam_gain"),
   tau_3db("tau_3dB"),
   pulse_gate("pulse_gate"),
   range_offset("range_offset"),
   frequency_offset("frequency_offset"),
   range_3db("range_3db"), 
   X0("X0"),
   dummImat("dummy"),
   dummDmat("dummy"),
   usable_area("usable_area"),   
   range_lower("range_lower"),
   range_upper("range_upper"),
   dop_lower("dop_lower"),
   dop_upper("dop_upper"),
   range_center("range_center"),
   dop_center("dop_center"),
   area_process_window("area_process_window"),
   avg_incidence("avg_incidence"),
   avg_backscatter("avg_backscatter"),
   avg_antenna("avg_antenna"),
   avg_azi("avg_azi"),
   avg_elev("avg_elev"),
   avg_area("avg_area"),
   lat("lat"), 
   lon("lon"), 
   range("range"),
   alongtrack("alongtrack"),
   crosstrack("crosstrack"),
   doppler("doppler"), 
   thetai("theta"), 
   center_area("center_area"),
   bin_ratiodB("bin_ratiodB"),
   normalized_bin_power("normalized_bin_power"),
   bin_antenna("bin_antenna"),
   bin_antennadB("bin_antennadB"),
   is_good_bin("is_good_bin"),
   is_good_bin_amb("is_good_bin_amb"),
   is_good_bin_nesigma0("is_good_bin_nesigma0"),
   is_good_bin_totalSNR("is_good_bin_totalSNR"),	
   is_good_bin_pn_to_amb("is_good_bin_pn_to_amb"),
   patch_sum("patch_sum"),    
   bin_power_main("bin_power_main"),
   bin_power_amb("bin_power_amb"),
   Return_signal("return_signal"),
   total_noise("total_noise"),
   center_patch_along("center_patch_along"),
   center_patch_cross("center_patch_cross"),
   read_data_(false),
   Nrange_patch_(0),
   Ndop_patch_(0),
   i_range_center_(0),
   j_dop_center_(0),
   i_patch_center_(0),
   Nrange_bin_(0),
   Ndop_bin_(0),
   epoch_time_(" ")

  {
  max_beam_gain.resize(5);
  tau_3db.resize(5);
  pulse_gate.resize(5);
  range_offset.resize(5);
  frequency_offset.resize(5);
  range_3db.resize(5); 
  X0.resize(5);
  dummImat.resize(5,5);
  dummDmat.resize(5,5);
  range_center.resize(5);
  dop_center.resize(5);
  area_process_window.resize(5);
  usable_area.resize(5); 
  }

//-----------------------------------------------
//Read amb file
//---------------------------------------------
void AmbiguityData::readAmbiguityData(ifstream& amb_data_file) throw(ErrorMessage)
  {

  double x;
  getline(amb_data_file,epoch_time_);
  amb_data_file>>x; sar_imaging_time_ = Uvar(x,"min");

  //------------------------------------------------
  //first read number of range/dop patches and range/dop cells
  //and resize all relevant arrays
  //---------------------------------------------------
  unsigned int i_input;
  amb_data_file >>i_input ; Nrange_patch_=i_input;
  amb_data_file>>i_input; Ndop_patch_= i_input;
  amb_data_file>>i_input;Nrange_bin_=i_input;
  amb_data_file>>i_input; Ndop_bin_=i_input;

  i_range_center_ = int(Nrange_patch_/2);
  j_dop_center_= int(Ndop_patch_/2);
  i_patch_center_ = int(Nrange_patch_ * Ndop_patch_/2);
  //-------------------------------------------
  //Resize and initialize Array Variables
  //-------------------------------------------   
  range_lower.resize(5,Nrange_patch_);
  range_upper.resize(5,Nrange_patch_);
  dop_lower.resize(5,Nrange_patch_);
  dop_upper.resize(5,Nrange_patch_);
  
  avg_incidence.resize(5,Nrange_patch_,Ndop_patch_);
  avg_backscatter.resize(5,Nrange_patch_,Ndop_patch_);
  avg_antenna.resize(5,Nrange_patch_,Ndop_patch_);
  avg_azi.resize(5,Nrange_patch_,Ndop_patch_);
  avg_elev.resize(5,Nrange_patch_,Ndop_patch_);
  avg_area.resize(5,Nrange_patch_,Ndop_patch_);
    
  lat.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);   
  lon.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);
  range.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);  
  doppler.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);    
  thetai.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);
  alongtrack.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);
  crosstrack.resize(5,Nrange_patch_*Ndop_patch_,Nrange_bin_*Ndop_bin_);
     
  bin_ratiodB.resize(5,Nrange_bin_,Ndop_bin_);
  patch_sum.resize(5,Nrange_bin_,Ndop_bin_);
  bin_power_main.resize(5,Nrange_bin_,Ndop_bin_);
  bin_power_amb.resize(5,Nrange_bin_,Ndop_bin_);
  normalized_bin_power.resize(5,Nrange_bin_,Ndop_bin_);
  center_area.resize(5,Nrange_bin_,Ndop_bin_);
  Return_signal.resize(5,Nrange_bin_,Ndop_bin_);
  total_noise.resize(5,Nrange_bin_,Ndop_bin_);
  center_patch_along.resize(5,Nrange_bin_,Ndop_bin_);
  center_patch_cross.resize(5,Nrange_bin_,Ndop_bin_);
  
  is_good_bin.resize(5,Nrange_bin_,Ndop_bin_);
  is_good_bin_amb.resize(5,Nrange_bin_,Ndop_bin_);
  is_good_bin_nesigma0.resize(5,Nrange_bin_,Ndop_bin_);
  is_good_bin_totalSNR.resize(5,Nrange_bin_,Ndop_bin_);
  is_good_bin_pn_to_amb.resize(5,Nrange_bin_,Ndop_bin_);
  bin_antenna.resize(5,Nrange_bin_,Ndop_bin_);
  bin_antennadB.resize(5,Nrange_bin_,Ndop_bin_);
    
  range_lower = Uvar(0.0,"km");
  range_upper = Uvar(0.0,"km");
  avg_incidence = Uvar(0.0,"deg");
  avg_azi = Uvar(0.0,"deg");
  avg_elev =Uvar(0.0,"deg");
  bin_power_main=Uvar(0.0,"kg km km/(s s s)");
  bin_power_amb=Uvar(0.0,"kg km km/(s s s)");
  normalized_bin_power=Uvar(0.0,"kg km km/(s s s)");
  center_area = Uvar(0.0,"km km");
  area_process_window=Uvar(0.0,"km km");
  


  for (unsigned int i_beam = 0; i_beam < 5; ++i_beam)
    {  
    amb_data_file>>i_input;
    amb_data_file>>max_beam_gain(i_beam);
    amb_data_file>>x;frequency=Uvar(x,"Hz");
    amb_data_file>>x;tau_p=Uvar(x,"us");
    amb_data_file>>x;Pt=coerce_base(x,"W");
    amb_data_file>>x;Tsys=Uvar(x,"K");
    amb_data_file>>x;BR=Uvar(x,"Hz");
    amb_data_file>>x;Bn=Uvar(x,"Hz");
    amb_data_file>>x;fadc=Uvar(x,"Hz");  
    amb_data_file>>x;radius=Uvar(x,"km");
    amb_data_file>>x;Pn=coerce_base(x,"W");    
    amb_data_file>>x;prf=Uvar(x,"Hz");
    amb_data_file>>x;min_bpd=Uvar(x,"ms");
    amb_data_file>>x;pbw=Uvar(x,"Hz");
    amb_data_file>>pulse_bandwidth_ratio;
    amb_data_file>>x;tau_3db(i_beam)=Uvar(x,"us");   
    amb_data_file>>x;pulse_gate(i_beam)=Uvar(x,"m");
    amb_data_file>>x;lambda=Uvar(x,"cm");
    amb_data_file>>x;receive_window=Uvar(x,"us");
    amb_data_file>>x;bpd=Uvar(x,"us");
    amb_data_file>>N_p;
    amb_data_file>>Ns_rw;
    amb_data_file>>duty_cycle;
    amb_data_file>>x;rtt=Uvar(x,"ms");
    amb_data_file>>muhleman_k1;
    amb_data_file>>muhleman_k2;
    amb_data_file>>x;range_offset(i_beam)=Uvar(x,"m");
    amb_data_file>>x;frequency_offset(i_beam)=Uvar(x,"Hz");
    amb_data_file>>x;range_3db(i_beam)=Uvar(x,"km");
    amb_data_file>>x;X0(i_beam)=Uvar(x,"kg km km/(s s s)");
    amb_data_file>>x;x_res=Uvar(x,"km");
    amb_data_file>>x;rg_res=Uvar(x,"km");

    amb_data_file>>x;usable_area(i_beam)= Uvar(x,"km km");
    //------------------------------------------
    //Read Array Variables
    //-----------------------------------------
    for (unsigned int i =0 ; i < Nrange_patch_;i++)
      {       
      amb_data_file>>x; range_lower(i_beam,i)=Uvar(x,"km");
      amb_data_file>>x; range_upper(i_beam,i)=Uvar(x,"km");
      }
    for (unsigned int i =0 ; i < Ndop_patch_;i++)
      {
      amb_data_file>>x; dop_lower(i_beam,i)=Uvar(x,"Hz");
      amb_data_file>>x; dop_upper(i_beam,i)=Uvar(x,"Hz");
      }
      
    for (unsigned int i =0 ; i < Nrange_patch_;i++)
    for (unsigned int j =0 ; j < Ndop_patch_;j++)
      {
      amb_data_file>>x;avg_incidence(i_beam,i,j) = Uvar(x,"deg");
      amb_data_file>>avg_backscatter(i_beam,i,j);
      amb_data_file>>avg_antenna(i_beam,i,j);
      amb_data_file>>x;avg_azi(i_beam,i,j) = Uvar(x,"deg");
      amb_data_file>>x;avg_elev(i_beam,i,j) = Uvar(x,"deg");      
      amb_data_file>>x;patch_sum(i_beam,i,j)= Uvar(x,"kg km km/(s s s)");
      amb_data_file>>x;avg_area(i_beam,i,j) = Uvar(x,"km km");
      }
       
    for (unsigned int i =0 ; i < Nrange_patch_;i++)
    for (unsigned int j = 0; j < Ndop_patch_;j++)
    for (unsigned int i_range = 0; i_range<Nrange_bin_; i_range++)
    for (unsigned int j_dop =0; j_dop<Ndop_bin_; j_dop++)
      {
      unsigned int patch_index = i*Ndop_patch_+j;
      unsigned int bin_index = i_range*Ndop_bin_+j_dop;	
      amb_data_file>>x;lat(i_beam,patch_index,bin_index)=Uvar(x,"deg");
      amb_data_file>>x;lon(i_beam,patch_index,bin_index)=Uvar(x,"deg");
      amb_data_file>>x;thetai(i_beam,patch_index,bin_index)=Uvar(x,"deg");
      amb_data_file>>x;range(i_beam,patch_index,bin_index)=Uvar(x,"km");
      amb_data_file>>x;doppler(i_beam,patch_index,bin_index)=Uvar(x,"Hz");  
      amb_data_file>>x;alongtrack(i_beam,patch_index,bin_index)=Uvar(x,"km");  
      amb_data_file>>x;crosstrack(i_beam,patch_index,bin_index)=Uvar(x,"km");
      }	      
		        
    for (unsigned int i_range = 0; i_range<Nrange_bin_; i_range++)
    for (unsigned int j_dop =0; j_dop<Ndop_bin_; j_dop++)
      { 
      int i1;
      
      amb_data_file>>x;
      bin_power_main(i_beam,i_range,j_dop)=Uvar(x,"kg km km/(s s s)");
      
      amb_data_file>>x;
      bin_power_amb(i_beam,i_range,j_dop)=Uvar(x,"kg km km/(s s s)");
      if (bin_power_main(i_beam,i_range,j_dop) == Uvar(0.0,"kg km km/(s s s)")
	  ||bin_power_amb(i_beam,i_range,j_dop) == Uvar(0.0,"kg km km/(s s s)"))
	{
	bin_ratiodB(i_beam,i_range,j_dop) = -50.0;
	}
      else
	{
	bin_ratiodB(i_beam,i_range,j_dop) = 10.0*
	  log(bin_power_main(i_beam,i_range,j_dop)
	      /bin_power_amb(i_beam,i_range,j_dop))/log(10.0);
	}
      
      amb_data_file>>x;
      normalized_bin_power(i_beam,i_range,j_dop)
	= Uvar(x,"kg km km/(s s s)");
      amb_data_file>>x;
      center_area(i_beam,i_range,j_dop)=Uvar(x,"km km");
      area_process_window(i_beam)+=center_area(i_beam,i_range,j_dop);
      amb_data_file>>x;
      bin_antenna(i_beam,i_range,j_dop)= x;
      if( bin_antenna(i_beam,i_range,j_dop)== 0.0)
	{
	bin_antennadB(i_beam,i_range,j_dop)= -50.0;
	}
      else
	{
	bin_antennadB(i_beam,i_range,j_dop)= 10.0 * log(x)/log(10.0);
	}
      amb_data_file>>i1;is_good_bin_amb(i_beam,i_range,j_dop)=i1;	
      amb_data_file>>i1;is_good_bin_nesigma0(i_beam,i_range,j_dop)=i1;	
      amb_data_file>>i1;is_good_bin_totalSNR(i_beam,i_range,j_dop)=i1;	
      amb_data_file>>i1;is_good_bin_pn_to_amb(i_beam,i_range,j_dop) = i1;
      amb_data_file>>i1;//is_good_bin(i_beam,i_range,j_dop)=i1;	
      }
    }    //Loop over all 5 beams  

  //--------------------------------
  //re-evaluate good area 
  //-------------------------------
  is_good_bin = 0;
  usable_area = Uvar(0,"km km");
  for (unsigned int i_beam = 0; i_beam < 5; ++i_beam)
  for (unsigned int i_range = 0; i_range <Nrange_bin_;++i_range)
  for (unsigned int j_dop=0; j_dop < Ndop_bin_; ++j_dop)
    {
      if (is_good_bin_amb(i_beam,i_range,j_dop)==1 &&
	  is_good_bin_nesigma0(i_beam,i_range,j_dop)==1)
	{
	  is_good_bin(i_beam,i_range,j_dop) = 1;
	 usable_area(i_beam) += center_area(i_beam,i_range,j_dop);
	}
    }

  //-----------------------------------------------------------------------
  //rearrange center patch along and cross into Nrange_bin * Ndop_bin arrays
  //-----------------------------------------------------------------------
  for (unsigned int i_beam = 0; i_beam<5; ++i_beam)
    {
    unsigned int i_patch = (Nrange_patch_ * Ndop_patch_)/2;
    for (unsigned int i_range = 0; i_range < Nrange_bin_; ++i_range)
    for (unsigned int j_dop = 0; j_dop < Ndop_bin_ ; ++j_dop)
      {
      unsigned int bin_index = i_range * Ndop_bin_ + j_dop;
      center_patch_along(i_beam,i_range,j_dop)=
	                     alongtrack(i_beam,i_patch,bin_index);
      center_patch_cross(i_beam,i_range,j_dop)=
                             crosstrack(i_beam,i_patch,bin_index);
      }
    }

  read_data_ = true;
  }

//----------------------------------------------------
//write summary into a file
//----------------------------------------------------
void AmbiguityData::displayAmbiguitySummary(ofstream& outputsummary) throw(ErrorMessage)
  {

    //----------------------------------------------
    //output file format
    //----------------------------------------------
    outputsummary<<"Epoch time "<<epoch_time_<<endl;
    outputsummary<<"sar_imaging_time in min"<<sar_imaging_time_.getInUnits("min")<<endl;
    outputsummary<<"Nrange_patch "<< Nrange_patch_<<endl;
    outputsummary<< "Ndop_patch "<<Ndop_patch_<<endl;
    outputsummary<<"Nrange_bin "<<Nrange_bin_<<endl;
    outputsummary<<"Ndop_bin "<< Ndop_bin_<<endl;
    outputsummary<<"max beam gain "<<max_beam_gain<<endl;
    outputsummary<<"frequency  "<<frequency<<endl;
    outputsummary<<"tau_p  "<<tau_p<<endl;
    outputsummary<<"Pt  "<<Pt<<endl;
    outputsummary<<"Tsys  "<<Tsys<<endl;
    outputsummary<<"Br  "<<BR<<endl;
    outputsummary<<"Bn  "<<Bn<<endl;
    outputsummary<<"fadc  "<<fadc<<endl;  
    outputsummary<<"radius  "<<radius<<endl;
    outputsummary<<"Pn  "<<Pn<<endl;    
    outputsummary<<"prf  "<<prf<<endl;
    outputsummary<<"min bpd "<<min_bpd<<endl;
    outputsummary<<"pbw  "<<pbw<<endl;
    outputsummary<<"pulse bandwidth ratio  "<<pulse_bandwidth_ratio<<endl;
    outputsummary<<"tau_3db "<<tau_3db<<endl;   
    outputsummary<<"pulse_gate  "<<pulse_gate<<endl;
    outputsummary<<"lambda  "<<lambda<<endl;
    outputsummary<<"receive widnow  "<<receive_window<<endl;
    outputsummary<<"bpd  "<<bpd<<endl;
    outputsummary<<"N_p  "<<N_p<<endl;
    outputsummary<<"Ns_wr  "<<Ns_rw<<endl;
    outputsummary<<"duty_cycle  "<<duty_cycle<<endl;
    outputsummary<<"rtt  "<<rtt<<endl;
    outputsummary<<"muhleman constant k1 "<<muhleman_k1<<endl;
    outputsummary<<"muhelman constant k2 "<<muhleman_k2<<endl;
    outputsummary<<"range offset  "<<range_offset<<endl;
    outputsummary<<"frequency offset  "<<frequency_offset<<endl;
    outputsummary<<"range_3db  "<<range_3db<<endl;
    outputsummary<<"X0  "<<X0<<endl;

    outputsummary.precision(5);
    outputsummary<<"CASSINI RADAR"<<endl;

    for (unsigned int i_beam = 0; i_beam<5;++i_beam)
      {
      outputsummary<<"beam number "<<i_beam +1 <<endl;
      outputsummary<<"dop_low/upper(KHz) range_lower/upper(km)"<<endl;
      for (unsigned int i = 0 ; i < Nrange_patch_;i++)
      for (unsigned int j = 0; j < Ndop_patch_; j++)
	{
	int ipatch = i - int(Nrange_patch_/2);
	int jpatch = j - int(Ndop_patch_/2);
	outputsummary<<setw(3)<<ipatch;
	outputsummary<<setw(3)<<jpatch;
	outputsummary<<setw(10)<<dop_lower(i_beam,j).getInUnits("KHz");
	outputsummary<<setw(10)<<dop_upper(i_beam,j).getInUnits("KHz");
	outputsummary<<setw(10)<<range_lower(i_beam,i).getInUnits("km");
	outputsummary<<setw(10)<<range_upper(i_beam,i).getInUnits("km")<<endl;
	}
      outputsummary<<endl;
      }

    outputsummary<<"Patch incident angle,backscatter azi and elev, and antenna gain, area"<<endl;
    for (unsigned int i_beam = 0 ; i_beam < 5; ++i_beam)
      {
      outputsummary<<"beam number "<<i_beam +1<<endl;
      for (unsigned int i = 0 ; i < Nrange_patch_;i++)
      for (unsigned int j = 0; j < Ndop_patch_; j++)
	{ 
	int ipatch = i - int(Nrange_patch_/2);
	int jpatch = j - int(Ndop_patch_/2);  
	double gain_ratio= avg_antenna(i_beam,i,j)/max_beam_gain(i_beam);
	if(gain_ratio == 0.0 ) gain_ratio = 1e-10;
	double gain_dB = 10.0 * log(gain_ratio)/log(10.0); 
	if(avg_backscatter(i_beam,i,j) == 0.0) 
	  avg_backscatter(i_beam,i,j) = 1e-10;
	double backscatter_dB = 10.0 
	  * log(avg_backscatter(i_beam,i,j))/log(10.0) ; 
	outputsummary<<setw(3)<<ipatch;
	outputsummary<<setw(3)<<jpatch;
	outputsummary<<setw(10)<<avg_incidence(i_beam,i,j).getInUnits("deg");
	outputsummary<<setw(10)<<backscatter_dB;
	outputsummary<<setw(10)<<avg_azi(i_beam,i,j).getInUnits("deg");
	outputsummary<<setw(10)<<avg_elev(i_beam,i,j).getInUnits("deg");
	outputsummary<<setw(10)<<gain_dB;
	outputsummary<<setw(10)<<avg_area(i_beam,i,j).getInUnits("km km")<<endl;
	}
      outputsummary<<endl;
      }
    
    //----------------
    //show amb level
    //---------------

    outputsummary<<"amb level "<<endl;
    for (unsigned int i_beam = 0 ; i_beam < 5; ++i_beam)
      {
      outputsummary<<"beam number "<<i_beam +1<<endl;
      for (unsigned int i = 0 ; i < Nrange_patch_;i++)
      for (unsigned int j = 0; j < Ndop_patch_; j++)
	{
	int ipatch = i - int(Nrange_patch_/2);
	int jpatch = j - int(Ndop_patch_/2);  
	if(patch_sum(i_beam,i,j) == Uvar(0.0,"kg km km/(s s s)"))
	  patch_sum(i_beam,i,j) = 1e-10 * 
            patch_sum(i_beam,int(Nrange_patch_/2),int(Ndop_patch_/2));
	double amb_level_ratio = 
            patch_sum(i_beam,i,j).getInUnits("kg km km/(s s s)");
	amb_level_ratio /= 
            patch_sum(i_beam,int(Nrange_patch_/2),
		      int(Ndop_patch_/2)).getInUnits("kg km km/(s s s)");
	double amb_level_ratiodB = 10.0 * log(amb_level_ratio)/log(10.0);
	outputsummary<<setw(3)<<ipatch<<setw(3)
		     <<jpatch<<setw(10)<<amb_level_ratiodB<<endl;
	}
      outputsummary<<endl;
      }
    //calculate each beam's area for three different conditions
    //(1) amb > + 14 dB
    //(2) Pn > Pamb
    //(3) ne sigmao < -10 dB
    Uvec area_amb("amb_area",5);
    Uvec area_pn_pamb("pn_pamb_area",5);
    Uvec area_nesigma0("ne_sigma0",5);
    area_amb = Uvar(0,"km km");
    area_pn_pamb = Uvar(0,"km km");
    area_nesigma0 = Uvar(0,"km km");

    for (unsigned int i_beam = 0 ; i_beam < 5; ++i_beam)
    for (unsigned int i_range=0; i_range<Nrange_bin_-1;++i_range)
    for (unsigned int j_dop=0; j_dop < Ndop_bin_-1;++j_dop)
      {
      if (is_good_bin_amb(i_beam,i_range,j_dop)==1) 
	area_amb(i_beam) += center_area(i_beam,i_range,j_dop);	
      if (is_good_bin_nesigma0(i_beam,i_range,j_dop)==1)
	area_nesigma0(i_beam)+= center_area(i_beam,i_range,j_dop);	
      if(is_good_bin_pn_to_amb(i_beam,i_range,j_dop)==1)
	area_pn_pamb(i_beam) +=center_area(i_beam,i_range,j_dop);      	
      }

    for (unsigned int i_beam = 0 ; i_beam < 5; ++i_beam)
      {
	//currently usable area does not satisfy Pn > Pamb
      outputsummary<<"beam "<<i_beam+1<<" usable area w/o pn> pamb "
		   <<usable_area(i_beam)<<endl;
      outputsummary<<"beam "<<i_beam+1<<" good amb ratio "
		   <<area_amb(i_beam)<<endl;
      outputsummary<<"beam "<<i_beam+1<<" area with pn> pamb "
		   <<area_pn_pamb(i_beam)<<endl;
      outputsummary<<"beam "<<i_beam+1<<" area with nesigma0 < -10 db "
		   <<area_nesigma0(i_beam)<<endl;
      outputsummary<<"beam "<<i_beam+1<<"process window area "
		   <<area_process_window(i_beam)<<endl;
      }
    outputsummary<<"total usable area "<<total_usable_area()<<endl;
    outputsummary<<"total usable area w/o Pn>Pamb condition"<<
      usable_area_wo_Pn_Pamb_requirement()<<endl;
  }


//----------------------------------------------------
//return number of range patches used for amb calculation
//----------------------------------------------------
unsigned int AmbiguityData::getNrangepatch() throw(ErrorMessage)
  {
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(Nrange_patch_);
  } 

//----------------------------------------------------
//return number of dop patches used for amb calculation
//----------------------------------------------------
unsigned int AmbiguityData::getNdoppatch() throw (ErrorMessage)
 {
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(Ndop_patch_);
  } 

//----------------------------------------------------
//return number of range bins inside processing window
//----------------------------------------------------
unsigned int AmbiguityData::getNrangebin() throw(ErrorMessage)
 {
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(Nrange_bin_);
  } 

//----------------------------------------------------
//return number of dop bins inside processing window
//----------------------------------------------------
unsigned int AmbiguityData::getNdopbin() throw (ErrorMessage)
 {
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(Ndop_bin_);
  } 

//----------------------------------------------------
//return epoch time
//----------------------------------------------------
string AmbiguityData::getEpochtime() throw(ErrorMessage)
 {
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(epoch_time_);
  } 

//----------------------------------------------------
//return sar imaging time w.r.t. epoch time
//----------------------------------------------------
Uvar AmbiguityData::getSarimagingtime() throw(ErrorMessage)
{
  if (read_data_ == false)
    {
    throw ErrorMessage("Data have not been read in");
    }
  return(sar_imaging_time_);
  } 


//----------------------------------------------------
//return lat and lon values inside beam's processing window
//----------------------------------------------------
void AmbiguityData::BeamSwathLonLat(const unsigned int& i_beam, Uvec& lon_track, Uvec& lat_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    unsigned int bin_index = i * Ndop_bin_ + j;
    lon_track(bin_index) = lon(i_beam,i_patch_center_,bin_index);
    lat_track(bin_index) = lat(i_beam,i_patch_center_,bin_index);
    }
  }

//----------------------------------------------------
//return cross and along track values inside beam's processing window
//----------------------------------------------------
void AmbiguityData::BeamSwathCrossAlong(const unsigned int& i_beam, Uvec& cross_track, Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    unsigned int bin_index = i * Ndop_bin_ + j;
    along_track(bin_index) = alongtrack(i_beam,i_patch_center_,bin_index);
    cross_track(bin_index) = crosstrack(i_beam,i_patch_center_,bin_index);
    }
  }

//----------------------------------------
//Return cross and along track values of patch number i for beam i_beam
//----------------------------------------
void AmbiguityData::CrossAlongPatches(const unsigned int& i_beam, 
			 const unsigned int& i_patch, 
			 Uvec& cross_track,
			 Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i=0 ;i <Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    unsigned int bin_index = i * Ndop_bin_ + j;
    cross_track(bin_index) = crosstrack(i_beam,i_patch,bin_index);
    along_track(bin_index) = alongtrack(i_beam,i_patch,bin_index);
    }
  }


//-------------------------------------
//return signal vs crosstrack
//------------------------------------
 void AmbiguityData::SignalvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& return_power,
			  Uvec& cross_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (alongtrack_index >= Ndop_bin_)
      {
      throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
      }

    for (unsigned int i = 0; i < Nrange_bin_;++i)
      {
      return_power(i) = bin_power_main(i_beam,i,alongtrack_index);
      cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
      }

  }


//-------------------------------------
//return ambiguity power vs crosstrack
//------------------------------------
 void AmbiguityData::AmbiguityvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& return_power,
			  Uvec& cross_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (alongtrack_index >= Ndop_bin_)
      {
      throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
      }

    for (unsigned int i = 0; i < Nrange_bin_;++i)
      {
      return_power(i) = bin_power_amb(i_beam,i,alongtrack_index);
      cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
      }

  }

//-------------------------------------
//return amb ratio vs crosstrack
//------------------------------------
void AmbiguityData::AmbratiodBvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& ratio_dB,
			  Uvec& cross_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (alongtrack_index >= Ndop_bin_)
      {
      throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
      }

    for (unsigned int i = 0; i < Nrange_bin_;++i)
      {
      ratio_dB(i) = bin_ratiodB(i_beam,i,alongtrack_index);
      cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
      }
  }

//--------------------------------
//return thermal SNR dB
//------------------------------
void AmbiguityData::ThermalSNRdBvsCrosstrack(const unsigned int& i_beam,
				const unsigned int& alongtrack_index,
				Uvec& ratio_dB,
				Uvec& cross_track) throw (ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (alongtrack_index >= Ndop_bin_)
      {
      throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
      }

    for (unsigned int i = 0; i < Nrange_bin_-1;++i)
      {
      double ratio= bin_power_main(i_beam,i,alongtrack_index).getInUnits("kg km km/(s s s)")/ Pn.getInUnits("kg km km/(s s s)");
      ratio *=(x_res*rg_res/center_area(i_beam,i,alongtrack_index)).getValue();
      if (ratio == 0.0)
	{
	  ratio_dB(i) = Uvar(-100," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      else
	{
	  ratio_dB(i) =10.0 * log(ratio)/log(10.0);
	  if (ratio_dB(i)<-100) ratio_dB(i) = -100.0;
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      }
  }
//-------------------------------------
//return antenna dB  vs crosstrack
//------------------------------------
void AmbiguityData::AntennadBvsCrosstrack(const unsigned int& i_beam,
			  const unsigned int& alongtrack_index,
			  Uvec& antenna_dB,
			  Uvec& cross_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (alongtrack_index >= Ndop_bin_)
      {
      throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
      }

    for (unsigned int i = 0; i < Nrange_bin_;++i)
      {
      antenna_dB(i) = bin_antennadB(i_beam,i,alongtrack_index);
      cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
      }
  }


//-------------------------------------
//return signal vs alongtrack
//------------------------------------
 void AmbiguityData::SignalvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& return_power,
			  Uvec& along_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (crosstrack_index >= Nrange_bin_)
      {
      throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
      }

    for (unsigned int i = 0; i < Ndop_bin_;++i)
      {
      return_power(i) = bin_power_main(i_beam,crosstrack_index,i);
      along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
      }
  }

//-------------------------------------
//return ambiguity power  vs alongtrack
//------------------------------------
 void AmbiguityData::AmbiguityvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& return_power,
			  Uvec& along_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (crosstrack_index >= Nrange_bin_)
      {
      throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
      }

    for (unsigned int i = 0; i < Ndop_bin_;++i)
      {
      return_power(i) = bin_power_amb(i_beam,crosstrack_index,i);
      along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
      }
  }


//-----------------------------------
//return amb ratio vs along track
//----------------------------------
void AmbiguityData::AmbratiodBvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& ratio_dB,
			  Uvec& along_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (crosstrack_index >= Nrange_bin_)
      {
      throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
      }
    for (unsigned int i = 0; i < Ndop_bin_;++i)
      {
      ratio_dB(i) = bin_ratiodB(i_beam,crosstrack_index,i);
      along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
      }

  }


//--------------------
//return thermal snr
//---------------------

void AmbiguityData::ThermalSNRdBvsAlongtrack(const unsigned int& i_beam,
				const unsigned int& crosstrack_index,
				Uvec& ratio_dB,
				Uvec& along_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (crosstrack_index >= Nrange_bin_)
      {
      throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
      }

    for (unsigned int i = 0; i < Ndop_bin_-1;++i)
      {
      double ratio=bin_power_main(i_beam,crosstrack_index,i).getInUnits("kg km km/(s s s)")/ Pn.getInUnits("kg km km/(s s s)");;
      ratio *=(x_res*rg_res/center_area(i_beam,crosstrack_index,i)).getValue();
      if (ratio == 0.0)
	{
	  ratio_dB(i) = Uvar(-100," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      else
	{
	  ratio_dB(i) =10.0 * log(ratio)/log(10.0);
	  if (ratio_dB(i) < -100) ratio_dB(i) = -100.0;
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      }
  }

//-----------------------------------
//return antenna dB vs along track
//----------------------------------
void AmbiguityData::AntennadBvsAlongtrack(const unsigned int& i_beam,
			  const unsigned int& crosstrack_index,
			  Uvec& antenna_dB,
			  Uvec& along_track) throw(ErrorMessage)
  {
    if(read_data_ == false)
      {
      throw ErrorMessage("Ambiguity data have not been read in");
      }
    if (crosstrack_index >= Nrange_bin_)
      {
      throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
      }
    for (unsigned int i = 0; i < Ndop_bin_;++i)
      {
      antenna_dB(i) = bin_antennadB(i_beam,crosstrack_index,i);
      along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
      }

  }






//-------------------------------------------
//Return locations of good bins
//-------------------------------------------
void AmbiguityData::Goodbins(const unsigned int& i_beam,Uvec& cross_track, Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    if (is_good_bin(i_beam,i,j) == 1)
      {
      unsigned int bin_index = i * Ndop_bin_ + j;
      along_track(bin_index) = alongtrack(i_beam,i_patch_center_,bin_index);
      cross_track(bin_index) = crosstrack(i_beam,i_patch_center_,bin_index);
      }
    }
  }

//-------------------------------------------
//Return locations of bins whose amb ration is > +14 dB
//-------------------------------------------
 void AmbiguityData::Goodbins_Signal_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    if (is_good_bin_amb(i_beam,i,j) == 1)
      {
      unsigned int bin_index = i * Ndop_bin_ + j;
      along_track(bin_index) = alongtrack(i_beam,i_patch_center_,bin_index);
      cross_track(bin_index) = crosstrack(i_beam,i_patch_center_,bin_index);
      }
    }
  }

//-------------------------------------------
//Return locations of bins whose Pn is larger than Pamb
//-------------------------------------------
 void AmbiguityData::Goodbins_Pn_to_Pamb(const unsigned int& i_beam, Uvec& cross_track,
			   Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    if (is_good_bin_pn_to_amb(i_beam,i,j) == 1)
      {
      unsigned int bin_index = i * Ndop_bin_ + j;
      along_track(bin_index) = alongtrack(i_beam,i_patch_center_,bin_index);
      cross_track(bin_index) = crosstrack(i_beam,i_patch_center_,bin_index);
      }
    }
  }



//-------------------------------------------
//Return locations of bins whose noise equivalent sigma0 is < -10 dB
//-------------------------------------------
 void AmbiguityData::Goodbins_Nesigma0(const unsigned int& i_beam, 
				       Uvec& cross_track,
				       Uvec& along_track) 
   throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }
  for (unsigned int i = 0; i < Nrange_bin_; ++i)
  for (unsigned int j = 0; j < Ndop_bin_; ++j)
    {
    if (is_good_bin_nesigma0(i_beam,i,j) == 1)
      {
      unsigned int bin_index = i * Ndop_bin_ + j;
      along_track(bin_index) = alongtrack(i_beam,i_patch_center_,bin_index);
      cross_track(bin_index) = crosstrack(i_beam,i_patch_center_,bin_index);
      }
    }
  }




//----------------------------------------------
//Return along and cross track values along the edges
// of beam process window
//------------------------------------------------
void AmbiguityData::Processwindowedges(const unsigned int& i_beam,Uvec& beam_edges_cross,Uvec& beam_edges_along) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }
  Uvec x_cross("x_cross",Nrange_bin_*Ndop_bin_);
  Uvec y_along("y_along",Nrange_bin_*Ndop_bin_);
  unsigned int k_count = 0;

  unsigned int i_range_lower = 0;
  unsigned int i_range_upper = Nrange_bin_-1; 
  unsigned int j_dop_lower = 0;
  unsigned int j_dop_upper = Ndop_bin_ - 1;

  for (unsigned int j = 0 ; j < Ndop_bin_; j++)
    {  
      x_cross(k_count) = center_patch_cross(i_beam,i_range_lower,j);
      y_along(k_count) = center_patch_along(i_beam,i_range_lower,j);
      k_count = k_count + 1;
    }

  for (unsigned int i = 0; i <Nrange_bin_;i++)
    {
      x_cross(k_count) = center_patch_cross(i_beam,i,j_dop_upper);
      y_along(k_count)= center_patch_along(i_beam,i,j_dop_upper);
      k_count = k_count + 1;
    }

  for (unsigned int j = 0 ; j < Ndop_bin_; j++)
    { 
      unsigned int j1 = Ndop_bin_-1 - j;
      x_cross(k_count) = center_patch_cross(i_beam,i_range_upper,j1);
      y_along(k_count)= center_patch_along(i_beam,i_range_upper,j1);
      k_count = k_count + 1;
    }

  for (unsigned int i = 0; i <Nrange_bin_;i++)
    {
      unsigned int i1 = Nrange_bin_ - 1 - i;
      x_cross(k_count) = center_patch_cross(i_beam,i1,j_dop_lower);
      y_along(k_count)= center_patch_along(i_beam,i1,j_dop_lower);
      k_count = k_count + 1;
    }
 
  beam_edges_cross.resize(k_count);
  beam_edges_along.resize(k_count);
  for (unsigned int i =0; i < k_count; ++i)
    {
    beam_edges_along(i) = y_along(i);
    beam_edges_cross(i) = x_cross(i);
    }
  } 

//-----------------------------------------------------------------
//Return alongtrack information whether it is a good bin or not?
//---------------------------------------------------------------
void AmbiguityData::Is_good_bin_Alongtrack(
			    const unsigned int& i_beam,
			    const unsigned int& i_indicator,
			    const unsigned int& crosstrack_index, 
			    Uvec& bin_indicator,
			    Uvec& along_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }
  if (crosstrack_index >= Nrange_bin_)
    {
    throw ErrorMessage("crosstrack index should be less than"+toStr(Nrange_bin_));
    }
  if (i_indicator < 1 || i_indicator > 5)
    throw ErrorMessage("there are only five options 1..5 ");
  bin_indicator=Uvar(0," ");
  switch(i_indicator)
    {
    case 1:
      for (unsigned int i = 0; i < Ndop_bin_;++i)
	{
	  bin_indicator=Uvar(is_good_bin_amb(i_beam,crosstrack_index,i)," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      break;
    case 2:
      for (unsigned int i = 0; i < Ndop_bin_;++i)
	{
	  bin_indicator
	    =Uvar(is_good_bin_nesigma0(i_beam,crosstrack_index,i)," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      break;
    case 3:
      for (unsigned int i = 0; i < Ndop_bin_;++i)
	{
	  bin_indicator
	    =Uvar(is_good_bin_totalSNR(i_beam,crosstrack_index,i)," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      break;
    case 4:
      for (unsigned int i = 0; i < Ndop_bin_;++i)
	{
	  bin_indicator
	    =Uvar(is_good_bin_pn_to_amb(i_beam,crosstrack_index,i)," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      break;
    case 5:
      for (unsigned int i = 0; i < Ndop_bin_;++i)
	{
	  bin_indicator=Uvar(is_good_bin_amb(i_beam,crosstrack_index,i)," ");
	  along_track(i)= center_patch_along(i_beam,crosstrack_index,i);
	}
      break;
    }
  }

//-----------------------------------------------------------------
//Return crosstrack information whether it is a good bin or not?
//---------------------------------------------------------------
void AmbiguityData::Is_good_bin_Crosstrack(
			    const unsigned int& i_beam,
			    const unsigned int& i_indicator,
			    const unsigned int& alongtrack_index, 
			    Uvec& bin_indicator,
			    Uvec& cross_track) throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }
  if (alongtrack_index >= Ndop_bin_)
    {
    throw ErrorMessage("alongtrack index should be less than"+toStr(Ndop_bin_));
    }
  if (i_indicator < 1 || i_indicator > 5)
    throw ErrorMessage("there are only five options 1..5 ");
  bin_indicator = Uvar(0," ");
  switch(i_indicator)
    {
    case 1:
      for (unsigned int i = 0; i < Nrange_bin_;++i)
	{
	  bin_indicator(i)
	    =Uvar(is_good_bin_amb(i_beam,i,alongtrack_index)," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      break;
    case 2:
      for (unsigned int i = 0; i < Nrange_bin_;++i)
	{
	  bin_indicator(i)
	    =Uvar(is_good_bin_nesigma0(i_beam,i,alongtrack_index)," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      break;
    case 3:
      for (unsigned int i = 0; i < Nrange_bin_;++i)
	{
	  bin_indicator(i)
	    =Uvar(is_good_bin_totalSNR(i_beam,i,alongtrack_index)," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      break;
    case 4:
      for (unsigned int i = 0; i < Nrange_bin_;++i)
	{
	  bin_indicator(i)
	    =Uvar(is_good_bin_pn_to_amb(i_beam,i,alongtrack_index)," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      break;
    case 5:
      for (unsigned int i = 0; i < Nrange_bin_;++i)
	{
	  bin_indicator(i)=Uvar(is_good_bin(i_beam,i,alongtrack_index)," ");
	  cross_track(i)= center_patch_cross(i_beam,i,alongtrack_index);
	}
      break;
    }
  }

//----------------------------------------------------
//calculate total usable area, excluding good areas covered by
//any two  beams
//----------------------------------------------------
Uvar AmbiguityData::total_usable_area() throw(ErrorMessage)
  {
  if(read_data_ == false)
    {
    throw ErrorMessage("Ambiguity data have not been read in");
    }

  //---------------
  //Calculate usable area covered by all 5 beams
  //,excluding overlapping regions
  //(1)count all usable area of beam 1
  //(2) for beam i > 1, check whether good area covered by beam i is also
  //  a good area covered by beam i-1.  If so, do not count it
  //--------------
  
  //calculate range and dop centers of each beam
  for (unsigned int i_beam = 0; i_beam < 5; ++i_beam)
    {
    range_center(i_beam) =(range_lower(i_beam,i_range_center_) 
			     + range_upper(i_beam,i_range_center_))/2.0;
    dop_center(i_beam) = (dop_lower(i_beam,j_dop_center_)
			    + dop_upper(i_beam,j_dop_center_))/2.0;
    }

  Uvar total_area("total_area",0,"km km");
  for (unsigned int i_range = 0; i_range < Nrange_bin_ ; ++i_range)
  for (unsigned int j_dop =0; j_dop < Ndop_bin_; ++j_dop)
    {
    if (is_good_bin(0,i_range,j_dop)==1)
      total_area += center_area(0,i_range,j_dop); 
    }

  for (unsigned int i_beam = 1; i_beam < 5; ++i_beam)
  for (unsigned int i_range = 0; i_range < Nrange_bin_ ; ++i_range)
  for (unsigned int j_dop =0; j_dop < Ndop_bin_; ++j_dop)
    {
    if (is_good_bin(i_beam,i_range,j_dop)==1)
      {
      Uvar range = range_center(i_beam) -0.5*pulse_gate(i_beam);
      range += double(i_range)/double(Nrange_bin_)*pulse_gate(i_beam);
      Uvar dop = dop_center(i_beam)-0.5*pbw;
      dop += double(j_dop)/double(Ndop_bin_)*pbw;

      //check whether range and dop values are outside beam-1 process window
      if (range <= range_upper(i_beam-1,i_range_center_)
	  && range >= range_lower(i_beam-1,i_range_center_) 
	  && dop <= dop_upper(i_beam-1,j_dop_center_)
	  && dop >= dop_lower(i_beam-1,j_dop_center_))
	{
	//find (range,dop) index in beam i_beam -1
	Uvar range_index;
	Uvar dop_index ;
	range_index= range-range_center(i_beam-1) + 0.5* pulse_gate(i_beam-1);
	range_index =range_index /  pulse_gate(i_beam-1);//0..1
	dop_index = dop - dop_center(i_beam-1) + 0.5*pbw;
	dop_index = dop_index / pbw; //0..1

	unsigned int i_range_beam_1 =int(
		   double(Nrange_bin_-1)*range_index.getValue());
	unsigned int j_dop_beam_1 =int(
		   double(Ndop_bin_-1)*dop_index.getValue());

	if (i_range_beam_1 >= Nrange_bin_ || j_dop_beam_1 >= Ndop_bin_)
	  {
	    cout<<"Warning:: i_range_beam_1 or j_dop_beam_1 is out of range"<<endl;
	    cout<<"range index and dop index "<<endl;
	    cout<<i_range_beam_1<<j_dop_beam_1<<endl;
	    cout<<"range range center dop dop center "<<endl;
	    cout<<range<<range_center(i_beam-1)
		<<dop<<dop_center(i_beam-1)<<endl;
	    cout<<"pulse gate and pbw "<<endl;
	    cout<<pulse_gate(i_beam-1)<<pbw<<endl;
	    cout<<"dop difference"<<dop - dop_center(i_beam-1) + 0.5*pbw<<endl;
	    cout<<"diff/pbw "<<
	      (dop - dop_center(i_beam-1) + 0.5*pbw)/pbw<<endl;
	    cout<<"dop boundary "<< dop_upper(i_beam-1,j_dop_center_)<<
	      dop_lower(i_beam-1,j_dop_center_)<<endl;
	  }

	if (is_good_bin(i_beam-1,i_range_beam_1,j_dop_beam_1) == 0)
	  total_area += center_area(i_beam,i_range,j_dop);
	}
      else
	{
	total_area += center_area(i_beam,i_range,j_dop);
	//good area outside beam i_beam -1
	}
      }
    }

  return(total_area);
  }


//----------------------------------------------------------
//Return usable area with one less restrict conditions: Pn is
// not a dominant noise source 
//----------------------------------------------------------
Uvar AmbiguityData::usable_area_wo_Pn_Pamb_requirement() throw (ErrorMessage)
  {
  if(read_data_ == false)
    {
      throw ErrorMessage("Ambiguity data have not been read in");
    }

  //---------------
  //Calculate usable area covered by all 5 beams
  //,excluding overlapping regions
  //(1)count all usable area of beam 1
  //(2) for beam i > 1, check whether good area covered by beam i is also
  //  a good area covered by beam i-1.  If so, do not count it
  //(3) here, we exclude a condition to satisfy usable_area_calculation
  //: Pn > Pamb
  // (4) To do (3), let's make a new index matrix for good bin indicator
  //--------------
  
  I3D new_good_bin("new_good_bin",5,Nrange_bin_,Ndop_bin_);
  new_good_bin = 0;
  for (unsigned int i_beam =0; i_beam <5; ++i_beam)
  for(unsigned int i_range = 0; i_range < Nrange_bin_;++i_range)
  for(unsigned int j_dop=0; j_dop < Ndop_bin_;++j_dop)
    {
      if (is_good_bin_amb(i_beam,i_range,j_dop) ==1 
	  && is_good_bin_nesigma0(i_beam,i_range,j_dop) ==1 )
	new_good_bin(i_beam,i_range,j_dop) =1;
    }

  //calculate range and dop centers of each beam
  for (unsigned int i_beam = 0; i_beam < 5; ++i_beam)
    {
    range_center(i_beam) =(range_lower(i_beam,i_range_center_) 
			     + range_upper(i_beam,i_range_center_))/2.0;
    dop_center(i_beam) = (dop_lower(i_beam,j_dop_center_)
			    + dop_upper(i_beam,j_dop_center_))/2.0;
    }

  Uvar total_area("total_area",0,"km km");
  for (unsigned int i_range = 0; i_range < Nrange_bin_ ; ++i_range)
  for (unsigned int j_dop =0; j_dop < Ndop_bin_; ++j_dop)
    {
    if (new_good_bin(0,i_range,j_dop)==1)
      total_area += center_area(0,i_range,j_dop); 
    }

  for (unsigned int i_beam = 1; i_beam < 5; ++i_beam)
  for (unsigned int i_range = 0; i_range < Nrange_bin_ ; ++i_range)
  for (unsigned int j_dop =0; j_dop < Ndop_bin_; ++j_dop)
    {
    if (new_good_bin(i_beam,i_range,j_dop)==1)
      {
      Uvar range = range_center(i_beam) -0.5*pulse_gate(i_beam);
      range += double(i_range)/double(Nrange_bin_)*pulse_gate(i_beam);
      Uvar dop = dop_center(i_beam)-0.5*pbw;
      dop += double(j_dop)/double(Ndop_bin_)*pbw;

      //check whether range and dop values are outside beam-1 process window
      if (range <= range_upper(i_beam-1,i_range_center_)
	  && range >= range_lower(i_beam-1,i_range_center_) 
	  && dop <= dop_upper(i_beam-1,j_dop_center_)
	  && dop >= dop_lower(i_beam-1,j_dop_center_))
	{
	//find (range,dop) index in beam i_beam -1
	Uvar range_index;
	Uvar dop_index ;
	range_index= range-range_center(i_beam-1) + 0.5* pulse_gate(i_beam-1);
	range_index =range_index /  pulse_gate(i_beam-1);//0..1
	dop_index = dop - dop_center(i_beam-1) + 0.5*pbw;
	dop_index = dop_index / pbw; //0..1

	unsigned int i_range_beam_1 =int(
		   double(Nrange_bin_-1)*range_index.getValue());
	unsigned int j_dop_beam_1 =int(
		   double(Ndop_bin_-1)*dop_index.getValue());

	if (i_range_beam_1 >= Nrange_bin_ || j_dop_beam_1 >= Ndop_bin_)
	  {
	    cout<<"Warning :: i_range_beam_1 or j_dop_beam_1 is out of range"<<endl;
	    cout<<"range index and dop index "<<endl;
	    cout<<i_range_beam_1<<j_dop_beam_1<<endl;
	    cout<<"range range center dop dop center "<<endl;
	    cout<<range<<range_center(i_beam-1)
		<<dop<<dop_center(i_beam-1)<<endl;
	    cout<<"pulse gate and pbw "<<endl;
	    cout<<pulse_gate(i_beam-1)<<pbw<<endl;
	    cout<<"dop difference"<<dop - dop_center(i_beam-1) + 0.5*pbw<<endl;
	    cout<<"diff/pbw "<<
	      (dop - dop_center(i_beam-1) + 0.5*pbw)/pbw<<endl;
	    cout<<"dop boundary "<< dop_upper(i_beam-1,j_dop_center_)<<
	      dop_lower(i_beam-1,j_dop_center_)<<endl;
	  }

	if (new_good_bin(i_beam-1,i_range_beam_1,j_dop_beam_1) == 0)
	  total_area += center_area(i_beam,i_range,j_dop);
	}
      else
	{
	total_area += center_area(i_beam,i_range,j_dop);
	//good area outside beam i_beam -1
	}
      }
    }
    return (total_area);
  }
