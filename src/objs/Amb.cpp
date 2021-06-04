//==========================================================
//
//
//Calculate amb for a given time t, state_vector (sc_state)
//==========================================================






//---------------
// Spice routines
//---------------

#include <SpiceUsr.h>

//---------------
// Other includes
//---------------

#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "Frame.h"
#include "TargetGeom.h"
#include "Units.h"
#include "Time.h"
#include "Amb.h"
#include "Constants.h"
#include "TemplateUtils.h"


//-------------------------------------------------------------
//Calculate ambiguity at t for Beam
// need to pass Beam information later
//-------------------------------------------------------------

void cal_ambiguity(const string& target_name,
		   Beam& beam,
		   const Frame& ftitan,
		   const StateVector& sc_state,
		   const Uvar& prf,
		   const Uvar& pulse_gate,
		   const Uvar& pbw,
		   const Uvar& lambda,
		   const Uvar& X0,
		   const Uvar& Pn,
		   const Uvar& x_res,
		   const Uvar& rg_res,
		   const double& Ni,
		   const double&  muhleman_k1,
		   const double&  muhleman_k2,
		   const unsigned int& Nrange_patch,
		   const unsigned int& Ndop_patch,
		   const unsigned int& Nrange_bin,
		   const unsigned int& Ndop_bin,
		   const Uvar& range_offset,
		   const Uvar& frequency_offset,
		   ofstream& outputfile) 
                  throw(Unit::UnitError,ErrorMessage)
  { 
    
    typedef Array2D<DirectionVector> Dirmat;    
    Uvar usable_area = Uvar(0);
    Time t = sc_state.time();
    unsigned int i_range_center = int (Nrange_patch/2);
    unsigned int j_dop_center = int(Ndop_patch/2);
    Uvar speed_light = Uvar("speed_light");       

    //get pri
    Uvar pri = 1.0/prf;

    //get max beam gain
    //double max_gain = beam.getMaxGain();
    Uvar azimuth_width = beam.getAzimuthWidthOneWay();
    Uvar elevation_width=beam.getElevationWidthOneWay();
    //------------------------------------------
    //First use Frame method to get range and dop
    //------------------------------------------
    unsigned int beamnum = beam.getBeamNumber();
    Frame fbeam("CASSINI_RADAR_" + toStr(beamnum),"Cassini");
    DirectionVector boresight("boresight",fbeam,t,0,0,1);
    Uvar range_F, dop_F, thetai_F,altitude_F;
    TargetGeom tg(t);
    tg.setState(sc_state);
    tg.setLookDirection(boresight);
    tg.setTarget(target_name,ftitan);
    Uvar radius = tg.radius();
    range_F = tg.range();
    dop_F = tg.doppler(lambda); 
    thetai_F = tg.incidenceAngle();
    altitude_F = tg.altitude();
    
    //-----------------------------------------
    //Then, use non-frame method(vector equations)
    // to calculate boresight range and dop
    //------------------------------------------
    
    PositionVector sc_pos = sc_state.position();
    sc_pos.representIn(ftitan);
    DirectionVector pos_dir = sc_pos;
    pos_dir.representIn(fbeam);   
    
    FloatVector velocity = sc_state.velocity();
    velocity.representIn(fbeam);
    DirectionVector velocity_direction = velocity;
   
    //---------
    //doppler
    //---------
    double v_b = dot(velocity_direction,boresight);
    Uvar dop0 = 2.0* v_b * velocity.magnitude()/lambda;
   
    //---------
    //range
    //---------
    double p_b = dot(pos_dir,boresight);
    Uvar in_sqrt = sc_pos.magnitude() *sc_pos.magnitude() * p_b * p_b;
    in_sqrt += radius * radius - sc_pos.magnitude()*sc_pos.magnitude();
    Uvar range0 = Uvar(0.0);
    if(in_sqrt.getValue() > 0.0) 
      {
	range0 = -sc_pos.magnitude() * p_b - sqrt(in_sqrt);
      }
    else
      {
	throw ErrorMessage("No intercepting point in the boresight direction");
	
      }  
   
  
    //-----------------------------------------
    //compare frame-method with non frame method
    //------------------------------------------
    double diff_range = range0.getInUnits("km") - range_F.getInUnits("km");
    double diff_dop   = dop0.getInUnits("Hz") - dop_F.getInUnits("Hz");
    
    if (fabs(diff_range)>0.001)
      {
	throw ErrorMessage("Frame method does not produce the same range");
      }
    if (fabs(diff_dop) > 0.001)
      {
	throw ErrorMessage("Frame method does not produce the same dop");
      }
    

    //----------------------------------------------------------
    //Determine range and dop centers for each patch
    //----------------------------------------------------------      
    Uvec range_center("range_center",Nrange_patch);
    Uvec dop_center("dop_center",Ndop_patch);
    Umat range_bin("range_bin",Nrange_patch,Nrange_bin);
    Umat dop_bin("dop_bin",Ndop_patch,Ndop_bin);
    range_center = Uvar(0.0);
    dop_center = Uvar(0.0);
    range_bin = Uvar(0.0);
    dop_bin = Uvar(0.0);
   

    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      range_center(i) = range0;
      range_center(i) += (double(i) -double(i_range_center))*speed_light*pri/2.0;
      range_center(i) -= range_offset;
      for (unsigned int i_range=0; i_range< Nrange_bin; i_range++)
      { 
	range_bin(i,i_range) = range_center(i) -0.5*pulse_gate;
	range_bin(i,i_range)+= double(i_range)/double(Nrange_bin-1)*pulse_gate;
      }
    }

    for (unsigned int j = 0; j < Ndop_patch; j++)
    {
      dop_center(j) = dop0 + (double(j)- double(j_dop_center))*prf;
      dop_center(j) -= frequency_offset;
      for (unsigned int j_dop=0; j_dop<Ndop_bin; j_dop++)
      {
	dop_bin(j,j_dop) = dop_center(j) -0.5 * pbw ;
	dop_bin(j,j_dop) += double(j_dop)/double(Ndop_bin-1)*pbw;
      }
    }

    //--------------------------------------------------------
    //Display range and dop settings on the screen
    //-------------------------------------------------------  
    /*
    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      cout<<"range center "<<range_center(i)<<endl ;
      if (i == i_range_center)
      {
	for (unsigned int i_range=0; i_range< Nrange_bin; i_range++)
	{ 
	  //cout<<"range bin "<<range_bin(i,i_range) <<endl;
	}
      }
    }

    for (unsigned int i = 0; i < Ndop_patch; i++)
    {
      cout<<"dop center "<<dop_center(i)<<endl;
      if (i == j_dop_center)
      {
	for (unsigned int i_dop=0; i_dop<Ndop_bin; i_dop++)
	{
	  //cout<<"dop_bin "<<dop_bin(i,i_dop)<<endl;
	}
      }
    }
    */



    Imat no_of_solution("no_solution",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat look_azimuth("look_azimuth_angle",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat look_elevation("look_elevation_angle",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  lat("bin_patch_lat",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  lon("bin_patch_lon",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  range("range",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  doppler("doppler",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  thetai("thetai",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat  alongtrack("alongtrack",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat crosstrack("crosstrack",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Umat bin_power_main("bin_power_main",Nrange_bin,Ndop_bin);
    Umat bin_power_amb("bin_power_amb",Nrange_bin,Ndop_bin);
    Umat normalized_bin_power("normalized_bin_power",Nrange_bin,Ndop_bin);
    Dmat  bin_backscatter("bin_backscatter",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dmat bin_antenna("bin_antenna",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dmat center_patch_antenna("center_bin_antenna",Nrange_bin,Ndop_bin);
    Umat patch_sum("patch_strength",Nrange_bin,Ndop_bin);
    Umat total_noise("total_nose_summed",Nrange_bin,Ndop_bin);
    Umat center_area("center_area",Nrange_bin,Ndop_bin);
    Umat area("area_element",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);
    Dirmat surface_intercept("surface_intercept",Nrange_patch*Ndop_patch,Nrange_bin*Ndop_bin);

   

    //Initialize variables
    no_of_solution= 0;
    look_azimuth=Uvar(0.0);
    look_elevation= Uvar(0.0);
    lat= Uvar(0.0);
    lon= Uvar(0.0);
    range= Uvar(0.0);
    doppler= Uvar(0.0);
    thetai= Uvar(0.0);
    alongtrack=Uvar(0.0);
    crosstrack=Uvar(0.0);
    bin_power_main= Uvar(0.0);
    bin_power_amb= Uvar(0.0);
    normalized_bin_power= Uvar(0.0);
    bin_backscatter= 0.0;
    bin_antenna= 0.0;
    center_patch_antenna= 0.0;
    patch_sum= Uvar(0.0);
    total_noise= Pn; // Pn = k * Tsys * Bn
    center_area= Uvar(0.0); 
    area= Uvar(0.0);

    //------------------------------------------------------------------
    //Find  intercepting points corresponding range and dop values
    //no_interception: store how many range and dop bins have a solution
    //------------------------------------------------------------------
  
    TargetGeom tg1(t);
    for (unsigned int i = 0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
	{
	unsigned int patch_index = i * Ndop_patch+j;   
	for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
	for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
	  {
	    unsigned int bin_index = i_range * Ndop_bin+j_dop;	
	    tg1.setState(sc_state);
	    tg1.setTarget(target_name,ftitan);
	    tg1.setRangeDoppler(range_bin(i,i_range), dop_bin(j,j_dop),lambda);
	    DirectionVector look = tg1.lookDirection();
	    if (tg1.foundSurfaceIntercept()==false)
	      {
	      no_of_solution(patch_index,bin_index) = 0;//no surface intercept point
	      }
	    else
	      {
	      no_of_solution(patch_index,bin_index) = 1;   
	      DirectionVector new_look= look.representIn(fbeam);
	      new_look.getAzimuthElevation(look_azimuth(patch_index,bin_index),
				     look_elevation(patch_index,bin_index));
	      lat(patch_index,bin_index)= tg1.lat();
	      lon(patch_index,bin_index) = tg1.lon();
	      range(patch_index,bin_index)= tg1.range();
	      doppler(patch_index,bin_index) =tg1.doppler(lambda);
	      thetai(patch_index,bin_index)= tg1.incidenceAngle(); 
	      tg1.interceptAlongCross(alongtrack(patch_index,bin_index),
				      crosstrack(patch_index,bin_index));
	      surface_intercept(patch_index,bin_index)=tg1.surfaceIntercept();
	      //DirectionVector mirror_look
	      } 
	    tg1.reset(t);
	  }//Loop over i_range and j_dop
	}//Loop over range_patch and dop_patch

    //    cout<<"Finish finding all surface intercept points"<<endl;
    //    cout<<"Start calculating area elements and radar return signal "<<endl;
    //------------------------------------------	
    //End of finding all range and dop solutions
    //if there is a surface intercept point, no_of_s
    //----------------------------------------------
    //****************************************
    //    Main Ambiguity Calculation
    //there is a solution for the given  range and dop value
    //range_bin(i,i_range) and dop_bin(j,j_dop)
    //let's try to get radar return signal strenth here
    //*******************************
       
    for (unsigned int i =0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
      {
      unsigned int patch_index = i * Ndop_patch+j;   
      for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
      for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
	{
	unsigned int bin_index = i_range * Ndop_bin+j_dop;	
	if(i_range < Nrange_bin -1 
	   && j_dop <  Ndop_bin -1 
	   && no_of_solution(patch_index,bin_index)==1)
	  { 
	  unsigned int bin_Up,  bin_Left,bin_Diag;
	  bin_Up = (i_range+1)*Ndop_bin+j_dop;	
	  bin_Left = i_range*Ndop_bin+(j_dop+1);
	  bin_Diag = (i_range+1)*Ndop_bin+j_dop+1;
	 
	  //------------------------------------
	  //find Area element covered by (range,dop)
	  //divide an area enclosed by two range and dop values
	  // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
	  //into two triangles and calclate their areas and add them
	  //together
	  //-------------------------------------
	  if (no_of_solution(patch_index,bin_Up) * 
	      no_of_solution(patch_index,bin_Diag)*
	      no_of_solution(patch_index,bin_Left) == 1)
	    {
	      //Here I will calculate surface area, beam , backscatter
	      //radar return.  Beware that I do not include 
	      //i_range = Nrang_bin-1 j_dop = Ndop_bin-1
	      // their contributions are counted when i_range= Nrange_bin-2
	      // and j_dop = Ndop_bin -2
	    Uvar side_a = radius * 
	      fabs((surface_intercept(patch_index,bin_index)
		    .angle(surface_intercept(patch_index,bin_Left))
		    .getInUnits("rad")));
	    Uvar side_b = radius * 
	      fabs((surface_intercept(patch_index,bin_index)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    Uvar side_c = radius * 
	      fabs((surface_intercept(patch_index,bin_Left)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    Uvar cos_c = (side_a*side_a + side_b*side_b - side_c*side_c)
	      /(2.0 * side_a*side_b);
	    Uvar angle_c = Uvar(acos(cos_c.getValue()),"rad");
	    
	    area(patch_index,bin_index) = 
	      0.5 * sin(angle_c) * side_a * side_b;
	    
	    Uvar side_d = radius * 
	      fabs((surface_intercept(patch_index,bin_Left)
		    .angle(surface_intercept(patch_index,bin_Diag))
		    .getInUnits("rad")));
	    Uvar side_e = radius * 
	      fabs((surface_intercept(patch_index,bin_Diag)
		    .angle(surface_intercept(patch_index,bin_Up))
		    .getInUnits("rad")));
	    cos_c = (side_d*side_d + side_e*side_e - side_c*side_c)
	      /(2.0 * side_d*side_e);
	    angle_c = Uvar(acos(cos_c.getValue()),"rad");
	    area(patch_index,bin_index) += 
	      0.5 * sin(angle_c) * side_d * side_e;	    
	    //--------
	    //averaged beam gain over bin_index, bin_Up,bin_Left,bin_Diag
	    //---------
	    
	    double beam_gain = beam.bilinear(
				      look_azimuth(patch_index,bin_index),
				    look_elevation(patch_index,bin_index));
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Up),
				    look_elevation(patch_index,bin_Up));		
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Left),
			            look_elevation(patch_index,bin_Left));
	    beam_gain +=beam.bilinear(look_azimuth(patch_index,bin_Diag),
				    look_elevation(patch_index,bin_Diag));
	    beam_gain /= 4.0;
	        
	    //--------------------
	    //averaged back scattering coeff
	    //--------------------
	    double sigma0_elem = muhleman_k1 
	      * cos(thetai(patch_index,bin_index))
	      / pow(sin(thetai(patch_index,bin_index))
		    +muhleman_k2*cos(thetai(patch_index,bin_index))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Up))
	      / pow(sin(thetai(patch_index,bin_Up))
		    +muhleman_k2*cos(thetai(patch_index,bin_Up))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Left))
	      / pow(sin(thetai(patch_index,bin_Left))
		    +muhleman_k2*cos(thetai(patch_index,bin_Left))
		    ,3.0);
	    sigma0_elem += muhleman_k1 
	      * cos(thetai(patch_index,bin_Diag))
	      / pow(sin(thetai(patch_index,bin_Diag))
		    +muhleman_k2*cos(thetai(patch_index,bin_Diag))
		    ,3.0);
	    sigma0_elem /= 4.0;
	   
	    
	    //-----------------------------
	    //store beam gain and sigm0
	    //----------------------------
	    bin_antenna(patch_index,bin_index) = beam_gain;
	    bin_backscatter(patch_index,bin_index) = sigma0_elem;
	    //---------------------------
	    //average range
	    //---------------------------
	    Uvar avg_range =(range_bin(i,i_range) +range_bin(i,i_range+1))/2.0;    
	    		
	    //---------------------------------------
	    //find X = G^2 * lambda^2 * sigma0 / R^4
	    //--------------------------------------
	    Uvar X = beam_gain * beam_gain * lambda * lambda* sigma0_elem
	      /pow(avg_range,4);  	
	  
	    if (i == i_range_center && j == j_dop_center)
	      {
	      bin_power_main(i_range,j_dop)=X*X0*area(patch_index,bin_index);  
	      normalized_bin_power(i_range,j_dop)= bin_power_main(i_range,j_dop);
	      normalized_bin_power(i_range,j_dop) /=sigma0_elem;  
	      patch_sum(i,j) += X*X0*area(patch_index,bin_index);  
	      center_patch_antenna(i_range,j_dop)= beam_gain;
	      center_area(i_range,j_dop)=area(patch_index,bin_index); 
	      
	      //Debugging routine
	      /*
	      cout<<"beam number range and dop "<<beamnum<<" "
		  <<i_range<<" "<<j_dop<<endl;
	      cout<<look_azimuth(patch_index,bin_index).getInUnits("deg")<<" "
	       <<look_elevation(patch_index,bin_index).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Up).getInUnits("deg")<<
	       " "<<look_elevation(patch_index,bin_Up).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Left).getInUnits("deg")<<" "
		<<look_elevation(patch_index,bin_Left).getInUnits("deg")<<endl;
	      cout<<look_azimuth(patch_index,bin_Diag).getInUnits("deg")<<" "
	       <<look_elevation(patch_index,bin_Diag).getInUnits("deg")<<endl;
	      cout<<"gain area sigma0 "<<10.0*log(beam_gain/max_gain)/log(10.0)
	       <<" " <<area(patch_index,bin_index)<<" "<<sigma0_elem<<endl<<endl;
	      */
	      // End of Debugging routine
	      }
	    else
	      {
	      bin_power_amb(i_range,j_dop) += X * X0 
		* area(patch_index,bin_index);
	      patch_sum(i,j) +=X * X0 
		* area(patch_index,bin_index);
	      }
	    }//if there is a solution at Up, Left, 
	     //Diag and (i_range != Nrange_bin -1) and (j_dop !=Ndob_bin-1)
	  }//if there is a solution   
	}//Loop over range/dop bin
      }//Loop over range/dop patch
      


   
      
    //cout<<"Starting counting usable area "<<endl;
   
    //-----------------------------------------------------------------
    //Determine whether each bin is good for imaging
    // criteria
    // (1) amb_ratio > 14 dB
    // (2) thermal noise equivalent sigma 0 < -10 dB
    // (3) total snr > thermal snr + 3db
    //-----------------------------------------------------------------
    Imat is_good_bin("is_this_a_good_bin",Nrange_bin,Ndop_bin);   
    Imat is_good_bin_amb("amb_ratio_requirements",Nrange_bin,Ndop_bin);   
    Imat is_good_bin_nesigma0("noise_equivalent_sigma0",Nrange_bin,Ndop_bin);  
    Imat is_good_bin_totalSNR("total_snr_to_theraml_snr",Nrange_bin,Ndop_bin);
    Imat is_good_bin_pn_to_amb("ration_of_pn_to_amb",Nrange_bin,Ndop_bin);
    is_good_bin = 0;
    is_good_bin_amb = 0;
    is_good_bin_nesigma0=0;
    is_good_bin_totalSNR = 0;
    is_good_bin_pn_to_amb = 0;

    for (unsigned int i_range = 0; i_range < Nrange_bin-1; i_range++)
    for (unsigned int j_dop = 0 ; j_dop < Ndop_bin-1; j_dop++)
    { 
      double bin_ratio;
      double bin_ratiodB;
      unsigned int patch_index = i_range_center * Ndop_patch+j_dop_center;
      unsigned int bin_index = i_range * Ndop_bin+j_dop;	

      if (no_of_solution(patch_index,bin_index) == 1)
	{
	  if (bin_power_main(i_range,j_dop)==Uvar(0.0))
	  {
	  bin_ratio = 1e-100;
	  }
	else
	  {
	  bin_ratio = 
	    bin_power_main(i_range,j_dop).getInUnits("kg m m/(s s s)");
	  bin_ratio /= 
	    bin_power_amb(i_range,j_dop).getInUnits("kg m m/(s s s)");	
	  }
	}
      else if (no_of_solution(patch_index,bin_index) == 0)
	{
	bin_ratio = 1e-100;	
	}
      else
	{
	throw ErrorMessage("no_solution data set has been contaminated");
	}
      bin_ratiodB = 10.0*  log(bin_ratio)/log(10.0);

      Uvar power_sigma0_normalized =  normalized_bin_power(i_range,j_dop) ;  
      power_sigma0_normalized *= x_res * rg_res /area(patch_index,bin_index);

      double ne0 = Pn.getInUnits("kg m m/(s s s)");
      ne0 /= power_sigma0_normalized.getInUnits("kg m m/(s s s)");
      double ne0_dB = 10.0 * log(ne0)/log(10.0);
      
      Uvar main_power_normalized= bin_power_main(i_range,j_dop);
      main_power_normalized *= x_res * rg_res /area(patch_index,bin_index);     
      double thermal_snr = main_power_normalized.getInUnits("kg m m/(s s s)")/
      	Pn.getInUnits("kg m m/(s s s)");
      double thermal_snrdB = 10.0 * log(thermal_snr)/log(10.0);

      Uvar amb_power_normalized = bin_power_amb(i_range,j_dop);
      amb_power_normalized *= x_res * rg_res /area(patch_index,bin_index);

      double total_snr = main_power_normalized.getInUnits("kg m m/(s s s)");
      total_snr /=(amb_power_normalized.getInUnits("kg m m/(s s s)")
		   +Pn.getInUnits("kg m m/(s s s)"));
      double total_snrdB = 10.0 * log(total_snr)/log(10.0);
      double pn_to_amb = Pn.getInUnits("kg m m/(s s s)")
                /amb_power_normalized.getInUnits("kg m m/(s s s)");

      if (bin_ratiodB >= 14.0) is_good_bin_amb(i_range,j_dop) = 1;
      if (ne0_dB <= -10.0)    is_good_bin_nesigma0(i_range,j_dop)=1;
      if ((total_snrdB + 3.0)>= thermal_snrdB) 
               is_good_bin_totalSNR(i_range,j_dop) = 1;
      if (pn_to_amb >= 1.0) is_good_bin_pn_to_amb(i_range,j_dop) = 1;

      if (is_good_bin_amb(i_range,j_dop)==1
	  && is_good_bin_nesigma0(i_range,j_dop)==1
	  && is_good_bin_totalSNR(i_range,j_dop)==1)
      {
	usable_area += area(patch_index,bin_index);
	is_good_bin(i_range,j_dop) = 1;
      }
    }

    
    //---------------------------------------------------------
    //Generate range and dop range for each patch for outputfile
    //---------------------------------------------------------
    Uvec range_lower("range_lower",Nrange_patch);
    Uvec range_upper("range_upper",Nrange_patch);
    Uvec dop_lower("dop_lower",Ndop_patch);
    Uvec dop_upper("dop_upper",Ndop_patch);

 
    for (unsigned int i = 0 ; i <Nrange_patch; i++)
      {
      range_lower(i) = range_center(i) -0.5*pulse_gate;
      range_upper(i)=range_center(i)+0.5*pulse_gate;
      }
    for (unsigned int j = 0; j < Ndop_patch; j++)
      {      
      dop_lower(j) = dop_center(j) -0.5 * pbw ;
      dop_upper(j) = dop_center(j) + 0.5*pbw;      
      }

    //---------------------------------------------------------
    //Find average incident angle, antenna_gain, backscatter
    // azi
    //---------------------------------------------------------
    Umat avg_incidence("avg_incidence",Nrange_patch,Ndop_patch);
    Dmat avg_backscatter("avg_backsctter",Nrange_patch,Ndop_patch);
    Dmat avg_antenna("avg_antenna_gain",Nrange_patch,Ndop_patch);
    Umat avg_azi("avg_azi_angle",Nrange_patch,Ndop_patch);
    Umat avg_elev("avg_elev",Nrange_patch,Ndop_patch);
    Umat avg_area("avg_area",Nrange_patch,Ndop_patch);

    avg_incidence = Uvar(0.0);
    avg_backscatter = 0.0;
    avg_antenna = 0.0;
    avg_azi = Uvar(0.0);
    avg_elev = Uvar(0.0);
    avg_area = Uvar(0.0);

    for (unsigned int i =0 ; i < Nrange_patch;i++)
    for (unsigned int j = 0; j < Ndop_patch;j++)
      {
      unsigned int patch_index = i*Ndop_patch+j;
      unsigned int bin_count = 0;
      for (unsigned int i_range = 0; i_range<Nrange_bin-1; i_range++)
      for (unsigned int j_dop =0; j_dop<Ndop_bin-1; j_dop++)
	{
        unsigned int bin_index = i_range*Ndop_bin+j_dop;	
	if (no_of_solution(patch_index,bin_index) == 1)
	  {
	  avg_incidence(i,j) += thetai(patch_index,bin_index);
	  avg_backscatter(i,j) += bin_backscatter(patch_index,bin_index);
	  avg_antenna(i,j) +=bin_antenna(patch_index,bin_index);
	  avg_azi(i,j) += look_azimuth(patch_index,bin_index);
	  avg_elev(i,j) += look_elevation(patch_index,bin_index);
	  avg_area(i,j) += area(patch_index,bin_index);
	  bin_count ++;
	  }
	}//Loop over range/dop_bin
      if (bin_count != 0)
	{
	avg_incidence(i,j) /= double(bin_count);
	avg_backscatter(i,j) /=double(bin_count);
	avg_antenna(i,j)/=double(bin_count);
	avg_azi(i,j) /= double(bin_count);
	avg_elev(i,j) /=double(bin_count);
	avg_area(i,j) /= double(bin_count);
	}
      }//Loop over range/dop_patch
  


    //------------------------------------------
    //Output data file format
    //Let's write down the output to the file
    //-------------------------------------------   
    outputfile<<usable_area.getInUnits("km km")<<endl;

    for (unsigned int i = 0 ; i <Nrange_patch; i++)
    {
      outputfile<<range_lower(i).getInUnits("km")<<endl; 
      outputfile<<range_upper(i).getInUnits("km")<<endl;
    }

    for (unsigned int j = 0; j < Ndop_patch; j++)
    {      
      outputfile<<dop_lower(j).getInUnits("Hz")<<endl;
      outputfile<<dop_upper(j).getInUnits("Hz")<<endl;  
    }

    for (unsigned int i = 0; i<Nrange_patch;i++)
    {
    for(unsigned int j = 0; j <Ndop_patch;j++)
    {
      outputfile<<avg_incidence(i,j).getInUnits("deg")<<endl;
      outputfile<<avg_backscatter(i,j)<<endl;
      outputfile<<avg_antenna(i,j)<<endl;
      outputfile<<avg_azi(i,j).getInUnits("deg")<<endl; 
      outputfile<<avg_elev(i,j).getInUnits("deg")<<endl;
      outputfile<<patch_sum(i,j).getInUnits("kg m m/(s s s)")<<endl;
      outputfile<<avg_area(i,j).getInUnits("km km")<<endl;
    }
    }

    for (unsigned int i =0 ; i < Nrange_patch;i++)
    {
    for (unsigned int j = 0; j < Ndop_patch;j++)
    {
    for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
    {
    for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
    {
      unsigned int patch_index = i*Ndop_patch+j;
      unsigned int bin_index = i_range*Ndop_bin+j_dop;	
      outputfile<<lat(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<lon(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<thetai(patch_index,bin_index).getInUnits("deg")<<endl;
      outputfile<<range(patch_index,bin_index).getInUnits("km")<<endl;
      outputfile<<doppler(patch_index,bin_index).getInUnits("Hz")<<endl;
      outputfile<<alongtrack(patch_index,bin_index).getInUnits("km")<<endl;
      outputfile<<crosstrack(patch_index,bin_index).getInUnits("km")<<endl;
      
    }//Loop over dop bin
    }//Loop over range bin
    }//Loop over dop patch
    }//Loop over range patch
        
    for (unsigned int i_range = 0; i_range<Nrange_bin; i_range++)
    {
    for (unsigned int j_dop =0; j_dop<Ndop_bin; j_dop++)
    {
      outputfile<<bin_power_main(i_range,j_dop).getInUnits("kg m m/(s s s)")<<endl;   
      outputfile<<bin_power_amb(i_range,j_dop).getInUnits("kg m m/(s s s)")<<endl;
      outputfile<<normalized_bin_power(i_range,j_dop).getInUnits("kg m m/(s s s)")<<endl;
      outputfile<<center_area(i_range,j_dop).getInUnits("km km")<<endl;
      outputfile<<center_patch_antenna(i_range,j_dop)<<endl;
      outputfile<<is_good_bin_amb(i_range,j_dop) <<endl;
      outputfile<<is_good_bin_nesigma0(i_range,j_dop)<<endl;
      outputfile<<is_good_bin_totalSNR(i_range,j_dop) <<endl;
      outputfile<<is_good_bin_pn_to_amb(i_range,j_dop)<<endl;
      outputfile<<is_good_bin(i_range,j_dop)<<endl;
    }
    }
  }//END of cal_amb


//------------------------------------------------
//Calculate amb table for BS
//-----------------------------------------------
void cal_ambiguity_geometry_table(const string& target_name,
				  Beam& beam,
				  const Frame& ftitan,
				  const StateVector& sc_state,
				  const Uvar& lambda,
				  const double&  muhleman_k1,
				  const double&  muhleman_k2,
				  const Uvec& dop_freq_values,
				  const Uvec& range_values,
				  Dmat& cell_one_way_gain,
				  Dmat& cell_radar_geom_factor,
				  Umat& cell_cross_track,
				  Umat& cell_along_track,
				  Umat& cell_incidenceangle,
				  Umat& cell_area,
				  Umat& cell_lat,
				  Umat& cell_lon,
				  Dmat& cell_sigma0,
				  Uvec& range_axis,
				  Uvec& dop_axis)
  throw(Unit::UnitError,ErrorMessage)
  {
  typedef Array1D<DirectionVector> Dirvec;    
  
  Time t = sc_state.time();
  Uvar speed_light = Uvar("speed_light");       


  //get max beam gain
  //double max_gain = beam.getMaxGain();
  Uvar azimuth_width = beam.getAzimuthWidthOneWay();
  Uvar elevation_width=beam.getElevationWidthOneWay();
  
  //------------------------------------------
  //First use Frame method to get range and dop
  //------------------------------------------
  unsigned int beamnum = beam.getBeamNumber();
  Frame fbeam("CASSINI_RADAR_" + toStr(beamnum),"Cassini");
  DirectionVector boresight("boresight",fbeam,t,0,0,1);
  Uvar range_F, dop_F, thetai_F,altitude_F;
  TargetGeom tg(t);
  tg.setState(sc_state);
  tg.setLookDirection(boresight);
  tg.setTarget(target_name,ftitan);
  Uvar radius = tg.radius();
  range_F = tg.range();
  dop_F = tg.doppler(lambda); 
  thetai_F = tg.incidenceAngle();
  altitude_F = tg.altitude();
  
  //-----------------------------------------
  //Then, use non-frame method(vector equations)
  // to calculate boresight range and dop
  //------------------------------------------  
  PositionVector sc_pos = sc_state.position();
  sc_pos.representIn(ftitan);
  DirectionVector pos_dir = sc_pos;
  pos_dir.representIn(fbeam);   
  
  FloatVector velocity = sc_state.velocity();
  velocity.representIn(fbeam);
  DirectionVector velocity_direction = velocity;
  
  //---------
  //doppler
  //---------
  double v_b = dot(velocity_direction,boresight);
  Uvar dop0 = 2.0* v_b * velocity.magnitude()/lambda;
  
  //---------
  //range
  //---------
  double p_b = dot(pos_dir,boresight);
  Uvar in_sqrt = sc_pos.magnitude() *sc_pos.magnitude() * p_b * p_b;
  in_sqrt += radius * radius - sc_pos.magnitude()*sc_pos.magnitude();
  Uvar range0 = Uvar(0.0);
  if(in_sqrt.getValue() > 0.0) 
    {
      range0 = -sc_pos.magnitude() * p_b - sqrt(in_sqrt);
    }
  else
    {
      throw ErrorMessage("No intercepting point in the boresight direction");
      
    }  
  
  //-----------------------------------------
  //compare frame-method with non frame method
  //------------------------------------------
  double diff_range = range0.getInUnits("km") - range_F.getInUnits("km");
  double diff_dop   = dop0.getInUnits("Hz") - dop_F.getInUnits("Hz");  
  if (fabs(diff_range)>0.001)
    {
      throw ErrorMessage("Frame method does not produce the same range");
    }
  if (fabs(diff_dop) > 0.001)
    {
      throw ErrorMessage("Frame method does not produce the same dop");
    }
  
  
  
 
  unsigned int Ngrid = dop_freq_values.size();
  if (Ngrid != range_values.size()) throw ErrorMessage("Need square matrix");

 
	

    Ivec no_of_solution("no_solution",Ngrid*Ngrid);
    Uvec look_azimuth("look_azimuth_angle",Ngrid*Ngrid);
    Uvec look_elevation("look_elevation_angle",Ngrid*Ngrid);
    Uvec  lat("bin_patch_lat",Ngrid*Ngrid);
    Uvec  lon("bin_patch_lon",Ngrid*Ngrid);
    Uvec  range("range",Ngrid*Ngrid);
    Uvec  doppler("doppler",Ngrid*Ngrid);
    Uvec  thetai("thetai",Ngrid*Ngrid);
    Uvec  alongtrack("alongtrack",Ngrid*Ngrid);
    Uvec crosstrack("crosstrack",Ngrid*Ngrid);
    Uvec area("area_element",Ngrid*Ngrid);
    Dirvec surface_intercept("surface_intercept",Ngrid*Ngrid);
    Dvec beam_gain("beam_gain",Ngrid*Ngrid);
    Dvec sigma0_elem("sigma0_elem",Ngrid*Ngrid);
    
    //Initialize variables
    no_of_solution= 0;
    look_azimuth=Uvar(0.0);
    look_elevation= Uvar(0.0);
    lat= Uvar(0.0);
    lon= Uvar(0.0);
    range= Uvar(0.0);
    doppler= Uvar(0.0);
    thetai= Uvar(0.0);
    alongtrack=Uvar(0.0);
    crosstrack=Uvar(0.0);
    area= Uvar(0.0);
    beam_gain= 0.0;
    sigma0_elem= 0.0;

   

    //------------------------------------------------------------------
    //Find  intercepting points corresponding range and dop values
    //no_interception: store how many range and dop bins have a solution
    //------------------------------------------------------------------  
    TargetGeom tg1(t); 
    for (unsigned int i_range = 0; i_range<Ngrid; ++i_range)
    for (unsigned int j_dop =0; j_dop<Ngrid; ++j_dop)
      {
      unsigned int bin_index = i_range * Ngrid + j_dop;	
      tg1.setState(sc_state);
      tg1.setTarget(target_name,ftitan);
      tg1.setRangeDoppler(range0+range_values(i_range), dop0+dop_freq_values(j_dop),lambda);
      DirectionVector look = tg1.lookDirection();
  
      if (tg1.foundSurfaceIntercept()==false)
	{
	no_of_solution(bin_index) = 0;//no surface intercept point
	lat(bin_index)= Uvar(0,"rad");
	lon(bin_index) =  Uvar(0,"rad");
	range(bin_index)= range0+range_values(i_range);
	doppler(bin_index) = dop0+dop_freq_values(j_dop);
	thetai(bin_index)= Uvar(0,"rad"); 
	alongtrack(bin_index)=Uvar(0,"km");
	crosstrack(bin_index)=Uvar(0,"km");
	surface_intercept(bin_index)=DirectionVector();
	beam_gain(bin_index) = 0.0;
	sigma0_elem(bin_index) = 0.0;
	}
      else
	{
	no_of_solution(bin_index) = 1;   
	DirectionVector new_look= look.representIn(fbeam);
	new_look.getAzimuthElevation(look_azimuth(bin_index),
				     look_elevation(bin_index));
	lat(bin_index)= tg1.lat();
	lon(bin_index) = tg1.lon();
	range(bin_index)= tg1.range();
	doppler(bin_index) =tg1.doppler(lambda);
	thetai(bin_index)= tg1.incidenceAngle(); 
	tg1.interceptAlongCross(alongtrack(bin_index),
				crosstrack(bin_index));
	surface_intercept(bin_index)=tg1.surfaceIntercept();
	beam_gain(bin_index) =  beam.bilinear(look_azimuth(bin_index),
					      look_elevation(bin_index));
	sigma0_elem(bin_index)= muhleman_k1 
	  * cos(thetai(bin_index))
	  / pow(sin(thetai(bin_index))
		+muhleman_k2*cos(thetai(bin_index))
		,3.0);
	} 
      tg1.reset(t);
      }//Loop over i_range and j_dop
 

    //    cout<<"Finish finding all surface intercept points"<<endl;
    //    cout<<"Start calculating area elements and radar return signal "<<endl;
    //------------------------------------------	
    //End of finding all range and dop solutions
    //if there is a surface intercept point, no_of_s
    //----------------------------------------------
    //****************************************
    //    Main Ambiguity Calculation
    //there is a solution for the given  range and dop value
    //range_bin(i,i_range) and dop_bin(j,j_dop)
    //let's try to get radar return signal strenth here
    //*******************************
    for (unsigned int i_range = 0; i_range<Ngrid-1; ++i_range)
    for (unsigned int j_dop =0; j_dop<Ngrid-1; ++j_dop)
      {
      unsigned int bin_index = i_range * Ngrid + j_dop;	
      unsigned int bin_Up,  bin_Left,bin_Diag;
      bin_Up = (i_range+1)*Ngrid+j_dop;	
      bin_Left = i_range*Ngrid+(j_dop+1);
      bin_Diag = (i_range+1)*Ngrid+j_dop+1;
      if (no_of_solution(bin_index)
	  *no_of_solution(bin_Up)
	  *no_of_solution(bin_Left)
	  *no_of_solution(bin_Diag) == 1)
	{
	//------------------------------------
	//find Area element covered by (range,dop)
	//divide an area enclosed by two range and dop values
        // bin_Up (range++), bin_Left(dop++), bin_Diag(range++,dop++)
        //into two triangles and calclate their areas and add them
        //together
        //-------------------------------------           
	Uvar side_a = radius * 
	  fabs((surface_intercept(bin_index)
		.angle(surface_intercept(bin_Left))
		.getInUnits("rad")));
	Uvar side_b = radius * 
	  fabs((surface_intercept(bin_index)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	Uvar side_c = radius * 
	  fabs((surface_intercept(bin_Left)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	Uvar cos_c = (side_a*side_a + side_b*side_b - side_c*side_c)
	  /(2.0 * side_a*side_b);
	Uvar angle_c = Uvar(acos(cos_c.getValue()),"rad");
	
	area(bin_index) = 0.5 * sin(angle_c) * side_a * side_b;
	
	Uvar side_d = radius * 
	  fabs((surface_intercept(bin_Left)
		.angle(surface_intercept(bin_Diag))
		.getInUnits("rad")));
	Uvar side_e = radius * 
	  fabs((surface_intercept(bin_Diag)
		.angle(surface_intercept(bin_Up))
		.getInUnits("rad")));
	cos_c = (side_d*side_d + side_e*side_e - side_c*side_c)
	  /(2.0 * side_d*side_e);
	angle_c = Uvar(acos(cos_c.getValue()),"rad");
	area(bin_index) += 0.5 * sin(angle_c) * side_d * side_e;	    
	
	
	//-----------------------------
	//Store results into Ngrid-1 X Ngrid-1 matrix form
	//-----------------------------    
	
	//----------------------------
	//cell area
	//---------------------------
	cell_area(i_range,j_dop) = area(bin_index);
	
	//--------------------
	//cell range and doppler
	//-----------------
	range_axis(i_range) = 
	  (range(bin_index)+range(bin_Up)
	   +range(bin_Left)+range(bin_Diag))/4.0;
	

	dop_axis(j_dop) = 
	  (doppler(bin_index)+doppler(bin_Up)
	   +doppler(bin_Left)+doppler(bin_Diag))/4.0;

	//-----------------------------
	//store beam gain and sigm0
	//----------------------------
	cell_one_way_gain(i_range,j_dop) = 
	  (beam_gain(bin_index)+beam_gain(bin_Up)
	   +beam_gain(bin_Left)+beam_gain(bin_Diag))/4.0;
	
	
	cell_sigma0(i_range,j_dop) = 
	  (sigma0_elem(bin_index) + sigma0_elem(bin_Up)
	   +sigma0_elem(bin_Left) + sigma0_elem(bin_Diag))/4.0;
	
	//-----------------------------------
	//incidenceangle
	//----------------------------------
	cell_incidenceangle(i_range,j_dop) = 
	  (thetai(bin_index) + thetai(bin_Up)
	   +thetai(bin_Left)+thetai(bin_Diag))/4.0;
	
	
	//------------------------------------------------------
	//crosstrack and alongtrack distances
	//------------------------------------------------------ 
	cell_cross_track(i_range,j_dop) =
	  ( crosstrack(bin_index)+crosstrack(bin_Up)
	    +crosstrack(bin_Left)+crosstrack(bin_Diag))/4.0;
	
	cell_along_track(i_range,j_dop) = 
	  (alongtrack(bin_index)+alongtrack(bin_Up)
	   +alongtrack(bin_Left)+alongtrack(bin_Diag))/4.0;
	
	//------------------------
	//lat and lon values
	//--------------------------
	cell_lat(i_range,j_dop)=
	  (lat(bin_index) + lat(bin_Up)
	   + lat(bin_Left)+ lat(bin_Diag))/4.0;
	
	cell_lon(i_range,j_dop)=
	  (lon(bin_index) + lon(bin_Up)
	   + lon(bin_Left)+ lon(bin_Diag))/4.0;
	
	//---------------------------------------
	//find X = G^2 * lambda^2 * sigma0 / R^4
	//--------------------------------------
	
	Uvar X = 
	  cell_one_way_gain(i_range,j_dop) * cell_one_way_gain(i_range,j_dop)
	  * lambda * lambda
	  *cell_sigma0(i_range,j_dop)
	  /pow(range_axis(i_range),4);  	     
	
	cell_radar_geom_factor(i_range,j_dop)= 
	  (X*cell_area(i_range,j_dop)).getValue();  
	}
      else
	{
	cell_area(i_range,j_dop) = Uvar(0,"km km");
	range_axis(i_range) = 
	  (range(bin_index)+range(bin_Up)
	   +range(bin_Left)+range(bin_Diag))/4.0;
	
	dop_axis(j_dop) = 
	  (doppler(bin_index)+doppler(bin_Up)
	   +doppler(bin_Left)+doppler(bin_Diag))/4.0;	  
	cell_one_way_gain(i_range,j_dop) = 0.0;
	cell_sigma0(i_range,j_dop) = 0.0;
	cell_incidenceangle(i_range,j_dop) = Uvar(0,"deg");
	cell_cross_track(i_range,j_dop) = Uvar(0,"km");
	cell_along_track(i_range,j_dop)=Uvar(0,"km");
	cell_lat(i_range,j_dop) = Uvar(0,"deg");
	cell_lon(i_range,j_dop) = Uvar(0,"deg");
	cell_radar_geom_factor(i_range,j_dop) = 0.0;
	}//no solution, return zero
      }//Loop over range/dop values   

    //------------------------
    //display range and dop axis on the screen
    //-----------------------
    //for (unsigned int i = 0; i < Ngrid-1;++i)
    //{
    //cout<<"range and dop "<< range_axis(i)<<" "<<dop_axis(i)<<endl;
    //}

  }




//-------------------------------------------------
// calculate sinc pattern beam gain 
//-------------------------------------------------
double sincgain(const Uvar& azim, 
		const Uvar& elev, 
		const Uvar& azimuth_width,
		const Uvar& elevation_width,
		const double& max_gain) throw(Unit::UnitError,ErrorMessage)
  {
    Uvar au,bv;
    double sau,sbv;
    
    au = 0.88 * azim;
    au /= azimuth_width;
    bv = 0.88* elev;
    bv /= elevation_width;
    sau = std::sinc(au.getValue()); 
    sbv = std::sinc(bv.getValue());
    return( sau * sau * sbv * sbv *max_gain); 
  }








