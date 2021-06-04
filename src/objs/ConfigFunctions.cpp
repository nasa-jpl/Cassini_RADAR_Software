//----------------------------------------------------------------------------
// ConfigFunctions.cpp
//
// This file contains definitions for functions that receive a config object
// and perform some kind of initialization action using some of the config
// parameters.  Results may be returned as arguments or return values.
//----------------------------------------------------------------------------

//----------------------
// Configuration Control
//----------------------

static const char rcs_id_configfunctions_c[] =
  "@(#) $Id: ConfigFunctions.cpp,v 11.5 2011/09/16 00:03:30 richw Exp $";

#include <string>
#include "Frame.h"
#include "Units.h"
#include "Time.h"
#include "TargetGeom.h"
#include "ConfigFunctions.h"
#include "TemplateUtils.h"


using std::cout;
using std::cout;
using std::cout;

/**** Subsumed in class Flyby

//--------------------------------------------------------------------------
// config_time_specification(cfg,epoch_time,start_time,end_time,interval)
//
// Read the config parameters in cfg which specify a time interval to deal
// with.  The results are returned as an epoch_time, and an epoch
// relative set of time steps from start_time to end_time stepping by
// interval.
//
// The following config keywords are processed by this function:
//
// epoch_selection (absolute,closest_flyby,designated_flyby,ckernel_valid)
// epoch_time (UTC time string)
// time_step
// start_time
// end_time
//
//--------------------------------------------------------------------------

void config_time_specification(const Config& cfg, Time& epoch_time,
  Uvar& start_time, Uvar& end_time, Uvar& interval)
  throw(ErrorMessage,Unit::UnitError)
  {
  string epoch_selection = cfg.str("epoch_selection");
  if (epoch_selection == "absolute")
    {
    //-------------------------------------------------------------------
    // User supplies epoch, start, end, and interval times directly.
    //-------------------------------------------------------------------

    string epoch_str = cfg.str("epoch_time");
    epoch_time.setUtc(epoch_str);
    interval = cfg["time_step"];
    start_time = cfg["start_time"];
    end_time = cfg["end_time"];
    }
  else if (epoch_selection == "closest_flyby")
    {
    //-------------------------------------------------------------------
    // User supplies an epoch_time close to the desired flyby, and epoch
    // relative start and end times and the interval.  The epoch of the
    // closest approach of the flyby closest to the input epoch_time is
    // located and used as the epoch_time.  In this case, epoch_time is
    // an input/output variable and is modified.
    //-------------------------------------------------------------------

    string epoch_str = cfg.str("epoch_time");
    epoch_time.setUtc(epoch_str);
    interval = cfg["time_step"];
    start_time = cfg["start_time"];
    end_time = cfg["end_time"];
    Uvar epoch_accuracy("epoch_accuracy");
    epoch_accuracy = cfg["epoch_accuracy"];

    //-------------------------------------------------------------
    // Locate nearest closest approach using golden section search.
    //-------------------------------------------------------------

    Uvar t_mid("t_mid");
    t_mid = epoch_time + Uvar(0.001,"s");
    Uvar t2("t2");
    Uvar final_time("final_time");
    Uvar lowest_altitude("lowest_altitude");
    TargetAltitude tga("Cassini","Titan");
    bracket_minimum(epoch_time,t_mid,t2,tga);
    golden_section_search(epoch_time,t_mid,t2,epoch_accuracy,tga,
      epoch_time,lowest_altitude);

    }
  else if (epoch_selection == "designated_flyby")
    {
    throw
      ErrorMessage("epoch_selection = designated_flyby not implemented yet");
    }
  else if (epoch_selection == "ckernel_valid")
    {
    throw ErrorMessage("epoch_selection = ckernel_valid not implemented yet");
    }
  else
    {
    throw ErrorMessage("Invalid epoch_selection parameter: " + epoch_selection);
    }

  }

**/
