//============================================================================
// ConfigFunctions.h
//
// This file contains declarations for functions that receive a config object
// and perform some kind of initialization action using some of the config
// parameters.  Results may be returned as arguments or return values.
//============================================================================

#ifndef ConfigFunctions_H
#define ConfigFunctions_H

/** Subsumed by class Flyby
void config_time_specification(const Config& cfg, Time& epoch_time,
  Uvar& start_time, Uvar& end_time, Uvar& interval)
  throw(ErrorMessage,Unit::UnitError);
**/

/**** Moved to Flyby.h

class TargetAltitude
  {
  public:

  TargetAltitude(const string& observer_name, const string& target_name)
    : observer_name_(observer_name), target_name_(target_name)
    { }

  Uvar operator()(const Uvar& t)
    {
    TargetGeom tg(t);
    tg.setState(observer_name_);
    tg.setTarget(target_name_);
    return(tg.altitude());
    }

  private:

  string observer_name_;
  string target_name_;
  };
***/

#endif
