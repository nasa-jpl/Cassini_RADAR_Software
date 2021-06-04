#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H
#include"Config.h"

class DebugInfo{
 public:
  DebugInfo(const string& routine_name);
  DebugInfo(int l){level=l;}
  static void config(Config& cfg);
  static void report(ofstream& fs);
  static void report();
  static void setLevel(const string& routine_name, int l);
  static ofstream file;
  static int global_level;
  static bool allWarnings;
  int level;

 private:

  typedef map<string,int> integers;
  static integers level_map_;
  

};
#endif
