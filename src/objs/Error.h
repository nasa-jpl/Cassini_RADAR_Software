#ifndef Error_H
#define Error_H

#include <string>

//-------------------------------------------------------------------------
// General error class that simply stores a message that can be set through
// a constructor.  This class provides an easy error handling system for
// situations where you just want to quit with a message (seems to be
// the case most of the time).
//-------------------------------------------------------------------------

class ErrorMessage
  {
  public:
  std::string msg;
  ErrorMessage() { msg = "Error exception"; }
  ErrorMessage(const char* emsg) { msg = emsg; }

  ErrorMessage(const std::string& emsg) { msg = emsg; }
  void throwMe() { throw *this; }
  };

#endif
