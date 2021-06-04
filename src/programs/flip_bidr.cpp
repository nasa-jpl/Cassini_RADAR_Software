//----------------------------------------------------------------------
// subsample_beammask:  Generate a reduced-resolution version of a BIDR
//                      beam mask product.
//
// Usage:  subsample_beammask [--ppd=N] inputfile outputfile
//
// where --ppd is the new resolution (pixels per degree) with legal
// values 2, 8, 32, 128, or 256; --ppd must also be less than or equal
// to the resolution of the input product.  If --ppd is omitted, the
// new resolution defaults to 128.  "inputfile" must be an existing BIDR
// beam mask product.  The reduced-resolution product is saved in
// "outputfile."  The reduced-resolution product is formed by
// inclusive-OR'ing the original pixels.  The program updates the
// relevant values in the PDS label to reflect the new resolution,
// including PRODUCT_ID.
//
// Note:  this version of the program acts upon the entire extent of the
// input product.  Specification of a subregion to process is only
// partially supported within the BIDRFile and BIDRFileUpdates classes.
//----------------------------------------------------------------------
#include <typeinfo>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include <getopt.h>
#include "BIDRFile.h"
#include "Error.h"
#include "Utils.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;

static const char rcs_id_subsample_beammask[] =
  "@(#) $Id: flip_bidr.cpp,v 11.6 2017/04/10 23:49:10 cveerama Exp $";

void usage(char *cmd)
{
  cerr << cmd << ": input_file output_file" << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  try
  {
    char *cmd = argv[0];
    char *bidr_input = NULL;
    char *bidr_output = NULL;


    // Parse command line options

    if (argc != 3)
    {
      usage(cmd);
    }
    int clidx=1;
    bidr_input = argv[clidx++];
    bidr_output = argv[clidx++];

    if (!fileExists(bidr_input))
    {
      ErrorMessage e("input file \"" + string(bidr_input) +
        "\" either not found or is unreadable.");
      e.throwMe();
    }

    BIDRFile bf_in(bidr_input,"r");
    bf_in.readHeader();
    if (bf_in.type() == INVALID) {
      ErrorMessage e(argv[0] + string(" cannot handle input") +
        string(" BIDR products of type INVALID."));
      e.throwMe();
    }
    BIDRFile bf_out(bidr_output,"w");

    BIDRFileUpdates updates(bidr_input, bidr_output, 0, 0);
    updates.flip();
    updates.writeBufferToOutput();
  }
  catch(ErrorMessage& e)
  {
    cerr << "Error: " << e.msg << endl;
  }
  catch(std::exception& ex)
  {
    cerr << "Error: " << typeid(ex).name() << ": " << ex.what() << endl;
  }
  catch(...)
  {
    cerr << "Exception caught" << endl;
  }
  return 0;
}


