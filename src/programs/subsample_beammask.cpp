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

#include <iostream>
#include <stdlib.h>
#include <typeinfo>
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
  "@(#) $Id: subsample_beammask.cpp,v 11.5 2011/09/16 00:03:31 richw Exp $";

void usage(char *cmd);

int main(int argc, char* argv[])
{
  try
  {
    int c;
    int ppd_new = DEFAULT_NEW_BEAMMASK_RESOLUTION;
    int option_index = 0;
    int flagged_ppd_new = 1; // Always check ppd_new even when the default is used
    // int line_start;
    // int line_end;
    // int pixel_start;
    // int pixel_end;
    extern int optind;
    extern char *optarg;
    char *cmd = argv[0];
    char *bidr_input = NULL;
    char *bidr_output = NULL;
    struct option long_options[] =
    {
      {"ppd", 1, 0, PPD_OPTION_MARKER},
      {0, 0, 0, 0}
    };

    // Parse command line options
    while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1)
    {
      switch (c)
      {
	case PPD_OPTION_MARKER:
	  ppd_new = atoi(optarg);
	  flagged_ppd_new++;
	break;
	default:
	  ;
	break;
      }
    }

    if (argc - optind != 2)
    {
      usage(cmd);
    }
    bidr_input = argv[optind++];
    bidr_output = argv[optind++];

    if (!fileExists(bidr_input))
    {
      ErrorMessage e("input file \"" + string(bidr_input) +
        "\" either not found or is unreadable.");
      e.throwMe();
    }

    BIDRFile bf_in(bidr_input,"r");
    bf_in.readHeader();
    if (bf_in.type() != BEAMMASK) {
      ErrorMessage e(argv[0] + string(" only handles input") +
        string(" BIDR products of type BEAMMASK."));
      e.throwMe();
    }
    BIDRFile bf_out(bidr_output,"w");

    BIDRFileUpdates updates(bidr_input, bidr_output, flagged_ppd_new, ppd_new);
    updates.average();
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

void usage(char *cmd)
{
  cerr << "Usage:  " << cmd <<
    " [--ppd=<new resolution in pixels/degree>]" <<
    " <input file>" <<
    " <output file>" <<
    endl;
  cerr << "where --ppd is 2, 8, 32, 128, or 256 AND";
  cerr << " less than or equal to the input product's resolution." <<
    endl;
  cerr << "If omitted, --ppd defaults to 128." <<
    endl;
  exit(1);
}
