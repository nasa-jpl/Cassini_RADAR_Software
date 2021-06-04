// ----------------------------------------------------------
//                      RMSS FUNCTIONS  
// ----------------------------------------------------------

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>


#include "rmss_functions.h"
//#include "classes.h"
//#include "Ieb.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::terminate;

// ----------------------------------------------------------
//                      VERSION
// ----------------------------------------------------------
        
void RMSS_FUNCTIONS_VERSION ()
{
cout << "RMSS_FUNCTIONS_VERSION: 03/19/2004...not used any more..." << endl;
}

 
// ----------------------------------------------------------
//                      GETFILENAME PROCESS
// ----------------------------------------------------------
int to_report (char *the_string)
{
static FILE *reportfile;
char report_extension [] = ".rmss_rpt" ;
//char pwd_directory [DIR_STRING_SIZE] ;
static char wreportname[PRINTSTRINGSIZE] ;
static int firsttime = TRUE;
static int printscreen = DISABLE ;
char filename [] = "default_4now";


  if (firsttime == TRUE ){
        firsttime = FALSE;
//      if ( (getenv ("RMSS_PRINTTOOL")) == "ENABLE" ) 
//              printscreen = ENABLE ; 
 
//        sprintf ( wreportname , "%s%s", getfilename (), report_extension );
        sprintf ( wreportname , "%s%s", the_string, report_extension );

        if ((reportfile = fopen (wreportname,  "w")) == NULL )
        {       printf ( "cannot open file. \n\n" );

                exit (0) ;
        }
        printf ( "\n!! opened 'to_report' file: '%s'\n", wreportname );
        fprintf (reportfile, "\nREPORT_VERSION: 03/19/2004"); 
  }


/* WRITE TO FILE or CLOSE*/


if (  !(strcmp (the_string, "CLOSE"))  ) {

        fclose (reportfile);
        printf ("\nClosing %s\n\n", wreportname);

} else {
         fprintf (reportfile, "\n%s", the_string);
         if (printscreen == ENABLE ) cout << the_string << endl ;
}


return 0;
}


//----------------------------------------------------








