PDS_VERSION_ID       = PDS3

RECORD_TYPE          = STREAM
INSTRUMENT_HOST_NAME = "CASSINI ORBITER"
INSTRUMENT_NAME      = "CASSINI RADAR"
OBJECT               = TEXT
  PUBLICATION_DATE   = XXXX
  NOTE               = "Description of Volume VOLUME ID GOES HERE"
END_OBJECT           = TEXT
END


Volume VOLUME ID GOES HERE:  VOLUME NAME GOES HERE

1.  Introduction

This volume contains Cassini Radar Burst Order Data Products (BODP)
and, if applicable, Basic Image Data Record (BIDR) products and 
Synthetic Aperture Radar Topograpgy (SARTopo) products for a single 
observation (data take).  Each product contains the results from
processing of the Cassini Radar instrument in passive (radiometer) 
and/or active (scatterometer, altimeter, synthetic aperture radar) 
modes. A given observation may contain data for only a subset of these
modes. The volume also contains product documentation and the set of
antenna pattern beam files and processing parameters that were used 
for these products.

The products on this volume are the exact products released by the
Cassini Project.  Supporting documentation and labels conform to the
Planetary Data System (PDS) Standards Reference, Version 3.6, Jet
Propulsion Laboratory (JPL) Document D-7669, Part 2, August 1, 2003; and
the Planetary Science Data Dictionary Document, Jet Propulsion Laboratory
(JPL) Document D-7116, Rev. E, August 28, 2002, Online Update Version
pdscat1r57, April 25, 2005.  The most recent copies of both documents are
available at the PDS web site, http://pds.jpl.nasa.gov .

For more information about the volume directory structure, refer to the
Volume Software Interface Specification for Cassini Radar Instrument Team
Data Products, JPL D-27890.  A copy of the document is located on this volume
as file VOLSIS.PDF in the DOCUMENT directory.


2.  Disk Format

Not applicable, as these products are being delivered to PDS electronically.
Each volume, however, has been constrained to a maximum size of 4.3 GB so
that it can fit on a single DVD-R disc.


3.  File Formats

Science Data File Format

Four types of BODP science data files have been defined.

  Short Burst Data Record (SBDR) Product
    Instrument telemetry and calibrated science data in burst order

  Long Burst Data Record (LBDR) Product
    Instrument telemetry, calibrated science data, and time-sampled
    echo data in burst order; only includes records for which radar
    was in active mode (scatterometer, altimeter, synthetic aperture
    radar)

  Altimetry Burst Data Record (ABDR) Product
    Instrument telemetry, calibrated science data, and altimeter profile
    data in burst order; only includes records for which radar was in
    altimeter mode

  Altimetry Burst Data Record Summary (ABDR-SUMMARY) Product
    Estimated surface height and related measurements derived from the
    average of all pulses in a burst, in burst order; only includes records
    for which radar was in altimeter mode

The SBDR product should be available for all observations.  LBDR products
will not be available if no active-mode data was collected, and ABDR and
ABDR-SUMMARY products will only be available if altimetry data was collected.

The Cassini BODP science data files are formatted using the "table
object" structure of the PDS standards.  Each row of the table consists
of a data record comprising radar and ancillary data for a single burst.
Rows occur in time order.  A column represents a particular data field
with a specific meaning and data type (i.e. integer, real, etc.)  Real
values are either 32-bit or 64-bit IEEE floating-point values with little
endian byte order (PC_REAL).

Cassini BIDR products result from the processing of synthetic aperture
radar (SAR) data into imagery (sigma_0 backscatter), incidence angle,
and related products.  The Cassini BIDR science data files are formatted
using the "image object" structure of the PDS standards.  The size and
coordinate system of the image is defined in the attached PDS label.
Image pixels may have several different meanings and data types.

Cassini SARTopo products contain surface height estimates (and related
measurements) of Titan in profiles between two SAR beams within a single 
SAR image (BIDR). The size of the product is defined in the attached PDS 
label. These products are formatted using the "table object" structure of
the PDS standards similar to the ABDR-SUMMARY products.

To conserve space, each ABDR, LBDR and BIDR product on this volume has been
compressed with the ZIP utility and stored in its own archive with the
extension ".ZIP".  Use the UNZIP utility to extract the original product
file from a ZIP archive.  The source code for ZIP and UNZIP is publicly
available, but both utilities come pre-installed on many types of Unix and
Linux systems as the commands "zip" and "unzip," respectively.  On those
systems, the command-line syntax for extracting a file with "unzip" is

    unzip file.ZIP -d output_directory

where "file.ZIP" is the ZIP archive and "output_directory" is the name of the
directory into which the uncompressed product will be written.  The output
directory must be writable by the user and be located on a file system with
sufficient free space to hold the uncompressed product.  The command

    unzip -l file.ZIP

will list the name and uncompressed file size (in bytes) of the file stored
in archive "file.ZIP."  The uncompressed file size is also stored as the
value of REQUIRED_STORAGE_BYTES in the ZIP archive's detached label, which
is located in the same directory as the ZIP archive on the volume.

For more information about the format and content of the BODP files, refer
to the Cassini Radar Burst Ordered Data Product Software Interface
Specification, JPL Document D-27891.  A copy of the document is located on
this volume as file BODPSIS.PDF in the DOCUMENT directory.

For more information about the format and content of the BIDR files and 
SARTopo files, refer to the Cassini Radar Basic Image Data Record Software
Interface Specification, JPL D-27889.  A copy of the document is located 
on this volume as file BIDRSIS.PDF in the DOCUMENT directory.

For a pointer to more information about the ZIP/UNZIP utilities, refer to
the software catalog SOFTWARE.CAT in the CATALOG directory on this volume.

Document File Format

Document files with the .TXT suffix exist in the root, CALIB, CATALOG,
DOCUMENT, EXTRAS, INDEX, and SOFTWARE directories and their subdirectories.
These files are ASCII files with attached PDS labels and with line lengths
that do not exceed 80 characters.  Document files have variable-length lines
with an end-of-line designator consisting of a carriage-return character
followed by a linefeed character.  This allows the files to be read by the
MacOS, Microsoft Windows, and Unix operating systems.

Tabular File Format

A tabular file (INDEX.TAB) exists in the INDEX directory.  This file is
an ASCII file formatted for direct reading into many database management
systems on various computers.  Commas separate all fields, and character
fields are enclosed in double quotation marks (").  Character fields
are left justified and padded on the right with spaces to keep quotation
marks in the same columns of successive records.  Number fields are right
justified and padded with spaces on the left to keep fixed-width fields.
The records are of fixed length, and the last two bytes of each record
contain the ASCII carriage-return and line-feed character.  This allows
a table to be treated as a fixed length record file on computers that
support this file type and as a text file with embedded line delimiters
on those that don't.

A detached PDS label describes INDEX.TAB.  The PDS label file name is the
same as the tabular file but with the extension ".LBL".

PDS Label Format

All data files in the Cassini Radar Instrument Team Archive Collection
have PDS labels.  A PDS label, whether attached to or detached from its
associated file, provides descriptive information about the associated
file.  The PDS label is an object-oriented structure consisting of sets of
"keyword = value" declarations.  The object to which the label refers
(e.g.  IMAGE, TABLE, etc.) denotes the form:

    ^object = location

in which the caret character (^, also called a pointer in this context)
indicates where to find the object.  In an attached label, the location
is an integer representing the starting record number of the object
(the first record in the file is record 1).  In a detached label, the
location denotes the name of the file containing the object, along with
the starting record or byte number, if there is more than one object
in the file.

For example:

    ^HEADER = ("F01.IMG",1)
    ^IMAGE = ("F01.IMG",1025 <BYTES>)

indicates that the IMAGE object begins at byte 1025 of the file F01.IMG,
in the same directory as the detached label file.  Below is a list of
the possible formats for the ^object definition.

    ^object = n
    ^object = n<BYTES>
    ^object = "filename.ext"
    ^object = ("filename.ext",n)
    ^object = ("[dirlist]filename.ext",n)
    ^object = ("filename.ext",n<BYTES>)
    ^object = ("[dirlist]filename.ext",n<BYTES>)

where:

"n" is the starting record or byte number of the object, counting from
the beginning of the file (record 1, byte 1),

"<BYTES>" indicates that the number given is in units of bytes,

"filename" is the up to 27-character, alphanumeric upper-case file name,

"ext" is the 3 character upper-case file extension,

"[dirlist]" is a period-delimited path-list of parent directories, in
upper case, that specifies the object file directory (used only when the
object is not in the same directory as the label file).  The list begins
at the directory level below the root directory of the CDROM.  "[dirlist]"
may be omitted when the object being described is located in the same
directory as the detached label.

Detached label files have variable length lines with an end-of-line
designator consisting of a carriage-return character followed by a
line-feed character.  This allows the files to be read by the MacOS,
Microsoft Windows, and Unix operating Systems.

Catalog File Format

Catalog files (suffix.CAT) exist in the root and CATALOG directories.
They are formatted in an object-oriented structure consisting of sets of
"keyword = value" declarations.


4.  Volume Contents

The files in this volume are organized in one top-level directory with
several subdirectories.  The following outline shows the structure and
content of these directories.  In the table, directory names are enclosed
in square brackets ([]), upper-case letters indicate an actual directory
or file name, and lower-case letters indicate the general form of a set
of directory or file names.

FILE/DIRECTORY                     DESCRIPTION

Top-level directory
   |
   |- AAREADME.TXT  Describes volume contents, organization and disk
   |                use (this file).
   |
   |- ERRATA.TXT    Cumulative listing of known errors and anomalies in
   |                the volumes.
   |
   |- VOLDESC.CAT   Description of the contents of this volume.
   |
   |- [INDEX]            Contains index tables for the data products on
   |    |                this disk.
   |    |
   |    |- INDXINFO.TXT  Description of the contents of this directory.
   |    |
   |    |- INDEX.TAB     Index table file for Cassini Radar Instrument Team
   |    |                data products.
   |    |
   |    |- INDEX.LBL     Detached label file describing contents of INDEX.TAB.
   |
   |- [DOCUMENT]                 Contains various documents that facilitate
   |    |                        the understanding of the datasets included
   |    |                        on this volume.
   |    |
   |    |- DOCINFO.TXT           Description of the contents of this
   |    |                        directory.
   |    |
   |    |- BIDRSIS.CSS           Cascading Style Sheet for the HTML version
   |    |                        of the BIDR SIS.
   |    |
   |    |- BIDRSIS.HTML          HTML version of the BIDR SIS.
   |    |
   |    |- BIDRSIS.LBL           Detached label describing the files that
   |    |                        comprise the BIDR SIS.
   |    |
   |    |- BIDRSIS.PDF           PDF version of the BIDR (Basic Image Data
   |    |                        Records) SIS.
   |    |
   |    |- BIDRSIS_FIGnn.GIF     Embedded figures for BIDRSIS.HTML, where nn
   |    |                        is the 2-digit figure number.
   |    |
   |    |- BODPSIS.HTML          HTML version of the BODP SIS.
   |    |
   |    |- BODPSIS.LBL           Detached label describing the files that
   |    |                        comprise the BODP SIS.
   |    |
   |    |- BODPSIS.PDF           PDF version of the BODP (Burst Ordered Data
   |    |                        Products) SIS.
   |    |
   |    |- BODPSIS_FIGnn.GIF     Embedded figures in GIF format for
   |    |                        BODPSIS.HTML, where nn is the 2-digit figure
   |    |                        number.
   |    |
   |    |- BODPSIS_FIGnn.JPG     Embedded figures in JPEG format for
   |    |                        BODPSIS.HTML, where nn is the 2-digit figure
   |    |                        number.
   |    |
   |    |- VOLSIS.CSS            Cascading Style Sheet for the HTML version
   |    |                        of the Volume SIS.
   |    |
   |    |- VOLSIS.HTML           HTML version of the Volume SIS.
   |    |
   |    |- VOLSIS.LBL            Detached label describing the files that
   |    |
   |    |- VOLSIS.PDF            PDF version of the Volume SIS.
   |
   |- [CATALOG]     This directory contains the catalog files.
   |    |
   |    |- CATINFO.TXT   Description of the contents of this directory.
   |    |
   |    |- ABDRDS.CAT    Information about the ABDR data sets.
   |    |                If no altimetry data was collected during the
   |    |                observation, then this catalog will not be
   |    |                present on the volume.
   |    |
   |    |- ABDRSUMDS.CAT Information about the ABDR-SUMMARY data sets.
   |    |                If no altimetry data was collected during the
   |    |                observation, then this catalog will not be
   |    |                present on the volume.
   |    |
   |    |- BIDRDS.CAT    Information about the BIDR data sets.
   |    |                If no SAR data was collected during the
   |    |                observation, then this catalog will not be
   |    |                present on the volume.
   |    |
   |    |- DSMAP.CAT     Information about the map projection used for
   |    |                the BIDR products; includes definition of oblique
   |    |                cylindrical projection.
   |    |
   |    |- INST.CAT      Information about the Cassini radar instrument.
   |    |
   |    |- INSTHOST.CAT  Information about the Cassini spacecraft.
   |    |
   |    |- LBDRDS.CAT    Information about the LBDR data sets.
   |    |                If no active-mode data was collected during the
   |    |                observation, then this catalog will not be
   |    |                present on the volume.
   |    |
   |    |- MISSION.CAT   Information about the Cassini mission.
   |    |
   |    |- PERSON.CAT    Information about Cassini team members responsible
   |    |                for generating the radar products.
   |    |
   |    |- PROJREF.CAT   Cassini Project references mentioned in other
   |    |                catalog files.
   |    |
   |    |- REF.CAT       Cassini Radar references mentioned in other catalog
   |    |                files.
   |    |
   |    |- STDRDS.CAT Information about the SARTopo data sets.
   |    |
   |    |- SBDRDS.CAT    Information about the SBDR data sets.
   |    |
   |    |- SOFTWARE.CAT  Information about software that is necessary or
   |    |                 helpful in accessing or processing the data on
   |    |                 this volume.
   |    |
   |    |- TARGETALBIORIX.CAT  Information about Cassini's target Albiorix.
   |    |
   |    |- TARGETABESTLA.CAT  Information about Cassini's target Bestla.
   |    |
   |    |- TARGETERRIAPUS.CAT  Information about Cassini's target Erriapus.
   |    |
   |    |- TARGETKIVIUQ.CAT  Information about Cassini's target Kiviuq.
   |    |
   |    |- TARGETPAALIAQ.CAT  Information about Cassini's target Paaliaq.
   |    |
   |    |- TARGETSIARNAQ.CAT  Information about Cassini's target Siarnaq.
   |    |
   |    |- TARGETSKATHI.CAT  Information about Cassini's target Skathi.
   |    |
   |    |- TARGETTARVOS.CAT  Information about Cassini's target Tarvos.
   |    |
   |    |- TARGETYMIR.CAT  Information about Cassini's target Ymir.
   |    |
   |    |- TARGETHYROKKIN.CAT Information about Cassini's target Hyrokkin.
   |    |
   |- [DATA]        This directory contains the data products, organized
   |    |           into subdirectories by product type.
   |    |
   |    |- [ABDR]   Directory containing ABDR data product files.
   |    |    |      If no altimetry data was collected during the
   |    |    |      observation, then this directory will not be
   |    |    |      present on the volume.
   |    |    |
   |    |    |- ABDR_yy_Dzzz_Vnn.ZIP   ZIP-compressed ABDR product file,
   |    |    |                         where
   |    |    |                         yy  = Radar mode of data
   |    |    |                         zzz = 3-digit radar observation counter
   |    |    |                         nn  = 2-digit version number
   |    |    |
   |    |    |- ABDR_yy_Dzzz_Vnn.LBL   Detached label for ABDR_yy_Dzzz_Vnn.ZIP
   |    |    |
   |    |    |- ABDR.FMT               Format file for the ABDR_TABLE object
   |    |                              referenced in the attached label for
   |    |                              an Altimeter Burst Data Record product.
   |    |
   |    |- [ABDR_SUMMARY]   Directory containing ABDR-SUMMARY data product
   |    |    |              files.  If no altimetry data was collected during
   |    |    |              the observation, then this directory will not be
   |    |    |              present on the volume.
   |    |    |
   |    |    |- ABDR_SUMMARY_yy_Dzzz_Vnn.CSV
   |    |    |              ABDR-SUMMARY product file,
   |    |    |              where
   |    |    |              yy  = Radar mode of data
   |    |    |              zzz = 3-digit radar observation counter
   |    |    |              nn  = 2-digit version number
   |    |    |
   |    |    |- ABDR_SUMMARY_yy_Dzzz_Vnn.LBL
   |    |    |              Detached label for ABDR_SUMMARY_yy_Dzzz_Vnn.CSV
   |    |    |
   |    |    |- ABDR_SUMMARY.FMT  Format file for the ABDR-SUMMARY product.
   |    |    |
   |    |- [LBDR]   Directory containing LBDR data product files.
   |    |    |      If no active-mode data was collected during the
   |    |    |      observation, then this directory will not be
   |    |    |      present on the volume.
   |    |    |
   |    |    |- LBDR_yy_Dzzz_Pm_Vnn.ZIP   ZIP-compressed LBDR product file,
   |    |    |                            where
   |    |    |                            yy  = Radar mode of data
   |    |    |                            zzz = 3-digit radar observation
   |    |    |                                  counter
   |    |    |                            m   = 1 or 2 for the two parts of a
   |    |    |                                  LBDR product whose total size
   |    |    |                                  is greater than or equal to 2
   |    |    |                                  GB.  If the original LBDR is
   |    |    |                                  smaller than 2 GB, then the
   |    |    |                                  "Pm_" is omitted from the
   |    |    |                                  filename.
   |    |    |                            nn  = 2-digit version number
   |    |    |
   |    |    |- LBDR_yy_Dzzz_Pm_Vnn.LBL   Detached label for
   |    |    |                            LBDR_yy_Dzzz_Pm_Vnn.ZIP.
   |    |    |
   |    |    |- LBDR_yy_Dzzz_Vnn.ZIP      ZIP-compressed LBDR product file
   |    |    |                            for a LBDR product whose total size
   |    |    |                            does not exceed 2 GB, where
   |    |    |                            yy  = Radar mode of data
   |    |    |                            zzz = 3-digit radar observation
   |    |    |                                  counter
   |    |    |                            nn  = 2-digit version number
   |    |    |
   |    |    |- LBDR_yy_Dzzz_Vnn.LBL      Detached label for
   |    |    |                            LBDR_yy_Dzzz_Vnn.ZIP.
   |    |    |
   |    |    |- LBDR.FMT               Format file for the LBDR_TABLE object
   |    |                              referenced in the attached label for
   |    |                              a Long Burst Data Record product.
   |    |
   |    |- [SBDR]   Directory containing SBDR data product files.
   |    |    |
   |    |    |- SBDR_yy_Dzzz_Vnn.TAB   SBDR product file, where
   |    |    |                         yy  = Radar mode of data
   |    |    |                         zzz = 3-digit radar observation counter
   |    |    |                         nn  = 2-digit version number
   |    |    |
   |    |    |- SBDR.FMT               Format file for the SBDR_TABLE object
   |    |                              referenced in the attached label for
   |    |                              a Short Burst Data Record product.
   |    |
   |    |- [BIDR]   Directory containing BIDR data product files.
   |    |    |      If no SAR data was collected during the
   |    |    |      observation, then this directory will not be
   |    |    |      present on the volume.
   |    |    |
   |    |    |- BIbcdeefggg_Dhhh_TiiiSjj_Vnn.ZIP    Zip-compressed BIDR
   |    |    |                                product file, where
   |    |    |                        b     = Kind and bit-type of data:
   |    |    |                            F = Primary dataset:  sigma0 with
   |    |    |                                correction for incidence angle
   |    |    |                                effects and noise subtracted,
   |    |    |                                in 32-bit floating-point format
   |    |    |                                (linear scale, not dB)
   |    |    |                            B = Primary dataset in unsigned
   |    |    |                                byte format (dB, normalized
   |    |    |                                to [0, 255])
   |    |    |                            D = Standard deviation of sigma0 
   |    |    |                                with noise subtraction but
   |    |    |                                without correction for incidence
   |    |    |                                angle, in 32-bit floating-point 
   |    |    |                                format (linear scale, not dB)   
   |    |    |                            S = Sigma0 with noise subtracted
   |    |    |                                but without correction for
   |    |    |                                incidence angle effects, in
   |    |    |                                32-bit floating-point format
   |    |    |                                (linear scale, not dB)
   |    |    |                            U = Sigma0 without either
   |    |    |                                correction for incidence angle
   |    |    |                                effects or noise subtraction,
   |    |    |                                in 32-bit floating-point format
   |    |    |                                (linear scale, not dB)
   |    |    |                            X = Noise Equivalent sigma0 with 
   |    |    |                                noise subtraction but without
   |    |    |                                correction for incidence angle
   |    |    |                                effects, in 32-bit floating-
   |    |    |                                point format (linear scale, 
   |    |    |                                not dB)
   |    |    |                            E = Incidence angle map, 32-bit
   |    |    |                                floating-point values in degrees
   |    |    |                            T = Latitude map, 32-bit floating-
   |    |    |                                point values in degrees
   |    |    |                            N = Longitude map, 32-bit floating-
   |    |    |                                point values in degrees
   |    |    |                            M = Beam mask map, 8-bit unsigned
   |    |    |                                byte values
   |    |    |                            L = Number of looks map, 32-bit
   |    |    |                                integer values
   |    |    |                        c     = Map projection:
   |    |    |                            Q = Oblique cylindrical
   |    |    |                        d     = Map resolution:
   |    |    |                            B =   2 pixels/degree
   |    |    |                            C =   4 pixels/degree
   |    |    |                            D =   8 pixels/degree
   |    |    |                            E =  16 pixels/degree
   |    |    |                            F =  32 pixels/degree
   |    |    |                            G =  64 pixels/degree
   |    |    |                            H = 128 pixels/degree
   |    |    |                            I = 256 pixels/degree
   |    |    |                        ee    = Absolute value of latitude at
   |    |    |                                center of file, rounded to
   |    |    |                                nearest degree
   |    |    |                        f     = Hemisphere of center of file:
   |    |    |                            N = Northern
   |    |    |                            S = Southern
   |    |    |                        ggg   = West longitude at center of
   |    |    |                                file, rounded to nearest degree
   |    |    |                        hhh   = 3-digit data take ID
   |    |    |                                (observation counter) from which
   |    |    |                                data are included, left-padded
   |    |    |                                with zeroes as needed
   |    |    |                        iii   = 3-digit index of Titan flyby
   |    |    |                                from which data are included,
   |    |    |                                left-padded with zeroes as
   |    |    |                                needed.  (Note:  the index
   |    |    |                                for the "Ta" flyby is "00A")
   |    |    |                        jj    = 2-digit segment number (00-99)
   |    |    |                        nn    = 2-digit version number (01-99)
   |    |    |
   |    |    |- BIbcdeefggg_Dhhh_TiiiSjj_Vnn.LBL    Detached label for
   |    |    |                                 BIbcdeefggg_Dhhh_Tiii_Vnn.ZIP
   |    |
   |    |- [STDR] Directory containing SARTopo product files.
   |    |
   |    |    |- SARTOPO_TaaaSbb_Bcc_Vvv_yymmdd.CSV
   |    |    |                      aaa     = Titan flyby number
   |    |    |                      bb      = SAR imaging segment number
   |    |    |                      cc      = Overlapping beam numbers used
   |    |    |                                to produce SARTopo
   |    |    |                      vv      = Version number 
   |    |    |                      yymmdd  = Year, Month, Day of SARTopo 
   |    |    |                                file creation  
   |    |    |
   |    |    |- SARTOPO_TaaaSbb_Bcc_Vvv_yymmdd.LBL   Detached label for
   |    |    |                            SARTOPO_TaaaSbb_Bcc_Vvv_yymmdd.CSV
   |    |    |
   |    |    |- STDR.FMT  Format file for the SARTOPO product.
   |    |    
   |- [CALIB]           This directory contains files that were used to
   |    |               calibrate the Cassini radar data.
   |    |
   |    |- CALINFO.TXT  Description of the contents of this directory.
   |    |
   |    |- [BEAMPAT]    Directory containing antenna gain pattern files
   |    |    |          for each beam and documentation for those files.
   |    |    |
   |    |    |- PATINFO.TXT    Description of the contents of this directory.
   |    |    |
   |    |    |- BEAMx_Vnn.PAT  Antenna gain pattern file, where
   |    |    |                 x  = Beam number (1 to 5)
   |    |    |                 nn = 2-digit version number
   |    |    |
   |    |    |- BEAMx_Vnn.LBL  Detached label for BEAMx_Vnn.PAT.
   |    |
   |    |- [CONFIG]     Directory containing processor configuration files
   |         |          and documentation for those files.
   |         |
   |         |- CFGINFO.TXT       Description of the contents of this
   |         |                    directory.
   |         |
   |         |- PCF_yyyy_Vnn.CFG  Processor configuration file for this
   |         |                    observation, where
   |         |                    yyyy = 4-digit radar observation counter
   |         |                    nn   = 2-digit version number
   |         |
   |         |- PCF_yyyy_Vnn.LBL  Detached label for PCF_yyyy_Vnn.CFG.
   |
   |- [EXTRAS]      This directory contains files that do not fit anywhere
   |    |           else.  A Cassini Radar Transition file should be present
   |    |           for each observation.  The directory may also contain
   |    |           memorandums, plots, or analyses pertinent to the data
   |    |           products being delivered.
   |    |
   |    |- CRT_yyy_Vnn.TAB        Cassini Radar Transition file, where
   |    |                         yyy = 3-digit radar observation counter
   |    |                         nn  = 2-digit version number
   |    |
   |    |- xx_SEQUENCE_DESIGN_MEMO.PDF  Memorandum (in PDF format) describing
   |    |                               command sequence design for
   |    |                               observation, where
   |    |                               xx = Observation ID (e.g., "TA"),
   |    |                                    expressed as a combination of
   |    |                                    either two or three upper-case
   |    |                                    letters and/or numerals.
   |    |
   |    |- TITAN_FLYBY_ID_TABLE.TXT    This file contains a cross reference
   |    |                              table of the four different 
   |    |                              identifying numbers for each Titan 
   |    |                              flyby in the Prime Mission. The four 
   |    |                              columns from left to right are: 
   |    |                              Titan flyby number, data take number,
   |    |                              sequence number, and rev number.
   |    |
   |    |- CPCK11FEB2008_SCI.TPC  SPICE planetary constants kernel with
   |                              short term Titan spin model relevant
   |                              for Titan passes up to T30. The spin
   |                              model was determined from Cassini RADAR
   |                              data.
   |
   |
   |- [SOFTWARE]         This directory contains information about software
        |                applications that are necessary or helpful in
        |                accessing or processing the data products on this
        |                volume.
        |
        |- SOFTINFO.TXT  Description of the contents of this directory,
                         including information on where to obtain the
                         software that was used to compress the LBDR and
                         (if present) the BIDR data products on this volume.


5.  Recommended DVD Drives and Driver Software

Not applicable for electronic delivery.


6.  Errata and Disclaimer

A list of anomalies and errors is maintained in the file ERRATA.TXT
at the root directory of this volume.

Although considerable care has gone into making this volume, errors
are both possible and likely.  Users of the data are advised to exercise
the same caution as they would when dealing with any other unknown data
set.

Reports of errors or difficulties would be appreciated.  Please contact
one of the persons listed herein.



7. Version Status

Ta-T30
Flybys were initially processed with the nominal International 
Astronomical Union (IAU) Titan spin model. Data products were submitted
to PDS with a PRODUCT_VERSION_ID 01. These submissions did not include 
SARTopo products nor the noise equivalent sigma0 & sigma0 standard 
deviation backplanes.

Flybys were later reprocessed with a non-synchronous rotation Titan
spin model (pck = cpck11Feb2008_Sci.tpc). Data products were submitted
to PDS with a PRODUCT_VERSION_ID 02. These submissions did not include
SARTopo products nor the noise equivalent sigma0 & sigma0 standard 
deviation of sigma0 backplanes.

Flybys were or soon will be reprocessed with a synchronous rotation 
Titan spin model (pck=cpck25Aug2011.tpc). Data products were or will be
submitted to PDS with a PRODUCT_VERSION_ID 03. These submissions will
include SARTopo products and the noise equivalent sigma0 & sigma0
standard deviation backplanes.


T36-T71
Flybys were initially processed with a non-synchronous rotation Titan 
spin model. Data products were submitted to PDS with a
PRODUCT_VERSION_ID 02. These submissions did not include SARTopo 
products nor the noise equivalent sigma0 & sigma0 standard deviation 
backplanes.

Flybys were later or soon will be reprocessed with a synchronous
rotation Titan spin model (pck=cpck25Aug2011.tpc). Data products
were or will be submitted to PDS with a PRODUCT_VERSION_ID 03. 
These submissions will include SARTopo products and the noise 
equivalent sigma0 & sigma0 standard deviation backplanes.

T77 and T104
Flybys were processed with a synchronous rotation Titan spin model.
Data products were submitted to PDS with a PRODUCT_VERSION_ID 02.
These submissions did not include SARTopo products nor the noise 
equivalent sigma0 & sigma0 standard deviation backplanes.

In order to include sigma0 & standard deviation of sigma0 backplanes,
the flybys were rerun through the processor to make a complete set
of backplanes. Data products were or will be submitted to PDS with a
PRODUCT_VERSION_ID 03. These submissions will include SARTopo products 
and the noise equivalent sigma0 & sigma0 standard deviation backplanes.

T108 and above
Flybys were processed with a synchronous rotation Titan spin model.
Data products were submitted to PDS with a PRODUCT_VERSION_ID 02. 
These submissions already include SARTopo products and the noise 
equivalent sigma0 & sigma0 standard deviation backplanes.

8.  Contact Information

For questions concerning this volume set, data products, and documentation,
contact

    Philip S. Callahan
    Jet Propulsion Laboratory
    M/S 300-319
    4800 Oak Grove Drive
    Pasadena, CA 91109-8099 USA
    Phone:  +1 818-354-4753
    E-mail:  Philip.S.Callahan@jpl.nasa.gov

    Bryan W. Stiles
    Jet Propulsion Laboratory
    M/S 300-319
    4800 Oak Grove Drive
    Pasadena, CA 91109-8099 USA
    Phone:  +1 818-354-5329
    E-mail:  Bryan.W.Stiles@jpl.nasa.gov

For questions concerning PDS standards and usage, contact

    Patricia A. Garcia
    United States Geological Survey
    2255 North Gemini Drive
    Flagstaff, AZ 86001
    Phone:  (928) 556-7246
    E-mail:  pgarcia@usgs.gov

