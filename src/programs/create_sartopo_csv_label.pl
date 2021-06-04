#!/usr/bin/perl -w
#
# Usage: create_sartopo_csv_label.pl --output=<out> --bidr=<bidr> --template=<template>  --csv=<csv> 
#
#
# create_sartopo_csv_label.pl . -output="SARTOPO_XXXXX.LBL" -bidr="BIXXXXXXXX" -template="STDR.FMT"  -csv="SARTOPO_XXXXXXX.csv"

use strict;
use File::Basename;
use Getopt::Long;
use File::stat;
use POSIX qw(strftime);

##-----------------------
## Parse the command line
##-----------------------

my $label;
my $bidr;
my $template;
my $fileCSV;

GetOptions("output=s"   => \$label,
           "bidr=s"     => \$bidr,
           "template=s" => \$template,
           "csv=s"      => \$fileCSV);

if (@ARGV != 0) {
    usage();
    exit 1;
}

if (!defined($label)) {
    print STDERR "Error:  output label must be specified via --output\n";
    usage();
    exit 1;
}

if (!defined($bidr)) {
    print STDERR "Error:  BIDR file must be specified via --bidr\n";
    usage();
    exit 1;
}

if (!defined($template)) {
    print STDERR "Error:  STDR label template must be specified via --template\n";
    usage();
    exit 1;
}

if (!defined($fileCSV)) {
    print STDERR "Error: STDR summary file must be specified via --csv\n";
    usage();
    exit 1;
}

my %keyword_updates;
my %keyword_update_spaces;

my $TMPFILE_SUFFIX = ".tmp.$$";
my $label_tmp = $label . $TMPFILE_SUFFIX;

$\ = "\r\n";

#clear file if already existing
open(DAT, ">$label");
close DAT;

getKeyword("PDS_VERSION_ID", $bidr, $label);

# FILE FORMATS AND LENGTH
open(DAT, ">>$label");
print DAT "";
print DAT "/*       FILE FORMATS AND LENGTH        */";
print DAT "";
close DAT;

getKeyword("RECORD_TYPE", $bidr, $label);
getKeyword("RECORD_BYTES", $bidr, $label);
getKeyword("FILE_RECORDS", $bidr, $label);

# PRODUCT DESCRIPTION   
open(DAT, ">>$label");
print DAT "";
print DAT "/*       PRODUCT DESCRIPTION     */";
print DAT "";
close DAT;
 
getKeyword("DATA_SET_ID", $bidr, $label);
getKeyword("DATA_SET_NAME", $bidr, $label);
my $data_set_id = getPDSFileScalarValue("DATA_SET_ID", $bidr);
my $data_set_name = getPDSFileScalarValue("DATA_SET_NAME", $bidr);
getKeyword("PRODUCER_INSTITUTION_NAME", $bidr, $label);
getKeyword("PRODUCER_ID", $bidr, $label);
getKeyword("PRODUCER_FULL_NAME", $bidr, $label);
getKeyword("PRODUCT_ID", $bidr, $label); 
my $product_id = getPDSFileScalarValue("PRODUCT_ID", $bidr);
getKeyword("PRODUCT_VERSION_ID", $bidr, $label);           
getKeyword("INSTRUMENT_HOST_NAME", $bidr, $label);        
getKeyword("INSTRUMENT_HOST_ID", $bidr, $label);          
getKeyword("INSTRUMENT_NAME", $bidr, $label);             
getKeyword("INSTRUMENT_ID", $bidr, $label);               
getKeyword("TARGET_NAME", $bidr, $label);                  
getKeyword("START_TIME", $bidr, $label);                 
getKeyword("STOP_TIME", $bidr, $label);                    
getKeyword("SPACECRAFT_CLOCK_START_COUNT", $bidr, $label); 
getKeyword("SPACECRAFT_CLOCK_STOP_COUNT", $bidr, $label);  
getKeyword("PRODUCT_CREATION_TIME", $bidr, $label);      
getKeyword("SOURCE_PRODUCT_ID", $bidr, $label);  
getKeyword("MISSION_NAME", $bidr, $label);                 
getKeyword("SOFTWARE_VERSION_ID", $bidr, $label);     
getKeyword("MAXIMUM_LATITUDE", $bidr, $label); 
getKeyword("MINIMUM_LATITUDE", $bidr, $label); 
getKeyword("EASTERNMOST_LONGITUDE", $bidr, $label); 
getKeyword("WESTERNMOST_LONGITUDE", $bidr, $label); 
getKeyword("LOOK_DIRECTION", $bidr, $label); 

my $max_lat = getPDSFileScalarValue("MAXIMUM_LATITUDE", $bidr);
my $min_lat = getPDSFileScalarValue("MINIMUM_LATITUDE", $bidr);
my $east_lon = getPDSFileScalarValue("EASTERNMOST_LONGITUDE", $bidr);
my $west_lon = getPDSFileScalarValue("WESTERNMOST_LONGITUDE", $bidr);
my $look_dir = getPDSFileScalarValue("LOOK_DIRECTION", $bidr);

$keyword_update_spaces{"MAXIMUM_LATITUDE"} = $max_lat;
$keyword_update_spaces{"MINIMUM_LATITUDE"} = $min_lat;
$keyword_update_spaces{"EASTERNMOST_LONGITUDE"} = $east_lon;
$keyword_update_spaces{"WESTERNMOST_LONGITUDE"} = $west_lon;
$keyword_update_spaces{"LOOK_DIRECTION"} = $look_dir;

removeKeywordsLeadingSpace($label, $label_tmp, %keyword_update_spaces);

     
getKeyword("DESCRIPTION", $bidr, $label);         
my $description = getPDSFileScalarValue("DESCRIPTION", $bidr);
getKeyword("PROCESSING_HISTORY_TEXT", $bidr, $label);  

# POINTERS TO START RECORDS OF OBJECTS IN FILE
open(DAT, ">>$label");
print DAT "";
print DAT "/*        POINTERS TO START RECORDS OF OBJECTS IN FILE */";
print DAT "";
print DAT "^SPREADSHEET                 = \"" . basename($fileCSV) . "\"";
print DAT "";
close DAT;

# Calculate the (theoretical) maximum length of each row in the STDR-CSV
# file from the field widths given in the template for the SPREADSHEET
# object definition.  The calculation also assumes that any field with a
# DATA_TYPE of CHARACTER (or EBCDIC_CHARACTER) is enclosed in double quotes.

open(FILE, $template) or
    die "Cannot open $template for reading:  $!, exiting";
my $terminator_orig = $/;
$/ = "\r\n";

my $field_count = 0;
my $max_row_length = 0;
my $character_fields = 0;
while (<FILE>)
{
    chomp;
    if (/^[\s]*OBJECT[\s]+=[\s]+FIELD[\s]*$/)
    {
        $field_count++;
    }
    if (/^[\s]*BYTES[\s]+=[\s]+(.*)[\s]*$/)
    {
        $max_row_length += $1;
    }
    if (/^[\s]*DATA_TYPE[\s]+=[\s]+.*CHARACTER[\s]*$/)
    {
        $character_fields++;
    }
}
$max_row_length += ($field_count - 1); # include the commas
$max_row_length += length($/); # and the record terminator
$max_row_length += 2 * $character_fields; # and two double quotes per character field

$/ = $terminator_orig;
close FILE;


# Update keywords in the label, which has not yet been concatenated to
# the template for the SPREADSHEET object definition.  Each file is being
# updated separately because the template also contains DESCRIPTION
# parameters that should not be modified.  (Currently, updateKeywords()
# modifies all instances of the specified keyword(s) in a file.)

$keyword_updates{"RECORD_TYPE"} = "STREAM";

$terminator_orig = $/;
$/ = "\r\n";
my $recSize = maxRecordSize($fileCSV);
if ($recSize==156){
    $keyword_updates{"RECORD_BYTES"} = $recSize;}
else{
    die "RECORD_BYTES is not 156 bytes. It is $recSize bytes. Check SARTOPO CSV files for extra characters at end of row\n";
}

my $nRows = numRows($fileCSV);
$/ = $terminator_orig;
$keyword_updates{"FILE_RECORDS"} = $nRows ;

$data_set_id =~ s/BIDR/STDR/;
$keyword_updates{"DATA_SET_ID"} = $data_set_id;
$data_set_name =~
    s/BIDR/SARTOPO/;
$keyword_updates{"DATA_SET_NAME"} = $data_set_name;
#$product_id =~ s/ABDR_04/SARTOPO_05/;
#$keyword_updates{"PRODUCT_ID"} = $product_id;

#$source_product_id = $bidr;
#$source_product_id =~ s/.IMG//;
$keyword_updates{"SOURCE_PRODUCT_ID"} = $product_id;

my $basename_CSV = basename($fileCSV);
$basename_CSV =~ s/\..*//; # removes filenname extension
$product_id =~ s/$product_id/$basename_CSV/;
$keyword_updates{"PRODUCT_ID"} = $product_id;

my $product_creation_time = GetFileCreationTime($fileCSV);
$keyword_updates{"PRODUCT_CREATION_TIME"} = $product_creation_time;

#$description =~
#    s/ALTIMETER BURST DATA RECORD/SARTOPO/;
#$keyword_updates{"DESCRIPTION"} = $description;

#$product_id =~ s/_P1_/_/;
#$keyword_updates{"PRODUCT_ID"} = $product_id;

#$product_id =~ s/_P2_/_/;
#$keyword_updates{"PRODUCT_ID"} = $product_id;


updateKeywords($label, $label_tmp, %keyword_updates);

# Update keywords in the template

undef %keyword_updates;
my $template_tmp = $template . $TMPFILE_SUFFIX;

$keyword_updates{"ROWS"} = $nRows ;
$keyword_updates{"ROW_BYTES"} = $max_row_length;
$keyword_updates{"FIELDS"} = $field_count;

updateKeywords($template, $template_tmp, %keyword_updates);

# Append the template to the preliminary label to form the final detached label

open(FILE, $template) or
    die "Cannot open $template for reading:  $!, exiting";
open(DAT, ">>$label") or
    die "Cannot open $label for reading:  $!, exiting";
$terminator_orig = $/;
$/ = "\r\n";

while (<FILE>)
{
    chomp;
    print DAT;
}

$/ = $terminator_orig;
close FILE;
close DAT;

# We're done
exit 0;

#----------------------------------------------------------------------
# getKeyword($keyword, $bidr, $label)
#
# Read the value of a PDS keyword from the given file and return the
# value. The first argument is the keyword that needs to be searched 
# for. The second argument is the name of the bidr file. The third
# argument is the name of the output label file.
#
#----------------------------------------------------------------------
sub getKeyword 
{
    my $keyword = shift;
    my $filename = shift;
    my $label = shift;
    my $field1;
    my $field2;
    my $field3;
    my $field4;
    my $line_terminator_orig = $\;

    $\ = "\r\n";
    open(FILE, $filename) or
        die "Cannot open $filename for reading:  $!, exiting";
    open(DAT, ">> $label") or
        die "Cannot open $label for reading:  $!, exiting";
    while (<FILE>)
    {
	($field1, $field4) = split '=';
	($field2) = split ' ', $field1;
        ($field3) = split '\r';
	if (defined($field2)) {
	    if ($field2 eq $keyword)
	    {
		print DAT "$field3";
		last;
	    }
	}
    }
    close FILE;
    close DAT;
    $\ = $line_terminator_orig;
    return $field4;
}

#----------------------------------------------------------------------
# updateKeywords($file, $file_tmp, %keyword_updates)
#
# Update the values of specified PDS keywords in a file.  The first
# argument is the original file, which should be a PDS text file or
# detached label.  The second argument is the name of a temporary file
# to use for the output.  The third argument is a hash whose keys are
# the PDS keywords to update and whose values are the new values of
# those keywords.  The original file is replaced by the temporary file
# after all the lines in the file have been read.
#
# The function returns the number of keywords updated.
#
# Note:  no warnings are printed for keywords that are not in the file.
#----------------------------------------------------------------------
sub updateKeywords 
{
    my $file = shift;
    my $file_tmp = shift;
    my %keyword_updates = @_ ;
    my $PLACEHOLDER = "XXXX";
    my $terminator_orig = $/;
  
    $/ = "\r\n";
    open IN, "<$file" or
	die "Cannot open $file for reading:  $!, exiting";
    open OUT, ">$file_tmp" or
	die "Cannot open $file_tmp for writing:  $!, exiting";
    my $count = 0;
    while (<IN>)
    {
	chomp;
	my $s = $_;
	foreach my $k (keys %keyword_updates) 
	{
	    if (/^([\s]*)$k([\s]+=[\s]+["])(.*)(["][\s]*.*)$/)
	    {
	        $s = "$1$k$2" . $keyword_updates{$k} . "$4";
		$count++;
	    }
	    elsif (/^([\s]*)$k([\s]+=[\s]+)(.*)([\s]*.*)$/)
	    {
	        $s = "$1$k$2" . $keyword_updates{$k} . "$4";
		$count++;
	    }
#	    my $found = s/$PLACEHOLDER/$keyword_updates{$k}/;
	}
	print OUT $s;
    }
    close IN;
    close OUT;

    rename $file_tmp, $file or
	die "Cannot save $file_tmp as $file:  $!, exiting";
    $/ = $terminator_orig;
    return $count;
}

#----------------------------------------------------------------------
# correctKeywords($file, $file_tmp, $keyReplace, $keyReplaced)
#
# Correct the values of specified PDS keywords in a file.  The first
# argument is the original file, which should be a PDS text file or
# detached label.  The second argument is the name of a temporary file
# to use for the output.  The third argument is a hash whose keys are
# the PDS keywords to correct and whose values are the new values of
# those keywords.  The fourth argument is the correct keyword value. 
# The original file is replaced by the temporary file after all the
# lines in the file have been read.
#
# The function returns the number of keywords updated.
#
# Note:  no warnings are printed for keywords that are not in the file.
#----------------------------------------------------------------------
sub correctKeywords 
{
    my $file = shift;
    my $file_tmp = shift;
    my $keyReplace = shift ;
    my $keyReplaced = shift;
    my $terminator_orig = $/;
  
    $/ = "\r\n";
    open IN, "<$file" or
	die "Cannot open $file for reading:  $!, exiting";
    open OUT, ">$file_tmp" or
	die "Cannot open $file_tmp for writing:  $!, exiting";
  
    while (<IN>)
    {
	chomp;
	my $found = s/$keyReplace/$keyReplaced/;
	print OUT;
    }
    close IN;
    close OUT;
    
    rename $file_tmp, $file or
	die "Cannot save $file_tmp as $file:  $!, exiting";
    $/ = $terminator_orig;
}


#----------------------------------------------------------------------
# maxRecordSize($file)
#
# Determine the size of the longest record(i.e. row) in bytes
# The first argument is the name of the abdr CSV file
# 
#----------------------------------------------------------------------
sub maxRecordSize
{
    my $file =shift;
    open(FILE,$file);
    my $max = 0;
    my $length = 0;
    while (<FILE>)
    {
	$length = length($_);
	if ($length > $max)
	{
	    $max = $length;
	}
    }
    close FILE;
    return $max;
}

#----------------------------------------------------------------------
# numRows($file)
#
# Determine the number of rows in a file CSV file.
# The first argument is the name of the abdr CSV file
# 
#----------------------------------------------------------------------
sub numRows
{
    my $file =shift;
    open(FILE,$file);
    my $nRow = 0;
    while (<FILE>)
    {
	$nRow = $nRow+1;
	
    }
    close FILE;
    return $nRow;
}


#----------------------------------------------------------------------
# getPDSFileScalarValue($keyword, $filename)
#
# Read the value of a PDS keyword from the given file and return the
# value after stripping any opening and closing quotes and units.
# Both string and numerical values are handled, but non-scalar values
# (i.e., sequences and arrays) are not.
#
# See comments below for keywordFound() in regards to handling of a
# leading ^.
#----------------------------------------------------------------------
sub getPDSFileScalarValue {
    my $keyword = shift;
    my $filename = shift;
    open IN, "<$filename" or
        die "Cannot open $filename for reading:  $!, exiting";
    while (<IN>) {
        chomp;
        my $value;
        if ($keyword =~ /^\136(.*)/) {
            # Strip the ^ from the keyword and match the ^ and keyword as
            # separate items.
            $keyword = $1;
            if (/[^\w]*\136$keyword[\s]*=[\s]*(.*)[\s]*\r/) {
                $value = $1;
                # Strip opening or closing quotes, if any
                $value =~ s/^["]//;
                $value =~ s/["]$//;
                # Strip units in < >, if any
                $value =~ s/[\s]*[<].*[>][\s]*//;
                close IN;
                return $value; 
            } 
        } else {
            if (/[^\w]*$keyword[\s]*=[\s]*(.*)[\s]*\r/) {
                $value = $1;
                # Strip opening or closing quotes, if any
                $value =~ s/^["]//;
                $value =~ s/["]$//;
                # Strip units in < >, if any
                $value =~ s/[\s]*[<].*[>][\s]*//;
                close IN;
                return $value;
            }
        }
    }
    close IN;
    return "";
}


#----------------------------------------------------------------------
# keywordFound($keyword, $line)
#
# Return true (1) if a PDS keyword was found on the input, 0 if not.
#
# PDS keywords that begin with a pointer (^, or octal 136) require
# special handling because ^ is a metacharacter in a Perl regular
# expression.
#----------------------------------------------------------------------
sub keywordFound {
    my $keyword = shift;
    my $line = shift;
    chomp $line;
    if ($keyword =~ /^\136(.*)/) {
	# Strip the ^ from the keyword and match the ^ and keyword as
	# separate items.
	$keyword = $1;
	if ($line =~ /[ ]*\136$keyword[ ]*=[ ]*.*[ ]*\r/) {
	    return 1;
	}
    } else {
	if ($line =~ /[ ]*$keyword[ ]*=[ ]*.*[ ]*\r/) {
	    return 1;
	}
    }
    return 0;
}

sub usage {
    print STDERR "Usage:  $0 --bidr=<BIDRfile> --csv=<STDR_CSV> " .
        "--output=<outputLabel> --template=<STDRlabelTemplate>\n";
}


# ----------------------
# Get File Creation Time
# ----------------------
sub GetFileCreationTime{
    my $input_file =shift;
    my $stat_epoch = stat($input_file)->mtime;
    my $time_tag = POSIX::strftime('%Y-%jT%H:%M:%S.000', localtime($stat_epoch)); 
    return $time_tag;
}


# -------------------------------------
# Remove spaces infront of Keyword Name
# -------------------------------------
sub removeKeywordsLeadingSpace 
{
    my $file = shift;
    my $file_tmp = shift;
    my %keyword_updates = @_ ;
    my $terminator_orig = $/;
  
    $/ = "\r\n";
    open IN, "<$file" or
	die "Cannot open $file for reading:  $!, exiting";
    open OUT, ">$file_tmp" or
	die "Cannot open $file_tmp for writing:  $!, exiting";
    my $count = 0;
    while (<IN>)
    {
	chomp;
	my $s = $_;
	foreach my $k (keys %keyword_updates) 
	{
	    if (/^([\s]*)$k([\s]+=[\s]+["])(.*)(["][\s]*.*)$/)
	    {
	        $s = "$k$1$2" . $keyword_updates{$k} . "$4";
		$count++;
	    }
	    elsif (/^([\s]*)$k([\s]+=[\s]+)(.*)([\s]*.*)$/)
	    {
	        $s = "$k$1$2" . $keyword_updates{$k} . "$4";
		$count++;
	    }
	}
	print OUT $s;
    }
    close IN;
    close OUT;

    rename $file_tmp, $file or
	die "Cannot save $file_tmp as $file:  $!, exiting";
    $/ = $terminator_orig;
    return $count;
}
