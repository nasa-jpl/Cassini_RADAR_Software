#! /bin/csh
# Author: R. Glenister
set filenew=filewithnoname
if( $#argv < 3 ) then
  echo "Usage:  $0 barestring1 barestring2 files"
  echo "        Substitutes barestring2 for barestring1 in files"
  exit
endif

# GET A GOOD DELIMITER (GOOD =  ONE NOT IN THE STRING )
#if( `expr "$1$2" : '.*&'` != 0 )then
#  echo "$0:  Will choke on symbol '&'"
#endif
if( `expr "$1$2" : '.*/'` == 0 )then
  set delim=/
else if( `expr "$1$2" : '.*+'` == 0 )then
  set delim=+
else if( `expr "$1$2" : '.*^'` == 0 )then
  set delim=^
else if( `expr "$1$2" : '.*@'` == 0 )then
  set delim=@
else if( `expr "$1$2" : '.*z'` == 0 )then
  set delim=z
else
  echo "$0 :  Can't find a good delimiter for the sed function."
  exit
endif
set global_flag=g
set iarg=2
while ($iarg<$#argv)
@ iarg++
  set file=$argv[$iarg]
# echo $iarg $file
  cat $file | sed "s$delim$1$delim$2$delim$global_flag" >! $filenew
  set sv=$status
  if( $sv == 0 )then
    \mv $filenew $file 
  else
    "$0: Bad exit status.  Took no action."
    rm -f $filenew
  endif
end
