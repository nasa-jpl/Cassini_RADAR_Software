#!/bin/csh

set file = $1
set dir = $file:h

mkdir $dir/headers
clip_header $file head body
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

\cp head head.old

set eastold = "`grep EASTERNMOST_LONGITUDE head | replace_ctrlM`" 
set westold = "`grep WESTERNMOST_LONGITUDE head | replace_ctrlM`" 
set minlatold = "`grep MINIMUM_LATITUDE head | replace_ctrlM`" 
set maxlatold = "`grep MAXIMUM_LATITUDE head | replace_ctrlM`" 
set seg = `grep PRODUCT_ID head | grep -v SOURCE | awk '{print $3}' | cut -d '_' -f 3`

#echo $seg

set eastnew = `awk '$1=="'$seg'" {print $4}' < ~/reproc/allbounds.txt `
set westnew = `awk '$1=="'$seg'" {print $5}' < ~/reproc/allbounds.txt `
set minlatnew = `awk '$1=="'$seg'" {print $2}' < ~/reproc/allbounds.txt `
set maxlatnew = `awk '$1=="'$seg'" {print $3}' < ~/reproc/allbounds.txt `

#echo $minlatnew
if ( $eastnew != "" ) then
set eastnew = `make_bounds_str "$eastold" "$eastnew" `
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

set westnew = `make_bounds_str "$westold" "$westnew" `
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

set minlatnew = `make_bounds_str "$minlatold" "$minlatnew" `
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

set maxlatnew = `make_bounds_str "$maxlatold" "$maxlatnew" `
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

substring.csh "$eastold" "$eastnew" head 
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif


substring.csh "$westold" "$westnew" head 
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

#echo substring.csh "$minlatold" "$minlatnew" head 
substring.csh "$minlatold" "$minlatnew" head 
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

substring.csh "$maxlatold" "$maxlatnew" head 
if($status) then
echo "Failed to fix_bounds for $file"
exit
endif

substring.csh '@' ' ' head
set sznew = ` wc -c head `
set sznew = $sznew[1]

set szold = ` wc -c head.old `
set szold = $szold[1]



if ( ! -e $dir/headers/${file:t}.head.badbounds ) then
cp head.old $dir/headers/${file:t}.head.badbounds
endif
if( $sznew != $szold ) then
  echo "Header size incorrect error for $file "
  exit
endif

cat head body > tmp
remove_nan tmp $file
\rm tmp
endif


\cp head $dir/headers/${file:t}.head
