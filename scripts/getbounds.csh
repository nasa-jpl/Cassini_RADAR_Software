#!/bin/csh

set file = $1

clip_header $file head

set flybyseg = `grep PRODUCT_ID head | grep -v SOURCE | awk '{print $3}' | cut -d '_' -f 3`
set bounds = `compute_bidr_bound_box $file`
set minlat = `grep MINIMUM_LATITUDE head | awk '{print $3}'`
set maxlat = `grep MAXIMUM_LATITUDE head | awk '{print $3}'`
set elon = `grep EASTERNMOST_LONGITUDE head | awk '{print $3}'`
set wlon = `grep WESTERNMOST_LONGITUDE head | awk '{print $3}'`
echo "$flybyseg $bounds $minlat $maxlat $elon $wlon" 
\rm head
