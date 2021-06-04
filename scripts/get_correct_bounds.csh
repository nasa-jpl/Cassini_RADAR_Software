#!/bin/csh

\rm allbounds.txt

foreach file  ( /home/ras/dat/tour/s*/t*/reproc1/{corramb,seg*}/*bidr*beammask* )

if ( $file:e == 'gz' ) then
gunzip $file
set file = $file:r
endif

if( $file:e != '128' ) then
getbounds.csh $file >> allbounds.txt
endif

gzip $file 
echo $file done
end


