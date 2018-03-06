#!/bin/tcsh

setenv code $GAaP/shapelets/kk/
setenv bigim $code/bigim/

# parameter 1 is input fits file name
# parameter 2 is input star catalogue name
# gaussianized output file (gpsf) is {input}.gpsf.fits
# gaussianized & tweaked output file (ggpsf) is {input}.ggpsf.fits

#echo ====================
#echo PSF-GAUSSIANIZING IMAGE $1

set im  = $1
if ! -e $im then
 echo ====  INPUT IMAGE $im NOT FOUND    ==============
 exit
endif

set psf = $2
if ! -e $psf then
 echo ====  PSF CATALOGUE $psf NOT FOUND    ==============
 exit
endif

set psfsh = $im.psf.sh
set psfmap = $im.psf.map
set kersh = $im.ker.sh
set kermap = $im.ker.map

echo 10 > orders.par; echo 4 >> orders.par

if (! -e $psfmap || -z $psfmap) then
# echo BUILDING SHAPELET-BASED PSF MAP
 \rm -f inimage.fits; ln -s $im inimage.fits
 sort -n -k 2 $psf | $bigim/psfcat2sherr_2 > $psfsh
 $bigim/fitpsfmap2 < $psfsh > $psfmap
 \cp -f psfsticks.ps $im.dirtypsfsticks.ps
 \cp -f psfstars.ps $im.dirtypsfstars.ps
 $code/showpsfmap < $psfmap
 \cp -f psfmap.ps $im.dirtypsfmap.ps
else
 echo ====  USING EXISTING PSF MAP $psfmap
endif

exit

if (! -e $kermap || -z $kermap) then
 echo BUILDING SHAPELET-BASED GAUSSIANIZATION KERNEL
 echo -0.8 > ker.in
 \rm -f inimage.fits; ln -s $im inimage.fits
 $code/gaussianizepsf4 < $psfmap > $kersh
 $bigim/fitkermap < $kersh > $kermap
 $code/showpsfmap < $kermap 
 \cp -f psfmap.ps $im.kermap.ps
 \cp -f ker.in gpsfsig.dat
else
 echo ====  USING EXISTING GAUSSIANIZATION KERNEL $kermap
 head -1 $kermap |awk '{print($2)}' > gpsfsig.dat
endif

set gim  = $im.gpsf.fits
if ! -e $gim then
 \rm -f inimage.fits; ln -s $im inimage.fits
 echo BUILDING GAUSSIANIZED IMAGE
 $bigim/imxshmap < $kermap
 \mv -f convolved.fits $gim
else
 echo ====  ALREADY HAVE GAUSSIANIZED IMAGE $gim
endif

set beta2 = `awk '{print($1*1.4)}' gpsfsig.dat`
echo SCALE RADIUS FOR PSF TWEAK IS $beta2

set gpsfsh = $im.gpsf.sh
set gpsfmap = $im.gpsf.map

if (! -e $gpsfmap || -z $gpsfmap) then
 echo COMPUTING SHAPELET-BASED PSF MAP FOR GAUSSIANIZED IMAGE
 \rm -f inimage.fits; ln -s $gim inimage.fits
 (echo $beta2; sort -n -k 2 $psf) | $bigim/psfcat2sherr_fixbeta > $gpsfsh
 #(echo 0 $beta2; sort -n -k 2 $psf) | $bigim/psfcat2sherr_fixbgbeta > $gpsfsh
 $bigim/fitpsfmap2 < $gpsfsh > $gpsfmap
 \cp -f psfsticks.ps $im.gpsfsticks.ps
 \cp -f psfstars.ps $im.gpsfstars.ps
 $code/showpsfmap < $gpsfmap 
 \cp -f psfmap.ps $im.gpsfmap.ps
else
 echo ====  USING EXISTING PSF MAP FROM GAUSSIANIZED IMAGE $gpsfmap
endif

set ggim  = $im.ggpsf.fits
if ! -e $ggim then
 echo TWEAKING GAUSSIANIZATION
 \rm -f inimage.fits; ln -s $gim inimage.fits
 (cat gpsfsig.dat $gpsfmap) | $code/psfmapmingaus > $im.delker.map
 $bigim/imxshmaptweak < $im.delker.map
 \mv -f convolved2.fits $ggim
else
 echo ====  ALREADY HAVE TWEAKED GAUSSIANIZED IMAGE $ggim
endif

exit























if ! -e ggpsf.i.psf.sh then
 echo COMPUTING SHAPELET-BASED PSF MAP FOR GAUSSIANIZED AND TWEAKED IMAGE
 \rm -f inimage.fits; ln -s $ggim inimage.fits
 (echo 0 $beta2; sort -n -k 2 i.psf.cat) | $bigim/psfcat2sherr_fixbgbeta > ggpsf.i.psf.sh
 $bigim/fitpsfmap2 < ggpsf.i.psf.sh > ggpsf.i.psf.map
 \cp -f psfsticks.ps ggpsfsticks.ps
 \cp -f psfstars.ps ggpsfstars.ps
 $code/showpsfmap < ggpsf.i.psf.map 
 \cp -f psfmap.ps ggpsfmap.ps
endif



if ! -e i.gal.g then
 echo COMPUTING SHAPELET BASED SHEARS
 \rm inimage.fits; ln -s $im inimage.fits
 \rm psf.map; ln -s i.psf.map psf.map
 sort -n -k 2 i.gal.cat | $bigim/cat2sherr  | sed -e 's/NaN/1./' > i.gal.sh
 $code/sortbybeta < i.gal.sh | $code/sh2gerrshift012 | sort -n -k 9 > i.gal.g
endif

# THIS TAKES 7hrs PER IMAGE!
if ! -e i.pixpsf.fits then
 echo BUILDING PIXEL-BASED PSF MAP
 echo 6 3 4 > pixmap.par
 \rm -f inimage.fits; ln -s $im inimage.fits
 $bigim/pixpsfmap < i.psf.cat >& pixpsfmap.log
 $bigim/pixpsfxshcpts8 
 \mv pixpsf.fits i.pixpsf.fits
 \mv psfstars.fits i.psfstars.fits
 \mv pixpsfshearcpts8.fits i.pixpsfshearcpts8.fits
endif

if ! -e gpsf.i.pixpsf.fits then
 echo BUILDING PIXEL-BASED PSF MAP FROM GAUSSIANIZED IMAGE
 echo 6 3 4 > pixmap.par
 \rm -f inimage.fits; ln -s $gim inimage.fits
 $bigim/pixpsfmap < i.psf.cat >& gpsf.pixpsfmap.log
 $bigim/pixpsfxshcpts8 
 \mv pixpsf.fits gpsf.i.pixpsf.fits
 \mv psfstars.fits gpsf.i.psfstars.fits
 \mv pixpsfshearcpts8.fits gpsf.i.pixpsfshearcpts8.fits
endif

if ! -e i.gal.gx then
 echo COMPUTING PIXEL-BASED SHEARS
 \rm pixpsfshearcpts8.fits; ln -s i.pixpsfshearcpts8.fits pixpsfshearcpts8.fits
 \rm inimage.fits; ln -s $im inimage.fits
 sort -n -k 2 i.gal.cat | $bigim/pix2g8 > i.gal.gx
endif

if ! -e gpsf.i.gal.g then
 echo COMPUTING GAUSSIANIZED SHAPELET BASED SHEARS
 \rm inimage.fits; ln -s $gim inimage.fits
 sort -n -k 2 i.psf.cat | $bigim/psfcat2sherr > gpsf.i.psf.sh
 $bigim/fitpsfmap2 < gpsf.i.psf.sh > gpsf.i.psf.map
 \cp -f psfsticks.ps gpsfsticks.ps
 \cp -f psfstars.ps gpsfstars.ps
 $code/showpsfmap < gpsf.i.psf.map 
 \cp -f psfmap.ps gpsfmap.ps
 \rm psf.map; ln -s gpsf.i.psf.map psf.map
 sort -n -k 2 i.gal.cat | $bigim/cat2sherr  | sed -e 's/NaN/1./' > gpsf.i.gal.sh
 $code/sortbybeta < gpsf.i.gal.sh | $code/sh2gerrshift012 | sort -n -k 9 > gpsf.i.gal.g
 \rm psf.map inimage.fits
endif

