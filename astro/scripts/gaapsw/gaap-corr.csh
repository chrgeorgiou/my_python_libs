#!/bin/tcsh

setenv code  $GAaP/shapelets/kk/
setenv bigim  $code/bigim

# calculate the Gaap photometry from an LDAC catalogue and a GPSF image.
# measure pixel correlation on pre-GPSF image and include it.
# 
# this script also analyses the pre-gaussianization image and incorporates the
# noise correlation it finds there into the gaussianized image's noise ACF.
#
# parameters:
#  $1 catalogue with apertures
#  $2 pre-GPSF image
#  $3 GPSF image
#  $4 weight image
#  $5 minimum aperture size (arcsec, added in quadrature to source size)
#  $6 maximum aperture size (arcsec)
#  $7 kernel map
#  $8 output GaaP file
#
# output files: 
#    A fits image of covariance from pre-GPSF image
#        ($2.cov.fits) and of a shapelet fit to it ($2.covsh.fits)
#    A shapelet map of the total ACF (orig * kernel ACF) $3.totacf.map
#        and a fits image of this ACF $3.totacf.map.fits

if ! ($# == 8) then
 echo Usage: gaap-corr.csh \<cat\> \<pre-GPSFimg\> \<GPSFimg\> \<wtimg\> \<minaper\> \<maxaper\> \<kermap\> \<out\>
 exit
endif

set script = $code/gpsf

set catsky = $1
set origimage = $2
set gpsfimage = $3
set wtimage = $4
set minaper = $5
set maxaper = $6
set kermap = $7
set gaapout = $8

# Check whether the input image exists
 if ! -e $gpsfimage then
  echo ":<(" GPSF image $gpsfimage missing ===
  goto done
 endif

# Format for $catsky aperture catalogue: text file with one line per aperture
# columns X_WORLD Y_WORLD MAG_AUTO MAGERR_AUTO FLUX_RADIUS A_WORLD B_WORLD THETA_WORLD SeqNr 

# Check whether there is an input catalogue.
 if ! -e $catsky then
  echo ":<(" No input aperture catalogue $catsky ===
  goto done
 endif

# Make input catalogue for Gaap photometry: needs X,Y,A",B",PAworld,ID
# NB PAworld uses opposite sign convention (N of W) from THETA_WORLD from SExtractor (N of E)
 awk '{print($1,$2)}' $catsky > radec.cat
 sky2xy  $gpsfimage @radec.cat | awk '{print($5,$6)}' > xy.cat
 paste xy.cat $catsky  | awk -v aper=$minaper -v maxaper=$maxaper '{a=$8*3600;b=$9*3600;a=sqrt(aper**2+a**2);if (a>maxaper) a=maxaper;b=sqrt(aper**2+b**2);if (b>maxaper) b=maxaper;print($1,$2,a,b,-$10,$11)}' > gaap.in

# Make kernel ACF map of GPSF image if needed
 set keracfmap =  `basename $kermap |sed -e 's/ker.map/keracf.map/' `
 if -e $keracfmap then
  echo +++ ALREADY HAVE KERNEL ACF MAP $keracfmap ===
 else
  # check there is a kernel map frame
  if ! -e $kermap then
   echo ":<(" KERNEL MAP FILE $kermap MISSING ===
   goto done
  endif
  echo === MAKING KERNEL ACF MAP $keracfmap ===
  if ! -e orders.par then
    head -2 $kermap | awk '{print($1)}' > orders.par
  endif
  $bigim/keracfmap < $kermap | $bigim/fitkermap > $keracfmap
 endif

# Estimate pre-GPSF pixel correlation and convolve with kernel ACF
# set totacfmap =  `basename $kermap |sed -e 's/ker.map/totacf.map/' `
# if ! -e $totacfmap then
#  ln -sf $origimage inimage.fits
#  ln -sf $wtimage weights.fits
#  echo NOMAP|$bigim/gapphot4 > /dev/null
#  set acfsigorig = `head -1 keracffitted.map|awk '{print($2)}'`
#  set acfsigconv = `head -1 $keracfmap |awk '{print($2)}'`
#  \mv -f cov.fits $origimage.cov.fits
#  \mv -f covsh.fits $origimage.covsh.fits
#  echo === APPROXIMATED NOISE ACF WITH A GAUSSIAN OF WIDTH $acfsigorig
#  (echo $acfsigorig $acfsigconv; cat $keracfmap) |$bigim/kermapxgauss > $totacfmap
#  \rm -f keracffitted.map keracffitted_pknorm.map
# echo === MADE TOTAL NOISE ACF MAP $totacfmap ===
# else
#  echo +++ ALREADY HAVE TOTAL NOISE ACF MAP $totacfmap ===
# endif
set totacfmap = $keracfmap
# Now do gapphot
 ln -sf $gpsfimage inimage.fits
 ln -sf $wtimage weights.fits
 # write GPSFSIG to gpsfsig.dat file in case GPSFSIG keyword not defined (Astrowise writes it to HISTORY kw)
# dfits $gpsfimage| grep GPSFSIG | sed 's/.*GPSFSIG =//'|awk '{print($1)}' > gpsfsig.dat
 (echo $totacfmap ; cat gaap.in) | $bigim/gapphot4 > $gaapout
 \mv -f keracf.fits $totacfmap.fits
 \rm -f gpsfsig.dat inimage.fits weights.fits xy.cat gaap.in radec.cat
 \rm -f bgnoise.dat noise-est.ps





done:

exit