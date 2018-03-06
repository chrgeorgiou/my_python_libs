#
#LINUX LIBRARIES FOR veluwe
libs = -L. -lshape -lutil -L ../../lib -lpgplot -lcfitsio -lX11 -O
#libs = -L. -lshape -lutil -L ~/lib -lpgplot -lcfitsio -lX11 -O
# OLD LIBRARIES:
#FC8
#libs = -L. -lshape -lutil -L ~/lib -lpgplot -lcfitsio -L/usr/X11R6/lib -lX11 -O
#
#libs = -L. -lshape -lutil -L/usr/X11R6/lib -lX11 /home/strw/kuijken/lib/pgplot/libpgplot.a -lcfitsio -O
#libs = -L. -lshape -lutil -L/usr/X11R6/lib -L/usr/lib -lX11 -lpgplot -lpng -lcfitsio -O
#
F77 = gfortran
veluwemakefile:
	cat makefile.veluwe makefile.common > makefile
linuxmakefile:
	cat makefile.linux makefile.common > makefile
macmakefile:
	cat makefile.macos makefile.common > makefile
linuxfc10makefile:
	cat makefile.linuxfc10 makefile.common > makefile
eemsmakefile:
	cat makefile.eems makefile.common > makefile


# mm is maximum shapelet order that will be used
# mmpix is maximum size of a fits image that can be analyzed (mmpix x mmpix)
# mmpoly is maximum polynomial order to be used for spatial variation of PSF
#
# Make sure the following derived parameters are consistent with these choices
# msh = (mm+1) * (mm+2) /2 +2
# mmpoly2=(mmpoly+1)*(mmpoly+2)/2
#
flg1 = -Dmm=16 -Dmsh=155 -Dmmpix=9000 -Dmmcat=1000000 -Dmmpoly=6 -Dmmpoly2=28
flg2 = -ffixed-line-length-0 -O -Wall -fbounds-check

.F.o: 
	$(F77) $(flg1) $(flg2) -c  -O $*.F

set: psfcat2sherr_2 fitpsfmap2 fitkermap imxshmap psfcat2sherr_fixbeta psfcat2sherr_fixbgbeta imxshmaptweak gapphot4 keracfmap kermapxgauss

libshape.a: code.o convolve.o fitshear.o operators.o plotshape.o \
            setpsfmat.o setroundshear.o setroundshearshift.o
	ar crs libshape.a code.o convolve.o fitshear.o operators.o plotshape.o\
            setpsfmat.o setroundshear.o setroundshearshift.o; ranlib libshape.a
libutil.a: fits.o mrq.o workhorse.o fitpoly.o readpars.o
	ar crs libutil.a fits.o mrq.o workhorse.o fitpoly.o readpars.o; ranlib libutil.a

testmain: testmain.o libshape.a libutil.a
	$(F77) testmain.o $(libs) -o testmain
testshear: testshear.o  libshape.a libutil.a
	$(F77) testshear.o $(libs) -o testshear
testsh2gerr: testsh2gerr.o libshape.a libutil.a
	     $(F77) testsh2gerr.o $(libs) -o testsh2gerr
testdecode: testdecode.o  libshape.a libutil.a
	$(F77) testdecode.o $(libs) -o testdecode
testortho: testortho.o  libshape.a libutil.a
	$(F77) testortho.o $(libs) -o testortho
testoper: testoper.o  libshape.a libutil.a
	$(F77) testoper.o $(libs) -o testoper
testconvolve: testconvolve.o  libshape.a libutil.a
	$(F77)  testconvolve.o $(libs) -o testconvolve
testround: testround.o  libshape.a libutil.a
	$(F77) testround.o $(libs) -o testround
testplot: testplot.o  libshape.a libutil.a
	$(F77)  testplot.o $(libs) -o testplot
testgaussians: testgaussians.o  libshape.a libutil.a
	$(F77) testgaussians.o $(libs) -o testgaussians
testmrq: testmrq.o libshape.a libutil.a
	$(F77) testmrq.o $(libs) -o testmrq
testpsfmat: testpsfmat.o  libshape.a libutil.a
	$(F77) testpsfmat.o $(libs) -o testpsfmat
testfitshear: testfitshear.o  libshape.a libutil.a
	$(F77) testfitshear.o $(libs) -o testfitshear
testfitshearmoffat: testfitshearmoffat.o  libshape.a libutil.a
	$(F77) testfitshearmoffat.o $(libs) -o testfitshearmoffat
thetest: thetest.o  libshape.a libutil.a
	$(F77) thetest.o $(libs) -o thetest
testreal: testreal.o  libshape.a libutil.a
	$(F77) testreal.o $(libs) -o testreal
testfits: testfits.o libshape.a libutil.a
	$(F77) testfits.o $(libs) -o testfits
testsherr: testsherr.o libshape.a libutil.a
	$(F77) testsherr.o $(libs) -o testsherr
testcentershape: testcentershape.o libshape.a libutil.a
	$(F77) testcentershape.o $(libs) -o testcentershape
testresiduals: testresiduals.o libshape.a libutil.a
	$(F77) testresiduals.o $(libs) -o testresiduals
testresiduals012: testresiduals012.o libshape.a libutil.a
	$(F77) testresiduals012.o $(libs) -o testresiduals012
testresiduals012_o2: testresiduals012_o2.o libshape.a libutil.a
	$(F77) testresiduals012_o2.o $(libs) -o testresiduals012_o2
testresiduals012conv: testresiduals012conv.o libshape.a libutil.a
	$(F77) testresiduals012conv.o $(libs) -o testresiduals012conv
testresiduals012mcbeta: testresiduals012mcbeta.o libshape.a libutil.a
	$(F77) testresiduals012mcbeta.o $(libs) -o testresiduals012mcbeta
testresiduals012mc: testresiduals012mc.o libshape.a libutil.a
	$(F77) testresiduals012mc.o $(libs) -o testresiduals012mc
testresiduals012_o2_mc: testresiduals012_o2_mc.o libshape.a libutil.a
	$(F77) testresiduals012_o2_mc.o $(libs) -o testresiduals012_o2_mc
testresiduals012_o3_mc: testresiduals012_o3_mc.o libshape.a libutil.a
	$(F77) testresiduals012_o3_mc.o $(libs) -o testresiduals012_o3_mc
testresiduals012mcnoise: testresiduals012mcnoise.o libshape.a libutil.a
	$(F77) testresiduals012mcnoise.o $(libs) -o testresiduals012mcnoise
test012: test012.o libshape.a libutil.a
	$(F77) test012.o  $(libs) -o test012
testfitmc: testfitmc.o libshape.a libutil.a
	$(F77) testfitmc.o $(libs) -o testfitmc 
testphotmc: testphotmc.o libshape.a libutil.a
	$(F77) testphotmc.o $(libs) -o testphotmc 
testsersic: testsersic.o libshape.a libutil.a
	$(F77) testsersic.o  $(libs) -o testsersic
testpsfcorphot: testpsfcorphot.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot.o psfcorphot.o $(libs) -o testpsfcorphot
testpsfcorphot_mc_noise: testpsfcorphot_mc_noise.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_noise.o psfcorphot.o $(libs) -o testpsfcorphot_mc_noise
testpsfcorphot_mc_new: testpsfcorphot_mc_new.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_new.o psfcorphot.o $(libs) -o testpsfcorphot_mc_new
testpsfcorphot_mc_new12: testpsfcorphot_mc_new12.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_new12.o psfcorphot.o $(libs) -o testpsfcorphot_mc_new12
testpsfcorphot_mc_new2: testpsfcorphot_mc_new2.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_new2.o psfcorphot.o $(libs) -o testpsfcorphot_mc_new2
testpsfcorphot_mc_new3: testpsfcorphot_mc_new3.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_new3.o psfcorphot.o $(libs) -o testpsfcorphot_mc_new3
testpsfcorphot_mc_noise_new: testpsfcorphot_mc_noise_new.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_noise_new.o psfcorphot.o $(libs) -o testpsfcorphot_mc_noise_new
testpsfcorphot_mc: testpsfcorphot_mc.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc.o psfcorphot.o $(libs) -o testpsfcorphot_mc
testpsfcorphot_mc_no13: testpsfcorphot_mc_no13.o psfcorphot.o libshape.a libutil.a
	$(F77) testpsfcorphot_mc_no13.o psfcorphot.o $(libs) -o testpsfcorphot_mc_no13
testcov: testcov.o libshape.a libutil.a
	$(F77) testcov.o $(libs) -o testcov

gethighordershearcoefs: gethighordershearcoefs.o round.o libshape.a libutil.a
	$(F77) gethighordershearcoefs.o round.o $(libs) -o gethighordershearcoefs
writehighordershearcoefs: writehighordershearcoefs.o round.o libshape.a libutil.a
	$(F77) writehighordershearcoefs.o round.o $(libs) -o writehighordershearcoefs
listpsfcorphot: listpsfcorphot.o psfcorphot.o libshape.a libutil.a
	$(F77) listpsfcorphot.o psfcorphot.o $(libs) -o listpsfcorphot
listpsfcorphot2: listpsfcorphot2.o psfcorphot.o libshape.a libutil.a
	$(F77) listpsfcorphot2.o psfcorphot.o $(libs) -o listpsfcorphot2
listpsfcorphot3: listpsfcorphot3.o psfcorphot.o libshape.a libutil.a
	$(F77) listpsfcorphot3.o psfcorphot.o $(libs) -o listpsfcorphot3
listpsfcorphot4: listpsfcorphot4.o psfcorphot.o libshape.a libutil.a
	$(F77) listpsfcorphot4.o psfcorphot.o $(libs) -o listpsfcorphot4
listpsfcorphot40: listpsfcorphot40.o psfcorphot.o libshape.a libutil.a
	$(F77) listpsfcorphot40.o psfcorphot.o $(libs) -o listpsfcorphot40
listfitbeta: listfitbeta.o libshape.a libutil.a
	$(F77) listfitbeta.o $(libs) -o listfitbeta
testhighorders: testhighorders.o libshape.a libutil.a
	$(F77) testhighorders.o $(libs) -o testhighorders

fakeimage: fakeimage.o libshape.a libutil.a
	$(F77) fakeimage.o $(libs) -o fakeimage
fakeimage2: fakeimage2.o libshape.a libutil.a
	$(F77) fakeimage2.o $(libs) -o fakeimage2
fakeimage3: fakeimage3.o libshape.a libutil.a
	$(F77) fakeimage3.o $(libs) -o fakeimage3
fakeimage4: fakeimage4.o libshape.a libutil.a
	$(F77) fakeimage4.o $(libs) -o fakeimage4
mcimage: mcimage.o libutil.a
	$(F77) mcimage.o $(libs) -o mcimage
mcimage_ell: mcimage_ell.o libshape.a libutil.a
	$(F77) mcimage_ell.o $(libs) -o mcimage_ell
mcimage_ell_step4: mcimage_ell_step4.o libshape.a libutil.a
	$(F77) mcimage_ell_step4.o $(libs) -o mcimage_ell_step4
mcimage_flx: mcimage_flx.o libshape.a libutil.a
	$(F77) mcimage_flx.o $(libs) -o mcimage_flx
mkimage_flx: mkimage_flx.o libshape.a libutil.a
	$(F77) mkimage_flx.o $(libs) -o mkimage_flx

fracpower: fracpower.o libshape.a libutil.a
	$(F77) fracpower.o $(libs) -o fracpower

selectpsfstars: selectpsfstars.o libshape.a libutil.a
	$(F77) selectpsfstars.o $(libs) -o selectpsfstars
selectpsfstars_auto: selectpsfstars_auto.o libshape.a libutil.a
	$(F77) selectpsfstars_auto.o $(libs) -o selectpsfstars_auto
psfcat2sh: psfcat2sh.o  libshape.a libutil.a
	$(F77) psfcat2sh.o $(libs) -o psfcat2sh
psfcat2sherr: psfcat2sherr.o  libshape.a libutil.a
	$(F77) psfcat2sherr.o $(libs) -o psfcat2sherr
psfcat2sherr_nobg: psfcat2sherr_nobg.o  libshape.a libutil.a
	$(F77) psfcat2sherr_nobg.o $(libs) -o psfcat2sherr_nobg
psfcat2sherr_nobg2: psfcat2sherr_nobg2.o  libshape.a libutil.a
	$(F77) psfcat2sherr_nobg2.o $(libs) -o psfcat2sherr_nobg2
psfcat2sherr_fixbeta: psfcat2sherr_fixbeta.o  libshape.a libutil.a
	$(F77) psfcat2sherr_fixbeta.o $(libs) -o psfcat2sherr_fixbeta
psfcat2sherr_fixbgbeta: psfcat2sherr_fixbgbeta.o  libshape.a libutil.a
	$(F77) psfcat2sherr_fixbgbeta.o $(libs) -o psfcat2sherr_fixbgbeta
psfcat2sherr_noctr: psfcat2sherr_noctr.o  libshape.a libutil.a
	$(F77) psfcat2sherr_noctr.o $(libs) -o psfcat2sherr_noctr
psfcat2sherr_2: psfcat2sherr_2.o  libshape.a libutil.a
	$(F77) psfcat2sherr_2.o $(libs) -o psfcat2sherr_2
fitpsfmap: fitpsfmap.o  libshape.a libutil.a
	$(F77) fitpsfmap.o $(libs) -o fitpsfmap
fitpsfmap0: fitpsfmap0.o  libshape.a libutil.a
	$(F77) fitpsfmap0.o $(libs) -o fitpsfmap0
fitpsfmap2: fitpsfmap2.o  libshape.a libutil.a
	$(F77) fitpsfmap2.o $(libs) -o fitpsfmap2
fitpsfmap3: fitpsfmap3.o  libshape.a libutil.a
	$(F77) fitpsfmap3.o $(libs) -o fitpsfmap3
fitkermap: fitkermap.o libshape.a libutil.a
	$(F77) fitkermap.o $(libs) -o fitkermap
fitkermaptwk: fitkermaptwk.o libshape.a libutil.a
	$(F77) fitkermaptwk.o $(libs) -o fitkermaptwk
keracfmap: keracfmap.o libshape.a libutil.a
	$(F77) keracfmap.o $(libs) -o keracfmap
kermapxgauss: kermapxgauss.o libshape.a libutil.a
	$(F77) kermapxgauss.o $(libs) -o kermapxgauss
kermapxkermap: kermapxkermap.o libshape.a libutil.a
	$(F77) kermapxkermap.o $(libs) -o kermapxkermap
sumkermaps: sumkermaps.o libshape.a libutil.a
	$(F77) sumkermaps.o $(libs) -o sumkermaps
showpsfmap: showpsfmap.o libshape.a libutil.a
	$(F77) showpsfmap.o $(libs) -o showpsfmap
showpsfmap2: showpsfmap2.o libshape.a libutil.a
	$(F77) showpsfmap2.o $(libs) -o showpsfmap2
showpsfmaptwk: showpsfmaptwk.o libshape.a libutil.a
	$(F77) showpsfmaptwk.o $(libs) -o showpsfmaptwk
showdxdy: showdxdy.o libshape.a libutil.a
	$(F77) showdxdy.o $(libs) -o showdxdy
psfcoefs: psfcoefs.o libshape.a 
	$(F77) psfcoefs.o $(libs) -o psfcoefs
psfmap2fits: psfmap2fits.o libshape.a libutil.a
	$(F77) psfmap2fits.o $(libs) -o psfmap2fits
pixpsfmap: pixpsfmap.o libshape.a libutil.a
	$(F77) pixpsfmap.o $(libs) -o pixpsfmap
pixpsfxshcpts: pixpsfxshcpts.o  libshape.a libutil.a
	$(F77) pixpsfxshcpts.o $(libs) -o pixpsfxshcpts
pixpsfxshcpts8: pixpsfxshcpts8.o  libshape.a libutil.a
	$(F77) pixpsfxshcpts8.o $(libs) -o pixpsfxshcpts8
pix2g: pix2g.o libshape.a libutil.a
	$(F77) pix2g.o $(libs) -o pix2g
pix2g8: pix2g8.o libshape.a libutil.a
	$(F77) pix2g8.o $(libs) -o pix2g8
pix2g-noagn: pix2g-noagn.o libshape.a libutil.a
	$(F77) pix2g-noagn.o $(libs) -o pix2g-noagn
carve: carve.o libutil.a
	$(F77) carve.o $(libs) -o carve
psfcat2gausker: psfcat2gausker.o libshape.a libutil.a
	$(F77)  psfcat2gausker.o $(libs) -o psfcat2gausker
psfcat2gauskerwithtweak: psfcat2gauskerwithtweak.o libshape.a libutil.a
	$(F77)  psfcat2gauskerwithtweak.o $(libs) -o psfcat2gauskerwithtweak

colplot: colplot.o
	$(F77) colplot.o $(libs) -o colplot

psf: selectpsfstars selectpsfstars_auto psfcat2sherr psfcat2sherr_2 psfcat2sherr_fixbeta psfcat2sherr_fixbgbeta fitpsfmap fitpsfmap0 fitpsfmap2 fitpsfmap3

cat2sh: cat2sh.o libshape.a libutil.a
	$(F77) cat2sh.o $(libs) -o cat2sh
cat2sherr: cat2sherr.o libshape.a libutil.a
	$(F77) cat2sherr.o $(libs) -o cat2sherr
cat2sherr_noctr: cat2sherr_noctr.o libshape.a libutil.a
	$(F77) cat2sherr_noctr.o $(libs) -o cat2sherr_noctr
cat2sherr_nobg: cat2sherr_nobg.o libshape.a libutil.a
	$(F77) cat2sherr_nobg.o $(libs) -o cat2sherr_nobg
cat2shcov: cat2shcov.o libshape.a libutil.a
	$(F77) cat2shcov.o $(libs) -o cat2shcov
sortbybeta: sortbybeta.o libutil.a
	$(F77) sortbybeta.o $(libs) -o sortbybeta
listbeta: listbeta.o libutil.a
	$(F77) listbeta.o $(libs) -o listbeta
sh2g: sh2g.o libshape.a libutil.a
	$(F77) sh2g.o $(libs) -o sh2g
sh2gerr: sh2gerr.o libshape.a libutil.a
	$(F77) sh2gerr.o $(libs) -o sh2gerr
sh2gerrshift: sh2gerrshift.o libshape.a libutil.a
	$(F77) sh2gerrshift.o $(libs) -o sh2gerrshift
sh2gerrshift012_o4: sh2gerrshift012_o4.o libshape.a libutil.a
	$(F77) sh2gerrshift012_o4.o $(libs) -o sh2gerrshift012_o4
sh2gerrshift_o4: sh2gerrshift_o4.o libshape.a libutil.a
	$(F77) sh2gerrshift_o4.o $(libs) -o sh2gerrshift_o4
sh2gerrshift012: sh2gerrshift012.o libshape.a libutil.a
	$(F77) sh2gerrshift012.o $(libs) -o sh2gerrshift012
sh2p: sh2p.o libshape.a libutil.a
	$(F77) sh2p.o $(libs) -o sh2p
imxshmap: imxshmap.o libshape.a libutil.a
	$(F77) imxshmap.o $(libs) -o imxshmap
imxshmap_rvdb: imxshmap_rvdb.o libshape.a libutil.a
	$(F77) imxshmap_rvdb.o $(libs) -o imxshmap_rvdb
imxshmap2: imxshmap2.o libshape.a libutil.a
	$(F77) imxshmap2.o $(libs) -o imxshmap2
imxshmap3: imxshmap3.o libshape.a libutil.a
	$(F77) imxshmap3.o $(libs) -o imxshmap3
imxshmapwithtweak: imxshmapwithtweak.o libshape.a libutil.a
	$(F77) imxshmapwithtweak.o $(libs) -o imxshmapwithtweak
trimwtmap: trimwtmap.o libshape.a libutil.a
	$(F77) trimwtmap.o $(libs) -o trimwtmap
growbadpix: growbadpix.o libshape.a libutil.a
	$(F77) growbadpix.o $(libs) -o growbadpix
imclip: imclip.o libshape.a libutil.a
	$(F77) imclip.o $(libs) -o imclip
imxshmaptweak: imxshmaptweak.o libshape.a libutil.a
	$(F77) imxshmaptweak.o $(libs) -o imxshmaptweak
isotropnoise: isotropnoise.o libshape.a libutil.a
	$(F77) isotropnoise.o $(libs) -o isotropnoise
checkgpsf: checkgpsf.o libshape.a libutil.a
	$(F77) checkgpsf.o $(libs) -o checkgpsf
checkgpsf2d: checkgpsf2d.o libshape.a libutil.a
	$(F77) checkgpsf2d.o $(libs) -o checkgpsf2d
checkkerneltwk: checkkerneltwk.o libshape.a libutil.a
	$(F77) checkkerneltwk.o $(libs) -o checkkerneltwk
reconstructimage: reconstructimage.o libshape.a libutil.a
	$(F77) reconstructimage.o $(libs) -o reconstructimage
shearstats: shearstats.o libshape.a libutil.a 
	$(F77) shearstats.o $(libs) -o shearstats
shearstats2: shearstats2.o libshape.a libutil.a 
	$(F77) shearstats2.o $(libs) -o shearstats2
shearstats2_ludo: shearstats2_ludo.o libshape.a libutil.a 
	$(F77) shearstats2_ludo.o $(libs) -o shearstats2_ludo
shearstatsmag: shearstatsmag.o libshape.a libutil.a 
	$(F77) shearstatsmag.o $(libs) -o shearstatsmag
shearstatsfw: shearstatsfw.o libshape.a libutil.a 
	$(F77) shearstatsfw.o $(libs) -o shearstatsfw
stackshear: stackshear.o libshape.a libutil.a 
	$(F77) stackshear.o $(libs) -o  stackshear
showsh: showsh.o libshape.a libutil.a
	$(F77) showsh.o $(libs) -o showsh
maskcat: maskcat.o libshape.a libutil.a
	$(F77) maskcat.o $(libs) -o maskcat
maskrawcat: maskrawcat.o libshape.a libutil.a
	$(F77) maskrawcat.o $(libs) -o maskrawcat
lensfitstars2sh: lensfitstars2sh.o  libshape.a libutil.a
	$(F77) lensfitstars2sh.o  $(libs) -o lensfitstars2sh
lensfitshmapresid: lensfitshmapresid.o  libshape.a libutil.a
	$(F77) lensfitshmapresid.o  $(libs) -o lensfitshmapresid
sh2ringtest: sh2ringtest.o  libshape.a libutil.a
	$(F77) sh2ringtest.o $(libs) -o sh2ringtest
stack4: stack4.o libutil.a
	$(F77) stack4.o $(libs) -o stack4
shmap2fwhm: shmap2fwhm.o 
	$(F77) shmap2fwhm.o -o shmap2fwhm

gapphot: gapphot.o libutil.a libshape.a
	$(F77) gapphot.o $(libs) -o gapphot
gapphot2: gapphot2.o libutil.a libshape.a
	$(F77) gapphot2.o $(libs) -o gapphot2
gapphot2a: gapphot2a.o libutil.a libshape.a
	$(F77) gapphot2a.o $(libs) -o gapphot2a
gapphot3: gapphot3.o libutil.a libshape.a
	$(F77) gapphot3.o $(libs) -o gapphot3
gapphot4: gapphot4.o libutil.a libshape.a
	$(F77) gapphot4.o $(libs) -o gapphot4
simpleapphot: simpleapphot.o libutil.a libshape.a
	$(F77) simpleapphot.o $(libs) -o simpleapphot

all: psf cat2sherr listpsfcorphot40 listfitbeta carve imxshmap maskcat maskrawcat pixpsfmap


clean:
	\rm *.o *.a *~ testmain testshear testdecode testmrq testortho testoper testconvolve testround testplot testgaussians testpsfmat testfitshear thetest testreal testfits test012 testsherr testresiduals testsersic selectpsfstars selectpsfstars_auto psfcat2sh fitpsfmap fitpsfmap0 showpsfmap cat2sh cat2sherr cat2shcov reconstructimage sh2g sh2gerr sh2gerrshift sh2gerrshift_o4 sh2gerrshift012 sh2p shearstats shearstats2 stackshear showsh fracpower fakeimage fakeimage2 fakeimage3 fakeimage4 testcentershape testfitmc testphotmc testresiduals012 testresiduals012conv testresiduals012mc testresiduals012mcnoise sortbybeta listpsfcorphot40 psfcat2sherr psfcat2sherr_fixbgbeta listfitbeta gethighordershearcoefs writehighordershearcoefs listbeta fitpsfmap fitpsfmap2 fitpsfmap3 carve imxshmap maskcat maskrawcat pixpsfmap gapphot gapphot2 gapphot2a gapphot3 gapphot4 psfcat2sherr_2 fitkermap checkgpsf keracfmap checkgpsf2d isotropnoise imxshmaptweak growbadpix kermapxgauss imxshmap2 imxshmap3 psfcat2sherr_fixbeta


backup:
	tar cfz backupdir/shpl-`date -Iseconds| sed s/://g`.tgz *.F *.csh *.smo Notes makefile* *readme* *.txt

towork:
	./tokuinder.csh
fromwork:
	./fromkuinder.csh

tests: testmain testshear testdecode testmrq testortho testoper testconvolve testround testplot testgaussians testpsfmat testfitshear thetest testreal testfits testsherr testcentershape fakeimage fakeimage2 fakeimage3 fakeimage4 test012 testresiduals testresiduals012 testresiduals012conv testresiduals012mc testresiduals012mcnoise testsersic testfitmc testphotmc
