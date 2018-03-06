#!/bin/tcsh
#
#run sextractor on fits image; extract star catalogue based on flux/fwhm plot
#

# point swpath variable at installation directory
set swpath = $GAaP
setenv PATH {$swpath}:$PATH

set im = $1
set psfcat = $2

if ! -e default.param then
 \cp -f $swpath/sexconfig/default.* .
endif

sex $im -CATALOG_NAME tmp1.cat -VERBOSE_TYPE QUIET
# find 20th brightest source
set n = `grep -v \# tmp1.cat |wc -l`
#echo Catalogue contains $n lines
set f20 = `grep -v \# tmp1.cat|sort -gr -k 3|head -20|tail -1|awk '{print($3)}'`
#echo Top 20th flux is $f20
# find median radius of sources within factor 30 of 20th flux
set fw1 = `awk -v f20=$f20 '$3<f20 && $3>f20/30 && $5>0.5 {print($5)}' tmp1.cat|median`
#echo Median radius of top decade in flux is $fw1
set fw2 = `awk -v f20=$f20 -v fw1=$fw1 '$3<f20 && $3>f20/100 && $5>0.5 && $5<1.5*fw1 {print($5)}' tmp1.cat|median`
#echo Next iteration: median radius of top decade in flux is $fw2
awk -v fw=$fw2 '$5>fw/2 && $5<fw+1 && $3>0' tmp1.cat | sort -gr -k 3 > tmp2.cat
# find top two decades in flux; remove brightest 10 pct
set f1 = `head -1 tmp2.cat|awk '{print($3)}'`
awk -v f1=$f1 '$3>f1/100' tmp2.cat > tmp3.cat
set n90 = `cat tmp3.cat|wc -l|awk '{print(int($1*0.9))}'`
tail -$n90 tmp3.cat > $psfcat
\cp tmp1.cat sextractor.cat
\rm tmp1.cat tmp2.cat tmp3.cat

