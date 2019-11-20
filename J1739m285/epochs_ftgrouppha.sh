#!/bin/bash

epochs=(1 2 3 4 5 6 7 8 9 10)

cd /Volumes/LaCie/AstroData/J1739-285/nustar

for x in "${epochs[@]}"
do
	cd products_epoch"$x"
	ftgrouppha infile=nu90501343002A01_sr.pha outfile=ftgroupedA.pha backfile=nu90501343002A01_bk.pha grouptype=opt respfile=nu90501343002A01_sr.rmf clobber=true >& ftgroupphaA.log
	ftgrouppha infile=nu90501343002B01_sr.pha outfile=ftgroupedB.pha backfile=nu90501343002B01_bk.pha grouptype=opt respfile=nu90501343002B01_sr.rmf clobber=true >& ftgroupphaB.log
	cd /Volumes/LaCie/AstroData/J1739-285/nustar
done