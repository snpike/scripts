#!/bin/bash

DS9=/Users/sean/ds9/ds9

echo Have you run setsas and assigned environment variable SAS_ODF and SAS_CCF'?' Press enter to continue or type quit to exit.
read temp
if [ "$temp" == q ] || [ "$temp" == quit ]
then
	echo Exiting.
	exit 0
else
	echo Continuing.
fi

echo Which instrument'? (Enter pn, mos1, mos2, or all)'
read inst
while [ ! "$inst" == pn ] && [ ! "$inst" == mos1 ] && [ ! "$inst" == mos2 ] && [ ! "$inst" == all ] 
do
	echo Please enter pn, mos1, mos2, or all
	read inst
done

if [ "$inst" == all ]
then
	dets=(mos1 mos2 pn)
else
	dets=("$inst")
fi

regs=(src bkg)

echo Filtering data.

evselect table=mos1.fits withfilteredset=yes expression='(PATTERN<=12)&&(PI in [200:12000])&&#XMMEA_EM' filteredset=mos1_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >& mos1_filt.log

evselect table=mos2.fits withfilteredset=yes expression='(PATTERN<=12)&&(PI in [200:12000])&&#XMMEA_EM' filteredset=mos2_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >& mos2_filt.log

evselect table=pn.fits withfilteredset=yes expression='(PATTERN<=4)&&(PI in [200:15000])&&(FLAG==0)&&#XMMEA_EP' filteredset=pn_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >& pn_filt.log

echo Making images.

for x in "${dets[@]}"
do
	echo Processing "$x" data.
	filt_tag=filt
	evselect table="$x"_filt.fits withimageset=yes imageset="$x"_filt_im.fits xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600 >& "$x"_filt_im.log
	echo Pick the source and background regions in ds9
	$DS9 "$x"_filt_im.fits -scale log
	echo Enter the source region file:'['"$x"'_filt_src.reg]'
	read src_reg
	if [ -z $src_reg ]
	then
		src_reg="$x"_filt_src.reg
	fi
	while [ ! -f "$src_reg" ]
	do
		echo I cannot find that file, try again:
		read src_reg
	done
	echo Enter the background region file:'['"$x"'_filt_bkg.reg]'
	read bkg_reg
	if [ -z $bkg_reg ]
	then
		bkg_reg="$x"_filt_bkg.reg
	fi
	while [ ! -f "$bkg_reg" ]
	do
		echo I cannot find that file, try again:
		read bkg_reg
	done
	srcphys=$(awk '/./ {line=$0}; END {print line}' "$src_reg")
	bkgphys=$(awk '/./ {line=$0}; END {print line}' "$bkg_reg")

	echo "$srcphys"
	echo "$bkgphys"

	evselect table="$x"_filt.fits withfilteredset=yes expression='((X,Y) in '"$srcphys"')' filteredset="$x"_filt_src.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >& "$x"_filt_src.log
	evselect table="$x"_filt.fits withfilteredset=yes expression='((X,Y) in '"$bkgphys"')' filteredset="$x"_filt_bkg.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >& "$x"_filt_bkg.log
	
	if [ "$x" == pn ]
	then
		evselect table="$x".fits withrateset=yes rateset="$x"_flareltcrv.fits maketimecolumn=yes timebinsize=100 makeratecolumn=yes expression='#XMMEA_EP&&(PI>10000&&PI<12000)&&(PATTERN==0)' >& "$x"_flareltcrv.log
	else
		evselect table="$x".fits withrateset=yes rateset="$x"_flareltcrv.fits maketimecolumn=yes timebinsize=100 makeratecolumn=yes expression='#XMMEA_EM&&(PI>10000)&&(PATTERN==0)' >& "$x"_flareltcrv.log
	fi

	echo Please check for flaring background.
	
	fv "$x"_flareltcrv.fits
	echo Do we need to filter'?'
	read resp 
	if [ "$resp" == yes ] || [ "$resp" == y ]
	then
		echo Above what "$x" rate should we filter'?'
		read max_rate
		tabgtigen table="$x"_flareltcrv.fits expression='RATE<='"$max_rate" gtiset="$x"_gti.fits >& "$x"_gti.log
		filt_tag=clean
		evselect table="$x"_filt_src.fits withfilteredset=yes filteredset="$x"_"$filt_tag"_src.fits keepfilteroutput=yes expression='gti('"$x"'_gti.fits,TIME)' >& "$x"_"$filt_tag"_src.fits
		evselect table="$x"_filt_bkg.fits withfilteredset=yes filteredset="$x"_"$filt_tag"_bkg.fits keepfilteroutput=yes expression='gti('"$x"'_gti.fits,TIME)' >& "$x"_"$filt_tag"_bkg.fits
	else
		echo Continuing
	fi

	epatplot set="$x"_"$filt_tag"_src.fits plotfile="$x"_epat.ps useplotfile=yes withbackgroundset=yes backgroundset="$x"_"$filt_tag"_bkg.fits >& "$x"_epat.log

	echo Now check for pileup

	open "$x"_epat.ps

	echo Do we need to abort'?'
	read resp
	if [ "$resp" == yes ] || [ "$resp" == y ]
	then
		echo Run this script again and try removing a small region in the center of the source.
		echo Exiting.
		exit 1
	else
		echo Continuing
	fi

	echo Running backscale, rmfgen, and arfgen.
	for y in "${regs[@]}"
	do
		if [ "$x" == pn ]
		then
			evselect table="$x"_"$filt_tag"_"$y".fits energycolumn=PI withspectrumset=yes spectrumset="$x"_"$filt_tag"_"$y".pha spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 >& "$x"_"$filt_tag"_"$y"_pha.log
		else
			evselect table="$x"_"$filt_tag"_"$y".fits energycolumn=PI withspectrumset=yes spectrumset="$x"_"$filt_tag"_"$y".pha spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 >& "$x"_"$filt_tag"_"$y"_pha.log
		fi

		backscale spectrumset="$x"_"$filt_tag"_"$y".pha badpixlocation="$x"_"$filt_tag"_"$y".fits >& "$x"_"$y"_backscale.log
	done

	rmfgen rmfset="$x"_rmf.fits spectrumset="$x"_"$filt_tag"_src.pha >& "$x"_rmf.log
	arfgen arfset="$x"_arf.fits withrmfset=yes rmfset="$x"_rmf.fits spectrumset="$x"_"$filt_tag"_src.pha withbadpixcorr=yes badpixlocation="$x"_"$filt_tag"_src.fits >& "$x"_arf.log

done

echo Done. Run grppha to make grouped spectral files.
