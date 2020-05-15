#!/bin/bash

DS9=/Users/sean/ds9/ds9

if [ -z $1 ]
then
	echo An OBSID is required
	exit 1
fi

OBSID="$1"
echo $OBSID

xrtpipeline indir="$PWD"/"$OBSID" outdir="$PWD"/"$OBSID"_pipe_out steminputs=sw"$OBSID" srcra=OBJECT srcdec=OBJECT clobber=yes >& "$PWD"/"$OBSID"_pipe.log

if [ ! -d "$PWD"/"$OBSID"_pipe_out ]
then
	echo xrtpipeline failed
	exit 1
fi

echo xrtpipeline done

cd "$OBSID"_pipe_out
echo Pick the source and background regions in ds9
$DS9 sw"$OBSID"xpcw3po_cl.evt -scale log
echo Enter the source region file:'[src.reg]'
read src_reg
if [ -z $src_reg ]
then
	src_reg=src.reg
fi
while [ ! -f "$src_reg" ]
do
	echo I cannot find that file, try again:
	read src_reg
done
echo Enter the background region file:'[bkg.reg]'
read bkg_reg
if [ -z $bkg_reg ]
then
	bkg_reg=bkg.reg
fi
while [ ! -f "$bkg_reg" ]
do
	echo I cannot find that file, try again:
	read bkg_reg
done

$DS9 sw"$OBSID"xpcw3po_cl.evt -scale log -regions "$src_reg" -regions select all -regions centroid -regions system physical -regions skyformat degrees -regions save src_phys.reg -exit

echo extraction regions generated

if [ -f xsel.xco ]
then
	rm xsel.xco
fi

touch xsel.xco

echo "$OBSID" > xsel.xco
echo read events sw"$OBSID"xpcw3po_cl.evt >> xsel.xco
echo ./ >> xsel.xco
echo yes >> xsel.xco
echo filter region "$src_reg" >> xsel.xco
echo extract spectrum >> xsel.xco
echo save spectrum src_pc.pha >> xsel.xco
echo clear region >> xsel.xco
echo filter region "$bkg_reg" >> xsel.xco
echo extract spectrum >> xsel.xco
echo save spectrum bkg_pc.pha >> xsel.xco
echo exit >> xsel.xco
echo no >> xsel.xco

echo xselect script generated

xselect < xsel.xco

echo spectra extracted via xselect

srcphys=$(awk '/circle/ {line=$1}; END {print line}' src_phys.reg)

srcx_temp=$(echo "$srcphys" | cut -d ',' -f 1)
srcx=${srcx_temp:7}
srcy=$(echo "$srcphys" | cut -d ',' -f 2)

echo srcx and srcy extracted
echo srcx: "$srcx"
echo srcy: "$srcy"

xrtmkarf phafile=./src_pc.pha srcx="$srcx" srcy="$srcy" outfile=./src_pc.arf psfflag=yes expofile=sw"$OBSID"xpcw3po_ex.img >& xrtmkarf.log

echo arf file generated

grppha src_pc.pha src_pc_"$OBSID"_10.grp comm="group min 10 & chkey BACKFILE bkg_pc.pha & chkey ANCRFILE src_pc.arf & chkey RESPFILE $CALDB/data/swift/xrt/cpf/rmf/swxpc0to12s6_20130101v014.rmf & exit" >& grppha10.log

echo grouped spectrum generated
echo Processing has finished. Double check log files to make sure that each step was successful.


