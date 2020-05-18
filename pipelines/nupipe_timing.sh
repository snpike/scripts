#!/bin/bash

echo Enter the OBSID
read OBSID
while [ -z $OBSID ] 
do
    echo Make sure you are in the right directory and enter the OBSID
    read OBSID
done

nupipeline indir=./"$OBSID" steminputs=nu"$OBSID" saamode=NONE tentacle=no outdir=./"$OBSID"_pipe_out/ clobber=yes cleancols=no >& pipe.log

if [ ! -d "$PWD"/"$OBSID"_pipe_out ]
then
    echo nupipeline failed
    exit 1
fi

echo nupipeline done

cd "$OBSID"_pipe_out
echo Pick the FPMA source and background regions in ds9 '(save in decimal form)'
$DS9 nu"$OBSID"A01_cl.evt -scale log

echo Pick the FPMB source and background regions in ds9 '(save in decimal form)'
$DS9 nu"$OBSID"B01_cl.evt -scale log


cd ..