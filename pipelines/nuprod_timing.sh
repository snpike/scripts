#!/bin/bash

echo Enter the OBSID
read OBSID
while [ -z $OBSID ] 
do
    echo Make sure you are in the right directory and enter the OBSID
    read OBSID
done


cd "$OBSID"_pipe_out
echo Enter the source region file:'[srcA.reg]'
read srcA
if [ -z $srcA ]
then
    srcA=srcA.reg
fi
while [ ! -f "$srcA" ]
do
    echo I cannot find that file, try again:
    read srcA
done
echo Enter the background region file:'[bkgA.reg]'
read bkgA
if [ -z $bkgA ]
then
    bkgA=bkgA.reg
fi
while [ ! -f "$bkgA" ]
do
    echo I cannot find that file, try again:
    read bkgA
done

srcA_str=$(awk '/circle/ {line=$1}; END {print line}' "$srcA")

srcA_temp=$(echo "$srcA_str" | cut -d ',' -f 1)
srcA_ra=${srcA_temp:7}
srcA_dec=$(echo "$srcA_str" | cut -d ',' -f 2)

echo Enter the source region file:'[srcB.reg]'
read srcB
if [ -z $srcB ]
then
    srcB=srcB.reg
fi
while [ ! -f "$srcB" ]
do
    echo I cannot find that file, try again:
    read srcB
done
echo Enter the background region file:'[bkgB.reg]'
read bkgB
if [ -z $bkgB ]
then
    bkgB=bkgB.reg
fi
while [ ! -f "$bkgB" ]
do
    echo I cannot find that file, try again:
    read bkgB
done

srcB_str=$(awk '/circle/ {line=$1}; END {print line}' "$srcB")

srcB_temp=$(echo "$srcB_str" | cut -d ',' -f 1)
srcB_ra=${srcB_temp:7}
srcB_dec=$(echo "$srcB_str" | cut -d ',' -f 2)

cd ..

echo Enter the path to the clockfile: 
read clock
while [ -z $clock ]
do
    echo I can\'t find that file, try again: 
    read clock
done


nuproducts indir=./"$OBSID"_pipe_out/ instrument=FPMA steminputs=nu"$OBSID" outdir=./"$OBSID"_products srcregionfile=./"$OBSID"_pipe_out/"$srcA" bkgregionfile=./"$OBSID"_pipe_out/"$bkgA" binsize=1 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" write_baryevtfile=yes clobber=yes >& "$OBSID"_nuproducts_A.log &
nuproducts indir=./"$OBSID"_pipe_out/ instrument=FPMB steminputs=nu"$OBSID" outdir=./"$OBSID"_products srcregionfile=./"$OBSID"_pipe_out/"$srcB" bkgregionfile=./"$OBSID"_pipe_out/"$bkgB" binsize=1 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" write_baryevtfile=yes clobber=yes >& "$OBSID"_nuproducts_B.log &

# nuproducts indir=./"$OBSID"_pipe_out/ instrument=FPMA steminputs=nu"$OBSID" outdir=./"$OBSID"_timing_products srcregionfile=./"$OBSID"_pipe_out/"$srcA" bkgregionfile=./"$OBSID"_pipe_out/"$bkgA" binsize=0.0009765625 barycorr=yes write_baryevtfile=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& "$OBSID"_nuproducts_A_timing.log &
# nuproducts indir=./"$OBSID"_pipe_out/ instrument=FPMB steminputs=nu"$OBSID" outdir=./"$OBSID"_timing_products srcregionfile=./"$OBSID"_pipe_out/"$srcB" bkgregionfile=./"$OBSID"_pipe_out/"$bkgB" binsize=0.0009765625 barycorr=yes write_baryevtfile=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& "$OBSID"_nuproducts_B_timing.log

# barycorr infile=./"$OBSID"_pipe_out/nu"$OBSID"A01_cl.evt outfile=./timing_products/nu"$OBSID"A01_cl_bc.evt orbitfiles=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" ra="$srcA_ra" dec="$srcA_dec" clobber=yes >& barycorr_A.log &
# barycorr infile=./"$OBSID"_pipe_out/nu"$OBSID"B01_cl.evt outfile=./timing_products/nu"$OBSID"B01_cl_bc.evt orbitfiles=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" ra="$srcB_ra" dec="$srcB_dec" clobber=yes >& barycorr_B.log &

# barycorr infile=./"$OBSID"/hk/nu"$OBSID"A_fpm.hk outfile=./timing_products/nu"$OBSID"A_fpm_bc.hk orbitfiles=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" ra="$srcA_ra" dec="$srcA_dec" clobber=yes >& barycorr_A_hk.log &
# barycorr infile=./"$OBSID"/hk/nu"$OBSID"B_fpm.hk outfile=./timing_products/nu"$OBSID"B_fpm_bc.hk orbitfiles=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" ra="$srcA_ra" dec="$srcA_dec" clobber=yes >& barycorr_B_hk.log &

