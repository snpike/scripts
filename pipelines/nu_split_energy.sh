#!/bin/bash

echo Enter the OBSID
read OBSID
while [ -z $OBSID ] 
do
    echo Make sure you are in the right directory and enter the OBSID
    read OBSID
done

echo Enter the low end of the energy range '(integer, in keV)'
read e_low
while [ $(( e_low -lt 3 )) ] || [ $(( e_low -gt 78 )) ]
do
    echo The energy range must lie between 3 and 78 keV
    read e_low
done

echo Enter the high end of the energy range '(integer, in keV)'
read e_hi
while [ $(( e_hi -lt 3 )) ] || [ $(( e_hi -gt 78 )) ]
do
    echo The energy range must lie between 3 and 78 keV
    read e_hi
done

echo Enter the bin size '(integer seconds)'
read bin_size

pi_low=$(( ((e_low*1000) - 1600)/40 ))
pi_low=${pi_low%.*}
echo $pi_low

pi_hi=$(( ((e_hi*1000) - 1600)/40 ))
pi_hi=${pi_hi%.*}
echo $pi_hi

srcA_str=$(awk '/circle/ {line=$1}; END {print line}' ./"$OBSID"_pipe_out/srcA.reg)

srcA_temp=$(echo "$srcA_str" | cut -d ',' -f 1)
srcA_ra=${srcA_temp:7}
srcA_dec=$(echo "$srcA_str" | cut -d ',' -f 2)

srcB_str=$(awk '/circle/ {line=$1}; END {print line}' ./"$OBSID"_pipe_out/srcB.reg)

srcB_temp=$(echo "$srcB_str" | cut -d ',' -f 1)
srcB_ra=${srcB_temp:7}
srcB_dec=$(echo "$srcB_str" | cut -d ',' -f 2)

echo Enter the path to the clockfile: 
read clock
while [ -z $clock ]
do
    echo I can\'t find that file, try again: 
    read clock
done

nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_"$e_low"to"$e_hi"keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize="$bin_size" pilow="$pi_low" barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" pihigh="$pi_hi" clobber=yes >& nuproductsA_"$e_low"to"$e_hi"keV.log &
nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_"$e_low"to"$e_hi"keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize="$bin_size" pilow="$pi_low" barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcB_ra" srcdec_barycorr="$srcB_dec" pihigh="$pi_hi" clobber=yes >& nuproductsB_"$e_low"to"$e_hi"keV.log &