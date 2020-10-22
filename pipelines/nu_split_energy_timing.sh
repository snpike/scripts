#!/bin/bash

echo Enter the OBSID
read OBSID
while [ -z $OBSID ] 
do
    echo Make sure you are in the right directory and enter the OBSID
    read OBSID
done

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


# 3-78 keV
nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_3to78keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=35 pihigh=1909 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb inarffile="$CALDB"/data/nustar/fpm/bcf/arf/nuA20100101v006.arf clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& log_files/nuproductsA_3to78keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_3to78keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=35 pihigh=1909 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsA_3to78keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_3to78keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize=0.5 pilow=35 pihigh=1909 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsB_3to78keV.log &

# # 3-6 keV
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_3to6keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=35 pihigh=110 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsA_3to6keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_3to6keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize=0.5 pilow=35 pihigh=110 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsB_3to6keV.log &

# # 6-12 keV
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_6to12keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=111 pihigh=260 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsA_6to12keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_6to12keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize=0.5 pilow=111 pihigh=260 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsB_6to12keV.log &

# # 12-25 keV
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_12to25keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=261 pihigh=585 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsA_12to25keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_12to25keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize=0.5 pilow=261 pihigh=585 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsB_12to25keV.log &

# # 25-50 keV
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMA steminputs=nu"$OBSID" outdir=./products_25to50keV srcregionfile=./"$OBSID"_pipe_out/srcA.reg bkgregionfile=./"$OBSID"_pipe_out/bkgA.reg binsize=0.5 pilow=586 pihigh=1210 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"A.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsA_25to50keV.log &
# nuproducts indir=./"$OBSID"_pipe_out clobber=yes instrument=FPMB steminputs=nu"$OBSID" outdir=./products_25to50keV srcregionfile=./"$OBSID"_pipe_out/srcB.reg bkgregionfile=./"$OBSID"_pipe_out/bkgB.reg binsize=0.5 pilow=586 pihigh=1210 barycorr=yes orbitfile=./"$OBSID"_pipe_out/nu"$OBSID"B.attorb clockfile="$clock" srcra_barycorr="$srcA_ra" srcdec_barycorr="$srcA_dec" clobber=yes >& nuproductsB_25to50keV.log &
