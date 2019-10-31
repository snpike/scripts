#!/bin/bash

epochs=(1 2 3 4 5 6 7 8 9 10)

cd /disk/lif2/spike/J1739m285

for x in "${epochs[@]}"
do
	nuproducts indir=./90501343002_pipe_out instrument=FPMA steminputs=nu90501343002 outdir=./products_epoch"$x" srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=0.5 usrgtifile=./90501343002_pipe_out/epoch_"$x".curs_gti >& nuproductsA_epoch"$x".log
	nuproducts indir=./90501343002_pipe_out instrument=FPMB steminputs=nu90501343002 outdir=./products_epoch"$x" srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=0.5 usrgtifile=./90501343002_pipe_out/epoch_"$x".curs_gti >& nuproductsB_epoch"$x".log
	cd products_epoch"$x"
	grppha nu90501343002A01_sr.pha groupedA_25.pha comm="bad 0-35 & bad 1910-4095 & group min 25 & exit" >& grpphaA.log
	grppha nu90501343002B01_sr.pha groupedB_25.pha comm="bad 0-35 & bad 1910-4095 & group min 25 & exit" >& grpphaB.log
	cd /disk/lif2/spike/J1739m285
done