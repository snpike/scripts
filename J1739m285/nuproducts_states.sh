#!/bin/bash

states=(soft_soft soft soft_mid mid hard_mid hard hard_hard)

cd /disk/lif2/spike/J1739m285

for x in "${states[@]}"
do
	nuproducts indir=./90501343002_pipe_out instrument=FPMA steminputs=nu90501343002 outdir=./products_state_"$x" srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg usrgtifile=./"$x".gti >& nuproductsA_state_"$x".log
	nuproducts indir=./90501343002_pipe_out instrument=FPMB steminputs=nu90501343002 outdir=./products_state_"$x" srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg usrgtifile=./"$x".gti >& nuproductsB_state_"$x".log
	cd products_state_"$x"
	ftgrouppha infile=nu90501343002A01_sr.pha backfile=nu90501343002A01_bk.pha outfile=ftgroupedA.pha grouptype=opt respfile=nu90501343002A01_sr.rmf clobber=true >& ftgrouppha_A.log
	ftgrouppha infile=nu90501343002B01_sr.pha backfile=nu90501343002B01_bk.pha outfile=ftgroupedB.pha grouptype=opt respfile=nu90501343002B01_sr.rmf clobber=true >& ftgrouppha_B.log
	cd /disk/lif2/spike/J1739m285
done