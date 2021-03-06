#!/bin/bash

cd /disk/lif2/spike/J1739m285

# # 3.0 to 10.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_3to10keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=35 pihigh=210 >& nuproductsA_3to10keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_3to10keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=35 pihigh=210 >& nuproductsB_3to10keV.log

# # 10.0 to 50.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_10to50keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=211 pihigh=1210 >& nuproductsA_10to50keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_10to50keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=211 pihigh=1210 >& nuproductsB_10to50keV.log

# # 10.0 to 25.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_10to25keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=211 pihigh=585 >& nuproductsA_10to25keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_10to25keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=211 pihigh=585 >& nuproductsB_10to25keV.log

# # 3.0 to 12.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_3to12keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=35 pihigh=260 >& nuproductsA_3to12keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_3to12keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=35 pihigh=260 >& nuproductsB_3to12keV.log

# # 12.0 to 50.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_12to50keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=261 pihigh=1210 >& nuproductsA_12to50keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_12to50keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=261 pihigh=1210 >& nuproductsB_12to50keV.log

# # 12.0 to 25.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_12to25keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=261 pihigh=585 >& nuproductsA_12to25keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_12to25keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=261 pihigh=585 >& nuproductsB_12to25keV.log

# # 25.0 to 50.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_25to50keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=586 pihigh=1210 >& nuproductsA_25to50keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_25to50keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=586 pihigh=1210 >& nuproductsB_25to50keV.log

# # 3.0 to 6.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_3to6keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=35 pihigh=110 >& nuproductsA_3to6keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_3to6keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=35 pihigh=110 >& nuproductsB_3to6keV.log

# # 6.0 to 12.0
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_6to12keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=111 pihigh=260 >& nuproductsA_6to12keV.log
# nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_6to12keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=111 pihigh=260 >& nuproductsB_6to12keV.log

# 3.0 to 78.0
nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMA steminputs=nu90501343002 outdir=./products_3to78keV srcregionfile=./90501343002_pipe_out/srcA.reg bkgregionfile=./90501343002_pipe_out/bkgA.reg binsize=128 pilow=35 pihigh=1909 >& nuproductsA_3to78keV.log
nuproducts indir=./90501343002_pipe_out clobber=yes instrument=FPMB steminputs=nu90501343002 outdir=./products_3to78keV srcregionfile=./90501343002_pipe_out/srcB.reg bkgregionfile=./90501343002_pipe_out/bkgB.reg binsize=128 pilow=35 pihigh=1909 >& nuproductsB_3to78keV.log
