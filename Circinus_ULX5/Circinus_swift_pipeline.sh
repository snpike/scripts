#!/bin/bash

OBSID=$1
# outfile="$1_pipe_out"
# logfile=

xrtpipeline indir=./$OBSID outdir=./$OBSID_pipe_out steminputs=sw$OBSID srcra=OBJECT srcdec=OBJECT >& $OBSID_pipe.log
cd $OBSID_pipe_out
ds9