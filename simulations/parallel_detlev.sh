#!/bin/bash

n_jobs=$1
n_sims=$2
job=1

while [ $job -le $n_jobs ]
do 
    xterm -hold -e "zsh ./line_detlev.sh $job $n_sims" &

    ((job++))

done
