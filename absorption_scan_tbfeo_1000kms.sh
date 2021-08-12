#!/bin/bash

# ^ This tells your terminal this is a bash script

# Specify where to keep the results
RESULTS="./tbfeo_simpl_diskbb_gabs_1000kms_cstat.csv"

# If results already exist, clobber them
if [ -e $RESULTS ]
then
	rm $RESULTS
fi
touch $RESULTS

# Start energy at 5 eV, end at 10 eV, take steps of 5 eV
# for cyclotron, you should change these values
centroid_eV=5000
linestep=5
linemax=10000

# While loop until we reach the max energy
while [ $centroid_eV -le $linemax ]
do

	# Translate the energy in eV to keV -- arithmetic in bash is confusing so I just manipulate the strings here
    centroid_keV=${centroid_eV: : -3}.${centroid_eV: -3: 3}

    # All of the echo commands are essentially writing an xspec script called "sim_feo.xco"
	echo "query yes" > sim_feo.xco
	echo "data none" >> sim_feo.xco
	echo "model none" >> sim_feo.xco

	# Load in your baseline model -- @bb_nthcomp.xcm for example
	echo "@tbfeo_simpl_diskbb.xcm" >> sim_feo.xco
	# You can fit again for good measure or comment this out, probably doesn't matter
	echo "fit" >> sim_feo.xco

	# Tell xspec to log the baseline info in the file "tbfeo_simpl_diskbb.log" (change the file name)
    echo "log tbfeo_simpl_diskbb.log" >> sim_feo.xco
    echo "show all" >> sim_feo.xco
    echo "log none" >> sim_feo.xco

    # Here I freeze the continuum parameters then add the gaussian absorption
    # For cyclotron, "add 2 cyclabs" instead.
    echo "freeze 1-80" >> sim_feo.xco
	echo "addc 2 gabs" >> sim_feo.xco
    
    # Initialize the parameters 
    # for fpma
	echo "$centroid_keV,-1" >> sim_feo.xco
	echo "=2/300." >> sim_feo.xco
	echo "1e-2" >> sim_feo.xco

	# and for fpmb (leaving it blank will automatically tie them to fpma values)
	echo " " >> sim_feo.xco
	echo " " >> sim_feo.xco
	echo " " >> sim_feo.xco

	# Also allow the continuum to move up and down a little to compensate by adding a multiplicative factor
    echo "addc 1 const" >> sim_feo.xco
    # fpma
	echo "1.0" >> sim_feo.xco
	# fpmb
	echo " " >> sim_feo.xco

	# Fit the new model and log the results in a different file, in this case "tbfeo_simpl_diskbb_gabs.log"
    echo "fit" >> sim_feo.xco
    echo "log tbfeo_simpl_diskbb_gabs.log" >> sim_feo.xco
    echo "show all" >> sim_feo.xco
    echo "log none" >> sim_feo.xco

    # Exit the xspec session
   	echo "exit" >> sim_feo.xco	

   	### At this point, we have not actually run xspec, we have only written an xspec script!!

   	# Delete results from the previous iteration of the loop 
    rm tbfeo_simpl_diskbb.log
    rm tbfeo_simpl_diskbb_gabs.log

    # Now pass the script to xspec
	xspec - sim_feo.xco	

		# Read the resulting c-statistics (0=baseline, 1=new)
        cstat0=`awk '/Total fit statistic/ {cstat=$4}; END {print cstat}' tbfeo_simpl_diskbb.log`
        cstat1=`awk '/Total fit statistic/ {cstat=$4}; END {print cstat}' tbfeo_simpl_diskbb_gabs.log`

        # Write the results into the $RESULTS file as well as the centroid energy
        echo $centroid_keV , $cstat0 , $cstat1 >> $RESULTS

    # Increment centroid energy
	((centroid_eV+=$linestep))

done

# Clean up remaining log files and scripts
rm tbfeo_simpl_diskbb.log
rm tbfeo_simpl_diskbb_gabs.log
rm sim_feo.xco