#!/bin/bash

RESULTS="./tbabs_simpl_diskbb_gauss_chisq.csv"

if [ -e $RESULTS ]
then
	rm $RESULTS
fi
touch $RESULTS

numsims=1000
numsteps=1000

dirnum=1
while [ $dirnum -le $numsims ]
do
	dir="/Volumes/Samsung_1TB/AstroData/Circinus_ULX5/nustar"

   	cd $dir
	echo "query yes" > sim.xco
	echo "statistic chi" >> sim.xco
	echo "cpd /null" >> sim.xco
	echo "data none" >> sim.xco
	echo "model none" >> sim.xco
	echo "@tbabs_diskpbb.xcm" >> sim.xco
	echo "fakeit" >> sim.xco
	echo "y" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "\n" >> sim.xco
	echo "ignore **:**-0.3,10.-**" >> sim.xco
	echo "fit" >> sim.xco
        echo "log tbabs_diskpbb_fak.log" >> sim.xco
        echo "show all" >> sim.xco
        echo "log none" >> sim.xco
	echo "addc 4 gauss" >> sim.xco
	echo "2.4 0.02 0.3 0.3 10.0 10.0" >> sim.xco
	echo "0.0001,-1" >> sim.xco
	echo "3.2e-6 8e-7 0 0 1 1" >> sim.xco
	echo " " >> sim.xco
	echo " " >> sim.xco
	echo " " >> sim.xco
	echo " " >> sim.xco
	echo " " >> sim.xco
	echo " " >> sim.xco
    echo "fit" >> sim.xco
    echo "log tbabs_diskpbb_gauss_fak.log" >> sim.xco
    echo "show all" >> sim.xco
    echo "log none" >> sim.xco
   	echo "exit" >> sim.xco	

	rm *fak*
	xspec - sim.xco	

        chisq0=`awk '/Fit statistic : Chi-Squared/ {chisq=$6}; END {print chisq}' tbabs_diskpbb_fak.log`
        chisq1=`awk '/Fit statistic : Chi-Squared/ {chisq=$6}; END {print chisq}' tbabs_diskpbb_gauss_fak.log`

        echo $dirnum , $chisq0 , $chisq1 >> $RESULTS


	((dirnum++))
done
