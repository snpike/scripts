#!/bin/bash

jobnum=$1
n_sims=$2

RESULTS="./tbfeo_diskbb_simpl_diskbb_gabs_1000kms_detlev_cstat_"$jobnum".csv"

echo $RESULTS

if [ -e $RESULTS ]
then
	rm $RESULTS
fi
touch $RESULTS

dirnum=1
while [ $dirnum -le $n_sims ]
do 
    echo "query yes" > sim_detlev_"$jobnum".xco
    echo "data none" >> sim_detlev_"$jobnum".xco
    echo "model none" >> sim_detlev_"$jobnum".xco
    echo "@tbfeo_diskbb_simpl_diskbb_all.xcm" >> sim_detlev_"$jobnum".xco
    
    echo "fakeit" >> sim_detlev_"$jobnum".xco

    echo "y" >> sim_detlev_"$jobnum".xco
    echo "fake$jobnum" >> sim_detlev_"$jobnum".xco

    echo " " >> sim_detlev_"$jobnum".xco
    echo "39730.9, 1.00000, 39730.9" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "39823.9, 1.00000, 39823.9" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "34409.5, 1.00000, 34407.2" >> sim_detlev_"$jobnum".xco

    echo " " >> sim_detlev_"$jobnum".xco
    echo "19056.7, 1.00000, 19056.7" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "18975.2, 1.00000, 18975.2" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "19974.2, 1.00000, 19964.9" >> sim_detlev_"$jobnum".xco

    echo " " >> sim_detlev_"$jobnum".xco
    echo "116690., 1.00000, 116690." >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "115493., 1.00000, 115493." >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "78342.9, 1.00000, 78342.7" >> sim_detlev_"$jobnum".xco

    echo " " >> sim_detlev_"$jobnum".xco
    echo "99160.5, 1.00000, 99160.5" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "63966.5, 1.00000, 63967.2" >> sim_detlev_"$jobnum".xco

    echo " " >> sim_detlev_"$jobnum".xco
    echo "47117.5, 1.00000, 47117.5" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "47339.7, 1.00000, 47339.7" >> sim_detlev_"$jobnum".xco
    echo " " >> sim_detlev_"$jobnum".xco
    echo "24242.3, 1.00000, 24241.9" >> sim_detlev_"$jobnum".xco
    
    echo "ignore **:**-0.3,10.-**" >> sim_detlev_"$jobnum".xco
    echo "fit" >> sim_detlev_"$jobnum".xco
    echo "save all tbfeo_diskbb_simpl_diskbb_fak_$jobnum.xcm" >> sim_detlev_"$jobnum".xco
    echo "exit" >> sim_detlev_"$jobnum".xco

    xspec - sim_detlev_"$jobnum".xco

    centroid_eV=300
    linestep=5
    linemax=10000

    while [ $centroid_eV -le $linemax ]
    do
    # dir="/Volumes/Samsung_1TB/AstroData/Circinus_ULX5/nustar/nustar_30002038006/30002038006_products"

        # cd $dir

        centroid_keV=${centroid_eV: : -3}.${centroid_eV: -3: 3}
        echo "query yes" > sim_detlev_"$jobnum".xco
        echo "data none" >> sim_detlev_"$jobnum".xco
        echo "model none" >> sim_detlev_"$jobnum".xco
        echo "@tbfeo_diskbb_simpl_diskbb_fak_$jobnum.xcm" >> sim_detlev_"$jobnum".xco
        # echo "fit" >> sim_detlev_"$jobnum".xco
            echo "log tbfeo_diskbb_simpl_diskbb_$jobnum.log" >> sim_detlev_"$jobnum".xco
            echo "show all" >> sim_detlev_"$jobnum".xco
            echo "log none" >> sim_detlev_"$jobnum".xco

        echo "freeze 1-168" >> sim_detlev_"$jobnum".xco
        echo "addc 2 gabs" >> sim_detlev_"$jobnum".xco
        
        echo "$centroid_keV,-1" >> sim_detlev_"$jobnum".xco
        echo "=2/300." >> sim_detlev_"$jobnum".xco
        echo "1e-2" >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo "addc 1 const" >> sim_detlev_"$jobnum".xco
        echo "1.0" >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco  

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco

        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco
        echo " " >> sim_detlev_"$jobnum".xco    

        echo "fit" >> sim_detlev_"$jobnum".xco
        echo "log tbfeo_diskbb_simpl_diskbb_gabs_$jobnum.log" >> sim_detlev_"$jobnum".xco
        echo "show all" >> sim_detlev_"$jobnum".xco
        echo "log none" >> sim_detlev_"$jobnum".xco
        echo "exit" >> sim_detlev_"$jobnum".xco  

        xspec - sim_detlev_"$jobnum".xco 

        cstat0=`awk '/Total fit statistic/ {cstat=$4}; END {print cstat}' tbfeo_diskbb_simpl_diskbb_"$jobnum".log`
        cstat1=`awk '/Total fit statistic/ {cstat=$4}; END {print cstat}' tbfeo_diskbb_simpl_diskbb_gabs_"$jobnum".log`

        rm tbfeo_diskbb_simpl_diskbb_"$jobnum".log
        rm tbfeo_diskbb_simpl_diskbb_gabs_"$jobnum".log

        echo $dirnum, $centroid_keV , $cstat0 , $cstat1 >> $RESULTS


        ((centroid_eV+=$linestep))

    done

    rm fake"$jobnum"*
    rm tbfeo_diskbb_simpl_diskbb_fak_"$jobnum".xcm
    rm sim_detlev_"$jobnum".xco

	((dirnum++))
done
