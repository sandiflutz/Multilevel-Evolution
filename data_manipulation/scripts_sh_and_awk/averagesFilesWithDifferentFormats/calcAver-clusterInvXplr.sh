#!/bin/bash

#Entries are: 
#	-1: path for the sample files (default is the path for the data_manipulation directory)
#	-2: path for the output files (default is the path for the data_manipulation directory)
#	-3: cost value in the simulations (default is cost=0.1)
#	-4: Mh if the calculations are for only 1 value (default is Mh=1, 50, 100 and 500)
#	-5: MT if  the calculations are for only 1 value (default is MT=0, 1, 2 and 3)

#Important Paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="$(cd ../ && pwd)/"
scriptspath="${datamanippath}/scripts_sh_and_awk/"

#awk scripts used
scriptAverage="average.awk"
scriptSpace="createSpaceWhenFirstColumnChange.awk"

#base file name
basesamplename="clusterInvXplr"

#entries
if [ -n "$1" ]
then
	samplepath="$1"
else
	samplepath="${datamanippath}/"
fi
if [ -n "$2" ]
then
	outputpath="$2"
else
	outputpath="${datamanippath}/"
fi
if [ -n "$3" ]
then
	cost="$3"
else
	cost="0.10"
fi

if [ -n "$4" ]
then
	mh="$4"
	if  [ -n "$5" ]
	then
		mt="$5"
		awk -f $scriptAverage ${samplepath}${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
		awk -f ${scriptspath}${scriptSpace} ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
		rm ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
		awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}averinv_aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
	else
		for mt in 0 1 2 3
		do
			awk -f $scriptAverage ${samplepath}${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
			awk -f ${scriptspath}${scriptSpace} ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
			rm ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
			awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}averinv_aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
		done
	fi
else

        for mt in 0 1 2 3
        do
                for mh in 1 50 100 500
                do
                        awk -f $scriptAverage ${samplepath}${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
                        awk -f ${scriptspath}${scriptSpace} ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
                        rm ${samplepath}tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
			awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' ${outputpath}aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat > ${outputpath}averinv_aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat
                done

        done
fi
