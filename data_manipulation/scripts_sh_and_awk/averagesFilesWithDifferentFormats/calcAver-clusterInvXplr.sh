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
gpfilespath="${datamanippath}gp_scripts/"
scriptspath="${datamanippath}/scripts_sh_and_awk/"
averscriptpath="${scriptspath}averagesFilesWithDifferentFormats/"
protocolpath="${scriptspath}protocols/"
awkfilterspath="${scriptspath}awkFilters/"
creategpscrpath="${scriptspath}createGpscripts/"
heatmapath="${gpfilespath}Files_clinvdistXplr/"
averwXplrpath="${datamanippath}averInv_Plr/"
plrpath="${averwXplrpath}allaverInvXplrfiles/"

#awk scripts used
scriptAverage=average.awk
scriptSpace=createSpaceWhenFirstColumnChange.awk

#base file name
basesamplename="clusterInvXplr"

#Entries
if [ -n "$1" ] ; then
	samplepath="$1"
else
	samplepath="${heatmapath}"
fi

if [ -n "$2" ] ; then
	outputpath="$2"
else
	outputpath="${heatmapath}"
fi

if [ -n "$3" ] ; then
	cost="$3"
else
	cost="0.10"
fi

if [ -n "$4" ] ; then
	mh="$4"

	if  [ -n "$5" ] ; then
		mt="$5"
		tmpfile="tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
		output="aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
		outputspinoff="averinv_${outputbase}"
		awk -f $averscriptpath$scriptAverage $samplepath${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > $outputpath$tmpfile
		awk -f $scriptspath$scriptSpace $outputpath$tmpfile > $outputpath$output
		rm ${outputpath}${tmpfile}
		awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' $outputpath$output > $outputpath$outputspinoff
	else
		for mt in 0 1 2 3
		do
			tmpfile="tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
			output="aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
			outputspinoff="averinv_${output}"
			:q
			awk -f $averscriptpath$scriptAverage $samplepath${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > $outputpath$tmpfile
			awk -f $scriptspath$scriptSpace $outputpath$tmpfile > $outputpath$output
			rm $outputpath$tmpfile
			awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' $outputpath$output > $outputpath$outputspinoff
		done
	fi
else

        for mt in 0 1 2 3
        do
                for mh in 1 50 100 500
                do
			tmpfile="tmpaver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
			output="aver${basesamplename}_cost${cost}_mh${mh}.0_MT${mt}.dat"
			outputspinoff="averinv_${output}"
			awk -f $averscriptpath$scriptAverage $samplepath${basesamplename}*cost${cost}*mh${mh}.0*MT${mt}_1* > $outputpath$tmpfile
			awk -f $scriptspath$scriptSpace $outputpath$tmpfile > $outputpath$output
			rm $outputpath$tmpfile
			awk 'BEGIN{printf("#Plr <w> sample\n");}{if($2==0.){printf("%f %f %d\n",$1,$4,$8);}}' ${outputpath}${output} > ${outputpath}${outputspinoff}
                done

        done
fi
