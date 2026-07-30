#!/bin/bash

#Important Paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
scriptspath="${datamanippath}scripts_sh_and_awk/"
averscriptpath="${scriptspath}averagesFilesWithDifferentFormats/"
protocolpath="${scriptspath}protocols/"
awkfilterspath="${scriptspath}awkFilters/"
creategpscrpath="${scriptspath}createGpscripts/"

#Scripts used

awkscript=createClFracXplr_for_a_InvBIN.awk
scriptpath="${awkfilterspath}"

#Entries 
#1: path for the data files used in the gnuplot script
#2: cost 
#3:investment bin (from 0.01 to 0.99, with a difference of 0.01 between consecutive bins)

if [ -n "$1" ] ; then
        aversamplepath="$1"
else
        aversamplepath=$heatmapath
fi

if [ -n "$2" ] ; then
        cost=$2
else
        cost="0.10"
fi

if [ -n "$3" ] ; then
	invbin="$3"
else
	invbin="0.01"
fi

#File base names

samplebasename="averclusterInvXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6"
outputbasename="clfracXplr_invbin${invbin}_cost${cost}_mu1e-2_mb1e-6"

#Actions

for mt in 1 2
do
	for mh in 1 50 100 500
	do
		samplename="${samplebasename}_mh${mh}.0_rmh1_MT${mt}.dat"
		output="${samplebasename}_mh${mh}.0_rmh1_MT${mt}.dat"
		awk -v INVBIN=${invbin} -f $scriptpath$awkscript $samplepath$samplename >  $samplepath$output
	done
done
