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
heatmapath="${gpfilespath}Files_clinvdistXplr/"

#awk script used

awkscript=createClFracXplr_for_a_InvBIN.awk
scriptpath="${awkfilterspath}"

#Entries (1: path for the data files used in the gnuplot script; 2: cost)

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

#File base names
samplebasename="averclusterInvXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6"
outputbasename="lowinv_clfracXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6"

for mt in 1 2
do
	for mh in 1 50 100 500
	do
		samplename="${samplebasename}_mh${mh}.0_rmh1_MT${mt}.dat"
		output="${outputbasename}_mh${mh}.0_rmh1_MT${mt}.dat"
		awk -v INVBIN=0.01 -f $scriptpath$awkscript $aversamplepath$samplename > $aversamplepath$output
	done
done
