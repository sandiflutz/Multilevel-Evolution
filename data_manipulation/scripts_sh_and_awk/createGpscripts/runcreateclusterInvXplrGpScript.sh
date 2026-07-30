#!/bin/bash

#Important paths
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

#Script used
script=createclusterInvXplrGpScript.sh
scriptpath="${creategpscrpath}"

#Entry

if [ -n "$1"  ]
then
	cost=$1
else
	cost="0.10"
fi

#Actions

for mt in 0 1 2 3
do
	for mh in 1 50 100 500
	do
		samplename="averclusterInvXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_mh${mh}.0_rmh1_MT${mt}.dat"
		output="clusterInvXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_mh${mh}.0_rmh1_MT${mt}"
		$( cd $creategpscrpath && ./"${script}" $samplename $output ${mt} ${mh} ${cost} )
	done
done
