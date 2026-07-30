#!/bin/bash

#Important Paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
heatmapath="${gpfilespath}Files_clinvdistXplr/"

#Entries (1: path for the gnuplot script; 2:cost)

if [ -n "$1" ] ; then
        gppath="$1"
else
        gppath="${heatmapath}"
fi

if [ -n "$2" ] ; then
        cost="$2"
else
        cost="0.10"
fi

#Actions
for mt in 0 1 2 3
do
	for mh in 1 50 100 500
	do
		gpscript="clusterInvXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_mh${mh}.0_rmh1_MT${mt}.gp"
		$( cd ${gppath} ; gnuplot ${gpscript} )
	done
done

$( cd $heatmapath ; export PATH=$PATH:/usr/bin/ fragmaster )
