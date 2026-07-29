#!/bin/bash

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
scriptspath="${datamanippath}/scripts_sh_and_awk/"
gpfilespath="${datamanippath}gp_scripts/"

costgppath="${gpfilespath}Files_wXcost/"
costscriptpath="${scriptspath}/CreateCostFile/"
awkscriptspath="${costscriptpath}"

samplepath="${datamanippath}averInv_MT0123/CI1/"

#awk scripts used

scriptawk1="average_sample_wXt.awk"
scriptawk2="average_columns.awk"

#File base names

samplebase="averInvXt"
allsamplesnamebase="sample_averInvXt"

if [ -n "$1" ]; then
	outputbase="$1"
else
	outputbase="costXw_SL_L100_Ty101_Kh500_mu1e-2_mb1e-6_mh1.0_rmh1"
fi
if [ -n "$2" ]; then
        plr="$2"
else
        plr="0.00"
fi

if [ -n "$3" ]; then
        mh="$3"
else
        mh="1"
fi

#Actions


for mt in 0 
do
	fullsamplepath="${samplepath}MT${mt}_Mh${mh}_plr${plr}/"
	output="${outputbase}_mh${mh}.0_plr0.${plr}0_MT${mt}.dat"
	> ${fullsamplepath}$output
	$( cd $fullsamplepath && 
		for dir in cost*
		do
			cost=$(basename "$dir" | sed 's/cost//')
			dir="${dir%/}"
			samplename="${allsamplesnamebase}_${cost}_mh${mh}_plr${plr}_MT${mt}.dat"
			echo "Calculating sample averages for MT=${mt} and cost=${cost}"
			awk -v param=$cost -v paramname=$dir -f ${awkscriptspath}$scriptawk1 $dir/${samplebase}*.dat > ${dir}/$samplename
			echo "Calculating time averages for MT=${mt} and ${cost}"
			awk -f ${awkscriptspath}$scriptawk2 ${dir}/$samplename >> ${output}
		done
	)
	cp ${fullsamplepath}${output} $costgppath
done
