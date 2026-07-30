#!/bin/bash

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
scriptspath="${datamanippath}scripts_sh_and_awk/"
protocolpath="${scriptspath}protocols/"
awkfilterspath="${scriptspath}awkFilters/"
averscriptpath="${scriptspath}averagesFilesWithDifferentFormats/"
creategpscrpath="${scriptspath}createGpscripts/"

pathimg="$(cd ${projectpath} && cd ../Figures/main_images/ && pwd)/"
pathimgdata="${pathimg}data_and_script_files/"
averwXplrpath="${datamanippath}averInv_Plr/"
plrpath="${averwXplrpath}allaverInvXplrfiles/"
heatmapath="${gpfilespath}Files_clinvdistXplr/"
plrgpfilespath="${gpfilespath}Files_mhXplrXw/"

#Entries

if [ -n "$1"  ]
then
	cost="$1"
else
	cost="0.10"
fi

#File base names

basefname="averavInvXplr_multipleMh_SL_L100_Ty101_cost${cost}_mu1e-2_mb1e-6_rmh1_MT"
datbasename="averavInvXplr_cost${cost}"
samplebasename1="avInvXplr"
samplebasename2="clusterInvXplr"

#Actions

echo "Running scpLab.sh"
./scpLab.sh ${averwXplrpath} ${samplebasename1}*
./scpLab.sh ${heatmapath} ${samplebasename2}*
echo ""
echo "Calculating averages for the heatmap files of the frac. of cluster for each investm. bin and Plr with calcAver-clusterInvXplr_tmp.sh"
$(cd $averscriptpath && ./calcAver-clusterInvXplr.sh $datamanippath $plrgpfilespath $cost )

echo "Bringing every file that has investment against Plr together"
cp ${averwXplrpath}avInvXplr* ${plrpath}
cp ${heatmapath}averinv_* ${plrpath}

echo "Calculating averages using all files that have investment as a function of Plr using calcAllAverages_wXplr.sh"

$(cd $averscriptpath && ./calcAllAverages_wXplr.sh $plrpath $plrpath $cost)

echo "Bringing output files to the current directory"
cp ${plrpath}${datbasename}* ${plrgpfilespath}

echo "Running gnuplot scripts, fragmaster and copying final files to directory with the main images"
for mt in 0 1 2 3
do
	fname=${basefname}${mt}
	$( cd $plrgpfilespath ; gnuplot ${fname}.gp )
	$( cd $plrgpfilespath ; export PATH=$PATH:/usr/bin/ fragmaster)
	cp ${plrgpfilespath}${fname}.gp ${plrgpfilespath}${fname}_fm ${plrgpfilespath}${datbasename}*MT${mt}.dat $pathimgdata
	cp ${plrgpfilespath}${fname}.eps $pathimg
done

