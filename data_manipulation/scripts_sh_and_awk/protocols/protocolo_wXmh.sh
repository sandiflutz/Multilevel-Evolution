#!/bin/bash

#Important paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
currentdir="$(pwd)/"
scriptpath="${datamanippath}scripts_sh_and_awk/averagesFilesWithDifferentFormats/"
pathimg="$(cd ${projectpath} && cd ../Figures/main_images/ && pwd)/"
pathimgdata="${pathimg}data_and_script_files/"

mhpath="${datamanippath}averInv_Mh/averInvXmh/allfiles/"
mhgpdatapath="${gpfilespath}Files_mhXplrXw/"

#entry is cost value (default is 0.1)
if [ -n "$1" ]
then
	cost=$1
else
	cost="0.10"
fi

#base names for files
basename="averwXmh_multiplotMT_multiplePlr_SL_L100_Ty101_cost${cost}_mu1e-2_mb1e-6_rmh1"
basedatname="averavInvXmh_cost${cost}"

#actions

samplebasename="avInvXmh"

echo "Running scpLab.sh"
./scpLab.sh ${mhpath} ${samplebasename}

echo ""
echo "Calculating averages"
$(cd $scriptpath && ./calcAverages_wXmh.sh $mhpath $mhpath $cost)
echo "Bringing output files to the current directory"
cp ${mhpath}averavInvXmh*.dat $mhgpdatapath

echo "Running gnuplot script and fragmaster"
echo ""
$( cd $mhgpdatapath ; gnuplot ${basename}.gp)
$( cd $mhgpdatapath ; export PATH=$PATH:/usr/bin/ fragmaster )

echo "Copying image files to the directory with the main images"
echo ""
cp ${mhgpdatapath}${basename}.gp ${mhgpdatapath}${basename}_fm ${mhgpdatapath}${basedatname}_plr0.*_MT*.dat $pathimgdata
cp ${mhgpdatapath}${basename}.eps $pathimg

