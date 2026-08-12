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
averwXplrpath="${datamanippath}averInv_Plr/"
plrpath="${averwXplrpath}allaverInvXplrfiles/"

#script used
averscript=calcAver-clusterInvXplr.sh
lowinvscript=createLowInvClFracFiles.sh
rungpscrip=rungp.sh
sendfilesscript=sendfiles-clusterInvXplr.sh

lowinvscriptpath=$protocolpath

#Entry is cost value (default is 0.1)

if [ -n "$1" ]
then
	cost="$1"
else
	cost="0.10"
fi

#File base names

samplebasename="clusterInvXplr"

#Actions

echo "Running scpLab.sh"
./scpLab.sh ${heatmapath} ${samplebasename}*

echo "Bringing heat map sample files from $datamanippath "
cp $datamanippath${samplebasename}* $heatmapath
echo "Calculating averages for the heatmap of the fraction of clusters for each average investment bin and fraction of long-range migration Plr"
$(cd $averscriptpath && ./"${averscript}" $heatmapath $heatmapath $cost )
echo ""
echo "Running gnuplot scripts"
./"${rungpscrip}" $heatmapath $cost

echo ""
echo "Creating Low Investment against Plr files"
$( cd $lowinvscriptpath && ./"${lowinvscript}" $heatmapath $cost )
echo ""
echo "Running gnuplot scripts"
gpscr_clInvXplr="clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_rmh1_MT12.gp"
gpscr_lowinvMT1="lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_rmh1_MT1.gp"
gpscr_lowinvMT2="lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_rmh1_MT2.gp"

$( cd $heatmapath && gnuplot $gpscr_clInvXplr ) 
$( cd $heatmapath && gnuplot $gpscr_lowinvMT1 ) 
$( cd $heatmapath && gnuplot $gpscr_lowinvMT2 )
echo "Running fragmaster"
$(  cd $heatmapath ; export PATH=$PATH:/usr/bin/ fragmaster )

echo "Copying final files to the main images directory"
./"${sendfilesscript}" $heatmapath $cost
