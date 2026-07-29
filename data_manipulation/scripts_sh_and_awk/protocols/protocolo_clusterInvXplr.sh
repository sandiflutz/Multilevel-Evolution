#!/bin/bash

#Important Paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="$(cd ../ && pwd)/"
scriptpath="${datamanippath}/scripts_sh_and_awk/averagesFilesWithDifferentFormats/"

#script used
script="./calcAver-clusterInvXplr.sh"

#entry is cost value (default is 0.1)

if [ -n "$1" ]
then
	cost="$1"
else
	cost="0.10"
fi

#actions
echo ""
echo "Calculating averages for the heatmap of the fraction of clusters for each average investment bin and fraction of long-range migration Plr"
$(cd $scriptpath && ./$script $datamanippath $currentdir ${cost})
echo ""
echo "Running gnuplot scripts"
./rungp.sh

echo ""
echo "Creating Low Investment against Plr files"
./createLowInvClFracFiles.sh
echo ""
echo "Running gnuplot scripts"
gnuplot clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1MT12.gp
gnuplot lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT1.gp
gnuplot lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT2.gp
echo ""
echo "Running fragmaster"
fragmaster

echo ""
echo "Copying final files to the main images directory"
./copyeps_to_figures.sh

