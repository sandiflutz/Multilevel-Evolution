#!/bin/bash

#Important paths

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
currentdir="$(pwd)/"
scriptspath="${datamanippath}scripts_sh_and_awk/"
averscriptpath="${scriptspath}averagesFilesWithDifferentFormats/"
protocolpath="${scriptspath}protocols/"
awkfilterspath="${scriptspath}awkFilters/"
creategpscrpath="${scriptspath}createGpscripts/"

pathimg="$(cd ${projectpath} && cd ../Figures/main_images/ && pwd)/"
pathimgdata="${pathimg}data_and_script_files/"
heatmapath="${gpfilespath}Files_clinvdistXplr/"

pathcp="${pathimg}clusterinvdistXplr/"
pathcpdata="${pathcp}data_and_script_files/"

#Entries

if [ -n "$1" ] ; then
	filespath="$1"
else
	filespath="${heatmapath}"
fi
if [ -n "$2" ] ; then
	cost=$2
else
	cost="0.10"
fi

#Base file names

basenameclinv="clusterInvXplr_"
basedatnameclinv="averclusterInvXplr_SL_L100_Ty101_Kh500"
basedatnameaverinv="averinv_averclusterInvXplr"
basenamelowinv="lowinv"
gpmultiplotname="clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost${cost}_mu1e-2_mb1e-6_rmh1_MT12"

#Actions

echo "Copying image files of the frac. of low investment clusters for cost=0.1 and MT=1 and 2"
cp $filespath${basenamelowinv}*.gp $filespath${basenamelowinv}*_fm $filespath${basenamelowinv}*.dat $pathcpdata
echo "Copying image files of the heatmap multiplot for cost=0.1 and MT=1 and 2"
cp $filespath${gpmultiplotname}_fm  $filespath${gpmultiplotname}.gp $pathcpdata
cp $filespath${gpmultiplotname}.eps $pathcp

for mt in 0 1 2 3
do
	for mh in 1 50 100 500
	do
		echo "Copying image files of the heatmap for cost=0.1, MT=${mt} and Mh=${mh}"
		cp $filespath${basefname}*mh${mh}.0*MT${mt}.eps $pathcp
		cp $filespath${basefname}*mh${mh}.0*MT${mt}.gp $pathcpdata
		cp ${filespath}${basefname}*mh${mh}.0*MT${mt}_fm $pathcpdata
		cp $filespath${basedatname}*mh${mh}.0*MT${mt}.dat $pathcpdata
		cp $filespath${basedatnamew}*mh${mh}.0*MT${mt}.dat $pathcpdata
	done
done
