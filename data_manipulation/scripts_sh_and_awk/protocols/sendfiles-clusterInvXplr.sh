#!/bin/bash

#Important paths

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
gpfilespath="$(cd ../ && pwd)/"
currentdir="$(pwd)/"
pathimg="$(cd ${projectpath} && cd ../Figures/main_images/ && pwd)/"
pathimgdata="${pathimg}data_and_script_files/"

pathcp="${pathimg}clusterinvdistXplr/"
pathcpdata="${pathcp}data_and_script_files/"

#Base file names

basefname="clusterInvXplr_"
basedatname="averclusterInvXplr_SL_L100_Ty101_Kh500"
basedatnamew="averinv_averclusterInvXplr"
echo ""

echo "Copying image files of the frac. of low investment clusters for cost=0.1 and MT=1 and 2"
cp lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT*.gp lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT*_fm lowinv_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_mh*.0_rmh1_MT*.dat $pathcpdata
cp lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT1.eps lowinv_multipleMh_clfracXplr_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1_MT2.eps $pathcp
echo "Copying image files of the heatmap multiplot for cost=0.1 and MT=1 and 2"
cp clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1MT12.gp clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1MT12_fm $pathcpdata
cp clusterInvXplr_MultipleMh_SL_L100_Ty101_Kh500_cost0.10_mu1e-2_mb1e-6_rmh1MT12.eps $pathcp

for mt in 0 1 2 3
do
	for mh in 1 50 100 500
	do
		echo ""
		echo "Copying image files of the heatmap for cost=0.1, MT=${MT} and Mh=${mh}"
		cp -v ${basefname}*mh${mh}.0*MT${mt}.eps $pathcp
		cp -v ${basefname}*mh${mh}.0*MT${mt}.gp ${basefname}*mh${mh}.0*MT${mt}_fm ${basedatname}*mh${mh}.0*MT${mt}.dat ${basedatnamew}*mh${mh}.0*MT${mt}.dat $pathcpdata
	done
done
