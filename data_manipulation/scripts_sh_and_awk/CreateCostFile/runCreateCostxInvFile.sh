#!/bin/bash

#Important paths

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
scriptspath="${datamanippath}/scripts_sh_and_awk/"

costscriptpath="${scriptspath}CreateCostFile/"

#Script used

script=createCostxInvFile.sh

#Output file base name

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
	
$( cd ${costscriptpath} && ./${script} $outputbase $plr $mh )
