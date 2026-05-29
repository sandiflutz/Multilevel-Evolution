network=$1 


if [[ "$network" == "CG" ]]
then
	name=costXw_CG_L100_Ty101_Kh500_Gh100_Bv1e-3_mu1e-2_mb1e-6.dat
	./createCostxInvFile.sh $name 1 3
	cp $name ~/docs/work/France/research/prog/project/data_manipulation/gp_scripts/
elif [[ "$network" == "SL" ]]
then
	name=costXw_SL_L100_Ty101_Kh500_Gh100_Bv1e-3_mu1e-2_mb1e-6_mh4_rmh1_plr0.00.dat
	./createCostxInvFile.sh $name 1 11
	cp  $name ~/docs/work/France/research/prog/project/data_manipulation/gp_scripts/
else
	echo "Error. Expected syntax is: ./runCreateCostxInvFile.sh <network type(=CG,SL)>"
fi
