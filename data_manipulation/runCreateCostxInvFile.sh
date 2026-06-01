parfile=$1 

read net lsize types kh gh bv mu mb mh rmh plr < <(grep -v '^#' $parfile | grep -v '^$')

if [[ $net -eq 0 ]]
then
	name=costXw_CG_L${lsize}_Ty${types}_Kh${kh}_Gh${gh}_Bv${bv}_mu${mu}_mb${mb}.dat
	./createCostxInvFile.sh $name 1 3
	cp $name ~/docs/work/France/research/prog/project/data_manipulation/gp_scripts/
elif [[ $net -eq 1 ]]
then
	name=costXw_SL_L${lsize}_Ty${types}_Kh${kh}_Gh${gh}_Bv${bv}_mu${mu}_mb${mb}_mh${mh}_rmh${rmh}_plr${plr}.dat
	./createCostxInvFile.sh $name 1 11
	cp  $name ~/docs/work/France/research/prog/project/data_manipulation/gp_scripts/
fi
