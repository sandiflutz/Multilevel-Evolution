#!/bin/bash
#parfile is the file with the main parameters used in the simulation that are going to the output file name
#core name is the firs part of the name of the sample files, which is the same for every file in the sample
#maxcol is to the side in which column of the sample file stop calculating average and standart deviation (since columns at the end may not relevant)


if [ -n "$1" ]
then
	parfile=$1
	read net lsize types kh gh bv mu mb mh rmh plr < <(grep -v '^#' $parfile | grep -v '^$')
	saida=avertranstimeXcost_L${lsize}_Ty${types}_Kh${kh}_Gh${gh}_Bv${bv}_mu${mu}_mb${mb}_plr${plr}.dat
else
	saida=avertranstimeXcost.dat
fi

if [ -n "$2" ]
then
	corename=$2
else
	corename=averInvXt
fi
if [ -n "$3" ]
then
	maxcol=$3
else
	maxcol=7
fi

> "${saida}"

for dir in cost*/
do
	cost=$(basename "$dir" | sed 's/cost//')
	dir="${dir%/}"
	samplefile=$(./makeSampleFile_ttrans.sh "${dir}" "${corename}")
       	awk -v param=$cost -v maxcol=$maxcol -f average_columns-transtimes.awk $dir/${samplefile} >> "${saida}"

done
