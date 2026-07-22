#!/bin/bash

if [ -n "$1" ]
then
	parfile=$1
	read net lsize types kh gh bv mu mb mh rmh plr < <(grep -v '^#' $parfile | grep -v '^$')
	saida=averslopeXcost_L${lsize}_Ty${types}_Kh${kh}_Gh${gh}_Bv${bv}_mu${mu}_mb${mb}_plr${plr}.dat
else
	saida=averslopeXcost.dat
fi

if [ -n "$2" ]
then
	corename=$2
else
	corename=averInvXt
fi

> "${saida}"

for dir in cost*/
do
	cost=$(basename "$dir" | sed 's/cost//')
	dir="${dir%/}"
	samplefile=$(./makeSampleFile_slope.sh "${dir}" "${corename}")
       	awk -v param=$cost -f average_columns.awk $dir/${samplefile} >> "${saida}"

done
