#!/bin/bash

script="average_using1ColumnAsId.awk"
if [ -n "$1"  ]
then
	samplepath=$1
else
	samplepath=""
fi

if [ -n "$2"  ]
then
	outputpath=$2
else
	outputpath=""
fi

if [ -n "$3" ]
then
	cost=$3
else
	cost=0.10
fi

declare -a plr=("0.00" "0.10" "0.25" "0.90")
arraylength=${#plr[@]}

for mt in 0 1 2 3
do

	for (( i=0; i<${arraylength}; i++ ));
	do
		awk -f $script ${samplepath}avInvXmh_*cost${cost}*plr${plr[$i]}*MT${mt}_1* > ${outputpath}averavInvXmh_cost${cost}_plr${plr[$i]}_MT${mt}.dat
	done
done
