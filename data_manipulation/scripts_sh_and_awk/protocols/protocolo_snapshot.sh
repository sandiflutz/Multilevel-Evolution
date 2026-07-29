#!/bin/bash
fps=$1
corefile=$2
output=$3
mt=$4
pathcp="/home/sandi/docs/work/France/research/prog/Figures/main_images/snapshots/MT${mt}/"
pathcpdatgp="${pathcp}dat_and_gp/"
gnuplot ${corefile}*.gp
./makegif.sh ${fps} ${corefile} ${output}
cp ${output}.gif ${pathcp}
if [[ -n "$5" ]]
then
	if [[ -n "$6" ]]
	then
		id=$6
		cp ${corefile}*${id}.gp ${corefile}*${id}.dat ${pathcpdatgp}
		cp ${corefile}*${id}.png ${pathcp}
	fi
	
	path=$5
	cp ${output}.gif ${path}
	mv ${corefile}*.gp ${corefile}*.dat ${path}/data/
fi

rm ${corefile}*.png

