#!/bin/bash

dir=$1
corename=$2

n=0
for file in ${dir}/${corename}*
do
	echo "$n"
	cost=${dir##*cost}
	cost=${cost%/}
	fname=${file##*/}
	fname=${fname%.dat}
	fimage=$(./createFilesForAnim_Transitions.sh "${fname}" ${cost})
	n=$((n+1))
done

(
    cd "$dir" || exit 1
	gnuplot *.gp
 	./makeavi.sh 1 "${fimage}"
)

rm ${dir}/*.gp
rm ${dir}/*.png
