#!/bin/bash

dir=$1
corename=$2
output_slope=sample_slope_${dir}.dat

> $dir/"$output_slope"

echo "$output_slope"


n=0
for file in ${dir}/${corename}*
do
	awk -v PRINT_AxT=0 -v PRINT_AMAX=1 -v samplewindow=100 -v file_index=$n -f findTransition_minsquares.awk "$file" >> $dir/${output_slope}
	n=$((n+1))        
done
