#!/bin/bash

dir=$1
corename=$2
output_ttrans=sample_transtimes_${dir}.dat

> $dir/"$output_ttrans"

echo "$output_ttrans"

n=0
for file in ${dir}/${corename}*
do
	awk -v PRINT_AxT=0 -v PRINT_AMAX=1 -v samplewindow=0 -v file_index=$n -f findTransition_minsquares.awk "$file" >> $dir/${output_ttrans}
	n=$((n+1))        
done
