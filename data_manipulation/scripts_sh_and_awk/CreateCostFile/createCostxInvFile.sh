#!/bin/bash

awk1=average_sample_wXt.awk
awk2=average_columns.awk
if [ -n "$1" ]
then
        saida=$1
else
	saida="outputfile.dat"
fi
> "$saida"

for dir in cost*/
do
    cost=$(basename "$dir" | sed 's/cost//')
    dir="${dir%/}"

    awk -v param=$cost -v paramname=$dir -f $awk1 $dir/averInvXt_*.dat > $dir/sample_averInvXt_${dir}.dat
    awk -f $awk2 $dir/sample_averInvXt_${dir}.dat >> $saida

done

