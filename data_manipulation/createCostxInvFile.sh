#!/bin/bash

awk1=average_sample_wXt.awk
awk2=average_time.awk
if [ -n "$1" ]
then
        saida=$1
else
	saida="costXw_L100_Ty101_Kh500_net1_Gh100_Bv1e-3_mu1e-2_mb1e-6_mh4_rmh1_plr0.00.dat"
fi
> "$saida"

if [ -n "$2" ]
then
        col1=$2
else
	col1=11
fi
if [ -n "$3" ]
then
        col2=$3
else
	col2=1
fi

for dir in cost*/
	do
    cost=$(basename "$dir" | sed 's/cost0.//')

    awk -v meascol=$col1 -v paramcol=$col2 -f $awk1 $dir/averInvXt_*.dat > $dir/sample_averInvXt_cost$cost.dat
    awk -f $awk2 $dir/sample_averInvXt_cost$cost.dat >> $saida

done

