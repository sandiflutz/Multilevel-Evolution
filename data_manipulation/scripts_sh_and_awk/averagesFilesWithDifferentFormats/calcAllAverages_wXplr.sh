#!/bin/bash

script=average_using1ColumnAsId.awk

if [ -n "$1" ]
then
	samplepath="$1"
else
	samplepath="../../gp_scripts/Files_mhXplrXw/allaverInvXplrfiles/ "
fi

if [ -n "$2" ]
then
	outputpath="$2"
else
	outputpath="../../gp_scripts/Files_mhXplrXw/allaverInvXplrfiles/ "
fi

if [ -n "$3" ]
then
	cost="$3"
else
	cost="0.10"
fi

for mh in 1 50 100 500 1000
do
	output=averavInvXplr_cost${cost}_mh${mh}.0_MT0.dat
	awk -v MULT=100 -v MAXCOLS=2 -v IDCOL=1 -f $script ${samplepath}/averinv_*cost${cost}*mh${mh}.0*MT0* > ${outputpath}/$output
done
for mh in 1 50 100 500 1000
do
	for mt in 1 2 3
	do
		output=averavInvXplr_cost${cost}_mh${mh}.0_MT${mt}.dat
		awk -v MULT=100 -v MAXCOLS=2 -f $script ${samplepath}/averinv_*cost${cost}*mh${mh}.0*MT${mt}* ${samplepath}/avInvXplr*cost${cost}*mh${mh}.0*MT${mt}* > ${outputpath}/$output
	done
done
awk -v MULT=100 -v MAXCOLS=2 -f $script ${samplepath}/avInvXplr*cost${cost}*mh5000.0*MT2* > ${outputpath}/averavInvXplr_cost${cost}_mh5000.0_MT2.dat

