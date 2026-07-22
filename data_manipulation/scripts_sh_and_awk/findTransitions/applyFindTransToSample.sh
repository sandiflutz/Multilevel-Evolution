#!/bin/bash

if [ -n "$1" ]
then
        basename=$1
	outfile=transtimes_${basename}.dat
	> "$outfile" 

        if [ -n "$2" ]
        then
                sampleinterval=$2
        else
                sampleinterval=1
        fi
        if [ -n "$3" ]
        then
                samplewindow=$3
        else
                samplewindow=200
        fi

	for file in ${basename}
	do
		awk -v sampleinterval=$sampleinterval -v samplewindow=$samplewindow -f findTransition_minsquares.awk "$file" >> "$outfile"
	done

fi
