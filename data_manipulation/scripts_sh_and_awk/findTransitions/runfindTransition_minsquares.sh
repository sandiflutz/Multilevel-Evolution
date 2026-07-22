#!/bin/bash

if [ -n "$1" ]
then
	fname=$1
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
	awk -v sampleinterval=$sampleinterval -v samplewindow=$samplewindow -f findTransition_minsquares.awk $fname
	#awk -v sampleinterval=$sampleinterval -v samplewindow=$samplewindow -f findTransition_minsquares.awk $fname > inclination_$fname

fi

