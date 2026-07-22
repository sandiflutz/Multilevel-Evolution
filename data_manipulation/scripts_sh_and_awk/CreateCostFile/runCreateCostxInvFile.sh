#!/bin/bash

if [ -n "$1" ]
then
	fname=$1 
else
	fname=outputfile.dat
fi
if [ -n "$2" ]
then
	fpath=$2 
else
	fpath="*"
fi
	
./createCostxInvFile.sh ${fname}
if [ -n "$2" ]
then
	fpath=$2 
	cp $fname $fpath
fi
