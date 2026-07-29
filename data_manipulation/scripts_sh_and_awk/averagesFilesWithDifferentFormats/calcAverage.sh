#!/bin/bash
if [ -n "$1" ]
then
        corename=$1
	awk -f average.awk $corename\_1*.dat > aver$corename\.dat
else
        echo 'Syntax: ./calcAverage.sh <core name>'
fi

