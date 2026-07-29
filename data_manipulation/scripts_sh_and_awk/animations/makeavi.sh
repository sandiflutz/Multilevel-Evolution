#!/bin/bash

if [ -n "$1" ]
then
	fps=$1
else
	fps=2
fi

if [ -n "$2" ]
then
	par=$2
fi

if [ -n "$3" ]
then
	target=$3.avi
else
	target='t.avi'
fi

type=png
mencoder "mf://*$par*.$type" -mf fps=$fps:type=$type -ovc lavc -lavcopts vcodec=mpeg4 -o $target
