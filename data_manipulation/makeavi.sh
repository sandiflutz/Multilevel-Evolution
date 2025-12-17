#!/bin/bash
if [ -n "$1" ]
then
	fps=$1
else
	fps=2
fi
type=png
par=$3
if [ -n "$2" ]
then
	target=$2.avi
else
	target='t.avi'
fi
mencoder "mf://*$par*.$type" -mf fps=$fps:type=$type -ovc lavc -lavcopts vcodec=mpeg4 -o $target
