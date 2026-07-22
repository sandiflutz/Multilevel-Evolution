#!/bin/bash
if [ -n "$1" ]
then
	fps=$1
else
	fps=2
fi
par=$2
if [ -n "$3" ]
then
	target=$3.gif
else
	target=anim$par.gif
fi
if [ -n "$4" ]
then
	scale=$4
else
	scale=500
fi
ffmpeg -framerate $fps -pattern_type glob -i "$par*.png" -vf "scale=$scale:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse"  -loop 0 "$target"
