#!/bin/bash

fname=$1
cost=$2

fttrans="cost${cost}/ttrans.dat"
fslope="cost${cost}/slope.dat"
corename=ttrans_${fname}
gpname="cost${cost}/${corename}.gp"

awk -v PRINT_AMAX=1 -f findTransition_minsquares.awk cost${cost}/${fname}.dat > "${fttrans}"
awk -v PRINT_AMAX=1 -v samplewindow=100 -f findTransition_minsquares.awk cost${cost}/${fname}.dat > "${fslope}"


t0=$(awk '{print $1}' "${fttrans}")
tf=$(awk '{print $2}' "${fttrans}")
amax=$(awk '{print $1}' "${fslope}")
bmax=$(awk '{print $2}' "${fslope}")

echo "set terminal pngcairo size 1080,810 enhanced font \"Helvetica,20\"" > "${gpname}"
echo "set output '${corename}.png'" >> "${gpname}"
echo "f(x)=a*x+b" >> "${gpname}"
echo "a=${amax}" >> "${gpname}"
echo "b=${bmax}" >> "${gpname}"
echo "set yr[0:1]" >> "${gpname}"
echo "set arrow 1 from ${t0},0 to ${t0},1 lt '-' lw 2 lc rgb 'black' nohead" >> "${gpname}"
echo "set arrow 2 from ${tf},0 to ${tf},1 lt '-' lw 2 lc rgb 'black' nohead" >> "${gpname}"
echo "plot '${fname}.dat' u 1:2 w l lw 2 lc rgb 'blue' t'',\\" >> "${gpname}"
echo "f(x) w l lw 2 lc rgb 'red' t''" >> "${gpname}"

#evince "${gpname}.eps" &

echo "${corename}"

rm "${fslope}"
rm "${fttrans}"
