#!/bin/bash
base_name=$1
fps=$2
gnuplot $base_name*.gp 
./makeavi.sh $fps $base_name animavi_$base_name
#rm $base_name*.gp $base_name*.dat $base_name*.png
