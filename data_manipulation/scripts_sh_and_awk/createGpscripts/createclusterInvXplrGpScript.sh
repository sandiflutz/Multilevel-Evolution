#!/bin/bash

#Important paths
projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"
gpfilespath="${datamanippath}gp_scripts/"
scriptspath="${datamanippath}scripts_sh_and_awk/"
averscriptpath="${scriptspath}averagesFilesWithDifferentFormats/"
protocolpath="${scriptspath}protocols/"
awkfilterspath="${scriptspath}awkFilters/"
creategpscrpath="${scriptspath}createGpscripts/"
heatmapath="${gpfilespath}Files_clinvdistXplr/"

filespath="${heatmapath}"

#Entrie are: 
#1- The data file name (without the extension) 
#2- The output file name (whithout the extension)
#3- MT (migration type) used in the simulations
#4- Mh (migration coefficient) used in the simulations
#5- cost used in the simulations


if [ -n "$1"  ]
then
	splotname=$1
	if [ -n "$2"  ] ; then
		output=$2
		MT=$3
		Mh=$4
		cost=$5
		if [[ $MT -eq 0 ]] ; then
			mt="homogeneous"
		elif [[ $MT -eq 1 ]] ; then
			mt="localcrowdness"
		elif [[ $MT -eq 2 ]] ; then
			mt="localinvestment"
		else
			mt="localcrowdnessinvestment"
		fi
		####gnuplot script#################
		echo "set term post eps enha color 18" > ${filespath}${output}.gp
		echo "set output'${filespath}${output}_fm.eps'" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "# --- Main plot settings ---" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "set xr[0:1]" >> ${filespath}${output}.gp
		echo "set yr[-0.001:0.99]" >> ${filespath}${output}.gp
		echo "set xtics 0.2" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp 
		echo "set pm3d map" >> ${filespath}${output}.gp
		printf '%s\n' 'set palette defined (0.0 "#000000",0.01 "#070551",0.02 "#0c0887",0.03 "#4b03a1",0.04 "#7d03a8",0.05 "#a82296",0.06 "#b32677",0.08 "#cb4679",0.1 "#d45266",0.11 "#e56b5d",0.12 "#f89441",0.13 "#fdc328",0.14 "#ffdb46",0.15 "#f0f921")' >> ${filespath}${output}.gp
		echo "set cbrange[0:0.15]" >> ${filespath}${output}.gp
		echo "set cbtics 0.05" >> ${filespath}${output}.gp
		echo "set cblabel '{/=15 clfrac}' offset 1" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "set xl'{/=25 Plr}'" >> ${filespath}${output}.gp
		echo "set yl'{/=25 AverageInvestment}' offset -1" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "set label 1 '{/=15 cost${cost}Mh${Mh}}' at 0.003,0.17 front" >> ${filespath}${output}.gp
		echo "set label 2 '{/=15 ${mt}}' at 0.002,0.1 front" >> ${filespath}${output}.gp
		printf '%s\n\n' 'set ytics ("0" 0, "0.2" 0.2, "0.4" 0.4, "0.6" 0.6,"0.8" 0.8,"1" 0.99)' >> ${filespath}${output}.gp
		echo "# --- Multiplot with inset ---" >> ${filespath}${output}.gp
		echo "set multiplot" >> ${filespath}${output}.gp
		echo "unset key" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "# --- Main plot ---" >> ${filespath}${output}.gp
		echo "splot '${filespath}${splotname}' u 1:2:3" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "# --- Inset (top-left) ---" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "set origin 0.1,0.32" >> ${filespath}${output}.gp
		echo "set size 0.4,0.4" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		echo "unset label" >> ${filespath}${output}.gp
		echo "unset pm3d" >> ${filespath}${output}.gp
		echo "unset colorbox" >> ${filespath}${output}.gp
		echo "unset cblabel" >> ${filespath}${output}.gp
		echo "unset arrow" >> ${filespath}${output}.gp
		echo"" >> ${filespath}${output}.gp
		printf '%s\n' 'set border lc rgb "white"' >> ${filespath}${output}.gp
		printf '%s\n' 'set tics textcolor rgb "white"' >> ${filespath}${output}.gp
		printf '%s\n' 'set xtics font ",10"' >> ${filespath}${output}.gp
		printf '%s\n' 'set ytics font ",10"' >> ${filespath}${output}.gp
		echo "set xr[0:1]" >> ${filespath}${output}.gp
		echo "set yr[0.45:1]" >> ${filespath}${output}.gp
		echo "set ytics 0.25" >> ${filespath}${output}.gp
		echo "set xtics 0.2" >> ${filespath}${output}.gp
		printf '%s' "set yl 'averageinvestment' offset 3,0.5 textcolor rgb " '"white"' >> ${filespath}${output}.gp
		echo "" >> ${filespath}${output}.gp
		printf '%s' "set xl 'plr' offset 0,0.8 textcolor rgb" '"white"' >> ${filespath}${output}.gp
		echo "" >> ${filespath}${output}.gp
		echo "plot '${filespath}averinv_${splotname}' u 1:2 w l lw 2 lc rgb 'red' t ''" >> ${filespath}${output}.gp
		echo "unset multiplot" >> ${filespath}${output}.gp
		####psfrag script#################
		printf '%s\n' '\psfrag{AverageInvestment}{\textcolor{black}{\hspace{-2mm}\large Average Investment}}' > ${filespath}${output}_fm
		printf '%s\n' '\psfrag{averageinvestment}{\textcolor{white}{\hspace{-2mm}\tiny Average Investment}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{Plr}{\textcolor{black}{\large $P_{lr}$}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{plr}{\textcolor{white}{\tiny $P_{lr}$}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{clfrac}{\textcolor{black}{\hspace{-13mm} Frac. of Clusters}}' >> ${filespath}${output}_fm
		printf '%s' '\psfrag{cost' "${cost}Mh${Mh}" '}{\fcolorbox{black}{white}{\tiny $\gamma=' "${cost}" '$, $M_h=' "${Mh}" '$}}' >> ${filespath}${output}_fm
		printf '\n%s\n' '\psfrag{homogeneous}{\fcolorbox{black}{white}{\tiny Mig.: Homogeneous}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{localcrowdness}{\fcolorbox{black}{white}{\tiny Mig.: Local Crowdness}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{localinvestment}{\fcolorbox{black}{white}{\tiny Mig.: Local Investmemt}}' >> ${filespath}${output}_fm
		printf '%s\n' '\psfrag{localcrowdnessinvestment}{\fcolorbox{black}{white}{\tiny Mig.: Local Crowd. and Inv.}}' >> ${filespath}${output}_fm
	else
		echo"use the syntax:" 
		echo"	./createclusterInvXplrGpScript.sh <file name for the splot map> <output name without extension> <mt> <mh> <cost>"
		echo"Make sure the files are in the current folder." 
		echo"" 
	fi

else
	echo"use the syntax:" 
	echo"	./createclusterInvXplrGpScript.sh <file name for the splot map> <output name without extension> <mt> <mh> <cost>"
	echo"" 
fi
