#!/bin/bash

if [ -n "$1" ]
then
        parfile=$1
        read net lsize types kh gh bv mu mb mh rmh plr < <(grep -v '^#' $parfile | grep -v '^$')
        saida=avertranstimeXcostL_${lsize}_Ty${types}_Kh${kh}_Gh${gh}_Bv${bv}_mu${mu}_mb${mb}.dat
else
        saida=avertranstimeXcost.dat
fi

if [ -n "$2" ]
then
	maxcol=$2
else
	maxcol=7
fi

echo "#1:cost 2:<t0> 3:stdt0 4:<tf> 5:stdtf 6:<a> 7:stda 8:<b> 9:stdb 10:<w0> 11:stdw0 12:<<wf>> 13:stdwf 14:<<w(t>tf)>t> 15:std(<w(t>tf)>t) 16:<tf-t0> 17:std(tf-t0)" > "$saida"

for dir in cost*/
do
	cost=$(basename "$dir" | sed 's/cost//')
	awk -v param=$cost -v maxcol=$maxcol -f average_columns-transtimes.awk ${dir}sample_transtimes_cost${cost}.dat >> "${saida}"
done
