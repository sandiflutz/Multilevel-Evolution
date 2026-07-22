#!/bin/bash

fname="runfindslope.txt"

par1=$(awk '{printf("%s",$1);}' "${fname}")
par2=$(awk '{printf("%s",$2);}' "${fname}")

./applyFindTransToSample_createCostxSlope.sh ${par1} ${par2}
