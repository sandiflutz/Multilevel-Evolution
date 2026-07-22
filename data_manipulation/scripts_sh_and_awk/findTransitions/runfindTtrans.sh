#!/bin/bash

fname="runfindttrans.txt"

par1=$(awk '{printf("%s",$1);}' "${fname}")
par2=$(awk '{printf("%s",$2);}' "${fname}")
par3=$(awk '{printf("%s",$3);}' "${fname}")

./applyFindTransToSample_createCostxTtrans.sh ${par1} ${par2} ${par3}
