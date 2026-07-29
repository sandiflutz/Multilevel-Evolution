#!/bin/bash

#Important Paths

projectpath="$(cd ../../../ && pwd)/"
datamanippath="$(cd ../../ && pwd)/"
currentdir="$(pwd)/"

address="alutz@sauron:/home/oar-jobs/alutz/project-multilevel/data_manipulation/"
outputaddress=$1
samplebasename=$2

scp ${address}${samplebasename}* ${outputaddress}

