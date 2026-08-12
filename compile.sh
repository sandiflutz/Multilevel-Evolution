#!/bin/bash

#Entrie:
#1=macro
#2=executable name

if [ -n "$1"  ]
then
	macro=$1
	if [ -n "$2"  ]
	then
		exe=$2
		make cleanobj
		make MACRO=${macro} EXEC=${exe}
		echo ""
		echo "To run the executable (recommended):"
		echo "	nice -20 ./${exe}"
		echo ""
	fi
fi
