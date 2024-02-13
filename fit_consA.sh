#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((i=1;i<=8;++i))
do
	for ((b=1;b<=1;++b))
	do 
		Rscript SIM_ConsA.R "$i" "$b"
	done
done