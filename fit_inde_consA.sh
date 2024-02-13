#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((i=1;i<=4;++i))
do
	for ((b=1;b<=1;++b))
	do 
		Rscript SIM_Inde_ConsA.R "$i" "$b"
	done
done