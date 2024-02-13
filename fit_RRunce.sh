#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((i=1;i<=5;++i))
do
	for ((b=1;b<=1;++b))
	do 
		Rscript SIM_RRunce.R "$i" "$b"
	done
done