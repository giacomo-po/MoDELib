#!/bin/bash

filename=$1;
meshID=$2;

#if [ ! -d ../T ]; then
#mkdir ../T
#fi
mv "$filename".1.ele T_"$meshID".txt;
sed '1d;$d' T_"$meshID".txt > T_temp.txt && mv T_temp.txt T_"$meshID".txt


#if [ ! -d ../N ]; then
#mkdir ../N
#fi
mv "$filename".1.node N_"$meshID".txt;
sed '1d;$d' N_"$meshID".txt > N_temp.txt && mv N_temp.txt N_"$meshID".txt
