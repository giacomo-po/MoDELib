#!/bin/bash

filename=$1;
meshID=$2;

if [ ! -d ../T ]; then
mkdir ../T
fi
mv "$filename".1.ele ../T/T_"$meshID".txt;
sed '1d;$d' ../T/T_"$meshID".txt > ../T/T_temp.txt && mv ../T/T_temp.txt ../T/T_"$meshID".txt


if [ ! -d ../N ]; then
mkdir ../N
fi
mv "$filename".1.node ../N/N_"$meshID".txt;
sed '1d;$d' ../N/N_"$meshID".txt > ../N/N_temp.txt && mv ../N/N_temp.txt ../N/N_"$meshID".txt