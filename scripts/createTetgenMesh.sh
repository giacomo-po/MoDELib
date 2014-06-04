#!/bin/bash

filename=$1;
meshID=$2;
maxVolume=$3;
#./tetgen -pznnq1.414a"$maxVolume" "$filename".stl;

# Run tetgen with the following options:
# -p = uses input file
# -z = numbers all output items starting from zero
# -F = suppresses output of .face and .edge file
# -q = a boundary conforming quality tetrahedral mesh is generated
# see also http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#sec35
./tetgen -pzFq1.414a"$maxVolume" "$filename".stl;

# Move and rename files
rm *.smesh;

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