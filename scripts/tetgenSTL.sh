#!/bin/bash

filename=$1;
#meshID=$2;
maxVolume=$2;
#./tetgen -pznnq1.414a"$maxVolume" "$filename".stl;

# Run tetgen with the following options:
# -p = uses input file
# -z = numbers all output items starting from zero
# -F = suppresses output of .face and .edge file
# -q = a boundary conforming quality tetrahedral mesh is generated
# see also http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#sec35
./tetgen -pqzFAa"$maxVolume" "$filename".stl;

# Move and rename files
rm *.smesh;
