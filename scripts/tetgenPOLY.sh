#!/bin/bash
# Run tetgen with the following options:
# -p = uses input file
# -q = a boundary conforming quality tetrahedral mesh is generated
# -z = numbers all output items starting from zero
# -F = suppresses output of .face and .edge file
# -a = applies a maximum tetrahedron volume constraint
# see also http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html#sec35
filename=$1;
#./tetgen -pqzFa "$filename".poly
./tetgen -pqzAFa "$filename".poly
