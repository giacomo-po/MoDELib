#!/bin/bash

filename=cylinder_2k_12k;
maxVolume=5000000;
./tetgen -pznnq1.414a"$maxVolume" "$filename".stl;

rm *.smesh;
mv "$filename".1.ele mesh.ele;
mv "$filename".1.node mesh.node;
mv "$filename".1.face mesh.face;
mv "$filename".1.neigh mesh.neigh;
../writeBCs
