#!/bin/bash

filename=cube_4000;
maxVolume=500000;
./tetgen -pznnq1.414a"$maxVolume" "$filename".stl;

rm *.smesh;
mv "$filename".1.ele mesh.ele;
mv "$filename".1.node mesh.node;
mv "$filename".1.face mesh.face;
mv "$filename".1.neigh mesh.neigh;
#writeBCs
