#!/bin/bash

filename=mesh;
meshID=0;
maxVolume=10000000;
../../../../scripts/createTetgenMesh.sh "$filename" "$meshID" "$maxVolume"
