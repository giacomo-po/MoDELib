#!/bin/bash

filename=mesh;
meshID=0;
maxVolume=1000000000;
../../../../scripts/createTetgenMesh.sh "$filename" "$meshID" "$maxVolume"