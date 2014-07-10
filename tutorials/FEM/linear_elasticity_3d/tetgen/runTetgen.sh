#!/bin/bash

filename=mesh;
meshID=6;
maxVolume=0.000001;
../../../../scripts/createTetgenMesh.sh "$filename" "$meshID" "$maxVolume"