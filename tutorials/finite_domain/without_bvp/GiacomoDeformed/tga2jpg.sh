#!/bin/bash
mkdir jpg

for i in `find tga -name image_\*.tga `
 do
 filename=$(basename "$i");
 echo converting "$filename";
# sips -s format jp2 "$i" --out jpg/"${filename%.tga}.jp2"
 convert "$i" jpg/"${filename%.tga}.jpg"
done
