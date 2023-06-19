#!/bin/bash
k=0;
for i in `find ../evl -name img_\*.png  |sort -k1.12 -n  `
# the option "-k1.x -n" means that the files are sorted using the x-th column in
# the filename. Note "../evl/img_" are 11 characters, so the 12-th is the ID
do
#    filename1=$(basename "$i");
    filename2=img_"$k".png;
    echo renaming $i into "$filename2";
    mv $i "$filename2";
    k=`expr $k + 1 `;
done
