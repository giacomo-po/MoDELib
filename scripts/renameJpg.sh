#!/bin/bash
k=0;
for i in `find . -name image_\*.jpg  |sort -k1.9 -n  `
# the option "-k1.9 -n" means that the files are sorted using the 9-th column in
# the filename, this because "./image_" are 8 characters, so the 9-th is the ID
do
    filename1=$(basename "$i");
    filename2=image_"$k".jpg;
    echo renaming "$filename1" into "$filename2";
    mv "$filename1" "$filename2";
    k=`expr $k + 1 `;
done
