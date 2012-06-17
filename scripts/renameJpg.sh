#!/bin/bash
#mkdir jpg
 k=1;
for i in `find jpg  -name image_\*.jpg  |sort -k1.11 -n  `
# the option "-k1.11 -n" means that the files are sorted using the 11-th column in the filename, which is the file id
 do

 filename1=$(basename "$i");
 filename2=image_"$k".jpg;
# echo renaming "$filename1" into "$filename2";
filename1WithPath=jpg/"$filename1";
filename2WithPath=jpg/"$filename2";
 echo renaming "$filename1WithPath" into "$filename2WithPath";

#  mv "jpg\$filename1" "jpg\$filename2";
 mv "$filename1WithPath" "$filename2WithPath";
  k=`expr $k + 1 `;
# sips -s format jp2 "$i" --out jpg/"${filename%.tga}.jp2"
# convert "$i" jpg/"${filename%.tga}.jpg"
done
