# requires imagemagick library
# put this file in your tga folder.
#!/bin/bash
mkdir ../jpg

for i in `find . -name image_\*.tga `
 do
 filename=$(basename "$i");
 echo converting "$filename";
 convert "$i" ../jpg/"${filename%.tga}.jpg"
done
