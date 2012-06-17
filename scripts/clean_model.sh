# execute with:
# sh clean_model.sh 

# create the ./trash directory
mkdir -p ./trash;

# find and permanently remove all files and folders ".AppleDouble" 
find . -name ".AppleDouble" -exec rm -r -i {} \;

# find and move to ./trash all files ".DS_Store"
find . -name ".DS_Store"    -exec mv {} ./trash \;

# find and move to ./trash all files "_2e_*" or ":2e_*"
find . -name _2e_\* -exec mv {} ./trash \;
find . -name :2e_\* -exec mv {} ./trash \;

# Check what survived
find . -name ".AppleDouble" -print;
find . -name ".DS_Store" -print;
find . -name _2e_\* -print;
find . -name :2e_\* -print;

