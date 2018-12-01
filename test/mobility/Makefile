MODEL_INCLUDE = ../../header_symlinks

program:
	g++ main.cpp -o bccMobility -std=c++14 -Ofast -msse4 -fopenmp -I./ -I/usr/local/include -I $(MODEL_INCLUDE)
	
clean:
		rm -f bccMobility
