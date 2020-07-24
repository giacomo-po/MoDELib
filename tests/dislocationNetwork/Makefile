MODEL_INCLUDE = ../../header_symlinks

program:
	g++ main.cpp -o loop -std=c++17 -fopenmp -O3 -I$(MODEL_INCLUDE) -I. -I/usr/local/include -Wall
	
clean:
		rm -f loop