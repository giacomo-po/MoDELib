MODEL_INCLUDE = ../../header_symlinks

program:
	g++ main.cpp -o mobility -std=c++17 -Ofast -msse4 -fopenmp -I./ -I/usr/local/include -I $(MODEL_INCLUDE)
	
clean:
		rm -f mobility
