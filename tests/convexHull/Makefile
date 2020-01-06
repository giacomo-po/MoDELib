MODEL_INCLUDE = ../../header_symlinks

program:
	g++ main.cpp -o hull -std=c++17 -Ofast -msse4 -I./ -I/usr/local/include -I $(MODEL_INCLUDE)
	
clean:
		rm -f hull
