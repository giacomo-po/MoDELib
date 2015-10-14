MODEL_INCLUDE = ../../

program:
	g++ main.cpp -o bench -std=c++11 -fopenmp -O3  -I./ -I/usr/local/include -I $(MODEL_INCLUDE)
	
clean:
		rm -f bench
