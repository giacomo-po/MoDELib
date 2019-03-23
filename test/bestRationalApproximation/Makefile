MODEL_INCLUDE = ../../

program:
	g++ main.cpp -o bra -std=c++11 -Ofast -msse4 -fopenmp -I./ -I/usr/local/include -I $(MODEL_INCLUDE)
	
clean:
		rm -f bra
