MODEL_INCLUDE = ../../

program:
	g++ main.cpp -o test -std=c++11 -Ofast -I$(MODEL_INCLUDE) -I/usr/local/include -I/opt/local/include 

#-L/opt/local/lib -lCGAL -lCGAL_Core -lgmp
	
clean:
	rm -f test output.txt
