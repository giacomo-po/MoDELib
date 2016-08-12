MODEL_INCLUDE = ../../
EIGEN_INCLUDE = /Users/Giacomo/Documents/eigen

##############################
IDIR += -I ./
IDIR += -I $(MODEL_INCLUDE)
IDIR += -I $(EIGEN_INCLUDE)

##############################
CFLAGS += -std=c++11  
CFLAGS += -O3 
CFLAGS += -fopenmp
CFLAGS += -Wall
CFLAGS += -Wextra

##############################
program:
	g++ main.cpp -o ppi $(CFLAGS) $(IDIR) 
	
clean:
		rm -f ppi
		