IDIR  = -I/usr/local/include/
IDIR += -I../../../


##########################################################
# COMPILER SETTINGS
##########################################################
CC	= g++
#CC	= icc

ifeq ($(CC),g++)
	CFLAGS	 = -O3
	CFLAGS	+= -fopenmp
	CFLAGS	+= -std=c++14
	CFLAGS  += -Wall
endif


ifeq ($(CC),icc)
	CFLAGS	 = -O3
#	CFLAGS	+= -openmp
	CFLAGS	+= -std=c++14
	CFLAGS  += -Wall
endif



##########################################################
# PROGRAM
##########################################################
microstructureGenerator: 
	$(CC) microstructureGenerator.cpp -o microstructureGenerator $(CFLAGS) $(IDIR)  
	
clean:
		rm -f microstructureGenerator
