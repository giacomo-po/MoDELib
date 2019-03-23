IDIR  = -I/usr/local/include/
IDIR += -I../../


##########################################################
# COMPILER SETTINGS
##########################################################
CC	= g++
#CC	= icc

ifeq ($(CC),g++)
	CFLAGS	 = -O3
#	CFLAGS	+= -fopenmp
	CFLAGS	+= -std=c++11
endif


ifeq ($(CC),icc)
	CFLAGS	 = -O3
#	CFLAGS	+= -openmp
	CFLAGS	+= -std=c++11
endif



##########################################################
# PROGRAM
##########################################################
PROGRAM =	test
OBJS 	= 	test


$(PROGRAM): 
	$(CC) main.cpp -o $(OBJS) $(CFLAGS) $(IDIR)  
	
clean:
		rm -f $(OBJS)

#	
#	