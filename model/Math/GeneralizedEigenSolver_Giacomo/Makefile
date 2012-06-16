IDIR  = -I/usr/local/include/
#IDIR += -I/Users/giacomo/Documents/UCLA/Research/C++/Model/
#IDIR += -I/home/giacomo/Documents/UCLA/Research/C++/Model/
IDIR += -I/home/giacomo/Dropbox/Model/


##########################################################
# COMPILER SETTINGS
##########################################################
#CC	= g++-mp-4.5
CC	= g++
#CC	= icc

ifeq ($(CC),g++)
	# optimization level
	CFLAGS	 = -O3
	# enable common warnings
	#CFLAGS	+= -Wall
	# ... and more warnings
	#CFLAGS	+= -Wextra
	# require strict ISO C and ISO C++
	#CFLAGS	+= -pedantic
	# enable openMP
	CFLAGS	+= -fopenmp
	CFLAGS	+= -std=gnu++0x
endif

ifeq ($(CC),g++-mp-4.5)
	# optimization level
	CFLAGS	 = -O3
	# enable common warnings
	CFLAGS	+= -Wall
	# ... and more warnings
	#CFLAGS	+= -Wextra
	# require strict ISO C and ISO C++
	#CFLAGS	+= -pedantic
	# enable openMP
	CFLAGS	+= -fopenmp
	#enable C++0x
	CFLAGS	+= -std=gnu++0x
endif

ifeq ($(CC),icc)
	# optimization level
	CFLAGS	 = -O3
	#CFLAGS	+= -fast
	# enable openMP
	CFLAGS	+= -openmp
endif

#CFLAGS += '-DCOMPLEX=std::complex<double>'

##########################################################
# OS SPECIFIC SETTINGS
##########################################################
#OS= $(shell uname -s)


LIBS  =  -L/usr/local/lib/
LIBS +=	 -lgsl
LIBS +=  -lgslcblas

	


##########################################################
# PROGRAM
##########################################################
PROGRAM = ges
OBJS 	= ges

$(PROGRAM): 
	$(CC) main.cpp -o $(OBJS)  $(CFLAGS)  $(IDIR) $(LIBS)
	
clean:
		rm -f $(OBJS)
