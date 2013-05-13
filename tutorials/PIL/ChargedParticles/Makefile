PIL_DIR = /Users/giacomo/Dropbox/PIL
EIGEN_DIR = /usr/local/include
BOOST_DIR = /usr/local/include
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(PIL_DIR)
IDIR += -I$(EIGEN_DIR)
IDIR += -I$(BOOST_DIR)

##########################################################
# COMPILER SETTINGS
##########################################################
#CC	= g++
CC	= /opt/local/bin/g++-mp-4.7
#CC	= icc
#CC = clang++


ifeq ($(CC), /opt/local/bin/g++-mp-4.7)

    # optimization level = fast
	CFLAGS  = -O2
    # auto-vectorization
    #CFLAGS += -ftree-vectorize
    # enable openMP
    #CFLAGS += -fopenmp
    # enable C++0x
    #CFLAGS += -std=c++0x
    CFLAGS += -std=c++11

    # enable warnings
	CFLAGS += -Wall
    # and more warnings
    CFLAGS += -Wextra

endif



ifeq ($(CC),icc)
	# optimization level
	CFLAGS	 = -O2
	CFLAGS	+= -fast

	# enable openMP
	CFLAGS	+= -openmp
	#enable C++0x
	CFLAGS	+= -std=c++0x

    # enable common warnings
	CFLAGS	+= -Wall

endif

##########################################################
# OS SPECIFIC SETTINGS
##########################################################
OS= $(shell uname -s)

#LIBS  =  -L/usr/local/lib/
#LIBS +=	 -lgsl
#LIBS +=  -lgslcblas

ifeq ($(OS),Darwin)
    #LIBS += -Wl,-no_pie
    #LIBS +=	-Wl,-allow_heap_execute
endif

ifeq ($(OS),Linux)
#	LIBS +=  -lm
endif



##########################################################
# MAKE PROGRAM
##########################################################
PROGRAM = particles
OBJS 	= particles

$(PROGRAM):
	$(CC) main.cpp -o $(OBJS)  $(CFLAGS)  $(IDIR) $(LIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)

