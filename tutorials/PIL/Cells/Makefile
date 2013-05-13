PIL_DIR = /u/home/campus/gpo/Desktop/PIL
#PIL_DIR = /home/giacomo/Dropbox/PIL
#EIGEN_DIR = /usr/local/include
EIGEN_DIR =/u/home/campus/gpo/local/include
#BOOST_DIR = /usr/local/include
BOOST_DIR =/u/home/campus/gpo/local/include
METIS_DIR =/u/home/campus/gpo/local/metis
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(PIL_DIR)
IDIR += -I$(EIGEN_DIR)
IDIR += -I$(BOOST_DIR)
IDIR += -I$(METIS_DIR)/include

##########################################################
# COMPILER SETTINGS
##########################################################
CC	= mpicxx
#CC	= g++
#CC	= /opt/local/bin/g++-mp-4.7
#CC	= icc
#CC = clang++


ifeq ($(CC), g++)

    # optimization level = fast
	CFLAGS  = -O2
    # auto-vectorization
    #CFLAGS += -ftree-vectorize
    # enable openMP
    #CFLAGS += -fopenmp
    # enable C++0x
    CFLAGS += -std=c++0x
    #CFLAGS += -std=c++11

    # enable warnings
	CFLAGS += -Wall
    # and more warnings
    CFLAGS += -Wextra

endif



ifeq ($(CC),mpicxx)
	# optimization level
	CFLAGS	 = -O2
	#CFLAGS	+= -fast

	# enable openMP
	#CFLAGS	+= -fopenmp
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

LIBS += -L$(METIS_DIR)/lib -lmetis

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
	$(CC) main.cpp -o $(OBJS) $(CFLAGS)  $(IDIR) $(LIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)

