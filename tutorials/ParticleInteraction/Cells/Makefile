USE_MPI = 1
USE_ICC = 0

#PIL_DIR = /u/home/campus/gpo/Desktop/PIL
#PIL_DIR = /home/giacomo/Dropbox/MODEL_MPI/model
PIL_DIR = /Users/giacomo/Dropbox/MODEL_MPI/model
MODEL_DIR = /Users/giacomo/Dropbox/MODEL_MPI
#EIGEN_DIR = /usr/local/include
EIGEN_DIR =/u/home/campus/gpo/local/include
#BOOST_DIR = /usr/local/include
BOOST_DIR =/u/home/campus/gpo/local/include
#METIS_DIR =/u/home/campus/gpo/local/metis
METIS_INCLUDE =/usr/local/include
METIS_LIB =/usr/local/lib

#MPI_INCLUDE = /usr/local/include
MPI_LIB = /usr/local/lib


IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(PIL_DIR)
IDIR += -I$(EIGEN_DIR)
IDIR += -I$(BOOST_DIR)
IDIR += -I$(MODEL_DIR)
#IDIR += -I$(MPI_INCLUDE)
IDIR += -I$(METIS_INCLUDE)

##########################################################
# COMPILER SETTINGS - DO NOT EDIT-
##########################################################
ifeq ($(USE_MPI), 0)
	ifeq ($(USE_ICC), 1)
		CC= icc
		# optimization level
		CFLAGS	 = -O3
		CFLAGS	+= -fast
		# enable openMP
		CFLAGS	+= -openmp
		#enable C++0x
		CFLAGS	+= -std=c++0x
	else
    		CC= g++
    		# optimization level = fast
		CFLAGS  += -O2
        	# enable openMP   
    		CFLAGS += -fopenmp
    		# enable C++0x
    		CFLAGS += -std=c++0x
	endif
else
    	ifeq ($(USE_ICC), 1)
		CC= icpc
		# optimization level
		CFLAGS += -O3
		CFLAGS += -fast
		# enable openMP
		CFLAGS += -openmp
		#enable C++0x
		CFLAGS	+= -std=c++0x
	else
    		CC= mpicxx
    		# optimization level
		CFLAGS  += -O2
        	# enable openMP   
    		CFLAGS += -fopenmp
    		# enable C++0x
    		CFLAGS += -std=c++0x        
	endif

	# Define macro
	CFLAGS += -D _MODEL_MPI_

	# Link to Metis partitioner
	LIBS += -L$(METIS_LIB) 
	LIBS +=	 -lmetis

endif

# enable common warnings… 
CFLAGS	+= -Wall
# … and more warnings       
CFLAGS += -Wextra


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

