#USE_MPI = 0
USE_ICC = 0


# MODEL LIBRARY #################################
 MODEL_INCLUDE = /Users/giacomo/Dropbox/Model_MPI
# MODEL_INCLUDE = /home/giacomo/Dropbox/Model

# MPI LIBRARY (PARALLEL ONLY) ###################
#MPI_INCLUDE = /usr/local/include
MPI_LIB = /usr/local/lib

# METIS LIBRARY (PARALLEL ONLY) #################
METIS_INCLUDE =/usr/local/include
METIS_LIB =/usr/local/lib

##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(MODEL_INCLUDE)
#IDIR += -I$(MPI_INCLUDE)
IDIR += -I$(METIS_INCLUDE)

##########################################################
# COMPILER SETTINGS - DO NOT EDIT-
##########################################################
all:
	@echo please use "make serial" or "make parallel"


##########################################################
# MAKE SERIAL
##########################################################
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

# enable common warnings… 
CFLAGS	+= -Wall
# … and more warnings       
CFLAGS += -Wextra

OBJS 	= DDserial

serial: 
	$(CC) main.cpp -o $(OBJS)  $(CFLAGS)  $(IDIR)


##########################################################
# MAKE PARALLEL
##########################################################
ifeq ($(USE_ICC), 1)
		MPICC= icpc
		# optimization level
		CFLAGS += -O3
		CFLAGS += -fast
		# enable openMP
		CFLAGS += -openmp
		#enable C++0x
		CFLAGS	+= -std=c++0x
else
    	MPICC = mpicxx
    	# optimization level
		CFLAGS  += -O2
        # enable openMP   
    	CFLAGS += -fopenmp
    	# enable C++0x
    	CFLAGS += -std=c++0x        
endif

# enable common warnings… 
CFLAGS	+= -Wall
# … and more warnings       
CFLAGS += -Wextra

# Link to Metis partitioner
LIBS += -L$(METIS_LIB) 
LIBS +=	 -lmetis

MPIOBJS = DDparallel

parallel:
	$(MPICC) main.cpp -o $(MPIOBJS) $(CFLAGS) -D _MODEL_DD_MPI_ $(IDIR) $(LIBS)


##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS) $(MPIOBJS)


##########################################################
# MAKE VIEWER
##########################################################	
OS= $(shell uname -s)

ifeq ($(OS),Darwin)
	VIEWERLIBS += -framework GLUT
	VIEWERLIBS += -framework OpenGL 
endif

ifeq ($(OS),Linux)
	VIEWERLIBS += -L/usr/lib/
	VIEWERLIBS += -lGL
	VIEWERLIBS += -lGLU
	VIEWERLIBS += -lglut
endif	

viewer: 
	$(CC) $(MODEL_INCLUDE)/model/DislocationDynamics/Visualization/DDviewer.cpp -o DDviewer -O3 $(IDIR)  $(VIEWERLIBS) 
