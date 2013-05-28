USE_MPI = 1
USE_ICC = 0


# MODEL LIBRARY #################################
 MODEL_INCLUDE = /Users/giacomo/Dropbox/Model_MPI
# MODEL_INCLUDE = /home/giacomo/Dropbox/Model

# MPI LIBRARY ###################################
#MPI_INCLUDE = /usr/local/include
MPI_LIB = /usr/local/lib

# METIS LIBRARY #################################
METIS_INCLUDE =/usr/local/include
METIS_LIB =/usr/local/lib

# METIS LIBRARY #################################
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
	CFLAGS += -D _MODEL_DD_MPI_

	# Link to Metis partitioner
	LIBS += -L$(METIS_LIB) 
	LIBS +=	 -lmetis

endif

# enable common warnings… 
CFLAGS	+= -Wall
# … and more warnings       
CFLAGS += -Wextra


##########################################################
# OS-SPECIFIC SETTINGS
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

##########################################################
# MAKE
##########################################################
PROGRAM = DDcode
OBJS 	= DDcode

$(PROGRAM): 
	$(CC) main.cpp -o $(OBJS)  $(CFLAGS)  $(IDIR) $(LIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)

##########################################################
# MAKE VIEWER
##########################################################	

viewer: 
	$(CC) $(MODEL_INCLUDE)/model/DislocationDynamics/Visualization/DDviewer.cpp -o DDviewer -O3 $(IDIR)  $(VIEWERLIBS) 

