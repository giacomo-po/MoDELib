 model_DIR = /Users/giacomo/Dropbox/Model
# model_DIR = /home/giacomo/Dropbox/Model
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(model_DIR)

##########################################################
# COMPILER SETTINGS
##########################################################
#CC	= /usr/local/bin/g++
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
    CFLAGS += -fopenmp
    # enable C++0x
    CFLAGS += -std=c++0x
    
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

ifeq ($(OS),Darwin)
	VIEWERLIBS =  -framework GLUT
	VIEWERLIBS += -framework OpenGL 
endif

ifeq ($(OS),Linux)
	VIEWERLIBS  =  -L/usr/lib/
	VIEWERLIBS +=	 -lGL
	VIEWERLIBS +=	 -lGLU
	VIEWERLIBS +=	 -lglut
endif	


##########################################################
# MAKE PROGRAM
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
	$(CC) $(model_DIR)/model/Dislocations/Visualization/DDviewer.cpp -o DDviewer -O3 $(IDIR)  $(VIEWERLIBS) 

##########################################################
# MAKE EMPTY
##########################################################
empty:
	sh emptyOutputs
