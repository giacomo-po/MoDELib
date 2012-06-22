model_DIR = /Users/giacomo/Dropbox/Model
IDIR  = -I/usr/local/include
IDIR  = -I./
IDIR += -I$(model_DIR)

##########################################################
# COMPILER SETTINGS
##########################################################
CC	= /usr/local/bin/g++
#CC	= g++
#CC	= icc

ifeq ($(CC),/usr/local/bin/g++)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS += -ftree-vectorize

	# enable warnings
	#CFLAGS	+= -Wall
	#CFLAGS	+= -Wextra 	# ... and more warnings
	#CFLAGS	+= -pedantic # require strict ISO C and ISO C++

	# enable openMP
	CFLAGS	+= -fopenmp

	#enable C++0x
	CFLAGS	+= -std=c++0x
	#CFLAGS	+= -std=gnu++0x
endif

ifeq ($(CC),g++)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS += -ftree-vectorize

	# enable warnings
	#CFLAGS	+= -Wall
	#CFLAGS	+= -Wextra 	# ... and more warnings
	#CFLAGS	+= -pedantic # require strict ISO C and ISO C++

	# enable openMP
	CFLAGS	+= -fopenmp

	#enable C++0x
	CFLAGS	+= -std=c++0x
	#CFLAGS	+= -std=gnu++0x
# enable profiling with gprof (run with: gprof DDcode)
#CFLAGS	+=-c 
#CFLAGS	+= -pg
endif

ifeq ($(CC),icc)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS	+= -fast
	# enable common warnings
	CFLAGS	+= -Wall
	# enable openMP
	CFLAGS	+= -openmp
	#enable C++0x
	CFLAGS	+= -std=c++0x
endif

#CFLAGS += '-DCOMPLEX=std::complex<double>'

##########################################################
# OS SPECIFIC SETTINGS
##########################################################
OS= $(shell uname -s)

LIBS  =  -L/usr/local/lib/
LIBS +=	 -lgsl
LIBS +=  -lgslcblas

ifeq ($(OS),Darwin)
	
endif

ifeq ($(OS),Linux)
	LIBS +=  -lm
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
