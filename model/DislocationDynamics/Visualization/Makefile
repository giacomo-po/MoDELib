IDIR  = -I/usr/local/include/
IDIR += -I/../../../


##########################################################
# COMPILER SETTINGS
##########################################################
#CC	= g++-mp-4.5
CC	= g++
#CC	= icc

ifeq ($(CC),g++)
	# optimization level
	CFLAGS	 = -O3
	# enable openMP
	CFLAGS	+= -fopenmp
endif


ifeq ($(CC),icc)
	# optimization level
	CFLAGS	 = -O3
	#CFLAGS	+= -fast
	# enable openMP
	CFLAGS	+= -openmp
endif

##########################################################
# OS SPECIFIC SETTINGS
##########################################################
OS= $(shell uname -s)

ifeq ($(OS),Darwin)
	LIBS =  -framework GLUT
	LIBS += -framework OpenGL 
#	LIBS += -framework Cocoa 
endif

ifeq ($(OS),Linux)
	LIBS  =  -L/usr/lib/
	LIBS +=	 -lGL
	LIBS +=	 -lGLU
	LIBS +=	 -lglut
endif

##########################################################
# PROGRAM
##########################################################
PROGRAM =	DDviewer
OBJS 	= 	DDviewer


$(PROGRAM): 
	$(CC) DDviewer.cpp -o $(OBJS) $(CFLAGS) $(IDIR)  $(LIBS) 
	
clean:
		rm -f $(OBJS)

#	
#	