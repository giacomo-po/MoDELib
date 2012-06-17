#mmdl_DIR = /home/giacomo/Documents/UCLA/Research/C++/Mmdl
#mmdl_DIR = /Users/giacomo/Documents/UCLA/Research/C++/Mmdl
#mmdl_DIR = /home/giacomo/Desktop/BackUp/Documents/UCLA/Research/C++/Mmdl
mmdl_DIR = /home/giacomo/Dropbox/Mmdl
#mmdl_DIR = /home/helios/a/mohamed3/Research/Indentation/Mmdl/

IDIR  = -I/usr/local/include
IDIR +=  -I./
IDIR += -I$(mmdl_DIR)


#IDIR += -I/Users/giacomo/Documents/UCLA/Research/C++/Mmdl/
#IDIR += -I/home/giacomo/Documents/UCLA/Research/C++/Mmdl/
#IDIR += -I/home/helios/a/mohamed3/Research/Libraries/include/
#IDIR += -I/home/mamdouh/Work/Indentation/Mmdl/


##########################################################
# COMPILER SETTINGS
##########################################################
#CC	= g++-mp-4.6
#CC	= g++-mp-4.5
#CC	= g++-mp-4.4

OS= $(shell uname -s)
ifeq ($(OS),Darwin)
	CC	= /usr/local/bin/g++
endif

ifeq ($(OS),Linux)
CC	= g++
endif

#CC	= icc

ifeq ($(CC),/usr/local/bin/g++)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS += -ftree-vectorize

	# enable common warnings
	#CFLAGS	+= -Wall
	# ... and more warnings
	#CFLAGS	+= -Wextra
	# require strict ISO C and ISO C++
	#CFLAGS	+= -pedantic
	# enable openMP
	CFLAGS	+= -fopenmp
	#enable C++0x
	#CFLAGS	+= -std=c++0x
	CFLAGS	+= -std=gnu++0x
endif

ifeq ($(CC),g++)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS += -ftree-vectorize

	# enable common warnings
	#CFLAGS	+= -Wall
	# ... and more warnings
	#CFLAGS	+= -Wextra
	# require strict ISO C and ISO C++
	#CFLAGS	+= -pedantic
	# enable openMP
	CFLAGS	+= -fopenmp
	#enable C++0x
	#CFLAGS	+= -std=c++0x
	CFLAGS	+= -std=gnu++0x
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
#OS= $(shell uname -s)

LIBS  =  -L/usr/local/lib/
LIBS +=	 -lgsl
LIBS +=  -lgslcblas
#LIBS +=	 -lsparse
#LIBS +=  -lspblas
#LIBS +=	 -lmv
#LIBS +=  -ltet

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
# MAKE EMPTY
##########################################################
empty:
	sh emptyOutputs

##########################################################
# MAKE VIEWER
##########################################################	
viewer: 
	g++ -O3 $(mmdl_DIR)/mmdl/Dislocations/Visualization/DDviewer.cpp -o DDviewer  $(IDIR)  $(VIEWERLIBS) 
