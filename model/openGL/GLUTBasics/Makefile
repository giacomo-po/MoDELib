 model_DIR = /Users/giacomo/Dropbox/Model
# model_DIR = /home/giacomo/Dropbox/Model
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(model_DIR)

##########################################################
# COMPILER SETTINGS
##########################################################
CC	= g++
CFLAGS  = -O2
	CFLAGS += -Wall
    CFLAGS += -Wextra


##########################################################
# OS SPECIFIC SETTINGS
##########################################################
OS= $(shell uname -s)

ifeq ($(OS),Darwin)
	LIBS =  -framework GLUT
	LIBS += -framework OpenGL 
endif

ifeq ($(OS),Linux)
	LIBS  =  -L/usr/lib/
	LIBS +=	 -lGL
	LIBS +=	 -lGLU
	LIBS +=	 -lglut
endif	


##########################################################
# MAKE PROGRAM
##########################################################
PROGRAM = DDviewer
OBJS 	= DDviewer

$(PROGRAM): 
	$(CC) DDviewer.cpp trackball.c SurfaceGeometry.c -o DDviewer -O3 $(CFLAGS) $(IDIR)  $(LIBS) 
##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)





