USE_ICC = 0

# MODEL LIBRARY #################################
 MODEL_INCLUDE =/Users/giacomo/Documents/MODEL
# MODEL_INCLUDE = /home/giacomo/Dropbox/Model

##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(MODEL_INCLUDE)

##########################################################
# COMPILER SETTINGS - DO NOT EDIT-
##########################################################
esc=`echo "\033"`
normalColor="${esc}[m\017"
greenColor="${esc}[0;32m"
redColor="${esc}[0;31m"
blueColor="${esc}[0;34m"

all:
	@echo ${redColor}please use \
	${greenColor}make serial \
	${redColor}or ${greenColor}make parallel \
	${redColor}or ${greenColor}make viewer \
	${redColor}or ${greenColor}make empty \
	${normalColor}


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
CFLAGS	+= -std=c++11
else
CC= g++
# optimization level = fast
CFLAGS  += -O3
# enable openMP
CFLAGS += -fopenmp
# enable C++0x
CFLAGS += -std=c++0x
endif

# enable common warnings… 
CFLAGS += -Wall
# … and more warnings       
CFLAGS += -Wextra

#disable model-assertions
#CFLAGS += -D _MODEL_NO_DEBUG_
#disable model-input-checks
#CFLAGS += -D _MODEL_NO_BAD_INPUT_CHECK_


OBJS 	= DDserial

serial: 
	@echo ${greenColor}making serial version... ${normalColor}
	$(CC) main.cpp -o $(OBJS)  $(CFLAGS)  $(IDIR)


##########################################################
# MAKE PARALLEL
##########################################################
ifeq ($(USE_ICC), 1)
mpiCC= icpc
# optimization level
mpiCFLAGS += -O3
mpiCFLAGS += -fast
# enable openMP
mpiCFLAGS += -openmp
#enable C++0x
mpiCFLAGS += -std=c++0x
else
mpiCC = mpicxx
# optimization level
mpiCFLAGS  += -O2
# enable openMP   
mpiCFLAGS += -fopenmp
# enable C++0x
mpiCFLAGS += -std=c++0x
endif

# enable common warnings… 
mpiCFLAGS += -Wall
# … and more warnings       
mpiCFLAGS += -Wextra

#disable model-assertions
mpiCFLAGS += -D _MODEL_NO_DEBUG_
#disable model-input-checks
mpiCFLAGS += -D _MODEL_NO_BAD_INPUT_CHECK_

# define _MODEL_DD_MPI_
mpiCFLAGS += -D _MODEL_DD_MPI_

mpiOBJS = DDparallel

parallel:
	@echo ${blueColor}making parallel version... ${normalColor}
	$(mpiCC) main.cpp -o $(mpiOBJS) $(mpiCFLAGS) $(IDIR) $(mpiLIBS)


##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS) $(mpiOBJS)


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

##########################################################
# MAKE EMPTY
##########################################################
#mkdirV=$(shell test -d V || mkdir V)
#mkdirE=$(shell test -d E || mkdir E)
Gexists=$(shell test -d G && echo 1)
Cexists=$(shell test -d C && echo 1)
Pexists=$(shell test -d P && echo 1)
Nexists=$(shell test -d N && echo 1)
Texists=$(shell test -d T && echo 1)
TGAexists=$(shell test -d tga && echo 1)
JPGexists=$(shell test -d jpg && echo 1)

empty:
	@find V/ -name V_\*.txt ! -name V_0.txt -exec rm {} \;
	@find V/ -name V_\*.bin ! -name V_0.bin -exec rm {} \;
	@find E/ -name E_\*.txt ! -name E_0.txt -exec rm {} \;
	@find E/ -name E_\*.bin ! -name E_0.bin -exec rm {} \;
ifeq ($(Gexists),1)
	@echo ${blueColor}emptying folder G ${normalColor}
	@find G/ -name G_\*.txt -exec rm {} \;
	@find G/ -name G_\*.bin -exec rm {} \;
endif
ifeq ($(Cexists),1)
	@echo ${blueColor}emptying folder C ${normalColor}
	@find C/ -name C_\*.txt -exec rm {} \;
	@find C/ -name C_\*.bin -exec rm {} \;
endif
ifeq ($(Pexists),1)
	@echo ${blueColor}emptying folder P ${normalColor}
	@find P/ -name P_\*.txt -exec rm {} \;
	@find P/ -name P_\*.bin -exec rm {} \;
endif	
ifeq ($(Nexists),1)
	@echo ${blueColor}emptying folder N ${normalColor}
	@find N/ -name N_\*.txt -exec rm {} \;
	@find N/ -name N_\*.bin -exec rm {} \;
endif	
ifeq ($(Texists),1)
	@echo ${blueColor}emptying folder T ${normalColor}
	@find T/ -name T_\*.txt -exec rm {} \;
	@find T/ -name T_\*.bin -exec rm {} \;
endif	
ifeq ($(TGAexists),1)
	@echo ${blueColor}emptying folder tga ${normalColor}
	@find tga/ -name image_\*.tga -exec rm {} \;
endif	
ifeq ($(JPGexists),1)
	@echo ${blueColor}emptying folder jpg ${normalColor}
	@find jpg/ -name image_\*.jpg -exec rm {} \;
endif
