USE_ICC = 0
usePARDISO = 0

##########################################################
# EDIT HERE
##########################################################
MODEL_INCLUDE =../../../
MKL_PATH=/opt/intel/composer_xe_2013.3.171/mkl

##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(MODEL_INCLUDE)

##########################################################
# COMPILER SETTINGS - DO NOT EDIT-
##########################################################
ifeq ($(USE_ICC), 1)
CC= icpc
mpiCC= icpc
# optimization level
CFLAGS	 = -O3
CFLAGS	+= -fast
# enable openMP
CFLAGS	+= -openmp 
#enable C++0x
CFLAGS	+= -std=c++0x
else
CC= g++
mpiCC = mpicxx
# optimization level = fast
#CFLAGS  += -O3
CFLAGS  += -Ofast
CFLAGS  += -msse4
#CFLAGS  += -mavx
#CFLAGS  += -march=native
# enable openMP
CFLAGS += -fopenmp
# enable C++0x
CFLAGS += -std=c++0x
endif
# enable common warnings… 
CFLAGS += -Wall
# … and more warnings       
CFLAGS += -Wextra
#CFLAGS += -ftime-report

#disable model-assertions
CFLAGS += -D _MODEL_NO_DEBUG_
#disable model-input-checks
CFLAGS += -D _MODEL_NO_BAD_INPUT_CHECK_

# Settings for PARDISO solver
ifeq ($(usePARDISO), 1)
CFLAGS += -D _MODEL_PARDISO_SOLVER_
IDIR += -I $(MKL_PATH)/include
LIBS += -L $(MKL_PATH)/lib
LIBS += -lmkl_intel_lp64 
LIBS += -lmkl_intel_thread 
LIBS += -lmkl_core
endif

##########################################################
# MAKE ALL
##########################################################

esc=`echo "\033"`
normalColor="${esc}[m\017"
greenColor="${esc}[0;32m"
redColor="${esc}[0;31m"
blueColor="${esc}[0;34m"

all:
	@echo ${redColor}please use \
	${greenColor}make DDomp \
	${redColor}or ${greenColor}make DDmpi \
	${redColor}or ${greenColor}make DDviewer \
	${redColor}or ${greenColor}make empty \
	${normalColor}


##########################################################
# MAKE DDOMP
##########################################################
OBJS 	= DDomp

$(OBJS): 
	@echo ${greenColor}making omp version... ${normalColor}
	$(CC) main.cpp -o  $(OBJS)  $(CFLAGS) $(LIBS) $(IDIR)


##########################################################
# MAKE DDMPI
##########################################################
#mpiCFLAGS += -D _MODEL_DD_MPI_
mpiCFLAGS += -D _MODEL_MPI_

mpiOBJS = DDmpi

$(mpiOBJS):
	@echo ${blueColor}making mpi version... ${normalColor}
	$(mpiCC) main.cpp -o $(mpiOBJS) $(CFLAGS) $(mpiCFLAGS) $(IDIR) $(LIBS) $(mpiLIBS)


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

DDviewer: 
	$(CC) $(MODEL_INCLUDE)/model/DislocationDynamics/Visualization/DDviewer.cpp -o DDviewer -O3 -std=c++11 $(IDIR)  $(VIEWERLIBS) 

##########################################################
# MAKE EMPTY
##########################################################
#mkdirV=$(shell test -d V || mkdir V)
#mkdirE=$(shell test -d E || mkdir E)
Fexists=$(shell test -d F && echo 1)
Gexists=$(shell test -d G && echo 1)
Cexists=$(shell test -d C && echo 1)
Pexists=$(shell test -d P && echo 1)
Dexists=$(shell test -d D && echo 1)
Sexists=$(shell test -d S && echo 1)
Uexists=$(shell test -d U && echo 1)
TGAexists=$(shell test -d tga && echo 1)
JPGexists=$(shell test -d jpg && echo 1)

empty:
	@find V/ -name V_\*.txt ! -name V_0.txt -exec rm {} \;
	@find V/ -name V_\*.bin ! -name V_0.bin -exec rm {} \;
	@find E/ -name E_\*.txt ! -name E_0.txt -exec rm {} \;
	@find E/ -name E_\*.bin ! -name E_0.bin -exec rm {} \;
ifeq ($(Fexists),1)
	@echo ${blueColor}emptying folder F ${normalColor}
	@find F/ -name F_\*.txt -exec rm {} \;
	@find F/ -name F_\*.bin -exec rm {} \;
endif
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
ifeq ($(Dexists),1)
	@echo ${blueColor}emptying folder D ${normalColor}
	@find D/ -name D_\*.txt -exec rm {} \;
	@find D/ -name D_\*.bin -exec rm {} \;
endif	
ifeq ($(Sexists),1)
	@echo ${blueColor}emptying folder S ${normalColor}
	@find S/ -name S_\*.txt -exec rm {} \;
	@find S/ -name S_\*.bin -exec rm {} \;
endif	
ifeq ($(Uexists),1)
	@echo ${blueColor}emptying folder U ${normalColor}
	@find U/ -name U_\*.txt -exec rm {} \;
	@find U/ -name U_\*.bin -exec rm {} \;
endif	
ifeq ($(TGAexists),1)
	@echo ${blueColor}emptying folder tga ${normalColor}
	@find tga/ -name image_\*.tga -exec rm {} \;
endif	
ifeq ($(JPGexists),1)
	@echo ${blueColor}emptying folder jpg ${normalColor}
	@find jpg/ -name image_\*.jpg -exec rm {} \;
endif
