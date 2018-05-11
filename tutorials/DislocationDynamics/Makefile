##########################################################
# EDIT HERE
##########################################################
USE_ICC = 0
usePARDISO = 0
MODEL_INCLUDE =../../..
EIGEN_INCLUDE =/usr/local/include
MKL_INCLUDE=/opt/intel/mkl/include
MKL_LIB=/opt/intel/mkl/lib

##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR += -I./
IDIR += -I$(MODEL_INCLUDE)
IDIR += -I$(EIGEN_INCLUDE)

##########################################################
# COMPILER SETTINGS - DO NOT EDIT-
##########################################################
OS= $(shell uname -s)

ifeq ($(USE_ICC), 1)
	CC= icpc
	mpiCC= icpc
	CFLAGS	 = -O3
	CFLAGS	+= -xhost
	#CFLAGS	+= -fast
	CFLAGS	+= -openmp 
	CFLAGS	+= -std=c++14
else
	CC= g++
	mpiCC = mpicxx
	CFLAGS += -Ofast
	CFLAGS += -fopenmp
	CFLAGS += -std=c++14

	ifeq ($(OS),Darwin)
		CFLAGS += -msse4
	endif

	ifeq ($(OS),Linux)
		CFLAGS += -msse4
#s		CFLAGS += -march=native
	endif

endif


# enable warnings
CFLAGS += -Wall 
#CFLAGS += -Wextra

############################################
# Compiler flags for MKL PARDISO
ifeq ($(usePARDISO), 1)
	CFLAGS += -D _MODEL_PARDISO_SOLVER_
	IDIR += -I $(MKL_INCLUDE)
	LIBS += -L $(MKL_LIB)

	ifeq ($(OS),Darwin)
		LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
	endif

	ifeq ($(OS),Linux)
		ifeq ($(USE_ICC), 1)
			LIBS += -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -Wl,--end-group 
		else
			LIBS += -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -Wl,--end-group 
		endif
	endif

endif

############################################
# Compiler flags for DDviewer (OpenGL flags)
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
# MAKE ALL
##########################################################
all:
	@echo please use \
	make DDomp \
	or make DDmpi \
	or make DDviewer \
	or make dipoleGenerator \
	or make frGenerator \
	or make DDconverter \
	or make empty \

##########################################################
# MAKE DDOMP
##########################################################
OBJS 	= DDomp

$(OBJS): 
	@echo making DDomp
	$(CC) main.cpp -o  $(OBJS)  $(CFLAGS) $(LIBS) $(IDIR)


##########################################################
# MAKE DDMPI
##########################################################
#mpiCFLAGS += -D _MODEL_DD_MPI_
mpiCFLAGS += -D _MODEL_MPI_

mpiOBJS = DDmpi

$(mpiOBJS):
	@echo making DDmpi
	$(mpiCC) main.cpp -o $(mpiOBJS) $(CFLAGS) $(mpiCFLAGS) $(IDIR) $(LIBS) $(mpiLIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
	rm -f $(OBJS) $(mpiOBJS)

##########################################################
# MAKE VIEWER
##########################################################
DDviewer: 
	$(CC) $(MODEL_INCLUDE)/model/DislocationDynamics/Visualization/DDviewer.cpp -o DDviewer -O3 -std=c++14 $(IDIR)  $(VIEWERLIBS) 

##########################################################
# MAKE DIPOLEGENERATOR
##########################################################
dipoleOBJS = dipoleGenerator

$(dipoleOBJS):
	@echo making dipoleGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/dipoleGenerator.cpp -o dipoleGenerator -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE DIPOLEGENERATOR
##########################################################
piOBJS = piGenerator

$(piOBJS):
	@echo making piGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/piGenerator.cpp -o piGenerator -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE FRANKREADGENERATOR
##########################################################
frOBJS = frGenerator

$(frOBJS):
	@echo making frGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/frGenerator.cpp -o frGenerator -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE straightGenerator
##########################################################
straightOBJS = straightGenerator

$(straightOBJS):
	@echo making straightGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/straightGenerator.cpp -o straightGenerator -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE junctionGenerator
##########################################################
junctionOBJS = junctionGenerator

$(junctionOBJS):
	@echo making junctionGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/junctionGenerator.cpp -o junctionGenerator -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE DDCONVERTER
##########################################################
converterOBJS = DDconverter

$(converterOBJS):
	@echo making DDconverter
	$(CC) $(MODEL_INCLUDE)/model/DislocationDynamics/DDconverter/DDconverter.cpp -o DDconverter -O3 -std=c++14 $(IDIR)

##########################################################
# MAKE EMPTY
##########################################################
#mkdirV=$(shell test -d V || mkdir V)
#mkdirE=$(shell test -d E || mkdir E)
Fexists=$(shell test -d F && echo 1)
Gexists=$(shell test -d G && echo 1)
Cexists=$(shell test -d C && echo 1)
Kexists=$(shell test -d K && echo 1)
Lexists=$(shell test -d L && echo 1)
Pexists=$(shell test -d P && echo 1)
Dexists=$(shell test -d D && echo 1)
Sexists=$(shell test -d S && echo 1)
Uexists=$(shell test -d U && echo 1)
Qexists=$(shell test -d Q && echo 1)
Zexists=$(shell test -d Z && echo 1)
TGAexists=$(shell test -d tga && echo 1)
JPGexists=$(shell test -d jpg && echo 1)

empty:
	@echo emptying folder V
	@find V/ -name V_\*.txt ! -name V_0.txt -exec rm {} +;
	@find V/ -name V_\*.bin -exec rm {} +;
	@echo emptying folder E
	@find E/ -name E_\*.txt ! -name E_0.txt -exec rm {} +;
	@find E/ -name E_\*.bin -exec rm {} +;
	@echo emptying folder L
	@find L/ -name L_\*.txt ! -name L_0.txt -exec rm {} +;
	@find L/ -name L_\*.bin -exec rm {} +;
ifeq ($(Fexists),1)
	@echo emptying folder F
	@find F/ -name F_\*.txt -exec rm {} +;
	@find F/ -name F_\*.bin -exec rm {} +;
endif
ifeq ($(Gexists),1)
	@echo emptying folder G
	@find G/ -name G_\*.txt -exec rm {} +;
	@find G/ -name G_\*.bin -exec rm {} +;
endif
ifeq ($(Cexists),1)
	@echo emptying folder C
	@find C/ -name C_\*.txt -exec rm {} +;
	@find C/ -name C_\*.bin -exec rm {} +;
endif
ifeq ($(Kexists),1)
	@echo emptying folder K
	@find K/ -name K_\*.txt -exec rm {} +;
	@find K/ -name K_\*.bin -exec rm {} +;
endif
ifeq ($(Pexists),1)
	@echo emptying folder P
	@find P/ -name P_\*.txt -exec rm {} +;
	@find P/ -name P_\*.bin -exec rm {} +;
endif	
ifeq ($(Dexists),1)
	@echo emptying folder D
	@find D/ -name D_\*.txt -exec rm {} +;
	@find D/ -name D_\*.bin -exec rm {} +;
endif	
ifeq ($(Sexists),1)
	@echo emptying folder S
	@find S/ -name S_\*.txt -exec rm {} +;
	@find S/ -name S_\*.bin -exec rm {} +;
endif	
ifeq ($(Uexists),1)
	@echo emptying folder U
	@find U/ -name U_\*.txt -exec rm {} +;
	@find U/ -name U_\*.bin -exec rm {} +;
endif
ifeq ($(Zexists),1)
	@echo emptying folder Z
	@find Z/ -name Z_\*.txt -exec rm {} +;
	@find Z/ -name Z_\*.bin -exec rm {} +;
endif
ifeq ($(Qexists),1)
	@echo emptying folder Q
	@find Q/ -name Q_\*.txt -exec rm {} +;
	@find Q/ -name Q_\*.bin -exec rm {} +;
endif	
ifeq ($(TGAexists),1)
	@echo emptying folder tga
	@find tga/ -name image_\*.tga -exec rm {} +;
endif	
ifeq ($(JPGexists),1)
	@echo emptying folder jpg
	@find jpg/ -name image_\*.jpg -exec rm {} +;
endif

##########################################################
# valgrind commands
# valgrind --vex-guest-max-insns=25 --track-origins=yes ./DDomp

