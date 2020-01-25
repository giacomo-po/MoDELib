##########################################################
# EDIT HERE
##########################################################
USE_ICC = 0
usePARDISO = 1
useGREATWHITE = 0
MODEL_INCLUDE =../../../..
GREATWHITE_INCLUDE= $(MODEL_INCLUDE)/../GreatWhite
EIGEN_INCLUDE =/usr/local/include
MKL_INCLUDE=/opt/intel/mkl/include


##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
include ../../buildSymlink.mk
IDIR += -I./
IDIR += -I$(all_header_dir)
#IDIR += -I$(MODEL_INCLUDE)
IDIR += -I$(EIGEN_INCLUDE)
IDIR += -I$(GREATWHITE_INCLUDE)

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
	CFLAGS	+= -qopenmp 
	CFLAGS	+= -std=c++17
else
	CC= g++
	mpiCC = mpicxx
	CFLAGS += -Ofast
	CFLAGS += -fopenmp
	CFLAGS += -std=c++17
	ifeq ($(OS),Darwin)
		CFLAGS += -msse4
	endif

	ifeq ($(OS),Linux)
#		CFLAGS += -msse4
		CFLAGS += -march=native
	endif

endif


# enable warnings
CFLAGS += -Wall 
CFLAGS += -Wextra

############################################
# Compiler flags for MKL PARDISO, see also
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
ifeq ($(usePARDISO), 1)
	CFLAGS += -D _MODEL_PARDISO_SOLVER_
	IDIR += -I $(MKL_INCLUDE)

	ifeq ($(OS),Darwin)
		MKL_LIB=/opt/intel/mkl/lib
		INTEL_LIB=/opt/intel/lib
		LIBS += -L $(MKL_LIB) -L $(INTEL_LIB)
#		LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread 
		LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl
# export DYLD_LIBRARY_PATH=/opt/intel/mkl/lib:$DYLD_LIBRARY_PATH
# export DYLD_LIBRARY_PATH=/opt/intel/lib:$DYLD_LIBRARY_PATH
	endif

	ifeq ($(OS),Linux)
		MKL_LIB=/opt/intel/mkl/lib/intel64
		LIBS += -L $(MKL_LIB)
		ifeq ($(USE_ICC), 1)
			LIBS += -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -Wl,--end-group 
		else
			LIBS += -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -Wl,--end-group 
		endif
# export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/opt/intel/lib/intel64:$LD_LIBRARY_PATH
	endif

endif

ifeq ($(useGREATWHITE), 1)
	CFLAGS += -D _MODEL_GREATWHITE_
endif


##########################################################
# MAKE ALL
##########################################################
all:
	@echo please use \
	make DDomp \
	or make DDmpi \
	or make microstructureGenerator \
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
# MAKE MICROSTRUCTUREGENERATOR
##########################################################
generatorOBJS = microstructureGenerator

$(generatorOBJS):
	@echo making microstructureGenerator
	$(CC) $(MODEL_INCLUDE)/tools/MicrostructureGenerators/microstructureGenerator.cpp -o microstructureGenerator -O3 -std=c++17 -Wall $(IDIR)


##########################################################
# MAKE DDCONVERTER
##########################################################
converterOBJS = DDconverter

$(converterOBJS):
	@echo making DDconverter
	$(CC) $(MODEL_INCLUDE)/tools/DDconverter/DDconverter.cpp -o DDconverter -O3 -std=c++17 -lstdc++fs $(IDIR)

##########################################################
# MAKE EMPTY
##########################################################
#mkdirV=$(shell test -d V || mkdir V)
#mkdirE=$(shell test -d E || mkdir E)
Eexists=$(shell test -d E && echo 1)
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
	@echo emptying folder evl
	@find evl/ -name evl_\*.txt ! -name evl_0.txt -exec rm {} +;
	@find evl/ -name evl_\*.bin ! -name evl_0.bin -exec rm {} +;
	@find evl/ -name ddAux_\*.txt -exec rm {} +;
	@find evl/ -name ddAux_\*.bin -exec rm {} +;
ifeq ($(Eexists),1)
	@echo emptying folder E
	@find E/ -name E_\*.txt ! -name E_0.txt -exec rm {} +;
	@find E/ -name E_\*.bin ! -name E_0.bin -exec rm {} +;
endif
ifeq ($(Fexists),1)
	@echo emptying folder F
	@find F/ -name F_\*.txt -exec rm {} +;
	@find F/ -name F_\*.bin -exec rm {} +;
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

