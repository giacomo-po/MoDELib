MODEL_INCLUDE = ../../../
usePARDISO = 0
MKL_PATH=/opt/intel/composer_xe_2013.3.171/mkl

##############################
IDIR += -I ./
IDIR += -I /usr/local/include
IDIR += -I $(MODEL_INCLUDE)

##############################
CFLAGS += -std=c++11  
CFLAGS += -O3 
CFLAGS += -fopenmp

##############################
LIBS += 

##############################
# Settings for PARDISO solver
ifeq ($(usePARDISO), 1)
CFLAGS += -D _MODEL_PARDISO_SOLVER_
IDIR += -I $(MKL_PATH)/include
LIBS += -L $(MKL_PATH)/lib
LIBS += -lmkl_intel_lp64 
LIBS += -lmkl_intel_thread 
LIBS += -lmkl_core
endif


program:
	g++ main.cpp -o fem $(CFLAGS) $(IDIR)  $(LIBS)
	
clean:
		rm -f fem
		
		
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2013.3.171/mkl/lib:$DYLD_LIBRARY_PATH
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2013.3.171/compiler/lib:$DYLD_LIBRARY_PATH
