MODEL_INCLUDE = ../../../
EIGEN_INCLUDE =/Users/giacomo/Documents/eigen
usePARDISO = 0
MKL_INCLUDE=/opt/intel/composer_xe_2015.2.132/mkl/include
MKL_LIB=/opt/intel/composer_xe_2015.2.132/mkl/lib

##############################
IDIR += -I ./
IDIR += -I $(MODEL_INCLUDE)
IDIR += -I $(EIGEN_INCLUDE)

##############################
CFLAGS += -std=c++11  
CFLAGS += -O3 
CFLAGS += -fopenmp
CFLAGS += -Wall
CFLAGS += -Wextra

##############################
LIBS += 

##############################
# Settings for PARDISO solver
ifeq ($(usePARDISO), 1)
CFLAGS += -D _MODEL_PARDISO_SOLVER_
IDIR += -I $(MKL_INCLUDE)
LIBS += -L $(MKL_LIB)
LIBS += -lmkl_intel_lp64 
LIBS += -lmkl_intel_thread 
LIBS += -lmkl_core
endif


program:
	g++ main.cpp -o fem $(CFLAGS) $(IDIR)  $(LIBS)
	
clean:
		rm -f fem
		
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.132/mkl/lib:$DYLD_LIBRARY_PATH
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.132/compiler/lib:$DYLD_LIBRARY_PATH

# Hoffman2
#export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
