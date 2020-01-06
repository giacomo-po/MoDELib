MKL_INCLUDE = /opt/intel/mkl/include

INTEL_LIB = /opt/intel/lib
MKL_LIB = /opt/intel/mkl/lib

#MKL_INCLUDE=/opt/intel/mkl/include
#MKL_LIB=/opt/intel/mkl/lib/intel64


IDIR += -I /usr/local/include
IDIR += -I $(MKL_INCLUDE)

LIBS += -L $(MKL_LIB)
LIBS += -L $(INTEL_LIB)

OS= $(shell uname -s)

ifeq ($(OS),Darwin)
	LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl
else
	LIBS+= -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -Wl,--end-group 
#-ldl -lpthread -lm
endif


program:
	g++ main.cpp -o test -std=c++11 -fopenmp -Ofast -msse4 $(IDIR) $(LIBS)
	
clean:
		rm -f test


#export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015.2.164/mkl/lib/intel64/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH =/opt/intel/composer_xe_2015.2.164/compiler/lib intel64/:$LD_LIBRARY_PATH