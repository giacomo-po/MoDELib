#MKLPATH=/opt/intel/composer_xe_2013.3.171/mkl


# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2013.3.171/mkl/lib:$DYLD_LIBRARY_PATH
# export DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2013.3.171/compiler/lib:$DYLD_LIBRARY_PATH



#LIBS += -L $(MKLPATH)/lib
#LIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

#-I $(MKLPATH)/include $(LIBS)

MODEL_INCLUDE = ../../../

program:
	g++ main.cpp -o fem -std=c++11  -O3 -fopenmp -I./ -I/usr/local/include -I $(MODEL_INCLUDE)  
	
clean:
		rm -f fem