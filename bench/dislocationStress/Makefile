##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR  = -I/usr/local/include
IDIR += -I./../../


serial:
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
CC= icpc
# optimization level
CFLAGS	 = -O2
CFLAGS	+= -fast
# enable openMP
CFLAGS	+= -openmp
#enable C++0x
CFLAGS	+= -std=c++0x
else
CC = g++
# optimization level = fast
CFLAGS  += -Ofast
CFLAGS  += -msse2
# enable openMP
CFLAGS += -fopenmp
# enable C++0x
CFLAGS += -std=c++0x
endif

# enable common warnings… 
CFLAGS += -Wall
# … and more warnings       
CFLAGS += -Wextra

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
mpiCFLAGS += -O2
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

