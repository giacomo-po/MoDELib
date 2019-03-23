USE_ICC = 0

# MODEL LIBRARY #################################
 MODEL_INCLUDE = /Users/giacomo/Documents/MODEL
# MODEL_INCLUDE = /home/giacomo/Dropbox/Model

# MPI LIBRARY (PARALLEL ONLY) ###################
#MPI_INCLUDE = /usr/local/include
MPI_LIB = /usr/local/lib

# METIS LIBRARY (PARALLEL ONLY) #################
#METIS_INCLUDE =/usr/local/include
#METIS_LIB =/usr/local/lib

##########################################################
# INCLUDE SETTINGS - DO NOT EDIT-
##########################################################
IDIR  = -I/usr/local/include
IDIR += -I./
IDIR += -I$(MODEL_INCLUDE)
#IDIR += -I$(MPI_INCLUDE)
#IDIR += -I$(METIS_INCLUDE)

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
CFLAGS	+= -std=c++0x
else
CC= g++
# optimization level = fast
CFLAGS  += -O3
# enable openMP
#CFLAGS += -fopenmp
# enable C++0x
CFLAGS += -std=c++0x
endif

# enable common warnings… 
CFLAGS += -Wall
# … and more warnings
CFLAGS += -Wextra

#disable model-assertions
CFLAGS += -D _MODEL_NO_DEBUG_
#disable model-input-checks
CFLAGS += -D _MODEL_NO_BAD_INPUT_CHECK_

OBJS 	= CPserial

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
mpiCFLAGS  += -O3
# enable openMP
#mpiCFLAGS += -fopenmp
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
mpiCFLAGS += -D _MODEL_MPI_

# Link to Metis partitioner
#mpiLIBS += -L$(METIS_LIB) 
#mpiLIBS +=	 -lmetis

mpiOBJS = CPparallel

parallel:
	@echo ${blueColor}making parallel version... ${normalColor}
	$(mpiCC) main.cpp -o $(mpiOBJS) $(mpiCFLAGS) $(IDIR) $(mpiLIBS)


##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS) $(mpiOBJS)


