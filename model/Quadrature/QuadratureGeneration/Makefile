# model_DIR = /Users/giacomo/Dropbox/Model
# model_DIR = /home/giacomo/Dropbox/Model
# IDIR  = -I/usr/local/include
# IDIR += -I./
# IDIR += -I$(model_DIR)

##########################################################
# COMPILER SETTINGS
##########################################################
CC	= g++
CFLAGS  = -O3
CFLAGS += -Wall
CFLAGS += -Wextra




##########################################################
# MAKE PROGRAM
##########################################################
PROGRAM = generateGaussLegendre
OBJS 	= generateGaussLegendre

$(PROGRAM): 
	$(CC) generateGaussLegendre.cpp -o $(OBJS)  $(CFLAGS)  -I./
#$(IDIR) $(LIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)
