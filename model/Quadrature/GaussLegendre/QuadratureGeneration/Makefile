IDIR += -I../../../../
IDIR += -I/usr/local/include/

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
	$(CC) generateGaussLegendre.cpp -o $(OBJS) -std=c++11 $(CFLAGS) $(IDIR)
#$(IDIR) $(LIBS)

##########################################################
# MAKE CLEAN
##########################################################
clean:
		rm -f $(OBJS)
