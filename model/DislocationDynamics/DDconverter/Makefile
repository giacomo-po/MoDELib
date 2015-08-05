IDIR  = -I/usr/local/include/
IDIR += -I/Users/giacomo/Documents/MODEL/


##########################################################
# COMPILER SETTINGS
##########################################################
CC	= g++
#CC	= icc

ifeq ($(CC),g++)
	# optimization level
	CFLAGS	 = -O3
	CFLAGS	+= -std=c++11
endif


ifeq ($(CC),icc)
	CFLAGS	 = -O3
	CFLAGS	+= -std=c++11
endif

##########################################################
# PROGRAM
##########################################################
PROGRAM =	DDconverter
OBJS 	= 	DDconverter


$(PROGRAM): 
	$(CC) DDconverter.cpp -o $(OBJS) $(CFLAGS) $(IDIR)  $(LIBS) 
	
clean:
		rm -f $(OBJS)

#	
#	