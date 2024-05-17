# Master Makefile for Sequence Utilities for Genome analysis
# Type make to build
include make.inc

O = o
SRC_DIR = ./

LIBRARY_DIR = -L./ 
INCLUDE = -I. 
CFLAGS = -g                 		                     # compiler switches to be applied to every module
OPTIM_SPEED = -O3             	  	                     # switches that give speed priority over size
OPTIM_SIZE = -O1              	  	                     # switches that give size priority over speed

OPTIONS = $(CFLAGS) $(INCLUDE)

ARCH     = ar
ARCHFLAGS= cr
RAN_LIB   = libranlib.a

all: sim_read_hash 


GENERATOR = \
    ranlib.o \
    com.o  \
    linpack.o



$(RAN_LIB) : $(GENERATOR)
	$(ARCH) $(ARCHFLAGS) $(RAN_LIB)  $(GENERATOR)                  
	ranlib $(RAN_LIB)





sim_read_hash:   $(RAN_LIB)
	upcxx  -o sim_read_hash -g sim_read_hash.cpp -L. -lranlib



%.o: %.c
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<





