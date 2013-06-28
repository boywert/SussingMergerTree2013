
CC = gcc
FLAG = -g -lm -fopenmp #-Wall
BASE = $(CURDIR)
SRC = $(CURDIR)/src
LIB = $(CURDIR)/lib
SRCALL = $(wildcard $(SRC)/*.c)
SO_FILES = $(LIB)/gadgetPyIO.so
OPT = 
OPT += -DREADPARTICLE 
#OPT += -DSUBFINDOUT
OPT += -DSTORESUBIDINMOSTBOUNDID
CC += $(OPT)

all: read $(SO_FILES) 

read:  $(SRCALL)
	$(CC) $(FLAG) $(SRCALL) -o read

$(LIB)/gadgetPyIO.so: $(SRC)/gadgetpyio-v1.0/setup.py $(SRC)/gadgetpyio-v1.0/gadgetmodule.c
	cd $(SRC)/gadgetpyio-v1.0 && python setup.py build &&  python setup2.py build && cp build/lib.linux-x86_64-2.6/*.so $(LIB) && cd $(BASE)
clean:
	rm -f read
	rm -f $(LIB)/*
