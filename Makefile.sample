
CC = gcc
FLAG = -g -lm -fopenmp #-Wall
BASE = $(CURDIR)
SRC = $(CURDIR)/src
LIB = $(CURDIR)/lib
#SRCALL = $(wildcard $(SRC)/*.c)
SRCALL = $(SRC)/allvars.h  $(SRC)/analysis.c  $(SRC)/read_catalogues.c  $(SRC)/read_output.c   $(SRC)/sorting.c  $(SRC)/utils.c
SO_FILES = $(LIB)/gadgetPyIO.so
OPT = 
#OPT += -DAQUARIUS000
OPT += -DREADPARTICLE 
#OPT += -DSUBFINDOUT
OPT += -DSTORESUBIDINMOSTBOUNDID
#OPT += -DHBTEXCLUSIVE
#OPT += -DRESETPARTICLECACHE
CC += $(OPT)

all: read removelowres $(SO_FILES) 

read:  $(SRCALL) $(SRC)/main.c
	$(CC) $(FLAG) $(SRCALL) $(SRC)/main.c -o read

removelowres:  $(SRCALL) $(SRC)/removelowres.c
	$(CC) $(FLAG) $(SRCALL) $(SRC)/removelowres.c -o removelowres

$(LIB)/gadgetPyIO.so: $(SRC)/gadgetpyio-v1.0/setup.py $(SRC)/gadgetpyio-v1.0/gadgetmodule.c
	cd $(SRC)/gadgetpyio-v1.0 && python setup.py build &&  python setup2.py build && cp build/lib.linux-x86_64-2.6/*.so $(LIB) && cd $(BASE)
clean:
	rm -f read
	rm -f $(LIB)/*
