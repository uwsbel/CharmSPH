SHELL := /bin/bash

# to be set appropiately
CHARMC         = $(CHARMDIR)/bin/charmc

OPTS            = -O3

DECL=
SUFFIX=

ifeq ($(FT), yes)
DECL=-DCMK_MEM_CHECKPOINT=1
FT=-syncft
endif

all: charmsph

charmsph: Main.o Cell.o Compute.o charmsph.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-language charm++ -o charmsph$(SUFFIX) Main.o Cell.o Compute.o

Main.o: Main.cc Main.h charmsph.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.cc

Cell.o: Cell.cc Cell.h charmsph.decl.h defs.h
	$(CHARMC) $(OPTS) -o Cell.o Cell.cc

charmsph.decl.h:	charmsph.ci
	$(CHARMC) -E charmsph.ci $(DECL)

Compute.o: Compute.cc Compute.h charmsph.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.cc

test: charmsph
	./charmrun +p4 ./charmsph 4 4 4 10 3 3 +balancer GreedyLB +LBDebug 1 ++local

clean:
	rm -f *.decl.h *.def.h *.o charmsph charmsph-ft charmsph.prj charmrun
