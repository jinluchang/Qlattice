SHELL=/bin/bash

CXX=mpic++ -std=c++0x -fopenmp -march=native -O2 -Wall

qlat_prefix=$$HOME/qlat-build/default

QLAT_INCLUDE=$(qlat_prefix)/include
QLAT_LIB=$(qlat_prefix)/lib

QLAT_CFLAGS=-std=c++0x -fno-strict-aliasing -Wfatal-errors
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/eigen3
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/qutils

QLAT_CXXFLAGS=$(QLAT_CFLAGS)

QLAT_LDFLAGS=$(QLAT_LIB)/libfftw3.a $(QLAT_LIB)/libz.a $(QLAT_LIB)/libcuba.a -lgsl -lgslcblas -lm


all: qlat.x

qlat.x: *.C
	-$(CXX) -o qlat.x $(QLAT_CXXFLAGS) *.C $(QLAT_LDFLAGS) 2>&1 | grep --color=always 'error:\|'

clean:
	rm *.o qlat.x
