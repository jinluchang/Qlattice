SHELL=/bin/bash

prefix=$$HOME/qlat-builds/0.1

CC=mpicc -O2 -Wall
CXX=mpic++ -O2 -Wall -std=c++0x

INCLUDE=$(prefix)/include
LIB=$(prefix)/lib

CFLAGS=-fopenmp -O2 -fno-strict-aliasing -Wall
CFLAGS+= -I$(INCLUDE)
CFLAGS+= -I$(INCLUDE)/eigen3

CXXFLAGS=$(CFLAGS)

LDFLAGS=-L$(LIB)
LDFLAGS+= -lgsl -lgslcblas -lm
LDFLAGS+= -lfftw3_omp -lfftw3
LDFLAGS+= -lhash-cpp

all: qlat.x

run: qlat.x
	. $(prefix)/setenv.sh ; time mpirun -x OMP_NUM_THREADS=8 --np 16 ./qlat.x
	make clean

qlat.x: *.C
	. $(qlat)/local/setenv.sh ; time make build
	[ -f $@ ]

build:
	$(CXX) -o qlat.x $(CXXFLAGS) *.C $(LDFLAGS) 2>&1 | grep --color=always 'error:\|' || true

clean:
	rm qlat.x || :

show-info:
	@echo CXX: $(CXX)
	@echo CXXFLAGS: $(CXXFLAGS)
	@echo LDFLAGS: $(LDFLAGS)
