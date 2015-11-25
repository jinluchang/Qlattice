SHELL=/bin/bash

lqps=$$HOME/lqps-builds/0.1

CC=mpicc
CXX=mpic++ -std=c++11

INCLUDE=$(lqps)/local/include
LIB=$(lqps)/local/lib

CFLAGS=-O2 -fno-strict-aliasing -Wall
CFLAGS+= -I$(INCLUDE)
CFLAGS+= -I$(INCLUDE)/eigen3

CXXFLAGS=$(CFLAGS)

LDFLAGS=-L$(LIB)
LDFLAGS+= -lfftw3_omp -lfftw3
LDFLAGS+= -llapacke -llapack -lblas -lgfortran
LDFLAGS+= -lgsl -lgslcblas -lm

all:
	. $(lqps)/local/setenv.sh ; make build

run: lqps.x
	. $(lqps)/local/setenv.sh ; time mpirun -x OMP_NUM_THREADS=2 -np 2 ./lqps.x

lqps.x: *.C
	. $(lqps)/local/setenv.sh ; make build

build:
	$(CXX) -o lqps.x $(CXXFLAGS) $(LDFLAGS) *.C

clean:
	rm lqps.x || :

show-info:
	@echo CXX: $(CXX)
	@echo CXXFLAGS: $(CXXFLAGS)
	@echo LDFLAGS: $(LDFLAGS)
