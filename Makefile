SHELL=/bin/bash

lqps=$$HOME/lqps-builds/0.1

CC=mpicc
CXX=mpic++ -std=c++11

INCLUDE=$(lqps)/local/include
LIB=$(lqps)/local/lib

CFLAGS=-fopenmp -O2 -fno-strict-aliasing -Wall
CFLAGS+= -I$(INCLUDE)
CFLAGS+= -I$(INCLUDE)/eigen3

CXXFLAGS=$(CFLAGS)

LDFLAGS=-L$(LIB)
LDFLAGS+= -llapacke -llapack -lblas -lgfortran
LDFLAGS+= -lgsl -lgslcblas -lm
LDFLAGS+= -lfftw3_omp -lfftw3

all: lqps.x

run: lqps.x
	. $(lqps)/local/setenv.sh ; time mpirun -x OMP_NUM_THREADS=2 -np 2 ./lqps.x

lqps.x: *.C
	. $(lqps)/local/setenv.sh ; time make build

build:
	$(CXX) -o lqps.x $(CXXFLAGS) *.C $(LDFLAGS)

clean:
	rm lqps.x || :

show-info:
	@echo CXX: $(CXX)
	@echo CXXFLAGS: $(CXXFLAGS)
	@echo LDFLAGS: $(LDFLAGS)
