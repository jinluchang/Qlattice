SHELL=/bin/bash

qlat=$$HOME/qlat-builds/0.1

CC=mpicc -O2 -Wall
CXX=mpic++ -O2 -Wall -std=c++11

INCLUDE=$(qlat)/local/include
LIB=$(qlat)/local/lib

CFLAGS=-fopenmp -O2 -fno-strict-aliasing -Wall
CFLAGS+= -I$(INCLUDE)
CFLAGS+= -I$(INCLUDE)/eigen3

CXXFLAGS=$(CFLAGS)

LDFLAGS=-L$(LIB)
LDFLAGS+= -lgsl -lgslcblas -lm
LDFLAGS+= -lfftw3_omp -lfftw3

all: qlat.x

run: qlat.x
	. $(qlat)/local/setenv.sh ; time mpirun -x OMP_NUM_THREADS=2 --np 16 ./qlat.x

qlat.x: *.C
	. $(qlat)/local/setenv.sh ; time make build

build:
	$(CXX) -o qlat.x $(CXXFLAGS) *.C $(LDFLAGS)

clean:
	rm qlat.x || :

show-info:
	@echo CXX: $(CXX)
	@echo CXXFLAGS: $(CXXFLAGS)
	@echo LDFLAGS: $(LDFLAGS)
