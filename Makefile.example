SHELL=/bin/bash

CXX=mpic++ -std=c++0x -fopenmp -march=native -O2 -Wall
# -DQLAT_NO_DEFAULT_TYPE
# -DSKIP_ASSERT
# -DOLD_CPP
# -DNO_OMP

qlat_prefix=$$HOME/qlat-build/1.0

QLAT_INCLUDE=$(qlat_prefix)/include
QLAT_LIB=$(qlat_prefix)/lib

QLAT_CFLAGS=-fno-strict-aliasing
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/eigen3
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/qutils

QLAT_CXXFLAGS=$(QLAT_CFLAGS)

QLAT_LDFLAGS=$(QLAT_LIB)/libfftw3.a $(QLAT_LIB)/libcuba.a $(QLAT_LIB)/libz.a -lm

all: qlat.x

run: qlat.x
	time mpirun -x OMP_NUM_THREADS=2 --np 2 ./qlat.x |& tee log
	# make clean

diff:
	-git diff | egrep -v Timer | egrep '^\+|^\-'

qlat.x: *.C
	time make build
	[ -f $@ ]

build:
	-$(CXX) -o qlat.x $(QLAT_CXXFLAGS) *.C $(QLAT_LDFLAGS) 2>&1 | grep --color=always 'error:\|'

clean:
	-rm qlat.x
	-rm -rf results*
	-rm -rf huge-data

dist-clean: clean
	-git checkout log

show-info:
	@echo CXX: $(CXX)
	@echo CXXFLAGS: $(QLAT_CXXFLAGS)
	@echo LDFLAGS: $(QLAT_LDFLAGS)
