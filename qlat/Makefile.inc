ifeq ($(strip $(QLAT_CXX)),)
	CXX=mpic++
else
	CXX=$(QLAT_CXX)
endif

ifeq ($(strip $(QLAT_CXXLD)),)
	CXXLD=mpic++ --shared
else
	CXXLD=$(QLAT_CXXLD)
endif

ifeq ($(strip $(QLAT_FLAGS)),)
	FLAGS=-std=c++11 -fopenmp -O2 -Wall -fno-strict-aliasing -fPIC
else
	FLAGS=$(QLAT_FLAGS)
endif

QLAT_INCLUDE=$(QLAT_PREFIX)/include
QLAT_LIB=$(QLAT_PREFIX)/lib

QLAT_CFLAGS=
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/eigen3
QLAT_CFLAGS+= -I$(QLAT_INCLUDE)/qutils

QLAT_CXXFLAGS=$(QLAT_CFLAGS)

QLAT_LDFLAGS=
QLAT_LDFLAGS+= -L$(QLAT_LIB)
QLAT_LDFLAGS+= -l:libfftw3.a -l:libfftw3f.a
QLAT_LDFLAGS+= -l:libcuba.a
QLAT_LDFLAGS+= -l:libz.a
QLAT_LDFLAGS+= -lm

CXXFLAGS=
CXXFLAGS+= $(QLAT_CXXFLAGS)
LDFLAGS=
LDFLAGS+= $(QLAT_LDFLAGS)
