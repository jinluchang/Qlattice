ifeq ($(strip $(QLAT_MPICXX)),)
	ifeq ($(strip $(MPICXX)),)
		MPICXX=mpic++ -std=c++11 -fopenmp -O2 -Wall -fno-strict-aliasing
	else
		MPICXX=$(MPICXX) -std=c++11 -fopenmp -O2 -Wall -fno-strict-aliasing
	endif
else
	MPICXX=$(QLAT_MPICXX)
endif

ifeq ($(strip $(QLAT_CXXFLAGS)),)
	ifeq ($(strip $(CXXFLAGS)),)
		CXXFLAGS=-fPIC
	endif
else
	CXXFLAGS=$(QLAT_CXXFLAGS)
endif

ifeq ($(strip $(QLAT_LDFLAGS)),)
	ifeq ($(strip $(LDFLAGS)),)
		LDFLAGS=--shared
	endif
else
	LDFLAGS=$(QLAT_LDFLAGS)
endif

ifeq ($(strip $(QLAT_LIBS)),)
	ifeq ($(strip $(LIBS)),)
		LIBS=
	endif
else
	LIBS=$(QLAT_LIBS)
endif

QLAT_INCLUDE=$(QLAT_PREFIX)/include
QLAT_LIB=$(QLAT_PREFIX)/lib

CXXFLAGS+= -I$(QLAT_INCLUDE)

LDFLAGS+= -L$(QLAT_LIB)

LIBS+= -lfftw3 -lfftw3f
LIBS+= -lcuba
LIBS+= -lz
LIBS+= -lm
