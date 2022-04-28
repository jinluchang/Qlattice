ifeq ($(strip $(QLAT_MPICXX)),)
	ifeq ($(strip $(MPICXX)),)
		MPICXX = mpic++ -std=c++11 -fopenmp -O2 -Wall -fno-strict-aliasing
	else
		MPICXX := $(MPICXX) -std=c++11 -fopenmp -O2 -Wall -fno-strict-aliasing
	endif
else
	MPICXX=$(QLAT_MPICXX)
endif

ifeq ($(strip $(QLAT_CXXFLAGS)),)
	CXXFLAGS +=
else
	CXXFLAGS = $(QLAT_CXXFLAGS)
endif

ifeq ($(strip $(QLAT_LDFLAGS)),)
	LDFLAGS +=
else
	LDFLAGS = $(QLAT_LDFLAGS)
endif

ifeq ($(strip $(QLAT_LIBS)),)
	LIBS +=
else
	LIBS = $(QLAT_LIBS)
endif

QLAT_INCLUDE = $(QLAT_PREFIX)/include
QLAT_LIB = $(QLAT_PREFIX)/lib

CXXFLAGS += -I$(QLAT_INCLUDE)

LDFLAGS += -L$(QLAT_LIB)

LIBS += -lfftw3_mpi -lfftw3f_mpi -lfftw3 -lfftw3f
LIBS += -lcuba
LIBS += -lgsl -lgslcblas
LIBS += -lz
LIBS += -lm
