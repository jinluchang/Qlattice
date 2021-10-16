SHELL=/bin/bash

ifeq ($(strip $(QLAT_PREFIX)),)
	QLAT_PREFIX=$(shell readlink -m "$$HOME/qlat-build/default")
endif

include $(QLAT_PREFIX)/include/qlat/Makefile.inc

CXXFLAGS+= $$(python3-config --includes)
# CXXFLAGS+= $$(python3 -c "import numpy; print(numpy.get_include())") 

NAME=cqlat.so

OBJS:=$(patsubst %.C, %.o, $(wildcard *.C))

all: build

build: $(OBJS)
	$(MPICXX) $(LDFLAGS) -o $(NAME) *.o $(LIBS)
	ls -lh $(NAME)
	mv $(NAME) ..

exports.h: *.C
	./update.sh

%.o : %.C exports.h
	$(MPICXX) $(CXXFLAGS) -c -o $@ $<
	@echo "$@" "$<"

clean:
	rm *.o
