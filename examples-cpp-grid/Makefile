SHELL=/usr/bin/env bash

mpi_options ?= --oversubscribe --bind-to none

tests = \
		grid-with-qlat

all: run-grid

compile-grid: $(patsubst %, %.compile, $(tests))

run-grid: $(patsubst %, %.run, $(tests))

clean: $(patsubst %, %.clean, $(tests))

compile-proj:
	-mkdir -p build
	-meson setup build >build/log.compile 2>&1
	-time meson compile -C build >>build/log.compile 2>&1

run-proj: compile-proj
	-mkdir build
	-touch build/log build/log.full
	-time meson compile -C build run
	-cat build/log | grep -a 'CHECK: \|INFO: \|WARNING' >build/log.txt
	-cat build/log.txt | grep -a 'CHECK: ' >build/log.check.txt.new
	-cat log | grep -a 'CHECK: ' >build/log.check.txt
	-if diff build/log.check.txt build/log.check.txt.new ; then echo passed ; else cat build/log.full ; fi
	-rm -rfv log
	-touch log
	-cp build/log.txt log
	-cat build/log.compile build/log.full >log.full

clean-proj:
	-rm -rf build log.full

%.compile:
	$(MAKE) -f ../Makefile -C $(patsubst %.compile, %, $@) compile-proj

%.run:
	$(MAKE) -f ../Makefile -C $(patsubst %.run, %, $@) run-proj

%.clean:
	$(MAKE) -f ../Makefile -C $(patsubst %.clean, %, $@) clean-proj
