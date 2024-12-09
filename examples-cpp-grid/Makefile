SHELL=/usr/bin/env bash

tests_grid = grid-with-qlat

all: run

all-grid: run-grid

compile-grid: $(patsubst %, %.compile, $(tests_grid))

run-grid: $(patsubst %, %.run, $(tests_grid))

clean: $(patsubst %, %.clean, $(tests_grid))

compile-proj:
	-mkdir -p build
	-meson setup build >build/log.compile 2>&1
	-time meson compile -C build >>build/log.compile 2>&1

run-proj: compile-proj
	-mkdir build
	-touch build/log build/log.full
	-time meson compile -C build run
	-cat build/log | grep '^CHECK: \|^INFO: ' >build/log.txt
	-cat build/log.txt | grep '^CHECK: ' >build/log.check.txt.new
	-cat log | grep '^CHECK: ' >build/log.check.txt
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
