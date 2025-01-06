SHELL=/usr/bin/env bash

mpi_options ?= --oversubscribe --bind-to none

tests = qlat-vector-utils \
		qlat-fft-tests \
		qlat-smear-tests \
		qcd-utils-tests \
		lbl-muon-part \
		rng-state-tests \
		field-rng-tests \
		benchmark \
		dist-io \
		selected-field \
		fields-io \
		fields-io-simple-demo \
		simple-1 \
		hmc \
		flowed-hmc \
		heatbath \
		template \
		dslash-tests \
		propagators \
		ldouble

all: run

compile: $(patsubst %, %.compile, $(tests))

run: $(patsubst %, %.run, $(tests))

clean: $(patsubst %, %.clean, $(tests))

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
