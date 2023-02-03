SHELL=/bin/bash

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
		simple-1 \
		hmc \
		flowed-hmc \
		heatbath \
		template \
		dslash-tests \
		propagators

tests_grid = grid-with-qlat

all: run

all-grid: run-grid

compile: $(patsubst %, %.compile, $(tests))

run: $(patsubst %, %.run, $(tests))

compile-grid: $(patsubst %, %.compile, $(tests_grid))

run-grid: $(patsubst %, %.run, $(tests_grid))

clean: $(patsubst %, %.clean, $(tests) $(tests_grid))

compile-proj:
	-meson setup build
	-time meson compile -C build

run-proj: compile-proj
	-mkdir build
	-touch build/log build/log.full
	-time meson compile -C build run
	-if diff build/log log | grep 'CHECK: ' ; then cat build/log.full ; else echo passed ; fi
	-cp build/log log
	-cp build/log.full log.full

clean-proj:
	-rm -rf build log.full

%.compile:
	$(MAKE) -f ../Makefile -C $(patsubst %.compile, %, $@) compile-proj

%.run:
	$(MAKE) -f ../Makefile -C $(patsubst %.run, %, $@) run-proj

%.clean:
	$(MAKE) -f ../Makefile -C $(patsubst %.clean, %, $@) clean-proj
