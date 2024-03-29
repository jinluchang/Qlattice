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
	-mkdir -p build
	-meson setup build >build/log.compile 2>&1
	-time meson compile -C build >>build/log.compile 2>&1

run-proj: compile-proj
	-mkdir build
	-touch build/log build/log.full
	-time meson compile -C build run
	-cat build/log | grep '^CHECK: ' >build/log.check
	-if diff build/log.check log ; then cat build/log.full ; else echo passed ; fi
	-rm -rfv log
	-cp build/log.check log
	-cat build/log.compile build/log.full >log.full

clean-proj:
	-rm -rf build log.full

%.compile:
	$(MAKE) -f ../Makefile -C $(patsubst %.compile, %, $@) compile-proj

%.run:
	$(MAKE) -f ../Makefile -C $(patsubst %.run, %, $@) run-proj

%.clean:
	$(MAKE) -f ../Makefile -C $(patsubst %.clean, %, $@) clean-proj
