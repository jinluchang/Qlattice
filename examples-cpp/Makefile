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
		heatbath \
		template \
		dslash-tests \
		propagators

all: run

compile: $(patsubst %, %.compile, $(tests))

run: $(patsubst %, %.run, $(tests))

clean: $(patsubst %, %.clean, $(tests))

compile-proj:
	meson setup build
	time ninja -C build

run-proj: compile-proj
	-rm log
	-time ninja -C build run
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
