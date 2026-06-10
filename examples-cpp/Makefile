SHELL=/usr/bin/env bash

# Run a single test — recommended:
#   ../nixpkgs/run-one-example-cpp.py <name>
#       e.g. ../nixpkgs/run-one-example-cpp.py simple-1
#       See --help for build variant options (--cuda, --cudasupport, --cu, --clang, --pypi).
#
#   The script sets up the environment from a qlat built by nix (via
#   nixpkgs/install-py-local-kernel-with-nix.sh, which creates the
#   ./result-py-local symlink) and then runs make with this Makefile.
#
# Run tests via make (uses whatever qlat is available in the environment):
#   make                          Run all C++ tests (same as make run).
#   make run                      Run all C++ tests.
#   make compile                  Compile all C++ tests.
#   make <name>.run               Run a single test (e.g. make simple-1.run).
#   make <name>.compile           Compile a single test.
#   make clean                    Remove all build artifacts.

mpi_options ?= --oversubscribe --bind-to none bash bind-gpu-qlat.sh

export mpi_options

tests = \
		vector-utils \
		fft-tests \
		smear-tests \
		gf-utils-tests \
		lbl-muon-part \
		rng-state-tests \
		field-rng-tests \
		mixed-size-cg \
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
		ldouble \
		latio

all: run

compile: $(patsubst %, %.compile, $(tests))

run: $(patsubst %, %.run, $(tests))

clean: $(patsubst %, %.clean, $(tests))

compile-proj:
	-mkdir -p build
	-meson setup build >build/log.compile 2>&1
	-cat build/log.compile
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
