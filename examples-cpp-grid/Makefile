SHELL=/usr/bin/env bash

# Run a single test — recommended:
#   ../nixpkgs/run-one-example-cpp-grid.py <name>
#       e.g. ../nixpkgs/run-one-example-cpp-grid.py grid-with-qlat
#       See --help for build variant options (--cuda, --cudasupport, --cu, --clang, --pypi).
#
#   The script sets up the environment from a qlat built by nix (via
#   nixpkgs/install-py-local-kernel-with-nix.sh, which creates the
#   ./result-py-local symlink) and then runs make with this Makefile.
#
# Run tests via make (uses whatever qlat is available in the environment):
#   make                          Run all Grid tests (same as make run-grid).
#   make run-grid                 Run all Grid tests.
#   make compile-grid             Compile all Grid tests.
#   make <name>.run               Run a single test (e.g. make grid-with-qlat.run).
#   make <name>.compile           Compile a single test.
#   make clean                    Remove all build artifacts.

mpi_options ?= --oversubscribe --bind-to none bash bind-gpu-qlat.sh

export mpi_options

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
