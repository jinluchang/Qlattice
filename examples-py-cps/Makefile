SHELL=/usr/bin/env bash

# Run a single test — recommended:
#   ../nixpkgs/run-one-example-py-cps.py <name>
#       e.g. ../nixpkgs/run-one-example-py-cps.py cps-qlat-io-test
#       See --help for build variant options (--cuda, --cudasupport, --cu, --clang, --pypi).
#
#   The script sets up the environment from a qlat built by nix (via
#   nixpkgs/install-py-local-kernel-with-nix.sh, which creates the
#   ./result-py-local symlink) and then runs make with this Makefile.
#
#   The script copies sources into ./tmp/examples-py-cps/ and runs make there.
#   Logs are in ./tmp/examples-py-cps/<name>.py.p/ — check that directory for
#   log.full.txt, log.txt, and log.check.txt when debugging a test run.
#
# Run tests via make (uses whatever qlat is available in the environment):
#   make                          Run all CPS tests (same as make run-cps).
#   make run-cps                  Run all CPS tests.
#   make <name>.log               Run a single test (e.g. make cps-qlat-io-test.log).
#   make clean                    Remove test build artifacts (*.p directories).

mpi_options ?= --oversubscribe --bind-to none bash bind-gpu-qlat.sh

export mpi_options

tests_cps = \
			cps-qlat-io-test.log \
			sparse-from-cps/log.txt

all: run-cps

run-cps: update-sources $(tests_cps)

update-sources:
	-for i in $(tests_cps) ; do touch "$${i%.log}.py" ; done
	make -C sparse-from-cps update-sources

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" --test -qmp-geom 1 1 1 2 >log.full.txt 2>&1 ; grep -a "CHECK: \|INFO: \|WARNING" log.full.txt >log.txt ; grep -a "CHECK: " log.txt >log.check.txt.new ; grep -a "CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v "$<".p/"$@".json.new "$@".json
	-cp -v "$<".p/log.txt "$@"

sparse-from-cps/log.txt: sparse-from-cps/run.py
	make -C sparse-from-cps

clean:
	-rm -rfv *.p
	make -C sparse-from-cps clean
