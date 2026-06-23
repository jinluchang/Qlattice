SHELL=/usr/bin/env bash

# Run a single test — recommended:
#   ../nixpkgs/run-one-example-py-gpt.py <name>
#       e.g. ../nixpkgs/run-one-example-py-gpt.py qtopo-measure-dwf
#       See --help for build variant options (--cuda, --cudasupport, --cu, --clang, --pypi).
#
#   The script sets up the environment from a qlat built by nix (via
#   nixpkgs/install-py-local-kernel-with-nix.sh, which creates the
#   ./result-py-local symlink) and then runs make with this Makefile.
#
#   The script copies sources into ./tmp/examples-py-gpt/ and runs make there.
#   Logs are in ./tmp/examples-py-gpt/<name>.py.p/ — check that directory for
#   log.full.txt, log.txt, and log.check.txt when debugging a test run.
#
# Run tests via make (uses whatever qlat is available in the environment):
#   make                          Run all GPT tests (same as make run-gpt).
#   make run-gpt                  Run all GPT tests.
#   make <name>.log               Run a single test (e.g. make qtopo-measure-dwf.log).
#   make clean                    Remove test build artifacts (*.p directories).
#
# Before running tests (make only):
#   1. Restore the reference .log and .log.json files from git (required for
#      the CHECK-line comparison that determines pass/fail):
#        git checkout -- *.log *.log.json
#
#   2. If a previous run left a *.p directory, run `make clean` first.
#
#   3. If the target .log already exists from a prior run, touch the .py source
#      file to trigger a rebuild (make considers .log up-to-date otherwise):
#        touch <name>.py
#
# How test verification works:
#   - Each test prints CHECK: lines during execution.
#   - After the run, these lines are extracted and diffed against the committed
#     reference .log file.
#   - A test passes when the CHECK: output matches exactly. A diff failure
#     means the output changed (or the reference .log is missing — run
#     `git checkout -- *.log` to restore it).

mpi_options ?= --oversubscribe --bind-to none bash bind-gpu-qlat.sh

export mpi_options

tests_gpt = \
			qtopo-measure-dwf.log \
			hlbl-muon-line.log \
			grid-qlat-io-test.log \
			gpt-qlat-auto-simple.log \
			gpt-qlat-data-gen-pi0-decay.log \
			gpt-qlat-data-gen-auto.log \
			gpt-qlat-data-gen-hlbl.log \
			gpt-qlat-data-gen-pipi-qed.log \
			gpt-qlat-data-gen-eta-c.log \
			gpt-qlat-data-gen-prop.log \
			gq-auto-checker.log \
			gpt-qlat-convert.log \
			gpt-qlat-sample-gauge-field.log \
			gpt-qlat-smear.log \
			gpt-qlat-mdwf.log \
			gpt-qlat-gfix-mdwf.log \
			gpt-qlat-free-invert.log \
			gpt-qlat-lanc.log \
			gpt-qlat-rbc-ukqcd-invert.log \
			gpt-qlat-madwf.log \
			docs-rbc_ukqcd.log

all: run-gpt

run-gpt: update-sources $(tests_gpt)

update-sources:
	-for i in $(tests_gpt) ; do touch "$${i%.log}.py" ; done

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" --test --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep -a "CHECK: \|INFO: \|WARNING" log.full.txt >log.txt ; grep -a "CHECK: " log.txt >log.check.txt.new ; grep -a "CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v "$<".p/"$@".json.new "$@".json
	-cp -v "$<".p/log.txt "$@"

clean:
	-rm -rfv *.p
