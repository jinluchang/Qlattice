SHELL=/usr/bin/env bash

mpi_options ?= --oversubscribe --bind-to none

tests_gpt = \
						hlbl-muon-line.log \
						grid-qlat-io-test.log \
						gpt-qlat-auto-simple.log \
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
						gpt-qlat-madwf.log

all: run-gpt

run-gpt: update-sources $(tests_gpt)

update-sources:
	-for i in $(tests_gpt) ; do touch "$${i%.log}.py" ; done

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" --test -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep -a "CHECK: \|INFO: \|WARNING" log.full.txt >log.txt ; grep -a "CHECK: " log.txt >log.check.txt.new ; grep -a "CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v "$<".p/"$@".json.new "$@".json
	-cp -v "$<".p/log.txt "$@"

clean:
	-rm -rfv *.p
