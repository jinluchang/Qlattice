SHELL=/usr/bin/env bash

tests_cps = \
			cps-qlat-io-test.log \
			sparse-from-cps/log.txt

all: run-cps

run-cps: update-sources $(tests_cps)

update-sources:
	-for i in $(tests_cps) ; do touch "$${i%.log}.py" ; done
	-touch qlat-ext/meson.build
	make -C sparse-from-cps update-sources

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep "^CHECK: \|^INFO: " log.full.txt >log.txt ; grep "^CHECK: " log.txt >log.check.txt.new ; grep "^CHECK: " ../"$@" >log.check.txt )
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
