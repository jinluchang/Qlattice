SHELL=/usr/bin/env bash

mpi_options ?= --oversubscribe --bind-to none

tests = \
		auto-cexprs.log \
		auto-cexprs-more.log \
		auto-cexprs-kpipi.log \
		auto-cexprs-kpipi-qed.log \
		auto-cexprs-uu.log \
		auto-cexprs-cy.log \
		auto-cexprs-more-cy.log \
		auto-cexprs-kpipi-cy.log \
		auto-cexprs-kpipi-qed-cy.log \
		auto-cexprs-uu-cy.log \
		auto-cexprs-distillation.log \
		auto-contract-01.log \
		auto-contract-02.log \
		auto-contract-03.log \
		auto-contract-04.log \
		auto-contract-05.log \
		auto-contract-06.log \
		auto-contract-07.log \
		mat-op-distillation.log \
		mat-op-qcd.log \
		qlat-id-node.log \
		set-rand.log \
		utils.log \
		benchmark-matrix-functions.log \
		make-sample-gauge-field.log \
		selected-convert-io.log \
		qar.log \
		lat-io.log \
		lat-io-test1.log \
		lat-io-test2.log \
		jackknife-random.log \
		jackknife-super.log \
		jackknife-random-bench.log \
		field-utils.log \
		field-analysis.log \
		field-multi-size.log \
		field-selection.log \
		qs-fsel-psel.log \
		selected-field.log \
		fields-io.log \
		hmc-pure-gauge.log \
		hmc-pure-gauge-flowed.log \
		hmc-pure-gauge-test.log \
		hmc-pure-gauge-flowed-test.log \
		hmc-pure-gauge-symmetric.log \
		hmc-pure-gauge-fa-sym.log \
		hmc-pions.log \
		qcd-utils.log \
		free-invert.log \
		q-fit-corr.log \
		q-fit-corr-2.log \
		q-hlt-reconstruction.log \
		mpi-utils.log \
		mat-mpi-test.log

all: run

run: update-sources $(tests) qlat-ext.log

update-sources:
	-for i in $(tests) ; do touch "$${i%.log}.py" ; done
	-touch qlat-ext/meson.build

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep -a "^CHECK: \|^INFO: \|WARNING" log.full.txt >log.txt ; grep -a "^CHECK: " log.txt >log.check.txt.new ; grep -a "^CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v "$<".p/"$@".json.new "$@".json
	-cp -v "$<".p/log.txt "$@"

qlat-ext.log: qlat-ext/meson.build
	-touch qlat-ext/meson.build
	-rm -rfv qlat-ext.p
	-mkdir -p qlat-ext.p
	-time meson setup qlat-ext.p/qlat-ext-build qlat-ext --prefix="$$PWD"/qlat-ext.p/qlat-ext-install -Dpython.platlibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages -Dpython.purelibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages
	-time meson compile -C qlat-ext.p/qlat-ext-build -j$(num_proc)
	-time meson install -C qlat-ext.p/qlat-ext-build
	-time ( cd qlat-ext.p ; PYTHONPATH="$$PWD/qlat-ext-install/lib/python3/dist-packages:$$PYTHONPATH" python3 qlat-ext-install/bin/run.py >log.full.txt 2>&1 ; grep -a "CHECK: " log.full.txt >log.check.txt.new ; grep -a "^CHECK: " ../"$@" >log.check.txt )
	-if diff qlat-ext.p/log.check.txt qlat-ext.p/log.check.txt.new ; then echo passed ; else cat qlat-ext.p/log.full.txt ; false ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v qlat-ext.p/log.check.txt.new "$@"

clean:
	-rm -rfv *.p
