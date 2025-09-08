SHELL=/usr/bin/env bash

mpi_options ?= --oversubscribe --bind-to none

tests = \
		qcalc-psel-closest-dis-sqr.log \
		qconvert-psel-txt-to-lati.log \
		qinst-map.log \
		qflow-scale.log \
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
		check-id-node.log \
		prop-smear.log \
		set-rand.log \
		utils.log \
		benchmark-matrix-functions.log \
		make-sample-gauge-field.log \
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
		fields-io.log \
		fsel-psel.log \
		psel-split.log \
		selected-convert-io.log \
		selected-field.log \
		selected-shuffle.log \
		test-qar.log \
		hmc-pure-gauge.log \
		hmc-pure-gauge-flowed.log \
		hmc-pure-gauge-test.log \
		hmc-pure-gauge-flowed-test.log \
		hmc-pure-gauge-symmetric.log \
		hmc-pure-gauge-fa-sym.log \
		hmc-pions.log \
		gf-utils.log \
		gf-flow-topo.log \
		free-invert.log \
		fit-corr.log \
		fit-corr-2.log \
		hlt-reconstruction.log \
		mpi-utils.log \
		mat-mpi-test.log

all: run

run: update-sources $(tests) py-qlat-ext.log

update-sources:
	-for i in $(tests) ; do touch "$${i%.log}.py" ; done
	-touch py-qlat-ext/meson.build

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 60m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" --test -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep -a "CHECK: \|INFO: \|WARNING" log.full.txt >log.txt ; grep -a "CHECK: " log.txt >log.check.txt.new ; grep -a "CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v "$<".p/"$@".json.new "$@".json
	-cp -v "$<".p/log.txt "$@"

py-qlat-ext.log: py-qlat-ext/meson.build
	-touch py-qlat-ext/meson.build
	-rm -rfv py-qlat-ext.p
	-mkdir -p py-qlat-ext.p
	-time meson setup py-qlat-ext.p/py-qlat-ext-build py-qlat-ext --prefix="$$PWD"/py-qlat-ext.p/py-qlat-ext-install -Dpython.platlibdir="$$PWD"/py-qlat-ext.p/py-qlat-ext-install/lib/python3/dist-packages -Dpython.purelibdir="$$PWD"/py-qlat-ext.p/py-qlat-ext-install/lib/python3/dist-packages
	-time meson compile -C py-qlat-ext.p/py-qlat-ext-build -j$(num_proc)
	-time meson install -C py-qlat-ext.p/py-qlat-ext-build
	-time ( cd py-qlat-ext.p ; PYTHONPATH="$$PWD/py-qlat-ext-install/lib/python3/dist-packages:$$PYTHONPATH" python3 py-qlat-ext-install/bin/run.py --test >log.full.txt 2>&1 ; grep -a "CHECK: " log.full.txt >log.check.txt.new ; grep -a "CHECK: " ../"$@" >log.check.txt )
	-if diff py-qlat-ext.p/log.check.txt py-qlat-ext.p/log.check.txt.new ; then echo passed ; else cat py-qlat-ext.p/log.full.txt ; false ; fi
	-rm -rfv "$@"
	-touch "$@"
	-cp -v py-qlat-ext.p/log.check.txt.new "$@"

clean:
	-rm -rfv *.p
