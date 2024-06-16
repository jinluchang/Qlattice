SHELL=/bin/bash

all: run

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
		mat-op-distillation.log \
		mat-op-qcd.log \
		qlat-id-node.log \
		set-rand.log \
		benchmark-matrix-functions.log \
		make-sample-gauge-field.log \
		qar.log \
		lat-io.log \
		lat-io-test1.log \
		lat-io-test2.log \
		jackknife-random.log \
		jackknife-super.log \
		jackknife-random-bench.log \
		field-utils.log \
		field-multi-size.log \
		field-selection.log \
		selected-field.log \
		fields-io.log \
		hmc-pure-gauge.log \
		hmc-pure-gauge-flowed.log \
		hmc-pure-gauge-test.log \
		hmc-pure-gauge-flowed-test.log \
		hmc-pure-gauge-symmetric.log \
		hmc-pions.log \
		qcd-utils.log \
		free-invert.log \
		q-fit-corr.log \
		q-hlt-reconstruction.log \
		mpi-utils.log \
		mat-mpi-test.log

tests_gpt = \
			hlbl-muon-line.log \
			grid-qlat-io-test.log \
			gpt-qlat-data-gen-auto.log \
			gpt-qlat-data-gen-hlbl.log \
			gpt-qlat-convert.log \
			gpt-qlat-sample-gauge-field.log \
			gpt-qlat-smear.log \
			gpt-qlat-mdwf.log \
			gpt-qlat-gfix-mdwf.log \
			gpt-qlat-free-invert.log \
			gpt-qlat-lanc.log \
			gpt-qlat-rbc-ukqcd-invert.log \
			gpt-qlat-madwf.log

tests_cps = \
			cps-qlat-io-test.log

all: run run-gpt run-cps

run: update-sources $(tests) qlat-ext.log

run-gpt: update-sources $(tests_gpt)

run-cps: update-sources $(tests_cps)

update-sources:
	-for i in $(tests) $(tests_gpt) $(tests_cps) ; do touch "$${i%.log}.py" ; done
	-touch qlat-ext/meson.build

%.log: %.py
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; fn_py="$<" ; fn_json="$${fn_py%.py}.log.json" ; cp -v ../"$$fn_py" ../"$$fn_json" . ; q_verbose=1 timeout -s KILL 30m mpiexec -n 2 $(mpi_options) python3 -m mpi4py ./"$$fn_py" -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep "^CHECK: \|^INFO: " log.full.txt >log.txt ; grep "^CHECK: " log.txt >log.check.txt.new ; grep "^CHECK: " ../"$@" >log.check.txt )
	-if diff "$<".p/log.check.txt "$<".p/log.check.txt.new ; then echo passed ; else tail -n 100 "$<".p/log.full.txt ; echo failed ; fi
	-rm -rfv "$@"
	-touch "$@"
	-mv -v "$<".p/"$@".json.new "$@".json
	-mv -v "$<".p/log.txt "$@"

qlat-ext.log: qlat-ext/meson.build
	-touch qlat-ext/meson.build
	-rm -rfv qlat-ext.p
	-mkdir -p qlat-ext.p
	-time meson setup qlat-ext.p/qlat-ext-build qlat-ext --prefix="$$PWD"/qlat-ext.p/qlat-ext-install -Dpython.platlibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages -Dpython.purelibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages
	-time meson compile -C qlat-ext.p/qlat-ext-build -j$(num_proc)
	-time meson install -C qlat-ext.p/qlat-ext-build
	-time ( cd qlat-ext.p ; PYTHONPATH="$$PWD/qlat-ext-install/lib/python3/dist-packages:$$PYTHONPATH" qlat-ext-install/bin/run.py >log.full.txt 2>&1 ; grep "CHECK: " log.full.txt >log.txt )
	-if diff qlat-ext.p/log.txt "$@" ; then echo passed ; else cat qlat-ext.p/log.full.txt ; false ; fi
	-rm -rfv "$@"
	-touch "$@"
	-mv -v qlat-ext.p/log.txt "$@"

clean:
	-rm -rfv *.p
