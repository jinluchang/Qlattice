SHELL=/bin/bash

all: run

tests = \
		auto-cexprs.log \
		auto-cexprs-more.log \
		auto-cexprs-kpipi.log \
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
		qplot.log \
		field-utils.log \
		field-selection.log \
		selected-field.log \
		fields-io.log \
		hmc-pure-gauge.log \
		hmc-pions.log \
		hmc-pure-gauge-flowed.log \
		qcd-utils.log \
		free-invert.log

tests_gpt = \
			grid-qlat-io-test.log \
			gpt-qlat-convert.log \
			gpt-qlat-sample-gauge-field.log \
			gpt-qlat-smear.log \
			gpt-qlat-mdwf.log \
			gpt-qlat-gfix-mdwf.log \
			gpt-qlat-free-invert.log \
			gpt-qlat-lanc.log \
			gpt-qlat-rbc-ukqcd-invert.log \
			gpt-qlat-madwf.log \
			gpt-qlat-data-gen-auto.log

tests_cps = \
			cps-qlat-io-test.log

all: run run-gpt run-cps

run: update-sources $(tests) qlat-ext.log

run-gpt: update-sources $(tests_gpt)

run-cps: update-sources $(tests_cps)

update-sources:
	-for i in $(tests) $(tests_gpt) $(tests_cps) ; do touch "$${i%.log}.py" ; done
	-touch qlat-ext/meson.build

%.log: %.py update-sources
	-rm -rfv "$<".p
	-mkdir -p "$<".p
	-time ( cd "$<".p ; cp ../"$<" . ; q_verbose=1 mpiexec -n 2 $(MPI_OPTIONS) ./"$<" -qmp-geom 1 1 1 2 --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >log.full.txt 2>&1 ; grep "CHECK: " log.full.txt >log.txt )
	-if diff "$<".p/log.txt "$@" ; then echo passed ; else cat "$<".p/log.full.txt ; false ; fi
	-mv -v "$<".p/log.txt "$@"

qlat-ext.log: qlat-ext/meson.build update-sources
	-rm -rfv qlat-ext.p
	-mkdir -p qlat-ext.p
	-time meson setup qlat-ext.p/qlat-ext-build qlat-ext --prefix="$$PWD"/qlat-ext.p/qlat-ext-install -Dpython.platlibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages -Dpython.purelibdir="$$PWD"/qlat-ext.p/qlat-ext-install/lib/python3/dist-packages
	-time meson compile -C qlat-ext.p/qlat-ext-build -j$(num_proc)
	-time meson install -C qlat-ext.p/qlat-ext-build
	-time ( cd qlat-ext.p ; PYTHONPATH="$$PWD/qlat-ext-install/lib/python3/dist-packages:$$PYTHONPATH" qlat-ext-install/bin/run.py >log.full.txt 2>&1 ; grep "CHECK: " log.full.txt >log.txt )
	-if diff "$<".p/log.txt "$@" ; then echo passed ; else cat "$<".p/log.full.txt ; false ; fi
	-mv -v "$<".p/log.txt "$@"

clean:
	-rm -rfv *.p
