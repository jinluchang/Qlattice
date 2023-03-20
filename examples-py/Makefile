SHELL=/bin/bash

all: run

tests = \
		qlat-id-node.log \
		set-rand.log \
		make-sample-gauge-field.log \
		qar.log \
		lat-io.log \
		lat-io-test1.log \
		field-utils.log \
		field-selection.log \
		selected-field.log \
		fields-io.log \
		hmc-pure-gauge.log \
		hmc-pions.log \
		hmc-pure-gauge-flowed.log \
		qcd-utils.log \
		auto-cexprs.log \
		auto-cexprs-more.log \
		auto-cexprs-kpipi.log \
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
			gpt-qlat-madwf.log

all: run run-gpt

run: update-sources $(tests) qlat-ext.log

run-gpt: update-sources $(tests_gpt)

update-sources:
	-for i in $(tests) $(tests_gpt) ; do touch "$${i%.log}.py" ; done
	-touch qlat-ext/meson.build

%.log: %.py
	-rm -rf results output_data cache
	-time ( q_verbose=1 mpirun --np 2 ./"$<" --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >"$@".full 2>&1 )
	-grep "CHECK: " "$@".full >"$@".new
	-if diff "$@".new "$@" ; then echo passed ; else cat "$@".full ; fi
	-mv -v "$@".new "$@"

qlat-ext.log: qlat-ext/meson.build
	-meson setup qlat-ext-build qlat-ext --prefix=$$PWD/qlat-ext-install -Dpython.platlibdir=$$PWD/qlat-ext-install/lib/python3/dist-packages -Dpython.purelibdir=$$PWD/qlat-ext-install/lib/python3/dist-packages
	-meson compile -C qlat-ext-build -j2
	-meson install -C qlat-ext-build
	-PYTHONPATH="$$PWD/qlat-ext-install/lib/python3/dist-packages:$$PYTHONPATH" qlat-ext-install/bin/run.py >"$@".full 2>&1
	-grep "CHECK: " "$@".full >"$@".new
	-if diff "$@".new "$@" ; then echo passed ; else cat "$@".full ; fi
	-mv -v "$@".new "$@"

clean:
	-rm -v *.log.full
	-rm -rf results output_data cache qlat-ext-build qlat-ext-install
