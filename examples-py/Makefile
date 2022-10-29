SHELL=/bin/bash

all: run

tests = auto-cexprs.log \
		auto-cexprs-kpipi.log \
		set-rand.log \
		make-sample-gauge-field.log \
		qar.log \
		lat-io.log \
		field-utils.log \
		qcd-utils.log \
		hmc-pure-gauge.log \
		field-selection.log \
		selected-field.log \
		fields-io.log \
		qlat-id-node.log \
		hmc-pions.log \
		free-invert.log

tests_gpt = gpt-qlat-convert.log \
			gpt-qlat-smear.log \
			gpt-qlat-sample-gauge-field.log \
			gpt-qlat-lanc.log \
			gpt-qlat-free-invert.log \
			gpt-qlat-mdwf.log \
			gpt-qlat-madwf.log \
			gpt-qlat-gfix-mdwf.log \
			gpt-qlat-rbc-ukqcd-invert.log

all: run-all

run-all: run $(tests_gpt)

run: update-sources $(tests) cqlat-ext.log

update-sources:
	-for i in $(tests) $(tests_gpt) ; do touch "$${i%.log}.py" ; done
	-touch cqlat-ext/meson.build

%.log: %.py
	-rm -rf results output_data cache
	-q_verbose=1 mpirun --np 2 ./"$<" --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 >"$@".full 2>&1
	-grep "CHECK: " "$@".full >"$@".new
	-if diff "$@".new "$@" ; then echo passed ; else cat "$@".full ; fi
	-mv -v "$@".new "$@"

cqlat-ext.log: cqlat-ext/meson.build
	-meson cqlat-ext-build cqlat-ext
	-cd cqlat-ext-build ; ninja -j2 ; q_verbose=1 python3 -c "import ctypes; import sys; flags = sys.getdlopenflags(); sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL); import cqlat_ext as ce ; ce.hello_world() ; import qlat as q ; q.timer_display()" >../"$@".full 2>&1
	-grep "CHECK: " "$@".full >"$@".new
	-if diff "$@".new "$@" ; then echo passed ; else cat "$@".full ; fi
	-mv -v "$@".new "$@"

clean:
	-rm -v *.log.full
	-rm -rf results output_data cache cqlat-ext-build
