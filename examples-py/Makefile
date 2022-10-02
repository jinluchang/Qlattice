SHELL=/bin/bash

all: run

tests = make-sample-gauge-field.log \
		qar.log \
		lat-io.log \
		field-utils.log \
		qcd-utils.log \
		hmc-pions.log \
		hmc-pure-gauge.log \
		field-selection.log \
		selected-field.log \
		fields-io.log \
		free-invert.log \
		qlat-id-node.log

tests_gpt = gpt-qlat-convert.log \
			gpt-qlat-smear.log \
			gpt-qlat-sample-gauge-field.log \
			gpt-qlat-lanc.log \
			gpt-qlat-free-invert.log \
			gpt-qlat-mdwf.log \
			gpt-qlat-madwf.log \
			gpt-qlat-gfix-mdwf.log \
			gpt-qlat-rbc-ukqcd-invert.log

run: $(tests) $(tests_gpt) cqlat-example 

%.log: %.py
	mpirun --np 2 ./"$<" --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 | tee "$@".full
	-rm -rfv results output_data
	-grep "CHECK: " "$@".full > "$@"

cqlat-example: cqlat-ext.log

cqlat-ext.log:
	meson cqlat-ext-build cqlat-ext
	cd cqlat-ext-build ; ninja -j2 ; python3 -c "import ctypes; import sys; flags = sys.getdlopenflags(); sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL); import cqlat_ext as ce ; ce.hello_world() ; import qlat as q ; q.timer_display()" | tee ../"$@".full
	-grep "CHECK: " "$@".full > "$@"

clean:
	-rm -v *.log.full
	-make -C cqlat clean
	-rm -rfv cqlat-ext-build

clean-logs: clean
	-touch $(tests_gpt)
	-rm -v cqlat.log $(tests)

clean-logs-gpt: clean-logs
	-rm -v $(tests_gpt) $(tests)
