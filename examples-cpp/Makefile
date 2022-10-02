SHELL=/bin/bash

tests = qcd-utils-tests \
		dslash-tests \
		lbl-muon-part \
		rng-state-tests \
		field-rng-tests \
		benchmark \
		dist-io \
		selected-field \
		fields-io \
		heatbath \
		propagators \
		hmc \
		simple-1 \
		template

all: run

compile:
	for i in $(tests) ; do make -C "$$i" compile & done ; time wait

run: compile
	-time for i in $(tests) ; do make -C "$$i" run ; done

clean:
	-time for i in $(tests) ; do make -C "$$i" clean ; done

clean-logs: clean
	-time for i in $(tests) ; do make -C "$$i" clean-logs ; done
