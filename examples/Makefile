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

all:
	for i in $(tests) ; do make -C "$$i" all & done ; time wait

diff:
	-git diff

run: all
	-time for i in $(tests) ; do make -C "$$i" run-test ; done
	make diff

clean:
	-time for i in $(tests) ; do make -C "$$i" clean ; done

dist-clean:
	-time for i in $(tests) ; do make -C "$$i" dist-clean ; done

clean-logs:
	-time for i in $(tests) ; do rm -v "$$i"/log ; done
