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
	for i in $(tests) ; do meson "$$i"-build "$$i" ; ninja -C "$$i"-build | tee -a log-build & done ; time wait

run: all
	-time for i in $(tests) ; do ninja -C "$$i"-build run ; cp "$$i"-build/log "$$i"/log ; done

clean:
	-rm -rfv *-build

clean-logs:
	-time for i in $(tests) ; do rm -v "$$i"/log ; done
