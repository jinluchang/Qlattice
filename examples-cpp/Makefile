SHELL=/bin/bash

tests = qlat-vector-utils \
		qlat-fft-tests \
		qlat-smear-tests \
		qcd-utils-tests \
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
	-time for i in $(tests) ; do make "$$i".compile ; done

run:
	-time for i in $(tests) ; do make "$$i".run ; done

clean:
	-time for i in $(tests) ; do make "$$i".clean ; done

%.compile:
	wd="$@" ; cd "$${wd%.compile}" && meson build
	wd="$@" ; cd "$${wd%.compile}" && time ninja -C build -j3

%.run: %.compile
	-wd="$@" ; cd "$${wd%.run}" && time ninja -C build run
	-wd="$@" ; cd "$${wd%.run}" && if diff build/log log | grep 'CHECK: ' ; then cat build/log.full ; else echo passed ; fi
	-wd="$@" ; cd "$${wd%.run}" && cp build/log log
	-wd="$@" ; cd "$${wd%.run}" && cp build/log.full log.full

%.clean:
	-wd="$@" ; cd "$${wd%.clean}" && rm -rf build log.full
