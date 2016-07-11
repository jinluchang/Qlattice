tests = lbl-muon-part \
		rng-state-tests

all:
	for i in $(tests) ; do make -C "$$i" all; done

run:
	for i in $(tests) ; do make -C "$$i" run ; done 2>&1 | tee log

clean:
	for i in $(tests) ; do make -C "$$i" clean ; done

