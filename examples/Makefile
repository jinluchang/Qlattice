tests = lbl-muon-part \
		rng-simple-tests

all:
	for i in $(tests) ; do make -C "$$i" all; done

run:
	for i in $(tests) ; do make -C "$$i" run ; done |& tee log

clean:
	for i in $(tests) ; do make -C "$$i" clean ; done

