tests = lbl-muon-part \
		rng-state-tests \
		field-rng-tests

all:
	time for i in $(tests) ; do make -C "$$i" all; done

run:
	time for i in $(tests) ; do make -C "$$i" run ; done 2>&1 | tee log

clean:
	time for i in $(tests) ; do make -C "$$i" clean ; done

