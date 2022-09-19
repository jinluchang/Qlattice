SHELL=/bin/bash

tests = crc32 \
		repartition \
		field-compare \
		gf-info \
		latio-show

all:
	for i in $(tests) ; do make -C "$$i" & done ; time wait

clean:
	for i in $(tests) ; do make -C "$$i" clean ; done ; time wait
