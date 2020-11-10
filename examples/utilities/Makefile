SHELL=/bin/bash

tests = crc32 \
		repartition \
		field-compare \
		latio-show

all:
	for i in $(tests) ; do make -C "$$i" & done ; time wait
