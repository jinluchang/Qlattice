SHELL=/bin/bash

ifeq ($(strip $(QLAT_PREFIX)),)
	QLAT_PREFIX=$(shell readlink -m "$$HOME/qlat-build/default")
endif

include $(QLAT_PREFIX)/pylib/Makefile.inc
