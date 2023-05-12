SHELL=/bin/bash

target= \
		mk-gf-gt \
		mk-sel \
		mk-topo \
		mk-lanc \
		mk-wsrc-prop \
		mk-fsel-self-loop \
		mk-psrc-prop \
		mk-smear-prop \
		mk-selected-data \
		auto-contract-meson-corr \
		auto-contract-meson-decay \
		auto-contract-0n2b-sd

src_path=~/Qlattice-cc/applications

.ONESHELL:
all:
	for t in $(target)
	do
	rsync -av --exclude=".*" $(src_path)/$${t}/ $${t}/
	make -C $${t} |& tee $${t}.log
	done
