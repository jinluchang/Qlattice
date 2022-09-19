SHELL=/bin/bash

target= \
		mk-gf-gt \
		mk-sel \
		mk-lanc \
		mk-wsrc-prop \
		mk-fsel-self-loop \
		mk-psrc-prop \
		mk-smear-prop \
		mk-topo \
		mk-selected-data

src_path=~/Qlattice-cc/pylib/apps

.ONESHELL:
all:
	for t in $(target)
	do
	rsync -av --exclude=".*" $(src_path)/$${t}/ $${t}/
	make -C $${t} |& tee $${t}.log
	done
