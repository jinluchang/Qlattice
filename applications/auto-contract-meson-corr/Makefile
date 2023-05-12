SHELL=/bin/bash

src_path=~/Qlattice-cc/applications/auto-contract-meson-corr

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 4
grid_options=--mpi_split 1.1.1.1 --mpi_split 1.1.1 --grouped 4

mpi_options=--bind-to none -mca btl_openib_allow_ib true
mpi_options=--bind-to none

run:
	. ~/qlat-build/default/setenv.sh ; q_verbose=6 q_num_threads=2 q_num_mp_processes=2 mpirun --np 4 --oversubscribe $(mpi_options) ./run.py |& tee log.txt
	-grep '^CHECK: ' log.txt > log.check.txt.new
	-diff log.check.txt.new log.check.txt
	-find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt.new
	-diff md5sum.txt.new md5sum.txt

update-check:
	-mv -v md5sum.txt.new md5sum.txt
	-mv -v log.check.txt.new log.check.txt

diff:
	for fn in * */* ; do if [ -f $(src_path)/$$fn ] ; then echo $$fn ; diff $(src_path)/$$fn $$fn ; fi ; done

sync:
	rsync -av --exclude=".*" $(src_path)/ ./

clean:
	-rm -rf locks cache __pycache__ log.txt *.new
	#-rm -rf results
