SHELL=/bin/bash

src_path=~/Qlattice-cc/applications/auto-contract-meson-decay

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 4
grid_options=--mpi_split 1.1.1.1 --mpi_split 1.1.1 --grouped 4

mpi_options=--bind-to none
mpi_options=--oversubscribe --bind-to none bind-gpu.sh
mpi_options=$(MPI_OPTIONS)

run:
	q_verbose=6 OMP_NUM_THREADS=2 mpiexec -n 4 $(mpi_options) ./run.py $(grid_options) |& tee log.txt
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
	-rm -rf locks __pycache__ log.txt *.new
	-rm -rf cache
	#-rm -rf results
