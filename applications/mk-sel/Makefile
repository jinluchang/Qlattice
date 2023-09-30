SHELL=/bin/bash

src_path=~/Qlattice-cc/applications/mk-sel

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=--mpi 1.2.1.2 --mpi 1.2.2 --mpi_split 1.1.1.1 --grouped 4

# mpi_options=--bind-to none -mca btl_openib_allow_ib true
# mpi_options=--oversubscribe --bind-to none bind-gpu.sh
mpi_options=$(mpi_options)

run:
	q_verbose=6 OMP_NUM_THREADS=2 mpiexec -n 4 $(mpi_options) ./run.py $(grid_options) |& tee log.txt
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

diff:
	for fn in * */* ; do if [ -f $(src_path)/$$fn ] ; then echo $$fn ; diff $(src_path)/$$fn $$fn ; fi ; done

sync:
	rsync -av --exclude=".*" $(src_path)/ ./

clean:
	-rm -rf locks __pycache__ log.txt *.new
	-rm -rf cache
	#-rm -rf results
