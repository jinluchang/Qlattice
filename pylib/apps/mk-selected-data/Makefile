SHELL=/bin/bash

src_path=~/Qlattice-cc/pylib/apps/mk-selected-data

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=

mpi_options=--bind-to none -mca btl_openib_allow_ib true
mpi_options=

run:
	q_verbose=1 mpirun --np 2 $(mpi_options) ./run.py --mpi 1.1.1.2 --mpi 1.1.2 $(grid_options) |& tee log.txt
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

diff:
	for fn in * */* ; do if [ -f $(src_path)/$$fn ] ; then echo $$fn ; diff $(src_path)/$$fn $$fn ; fi ; done

sync:
	rsync -av --exclude=".*" $(src_path)/ ./

clean:
	-rm -rf locks
	#-rm -rf results
