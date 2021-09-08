SHELL=/bin/bash

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=

mpi_options=
mpi_options=--bind-to none -mca btl_openib_allow_ib true

run:
	q_verbose=1 mpirun --np 2 $(mpi_options) ./run.py --mpi 1.1.1.2 --mpi 1.1.2 $(grid_options) |& tee log.txt
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

diff:
	diff ~/Qlattice-cc/pylib/apps/mk-gf-gt/run.py run.py
	diff ~/Qlattice-cc/pylib/apps/mk-gf-gt/Makefile Makefile

sync:
	cp -v ~/Qlattice-cc/pylib/apps/mk-gf-gt/run.py .
	cp -v ~/Qlattice-cc/pylib/apps/mk-gf-gt/Makefile .

clean:
	-rm -rf locks results
