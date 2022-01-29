SHELL=/bin/bash

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=

mpi_options=
mpi_options=--bind-to none -mca btl_openib_allow_ib true

run:
	. ~/qlat-build/default/setenv.sh ; q_verbose=1 mpirun --np 2 $(mpi_options) ./run.py --mpi 1.1.1.2 --mpi 1.1.2 $(grid_options) |& tee log.txt
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

src_path=~/Qlattice-cc/pylib/apps/mk-gf-gt

diff:
	diff $(src_path)/run.py run.py
	diff $(src_path)/Makefile Makefile

sync:
	cp -v $(src_path)/run.py .
	cp -v $(src_path)/Makefile .

clean:
	-rm -rf locks results
