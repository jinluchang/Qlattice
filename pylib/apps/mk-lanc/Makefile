SHELL=/bin/bash

grid_options=--dslash-asm --shm-hugepages --shm 4050
grid_options=

mpi_options=
mpi_options=--bind-to none -mca btl_openib_allow_ib true

run:
	mpirun --np 2 $(mpi_options) ./run.py --mpi 1.1.1.2 $(grid_options)
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

clean:
	-rm -rf locks results
