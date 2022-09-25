SHELL=/bin/bash

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=

mpi_options=--bind-to none -mca btl_openib_allow_ib true
mpi_options=

mpirun=mpirun --np 2
mpirun=

run:
	q_verbose=1 q_time_limit=2000000 q_num_threads=1 $(mpirun) $(mpi_options) ./run.py $(grid_options) |& tee log.txt

src_path=~/Qlattice-cc/applications/misc/fill-wsnk-prop

diff:
	-diff $(src_path)/run.py run.py
	-diff $(src_path)/Makefile Makefile

sync:
	cp -v $(src_path)/run.py .
	cp -v $(src_path)/Makefile .

clean:
	-rm -rf locks
	# -rm -rf results
