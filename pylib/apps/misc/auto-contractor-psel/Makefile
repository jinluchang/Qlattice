SHELL=/bin/bash

grid_options=--dslash-asm --shm-hugepages --shm 4050 --mpi_split 1.1.1.1 --grouped 2
grid_options=

mpi_options=--bind-to none -mca btl_openib_allow_ib true
mpi_options=

run:
	q_verbose=1 mpirun --np 1 $(mpi_options) ./run.py --mpi 1.1.1.1 $(grid_options) |& tee log.txt

src_path=~/Qlattice-cc/pylib/apps/misc/auto-contractor-psel

diff:
	-diff $(src_path)/run.py run.py
	-diff $(src_path)/jobs.py jobs.py
	-diff $(src_path)/cexpr.py cexpr.py
	-diff $(src_path)/Makefile Makefile

sync:
	cp -v $(src_path)/run.py .
	cp -v $(src_path)/jobs.py .
	cp -v $(src_path)/cexpr.py .
	cp -v $(src_path)/Makefile .

clean:
	-rm -rf locks
	#-rm -rf results
