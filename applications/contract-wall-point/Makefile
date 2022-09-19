SHELL=/bin/bash

all:
	time make -C cqlat -j16

run:
	mpirun --np 2 ./run.py

clean:
	-make clean -C cqlat
	-rm -rf cqlat.so
	-rm -rf __pycache__

sync:
	rsync -av ~/Qlattice-cc/pylib/apps/contract-wall-point/* .

diff:
	-diff -r ~/Qlattice-cc/pylib/apps/contract-wall-point/cqlat cqlat
	-diff -r ~/Qlattice-cc/pylib/apps/contract-wall-point/run.py run.py
	-diff -r ~/Qlattice-cc/pylib/apps/contract-wall-point/qlat_ext.py qlat_ext.py
	-diff -r ~/Qlattice-cc/pylib/apps/contract-wall-point/Makefile Makefile
	-diff -r ~/Qlattice-cc/pylib/apps/contract-wall-point/README.md README.md
