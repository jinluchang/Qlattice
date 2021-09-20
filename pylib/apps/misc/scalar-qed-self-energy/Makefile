SHELL=/bin/bash

all:
	time make -C cqlat -j16

run:
	mpirun --np 2 ./run.py

clean:
	-make clean -C cqlat
	-rm -rf cqlat.so
	-rm -rf __pycache__

src_path=~/Qlattice-cc/pylib/apps/misc/scalar-qed-self-energy

sync:
	rsync -av $(src_path)/* .

diff:
	-diff -r $(src_path)/cqlat cqlat
	-diff -r $(src_path)/qlat_ext.py qlat_ext.py
	-diff -r $(src_path)/run.py run.py
	-diff -r $(src_path)/Makefile Makefile
	-diff -r $(src_path)/note.md note.md
