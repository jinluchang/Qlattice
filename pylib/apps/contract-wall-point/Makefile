SHELL=/bin/bash

all:
	time make -C cqlat -j16

run:
	mpirun --np 2 ./run.py

clean:
	-make clean -C cqlat
	-rm -rf cqlat.so
	-rm -rf __pycache__
