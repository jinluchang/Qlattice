SHELL=/bin/bash

all:
	time make -C qcontract -j 16
	mpirun --np 2 ./run.py

clean:
	-make clean -C qcontract
	-rm -rf qcontract.so
	-rm -rf __pycache__
