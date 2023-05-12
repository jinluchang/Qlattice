SHELL=/bin/bash

all:

run:
	mpirun --np 2 ./run.py

clean:
	-rm -rf __pycache__

sync:
	rsync -av ~/Qlattice-cc/applications/summary-neutrino-hadron-scattering/* .

diff:
	-diff -r ~/Qlattice-cc/applications/summary-neutrino-hadron-scattering/run.py run.py
	-diff -r ~/Qlattice-cc/applications/summary-neutrino-hadron-scattering/Makefile Makefile
	-diff -r ~/Qlattice-cc/applications/summary-neutrino-hadron-scattering/README.md README.md

