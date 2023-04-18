#!/bin/bash

pwd
q_verbose=10 OMP_NUM_THREADS=2 time mpiexec $MPI_OPTIONS -n 4 ./qlat.x >log.out 2>log.err || q_verbose=10 OMP_NUM_THREADS=2 time mpiexec $MPI_OPTIONS -n 4 --oversubscribe ./qlat.x >log.out 2>log.err
cat log.out | grep -v '^Timer\|^check_status:\|^display_geo_node : id_node =' >log
cat log.out log.err > log.full
