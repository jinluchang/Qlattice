#!/bin/bash

pwd
q_verbose=10 OMP_NUM_THREADS=2 time mpiexec -n 2 $MPI_OPTIONS ./qlat.x >log.out 2>log.err || q_verbose=10 OMP_NUM_THREADS=2 time mpiexec -n 2 --oversubscribe $MPI_OPTIONS ./qlat.x >log.out 2>log.err
cat log.out | grep -v '^Grid :\|^Timer\|^check_status:\|^display_geometry_node : id_node =' >log
cat log.out log.err > log.full
