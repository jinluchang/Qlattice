#!/usr/bin/env bash

pwd
q_verbose=10 OMP_NUM_THREADS=2 time timeout -s KILL 30m mpiexec -n 4 $mpi_options ./qlat.x >log.out 2>log.err
cat log.out | grep -v '^Grid :\|^Timer\|^check_status:\|^display_geometry_node : id_node =' >log
cat log.out log.err > log.full
