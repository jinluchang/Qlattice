#!/bin/bash

pwd
q_verbose=10 time mpirun --oversubscribe -x OMP_NUM_THREADS=2 --np 4 ./qlat.x >log.out 2>log.err
cat log.out | grep -v '^Timer\|^check_status:\|^display_geo_node : id_node =' >log
cat log.out log.err > log.full
