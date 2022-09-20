#!/bin/bash

pwd
ls
q_verbose=10 time mpirun --oversubscribe -x OMP_NUM_THREADS=2 --np 2 ./qlat.x | tee log.full
cat log.full | grep -v '^Timer\|^check_status:\|^display_geo_node : id_node =' > log
