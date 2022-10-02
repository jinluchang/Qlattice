#!/bin/bash

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 720

export q_end_time=$(($(date +%s) + 720 * 60))
export q_num_threads=12
export q_verbose=1

make run
