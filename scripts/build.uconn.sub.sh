#!/bin/bash

#SBATCH -n 4           # Job num of cors
#SBATCH -N 1           # Job num of nodes
#SBATCH -t 480         # Job time in mins
#SBATCH -p generalsky  # Job node type

./build.sh uconn

./scripts/qlat-examples-py.sh

./scripts/qlat-examples-cpp.sh

./scripts/qlat-examples-py-gpt.sh
