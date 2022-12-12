#!/bin/bash

mkdir -p "$HOME/qlat-build"

sbatch -o "$HOME/qlat-build"/slurm-%A_%a.out scripts/build.uconn.sub.sh
