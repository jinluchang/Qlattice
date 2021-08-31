#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - setenv.sh >"$prefix/setenv.sh" << EOF
prefix="$prefix"

export num_proc=4
export PYTHONPATH=

module purge

module add DefApps
module add cuda
module add gcc
module add gmp
module add hdf5
module load python/3.8.10

EOF

echo "!!!! $name build !!!!"

rm -rf $temp_dir
