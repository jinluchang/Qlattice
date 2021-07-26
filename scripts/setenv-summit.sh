#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - setenv.sh >"$prefix/setenv.sh" << EOF
prefix="$prefix"

export num_proc=8
export PYTHONPATH=

module purge

module add DefApps
module add cuda
module add gcc
module add python
module add hdf5

EOF

echo "!!!! $name build !!!!"

rm -rf $temp_dir
