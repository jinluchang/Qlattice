#!/bin/bash

. conf.sh

name=setenv

echo "!!!! build $name !!!!"

mkdir -p "$prefix/bin"

cat - setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=32
fi
export PYTHONPATH=
module purge
module add python/3.8-anaconda-2020-11
module add gcc/9.3.0
module add openmpi/3.1.1-gnu
module list
export CC=CC
export CXX=CXX
export MPICC=MPICC
export MPICXX=MPICXX
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir
