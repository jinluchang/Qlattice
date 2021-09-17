#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"

if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=4
fi

export PYTHONPATH=

module purge

module add DefApps
module add cuda/10.1.243
module add gcc/7.5.0
module add gmp
module add hdf5
module add python/3.8.10

module list

export qlat_cxx="nvcc -x cu"
export qlat_cxxld="nvcc -link"
export qlat_flags="-w --shared -std=c++14 -O3 -ccbin mpicxx -gencode arch=compute_70,code=sm_70 -Xcompiler -fPIC -Xcompiler -fopenmp -Xcompiler -fno-strict-aliasing --expt-extended-lambda --expt-relaxed-constexpr"

EOF

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir
