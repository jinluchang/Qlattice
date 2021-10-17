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
module list
export QLAT_CXX="nvcc -w -std=c++14 -O3 -ccbin mpicxx -Xcompiler -fopenmp -Xcompiler -fno-strict-aliasing --expt-extended-lambda --expt-relaxed-constexpr -gencode arch=compute_70,code=sm_70"
export QLAT_CXXFLAGS="-x cu -Xcompiler -fPIC "
export QLAT_LDFLAGS="-link --shared"
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir
