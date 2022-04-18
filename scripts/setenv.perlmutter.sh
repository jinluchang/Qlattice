#!/bin/bash

. scripts/conf.sh

name=setenv

{

echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - scripts/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=16
fi
export PYTHONPATH=
module purge
module load  craype-network-ofi
module load  PrgEnv-gnu
module load  cudatoolkit
export CRAY_ACCEL_TARGET=nvidia80
export CRAY_CPU_TARGET=x86-64
module load  cmake
module list
export QLAT_MPICXX="nvcc -w -std=c++14 -O3 -ccbin CC -Xcompiler -fopenmp -Xcompiler -fno-strict-aliasing --expt-extended-lambda --expt-relaxed-constexpr -gencode arch=compute_80,code=sm_80"
export QLAT_CXXFLAGS="-x cu -Xcompiler -fPIC "
export QLAT_LDFLAGS="-link --shared -lcudart -lcufft "
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
