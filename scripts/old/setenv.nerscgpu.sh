#!/usr/bin/env bash

. scripts/res/conf.sh

name=setenv

mkdir -p $prefix

{

echo "!!!! build $name !!!!"

cat - scripts/res/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=4
fi
export PYTHONPATH=
module purge
module load cgpu
module load gcc/8.3.0
module load openmpi/4.0.3
module load cuda/11.1.1
module load python/3.8-anaconda-2020.11
module list
export QLAT_MPICXX="nvcc -w -std=c++14 -O3 -ccbin mpicxx -Xcompiler -fopenmp -Xcompiler -fno-strict-aliasing --expt-extended-lambda --expt-relaxed-constexpr -gencode arch=compute_70,code=sm_70"
export QLAT_CXXFLAGS="-x cu -Xcompiler -fPIC"
export QLAT_LDFLAGS="-link --shared -lcudart -lcufft "
EOF

./scripts/setup-scripts.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name-build.txt
