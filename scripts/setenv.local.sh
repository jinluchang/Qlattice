#!/bin/bash

. scripts/conf.sh

name=setenv

{

echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - scripts/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=4
fi
export PYTHONPATH=
module load python3
module list
export CC="mpicc"
export CXX="mpicxx"
export MPICC="mpicc"
export MPICXX="mpicxx"
export USE_COMPILER="gcc"
export QLAT_CC="mpicc -std=c++17 -fopenmp -O2 -xhost -Wall"
export QLAT_CXX="mpicxx -std=c++17 -fopenmp -O2 -xhost -Wall"
export QLAT_MPICC="mpicc -std=c++17 -fopenmp -O2 -xhost -Wall"
export QLAT_MPICXX="mpicxx -std=c++17 -fopenmp -O2 -xhost -Wall "
export QLAT_CXXFLAGS="-fPIC  -I$prefix/include/"
export QLAT_LDFLAGS="--shared"
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
