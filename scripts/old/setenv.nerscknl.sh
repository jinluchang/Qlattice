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
module load python3
module list
export CC="CC"
export CXX="CC"
export MPICC="cc"
export MPICXX="CC"
export USE_COMPILER="cc"
export QLAT_CC="CC -std=c++11 -qopenmp -O2 -xhost -Wall"
export QLAT_CXX="CC -std=c++11 -qopenmp -O2 -xhost -Wall"
export QLAT_MPICC="CC -std=c++11 -qopenmp -O2 -xhost -Wall"
export QLAT_MPICXX="CC -std=c++11 -qopenmp -O2 -xhost -Wall"
export QLAT_CXXFLAGS="-fPIC "
export QLAT_LDFLAGS="--shared"
EOF

./scripts/setup-scripts.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name-build.txt
