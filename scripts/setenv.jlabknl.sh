#!/bin/bash

. scripts/res/conf.sh

name=setenv

mkdir -p "$prefix"

{

echo "!!!! build $name !!!!"

cat - scripts/res/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=16
fi
export PYTHONPATH=
module purge
module add gcc-7.1.0
# module add openmpi
# load intel libraries
source /dist/intel/parallel_studio_xe/parallel_studio_xe/psxevars.sh intel64
module list
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
EOF

./scripts/setup-scripts.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name-build.txt
