#!/bin/bash

. scripts/conf.sh

name=setenv

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix"
cat - scripts/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=6
fi
export PYTHONPATH=
module purge
module add gcc-7.1.0
# load intel libraries
source /dist/intel/parallel_studio_xe/parallel_studio_xe/psxevars.sh intel64
module list
if [ -z "\$CFLAGS" ] ; then
    export CFLAGS="--wrapper-remove-arg='-cc=gcc' --wrapper-remove-arg='-cc=clang'"
fi
if [ -z "\$CXXFLAGS" ] ; then
    export CXXFLAGS="--wrapper-remove-arg='-cxx=g++' --wrapper-remove-arg='-cxx=clang++'"
fi
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
