#!/bin/bash

. conf.sh

name=setenv

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix"

cat - setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=8
fi
export PYTHONPATH=
CC_OLD="\$CC"
CXX_OLD="\$CXX"
module purge
if [ "\$CXX_OLD" = g++ ] ; then
    module add gcc/7.5.0
else
    module add intel/2019.5.281
    module add hdf5/1.10.5
    module add gsl/2.5
fi
module add cmake/3.14.5
module list
export CC="\$CC_OLD"
export CXX="\$CXX_OLD"
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
