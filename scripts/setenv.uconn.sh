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
    num_proc=2
fi
module purge
module add modules
module add pre-module
module add post-module
module add vim/8.1
module add gcc/9.2.0
module add mpi/openmpi/4.0.3
module add sqlite/3.30.1
module add libffi/3.2.1
module add tcl/8.6.6.8606
module add bzip2/1.0.6
module add lzma/4.32.7
module add python/3.9.1-sqlite3-rhel7.7-vmaf
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
