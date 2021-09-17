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
module add gcc
module add gmp
module add hdf5
module add python/3.8.10

module list

EOF

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir
