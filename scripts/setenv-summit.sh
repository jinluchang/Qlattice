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
module add cuda/11.4.0
module add gcc
module add gmp
module add hdf5
module load python/3.8.10

EOF

echo "!!!! $name build !!!!"

rm -rf $temp_dir
