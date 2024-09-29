#!/usr/bin/env bash

name=setenv-jlabknl

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    num_proc=16
fi
if [ -z "\$temp_dir" ] ; then
    temp_dir=/dev/shm/$(whoami)/temp
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

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
