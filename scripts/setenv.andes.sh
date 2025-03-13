#!/usr/bin/env bash

name=setenv-andes

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=8
fi
# module purge
# module load gcc/14.2.0
# module load openmpi/4.1.2
module list
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
