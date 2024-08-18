#!/usr/bin/env bash

name=setenv-bnlknl

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=12
fi
export PYTHONPATH=
CC_OLD="\$CC"
CXX_OLD="\$CXX"
module purge
module add gcc/9.3.0
module add intel/psxe2020
# module add anaconda3/2019.03-py3.7
# module add openmpi/1.10.4-gcc
# load intel libraries
# source /hpcgpfs01/software/Intel/psxe2019/bin/compilervars.sh -arch intel64
# source /hpcgpfs01/software/Intel/psxe2020/bin/compilervars.sh -arch intel64
# export INTEL_LICENSE_FILE=/hpcgpfs01/software/Intel/psxe2018.u1/licenses
module list
export CC="\$CC_OLD"
export CXX="\$CXX_OLD"
if [ -z "\$CC" ] ; then
    unset CC
fi
if [ -z "\$CXX" ] ; then
    unset CXX
fi
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
