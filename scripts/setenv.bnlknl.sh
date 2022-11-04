#!/bin/bash

. scripts/conf.sh

name=setenv

mkdir -p "$prefix"

{

echo "!!!! build $name !!!!"

cat - scripts/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=12
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
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name-build.txt
