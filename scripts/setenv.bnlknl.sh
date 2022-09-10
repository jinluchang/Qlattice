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
    num_proc=8
fi
export PYTHONPATH=
CC_OLD="\$CC"
CXX_OLD="\$CXX"
module purge
module add anaconda3/2019.03-py3.7
module add gcc/9.3.0
# load intel libraries
# source /hpcgpfs01/software/Intel/psxe2019/bin/compilervars.sh -arch intel64
source /hpcgpfs01/software/Intel/psxe2020/bin/compilervars.sh -arch intel64
export INTEL_LICENSE_FILE=/hpcgpfs01/software/Intel/psxe2018.u1/licenses
module list
export CC="\$CC_OLD"
export CXX="\$CXX_OLD"
export USE_COMPILER=intel
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
