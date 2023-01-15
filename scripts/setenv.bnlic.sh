#!/bin/bash

. scripts/res/conf.sh

name=setenv

mkdir -p "$prefix"

{

echo "!!!! build $name !!!!"

cat - scripts/res/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=8
fi
export PYTHONPATH=
module purge
module add gcc/9.3.0
module add openmpi/3.1.1-gnu
module list
EOF

./scripts/setup-scripts.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name-build.txt
