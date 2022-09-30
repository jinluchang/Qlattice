#!/bin/bash

. scripts/conf.sh

name=setenv

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix"

cat - scripts/setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"
export prefix="$prefix"
if [ -z "\$num_proc" ] ; then
    num_proc=2
fi
module ()
{
    eval \$(/apps2/Modules/\$MODULE_VERSION/bin/modulecmd bash \$*)
}
module purge
module add modules
module add pre-module
module add post-module
module add vim/8.1
module add git/2.27.0
module add gcc/9.2.0
module add mpi/openmpi/4.0.3
module list
EOF

./scripts/compiler-wrappers.sh

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-build.txt
