#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p "$prefix"
cat - setenv.sh >"$prefix/setenv.sh" << EOF
echo "Sourcing '$prefix/setenv.sh'"

if [ -z "\$prefix" ] ; then
    prefix="$prefix"
fi
if [ -z "\$num_proc" ] ; then
    num_proc=6
fi

export PYTHONPATH=

module purge

module add anaconda3-latest
module add gcc-7.2.0

# load intel libraries
source /dist/intel/parallel_studio_xe/parallel_studio_xe/psxevars.sh intel64

module list

EOF

mkdir -p "$prefix/bin"
cat - >"$prefix/bin/mpic++" << EOF
#!/bin/bash
mpiicpc "\$@"
EOF

chmod +x "$prefix/bin/mpic++"

. "$prefix/setenv.sh" >"$prefix/log.setenv.txt" 2>&1

echo "!!!! $name build !!!!"

rm -rf $temp_dir
