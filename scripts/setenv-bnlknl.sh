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
    num_proc=32
fi

export PYTHONPATH=

module purge

module add anaconda3/2019.03-py3.7
module add gcc/9.3.0

# load intel libraries
# source /hpcgpfs01/software/Intel/psxe2019/bin/compilervars.sh -arch intel64
source /hpcgpfs01/software/Intel/psxe2020/bin/compilervars.sh -arch intel64
export INTEL_LICENSE_FILE=/hpcgpfs01/software/Intel/psxe2018.u1/licenses

EOF

mkdir -p "$prefix/bin"
cat - >"$prefix/bin/mpic++" << EOF
#!/bin/bash
mpiicpc "\$@"
EOF

chmod +x "$prefix/bin/mpic++"

echo "!!!! $name build !!!!"

rm -rf $temp_dir
