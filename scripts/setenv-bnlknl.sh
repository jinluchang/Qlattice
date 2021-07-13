#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - setenv.sh >"$prefix/setenv.sh" << EOF
prefix="$prefix"

export num_proc=32
export PYTHONPATH=

module purge
module load anaconda3/2019.03-py3.7

# load intel libraries
# source /hpcgpfs01/software/Intel/psxe2019/bin/compilervars.sh -arch intel64
source /hpcgpfs01/software/Intel/psxe2020/bin/compilervars.sh -arch intel64
export INTEL_LICENSE_FILE=/hpcgpfs01/software/Intel/psxe2018.u1/licenses

# module add gcc/9.3.0

EOF

echo "!!!! $name build !!!!"

rm -rf $temp_dir
