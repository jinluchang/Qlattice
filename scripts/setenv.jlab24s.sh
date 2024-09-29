#!/usr/bin/env bash

name=setenv-jlab24s

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    num_proc=32
fi
if [ -z "\$temp_dir" ] ; then
    temp_dir=/dev/shm/$(whoami)/temp
fi
unset PYTHONPATH
unset PYTHONHOME
#
source /etc/profile.d/modules.sh
module use /qcd/dist/el9/modulefiles
module load oneapi/2023.0.0.25537
#
export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER_DUMP=1
export I_MPI_OFI_PROVIDER=psm3
#
module list
#
export TMPDIR=/dev/shm/$(whoami)/tmp
mkdir -p "\$TMPDIR"
if [ -z "\$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
