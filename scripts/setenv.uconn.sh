#!/usr/bin/env bash

name=setenv-uconn

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=8
fi
if [ -z "\$temp_dir" ] ; then
    temp_dir=/dev/shm/$(whoami)/temp
fi
echo "May need to run 'source /etc/profile' if your bashrc does not include it."
# source /etc/profile # cannot run it here as it requires interactive bash.
module purge
module add gcc/12.2.0
module add openmpi/5.0.6
module add cuda/12.3
module list
#
export OMPI_MCA_pml=ucx
#
# ``OMPI_MCA_pml=ucx'' does not work for mpi4py (need to turn of threading)
export MPI4PY_RC_THREADS=0
#
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
