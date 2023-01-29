#!/bin/bash

name=setenv.uconn.sh

source qcore/set-prefix.sh

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=8
fi
source /etc/profile
module purge
module add openmpi/4.1.4
module list
EOF

    #

    "$wd"/qcore/bin/mk-setenv-dir.sh --keep

    echo "!!!! $name build !!!!"

} } 2>&1 | tee $prefix/log.$name.txt
