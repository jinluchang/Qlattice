#!/usr/bin/env bash

name=setenv-default

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=2
fi
if [ -z "\$temp_dir" ] ; then
    if [ -d /dev/shm ] ; then
        temp_dir=/dev/shm/$(whoami)/temp
    fi
fi
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
