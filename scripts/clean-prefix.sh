#!/bin/bash

name=clean-prefix

source qcore/set-prefix.sh

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

    rmdir "$prefix" || true

    if [ -e "$prefix" ] ; then
        echo "$prefix already exist, continue to build will erase all its contents."
        echo "Use ./scripts/qlat.sh to build Qlat only."
        echo "Ctrl-C to stop."
        for i in {10..0} ; do
            echo -n "$i "
            sleep 1;
        done
        echo
        prefix_tmp=$(mktemp -d $prefix.tmp.XXXXX)
        mv "$prefix" "$prefix_tmp"
        rm -rf "$prefix_tmp" || true
    fi

    mkdir -p "$prefix"

    echo "!!!! $name build !!!!"

} }
