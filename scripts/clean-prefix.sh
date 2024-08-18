#!/usr/bin/env bash

name=clean-prefix

source qcore/set-prefix.sh

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

    rmdir "$prefix" || true

    if [ -e "$prefix" ] ; then
        echo "Use ./scripts/qlat.sh to build Qlat only."
        echo "$prefix already exist, continue to build will erase all its contents."
        while : ; do
            echo "REALLY Want to continue? Crtl-C to interupt, type 'yes' to continue"
            prompt=""
            read prompt
            if [ "$prompt" == "yes" ] ; then
                break
            fi
        done
        if [ "$prompt" == "yes" ] ; then
            prefix_tmp=$(mktemp -d $prefix.tmp.XXXXX)
            mv "$prefix" "$prefix_tmp"
            rm -rf "$prefix_tmp" || true
        fi
    fi

    mkdir -p "$prefix"

    echo "!!!! $name build !!!!"

} }
