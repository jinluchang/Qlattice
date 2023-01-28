#!/bin/bash

name=qcore

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    mkdir -p "$prefix/bin"

    cp -rpv qcore/bin/* "$prefix"/bin/

    cp -rpv qcore/setenv.sh "$prefix"/

    "$wd"/qcore/bin/mk-setenv.sh --keep

    echo "!!!! $name build !!!!"

} } 2>&1 | tee "$prefix/log.$name.txt"
