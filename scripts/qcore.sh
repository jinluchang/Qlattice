#!/bin/bash

source qcore/set-prefix.sh qcore

name=qcore

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

    mkdir -p "$prefix/bin"

    cp -rpv qcore/bin/* "$prefix"/bin/

    echo "!!!! $name build !!!!"

} } 2>&1 | tee "$prefix/log.$name.txt"
