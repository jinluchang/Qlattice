#!/bin/bash

name=setenv-default

source qcore/set-prefix.sh

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh

    "$wd"/qcore/bin/mk-setenv-dir.sh

    echo "!!!! $name build !!!!"

} } 2>&1 | tee $prefix/log.$name.txt
