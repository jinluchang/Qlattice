#!/bin/bash

name=setenv.default.sh

source qcore/set-prefix.sh

{ time {

    echo "!!!! build $name !!!!"

    qcore/bin/mk-setenv-dir.py

    echo "!!!! $name build !!!!"

} } 2>&1 | tee $prefix/log.$name.txt
