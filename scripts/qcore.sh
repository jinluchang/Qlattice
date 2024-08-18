#!/usr/bin/env bash

name=qcore

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rsync -av --delete qcore/bin "$prefix"/

    cp -rpv qcore/setenv.sh "$prefix"/

    "$wd"/qcore/bin/mk-setenv.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee "$prefix/log.$name.txt"
