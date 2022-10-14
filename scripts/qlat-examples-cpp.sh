#!/bin/bash

. scripts/conf.sh

name=qlat-examples-cpp

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/qlat-examples"
    mkdir -p "$build"
    cd "$build"

    rsync -av --delete "$wd"/examples-cpp "$build"/

    make -C examples-cpp clean-logs

    q_verbose=1 make -C examples-cpp run || true

    cd "$wd"/examples-cpp

    for log in */log ; do
        echo diff "$build/examples-cpp/$log" "$log"
        diff "$build/examples-cpp/$log" "$log" || true
    done

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
