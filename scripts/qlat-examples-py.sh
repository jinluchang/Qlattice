#!/bin/bash

. scripts/conf.sh

name=qlat-examples-py

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/qlat-examples"
    mkdir -p "$build"

    rsync -a --delete "$wd"/examples-py "$build"/

    q_verbose=1 make -C "$build"/examples-py run || true

    cd "$wd"

    for log in *.log ; do
        echo diff "$build/$log" "$log"
        diff "$build/$log" "$log" || cat "$build/$log".full || true
    done

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
