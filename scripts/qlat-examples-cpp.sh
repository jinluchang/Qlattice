#!/bin/bash

. scripts/conf.sh

name=qlat-examples-cpp

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/qlat-examples"
    mkdir -p "$build"

    rsync -a --delete "$wd"/examples-cpp "$build"/

    q_verbose=1 make -C "$build"/examples-cpp compile -j $num_proc || true

    q_verbose=1 make -C "$build"/examples-cpp run || true

    cd "$wd"

    for log in examples-cpp/*/log ; do
        echo diff "$build/$log" "$log"
        diff "$build/$log" "$log" | grep 'CHECK: ' && cat "$build/$log".full || true
        cp -rpv "$build/$log" "$log" || true
    done

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
