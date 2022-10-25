#!/bin/bash

. scripts/conf.sh

name=qlat-examples-py

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/qlat-examples"
    mkdir -p "$build"
    cd "$build"

    rsync -av --delete "$wd"/examples-py "$build"/

    make -C examples-py clean-logs

    q_verbose=1 make -C examples-py || true

    cd "$wd"/examples-py

    for log in *.log ; do
        echo diff "$build/examples-py/$log" "$log"
        diff "$build/examples-py/$log" "$log" || cat "$build/examples-py/$log".full
    done

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
