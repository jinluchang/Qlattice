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

    # rm -v examples-py/auto-cexprs.log
    make -C examples-py clean-logs

    q_verbose=1 make -C examples-py

    cd "$wd"/examples-py

    for log in *.log ; do
        echo diff "$build/examples-py/$log" "$log"
        diff "$build/examples-py/$log" "$log" || true
    done

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt
