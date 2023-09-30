#!/bin/bash

name=qlat-examples-py-cps

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run rsync -a --delete "$wd"/examples-py "$prefix"/

    export mpi_options="--oversubscribe $mpi_options"

    q_verbose=1 time-run make -C "$prefix"/examples-py update-sources || true
    q_verbose=1 time-run make -C "$prefix"/examples-py run-cps -j "$num_test" || true

    cd "$wd"

    for log in examples-py/*.log ; do
        echo diff "$prefix/$log" "$log"
        diff "$prefix/$log" "$log" | grep 'CHECK: ' && ( echo "$log" ; cat "$prefix/$log" || true )
        diff "$prefix/$log" "$log" >/dev/null 2>&1 || ( cp -rpv "$prefix/$log" "$log".new || true )
    done

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
