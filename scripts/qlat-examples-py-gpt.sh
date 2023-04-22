#!/bin/bash

name=qlat-examples-py-gpt

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run rsync -a --delete "$wd"/examples-py "$prefix"/

    q_verbose=1 time-run make -C "$prefix"/examples-py run-gpt -j "$num_test" || true

    cd "$wd"

    for log in examples-py/*.log ; do
        echo diff "$prefix/$log" "$log"
        diff "$prefix/$log" "$log" || cat "$prefix/$log".full || true
        cp -rpv "$prefix/$log" "$log" || true
    done

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
