#!/usr/bin/env bash

name=qlat-examples-py-cps

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run rsync -a --delete "$wd"/examples-py-cps "$prefix"/

    export mpi_options="--oversubscribe $mpi_options"

    q_verbose=1 time-run make -C "$prefix"/examples-py-cps update-sources || true
    q_verbose=1 time-run make -C "$prefix"/examples-py-cps run-cps -j "$num_test" || true

    cd "$wd"

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee $prefix/log.$name.txt
