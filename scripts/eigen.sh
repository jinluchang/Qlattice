#!/usr/bin/env bash

name=eigen

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$src_dir"
    cd "$src_dir"
    time-run tar xjf "$distfiles/$name"-*

    time-run rsync -a --delete "$src_dir/$name"-*/{Eigen,signature_of_eigen3_matrix_library,unsupported} "$prefix"/include/

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee "$prefix/log.$name.txt"
