#!/usr/bin/env bash

name=Hadrons

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$src_dir" || true
    time-run rsync -a --delete $distfiles/$name "$src_dir"/
    cd "$src_dir/$name"

    mkdir build
    cd build
    time-run ../configure \
        --with-grid="$(grid-config --prefix)" \
        --prefix="$prefix" \

    time-run make -j "$num_proc"
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
