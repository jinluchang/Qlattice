#!/usr/bin/env bash

name=gpt

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$src_dir" || true
    time-run rsync -a --delete "$distfiles/$name" "$src_dir"
    cd "$src_dir/$name"

    cd lib/cgpt

    echo "BASIS_SIZE(4)" > lib/basis_size.h
    echo "BASIS_SIZE(10)" >> lib/basis_size.h
    echo "BASIS_SIZE(25)" >> lib/basis_size.h

    echo "SPIN(4)" > lib/spin_color.h
    echo "COLOR(3)" >> lib/spin_color.h
    echo "COLOR(2)" >> lib/spin_color.h
    echo "COLOR(1)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,3)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,2)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,1)" >> lib/spin_color.h

    time-run ./clean

    time-run ./make %grid-config "$num_proc"

    rm -rfv "$prefix"/lib
    rm -rfv "$prefix"/src
    mkdir -pv "$prefix"/lib/python3/dist-packages
    mkdir -pv "$prefix"/src
    time-run rsync -a --delete "$src_dir/$name"/lib/gpt "$prefix"/lib/python3/dist-packages/
    time-run rsync -a --delete "$src_dir/$name"/lib/cgpt/build/cgpt.so "$prefix"/lib/python3/dist-packages/
    time-run rsync -a --delete "$src_dir/$name"/tests "$prefix"/src/
    time-run rsync -a --delete "$src_dir/$name"/applications "$prefix"/src/
    time-run rsync -a --delete "$src_dir/$name"/benchmarks "$prefix"/src/
    time-run rsync -a --delete "$src_dir/$name"/documentation "$prefix"/src/
    time-run rsync -a --delete "$src_dir/$name"/docker "$prefix"/src/

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
