#!/bin/bash

name=Hadrons

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    mkdir -p "$prefix"/src || true

    rsync -av --delete $distfiles/$name/ "$prefix"/src

    cd "$prefix/src"

    mkdir build

    cd build

    ../configure \
        --with-grid="$prefix/../Grid" \
        --prefix="$prefix" \

    make -j "$num_proc"
    make install

    cd $wd

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

} } 2>&1 | tee $prefix/log.$name.txt
