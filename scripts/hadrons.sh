#!/bin/bash

name=Hadrons

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    mkdir -p "$prefix"/src || true

    time-run rsync -a --delete $distfiles/$name/ "$prefix"/src

    cd "$prefix/src"

    mkdir build

    cd build

    time-run ../configure \
        --with-grid="$(find-library.py libGrid.a)" \
        --prefix="$prefix" \

    time-run make -j "$num_proc"
    time-run make install

    cd $wd

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

} } 2>&1 | tee $prefix/log.$name.txt
