#!/bin/bash

. scripts/res/conf.sh

name=Hadrons

{

    time {

    echo "!!!! build $name !!!!"

    mkdir -p "$prefix"/$name || true

    rsync -av --delete $distfiles/$name/ "$prefix"/$name/

    cd "$prefix/$name"

    mkdir build

    cd build

    ../configure \
        --with-grid="$prefix/grid-paboyle"

    make -j "$num_proc"

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
