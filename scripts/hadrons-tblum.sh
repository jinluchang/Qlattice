#!/bin/bash

name=Hadrons-tblum

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    source qcore/conf.sh ..

    mkdir -p "$prefix"/src || true

    debug rsync -av --delete $distfiles/$name/ "$prefix"/src

    cd "$prefix/src"

    if which qlat-include >/dev/null 2>&1 ; then
        for v in $(qlat-include) ; do
            export CPATH="$v":"$CPATH"
        done
    fi

    mkdir build

    cd build

    debug ../configure \
        --with-grid="$prefix/../Grid-tblum" \
        --prefix="$prefix"

    make -j "$num_proc"
    make install

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

} } 2>&1 | tee $prefix/log.$name.txt
