#!/bin/bash

. scripts/res/conf.sh

name=Hadrons-tblum

{

    time {

    echo "!!!! build $name !!!!"

    mkdir -p "$prefix"/$name || true

    rsync -av --delete $distfiles/$name/ "$prefix"/$name/

    cd "$prefix/$name"

    if which qlat-include >/dev/null 2>&1 ; then
        for v in $(qlat-include) ; do
            export CPATH="$v":"$CPATH"
        done
        if which organize-colon-list.py >/dev/null 2>&1 ; then
            export CPATH="$(organize-colon-list.py "$CPATH")"
        fi
    fi

    mkdir build

    cd build

    ../configure \
        --with-grid="$prefix/grid-tblum" \
        --prefix="$prefix/hadrons-tblum"

    make -j "$num_proc"
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
