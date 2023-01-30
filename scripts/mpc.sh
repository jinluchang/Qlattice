#!/bin/bash

name=mpc
source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir || true
    mkdir -p $src_dir || true
    cd $src_dir
    time-run tar xaf $distfiles/$name-*

    rm -rf $build_dir || true
    mkdir -p $build_dir || true
    cd $build_dir

    opts=""
    if [ -n "$(find-library.py libgmp.a)" ] ; then
        opts+=" --with-gmp=$(find-library.py libgmp.a)"
    fi
    if [ -n "$(find-library.py libmpfr.a)" ] ; then
        opts+=" --with-mpfr=$(find-library.py libmpfr.a)"
    fi

    time-run $src_dir/$name-*/configure \
        $opts \
        --prefix=$prefix

    make -j$num_proc
    make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
