#!/usr/bin/env bash

name=gcc

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
    if [ -n "$(find-library.sh libgmp.a)" ] ; then
        opts+=" --with-gmp=$(find-library.sh libgmp.a)"
    fi
    if [ -n "$(find-library.sh libmpfr.a)" ] ; then
        opts+=" --with-mpfr=$(find-library.sh libmpfr.a)"
    fi
    if [ -n "$(find-library.sh libmpc.a)" ] ; then
        opts+=" --with-mpc=$(find-library.sh libmpc.a)"
    fi
    if [ -n "$(find-library.sh libisl.a)" ] ; then
        opts+=" --with-isl=$(find-library.sh libisl.a)"
    fi

    time-run $src_dir/$name-*/configure \
        --prefix=$prefix \
        $opts \
        --disable-multilib

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
