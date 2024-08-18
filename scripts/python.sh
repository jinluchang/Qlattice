#!/usr/bin/env bash

name=Python

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

    if [ -f "$prefix/../libffi/lib64/libffi.a" ] ; then
        export LDFLAGS="$LDFLAGS -L$prefix/../libffi/lib64"
        export LIBS="$LIBS -lffi"
    elif [ -f "$prefix/../libffi/lib/libffi.a" ] ; then
        export LDFLAGS="$LDFLAGS -L$prefix/../libffi/lib"
        export LIBS="$LIBS -lffi"
    fi

    time-run $src_dir/$name-*/configure \
        --enable-optimizations --with-lto \
        --enable-shared \
        --prefix=$prefix

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
