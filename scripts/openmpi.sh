#!/usr/bin/env bash

name=openmpi

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
    if [ -n "$CUDA_ROOT" ] ; then
        opts+=" --with-cuda=${CUDA_ROOT}"
    fi
    if [ -n "$CUDA_LIBDIR" ] ; then
        opts+=" --with-cuda-libdir=${CUDA_LIBDIR}"
    fi

    time-run $src_dir/$name-*/configure \
        --prefix=$prefix \
        $opts \
        --with-hwloc=internal

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
