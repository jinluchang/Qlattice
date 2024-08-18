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
    if which nvcc >/dev/null 2>&1 ; then
        NVCC_PATH=$(which nvcc)
        CUDA_ROOT=${NVCC_PATH%/bin/nvcc}
        opts+=" --with-cuda=${CUDA_ROOT}"
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
