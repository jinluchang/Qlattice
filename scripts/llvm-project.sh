#!/bin/bash

name=llvm-project

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p $src_dir
    cd $src_dir
    debug tar xaf $distfiles/$name-*.xz

    cd $name-*
    mkdir -p build

    debug cmake \
        -S llvm -B build -G Ninja \
        -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;openmp;lld" \
        -DLLVM_ENABLE_RUNTIMES="compiler-rt;libc;libcxx;libcxxabi;libunwind;openmp" \
        -DCMAKE_INSTALL_PREFIX="$prefix" \
        -DCMAKE_PREFIX_PATH="$prefix" \
        -DCMAKE_BUILD_TYPE=Release \
        -DLLVM_ENABLE_FFI=On \
        -DLLVM_PARALLEL_LINK_JOBS=2 \
        -DLLVM_PARALLEL_COMPILE_JOBS="$num_proc" \
        -Wno-dev

    cd build
    ninja
    ninja install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
