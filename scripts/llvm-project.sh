#!/bin/bash

. conf.sh

name=llvm-project
echo "!!!! build $name !!!!"

mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*.xz

mkdir -p $build_dir
cd $build_dir

cmake -DCMAKE_INSTALL_PREFIX:PATH=$prefix \
    -DLLVM_ENABLE_PROJECTS=all \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    $src_dir/$name-*/llvm

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
