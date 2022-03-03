#!/bin/bash

. scripts/conf.sh

name=llvm-project

{

echo "!!!! build $name !!!!"

mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*.xz

cd $name-*
mkdir -p build

cmake \
    -S llvm -B build -G Ninja \
    -DLLVM_ENABLE_PROJECTS=all \
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

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-all.txt
