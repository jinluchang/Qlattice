#!/bin/bash

. conf.sh

name=openmpi

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xaf $distfiles/$name-*

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir

$src_dir/$name-*/configure \
    --prefix=$prefix \
    --with-hwloc=internal
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
