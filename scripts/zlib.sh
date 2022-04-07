#!/bin/bash

. scripts/conf.sh

name=zlib

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xaf $distfiles/$name-*

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir

CFLAGS="-O3 -fPIC"

$src_dir/$name-*/configure \
    --prefix=$prefix
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
