#!/bin/bash

. conf.sh

name=eigen
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir
cmake \
 -DCMAKE_INSTALL_PREFIX=$prefix \
 $src_dir/$name-*
make -j$num_proc
make install

cd $wd

echo "!!!! $name build !!!!"

rm -rf $temp_dir
