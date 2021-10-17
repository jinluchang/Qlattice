#!/bin/bash

. conf.sh

name=perl
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*.tar.*

cd $name-*
./Configure \
    -des \
    -Dprefix=$prefix
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
