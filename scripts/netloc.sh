#!/bin/bash

. conf.sh

name=netloc
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xjf $distfiles/$name-*.bz2

cd $name-*
./configure \
    --prefix=$prefix \
    --with-jansson=$prefix
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
