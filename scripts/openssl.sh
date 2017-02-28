#!/bin/bash

. conf.sh

name=openssl
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*.tar.*

cd $name-*
if [ ppc64 = "$(uname -m)" ] ; then
    ./Configure linux-ppc64 \
        --prefix=$prefix \
        shared
else
    ./config \
        --prefix=$prefix \
        shared
fi
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
