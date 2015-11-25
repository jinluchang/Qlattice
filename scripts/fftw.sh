#!/bin/bash

. conf.sh

name=fftw

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xzf $distfiles/$name-*.tar.gz

cd $name-*

if [ $arch = bgq ] ; then
    ./configure \
        --prefix=$prefix \
        --enable-openmp
else
    ./configure \
        --prefix=$prefix \
        --enable-openmp
fi

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
