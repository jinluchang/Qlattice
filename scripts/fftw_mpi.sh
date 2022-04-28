#!/bin/bash

. scripts/conf.sh

name=fftw

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xzf $distfiles/$name-*.tar.gz

cd $name-*

export CFLAGS="$CFLAGS -fPIC"
export CXXFLAGS="$CXXFLAGS -fPIC"

./configure \
    --prefix=$prefix \
    --enable-mpi \
#    MPICC=cc \
#     --enable-openmp

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name-mpi.txt
