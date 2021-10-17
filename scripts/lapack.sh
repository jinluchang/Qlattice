#!/bin/bash

. conf.sh

name=lapack
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*

cd $name-*
cmake \
  -DCMAKE_INSTALL_PREFIX:PATH=$prefix \
  -DLAPACKE:BOOL=ON \
  -DLAPACKE_WITH_TMG:BOOL=ON \
  .
make -j$num_proc
make install

cd $wd

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
