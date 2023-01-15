#!/bin/bash

. scripts/res/conf.sh

name=lapack

{

echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir

cmake $src_dir/$name-* \
  -DCMAKE_INSTALL_PREFIX=$prefix \
  -DLAPACKE=ON \
  -DLAPACKE_WITH_TMG=ON \
  .
make -j$num_proc
make install

cd $wd

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
