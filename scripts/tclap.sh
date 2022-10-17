#!/bin/bash

. scripts/conf.sh

name=tclap

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xaf $distfiles/$name-*

cd $src_dir/$name-*

./configure \
    --prefix=$prefix
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
