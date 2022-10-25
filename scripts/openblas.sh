#!/bin/bash

. scripts/conf.sh

name=OpenBLAS

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xaf $distfiles/$name-*

cd $src_dir/$name-*

make MAKE_NB_JOBS=$num_proc
make install PREFIX=$prefix

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
