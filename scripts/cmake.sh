#!/bin/bash

. scripts/conf.sh

name=cmake

{

echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
tar xaf $distfiles/$name-*.tar.*

cd $name-*
./bootstrap \
    --prefix=$prefix \
    --parallel=$num_proc \
    -- \
    -DCMAKE_BUILD_TYPE:STRING=Release
make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
