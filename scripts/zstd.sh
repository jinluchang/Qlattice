#!/bin/bash

. scripts/res/conf.sh

name=zstd

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xaf $distfiles/$name-*

rm -rf $build_dir || true
mkdir -p $build_dir || true

cd $src_dir/$name-*/build/meson

meson setup $build_dir \
    -Dbin_programs=true -Dbin_contrib=true \
    --prefix=$prefix

meson compile -C $build_dir
meson install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
