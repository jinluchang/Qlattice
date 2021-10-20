#!/bin/bash

. conf.sh

name=c-lime

{

echo "!!!! build $name !!!!"

rm -rf $src_dir || true
mkdir -p $src_dir || true
cd $src_dir
tar xzf $distfiles/$name.tar.gz

cd *"$name"*
./autogen.sh
cd ..

rm -rf $build_dir || true
mkdir -p $build_dir || true
cd $build_dir

export CFLAGS="$CFLAGS -fPIC"
export CXXFLAGS="$CXXFLAGS -fPIC"

"$src_dir"/*"$name"*/configure \
    --prefix="$prefix"

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
} |& tee $prefix/log.$name.txt
