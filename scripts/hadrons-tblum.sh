#!/bin/bash

. conf.sh

name=Hadrons-tblum

{

echo "!!!! build $name !!!!"

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/$name || true

rsync -av --delete $distfiles/$name/ "$prefix"/$name/

cd "$prefix/$name"

mkdir build

cd build

../configure \
    --with-grid="$prefix/grid-tblum"

make -j "$num_proc"

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
