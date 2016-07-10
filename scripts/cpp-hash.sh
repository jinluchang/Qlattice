#!/bin/bash

. conf.sh

name=cpp-hash
echo "!!!! build $name !!!!"

rm -rf $src_dir
mkdir -p $src_dir
cd $src_dir
rsync -avz $distfiles/$name .

cd $name
prefix=$prefix make

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
