#!/bin/bash

. conf.sh

name=qlat
echo "!!!! build $name !!!!"

rm -rf $prefix/include/qutils
rm -rf $prefix/include/qlat
mkdir -p $prefix/include/qutils
mkdir -p $prefix/include/qlat
cp -pv qutils/*.h $prefix/include/qutils
cp -pv qlat/*.h $prefix/include/qlat
cp -pv qlat-setup.h $prefix/include

echo "!!!! $name build !!!!"

rm -rf $temp_dir
