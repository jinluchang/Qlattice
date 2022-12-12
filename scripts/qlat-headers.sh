#!/bin/bash

. scripts/conf.sh

name=qlat-headers

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include
rm -rfv $prefix/include/qlat-utils
rm -rfv $prefix/include/qlat
rm -rfv $prefix/include/qlat-setup.h
cp -rpv qlat-utils/include/qlat-utils $prefix/include/
cp -rpv qlat/include/qlat $prefix/include/
cp -rpv qlat/include/qlat/qlat-setup.h $prefix/include/

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
