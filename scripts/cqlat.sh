#!/bin/bash

. conf.sh

name=cqlat
echo "!!!! build $name !!!!"

rm -rf $prefix/pylib/cqlat
mkdir -p $prefix/pylib
cp -rpv pylib/cqlat $prefix/pylib/

make -C $prefix/pylib/cqlat

echo "!!!! $name build !!!!"

rm -rf $temp_dir
