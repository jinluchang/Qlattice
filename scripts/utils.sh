#!/bin/bash

. conf.sh

name=qutils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/$name
cp -v utils/* $prefix/include/$name

echo "!!!! $name built !!!!"

rm -rf $temp_dir
