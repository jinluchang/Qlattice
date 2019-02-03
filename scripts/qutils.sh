#!/bin/bash

. conf.sh

name=qutils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/$name
cp -v qutils/* $prefix/include/$name

echo "!!!! $name built !!!!"

rm -rf $temp_dir
