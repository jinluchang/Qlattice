#!/bin/bash

. conf.sh

name=timer
echo "!!!! build $name !!!!"

mkdir -p $prefix/include
cp -v $distfiles/timer.h $prefix/include/

echo "!!!! $name build !!!!"

rm -rf $temp_dir
