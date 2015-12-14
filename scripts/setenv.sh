#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
echo -e "prefix=$prefix\n" | cat - setenv.sh >$prefix/setenv.sh

echo "!!!! $name build !!!!"

rm -rf $temp_dir
