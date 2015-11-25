#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cp -v setenv.sh $prefix

echo "!!!! $name build !!!!"

rm -rf $temp_dir
