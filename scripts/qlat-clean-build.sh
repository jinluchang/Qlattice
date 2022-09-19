#!/bin/bash

. scripts/conf.sh

name=qlat-clean-build

{

echo "!!!! build $name !!!!"

rm -rfv "$prefix/build-qlat"*

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
