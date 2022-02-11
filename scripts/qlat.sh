#!/bin/bash

./scripts/qlat-header.sh

. conf.sh

name=qlat

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/pylib
rm -rfv $prefix/pylib/cqlat
cp -rpv pylib/cqlat $prefix/pylib/

time make -C $prefix/pylib/cqlat -j $num_proc

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
