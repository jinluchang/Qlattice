#!/bin/bash

. conf.sh

name=cqlat
echo "!!!! build $name !!!!"

( cd ./pylib/cqlat ; ./update.sh )

rm -rfv $prefix/pylib/cqlat
mkdir -pv $prefix/pylib
cp -rpv pylib/cqlat $prefix/pylib/

time make -C $prefix/pylib/cqlat -j $num_proc

echo "!!!! $name build !!!!"

rm -rf $temp_dir
