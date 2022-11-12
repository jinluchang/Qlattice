#!/bin/bash

. scripts/conf.sh

name=python-packages

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts psutil
pip3 install $opts numpy
pip3 install $opts sympy
pip3 install $opts scipy

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
