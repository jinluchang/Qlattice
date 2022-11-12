#!/bin/bash

. scripts/conf.sh

name=python-meson-py

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts meson-python

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
