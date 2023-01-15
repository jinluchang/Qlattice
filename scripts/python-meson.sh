#!/bin/bash

. scripts/res/conf.sh

name=python-meson

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

time pip3 install $opts meson
time pip3 install $opts ninja

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
