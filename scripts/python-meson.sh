#!/bin/bash

. scripts/conf.sh

name=python-meson

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts meson

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
