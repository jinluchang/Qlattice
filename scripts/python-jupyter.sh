#!/bin/bash

. scripts/conf.sh

name=python-jupyter

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

# opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts --upgrade pip
pip3 install notebook
pip3 install jupyterlab
pip3 install matplotlib
pip3 install lz4

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
