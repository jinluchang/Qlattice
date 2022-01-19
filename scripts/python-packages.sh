#!/bin/bash

. conf.sh

name=python-packages

{

echo "!!!! build $name !!!!"

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts --upgrade pip
pip3 install $opts wheel
pip3 install $opts numpy
pip3 install $opts sympy
pip3 install $opts scipy

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
