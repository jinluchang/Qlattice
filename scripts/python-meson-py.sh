#!/bin/bash

. scripts/conf.sh

name=python-meson-py

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

opts="--verbose --no-index --no-build-isolation --no-cache-dir -f $distfiles/python-packages"

pip3 install $opts flit_core
pip3 install $opts tomli
pip3 install $opts packaging
pip3 install $opts pep517
pip3 install $opts build
pip3 install $opts pyproject-metadata
pip3 install $opts setuptools_scm
pip3 install $opts scikit-build
pip3 install $opts meson-python

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
