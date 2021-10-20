#!/bin/bash

. conf.sh

name=python-packages

{

echo "!!!! build $name !!!!"

pip3 install --verbose --no-index -f $distfiles/python-packages Cython
pip3 install --verbose --no-index -f $distfiles/python-packages wheel
pip3 install --verbose --no-index -f $distfiles/python-packages numpy
pip3 install --verbose --no-index -f $distfiles/python-packages mpmath
pip3 install --verbose --no-index -f $distfiles/python-packages sympy

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
