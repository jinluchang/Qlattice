#!/bin/bash

. scripts/res/conf.sh

name=python-jupyter

{

echo "!!!! build $name !!!!"

find ~/.cache/pip/wheels -type f || true

# rm -rfv ~/.cache/pip/wheels || true

# opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages"

time pip3 install notebook
time pip3 install jupyterlab
time pip3 install jupyterhub
time pip3 install pandas
time pip3 install plotly
time pip3 install seaborn
time pip3 install matplotlib
time pip3 install lz4

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
