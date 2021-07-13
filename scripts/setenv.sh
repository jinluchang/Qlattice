#!/bin/bash

. conf.sh

name=setenv
echo "!!!! build $name !!!!"

mkdir -p $prefix
cat - setenv.sh >"$prefix/setenv.sh" << EOF
prefix="$prefix"

EOF

echo "!!!! $name build !!!!"

rm -rf $temp_dir
