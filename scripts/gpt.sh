#!/bin/bash

. conf.sh

name=gpt

echo "!!!! build $name !!!!"

rm -rf "$prefix"/gpt || true
mkdir -p "$prefix" || true

cd "$prefix"
tar xaf $distfiles/gpt.tar.gz

mv waterret-gpt-* gpt
cd gpt/lib/cgpt
./make "$prefix"/Grid/build "$num_proc"

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
