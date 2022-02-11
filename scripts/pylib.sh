#!/bin/bash

./scripts/qlat-header.sh

. conf.sh

name=pylib

{

echo "!!!! build $name !!!!"

mkdir -pv "$prefix"/gpt/lib/gpt || true
rsync -av --delete $distfiles/gpt/lib/gpt/ "$prefix"/gpt/lib/gpt/

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
