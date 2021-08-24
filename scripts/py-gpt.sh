#!/bin/bash

. conf.sh

name=py-gpt

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/gpt/lib/gpt || true

rsync -av --delete $distfiles/gpt/lib/gpt/ "$prefix"/gpt/lib/gpt/

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
