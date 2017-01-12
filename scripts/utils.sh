#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/utils
cp -v $distfiles/timer.h $prefix/include/utils
cp -v $distfiles/show.h $prefix/include/utils
cp -v $distfiles/sha256.h $prefix/include/utils
cp -v $distfiles/rng-state.h $prefix/include/utils
cp -v $distfiles/sprng-sha256.h $prefix/include/utils

echo "!!!! $name build !!!!"

rm -rf $temp_dir
