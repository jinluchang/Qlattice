#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include
cp -v $distfiles/timer.h $prefix/include/
cp -v $distfiles/show.h $prefix/include/
cp -v $distfiles/sha256.h $prefix/include/
cp -v $distfiles/rng-state.h $prefix/include/
cp -v $distfiles/sprng-sha256.h $prefix/include/

echo "!!!! $name build !!!!"

rm -rf $temp_dir
