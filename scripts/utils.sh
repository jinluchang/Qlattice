#!/bin/bash

. conf.sh

name=utils
echo "!!!! build $name !!!!"

mkdir -p $prefix/include/qlat
cp -v $distfiles/timer.h $prefix/include/qlat
cp -v $distfiles/show.h $prefix/include/qlat
cp -v $distfiles/sha256.h $prefix/include/qlat
cp -v $distfiles/rng-state.h $prefix/include/qlat
cp -v $distfiles/sprng-sha256.h $prefix/include/qlat

echo "!!!! $name build !!!!"

rm -rf $temp_dir
