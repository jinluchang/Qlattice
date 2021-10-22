#!/bin/bash

. conf.sh

name=dist-update-hash

{

mkdir -p $distfiles

cd $distfiles

sha256sum *.tar.* > sha256sums.txt

echo >> sha256sums.txt

sha256sum python-packages/*.* >> sha256sums.txt

echo >> sha256sums.txt

for fn in * ; do
    if [ -e "$fn"/.git ] ; then
        echo -n "$fn: " >> sha256sums.txt
        ( cd "$fn" ; git rev-parse HEAD >> ../sha256sums.txt )
    fi
done

cat sha256sums.txt

} |& tee $prefix/log.$name.txt
