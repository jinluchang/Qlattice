#!/bin/bash

. conf.sh

name=clean-prefix
echo "!!!! build $name !!!!"

if [ -e "$prefix" ] ; then
    echo "$prefix already exist, continue to build will erase all its contents."
    echo "Use ./scripts/qlat.sh to build Qlat only."
    echo "Ctrl-C to stop."
    for i in {10..0} ; do
        echo -n "$i "
        sleep 1;
    done
    echo
fi

prefix_tmp=$(mktemp -d $prefix.tmp.XXXXX)
mv $prefix $prefix_tmp
rm -rf $prefix_tmp || true

mkdir -p $prefix || true

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
