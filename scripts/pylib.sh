#!/bin/bash

. conf.sh

name=pylib

echo "!!!! build $name !!!!"

mkdir -pv $prefix/pylib
rm -rfv $prefix/pylib/qlat
rm -rfv $prefix/pylib/rbc_ukqcd_params
rm -rfv $prefix/pylib/auto_contractor
rm -rfv $prefix/pylib/*.py
rm -rfv $prefix/pylib/Makefile.inc
cp -rpv pylib/qlat $prefix/pylib/
cp -rpv pylib/rbc_ukqcd_params $prefix/pylib/
cp -rpv pylib/auto_contractor $prefix/pylib/
cp -pv pylib/*.py $prefix/pylib/
cp -pv pylib/Makefile.inc $prefix/pylib/

mkdir -pv "$prefix"/gpt/lib/gpt || true
rsync -av --delete $distfiles/gpt/lib/gpt/ "$prefix"/gpt/lib/gpt/

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
