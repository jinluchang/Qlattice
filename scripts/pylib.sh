#!/bin/bash

. conf.sh

name=pylib

echo "!!!! build $name !!!!"

rm -rfv $prefix/pylib/qlat
rm -rfv $prefix/pylib/rbc_ukqcd_params
rm -rfv $prefix/pylib/auto_contractor
rm -rfv $prefix/pylib/*.py
mkdir -pv $prefix/pylib/qlat
cp -pv pylib/qlat/* $prefix/pylib/qlat
cp -rpv pylib/rbc_ukqcd_params $prefix/pylib/
cp -rpv pylib/auto_contractor $prefix/pylib/
cp -pv pylib/*.py $prefix/pylib/

mkdir -pv "$prefix"/gpt/lib/gpt || true
rsync -av --delete $distfiles/gpt/lib/gpt/ "$prefix"/gpt/lib/gpt/

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
