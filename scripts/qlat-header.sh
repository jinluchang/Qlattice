#!/bin/bash

. conf.sh

name=qlat-header

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include
rm -rfv $prefix/include/qutils
rm -rfv $prefix/include/qlat
cp -rpv qutils $prefix/include/
cp -rpv qlat $prefix/include/
cp -pv qlat-setup.h $prefix/include

mkdir -pv $prefix/pylib
rm -rfv $prefix/pylib/qlat
rm -rfv $prefix/pylib/rbc_ukqcd_params
rm -rfv $prefix/pylib/auto_contractor
rm -rfv $prefix/pylib/*.py
cp -rpv pylib/qlat $prefix/pylib/
cp -rpv pylib/rbc_ukqcd_params $prefix/pylib/
cp -rpv pylib/auto_contractor $prefix/pylib/
cp -pv pylib/*.py $prefix/pylib/
cp -pv pylib/Makefile.inc $prefix/pylib/

( cd ./pylib/cqlat ; ./update.sh )

echo "!!!! $name build !!!!"

rm -rf $temp_dir
