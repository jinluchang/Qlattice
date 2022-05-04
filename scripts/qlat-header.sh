#!/bin/bash

. scripts/conf.sh

name=qlat-header

{

echo "!!!! build $name !!!!"

mkdir -pv $prefix/include
rm -rfv $prefix/include/qutils
rm -rfv $prefix/include/qlat
cp -rpv qutils $prefix/include/
cp -rpv qlat $prefix/include/
cp -rpv qlat/qlat-setup.h $prefix/include/

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

cp pylib/apps/lat-io-glimpse/lat-io-glimpse $prefix/bin/
cp pylib/apps/lat-io-diff/lat-io-diff $prefix/bin/
cp pylib/apps/pickle-glimpse/pickle-glimpse $prefix/bin/
cp pylib/apps/gauge-fix-coulomb/gauge-fix-coulomb $prefix/bin/
cp pylib/apps/qlat-convert/qlat-convert $prefix/bin/
cp pylib/apps/topo-measure/topo-measure $prefix/bin/
cp pylib/apps/eigen-system-repartition/eigen-system-repartition $prefix/bin/
cp pylib/apps/eigen-system-checksum/eigen-system-checksum $prefix/bin/
cp pylib/apps/fields-checksum/fields-checksum $prefix/bin/
cp pylib/apps/crc32/crc32 $prefix/bin/
cp pylib/apps/qar/qar $prefix/bin/

( cd ./pylib/cqlat ; ./update.sh )

echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
