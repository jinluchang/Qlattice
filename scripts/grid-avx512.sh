#!/bin/bash

. conf.sh

name=Grid

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/Grid || true

rsync -av --delete $distfiles/waterret-Grid-*/ "$prefix"/Grid/

cd "$prefix"/Grid

INITDIR="$(pwd)"
rm -v "${INITDIR}/Eigen/Eigen/unsupported"
rm -v "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

mkdir build
cd build
../configure --enable-simd=AVX512 --enable-alloc-align=4k --enable-comms=mpi-auto \
    --with-lime="$prefix" --prefix="$prefix" \
    CXXFLAGS=-fPIC

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
