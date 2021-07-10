#!/bin/bash

echo "Need to run './scripts/mpfr.sh' first"

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
../configure --enable-comms=mpi \
    --enable-simd=GPU \
    --enable-accelerator=cuda \
    --enable-unified=no \
    --with-lime="$prefix" \
    --with-mpfr="$prefix" \
    --prefix=$prefix \
    CXX=nvcc \
    CXXFLAGS="-Xcompiler -fPIC -ccbin mpicxx -gencode arch=compute_70,code=sm_70 -std=c++11" 

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
