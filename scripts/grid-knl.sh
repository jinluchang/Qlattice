#!/bin/bash

. conf.sh

name=Grid

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/$name || true

rsync -av --delete $distfiles/$name/ "$prefix"/$name/

cd "$prefix/$name"

INITDIR="$(pwd)"
rm -v "${INITDIR}/Eigen/Eigen/unsupported"
rm -v "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

mkdir build
cd build
../configure \
    --enable-simd=KNL \
    --enable-alloc-align=4k \
    --enable-comms=mpi-auto \
    --enable-mkl \
    --enable-shm=shmget \
    --enable-shmpath=/dev/hugepages \
    --enable-gparity=no \
    --with-lime="$prefix" \
    --with-fftw="$prefix" \
    --prefix="$prefix" \
    CXXFLAGS=-fPIC \
    CXX=icpc \
    MPICXX=mpiicpc

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
