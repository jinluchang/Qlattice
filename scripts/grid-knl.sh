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
../configure \
    --enable-simd=KNL \
    --enable-mkl \
    --enable-alloc-align=2MB \
    --enable-comms=mpi-auto \
    --enable-shm=shmget \
    --enable-shmpath=/dev/hugepages \
    --with-lime="$prefix" \
    --prefix="$prefix" \
    CXXFLAGS=-fPIC \
    CXX=icpc \
    MPICXX=mpiicpc

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
