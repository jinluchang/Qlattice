#!/bin/bash

. conf.sh

name=Grid

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/$name || true

rsync -av --delete $distfiles/$name-lehner/ "$prefix"/$name/

cd "$prefix/$name"

git checkout c50f27e68bd4b3e4fb6a1da00aebe224f0a0bc23
echo '-- generating Make.inc files...'
./scripts/filelist
echo '-- generating configure script...'
autoreconf -fvi

INITDIR="$(pwd)"
rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
rm -rfv "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

export CXXFLAGS="$CXXFLAGS -fPIC"

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
    --prefix="$prefix"

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true
