#!/bin/bash

echo "Need to run './scripts/mpfr.sh' first"

. conf.sh

name=Grid

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/$name || true

rsync -av --delete $distfiles/$name/ "$prefix"/$name/

cd "$prefix/$name"

INITDIR="$(pwd)"
rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
rm -rfv "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

mkdir build
cd build
../configure \
    --enable-simd=GPU \
    --enable-alloc-align=4k \
    --enable-comms=mpi \
    --enable-unified=no \
    --enable-accelerator=cuda \
    --enable-accelerator-cshift \
    --enable-gparity=no \
    --with-lime="$prefix" \
    --with-fftw="$prefix" \
    --with-mpfr="$prefix" \
    --prefix=$prefix \
    CXX=nvcc \
    CXXFLAGS="-Xcompiler -fPIC -ccbin mpicxx -gencode arch=compute_70,code=sm_70 -std=c++14" \
    LDFLAGS="-Xcompiler -fopenmp"

    # --enable-shm=shmopen \

make -j$num_proc
make install

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir
