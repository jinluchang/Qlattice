#!/usr/bin/env bash

. scripts/res/conf.sh

name=Grid

{

    time {

    echo "!!!! build $name !!!!"

    mkdir -p "$prefix"/$name || true

    # rsync -av --delete $distfiles/$name/ "$prefix"/$name/
    rsync -av --delete $distfiles/$name-lehner/ "$prefix"/$name/

    cd "$prefix/$name"

    INITDIR="$(pwd)"
    rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
    rm -rfv "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

    export CC=
    export CXX=nvcc
    export CFLAGS=
    export CXXFLAGS="-Xcompiler -fPIC -ccbin mpicxx -gencode arch=compute_70,code=sm_70 -std=c++14"
    export LDFLAGS="-Xcompiler -fopenmp"
    export LIBS=
    export MPICXX=
    export MPICC=

    mkdir build
    cd build
    ../configure \
        --enable-simd=GPU \
        --enable-gen-simd-width=32 \
        --enable-alloc-align=4k \
        --enable-comms=mpi \
        --enable-unified=no \
        --enable-accelerator=cuda \
        --enable-accelerator-cshift \
        --enable-gparity=no \
        --disable-fermion-reps \
        --with-lime="$prefix" \
        --with-fftw="$prefix" \
        --with-gmp="$prefix" \
        --with-mpfr="$prefix" \
        --with-hdf5="$prefix" \
        --with-openssl="$prefix" \
        --prefix=$prefix

    # --enable-shm=nvlink \
    # --enable-setdevice \

    # --enable-shm=shmopen \

    make -j$num_proc
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
