#!/bin/bash

. scripts/conf.sh

name=Grid-paboyle

{

    time {

    echo "!!!! build $name !!!!"

    mkdir -p "$prefix"/$name || true

    rsync -av --delete $distfiles/$name/ "$prefix"/$name/

    cd "$prefix/$name"

    INITDIR="$(pwd)"
    rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
    rm -rfv "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

    export CXXFLAGS="$CXXFLAGS -fPIC"

    mkdir build
    cd build
    ../configure \
        --enable-simd=AVX2 \
        --enable-alloc-align=4k \
        --enable-comms=mpi-auto \
        --enable-gparity=no \
        --with-lime="$prefix" \
        --with-fftw="$prefix" \
        --prefix="$prefix/grid-paboyle"

    make -j$num_proc
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} |& tee $prefix/log.$name.txt