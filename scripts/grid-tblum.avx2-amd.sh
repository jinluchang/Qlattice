#!/bin/bash

. scripts/res/conf.sh

name=Grid-tblum

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

    export CXXFLAGS="$CXXFLAGS --wrapper-remove-arg=-xcore-avx2 -DUSE_QLATTICE"

    if which qlat-include >/dev/null 2>&1 ; then
        for v in $(qlat-include) ; do
            export CPATH="$v":"$CPATH"
        done
        if which organize-colon-list.py >/dev/null 2>&1 ; then
            export CPATH="$(organize-colon-list.py "$CPATH")"
        fi
    fi

    mkdir build
    cd build
    ../configure \
        --enable-simd=AVX2 \
        --enable-alloc-align=4k \
        --enable-comms=mpi-auto \
        --enable-gparity=no \
        --with-lime="$prefix" \
        --with-fftw="$prefix" \
        --with-hdf5="$prefix" \
        --prefix="$prefix/grid-tblum"

    make -j$num_proc
    make install

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
