#!/bin/bash

name=Grid-clehner

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$src_dir" || true
    time-run rsync -a --delete $distfiles/$name "$src_dir"/
    cd "$src_dir/$name"

    INITDIR="$(pwd)"
    rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
    rm -rfv "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

    export CXXFLAGS="$CXXFLAGS -fPIC -w -Wno-psabi"

    opts=""
    if [ -n "$(find-library.py libgmp.a)" ] ; then
        opts+=" --with-gmp=$(find-library.py libgmp.a)"
    fi
    if [ -n "$(find-library.py libmpfr.a)" ] ; then
        opts+=" --with-mpfr=$(find-library.py libmpfr.a)"
    fi
    if [ -n "$(find-library.py libfftw3.a)" ] ; then
        opts+=" --with-fftw=$(find-library.py libfftw3.a)"
    fi
    if [ -n "$(find-library.py liblime.a)" ] ; then
        opts+=" --with-lime=$(find-library.py liblime.a)"
    fi
    if [ -n "$(find-library.py libcrypto.a)" ] ; then
        opts+=" --with-openssl=$(find-library.py libcrypto.a)"
    fi
    if [ -n "$(find-library.py libhdf5_hl_cpp.a)" ] ; then
        opts+=" --with-hdf5=$(find-library.py libhdf5_hl_cpp.a)"
    fi

    if [ -n "$(find-library.py libz.a)" ] ; then
        LDFLAGS+=" -L$(find-library.py libz.a)/lib"
        CXXFLAGS+=" -I$(find-library.py libz.a)/include"
        export LDFLAGS
        export CXXFLAGS
    fi

    mkdir build
    cd build
    time-run ../configure \
        --enable-simd=AVX2 \
        --enable-alloc-align=4k \
        --enable-comms=mpi-auto \
        --enable-gparity=no \
        $opts \
        --prefix="$prefix"

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee "$prefix/log.$name.txt"
