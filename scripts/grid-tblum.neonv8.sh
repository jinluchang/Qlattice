#!/bin/bash

name=Grid-tblum

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

    export CXXFLAGS="$CXXFLAGS -fPIC -DUSE_QLATTICE -w -Wno-psabi"

    opts=""
    if [ -n "$(find-library.sh libgmp.a)" ] ; then
        opts+=" --with-gmp=$(find-library.sh libgmp.a)"
    fi
    if [ -n "$(find-library.sh libmpfr.a)" ] ; then
        opts+=" --with-mpfr=$(find-library.sh libmpfr.a)"
    fi
    if [ -n "$(find-library.sh libfftw3.a)" ] ; then
        opts+=" --with-fftw=$(find-library.sh libfftw3.a)"
    fi
    if [ -n "$(find-library.sh liblime.a)" ] ; then
        opts+=" --with-lime=$(find-library.sh liblime.a)"
    fi
    if [ -n "$(find-library.sh libcrypto.a)" ] ; then
        opts+=" --with-openssl=$(find-library.sh libcrypto.a)"
    fi
    if [ -n "$(find-library.sh libhdf5_hl_cpp.a)" ] ; then
        opts+=" --with-hdf5=$(find-library.sh libhdf5_hl_cpp.a)"
    fi

    if [ -n "$(find-library.sh libz.a)" ] ; then
        export LDFLAGS="$LDFLAGS -L$(find-library.sh libz.a)/lib"
        export CXXFLAGS="$CXXFLAGS -I$(find-library.sh libz.a)/include"
    fi

    if which qlat-config >/dev/null 2>&1 ; then
        export LIBS="$LIBS $(qlat-config --libs)"
        export LDFLAGS="$LDFLAGS $(qlat-config --ldflags)"
        export CXXFLAGS="$CXXFLAGS $(qlat-config --cxxflags)"
    fi

    mkdir build
    cd build
    time-run ../configure \
        --enable-simd=NEONv8 \
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
