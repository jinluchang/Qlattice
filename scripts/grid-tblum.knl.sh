#!/bin/bash

name=Grid-tblum

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$prefix"/src || true

    rsync -av --delete $distfiles/$name/ "$prefix"/src

    cd "$prefix/src"

    INITDIR="$(pwd)"
    rm -rfv "${INITDIR}/Eigen/Eigen/unsupported"
    rm -rfv "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
    ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"

    export CXXFLAGS="$CXXFLAGS -DUSE_QLATTICE"

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
        --enable-simd=KNL \
        --enable-alloc-align=4k \
        --enable-comms=mpi-auto \
        --enable-mkl \
        --enable-shm=shmget \
        --enable-shmpath=/dev/hugepages \
        --enable-gparity=no \
        --prefix="$prefix"

    make -j$num_proc
    make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
