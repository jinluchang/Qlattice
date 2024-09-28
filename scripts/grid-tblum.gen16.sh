#!/usr/bin/env bash

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
    if [ -n "$(find-header.sh Eigen/Eigen)" ] ; then
        eigen_path="$(find-header.sh Eigen/Eigen)"
        rsync -av --delete "$eigen_path/include/Eigen/" ${INITDIR}/Grid/Eigen/
        rsync -av --delete "$eigen_path/include/unsupported/Eigen/" ${INITDIR}/Grid/Eigen/unsupported/
        cd ${INITDIR}/Grid
        echo 'eigen_files =\' > ${INITDIR}/Grid/Eigen.inc
        find -L Eigen -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> ${INITDIR}/Grid/Eigen.inc
        cd ${INITDIR}
    else
        ln -vs "${INITDIR}/Eigen/Eigen" "${INITDIR}/Grid/Eigen"
        ln -vs "${INITDIR}/Eigen/unsupported/Eigen" "${INITDIR}/Grid/Eigen/unsupported"
    fi

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

    if python3 -m qlat >/dev/null 2>&1 ; then
        export LIBS="$LIBS $(python3 -m qlat qlat-config --libs)"
        export LDFLAGS="$LDFLAGS $(python3 -m qlat qlat-config --ldflags)"
        export CXXFLAGS="$CXXFLAGS $(python3 -m qlat qlat-config --cxxflags)"
        export LD_LIBRARY_PATH="$(python3 -m qlat qlat-config --LD_LIBRARY_PATH)"
    fi

    mkdir build
    cd build
    time-run ../configure \
        --enable-simd=GEN \
        --enable-gen-simd-width=16 \
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
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee "$prefix/log.$name.txt"
