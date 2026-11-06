#!/usr/bin/env bash

name=cmake

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf $src_dir
    mkdir -p $src_dir
    cd $src_dir
    time-run tar xaf $distfiles/$name-*.tar.*

    cd $name-*
    # cmcppdap uses uint32_t in its public header without including its definition.
    # Include it explicitly so builds do not depend on transitive standard headers.
    if ! grep -q '^#include <cstdint>$' Utilities/cmcppdap/include/dap/network.h; then
        sed -i '/^#include <memory>$/a #include <cstdint>' \
            Utilities/cmcppdap/include/dap/network.h
    fi
    time-run ./bootstrap \
        --prefix=$prefix \
        --parallel=$num_proc \
        -- \
        -DCMAKE_BUILD_TYPE:STRING=Release

    time-run make -j$num_proc
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
