#!/bin/bash

. scripts/conf.sh

name=qlat-examples-py

{

    time {

    echo "!!!! build $name !!!!"

    build="$prefix/qlat-examples"
    mkdir -p "$build"

    rsync -a --delete "$wd"/examples-py "$build"/

    if [ -n "$QLAT_MPICXX" ] ; then
        export CXX="$QLAT_MPICXX"
        export MPICXX="$QLAT_MPICXX"
    fi
    if [ -n "$QLAT_CXXFLAGS" ] ; then
        export CXXFLAGS="$QLAT_CXXFLAGS"
    fi
    if [ -n "$QLAT_LDFLAGS" ] ; then
        export LDFLAGS="$QLAT_LDFLAGS"
    fi
    if [ -n "$QLAT_LIBS" ] ; then
        export LIBS="$QLAT_LIBS"
    fi

    time q_verbose=1 make -C "$build"/examples-py run || true

    cd "$wd"

    for log in examples-py/*.log ; do
        echo diff "$build/$log" "$log"
        diff "$build/$log" "$log" || cat "$build/$log".full || true
        cp -rpv "$build/$log" "$log" || true
    done

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
