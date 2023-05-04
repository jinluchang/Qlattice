#!/bin/bash

name=qlat-examples-cpp

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run rsync -a --delete "$wd"/examples-cpp "$prefix"/

    if [ -n "$QLAT_MPICXX" ] ; then
        export MPICXX="$QLAT_MPICXX"
    fi
    export CXX="$MPICXX"
    if [ -n "$QLAT_CXXFLAGS" ] ; then
        export CXXFLAGS="$QLAT_CXXFLAGS"
    fi
    if [ -n "$QLAT_LDFLAGS" ] ; then
        export LDFLAGS="$QLAT_LDFLAGS"
    fi
    if [ -n "$QLAT_LIBS" ] ; then
        export LIBS="$QLAT_LIBS"
    fi

    # q_verbose=1 time-run make -C "$prefix"/examples-cpp compile -j$num_proc || true

    q_verbose=1 time-run make -C "$prefix"/examples-cpp run -j"$num_test" || true

    cd "$wd"

    for log in examples-cpp/*/log ; do
        echo diff "$prefix/$log" "$log"
        diff "$prefix/$log" "$log" | grep 'CHECK: ' && ( echo "$log" ; cat "$prefix/$log" || true )
        diff "$prefix/$log" "$log" >/dev/null 2>&1 || ( cp -rpv "$prefix/$log" "$log".new || true )
    done

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
} } 2>&1 | tee $prefix/log.$name.txt
