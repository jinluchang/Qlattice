#!/bin/bash

. scripts/conf.sh

name=qlat-packages

{

    time {

    echo "!!!! build $name !!!!"

    find ~/.cache/pip/wheels -type f || true

    build="$prefix/qlat-packages"
    rm -rfv "$build" || true
    mkdir -p "$build"

    cd "$build"

    # opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages -f $build"
    opts="--verbose -f $build"

    pip3 uninstall -y qlat-utils qlat qlat-grid || true

    python3 -m build -ns -o "$build" "$wd"/qlat-utils
    pip3 install $opts qlat-utils
    python3 -m build -ns -o "$build" "$wd"/qlat
    pip3 install $opts qlat
    if grid-config --prefix >/dev/null 2>&1 ; then
        python3 -m build -ns -o "$build" "$wd"/qlat-grid
        pip3 install $opts qlat-grid
    fi

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
