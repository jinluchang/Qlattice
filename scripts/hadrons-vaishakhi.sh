#!/usr/bin/env bash

name=Hadrons-vaishakhiN

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$src_dir" || true
    time-run rsync -a --delete $distfiles/$name "$src_dir"/
    cd "$src_dir/$name"

    export CXXFLAGS="$CXXFLAGS -fpermissive -DUSE_QLATTICE  -DHADRONS_DEFAULT_LANCZOS_NBASIS=384"

    if python3 -m qlat >/dev/null 2>&1 ; then
        export LD_LIBRARY_PATH="$(python3 -m qlat qlat-config --LD_LIBRARY_PATH)"
    fi

    mkdir build
    cd build
    time-run ../configure \
        --with-grid="$(grid-config --prefix)" \
        --prefix="$prefix"

    time-run make -j "$num_proc"
    time-run make install

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
