#!/bin/bash

name=gpt

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    mkdir -p "$prefix"/src || true

    debug rsync -a --delete $distfiles/$name/ "$prefix"/src/

    cd "$prefix/src"

    cd lib/cgpt

    echo "BASIS_SIZE(4)" > lib/basis_size.h
    echo "BASIS_SIZE(10)" >> lib/basis_size.h
    echo "BASIS_SIZE(30)" >> lib/basis_size.h
    echo "BASIS_SIZE(50)" >> lib/basis_size.h

    echo "SPIN(4)" > lib/spin_color.h
    echo "COLOR(3)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,3)" >> lib/spin_color.h

    ./clean

    ./make "$prefix"/../Grid-clehner/src/build "$num_proc"

    rm -rfv "$prefix"/lib
    mkdir -p "$prefix"/lib/python3/dist-packages
    ln -s "$prefix"/src/lib/gpt "$prefix"/lib/python3/dist-packages/gpt
    ln -s "$prefix"/src/lib/cgpt/build/cgpt.so "$prefix"/lib/python3/dist-packages/cgpt.so

    mk-setenv.sh
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
