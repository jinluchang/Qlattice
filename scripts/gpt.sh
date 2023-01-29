#!/bin/bash

name=gpt

source qcore/set-prefix.sh $name

{ time {

    echo "!!!! build $name !!!!"

    mkdir -p "$prefix"/$name || true

    rsync -av --delete $distfiles/$name/ "$prefix"/src/

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

    cd "$wd"

    mk-setenv.sh

    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

} } 2>&1 | tee $prefix/log.$name.txt
