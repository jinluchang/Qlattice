#!/bin/bash

. scripts/conf.sh

name=dist-update-hash

{

    ./scripts/update-sources.sh

    (

    distfiles="$wd/distfiles"

    mkdir -p "$distfiles"
    cd "$distfiles"

    sha256sum *.tar.* | sort > sha256sums.txt

    echo >> sha256sums.txt

    sha256sum python-packages/*.* | sort >> sha256sums.txt

    echo >> sha256sums.txt

    for fn in * ; do
        if [ -e "$fn"/.git ] ; then
            echo -n "$fn: "
            ( cd "$fn" ; git rev-parse HEAD )
        fi
    done | sort >> sha256sums.txt

    cat sha256sums.txt

)

} 2>&1 | tee $prefix/log.$name.txt
