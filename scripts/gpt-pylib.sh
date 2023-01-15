#!/bin/bash

. scripts/res/conf.sh

name=gpt-pylib

{

    time {

    echo "!!!! build $name !!!!"

    mkdir -pv "$prefix"/gpt/lib/gpt || true

    rsync -av --delete $distfiles/gpt/lib/gpt/ "$prefix"/gpt/lib/gpt/

    cd $wd
    echo "!!!! $name build !!!!"

    rm -rf $temp_dir || true

}

} 2>&1 | tee $prefix/log.$name.txt
