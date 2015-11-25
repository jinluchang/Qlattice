#!/bin/bash

. setenv.sh

set -e

if [ -z $num_proc ] ; then
    num_proc=16
fi

wd=$(pwd)
distfiles=$wd/distfiles
temp_dir=/dev/shm/$(whoami)/temp/lqps-build
src_dir=$temp_dir/src
build_dir=$temp_dir/build

rm -rf $temp_dir || true
