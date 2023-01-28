#!/bin/bash

# Usage:
# source qcore/set-prefix.sh subdir
# source qcore/conf.sh

# $prefix controls the installation directory

set -e

wd="$(pwd)"
distfiles="$wd/distfiles"

if [ -z "$temp_dir" ] ; then
   temp_dir="$prefix/temp"
fi

src_dir="$temp_dir/src"
build_dir="$temp_dir/build"

if [ -e "$temp_dir" ] ; then
    temp_dir_tmp="$(mktemp -d "$temp_dir".tmp.XXXXX)"
    mv "$temp_dir" "$temp_dir_tmp"
    rm -rf "$temp_dir".tmp.* || true
fi
