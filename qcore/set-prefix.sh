#!/usr/bin/env bash

# Usage:
# source qcore/set-prefix.sh subdir

# prefix=$prefix/subdir

subdir="$1"

prefix_default="$HOME/qlat-build/default"

if [ -z "$prefix" ] ; then
    prefix="$prefix_default"
fi

if [ -z "$subdir" ] ; then
    export prefix="$prefix"
else
    export prefix="$prefix"/"$subdir"
fi

mkdir -p "$prefix"
