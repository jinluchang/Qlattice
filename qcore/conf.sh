#!/usr/bin/env bash

# Usage:

# source qcore/conf.sh

# source qcore/set-prefix.sh subdir
# source qcore/conf.sh ..

# source qcore/set-prefix.sh
# source qcore/conf.sh .

# $prefix controls the installation directory
# $wd will be set to the current directory which should be the root directory of the repository

set -e

date

export wd="$(pwd)"

if ! [ -f "$wd/qcore/conf.sh" -a -f "$wd/qcore/set-prefix.sh" ] ; then
    echo "Need to run the scripts in the root directory of the repository."
    echo "Currently, wd='$wd'"
    return 1
fi

export LC_ALL="C"

func() {

local v

for v in "$@" ; do
    if [ -f "$prefix/$v"/setenv.sh ] ; then
        echo "Loading:" "$prefix/$v"/setenv.sh
        source "$prefix/$v"/setenv.sh
        echo "Loaded: " "$prefix/$v"/setenv.sh
    fi
done

}

func "$@"

if "$wd"/qcore/bin/python-check-version.py >/dev/null 2>&1 ; then
    "$wd"/qcore/bin/show-env.py
fi

export

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
