#!/bin/bash

# This script would install Qlattice and all its dependencies in ``$prefix''.
# Change $prefix environment variable if you need to install in other directories.
#
# CAUTION! This script could remove files in $prefix silently. Do not put anything
# important there. You have been warned.
#
# Authored by Luchang Jin

echo "Need to run ./scripts/download.sh to download all needed packages."

target="$1"

if [ -z "$target" ] ; then
    target=default
fi

echo "target is $target"

if [ -f "scripts/build.$target.sh" ] ; then
    :
else
    echo "No 'scripts/build.$target.sh' to build '$target'."
    exit
fi

if [ -e "$prefix" ] ; then
    echo "$prefix already exist, continue to build will erase all its contents."
    echo "Use ./scripts/qlat.sh to build Qlat only."
    echo "Ctrl-C to stop."
    for i in {10..0} ; do
        echo -n "$i "
        sleep 1;
    done
    echo
fi

rm -rf "$prefix" || true
mkdir -p "$prefix"

if [ -f "scripts/build.$target.sh" ] ; then
    "scripts/build.$target.sh"
fi

rm -rf $temp_dir
