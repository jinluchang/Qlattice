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

scripts/clean-prefix.sh

. scripts/res/conf.sh

mkdir -p "$prefix"
touch "$prefix"/build-target="$target".txt

scripts/dist-update-hash.sh
scripts/build."$target".sh
