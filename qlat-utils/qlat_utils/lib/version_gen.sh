#!/bin/bash

src_file="$1"
output="$2"

src_dir="$(dirname "$src_file")"

wd="$(pwd)"

cd "$src_dir"

if git describe --tags ; then
    version="$(git describe --tags)"
else
    version="v0.41-current" # default version
fi

cd "$wd"

sed "s/VERSION/$version/" "$src_file" >"$output"
