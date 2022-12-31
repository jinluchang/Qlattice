#!/bin/bash

input="$1"
output="$2"

shift 2

cat "$input" >"$output"
echo >>"$output"
for fn in "$@"; do
    echo "include \"$fn\"" >>"$output"
done
