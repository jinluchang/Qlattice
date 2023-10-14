#!/bin/bash

input="$1"
output="$2"

shift 2

cat "$input" >"$output"
echo >>"$output"
for name in "$@"; do
    echo "from .c import Field""$name" >>"$output"
    echo "from .c import SelectedField""$name" >>"$output"
    echo "from .c import SelectedPoints""$name" >>"$output"
done
