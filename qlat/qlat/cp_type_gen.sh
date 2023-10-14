#!/bin/bash

type_name="$1"
input="$2"
output="$3"

sed "s/TYPENAME/$type_name/g" "$input" >"$output"
