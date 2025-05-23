#!/usr/bin/env bash

# Usage:
# find-header.sh header.h
# return first prefix of the path in CPATH and C_INCLUDE_PATH (endswith /include) that contains the header.h file.

header_name="$1"

value="$CPATH"
if [ -n "$value" ] ; then
    IFS=':' read -a vs <<< "$value"
    for v in "${vs[@]}" ; do
        if [ -f "$v"/"$header_name" ] ; then
            echo "$(dirname "$v")"
            exit 0
        fi
    done
fi

value="$C_INCLUDE_PATH"
if [ -n "$value" ] ; then
    IFS=':' read -a vs <<< "$value"
    for v in "${vs[@]}" ; do
        if [ -f "$v"/"$header_name" ] ; then
            echo "$(dirname "$v")"
            exit 0
        fi
    done
fi
