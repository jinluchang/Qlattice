#!/usr/bin/env bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ intel = "$USE_COMPILER" ] ; then
    $run icpc "$@"
elif [ gcc = "$USE_COMPILER" ] ; then
    $run g++ "$@"
elif [ clang = "$USE_COMPILER" ] ; then
    $run clang++ "$@"
elif which icpc >/dev/null 2>&1 ; then
    $run icpc "$@"
elif which clang++ >/dev/null 2>&1 ; then
    $run clang++ "$@"
else
    $run g++ "$@"
fi
