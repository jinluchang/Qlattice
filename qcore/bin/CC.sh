#!/usr/bin/env bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ gcc = "$USE_COMPILER" ] ; then
    $run gcc "$@"
elif [ intel = "$USE_COMPILER" ] ; then
    $run icc "$@"
elif [ clang = "$USE_COMPILER" ] ; then
    $run clang "$@"
elif which icc >/dev/null 2>&1 ; then
    $run icc "$@"
elif which clang >/dev/null 2>&1 ; then
    $run clang "$@"
else
    $run gcc "$@"
fi
