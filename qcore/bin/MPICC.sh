#!/bin/bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicc "$@"
elif which mpiicc >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        $run mpiicc -cc=gcc "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        $run mpiicc -cc=clang "$@"
    elif which clang >/dev/null 2>&1 ; then
        $run mpiicc -cc=clang "$@"
    elif which gcc >/dev/null 2>&1 ; then
        $run mpiicc -cc=gcc "$@"
    else
        $run mpiicc "$@"
    fi
else
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CC=gcc $run mpicc "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CC=clang $run mpicc "$@"
    elif which clang >/dev/null 2>&1 ; then
        OMPI_CC=clang $run mpicc "$@"
    elif which gcc >/dev/null 2>&1 ; then
        OMPI_CC=gcc $run mpicc "$@"
    else
        $run mpicc "$@"
    fi
fi
