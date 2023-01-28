#!/bin/bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicpc "$@"
elif which mpiicpc >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        $run mpiicpc -cxx=g++ "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        $run mpiicpc -cxx=clang++ "$@"
    elif which clang++ >/dev/null 2>&1 ; then
        $run mpiicpc -cxx=clang++ "$@"
    elif which g++ >/dev/null 2>&1 ; then
        $run mpiicpc -cxx=g++ "$@"
    else
        $run mpiicpc "$@"
    fi
elif which mpicxx >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CXX=g++ $run mpicxx "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CXX=clang++ $run mpicxx "$@"
    elif which clang++ >/dev/null 2>&1 ; then
        OMPI_CXX=clang++ $run mpicxx "$@"
    elif which g++ >/dev/null 2>&1 ; then
        OMPI_CXX=g++ $run mpicxx "$@"
    else
        $run mpicxx "$@"
    fi
else
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CXX=g++ $run mpic++ "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CXX=clang++ $run mpic++ "$@"
    elif which clang++ >/dev/null 2>&1 ; then
        OMPI_CXX=clang++ $run mpic++ "$@"
    elif which g++ >/dev/null 2>&1 ; then
        OMPI_CXX=g++ $run mpic++ "$@"
    else
        $run mpic++ "$@"
    fi
fi
