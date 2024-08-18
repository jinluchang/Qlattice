#!/usr/bin/env bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicpc "$@"
elif which mpiicpc >/dev/null 2>&1 ; then
    $run mpiicpc -cxx=CXX.sh "$@"
elif which mpicxx >/dev/null 2>&1 ; then
    OMPI_CXX=CXX.sh $run mpicxx "$@"
else
    OMPI_CXX=CXX.sh $run mpic++ "$@"
fi
