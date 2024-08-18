#!/usr/bin/env bash

if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi

if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicc "$@"
elif which mpiicc >/dev/null 2>&1 ; then
    $run mpiicc -cc=CC.sh "$@"
else
    OMPI_CC=CC.sh $run mpicc "$@"
fi
