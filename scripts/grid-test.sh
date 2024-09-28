#!/usr/bin/env bash

name=grid-test

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    cd "$prefix"

    if [ -z "$grid_options" ] ; then
        # grid_options="--dslash-asm --shm-hugepages --shm 4050"
        # grid_options="--dslash-asm"
        grid_options=""
    fi

    if [ -z "$geo_options" ] ; then
        geo_options="--grid 16.16.16.16 --Ls 12 --N 10"
    fi

    if [ -z "$OMP_NUM_THREADS" ] ; then
        export OMP_NUM_THREADS=4
    fi

    if [ -z "$mpiexec" ] ; then
        mpiexec=""
    fi

    if [ -z "$mpi_options" ] ; then
        mpi_options="--mpi 1.1.1.1"
    fi

    echo "OMP_NUM_THREADS=\"$OMP_NUM_THREADS\""
    echo "mpiexec=\"$mpiexec\""
    echo "mpi_options=\"$mpi_options\""
    echo "geo_options=\"$geo_options\""
    echo "grid_options=\"$grid_options\""

    grid_prefix="$(grid-config --prefix)"

    time-run $mpiexec $grid_prefix/bin/Benchmark_dwf_fp32 $mpi_options $grid_options $geo_options

    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee $prefix/log.$name.txt
